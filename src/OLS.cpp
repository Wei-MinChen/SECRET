////////////////////////////////////////////////////////////// ////////
// Linear.cpp
// (c) 2023-2026 Wei-Min Chen
//
// This file is distributed as part of the SECRET source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile SECRET.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Oct 26, 2026

#include <time.h>
#include <math.h>
#include <string.h>
#include "Random.h"
#include "MathSVD.h"
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#include "Penalized.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

extern void localPCA(Matrix &risks, int dimPC, Vector &Values, Matrix &PCs, int coreCount);
extern int plotRSQ(const char *input, const char *traitName, const char *rpath, String &msg);
extern int plotBestHeatmap(const char *input, const char *rpath, String &msg);
extern int plotPGSCorr(const char *input, const char *statfile, const char *traitName, const char *rpath, String &msg);
extern int plotPGSHeatmap(const char *input, const char *rpath, String &msg);
extern int plotCalibration(const char *input, const char *rpath, String &msg);
extern void WritePGSCorrelations(const char *outfile0, Matrix &risks, StringArray &rsNames, IntArray &topIndex,
  IntArray &samples, IntArray &status, int coreCount, bool &caseControlDone, double pruneR2);
extern void ScreenLD1(Matrix &risks, double ldT, IntArray pairList[], int coreCount);
extern void ScreenLD2(Matrix &risks, double ldT, IntArray pairList[], Vector pairLD[], Vector &rsSD, int coreCount);
extern void MakeBinaryLD(int rsCount, IntArray pairList[], unsigned long long int **LD, int coreCount);
extern void MakeBinaryLD(Matrix &risks, double ldT, unsigned long long int **LD, Vector &rsSD, int coreCount);
extern void MakePrunedOrder(unsigned long long int **LD, IntArray &order);
extern void LDSummary(unsigned long long int **LD, int rsCount, int wordCount, int coreCount);

// allY  (idCount x traitCount): the phenotypes (quantitative traits first, then affections)
// allX  (idCount x covariateCount): the raw covariate columns, one row per sample.  Missing
//       values stay _NAN_; --logistic uses them (see MakeVariable), --linear does not.
// masks: bit t of masks[i] is set when sample i has a usable phenotype for trait t
void MakePhenos(Matrix &risks, Matrix &alltraits, StringArray &alltraitNames, IntArray &affections,
  IntArray &quantraits, IntArray &covariates, Matrix &allY, Matrix &allX,
  unsigned long long int **masks, int coreCount)
{
  int affectionCount = affections.Length();
  int quanCount = quantraits.Length();
  int idCount = alltraits.rows;
  int covariateCount = covariates.Length();
  allY.Zero();
  allX.Dimension(idCount, covariateCount);
  for(int i = 0; i < idCount; i++)
    for(int c = 0; c < covariateCount; c++)
      allX[i][c] = alltraits[i][covariates[c]];
  int twordCount = (affectionCount+quanCount-1)/64+1;
  for(int i = 0; i < idCount; i++)
    for(int w = 0; w < twordCount; w++)
      masks[i][w] = 0;
  if(quanCount){
    if (covariateCount)
      AdjustCovariate(risks, alltraits, alltraitNames, quantraits, covariates, allY, masks, coreCount);
    else
      AdjustNoCovariate(risks, alltraits, alltraitNames, quantraits, allY, masks, coreCount);
  }
  if(affectionCount){
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for(int i = 0; i < idCount; i++)
      for (int a = 0; a < affectionCount; a++) {
        int trait = affections[a];
        if(alltraits[i][trait] != _NAN_ && alltraits[i][trait] != 0 && risks[0][i] != _NAN_){
          allY[i][quanCount+a] = alltraits[i][trait];
          masks[i][(quanCount+a)>>6] |= ((unsigned long long int)0x1<<((quanCount+a)&0x3F));
        }
      } // End of a loop
  }
}

void popcountbycol(unsigned long long int **table, int idCount, IntArray &counts)
{
  int colCount = counts.Length();
  if(colCount==1){
    int sum = 0;
    for(int i = 0; i < idCount; i++)
      sum += table[i][0];
    counts[0] = sum;
    return;
  }
  int wordCount = (colCount-1)/64+1;
  const int IBLOCKSIZE=255;
  const int WBLOCKSIZE=32;
  unsigned long long int words[256];
  counts.Zero();
  for(int w = 0; w < wordCount; w += WBLOCKSIZE){
    int wMax = w < wordCount-WBLOCKSIZE? w+WBLOCKSIZE: wordCount;
    for(int i = 0; i < idCount; i += IBLOCKSIZE){
      int iMax = i<idCount-IBLOCKSIZE? i+IBLOCKSIZE: idCount;
      for(int s = 0; s < 8; s++)
        words[s] = 0;
      for(int ii = i; ii < iMax; ii++){
        int base = 0;
        for(int ww = w; ww < wMax; ww++, base+=8){
          unsigned long long int word = table[ii][ww];
          for(int s = 0; s < 8; s++, base++,word>>=1)
            words[base] += (word&0x0101010101010101);
        } // End of ww loop
      } // End of ii loop
      for(int ww = w; ww < wMax; ww++)
        for(int s = 0; s < 8; s++){
          unsigned long long int word = words[(ww-w)*8+s];
          for(int b = 0; b < 8; b++, word>>=8)
            counts[ww*64+b*8+s] += (word&0xFF);
        } // End of s loop
    } // End of i loop
  } // End of w loop
}

void MatrixLinearScan(Matrix &risks, Matrix &allY, unsigned long long int **masks,
  Matrix &GT, Matrix &sigmaG, Matrix &sigmaT, IntArray &SS, int coreCount)
{
  int idCount = allY.rows;
  int traitCount = allY.cols;
  int rsCount = risks.rows;
  Matrix GY(rsCount, traitCount), sumG(rsCount, traitCount), GG(rsCount, traitCount);
  const int TBLOCK=64;
  const int IBLOCK=256;
  const int RBLOCK=16;
  double V[3][RBLOCK][TBLOCK]; // 24KB
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(V)
#endif
  for(int r = 0; r < rsCount; r+=RBLOCK){
    int rMax = r<rsCount-RBLOCK? RBLOCK: rsCount-r;
    for(int t = 0; t < traitCount; t+=TBLOCK){
      int tMax = t < traitCount-TBLOCK? TBLOCK: traitCount-t;
      int t64 = (t>>6);
      memset(&V[0][0][0], 0, RBLOCK*TBLOCK*3*sizeof(double));
      for(int i = 0; i < idCount; i+=IBLOCK){ // CPU CACHE BLOCK
        int iMax = i < idCount-IBLOCK? i+IBLOCK: idCount;
        for(int rr = 0; rr < rMax; rr++)
          for(int ii = i; ii < iMax; ii++){
            double risk = risks[r+rr][ii];
            double risk2 = risk * risk;
            unsigned long long int word = masks[ii][t64];
            for(int tt = 0; tt < tMax; tt++, word >>= 1){
              int flag = (word & 0x1);
              V[0][rr][tt] += risk * allY[ii][t+tt];
              V[1][rr][tt] += risk * flag;
              V[2][rr][tt] += risk2 * flag;
            } // End of tt loop
          } // End of ii loop
      } // End of i loop CPU CACHE BLOCK
      int memSize = tMax * sizeof(double);
      for(int rr = 0; rr < rMax; rr++){
        memcpy(&GY[r+rr][t], &V[0][rr][0], memSize);
        memcpy(&sumG[r+rr][t], &V[1][rr][0], memSize);
        memcpy(&GG[r+rr][t], &V[2][rr][0], memSize);
      } // End of rr loop
    } // End of t loop
  } // End of OMP r loop
  popcountbycol(masks, idCount, SS);
  double sumY, YY;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(sumY, YY)
#endif
  for(int t = 0; t < traitCount; t++){
    sumY = YY = 0.0;
    for(int i = 0; i < idCount; i++){
      double y = allY[i][t];
      sumY += y;
      YY += y*y;
    } // End of i loop
    double st = (YY - sumY * sumY / SS[t]) / SS[t];
    for (int r = 0; r < rsCount; r++) {
      GT[t][r] = GY[r][t] - sumG[r][t] * sumY / SS[t];
      sigmaG[t][r] = (GG[r][t] - sumG[r][t] * sumG[r][t] / SS[t]) / SS[t];
      sigmaT[t][r] = st;
    } // End of r loop
  } // End of OMP t loop
}

void MatrixLinearScan_Sum(Matrix &risks, Matrix &allY, unsigned long long int **masks,
  Matrix &GY, Matrix &GG, Matrix &G, Vector &YY, Vector &Y, IntArray &SS, int coreCount)
{
  int idCount = allY.rows;
  int traitCount = allY.cols;
  int rsCount = risks.rows;
  Matrix gy(rsCount, traitCount), sumg(rsCount, traitCount), gg(rsCount, traitCount);
  const int TBLOCK=64;
  const int IBLOCK=256;
  const int RBLOCK=16;
  double V[3][RBLOCK][TBLOCK]; // 24KB
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(V)
#endif
  for(int r = 0; r < rsCount; r+=RBLOCK){
    int rMax = r<rsCount-RBLOCK? RBLOCK: rsCount-r;
    for(int t = 0; t < traitCount; t+=TBLOCK){
      int tMax = t < traitCount-TBLOCK? TBLOCK: traitCount-t;
      int t64 = (t>>6);
      memset(&V[0][0][0], 0, RBLOCK*TBLOCK*3*sizeof(double));
      for(int i = 0; i < idCount; i+=IBLOCK){ // CPU CACHE BLOCK
        int iMax = i < idCount-IBLOCK? i+IBLOCK: idCount;
        for(int rr = 0; rr < rMax; rr++)
          for(int ii = i; ii < iMax; ii++){
            double risk = risks[r+rr][ii];
            double risk2 = risk * risk;
            unsigned long long int word = masks[ii][t64];
            for(int tt = 0; tt < tMax; tt++, word >>= 1){
              int flag = (word & 0x1);
              V[0][rr][tt] += risk * allY[ii][t+tt];
              V[1][rr][tt] += risk * flag;
              V[2][rr][tt] += risk2 * flag;
            } // End of tt loop
          } // End of ii loop
      } // End of i loop CPU CACHE BLOCK
      int memSize = tMax * sizeof(double);
      for(int rr = 0; rr < rMax; rr++){
        memcpy(&gy[r+rr][t], &V[0][rr][0], memSize);
        memcpy(&sumg[r+rr][t], &V[1][rr][0], memSize);
        memcpy(&gg[r+rr][t], &V[2][rr][0], memSize);
      } // End of rr loop
    } // End of t loop
  } // End of OMP r loop
  popcountbycol(masks, idCount, SS);
  double sumy, yy;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(sumy, yy)
#endif
  for(int t = 0; t < traitCount; t++){
    sumy = yy = 0.0;
    for(int i = 0; i < idCount; i++){
      double y = allY[i][t];
      sumy += y;
      yy += y*y;
    } // End of i loop
    YY[t] = yy;
    Y[t] = sumy;
    for (int r = 0; r < rsCount; r++) {
      GY[t][r] = gy[r][t];
      GG[t][r] = gg[r][t];
      G[t][r] = sumg[r][t];
    } // End of r loop
  } // End of OMP t loop
}

double GenomicControlNumber(Vector &Z, IntArray &sortedrsIndex)
{
    int rsCount = Z.Length();
    IntArray validrs(0);
    Vector oneChisq(0);
    for (int r = 0; r < rsCount; r++) {
        if (Z[r] == _NAN_) continue;
        oneChisq.Push(-Z[r] * Z[r]);
        validrs.Push(r);
    }   // End of r loop
    int validCount = validrs.Length();
    if(validCount==0) return _NAN_;
    QuickIndex idx;
    idx.Index(oneChisq);
    sortedrsIndex.Dimension(validCount);
    for(int r = 0; r < validCount; r++)
      sortedrsIndex[r] = validrs[idx[r]];
    return -oneChisq[idx[validCount / 2]] / 0.4549364;
}

void WriteLinear(String &outfile, IntArray &sortedrsIndex, StringArray &rsName,
  StringArray &annotation, Vector &GT, Vector &sigmaG, Vector &sigmaT, int ss)
{
  int rsCount = sortedrsIndex.Length();
  if (rsCount == 0) return;
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  int annotCount = annotation.Length();
  fprintf(fp, "PGS\tN\tBetaS\tBetaT\tZ\tPvalue\tRsq");
  WriteCoverageHeader(fp);
  if (annotCount)
        fprintf(fp, "\tAnnotation");
  for (int r = 0; r < rsCount; r++) {
        int rs = sortedrsIndex[r];
        double statistic = GT[rs] / sqrt(sigmaT[rs] * sigmaG[rs] * ss);
        double pvalue = ndist(fabs(statistic)) * 2;
        fprintf(fp, "\n%s\t%d\t%.3lf\t%.3lf\t%.3lf\t%.2G\t%.5lf",
            (const char*)rsName[rs], ss, GT[rs] / sigmaT[rs] / ss, GT[rs] / sigmaG[rs] / ss,
            statistic, pvalue, statistic * statistic / ss);
        WriteCoverageFields(fp, rs);
        if (annotCount)
            fprintf(fp, "\t%s", (const char*)annotation[rs]);
    }   // End of r loop
  fprintf(fp, "\n");
  fclose(fp);
}

void WriteScoreWS(String &outfile, StringArray &topID, const char *traitName,
  Vector &toptrait, StringArray &toprsNames, Matrix &toprisks, Vector &wss)
{
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  int idCount = toprisks.cols;
  int rsCount = toprisks.rows;
  fprintf(fp, "ID\t%s", traitName);
  for (int r = 0; r < rsCount; r++)
        fprintf(fp, "\t%s", (const char*)toprsNames[r]);
  fprintf(fp, "\tWTSUM");
  for (int i = 0; i < idCount; i++) {
        fprintf(fp, "\n%s\t%lf", (const char*)topID[i], toptrait[i]);
        for (int r = 0; r < rsCount; r++)
            fprintf(fp, "\t%.3lf", toprisks[r][i]);
        fprintf(fp, "\t%.4lf", wss[i]);
  }   // End of i loop
  fprintf(fp, "\n");
  fclose(fp);
}

void WriteTable(String &outfile, Matrix &table, StringArray &header, StringArray &col0)
{
  int colCount = table.cols;
  int rowCount = table.rows;
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  fprintf(fp, "%s", (const char*)header[0]);
  for (int c = 0; c < colCount; c++)
    fprintf(fp, "\t%s", (const char*)header[c+1]);
  for (int r = 0; r < rowCount; r++) {
    fprintf(fp, "\n%s", (const char*)col0[r]);
    for (int c = 0; c < colCount; c++)
      if(table[r][c]==_NAN_)
        fprintf(fp, "\tNA");
      else
        fprintf(fp, "\t%.5lf", table[r][c]);
  } // End of r loop
  fprintf(fp, "\n");
  fclose(fp);
}

void WriteRotatedBinaryTable(String &outfile, unsigned long long int **table, StringArray &header, StringArray &col0)
{
  int colCount = header.Length()-1;
  int rowCount = col0.Length();
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  fprintf(fp, "%s", (const char*)header[0]);
  for (int c = 0; c < colCount; c++)
    fprintf(fp, "\t%s", (const char*)header[c+1]);
  for(int r = 0; r < rowCount; r++){
    fprintf(fp, "\n%s", (const char*)col0[r]);
    for(int c = 0; c < colCount; c++)
      fprintf(fp, "\t%d", int((table[c][r>>6]>>(r&0x3F))&0x1));
  } // End of r loop
  fprintf(fp, "\n");
  fclose(fp);
}

void WriteInfoWS(String &outfile, StringArray &toprsNames, StringArray &topannotation, Vector &weights)
{
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  int annotationCount = topannotation.Length();
  int toprsCount = weights.Length();
  fprintf(fp, "PGS\tWTSUM");
  if (annotationCount) fprintf(fp, "\tAnnotation");
  for (int r = 0; r < toprsCount; r++) {
        fprintf(fp, "\n%s\t%.4lf", (const char*)toprsNames[r], weights[r]);
        if (annotationCount) fprintf(fp, "\t%s", (const char*)topannotation[r]);
  }   // End of r loop
  fprintf(fp, "\n");
  fclose(fp);
}
/*
void WriteModels(String &outfile, StringArray &toprsNames, StringArray &topannotation, Matrix &weights)
{
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  int annotationCount = topannotation.Length();
  int toprsCount = weights.Length();
  fprintf(fp, "PGS\tWTSUM");
  if (annotationCount) fprintf(fp, "\tAnnotation");
  for (int r = 0; r < toprsCount; r++) {
        fprintf(fp, "\n%s\t%.4lf", (const char*)toprsNames[r], weights[r]);
        if (annotationCount) fprintf(fp, "\t%s", (const char*)topannotation[r]);
  }   // End of r loop
  fprintf(fp, "\n");
  fclose(fp);
}
*/

double ComputeRSQ(Vector &X, Vector &Y)
{
  int N = X.Length();
  double lxy = 0.0, meanX = 0.0, meanY = 0.0, varX = 0.0, varY = 0.0;
  for (int i = 0; i < N; i++) {
    double x = X[i];
    double y = Y[i];
    lxy += x * y;
    varX += x * x;
    varY += y * y;
    meanX += x;
    meanY += y;
  } // End of i loop
  meanX /= N;
  meanY /= N;
  varY = varY / N - meanY * meanY;
  varX = varX / N - meanX * meanX;
  double rsq = lxy / N - meanX * meanY;
  rsq = rsq * rsq / (varX * varY);
  return rsq;
}

double ComputeRSQ(Vector &X, Vector &Y, IntArray &mask)
{
  int N = X.Length();
  int ss = mask.Sum();
  double lxy = 0.0, meanX = 0.0, meanY = 0.0, varX = 0.0, varY = 0.0;
  for (int i = 0; i < N; i++) {
    double x = X[i];
    double y = Y[i];
    int flag = mask[i];
    lxy += x * y;
    varX += x * x * flag;
    varY += y * y;
    meanX += x * flag;
    meanY += y;
  } // End of i loop
  meanX /= ss;
  meanY /= ss;
  varY = varY / ss - meanY * meanY;
  varX = varX / ss - meanX * meanX;
  double rsq = lxy / ss - meanX * meanY;
  rsq = rsq * rsq / (varX * varY);
  return rsq;
}

double MakeWeightedSumRS(Matrix &toprisks, Vector &toptrait, Vector &betaWS, Vector &rs_WS)
{
    int rsCount = toprisks.rows;
    if (rsCount == 0) return _NAN_;
    int idCount = toprisks.cols;
    for (int i = 0; i < idCount; i++) {
        double sum = 0.0;
        for (int r = 0; r < rsCount; r++)
            sum += toprisks[r][i] * betaWS[r];
        rs_WS[i] = sum;
    }   // End of i loop
    return ComputeRSQ(rs_WS, toptrait);
}

void NowriteScan(Matrix &risks, Matrix &allY, unsigned long long int **masks, unsigned long long int **LD,
  Matrix &RsqBest, IntArray &indexBest, int toprsCount, Vector &RsqWTSUM, int coreCount)
{
  int idCount = allY.rows;
  int traitCount = allY.cols;
  int rsCount = risks.rows;
  Matrix GTs(traitCount, rsCount);
  Matrix sigmaGs(traitCount, rsCount);
  Matrix sigmaTs(traitCount, rsCount);
  IntArray SSs(traitCount);
  MatrixLinearScan(risks, allY, masks, GTs, sigmaGs, sigmaTs, SSs, coreCount);
//  int toprsCount0 = (toprsCount << 4);
//  if(toprsCount0 > rsCount) toprsCount0 = rsCount;
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
 {
#endif
  IntArray sortedrsIndex0(rsCount), sortedrsIndex(rsCount), mask(idCount);
  Vector Z(rsCount), rs_WS(idCount), Y(idCount);
  int modelrsCount;
#ifdef _OPENMP
#pragma omp for
#endif
  for (int t = 0; t < traitCount; t++) {
    Z.Set(_NAN_);
    for (int r = 0; r < rsCount; r++) {
      if (sigmaGs[t][r] < 1E-100 || sigmaTs[t][r] < 1E-100) continue;
      Z[r] = GTs[t][r] / sqrt(sigmaGs[t][r] * sigmaTs[t][r] * SSs[t]);
    }   // End of r loop
//    sortedrsIndex0.Dimension(rsCount);
    GenomicControlNumber(Z, sortedrsIndex0);
    indexBest[t] = sortedrsIndex0[0];
//    RsqBest[t] = Z[sortedrsIndex[0]] * Z[sortedrsIndex[0]] / SSs[t];
    sortedrsIndex = sortedrsIndex0;
    sortedrsIndex.Dimension(toprsCount);
    MakePrunedOrder(LD, sortedrsIndex);
    int modelrsCount = sortedrsIndex.Length();
    rs_WS.Zero();
    for(int r = 0; r < modelrsCount; r++){
      int rr = sortedrsIndex0[sortedrsIndex[r]];
      double beta = GTs[t][rr] / sigmaGs[t][rr] / SSs[t];
      for(int i = 0; i < idCount; i++)
        rs_WS[i] += risks[rr][i] * beta;
    } // End of r loop
    for(int i = 0; i < idCount; i++){
      mask[i] = ((masks[i][t>>6]>>(t&0x3F))&0x1);
      Y[i] = allY[i][t];
    } // End of i loop
    RsqWTSUM[t] = ComputeRSQ(rs_WS, Y, mask);
  } // End of OMP t loop
#ifdef _OPENMP
 } // extra bracket for omp
#endif
 for(int i = 0; i < traitCount; i++){
   int r = indexBest[i];
   for(int t = 0; t < traitCount; t++){
      double z = GTs[t][r] / sqrt(sigmaGs[t][r] * sigmaTs[t][r] * SSs[t]);
      RsqBest[i][t] = z * z / SSs[t];
   } // End of t loop
 } // End of i loop
}

void WithwriteScan(Matrix &risks, Matrix &allY, unsigned long long int **masks, unsigned long long int **LD,
  Matrix &RsqBest, IntArray &indexBest, int toprsCount, Vector &RsqWTSUM,
  IntArray *sortedrsIndex0, IntArray *sortedrsIndex, Matrix &GTs, Matrix &sigmaGs, Matrix &sigmaTs, int coreCount)
{
  int idCount = allY.rows;
  int traitCount = allY.cols;
  int rsCount = risks.rows;
  IntArray SSs(traitCount);
//************* Genome scan starts ************//
  MatrixLinearScan(risks, allY, masks, GTs, sigmaGs, sigmaTs, SSs, coreCount);
//************* Genome scan ends **************//
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
 {
#endif
  IntArray mask(idCount);
  Vector Z(rsCount), betaWS(toprsCount), rs_WS(idCount), Y(idCount);
#ifdef _OPENMP
#pragma omp for
#endif
  for (int t = 0; t < traitCount; t++) {
    Z.Set(_NAN_);
    for (int r = 0; r < rsCount; r++) {
      if (sigmaGs[t][r] < 1E-100 || sigmaTs[t][r] < 1E-100) continue;
      Z[r] = GTs[t][r] / sqrt(sigmaGs[t][r] * sigmaTs[t][r] * SSs[t]);
    }   // End of r loop
    GenomicControlNumber(Z, sortedrsIndex0[t]);
    indexBest[t] = sortedrsIndex0[t][0];
    sortedrsIndex[t] = sortedrsIndex0[t];
    sortedrsIndex[t].Dimension(toprsCount);
    MakePrunedOrder(LD, sortedrsIndex[t]);
    int modelrsCount = sortedrsIndex[t].Length();
    betaWS.Dimension(modelrsCount);
    for (int r = 0; r < modelrsCount; r++){
      int rr = sortedrsIndex0[t][sortedrsIndex[t][r]];
      betaWS[r] = GTs[t][rr] / sigmaGs[t][rr] / SSs[t];
    } // End of r loop
    rs_WS.Zero();
    for(int r = 0; r < modelrsCount; r++){
      int rr = sortedrsIndex0[t][sortedrsIndex[t][r]];
      double beta = betaWS[r];
      for(int i = 0; i < idCount; i++)
        rs_WS[i] += risks[rr][i] * beta;
    } // End of r loop
    for(int i = 0; i < idCount; i++){
      mask[i] = ((masks[i][t>>6]>>(t&0x3F))&0x1);
      Y[i] = allY[i][t];
    } // End of i loop
    RsqWTSUM[t] = ComputeRSQ(rs_WS, Y, mask);
  } // End of OMP t loop
#ifdef _OPENMP
 } // extra bracket for omp
#endif
  for(int i = 0; i < traitCount; i++){
    int r = indexBest[i];
    for(int t = 0; t < traitCount; t++){
      double z = GTs[t][r] / sqrt(sigmaGs[t][r] * sigmaTs[t][r] * SSs[t]);
      RsqBest[i][t] = z * z / SSs[t];
    } // End of t loop
  } // End of i loop
}

void ScanCV1Time(Matrix &risks, Matrix &allY, unsigned long long int **masks,
  unsigned long long int **LD, Vector &RsqBest, int toprsCount, Vector &RsqWS,
  unsigned long long int *testFlag, Matrix &GY0, Matrix &GG0, Matrix &G0,
  Vector &YY0, Vector &Y0, IntArray &SS0, int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
  int traitCount = RsqWS.Length();
  int twordCount = (traitCount-1)/64+1;
  IntArray tIndex(0);
  for(int i = 0; i < idCount; i+=64){
    int iMax = i < idCount - 64? i+64: idCount;
    unsigned long long int word = testFlag[i>>6];
    for(int ii = i; ii < iMax; ii++, word>>=1)
      if(word & 0x1) // set 0 for test set
        tIndex.Push(ii);
  } // End of i loop
  int testSize = tIndex.Length();
  unsigned long long int **tmasks = new unsigned long long int *[testSize];
  for(int i = 0; i < testSize; i++)
    tmasks[i] = new unsigned long long int[twordCount];
  Matrix trisks(rsCount, testSize), tallY(testSize, traitCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int i = 0; i < testSize; i++){
    int ii = tIndex[i];
    memcpy(&tallY[i][0], &allY[ii][0], traitCount*sizeof(double));
    memcpy(&tmasks[i][0], &masks[ii][0], twordCount*sizeof(unsigned long long int));
    for(int r = 0; r < rsCount; r++)
      trisks[r][i] = risks[r][ii];
  } // End of OMP i loop
  Matrix GY(traitCount, rsCount), GG(traitCount, rsCount), SG(traitCount, rsCount);
  Vector YY(traitCount), SY(traitCount);
  IntArray SS(traitCount);
  MatrixLinearScan_Sum(trisks, tallY, tmasks, GY, GG, SG, YY, SY, SS, coreCount);
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
{
#endif
  Vector Z(rsCount), Y(testSize), rs_WS(testSize), betaWS, sigmaGs(rsCount), GTs(rsCount);
  IntArray sortedrsIndex0(rsCount), sortedrsIndex(toprsCount), mask(testSize);
#ifdef _OPENMP
#pragma omp for
#endif
  for(int t = 0; t < traitCount; t++) {
    int SSs = SS0[t] - SS[t];
    double yy = YY0[t] - YY[t];
    double y = Y0[t] - SY[t];
    double sigmaTs = (yy - y*y/SSs)/SSs;
    for(int r = 0; r < rsCount; r++){
      GTs[r] = GY0[t][r] - GY[t][r] - (G0[t][r]-SG[t][r]) * y / SSs;
      sigmaGs[r] = (GG0[t][r] - GG[t][r] - (G0[t][r]-SG[t][r])*(G0[t][r]-SG[t][r])/SSs) / SSs;
    } // End of r loop
    if (sigmaTs < 1E-100) continue;
    Z.Set(_NAN_);
    for (int r = 0; r < rsCount; r++)
      if (sigmaGs[r] > 1E-100)
        Z[r] = GTs[r] / sqrt(sigmaGs[r] * sigmaTs * SSs);
    GenomicControlNumber(Z, sortedrsIndex0);
    int tW = (t>>6);
    int tS = (t&0x3F);
    for(int i = 0; i < testSize; i++){
      mask[i] = (tmasks[i][tW]>>tS)&0x1;
      Y[i] = tallY[i][t];
    } // End of i loop
    RsqBest[t] = ComputeRSQ(trisks[sortedrsIndex0[0]], Y, mask);
    sortedrsIndex = sortedrsIndex0;
    sortedrsIndex.Dimension(toprsCount);
    MakePrunedOrder(LD, sortedrsIndex);
    int modelrsCount = sortedrsIndex.Length();
    betaWS.Dimension(modelrsCount);
    for (int r = 0; r < modelrsCount; r++){
      int rr = sortedrsIndex0[sortedrsIndex[r]];
      betaWS[r] = GTs[rr] / sigmaGs[rr] / SSs;
    } // End of r loop
    rs_WS.Zero();
    for(int r = 0; r < modelrsCount; r++){
      double beta = betaWS[r];
      int rr = sortedrsIndex0[sortedrsIndex[r]];
      for(int i = 0; i < testSize; i++)
        rs_WS[i] += trisks[rr][i] * beta;
    } // End of r loop
    RsqWS[t] = ComputeRSQ(rs_WS, Y, mask);
  } // End of t loop
#ifdef _OPENMP
 } // extra bracket for omp
#endif
  for(int i = 0; i < testSize; i++)
    delete []tmasks[i];
  delete []tmasks;
}

void ScanCV(Matrix &risks, Matrix &allY, unsigned long long int **masks,
  unsigned long long int **LD, int fold, Matrix &RsqBest, int toprsCount,
  Matrix &RsqWS, unsigned long long int **testFlags, int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
  int traitCount = allY.cols;
  int repeat = RsqBest.rows;
  Random ran;
  IntArray testFlag(idCount);
  int imemsize = ((idCount-1)/64+1)*sizeof(unsigned long long int);
  int testsetCount = idCount * 1.0 / fold;
  for(int run = 0; run < repeat; run++){
    ran.Choose(testFlag, idCount, testsetCount);
    memset(&testFlags[run][0], 0, imemsize);
    for(int i = 0; i < idCount; i+=64){
      int iMax = i<idCount-64? i+64: idCount;
      int w = (i>>6);
      unsigned long long int word0 = 0x1;
      for(int ii = i; ii < iMax; ii++, word0 <<= 1)
        if(testFlag[ii]) testFlags[run][w] |= word0;
    } // End of i loop
  } // End of run loop
  Matrix GY0(traitCount, rsCount), GG0(traitCount, rsCount), G0(traitCount, rsCount);
  Vector YY0(traitCount), Y0(traitCount);
  IntArray SS0(traitCount);
  MatrixLinearScan_Sum(risks, allY, masks, GY0, GG0, G0, YY0, Y0, SS0, coreCount);
  for(int run = 0; run < repeat; run++){
    printf(".");fflush(stdout);
    ScanCV1Time(risks, allY, masks, LD, RsqBest[run], toprsCount, RsqWS[run],
      testFlags[run], GY0, GG0, G0, YY0, Y0, SS0, coreCount);
  } // End of run loop
}

void ScanCV1Time(Matrix &risks, Matrix &allY, unsigned long long int **masks,
  unsigned long long int ***LD, int r2Count, IntArray &toprsCount,
  Matrix *RsqWS, unsigned long long int *testFlag, Matrix &GY0, Matrix &GG0,
  Matrix &G0, Vector &YY0, Vector &Y0, IntArray &SS0, int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
  int traitCount = allY.cols;
  int parCount = toprsCount.Length();
  int twordCount = (traitCount-1)/64+1;
  IntArray tIndex(0);
  for(int i = 0; i < idCount; i+=64){
    int iMax = i < idCount - 64? i+64: idCount;
    unsigned long long int word = testFlag[i>>6];
    for(int ii = i; ii < iMax; ii++, word>>=1)
      if(word & 0x1) // set 0 for test set
        tIndex.Push(ii);
  } // End of i loop
  int testSize = tIndex.Length();
  unsigned long long int **tmasks = new unsigned long long int *[testSize];
  for(int i = 0; i < testSize; i++)
    tmasks[i] = new unsigned long long int[twordCount];
  Matrix trisks(rsCount, testSize), tallY(testSize, traitCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int i = 0; i < testSize; i++){
    int ii = tIndex[i];
    memcpy(&tallY[i][0], &allY[ii][0], traitCount*sizeof(double));
    memcpy(&tmasks[i][0], &masks[ii][0], twordCount*sizeof(unsigned long long int));
    for(int r = 0; r < rsCount; r++)
      trisks[r][i] = risks[r][ii];
  } // End of OMP i loop
  Matrix GY(traitCount, rsCount), GG(traitCount, rsCount), SG(traitCount, rsCount);
  Vector YY(traitCount), SY(traitCount);
  IntArray SS(traitCount);
  MatrixLinearScan_Sum(trisks, tallY, tmasks, GY, GG, SG, YY, SY, SS, coreCount);
  IntArray toprsCount0(parCount);
  for(int n = 0; n < parCount; n++)
    toprsCount0[n] = (toprsCount[n] << 4) > rsCount? rsCount: (toprsCount[n] << 4);
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount) private(toprsCount)
{
#endif
  Vector Z(rsCount), Y(testSize), rs_WS(testSize), betaWS, GTs(rsCount), sigmaGs(rsCount);
  IntArray sortedrsIndex(rsCount), mask(testSize);
#ifdef _OPENMP
#pragma omp for
#endif
  for(int t = 0; t < traitCount; t++) {
    int SSs = SS0[t] - SS[t];
    double yy = YY0[t] - YY[t];
    double y = Y0[t] - SY[t];
    double sigmaTs = (yy - y*y/SSs)/SSs;
    for(int r = 0; r < rsCount; r++){
      GTs[r] = GY0[t][r] - GY[t][r] - (G0[t][r]-SG[t][r]) * y / SSs;
      sigmaGs[r] = (GG0[t][r] - GG[t][r] - (G0[t][r]-SG[t][r])*(G0[t][r]-SG[t][r])/SSs) / SSs;
    } // End of r loop
    if (sigmaTs < 1E-100) continue;
    Z.Set(_NAN_);
    for (int r = 0; r < rsCount; r++)
      if (sigmaGs[r] > 1E-100)
        Z[r] = GTs[r] / sqrt(sigmaGs[r] * sigmaTs * SSs);
    GenomicControlNumber(Z, sortedrsIndex);
    int tW = (t>>6);
    int tS = (t&0x3F);
    for(int i = 0; i < testSize; i++){
      mask[i] = (tmasks[i][tW]>>tS)&0x1;
      Y[i] = tallY[i][t];
    } // End of i loop
    IntArray sortedrsIndex0;
    for(int n = 0; n < parCount; n++){
      sortedrsIndex0.Dimension(toprsCount0[n]);
      for(int r = 0; r < toprsCount0[n]; r++)
        sortedrsIndex0[r] = sortedrsIndex[r];
      for(int c = 0; c < r2Count; c++){
        MakePrunedOrder(LD[c], sortedrsIndex0);
        int modelrsCount = sortedrsIndex0.Length();
        betaWS.Dimension(modelrsCount);
        for (int r = 0; r < modelrsCount; r++){
          int rr = sortedrsIndex[sortedrsIndex0[r]];
          betaWS[r] = GTs[rr] / sigmaGs[rr] / SSs;
        } // End of r loop
        rs_WS.Zero();
        for(int r = 0; r < modelrsCount; r++){
          int rr = sortedrsIndex[sortedrsIndex0[r]];
          double beta = betaWS[r];
          for(int i = 0; i < testSize; i++)
          rs_WS[i] += trisks[rr][i] * beta;
        } // End of r loop
        RsqWS[t][n][c] = ComputeRSQ(rs_WS, Y, mask);
      } // End of c loop
    } // End of n loop
  } // End of t loop
#ifdef _OPENMP
 } // extra bracket for omp
#endif
  for(int i = 0; i < testSize; i++)
    delete []tmasks[i];
  delete []tmasks;
}

void AutoCV(Matrix &risks, Matrix &allY, unsigned long long int **masks,
    int ntop, double ldT, int fold, int repeat, IntArray &maxPar, Vector &maxR2,
    Vector &maxRsq, Vector &maxVar, int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
  int traitCount = allY.cols;
  int MAXR2COUNT;
  Vector r2;
  if(ldT==_NAN_){
    MAXR2COUNT = 6;
    r2.Dimension(MAXR2COUNT);
    for(int k = 0; k < MAXR2COUNT; k++)
      r2[k] = (k+1)*0.1;
  }else{
    MAXR2COUNT = 1;
    r2.Dimension(MAXR2COUNT);
    r2[0] = ldT;
  }
  time_t timer;
  struct tm *tblock;
  timer = time(NULL);
  tblock = localtime(&timer);
  printf("Computing PGS correlations starts at %s", asctime(tblock));
  unsigned long long int ***LD = new unsigned long long int **[MAXR2COUNT];
  int rwordCount = (rsCount-1) / 64 + 1;
  IntArray *pairList = new IntArray[coreCount];
  Vector *pairLD = new Vector[coreCount];
  ScreenLD1(risks, r2[0], pairList, coreCount);
  Vector rsSD(rsCount); rsSD.Zero();
  ScreenLD2(risks, r2[0], pairList, pairLD, rsSD, coreCount);
  for(int c = 0; c < MAXR2COUNT; c++){
    printf("R2 = %.3lf\n", r2[c]);
    LD[c] = new unsigned long long int *[rsCount];
    for(int r = 0; r < rsCount; r++){
      LD[c][r] = new unsigned long long int[rwordCount];
      memset(&LD[c][r][0], 0, rwordCount*sizeof(unsigned long long int));
    } // End of r loop
    if(c>0)
      for(int core = 0; core < coreCount; core++){
        int newpairCount = 0;
        int pairCount = pairList[core].Length();
        for(int p = 0; p < pairCount; p++)
          if(pairLD[core][p]*pairLD[core][p] >= r2[c]){
            pairList[core][newpairCount] = pairList[core][p];
            pairLD[core][newpairCount++] = pairLD[core][p];
          } // End of p loop
        pairList[core].Dimension(newpairCount);
        pairLD[core].Dimension(newpairCount);
      } // End of core loop
    MakeBinaryLD(rsCount, pairList, LD[c], coreCount);
    LDSummary(LD[c], rsCount, rwordCount, coreCount);
  } // End of c loop
  delete []pairList;
  delete []pairLD;
  timer = time(NULL);
  tblock = localtime(&timer);
  printf("Computing PGS correlations ends at %s", asctime(tblock));
  Random ran;
  int iwordCount = (idCount-1)/64+1;
  unsigned long long int **testFlags = new unsigned long long int *[repeat];
  for(int run = 0; run < repeat; run++)
    testFlags[run] = new unsigned long long int[iwordCount];
  IntArray testFlag(idCount);
  int imemsize = ((idCount-1)/64+1)*sizeof(unsigned long long int);
  int testsetCount = idCount * 1.0 / fold;
  for(int run = 0; run < repeat; run++){
    ran.Choose(testFlag, idCount, testsetCount);
    memset(&testFlags[run][0], 0, imemsize);
    for(int i = 0; i < idCount; i+=64){
      int iMax = i<idCount-64? i+64: idCount;
      int w = (i>>6);
      unsigned long long int word0 = 0x1;
      for(int ii = i; ii < iMax; ii++, word0 <<= 1)
        if(testFlag[ii]) testFlags[run][w] |= word0;
    } // End of i loop
  } // End of run loop
  int MAXTOPCOUNT;
  IntArray toprsCounts;
  if(ntop==0){
    MAXTOPCOUNT=6;
    toprsCounts.Dimension(MAXTOPCOUNT);
    toprsCounts[0] = 1;
    toprsCounts[1] = 3;
    for(int k = 2; k < MAXTOPCOUNT; k++)
      toprsCounts[k] = k*30 - 50;
  }else{
    MAXTOPCOUNT=1;
    toprsCounts.Dimension(MAXTOPCOUNT);
    toprsCounts[0] = ntop;
  }
  Matrix *means = new Matrix[traitCount];
  Matrix *vars = new Matrix[traitCount];
  Matrix *RsqWS = new Matrix[traitCount];
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int t = 0; t < traitCount; t++){
    means[t].Dimension(MAXTOPCOUNT, MAXR2COUNT); means[t].Zero();
    vars[t].Dimension(MAXTOPCOUNT, MAXR2COUNT); vars[t].Zero();
    RsqWS[t].Dimension(MAXTOPCOUNT, MAXR2COUNT);
  } // End of t loop
  printf("\nCross-validation starts"); fflush(stdout);
  Matrix GY0(traitCount, rsCount), GG0(traitCount, rsCount), G0(traitCount, rsCount);
  Vector YY0(traitCount), Y0(traitCount);
  IntArray SS0(traitCount);
  MatrixLinearScan_Sum(risks, allY, masks, GY0, GG0, G0, YY0, Y0, SS0, coreCount);
  for(int run = 0; run < repeat; run++){
    printf(".");fflush(stdout);
    ScanCV1Time(risks, allY, masks, LD, MAXR2COUNT, toprsCounts, RsqWS,
      testFlags[run], GY0, GG0, G0, YY0, Y0, SS0, coreCount);
    for(int t = 0; t < traitCount; t++)
      for(int k = 0; k < MAXTOPCOUNT; k++)
        for(int c = 0; c < MAXR2COUNT; c++){
          means[t][k][c] += RsqWS[t][k][c];
          vars[t][k][c] += RsqWS[t][k][c] * RsqWS[t][k][c];
        } // End of c loop
  } // End of run loop
  delete []RsqWS;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int t = 0; t < traitCount; t++){
    double max = -0.1;
    int indexN = -1;
    int indexC = -1;
    for(int k = 0; k < MAXTOPCOUNT; k++)
      for(int c = 0; c < MAXR2COUNT; c++)
        if(means[t][k][c] > max){
          max = means[t][k][c];
          indexN = k;
          indexC = c;
        }
    maxRsq[t] = max / repeat;
    maxPar[t] = toprsCounts[indexN];
    maxR2[t] = r2[indexC];
    maxVar[t] = vars[t][indexN][indexC]/repeat - maxRsq[t] * maxRsq[t];
  } // End of OMP t loop
  delete []means;
  delete []vars;
  for(int c = 0; c < MAXR2COUNT; c++){
    for(int r = 0; r < rsCount; r++)
      delete []LD[c][r];
    delete []LD[c];
  } // End of c loop
  delete []LD;
}

void Linear_ReadModel(StringArray &pgsName, Matrix &coef, StringArray &rsNames, Matrix &risks,
  Matrix &allY, unsigned long long int **masks, Matrix &rsWS, Matrix &RsqWS)
{
  int idCount = allY.rows;
  int traitCount = allY.cols;
  int modelrsCount = coef.rows;
  int methodCount = coef.cols;
  rsWS.Dimension(methodCount, idCount); rsWS.Zero();
  for(int m = 0; m < methodCount; m++)
    for(int r = 0; r < modelrsCount; r++){
      int rr = rsNames.Find(pgsName[r]);
      if(rr==-1) {printf("PGS %s does not exist\n", (const char*)pgsName[r]); continue;}
      double beta = coef[r][m];
      for(int i = 0; i < idCount; i++)
        rsWS[m][i] += risks[rr][i] * beta;
    } // End of r loop
  IntArray mask(idCount);
  Vector Y(idCount);
  RsqWS.Dimension(traitCount, methodCount);
  for(int t = 0; t < traitCount; t++){
    for(int i = 0; i < idCount; i++){
      mask[i] = ((masks[i][t>>6]>>(t&0x3F))&0x1);
      Y[i] = allY[i][t];
    } // End of i loop
    for(int m = 0; m < methodCount; m++)
      RsqWS[t][m] = ComputeRSQ(rsWS[m], Y, mask);
  } // End of r loop
}

void Linear_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *modelfile,
  const char *prefix, const char *covariateList, bool rplotFlag,
  bool cvFlag, int nTop, double ldT, int fold, int repeat, bool nowriteFlag, int coreCount)
{
  printf("\nOptions in effect:\n");
  printf("\t--linear\n");
  if(cvFlag) printf("\t--cv\n");
  printf("\t--pheno %s\n", phefile);
  if (gCovarFile[0]) printf("\t--covarfile %s\n", gCovarFile);
  if (gTraitList[0]) printf("\t--trait %s\n", gTraitList);
  if(covariateList[0]) printf("\t--covariate %s\n", covariateList);
  printf("\t--score %s\n", rsfile);
  if(dictfile[0]) printf("\t--diction %s\n", dictfile);
  if(modelfile[0]) printf("\t--readmodel %s\n", modelfile);
  if(nTop) printf("\t--ntop %d\n", nTop);
  if(fold) printf("\t--fold %d\n", fold);
  if(repeat) printf("\t--repeat %d\n", repeat);
  if (ldT != _NAN_ && ldT < 1 && ldT >= 0)
    printf("\t--pruner2 %.3lf\n", ldT);
  else ldT = _NAN_;
  if (rplotFlag) printf("\t--rplot\n");
  if(nowriteFlag) printf("\t--nowrite\n");
  if (coreCount) printf("\t--cpus %d\n", coreCount);
  if (prefix) printf("\t--prefix %s\n", prefix);
  printf("\n");

  StringArray sampleID, alltraitNames, rsNames, annotation;
  IntArray covariates, affections, quantraits;
  Matrix alltraits, risks;
  SECRET_Read(phefile, dictfile, rsfile, modelfile, prefix, covariateList,
    sampleID, alltraitNames, alltraits, covariates, affections, quantraits,
    rsNames, risks, annotation, coreCount);
  int rsCount = rsNames.Length();
  int idCount = sampleID.Length();
  if (idCount == 0) return;
  if(ldT==_NAN_){
    ldT = 0.2;
    printf("R2 to prune correlated PGSs for WTSUM is missing. --pruner2=%.1lf is set\n", ldT);
  }
  if (!nTop){
    nTop = rsCount<30? rsCount: 30;
    printf("The number of PGSs for a combined PGS (ntop) is missing. --ntop=%d is set\n", nTop);
  }
  int annotationCount = annotation.Length();
  IntArray traits = quantraits;
  traits.Stack(affections);
  int traitCount = traits.Length();
  if(traitCount == 0) { printf("Traits not found.\n"); return; }
  Matrix allY(idCount, traitCount);
  int twordCount = (traitCount-1)/64+1;
  unsigned long long int **masks = new unsigned long long int *[idCount];
  for(int i = 0; i < idCount; i++)
    masks[i] = new unsigned long long int [twordCount];
  Matrix allX(idCount, covariates.Length());
  MakePhenos(risks, alltraits, alltraitNames, affections, quantraits, covariates, allY, allX, masks, coreCount);
  String msg(256);
  if(modelfile[0]){
    Matrix coef;
    StringArray methodName, pgsName;
    ReadMatrix(modelfile, methodName, pgsName, coef, coreCount);
    Matrix rsWS, RsqWS;
    Linear_ReadModel(pgsName, coef, rsNames, risks, allY, masks, rsWS, RsqWS);
    printf("\nSummary of prediction Rsq for %d Traits\n%20s", traitCount, "TRAIT");
    int methodCount = coef.cols;
    for(int m = 0; m < methodCount; m++)
      printf(" %10s", (const char*)methodName[m]);
    printf("\n");
    for (int t = 0; t < traitCount; t++){
      printf("%20s", (const char*)alltraitNames[traits[t]]);
      for(int m = 0; m < methodCount; m++)
        printf(" %10.3lf", RsqWS[t][m]);
      printf("\n");
    } // End of t loop
    if(rplotFlag){
      IntArray mask(idCount);
      Vector Y(idCount);
      for(int t = 0; t < traitCount; t++){
        for(int i = 0; i < idCount; i++){
          mask[i] = ((masks[i][t>>6]>>(t&0x3F))&0x1);
          Y[i] = allY[i][t];
        } // End of i loop
        String outfile0 = prefix;
        outfile0.Add("_");
        outfile0.Add(alltraitNames[traits[t]]);
        Matrix pgsTable(idCount, methodCount+1); pgsTable.Zero();
        for(int i = 0; i < idCount; i++){
          pgsTable[i][0] = mask[i]? Y[i]: _NAN_;
          for(int m = 0; m < methodCount; m++)
            pgsTable[i][m+1] = rsWS[m][i];
        } // End of i loop
        StringArray header(methodCount+2);
        header[0] = "ID";
        header[1] = alltraitNames[traits[t]];
        for(int m = 0; m < methodCount; m++)
          header[m+2] = methodName[m];
        StringArray col0(idCount);
        for(int i = 0; i < idCount; i++)
          col0[i] = sampleID[i];
        String outfile = outfile0; outfile.Add("_combinedPGS.txt");
        WriteTable(outfile, pgsTable, header, col0);
        printf("  %d polygenic scores for trait %s write to %s\n",
          methodCount, (const char*)alltraitNames[traits[t]], (const char*)outfile);
        plotCalibration(outfile0, "", msg);
        printf("  %s", (const char*)msg);
      } // End of t loop
    }
    return;
  }
  int wordCount = (rsCount-1) / 64 + 1;
  unsigned long long int **LD = new unsigned long long int *[rsCount];
  for(int r = 0; r < rsCount; r++){
    LD[r] = new unsigned long long int[wordCount];
    memset(&LD[r][0], 0, wordCount*sizeof(unsigned long long int));
  } // End of r loop
  Vector rsSD(rsCount);
  if(ldT != _NAN_){
    time_t timer;
    struct tm *tblock;
    timer = time(NULL);
    tblock = localtime(&timer);
    printf("Computing PGS correlations starts at %s", asctime(tblock)); fflush(stdout);
    MakeBinaryLD(risks, ldT, LD, rsSD, coreCount);
    LDSummary(LD, rsCount, wordCount, coreCount);
    timer = time(NULL);
    tblock = localtime(&timer);
    printf("Computing PGS correlations ends at %s", asctime(tblock)); fflush(stdout);
  }
  int toprsCount = nTop;
  IntArray SS(traitCount); // Sample size
  popcountbycol(masks, idCount, SS);
  Matrix RsqBest(traitCount, traitCount);
  Vector RsqWS(traitCount);
  IntArray indexBest(traitCount);
  IntArray *sortedrsIndex0 = new IntArray[traitCount];
  IntArray *sortedrsIndex = new IntArray[traitCount];
  for(int t = 0; t < traitCount; t++){
    sortedrsIndex0[t].Dimension(rsCount);
    sortedrsIndex[t].Dimension(rsCount);
  } // End of t loop
  Matrix GTs(traitCount, rsCount);
  Matrix sigmaGs(traitCount, rsCount);
  Matrix sigmaTs(traitCount, rsCount);
  WithwriteScan(risks, allY, masks, LD, RsqBest, indexBest, toprsCount, RsqWS,
    sortedrsIndex0, sortedrsIndex, GTs, sigmaGs, sigmaTs, coreCount);

  StringArray header(traitCount+1), col0(traitCount);
  if(rplotFlag && traitCount>1) { // multiple trait only
    header[0] = "PGS_Best";
    for(int t = 0; t < traitCount; t++){
      header[t+1] = alltraitNames[traits[t]];
      col0[t] = (const char*)rsNames[indexBest[t]];
    } // End of t loop
    String outfile = prefix; outfile.Add("_BestPGS_Rsq.txt");
    WriteTable(outfile, RsqBest, header, col0);
    printf("Rsqs for Best PGSs from %d traits saved in file %s\n", traitCount, (const char*)outfile);
    outfile = prefix; outfile.Add("_BestPGS_Rsq");
    plotBestHeatmap(outfile, "", msg);
    printf("%s", (const char*)msg);
  } // End of BestPGS for multiple traits
  Matrix rsqs[2];
  for(int k = 0; k < 2; k++)
    rsqs[k].Dimension(traitCount, 3);
  IntArray pheno;
  Vector Y, meanBest(traitCount), meanWS(traitCount);
  if(cvFlag){ // Cross-validation
    if(fold==0) fold=10;
    if(repeat==0) repeat=10;
    int iwordCount = (idCount-1)/64+1;
    unsigned long long int **testFlags = new unsigned long long int *[repeat];
    for(int run = 0; run < repeat; run++)
      testFlags[run] = new unsigned long long int[iwordCount];
    Matrix RsqBest0(repeat, traitCount), RsqWS0(repeat, traitCount);
    printf("Cross-validation start");
    ScanCV(risks, allY, masks, LD, fold, RsqBest0, toprsCount, RsqWS0, testFlags, coreCount);
    printf("done.\n");
    Vector V(repeat);
    for(int t = 0; t < traitCount; t++){
      for(int w = 0; w < repeat; w++)
        V[w] = RsqBest0[w][t];
      meanBest[t] = V.Average();
      for(int w = 0; w < repeat; w++)
        V[w] = RsqWS0[w][t];
      meanWS[t] = V.Average();
    } // End of t loop
  }

    printf("\nSummary of %d Trait Prediction\n", traitCount);
    printf("%20s%10s%10s%10s  PGS_Best\n", "TRAIT", "N", "RSQ_Best", "RSQ_WTSUM");
    for (int t = 0; t < traitCount; t++) {
      printf("%20s%10d%10.3lf%10.3lf  %s", (const char*)alltraitNames[traits[t]],
        SS[t], RsqBest[t][t], RsqWS[t], (const char*)rsNames[indexBest[t]]);
      if(annotationCount) printf(":%s", (const char*)annotation[indexBest[t]]);
      printf("\n");
    } // End of t loop

  // the penalized models now always run, so the CV summary is printed
  // whenever --cv is given (it used to be suppressed by --net)
  if(cvFlag){
    printf("\nSummary of %d Trait Prediction with CV\n", traitCount);
    printf("%20s%10s%10s%10s\n", "TRAIT", "N", "RSQ_Best", "RSQ_WTSUM");
    for(int t = 0; t < traitCount; t++)
      printf("%20s%10d%10.3lf%10.3lf\n", (const char*)alltraitNames[traits[t]], SS[t], meanBest[t], meanWS[t]);
  }
  fflush(stdout);
  for(int t = 0; t < traitCount; t++){
    printf("Trait %s...\n", (const char*)alltraitNames[traits[t]]);
    int tw = (t>>6);
    int tb = (t&0x3F);
    Y.Dimension(0);
    pheno.Dimension(0);
    for(int i = 0; i < idCount; i++)
      if((masks[i][tw]>>tb) & 0x1) {
        Y.Push(allY[i][t]);
        pheno.Push(i);
      }// End of i loop
    int totalCount = Y.Length();
    String outfile0 = prefix;
    outfile0.Add("_");
    outfile0.Add(alltraitNames[traits[t]]);
    String outfile = outfile0; outfile.Add(".linear");
    WriteLinear(outfile, sortedrsIndex0[t], rsNames, annotation, GTs[t], sigmaGs[t], sigmaTs[t], SS[t]);
    printf("  PGS scan statistics write to %s\n", (const char*)outfile);
    {
      // correlations among the top PGSs (the case/control analysis is --logistic only)
      int topCount = toprsCount < sortedrsIndex0[t].Length()? toprsCount: sortedrsIndex0[t].Length();
      IntArray topIndex(topCount), corrStatus(0);
      for(int r = 0; r < topCount; r++)
        topIndex[r] = sortedrsIndex0[t][r];
      bool caseControlDone = false;
      WritePGSCorrelations(outfile0, risks, rsNames, topIndex, pheno, corrStatus, coreCount, caseControlDone, ldT);
      if(rplotFlag){
        plotPGSCorr(outfile0, outfile, alltraitNames[traits[t]], "", msg);
        printf("  %s", (const char*)msg);
      }
    }

    int modelrsCount = sortedrsIndex[t].Length();
    Vector beta(modelrsCount);
    Vector rs_WS(totalCount); rs_WS.Zero();
    for(int r = 0; r < modelrsCount; r++){
        int rr = sortedrsIndex0[t][sortedrsIndex[t][r]];
        beta[r] = GTs[t][rr] / sigmaGs[t][rr] / SS[t];
        for(int i = 0; i < totalCount; i++)
          rs_WS[i] += risks[rr][pheno[i]] * beta[r];
    } // End of r loop
    Matrix pgsTable(totalCount, toprsCount+2); pgsTable.Zero();
    for(int i = 0; i < totalCount; i++){
        pgsTable[i][0] = Y[i];
        for(int r = 0; r < toprsCount; r++)
          pgsTable[i][r+1] = risks[sortedrsIndex0[t][r]][pheno[i]];
        pgsTable[i][toprsCount+1] = rs_WS[i];
    } // End of i loop
    header.Dimension(toprsCount+3);
    header[0] = "ID";
    header[1] = alltraitNames[traits[t]];
    for(int r = 0; r < toprsCount; r++)
        header[r+2] = rsNames[sortedrsIndex0[t][r]];
    header[toprsCount+2] = "WTSUM";
    col0.Dimension(totalCount);
    for(int i = 0; i < totalCount; i++)
        col0[i] = sampleID[pheno[i]];
    outfile = outfile0; outfile.Add("_PGS.txt");
    WriteTable(outfile, pgsTable, header, col0);
    printf("  Top %d polygenic scores for trait %s write to %s\n",
        toprsCount, (const char*)alltraitNames[traits[t]], (const char*)outfile);
    // Ridge / elastic net / lasso are fitted in-process by Penalized.cpp.
    // plotRSQ() below still runs R for the plots, but it no longer writes
    // the tmp_*_net_*.txt files, so these results are not overwritten.
    PenalizedRegressionR2(outfile0, alltraitNames[traits[t]], repeat, coreCount, msg);
    printf("  %s", (const char*)msg);
    if(rplotFlag){
      // plots only; the penalized numbers were computed in C++ just above
      plotRSQ(outfile0, alltraitNames[traits[t]], "", msg);
      printf("  %s", (const char*)msg);
    }
    {
      StringArray lines;
      outfile = "tmp_"; outfile.Add(outfile0); outfile.Add("_net_rsq.txt");
      ReadLines(outfile, lines, false);
      StringArray tokens;
      for(int m = 0; m < 3; m++){
        tokens.Clear();
        tokens.AddTokens(lines[m], '\t');
        for(int k = 0; k < 2; k++)
          rsqs[k][t][m] = tokens[k].AsDouble();
      } // End of m loop
    }
    // WriteModels
    int modelCount = 5;
    Matrix weights(toprsCount, modelCount); weights.Zero();
    // prepare weights here
    header.Dimension(modelCount+1);
    header[0] = "PGS_NAME";
    header[1] = "BEST";
    header[2] = "WTSUM";
    header[3] = "RIDGE";
    header[4] = "NET";
    header[5] = "LASSO";
    for(int r = 0; r < toprsCount; r++)
      col0[r] = rsNames[sortedrsIndex0[t][r]];

    if(rsSD[sortedrsIndex0[t][sortedrsIndex[t][0]]] < 1E-100) error("rs SD not initiated.");
    weights[sortedrsIndex[t][0]][0] = 1.0;
    for(int r = 0; r < modelrsCount; r++)
      weights[sortedrsIndex[t][r]][1] = beta[r] / rsSD[sortedrsIndex0[t][sortedrsIndex[t][r]]];
    Matrix netCoef(toprsCount+1, 3);
    StringArray tmpcolNames;
    outfile = "tmp_"; outfile.Add(outfile0); outfile.Add("_net_coef.txt");
    ReadMatrix(outfile, tmpcolNames, netCoef, coreCount);
    for(int r = 0; r < toprsCount; r++)
      for(int c = 0; c < 3; c++)
        weights[r][c+2] = netCoef[r+1][c] / rsSD[sortedrsIndex0[t][r]];
    outfile = outfile0; outfile.Add("_PGS.mod");
    WriteTable(outfile, weights, header, col0);
    printf("  %d PGS models for trait %s write to %s\n",
      modelCount, (const char*)alltraitNames[traits[t]], (const char*)outfile);
  } // End of t loop
  delete []sortedrsIndex0;
  delete []sortedrsIndex;
  printf("\nSummary of %d Trait Prediction (without CV)\n", traitCount);
  printf("%20s%10s%10s%10s%10s%10s%10s\n",
    "TRAIT", "N", "RSQ_Best", "RSQ_WTSUM", "RSQ_Ridge", "RSQ_Net", "RSQ_Lasso");
  for (int t = 0; t < traitCount; t++)
    printf("%20s%10d%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
    (const char*)alltraitNames[traits[t]], SS[t], RsqBest[t][t], RsqWS[t],
    rsqs[0][t][0], rsqs[0][t][1], rsqs[0][t][2]);
  if(cvFlag){
    printf("\nSummary of %d Trait Prediction with CV\n", traitCount);
    printf("%20s%10s%10s%10s%10s%10s%10s\n",
      "TRAIT", "N", "RSQ_Best", "RSQ_WTSUM", "RSQ_Ridge", "RSQ_Net", "RSQ_Lasso");
    for(int t = 0; t < traitCount; t++)
      printf("%20s%10d%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
      (const char*)alltraitNames[traits[t]], SS[t], meanBest[t], meanWS[t],
      rsqs[1][t][0], rsqs[1][t][1], rsqs[1][t][2]);
  }
  for(int i = 0; i < idCount; i++)
    delete []masks[i];
  delete []masks;
  for(int r = 0; r < rsCount; r++)
    delete []LD[r];
  delete []LD;
  return;
}


void WriteRS(String &outfile, StringArray &topID, const char *traitName, Vector &toptrait, StringArray &toprsNames, Matrix &toprisks, Matrix &PCs, Vector &pcars)
{
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  int idCount = toprisks.rows;
  int rsCount = toprisks.cols;
  int dimPC = PCs.rows;
  fprintf(fp, "ID\t%s", traitName);
  for (int r = 0; r < rsCount; r++)
        fprintf(fp, "\t%s", (const char*)toprsNames[r]);
  for (int r = 0; r < dimPC; r++)
        fprintf(fp, "\tPC%d", r + 1);
  fprintf(fp, "\tPCARS");
  for (int i = 0; i < idCount; i++) {
        fprintf(fp, "\n%s\t%lf", (const char*)topID[i], toptrait[i]);
        for (int r = 0; r < rsCount; r++)
            fprintf(fp, "\t%.3lf", toprisks[i][r]);
        for (int r = 0; r < dimPC; r++)
            fprintf(fp, "\t%.4lf", PCs[r][i]);
        fprintf(fp, "\t%.4lf", pcars[i]);
  }   // End of i loop
  fprintf(fp, "\n");
  fclose(fp);
}


  /*
  if(cvFlag && nTop && ldT!=_NAN_){ // Cross-validation
    if(fold==0) fold=8;
    if(repeat==0) repeat=10;
    int iwordCount = (idCount-1)/64+1;
    unsigned long long int **testFlags = new unsigned long long int *[repeat];
    for(int run = 0; run < repeat; run++)
      testFlags[run] = new unsigned long long int[iwordCount];
    Matrix RsqBest(repeat, traitCount), RsqWS(repeat, traitCount);
    printf("Cross-validation start");
    ScanCV(risks, allY, masks, LD, fold, RsqBest, toprsCount, RsqWS, testFlags, coreCount);
    Vector meanBest(traitCount), sdBest(traitCount), meanWS(traitCount), sdWS(traitCount), V(repeat);
    for(int t = 0; t < traitCount; t++){
      for(int w = 0; w < repeat; w++)
        V[w] = RsqBest[w][t];
      meanBest[t] = V.Average();
      sdBest[t] = sqrt(V.Var());
      for(int w = 0; w < repeat; w++)
        V[w] = RsqWS[w][t];
      meanWS[t] = V.Average();
      sdWS[t] = sqrt(V.Var());
    } // End of t loop
    printf("\nSummary of %d Trait Prediction\n", traitCount);
    printf("%20s%10s%13s%13s%13s%13s\n", "TRAIT", "N_Total", "RSQBEST_MEAN", "RSQBEST_SD", "RSQWS_MEAN", "RSQWS_SD");
    for (int t = 0; t < traitCount; t++) {
      printf("%20s%10d%13.3lf%13.3lf%13.3lf%13.3lf",
      (const char*)alltraitNames[traits[t]], SS[t], meanBest[t], sdBest[t], meanWS[t], sdWS[t]);
      printf("\n");
    } // End of t loop
    if(!nowriteFlag){
      String outfile = prefix; outfile.Add("_CV_testset.txt");
      StringArray header(repeat+1);
      header[0] = "SampleID";
      String tS;
      for(int c = 0; c < repeat; c++){
        header[c+1] = "Run";
        tS = c+1;
        header[c+1].Add(tS);
      } // End of c loop
      WriteRotatedBinaryTable(outfile, testFlags, header, sampleID);
      printf("Masks of test set samples from %d cross-validations saved in file %s\n", repeat, (const char*)outfile);
      outfile = prefix; outfile.Add("_CV_RsqWTSUM.txt");
      header.Dimension(traitCount+1);
      header[0] = "Run";
      for(int t = 0; t < traitCount; t++)
        header[t+1] = alltraitNames[traits[t]];
      StringArray col0(repeat);
      for(int r = 0; r < repeat; r++)
        col0[r] = r+1;
      WriteTable(outfile, RsqWS, header, col0);
      printf("Rsq_WTSUM from %d cross-validations saved in file %s\n", repeat, (const char*)outfile);
      outfile = prefix; outfile.Add("_CV_RsqBest.txt");
      WriteTable(outfile, RsqBest, header, col0);
      printf("Rsq_Best from %d cross-validations saved in file %s\n", repeat, (const char*)outfile);
    }
    for(int run = 0; run < repeat; run++)
      delete[]testFlags[run];
    delete[]testFlags;

    const int PGSCOUNT=50;
    Vector Y;
    IntArray pheno;
    Matrix pgsQuantile(30, PGSCOUNT);
    Matrix rsqs[2];
    for(int k = 0; k < 2; k++)
      rsqs[k].Dimension(traitCount, 3);
    for(int t = 0; t < traitCount; t++){
      printf("Trait %s...\n", (const char*)alltraitNames[traits[t]]);
      int tw = (t>>6);
      int tb = (t&0x3F);
      Y.Dimension(0);
      pheno.Dimension(0);
      for(int i = 0; i < idCount; i++)
        if((masks[i][tw]>>tb) & 0x1) {
          Y.Push(allY[i][t]);
          pheno.Push(i);
        }// End of i loop
      int totalCount = Y.Length();
      String outfile0 = prefix;
      outfile0.Add("_");
      outfile0.Add(alltraitNames[traits[t]]);
      outfile = outfile0; outfile.Add(".linear");
      WriteLinear(outfile, sortedrsIndex0[t], rsNames, annotation, GTs[t], sigmaGs[t], sigmaTs[t], SS[t]);
      printf("  PGS scan statistics write to %s\n", (const char*)outfile);

      if(netFlag){
        Vector rs_WS(totalCount); rs_WS.Zero();
        int modelrsCount = sortedrsIndex[t].Length();
        for(int r = 0; r < modelrsCount; r++){
          int rr = sortedrsIndex0[t][sortedrsIndex[t][r]];
          double beta = GTs[t][rr] / sigmaGs[t][rr] / SS[t];
          for(int i = 0; i < totalCount; i++)
            rs_WS[i] += risks[rr][pheno[i]] * beta;
        } // End of r loop
        Matrix pgsTable(totalCount, toprsCount+2); pgsTable.Zero();
        for(int i = 0; i < totalCount; i++){
          pgsTable[i][0] = Y[i];
          for(int r = 0; r < toprsCount; r++)
            pgsTable[i][r+1] = risks[sortedrsIndex0[t][r]][pheno[i]];
          pgsTable[i][toprsCount+1] = rs_WS[i];
        } // End of i loop
        header.Dimension(toprsCount+3);
        header[0] = "ID";
        header[1] = alltraitNames[traits[t]];
        for(int r = 0; r < toprsCount; r++)
          header[r+2] = rsNames[sortedrsIndex0[t][r]];
        header[toprsCount+2] = "WTSUM";
        col0.Dimension(totalCount);
        for(int i = 0; i < totalCount; i++)
          col0[i] = sampleID[pheno[i]];
        outfile = outfile0; outfile.Add("_PGS.txt");
        WriteTable(outfile, pgsTable, header, col0);
        printf("  Top %d polygenic scores for trait %s write to %s\n",
          toprsCount, (const char*)alltraitNames[traits[t]], (const char*)outfile);
        // Ridge / elastic net / lasso are now fitted in-process by
        // Penalized.cpp, so --net no longer needs R or the glmnet package.
        // The two tmp_* files keep the same names and format, so the
        // reader just below is unchanged.
        PenalizedRegressionR2(outfile0, alltraitNames[traits[t]], repeat, coreCount, msg);
        printf("  %s", (const char*)msg);

        StringArray lines;
        outfile = "tmp_"; outfile.Add(outfile0); outfile.Add("_net_rsq.txt");
        ReadLines(outfile, lines);
        StringArray tokens;
        for(int m = 0; m < 3; m++){
          tokens.Clear();
          tokens.AddTokens(lines[m], '\t');
          for(int k = 0; k < 2; k++)
            rsqs[k][t][m] = tokens[k].AsDouble();
        } // End of m loop
      }
    } // End of t loop
    printf("\nSummary of %d Trait Prediction\n", traitCount);
    printf("%20s%10s%10s%10s%10s%10s%10s%13s%13s%13s\n",
      "TRAIT", "N", "RSQ_Best", "RSQ_WTSUM", "RSQ_Ridge", "RSQ_Net",
      "RSQ_Lasso", "RSQ_CV_Ridge", "RSQ_CV_Net", "RSQ_CV_Lasso" );
    for (int t = 0; t < traitCount; t++)
      printf("%20s%10d%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%13.3lf%13.3lf%13.3lf\n",
      (const char*)alltraitNames[traits[t]], SS[t], RsqBest[t][t], RsqWS[t],
      rsqs[0][t][0], rsqs[0][t][1], rsqs[0][t][2], rsqs[1][t][0], rsqs[1][t][1], rsqs[1][t][2]);


  }else if(cvFlag){ // autoCV
    if(fold==0) fold=8;
    if(repeat==0) repeat=20;
    IntArray maxPar(traitCount);
    Vector maxR2(traitCount), maxRsq(traitCount), maxVar(traitCount);
    AutoCV(risks, allY, masks, nTop, ldT, fold, repeat, maxPar, maxR2, maxRsq, maxVar, coreCount);
    printf("\nSummary of %d Trait Prediction with auto cross-validation\n", traitCount);
    printf("%20s%10s%10s%10s%13s%13s\n", "TRAIT", "N_Total", "Ntop", "R2", "RsqWS_Max", "RsqWS_SD");
    for (int t = 0; t < traitCount; t++)
      printf("%20s%10d%10d%10.3lf%13.3lf%13.3lf\n",
      (const char*)alltraitNames[traits[t]], SS[t], maxPar[t], maxR2[t], maxRsq[t], sqrt(maxVar[t]));
  }else if((!netFlag) && nowriteFlag){ // no Cross-Validation
    Matrix RsqBest(traitCount, traitCount);
    Vector RsqWS(traitCount);
    IntArray indexBest(traitCount);
    NowriteScan(risks, allY, masks, LD, RsqBest, indexBest, toprsCount, RsqWS, coreCount);
    printf("\nSummary of %d Trait Prediction\n", traitCount);
    printf("%20s%10s%10s%10s  PGS_Best\n", "TRAIT", "N", "RSQ_Best", "RSQ_WTSUM");
    for (int t = 0; t < traitCount; t++) {
      printf("%20s%10d%10.3lf%10.3lf  %s",
      (const char*)alltraitNames[traits[t]], SS[t], RsqBest[t][t], RsqWS[t], (const char*)rsNames[indexBest[t]]);
      if(annotationCount) printf(":%s", (const char*)annotation[indexBest[t]]);
      printf("\n");
    } // End of t loop
  }else{ // No CV, write
  */


