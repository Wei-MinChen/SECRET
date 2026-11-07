////////////////////////////////////////////////////////////// ////////
// RiskPrediction.cpp
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
// Oct 2, 2026

#include <math.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#include "Random.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

extern int plotPGSCorrList(const char *outprefix, const char *statfile, StringArray &inputs, StringArray &titles, const char *rpath, String &msg);
extern int plotAUCs(const char *input, const char *output, const char *rpath);
extern int plotTuneTraitAUCs(const char *prefix, const char *tuneTrait, const char *rpath);
extern void WritePGSCorrelations(const char *outfile0, Matrix &risks, StringArray &rsNames, IntArray &topIndex,
  IntArray &samples, IntArray &status, int coreCount, bool &caseControlDone, double pruneR2);

void ComputeAUCs(Matrix &risks, unsigned long long int **btraits, Matrix &AUC, int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
  int affectionCount = AUC.cols;
  const int IBLOCKSIZE = 255;
  int blockCount = (idCount-1)/IBLOCKSIZE+1;
  int AFFBLOCKSIZE = affectionCount < 32? affectionCount: 32;
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
{
#endif
  unsigned char thisblockCount[2][32];
  unsigned char *pCount[2][32], *perblockCount[2][32]; // 2^14 for N=64K
  IntArray perblockCount4[2];
  for(int k = 0; k < 2; k++){
    for(int a = 0; a < AFFBLOCKSIZE; a++){
      pCount[k][a] = new unsigned char[idCount+1];
      perblockCount[k][a] = new unsigned char[blockCount];
    } // End of a loop
    perblockCount4[k].Dimension(blockCount+1);
    perblockCount4[k].Zero();
  } // End of k loop
  Vector sensitivity(idCount+1), specificity(idCount+1);
  double auc[32];
#ifdef _OPENMP
#pragma omp for
#endif
  for (int r = 0; r < rsCount; r++) {
    QuickIndex rIdx(risks[r]);
    for(int a = 0; a < affectionCount; a += 32){
      int aMax = a < affectionCount-32? 32: affectionCount-a;
      for(int cc = 0; cc < 2; cc++)
        for(int aa = 0; aa < aMax; aa++)
          memset(&pCount[cc][aa][0], 0, idCount+1);
      int w = (a>>5);
      for(int i = 0; i < idCount; i+=IBLOCKSIZE){ // CPU CACHE BLOCK
        int iMax = i < idCount-IBLOCKSIZE? i+IBLOCKSIZE: idCount;
        for(int cc = 0; cc < 2; cc++)
          memset(&thisblockCount[cc][0], 0, AFFBLOCKSIZE);
        for(int ii = i; ii < iMax; ii++){
          unsigned long long int word = btraits[rIdx[ii]][w];
          for(int aa = 0; aa < aMax; aa++, word >>= 2)
            if(word & 0x2){
              int cc = (word & 0x1);
              for(int cut = i; cut <= ii; cut++)
                pCount[cc][aa][cut] ++;
              thisblockCount[cc][aa]++;
            } // End of aa loop
        } // End of ii loop
        int block = i/IBLOCKSIZE;
        for(int cc = 0; cc < 2; cc++)
          for(int aa = 0; aa < aMax; aa++)
            perblockCount[cc][aa][block] = thisblockCount[cc][aa];
      } // End of CPU CACHE BLOCK
      for(int aa = 0; aa < aMax; aa++){
        for(int cc = 0; cc < 2; cc ++)
          for(int bb = 0; bb < blockCount; bb++){
              int count = 0;
              for(int block = bb; block < blockCount; block++)
                count += perblockCount[cc][aa][block];
              perblockCount4[cc][bb] = count;
            } // End of bb loop
        if(perblockCount4[0][0])
          for(int cut = 0; cut < idCount+1; cut++)
            specificity[cut] = 1 - (perblockCount4[0][cut/IBLOCKSIZE+1]+pCount[0][aa][cut]) * 1.0 / perblockCount4[0][0];
        else memset(&specificity[0], 1, (idCount+1)*sizeof(double));
        if(perblockCount4[1][0])
          for(int cut = 0; cut < idCount+1; cut++)
            sensitivity[cut] = (perblockCount4[1][cut/IBLOCKSIZE+1]+pCount[1][aa][cut]) * 1.0 / perblockCount4[1][0];
        else memset(&sensitivity[0], 1, (idCount+1)*sizeof(double));
        // trapezoids between tie-group boundaries only, so tied scores count as 0.5
        // and the AUC no longer depends on the order of samples with identical scores
        double tauc = 0.0;
        int prevCut = 0;
        for (int cut = 1; cut <= idCount; cut++)
          if (cut == idCount || risks[r][rIdx[cut-1]] != risks[r][rIdx[cut]]) {
            tauc += (specificity[cut] - specificity[prevCut]) * (sensitivity[prevCut] + sensitivity[cut]);
            prevCut = cut;
          }
        auc[aa] = tauc*0.5;
      } // End of aa loop
      memcpy(&AUC[r][a], &auc[0], aMax*sizeof(double));
    } // End of a loop
  } // End of OMP r loop
  for(int k = 0; k < 2; k++)
    for(int aa = 0; aa < AFFBLOCKSIZE; aa++)
      delete []pCount[k][aa];
#ifdef _OPENMP
 } // extra bracket for omp
#endif
}

void maketestFlags(int repeat, int fold, int idCount, unsigned long long int **testFlags)
{
  IntArray testFlag(idCount);
  IntArray *testFlag64 = new IntArray[idCount];
  for(int i = 0; i < idCount; i++)
    testFlag64[i].Dimension(64);
  Random ran;
  int testsetCount = idCount * 1.0 / fold;
  for(int run = 0; run < repeat; run += 64){
    int rMax = run<repeat-64? 64: repeat-run;
    int w = (run>>6);
    for(int r = 0; r < rMax; r++){
      ran.Choose(testFlag, idCount, testsetCount);
      for(int i = 0; i < idCount; i++)
        testFlag64[i][r] = testFlag[i];
    } // End of r loop
    for(int i = 0; i < idCount; i++){
      unsigned long long int word = 0;
      for(int r = 0; r < rMax; r++)
        word |= ((unsigned long long int)testFlag64[i][r]<<r);
      testFlags[i][w] = word;
    } // End of i loop
  } // End of run loop
  delete []testFlag64;
}

void ComputeAUCsCV(Matrix &risks, unsigned long long int **btraits, int repeat, int fold,
  Vector &meanAUC, Vector &sdAUC, Vector &freqValidated, IntArray &validatedIndex, int coreCount)
{
  if(fold == 0) fold = 8;
  if(repeat == 0) repeat = 20;
  int rsCount = risks.rows;
  int idCount = risks.cols;
  int affectionCount = meanAUC.Length();
  const int IBLOCKSIZE = 255;
  int blockCount = (idCount-1)/IBLOCKSIZE+1;
  Matrix AUC(rsCount, affectionCount*repeat*2);
  int wordCount = (repeat-1)/64+1;
  unsigned long long int **testFlags = new unsigned long long int*[idCount];
  for(int i = 0; i < idCount; i++)
    testFlags[i] = new unsigned long long int[wordCount];
  maketestFlags(repeat, fold, idCount, testFlags);
  int RUNBLOCKSIZE = repeat < 64? repeat: 64;
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
{
#endif
  unsigned char thisblockCount[2][2][64];
  unsigned char *pCount[2][2][64], *perblockCount[2][2][64]; // 2^16 for N=64K
  IntArray perblockCount4[2];
  for(int cc = 0; cc < 2; cc++){
    for(int k = 0; k < 2; k++)
      for(int ww = 0; ww < RUNBLOCKSIZE; ww++){
        pCount[cc][k][ww] = new unsigned char[idCount+1];
        perblockCount[cc][k][ww] = new unsigned char[blockCount];
      } // End of ww loop
    perblockCount4[cc].Dimension(blockCount+1);
    perblockCount4[cc].Zero();
  } // End of cc loop
  Vector sensitivity(idCount+1), specificity(idCount+1);
  double auc[2][64];
  int iwordCount = (idCount-1)/32+1;
  unsigned long long int *btrait = new unsigned long long int[iwordCount];
#ifdef _OPENMP
#pragma omp for
#endif
  for (int r = 0; r < rsCount; r++) {
    QuickIndex rIdx(risks[r]);
    for (int a = 0; a < affectionCount; a++){
      int aw = (a>>5);
      int as = ((a&0x1F)<<1);
      for(int i = 0; i < idCount; i+=32){
        int iMax = i < idCount-32? i+32: idCount;
        int iw = (i>>5);
        btrait[iw] = 0;
        for(int ii = i; ii < iMax; ii++)
          btrait[iw] |= (((btraits[rIdx[ii]][aw]>>as) & 0x3)<<((ii-i)<<1));
      } // End of i loop
      for(int run = 0; run < repeat; run+=64){
        int wMax = run < repeat-64? 64: repeat-run;
        for(int cc = 0; cc < 2; cc++)
          for(int k = 0; k < 2; k++)
            for(int ww = 0; ww < wMax; ww++)
              memset(&pCount[cc][k][ww][0], 0, idCount+1);
        int w = (run >> 6);
        for(int i = 0; i < idCount; i+=IBLOCKSIZE){ // CPU CACHE BLOCK
          for(int cc = 0; cc < 2; cc++)
            for(int k = 0; k < 2; k++)
              memset(&thisblockCount[cc][k][0], 0, RUNBLOCKSIZE);
          int iMax = i < idCount - IBLOCKSIZE? i+IBLOCKSIZE: idCount;
          for(int ii = i; ii < iMax; ii++){ // 64KB
            char pheno = (btrait[ii>>5]>>((ii&0x1F)<<1));
            if(pheno&0x2){ // non-missing
              int cc = (pheno&0x1);
              unsigned long long int word = testFlags[rIdx[ii]][w];
              for(int ww = 0; ww < wMax; ww++, word >>= 1)
                if((word & 0x1)==0){ // training set
                  for(int cut = i; cut <= ii; cut++)
                    pCount[cc][0][ww][cut]++;
                  thisblockCount[cc][0][ww]++;
                }else{ // test set
                  for(int cut = i; cut <= ii; cut++)
                    pCount[cc][1][ww][cut]++;
                  thisblockCount[cc][1][ww]++;
                } // End of ww loop
            } // End of if pheno
          } // End of ii loop
          int block = i/IBLOCKSIZE;
          for(int cc = 0; cc < 2; cc++)
            for(int k = 0; k < 2; k++)
              for(int ww = 0; ww < wMax; ww++)
                perblockCount[cc][k][ww][block] = thisblockCount[cc][k][ww];
        } // End of CPU CACHE BLOCK
        for(int k = 0; k < 2; k++)
          for(int ww = 0; ww < wMax; ww++){
            for(int cc = 0; cc < 2; cc ++)
              for(int bb = 0; bb < blockCount; bb++){
                int count = 0;
                for(int block = bb; block < blockCount; block++)
                  count += perblockCount[cc][k][ww][block];
                perblockCount4[cc][bb] = count;
              } // End of bb loop
            if(perblockCount4[0][0])
              for(int cut = 0; cut < idCount+1; cut++)
                specificity[cut] = 1-(perblockCount4[0][cut/IBLOCKSIZE+1]+pCount[0][k][ww][cut])*1.0/ perblockCount4[0][0];
            else
              memset(&specificity[0], 1, (idCount+1)*sizeof(double));
            if(perblockCount4[1][0])
              for(int cut = 0; cut < idCount+1; cut++)
                sensitivity[cut] = (perblockCount4[1][cut/IBLOCKSIZE+1]+pCount[1][k][ww][cut])*1.0 / perblockCount4[1][0];
            else
              memset(&sensitivity[0], 1, (idCount+1)*sizeof(double));
            double tauc = (specificity[1] - specificity[0])*sensitivity[0];
            for(int cut = 0; cut < idCount-1; cut++)
              tauc += (specificity[cut+2] - specificity[cut]) * sensitivity[cut+1];
            tauc += (1 - specificity[idCount-1])*sensitivity[idCount];
            auc[k][ww] = tauc * 0.5;
          } // End of k&ww loop
        for(int k = 0; k < 2; k++)
          memcpy(&AUC[r][a*repeat*2+k*repeat+run], &auc[k][0], wMax*sizeof(double));
      } // End of run loop
    } // End of a loop
  } // End of OMP r loop
  for(int cc = 0; cc < 2; cc++)
    for(int k = 0; k < 2; k++){
      for(int ww = 0; ww < RUNBLOCKSIZE; ww++){
        delete []pCount[cc][k][ww];
        delete []perblockCount[cc][k][ww];
      } // End of ww loop
    } // End of k loop
  delete []btrait;
#ifdef _OPENMP
 } // extra bracket for omp
#endif
  for(int i = 0; i < idCount; i++)
    delete []testFlags[i];
  delete []testFlags;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int a = 0; a < affectionCount; a++){
    Vector bestAUC(repeat);
    IntArray bestIndex(repeat);
    for(int run = 0; run < repeat; run++){
      double max = -1;
      int maxIndex = -1;
      int pos = a*repeat*2 + run;
      for(int r = 0; r < rsCount; r++)
        if(AUC[r][pos] > max){
          max = AUC[r][pos];
          maxIndex = r;
        } // End of r loop
      bestAUC[run] = AUC[maxIndex][pos+repeat];
      bestIndex[run] = maxIndex;
    } // End of run loop
    meanAUC[a] = bestAUC.Average();
    sdAUC[a] = sqrt(bestAUC.Var());
    QuickIndex rIdx(bestIndex);
    validatedIndex[a] = bestIndex[rIdx[repeat/2]];
    freqValidated[a] = 0.0;
    for(int run = 0; run < repeat; run++)
      if(bestIndex[run] == validatedIndex[a])
        freqValidated[a] ++;
    freqValidated[a] /= repeat;
  } // End of OMP a loop
}

// ordering of (-rounded AUC, PGS) pairs by AUC (a function, not a lambda: C++ Builder 6 has no lambdas)
static bool RankedLess(const std::pair<double, int> &x, const std::pair<double, int> &y)
{
  return x.first < y.first;
}

// The PGSs of column d of the AUC table with the largest AUC first: ranked gets (-AUC rounded to the
// 5 digits written to the .auc file, PGS), ascending, so that PGSs with the same rounded AUC stay in
// the order of the PGSs.  invalid[r] != 0 (a PGS that is excluded) and a missing AUC are left out.
// Both the .auc file and the correlations of the top PGSs use it, so they rank the PGSs the same way.
static void RankPGSByAUC(Matrix &auc, int d, const std::vector<char> &invalid,
  std::vector<std::pair<double, int> > &ranked)
{
  const int rsCount = auc.rows;
  ranked.clear();
  for (int r = 0; r < rsCount; r++) {
    double a = auc[r][d];
    if (invalid[r] || a == _NAN_ || a != a) continue;
    ranked.push_back(std::make_pair(-floor(a * 1e5 + 0.5), r));
  }   // End of r loop
  std::stable_sort(ranked.begin(), ranked.end(), RankedLess);
}

// The .auc file of one trait, column d of the AUC table: one row for each PGS, the largest AUC first
// (ties in the order of the PGSs).  The PGSs in invalidPGS (excluded) are not written; a PGS without
// an AUC is written last, with NA.  Columns: PGS, the trait (its AUC), SNPused and Coverage (only with
// a coverage dictionary), Annotation (only with a dictionary).
void WriteAUC(FILE *fp, const char *traitName, int d, StringArray &rsNames,
  StringArray &annotation, IntArray &invalidPGS, Matrix &auc)
{
    int rsCount = rsNames.Length();
    int annotCount = annotation.Length();
    std::vector<char> invalid(rsCount, 0);
    for (int k = 0; k < invalidPGS.Length(); k++)
        invalid[invalidPGS[k]] = 1;
    std::vector<std::pair<double, int> > ranked;
    RankPGSByAUC(auc, d, invalid, ranked);
    std::vector<int> order;                       // the PGSs as they are written
    order.reserve(rsCount);
    for (size_t k = 0; k < ranked.size(); k++)
        order.push_back(ranked[k].second);
    for (int r = 0; r < rsCount; r++) {
        double a = auc[r][d];
        if (!invalid[r] && (a == _NAN_ || a != a)) order.push_back(r);
    }   // End of r loop: the PGSs without an AUC
    char buffer[0x20000];
    int pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "PGS\t%s", traitName);
    pbuffer += FormatCoverageHeader(&buffer[pbuffer]);    // SNPused, Coverage: only with a coverage dictionary
    if (annotCount)
        pbuffer += sprintf(&buffer[pbuffer], "\tAnnotation");
    for (size_t k = 0; k < order.size(); k++) {
        int r = order[k];
        double a = auc[r][d];
        pbuffer += sprintf(&buffer[pbuffer], "\n%s", (const char*)rsNames[r]);
        if (a == _NAN_ || a != a)
            pbuffer += sprintf(&buffer[pbuffer], "\tNA");
        else
            pbuffer += sprintf(&buffer[pbuffer], "\t%.5lf", a);
        pbuffer += FormatCoverageFields(&buffer[pbuffer], r);
        if (annotCount)
            pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)annotation[r]);
        if (pbuffer > 0xFFFF) {  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
        }
    }   // End of k loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
}

// Two bits per trait for ComputeAUCs: bit 1 = the trait is observed, bit 0 = affected (the
// traits are coded 1 = unaffected, 2 = affected).  columns are the positions in affections of
// the traits, in the order of the bits; a sample with sampleFlag[i] == 0 has no trait observed,
// so it is left out of the scan (sampleFlag NULL: every sample is used).
static unsigned long long int **PackTraits(Matrix &alltraits, IntArray &affections, IntArray &columns,
  IntArray *sampleFlag, int idCount, int coreCount)
{
  const int colCount = columns.Length();
  const int tDim = (colCount - 1) / 32 + 1;
  unsigned long long int **packed = new unsigned long long int *[idCount];
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for (int i = 0; i < idCount; i++) {
    packed[i] = new unsigned long long int[tDim];
    memset(&packed[i][0], 0, tDim * sizeof(unsigned long long int));
    if (sampleFlag != NULL && (*sampleFlag)[i] == 0) continue;
    for (int c = 0; c < colCount; c++) {
      double v = alltraits[i][affections[columns[c]]];
      if (v != _NAN_)
        packed[i][c >> 5] |= ((unsigned long long int)(v + 1.5) << ((c & 0x1F) << 1));
    }   // End of c loop
  }   // End of i loop
  return packed;
}

static void FreePackedTraits(unsigned long long int **packed, int idCount)
{
  for (int i = 0; i < idCount; i++)
    delete []packed[i];
  delete []packed;
}

// --auc --rplot: correlations among the top 30 PGSs (by AUC) of one trait, column d of the AUC table,
// for the heatmap of the top PGSs made by plotAUCs and for the PGS correlation plots.  The PGSs are
// ranked as in the .auc file (RankPGSByAUC), and the correlations use the samples that have the
// trait.  tprefix is <prefix>_<trait> (for a stratum <prefix>_<trait>_<strat>_<1 or 2>); the files are
//     <tprefix>_AUC_PGScorr.txt (with _case, _control and _diff) and <tprefix>_AUC_PGScorrplot.pdf
// d is also the position of the trait in btraits.
static void WriteAUCCorrelations(const char *tprefix, const char *traitName, int d, Matrix &auc,
  IntArray &invalidPGS, Matrix &risks, StringArray &rsNames, unsigned long long int **btraits,
  const char *rpath, int coreCount)
{
  const int rsCount = auc.rows, idCount = risks.cols;
  std::vector<char> invalid(rsCount, 0);
  for (int k = 0; k < invalidPGS.Length(); k++)
    invalid[invalidPGS[k]] = 1;
  std::vector<std::pair<double, int> > ranked;   // (-rounded AUC, PGS), sorted ascending
  RankPGSByAUC(auc, d, invalid, ranked);
  int top = (int) ranked.size() < 30 ? (int) ranked.size() : 30;
  if (top < 2) return;
  IntArray topIndex(top);
  for (int k = 0; k < top; k++)
    topIndex[k] = ranked[k].second;
  IntArray samples(0);
  for (int i = 0; i < idCount; i++)
    if ((btraits[i][d >> 5] >> ((d & 0x1F) << 1)) & 0x2)   // as ComputeAUCs: bit 1 = the trait is observed
      samples.Push(i);
  if (samples.Length() < 10) return;
  String outfile0(tprefix);
  outfile0.Add("_AUC");
  IntArray noStatus(0);
  bool caseControlDone = false;
  WritePGSCorrelations(outfile0, risks, rsNames, topIndex, samples, noStatus, coreCount, caseControlDone, _NAN_);
  // the correlation plots of the trait; the .auc file gives the annotation of the PGSs
  String statfile(tprefix), msg(256);
  statfile.Add(".auc");
  StringArray corrInputs, corrTitles;
  corrInputs.Push(outfile0);
  corrTitles.Push(traitName);
  plotPGSCorrList(outfile0, statfile, corrInputs, corrTitles, rpath, msg);
  printf("%s", (const char*)msg);
}

// --auc: everything that is written for one trait, column d of the AUC table (and of btraits):
//     <tprefix>.auc, the PGSs with the largest AUC first, and with --rplot the correlations of the
//     top PGSs and the figures (<tprefix>_aucplot, <tprefix>_AUC_PGScorrplot)
// tprefix is <prefix>_<trait>, or <prefix>_<trait>_<strat>_<1 or 2> for a stratum, and where is the
// text for the message (empty, or " in <strat>=<1 or 2>").
static void WriteAUCTrait(const char *tprefix, const char *traitName, int d, Matrix &auc,
  IntArray &invalidPGS, StringArray &rsNames, StringArray &annotation, Matrix &risks,
  unsigned long long int **btraits, const char *where, bool rplotFlag, const char *rpath, int coreCount)
{
  String outfile(tprefix);
  outfile.Add(".auc");
  FILE *fp = fopen(outfile, "wb");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
  WriteAUC(fp, traitName, d, rsNames, annotation, invalidPGS, auc);
  fclose(fp);
  printf("AUCs for %d PGS scores of trait %s%s write to %s\n", rsNames.Length() - invalidPGS.Length(),
    traitName, where, (const char*)outfile);
  if (rplotFlag) {
    WriteAUCCorrelations(tprefix, traitName, d, auc, invalidPGS, risks, rsNames, btraits, rpath, coreCount);
    plotAUCs(outfile, tprefix, rpath);
  }
}

void predictionAUC(const char *phefile, const char *dictfile, const char *rsfile,
  bool cvFlag, int fold, int repeat,
  const char *pcarsfile, const char *prefix, const char *tuneTrait,
  const char *covariateList, bool rplotFlag, const char *rpath, int coreCount)
{
  printf("\nOptions in effect:\n");
  printf("\t--auc\n");
  if(cvFlag) printf("\t--cv\n");
  printf("\t--pheno %s\n", phefile);
  if (gCovarFile[0]) printf("\t--covarfile %s\n", gCovarFile);
  if (gTraitList[0]) printf("\t--trait %s\n", gTraitList);
  if (covariateList[0])
    printf("\t--covariate %s\n", covariateList);
  if (tuneTrait[0])
    printf("\t--strat %s\n", tuneTrait);
  printf("\t--score %s\n", rsfile);
  if(fold) printf("\t--fold %d\n", fold);
  if(repeat) printf("\t--repeat %d\n", repeat);
  if(pcarsfile[0])
    printf("\t--readpcars %s\n", pcarsfile);
  if(dictfile[0])
    printf("\t--diction %s\n", dictfile);
  if(rplotFlag)
    printf("\t--rplot\n");
  if(coreCount)
    printf("\t--cpus %d\n", coreCount);
  if(prefix)
    printf("\t--prefix %s\n", prefix);
  printf("\n");

  StringArray sampleID, alltraitNames, rsNames, annotation;
  IntArray covariates, affections, quantraits;
  Matrix alltraits, risks;
  SECRET_Read(phefile, dictfile, rsfile, pcarsfile, prefix, covariateList,
    sampleID, alltraitNames, alltraits, covariates, affections, quantraits,
    rsNames, risks, annotation, coreCount);
  int alltraitCount = alltraitNames.Length();
  int idCount = sampleID.Length();
  if (idCount == 0) return;
  int rsCount = rsNames.Length();
  int affectionCount = affections.Length();
  printf("Computing AUCs..."); fflush(stdout);
  unsigned long long int **btraits = new unsigned long long int *[idCount];
  StringArray affectionNames(0);
  for (int a = 0; a < affectionCount; a++)
    affectionNames.Push(alltraitNames[affections[a]]);
  int tuneIndex = affectionNames.Find(tuneTrait);
  int tuneIndex0 = tuneIndex > -1 ? affections[tuneIndex] : -1;
  if (tuneTrait[0] != '\0')   // stratum files of an earlier run must not be mistaken for this run's
    for (int g = 1; g <= 2; g++) {
      String oldFile(prefix), number;
      number = g;
      oldFile.Add("_"); oldFile.Add(tuneTrait); oldFile.Add("_"); oldFile.Add(number); oldFile.Add(".auc");
      remove(oldFile);
    }   // End of g loop
  if (tuneIndex < 0 && tuneTrait[0] != '\0')
    printf("The --strat trait %s is not among the dichotomous traits in %s.\n", tuneTrait, phefile);

  IntArray invalidPGS(0);

  if (tuneIndex > -1) {
    // --strat: the AUC scan is run on all samples, and then separately in each stratum of the
    // --strat trait (its value 1 and its value 2).  The files of every trait are
    //   <prefix>_<trait>.auc               all samples, every dichotomous trait but the --strat trait
    //   <prefix>_<trait>_<strat>_1.auc     the samples with --strat trait = 1, the traits that can be
    //   <prefix>_<trait>_<strat>_2.auc     evaluated in both strata
    // each with the PGSs of the largest AUC first.
    // A PGS that predicts the --strat trait itself (AUC > 0.65 in all samples) is excluded
    // from all of them.
    const int MINSTRATUM = 10;            // cases and controls needed in each stratum
    IntArray columns(0);                  // the traits of the tables; the --strat trait comes last
    for (int a = 0; a < affectionCount; a++)
      if (a != tuneIndex) columns.Push(a);
    const int tableCount = columns.Length();
    if (tableCount == 0) {
      printf("done.\nThe --strat trait %s is the only dichotomous trait, so there is nothing to compute.\n", tuneTrait);
      for (int i = 0; i < idCount; i++) btraits[i] = NULL;
      delete []btraits;
      return;
    }
    columns.Push(tuneIndex);
    unsigned long long int **btraitsAll = PackTraits(alltraits, affections, columns, NULL, idCount, coreCount);
    Matrix aucAll(rsCount, tableCount + 1);
    ComputeAUCs(risks, btraitsAll, aucAll, coreCount);
    printf("done.\n"); fflush(stdout);
    for (int r = 0; r < rsCount; r++)
      if (aucAll[r][tableCount] > 0.65) invalidPGS.Push(r);
    int invalidPGSCount = invalidPGS.Length();
    if (invalidPGSCount) {
      printf("The following PGS scores are excluded due to large AUC for the --strat trait %s:\n", tuneTrait);
      for (int r = 0; r < invalidPGSCount; r++) {
        printf("\t%s\t%.3lf", (const char*)rsNames[invalidPGS[r]], aucAll[invalidPGS[r]][tableCount]);
        if (dictfile[0]) printf("\t%s", (const char*)annotation[invalidPGS[r]]);
        printf("\n");
      }   // End of r loop
    }
    printf("Trait %s is used for the stratification and is not in the AUC tables.\n", tuneTrait);
    StringArray names(tableCount);
    Matrix aucTable(rsCount, tableCount);
    for (int c = 0; c < tableCount; c++) {
      names[c] = affectionNames[columns[c]];
      for (int r = 0; r < rsCount; r++)
        aucTable[r][c] = aucAll[r][c];
    }   // End of c loop
    for (int c = 0; c < tableCount; c++) {
      String tprefix(prefix);
      tprefix.Add("_");
      tprefix.Add(names[c]);
      WriteAUCTrait(tprefix, names[c], c, aucTable, invalidPGS, rsNames, annotation, risks, btraitsAll,
        "", rplotFlag, rpath, coreCount);
    }   // End of c loop
    FreePackedTraits(btraitsAll, idCount);

    // the strata: which traits have enough cases and controls in both
    IntArray inStratum[2];
    int stratumSize[2] = {0, 0};
    for (int g = 0; g < 2; g++) {
      inStratum[g].Dimension(idCount);
      inStratum[g].Zero();
      for (int i = 0; i < idCount; i++) {
        double v = alltraits[i][tuneIndex0];
        if (v != _NAN_ && fabs(v - (g + 1)) < 0.5) { inStratum[g][i] = 1; stratumSize[g]++; }
      }   // End of i loop
    }   // End of g loop
    printf("AUC scan in each stratum of %s: %s=1 (%d samples) and %s=2 (%d samples)\n",
      tuneTrait, tuneTrait, stratumSize[0], tuneTrait, stratumSize[1]);
    IntArray keep(0);                     // positions in names of the traits evaluated in both strata
    for (int c = 0; c < tableCount; c++) {
      bool ok = true;
      for (int g = 0; g < 2; g++) {
        int cases = 0, controls = 0;
        for (int i = 0; i < idCount; i++) {
          double v = alltraits[i][affections[columns[c]]];
          if (!inStratum[g][i] || v == _NAN_) continue;
          if (v > 1.5) cases++; else controls++;
        }   // End of i loop
        if (cases < MINSTRATUM || controls < MINSTRATUM) {
          printf("  trait %s is skipped: %s=%d has %d cases and %d controls (at least %d each are needed)\n",
            (const char*)names[c], tuneTrait, g + 1, cases, controls, MINSTRATUM);
          ok = false;
          break;
        }
      }   // End of g loop
      if (ok) keep.Push(c);
    }   // End of c loop
    int keepCount = keep.Length();
    if (keepCount == 0) {
      printf("No trait has enough cases and controls in both strata; the stratum-specific AUC files are not made.\n");
    } else {
      StringArray keepNames(keepCount);
      IntArray keepColumns(keepCount);
      for (int k = 0; k < keepCount; k++) {
        keepNames[k] = names[keep[k]];
        keepColumns[k] = columns[keep[k]];
      }   // End of k loop
      for (int g = 0; g < 2; g++) {
        unsigned long long int **btraitsG = PackTraits(alltraits, affections, keepColumns, &inStratum[g], idCount, coreCount);
        Matrix aucG(rsCount, keepCount);
        ComputeAUCs(risks, btraitsG, aucG, coreCount);
        String stratumNumber, stratumLabel(" in ");
        stratumNumber = g + 1;
        stratumLabel.Add(tuneTrait);
        stratumLabel.Add("=");
        stratumLabel.Add(stratumNumber);
        for (int k = 0; k < keepCount; k++) {
          String stratumPrefix(prefix);
          stratumPrefix.Add("_");
          stratumPrefix.Add(keepNames[k]);
          stratumPrefix.Add("_");
          stratumPrefix.Add(tuneTrait);
          stratumPrefix.Add("_");
          stratumPrefix.Add(stratumNumber);
          WriteAUCTrait(stratumPrefix, keepNames[k], k, aucG, invalidPGS, rsNames, annotation, risks, btraitsG,
            stratumLabel, rplotFlag, rpath, coreCount);
        }   // End of k loop
        FreePackedTraits(btraitsG, idCount);
      }   // End of g loop
      if (rplotFlag)       // the AUC of each trait in the two strata, from <prefix>_<trait>_<strat>_1.auc and _2.auc
        for (int k = 0; k < keepCount; k++) {
          String tprefix(prefix);
          tprefix.Add("_");
          tprefix.Add(keepNames[k]);
          plotTuneTraitAUCs(tprefix, tuneTrait, rpath);
        }   // End of k loop
    }
    for (int i = 0; i < idCount; i++) btraits[i] = NULL;   // the branch packed its own
    delete []btraits;
    return;
  }else{ // No tune traits

    int tDim = (affectionCount-1)/32+1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int i = 0; i < idCount; i++){
      btraits[i] = new unsigned long long int[tDim];
      memset(&btraits[i][0], 0, tDim*sizeof(unsigned long long int));
      for(int a = 0; a < affectionCount; a++)
        if(alltraits[i][affections[a]] != _NAN_)
          btraits[i][a>>5] |= ((unsigned long long int)(alltraits[i][affections[a]] + 1.5) << ((a&0x1F)<<1));
    } // End of i loop
    Matrix auc(rsCount, affectionCount);
    ComputeAUCs(risks, btraits, auc, coreCount);
    printf("done.\n"); fflush(stdout);
    for (int a = 0; a < affectionCount; a++) {
      String tprefix(prefix);
      tprefix.Add("_");
      tprefix.Add(affectionNames[a]);
      WriteAUCTrait(tprefix, affectionNames[a], a, auc, invalidPGS, rsNames, annotation, risks, btraits,
        "", rplotFlag, rpath, coreCount);
    }   // End of a loop

    // CV below
    Vector meanAUC(affectionCount), sdAUC(affectionCount), freqValidated(affectionCount);
    IntArray validatedIndex(affectionCount);
    Matrix aucTrain(rsCount, affectionCount), sdTrain(rsCount, affectionCount),
      aucTest(rsCount, affectionCount), sdTest(rsCount, affectionCount);
    ComputeAUCsCV(risks, btraits, repeat, fold, meanAUC, sdAUC, freqValidated, validatedIndex, coreCount);
    int annotCount = annotation.Length();
    printf("\n%15s %7s %7s %10s %11s", "Trait", "AUC_Ave", "AUC_SD", "PGS_Likely", "Consistency");
    if(annotCount)
      printf(" %20s", "Annotation");
    printf("\n");
    for(int a = 0; a < affectionCount; a++){
      printf("%15s %7.3lf %7.3lf %10s %11.3lf",
        (const char*)alltraitNames[affections[a]], meanAUC[a], sdAUC[a],
        (const char*)rsNames[validatedIndex[a]], freqValidated[a]);
      if(annotCount)
        printf(" %20s", (const char*)annotation[validatedIndex[a]]);
      printf("\n");
    } // End of a loop
  }
  for(int i = 0; i < idCount; i++)
    delete []btraits[i];
  delete []btraits;
}

  /*
  if(cvFlag){
    Vector meanAUC(affectionCount), sdAUC(affectionCount), freqValidated(affectionCount);
    IntArray validatedIndex(affectionCount);
    Matrix aucTrain(rsCount, affectionCount), sdTrain(rsCount, affectionCount),
      aucTest(rsCount, affectionCount), sdTest(rsCount, affectionCount);
    ComputeAUCsCV(risks, btraits, repeat, fold, meanAUC, sdAUC, freqValidated, validatedIndex, coreCount);

    int annotCount = annotation.Length();
    printf("\n%15s %7s %7s %10s %11s", "Trait", "AUC_Ave", "AUC_SD", "PGS_Likely", "Consistency");
    if(annotCount)
      printf(" %20s", "Annotation");
    printf("\n");
    for(int a = 0; a < affectionCount; a++){
      printf("%15s %7.3lf %7.3lf %10s %11.3lf",
        (const char*)alltraitNames[affections[a]], meanAUC[a], sdAUC[a],
        (const char*)rsNames[validatedIndex[a]], freqValidated[a]);
      if(annotCount)
        printf(" %20s", (const char*)annotation[validatedIndex[a]]);
      printf("\n");
    } // End of a loop
  }else{



    Matrix auc(rsCount, affectionCount), meanRS1, meanRS2;
    ComputeAUCs(risks, btraits, auc, coreCount);
    printf("done.\n"); fflush(stdout);

    IntArray invalidPGS(0);
    StringArray affectionNames(0);
    for (int a = 0; a < affectionCount; a++)
        affectionNames.Push(alltraitNames[affections[a]]);
    int tuneTraitIndex = affectionNames.Find(tuneTrait);
    int tuneTraitIndex0 = tuneTraitIndex > -1 ? affections[tuneTraitIndex] : -1;
    if (tuneTraitIndex < 0 && tuneTrait[0] != '\0')
        printf("The --strat trait %s is not among the dichotomous traits in %s.\n", tuneTrait, phefile);
    if (tuneTraitIndex > -1) {
        for (int r = 0; r < rsCount; r++)
            if (auc[r][tuneTraitIndex] > 0.65) invalidPGS.Push(r);
        int invalidPGSCount = invalidPGS.Length();
        if (invalidPGSCount) {
            printf("The following PGS scores are excluded due to large AUC for the --strat trait %s:\n", tuneTrait);
            for (int r = 0; r < invalidPGSCount; r++) {
                printf("\t%s\t%.3lf", (const char*)rsNames[invalidPGS[r]], auc[invalidPGS[r]][tuneTraitIndex]);
                if (dictfile[0]) printf("\t%s", (const char*)annotation[invalidPGS[r]]);
                printf("\n");
            }   // End of r loop
        }
        printf("Trait %s removed from AUC calculation.\n", (const char*)alltraitNames[affections[tuneTraitIndex]]);
        affections.Delete(tuneTraitIndex);
        affectionCount = affections.Length();
        auc.DeleteColumn(tuneTraitIndex);
    }
    String outfile;
    outfile.Copy(prefix);
    outfile.Add(".auc");
    FILE *fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    WriteAUC(fp, alltraitNames, affections, rsNames, annotation, invalidPGS, auc);
    fclose(fp);
    printf("AUCs for %d PGS scores and %d dichotomous traits write to %s\n", rsCount, affectionCount, (const char*)outfile);
    delete []btraits;
//*
    if(tuneTraitIndex > -1) {
      Matrix newrisks[2];
      Matrix newtrait[2];
      for (int g = 0; g < 2; g++) {
        newrisks[g].Dimension(rsCount, idCount);
        newalltraits[g].Dimension(idCount, affectionCount);
      } // End of g loop
      int index[2];
      for (int g = 0; g < 2; g++) {
        newrisks[g].Set(_NAN_);
        index[0] = index[1] = 0;
      }   // End of a loop
      for (int i = 0; i < idCount; i++) {
        int g = alltraits[tuneTraitIndex0][i] - 1;
        if (g == 0 || g == 1) {
                for (int t = 0; t < alltraitCount; t++)
                    newalltraits[g][t][index[g]] = alltraits[t][i];
                for(int r = 0; r < rsCount; r++)
                    newrisks[g][r][index[g]] = risks[r][i];
                index[g]++;
        }
      }   // End of i loop
        for (int g = 0; g < 2; g++)
            newrisks[g].Dimension(rsCount, index[g]);
        printf("Computing AUCs stratified by %s.\n", tuneTrait);
        String outfile;
        for (int g = 0; g < 2; g++) {
//            ComputeAUCs(newrisks[g], newalltraits[g], affections, auc, meanRS1, meanRS2, coreCount);
            outfile.Copy(prefix);
            outfile.Add("_");
            outfile.Add(tuneTrait);
            outfile.Add(g?"_2":"_1");
            outfile.Add(".auc");
            fp = fopen(outfile, "wb");
            if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
            WriteAUC(fp, alltraitNames, affections, rsNames, annotation, invalidPGS, auc);
            fclose(fp);
            printf("AUCs for %d risk scores and %d dichotomous traits write to %s\n",
                rsCount, affectionCount, (const char*)outfile);
        }   // End of g loop
 //
  }
}
   */


/*
    outfile.Copy(prefix);
    outfile.Add(".stat");
    fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    char buffer[0x20000];
    int pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "PGS");
    for (int a = 0; a < affectionCount; a++)
        pbuffer += sprintf(&buffer[pbuffer], "\t%s_2,1", (const char*)alltraitNames[affections[a]]);
    for (int r = 0; r < rsCount; r++) {
        pbuffer += sprintf(&buffer[pbuffer], "\n%s", (const char*)rsNames[r]);
        for (int a = 0; a < affectionCount; a++)
        pbuffer += sprintf(&buffer[pbuffer], "\t%.3lf,%.3lf", meanRS2[r][a], meanRS1[r][a]);
        if (pbuffer > 0xFFFF) {  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
        }
    }   // End of r loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
    fclose(fp);
    printf("Polygenic score averages by trait (%d scores) write to %s\n", rsCount, (const char*)outfile);
    }
    */


        /*
    String outfile;
    outfile.Copy(prefix);
    outfile.Add(".auc");
    IntArray invalidPGS(0);
    FILE *fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    WriteAUCsCV(fp, alltraitNames, affections, rsNames, annotation, invalidPGS,
      aucTrain, sdTrain, aucTest, sdTest);
    fclose(fp);
    printf("AUCs for %d PGS scores and %d dichotomous traits write to %s\n", rsCount, affectionCount, (const char*)outfile);
*/

/*
    StringArray familyID, sampleID, alltraitNames;
    Matrix alltraits;   // Step 1: Read in traits
    ReadTraits(phefile, alltraitNames, familyID, sampleID, alltraits, coreCount);
    int alltraitCount = alltraitNames.Length();
    int idCount = sampleID.Length();
    IntArray covariates(0);
    StringArray covList;
    if (covariateList != "")    // Step 2: Make covariates
        MakeCovariates(covariateList, alltraitNames, covList, covariates);
    int covariateCount = covariates.Length();
    IntArray affections(0), quantraits(0);  // Step 3: Make traits
    MakeTraits(alltraits, alltraitNames, covariates, affections, quantraits);
    int affectionCount = affections.Length();
    if(affectionCount==0) {printf("Dichotomous traits not found.\n");return;}
    StringArray rsNames;
    Matrix risks;   // Step 4: Read in risk scores
    ReadRisks(rsfile, rsNames, familyID, sampleID, risks, coreCount);
    int rsCount = rsNames.Length();
    printf("%d risk scores for %d samples loaded.\n", rsCount, idCount);
    StringArray annotation(0);
    if (dictfile[0]) {  // Step 5: Read in risk score annotations
        ReadKeys(dictfile, rsNames, annotation);
        printf("Data dictionary loaded.\n");
    }
    int annotationCount = annotation.Length();
    if (covList.Find("PC1") > -1)   // Step 6: Adjust PCs
        AdjustPC(risks, alltraits, alltraitNames, covList, coreCount);
        */




double ComputeAUC(Vector & risks, IntArray & diseaseStatus, Vector & stat)
{
   Vector temprisk[2];
   for(int a = 0; a < 2; a++)
      temprisk[a].Dimension(0);
   double minRisk = 10000, maxRisk = -10000;
   int totalCount = risks.Length();
   for(int i = 0; i < totalCount; i++){
      int aff = diseaseStatus[i];
      if(aff != 2 && aff != 1 || risks[i] == _NAN_) continue;
      temprisk[aff-1].Push(risks[i]);
      if (risks[i] > maxRisk) maxRisk = risks[i];
      if (risks[i] < minRisk) minRisk = risks[i];
   }    // End of i loop
   stat.Dimension(2); stat.Zero();
   if (maxRisk - minRisk < 0.0001) return 0.5;
   if (temprisk[0].Length() == 0 || temprisk[1].Length() == 0) return 0.5;
   stat[0] = temprisk[0].Average();
   stat[1] = temprisk[1].Average();
   const int TOTALCUT = 10000;
   int Pcount[2][TOTALCUT], Ncount[2][TOTALCUT];
   for (int t = 0; t < TOTALCUT; t++)
       Pcount[0][t] = Pcount[1][t] = Ncount[0][t] = Ncount[1][t] = 0;
   double DELTA = TOTALCUT / (maxRisk - minRisk);
   QuickIndex riskIdx;
   for (int a = 0; a < 2; a++) {
       riskIdx.Index(temprisk[a]);
       int affCount = temprisk[a].Length();
       int oldT = 0, newT = 0;
       for (int i = 0; i < affCount; i++) {
           newT = int((temprisk[a][riskIdx[i]] - minRisk) * DELTA);
           if (newT >= TOTALCUT) newT = TOTALCUT - 1;   // the maximum value maps to TOTALCUT
           int k = affCount - i;
           for (int t = oldT; t < newT; t++) {
               Pcount[a][t] = i;
               Ncount[a][t] = k;
           }  // End of t loop
           Pcount[a][newT] = i + 1;
           Ncount[a][newT] = k - 1;
           oldT = newT + 1;
       }  // End of i loop
       for (int t = oldT; t < TOTALCUT; t++) {
           Pcount[a][t] = affCount;
           Ncount[a][t] = 0;
       } // End of t loop
   }  // End of a loop
   double sensitivity[TOTALCUT], specificity[TOTALCUT];
   for (int t = 0; t < TOTALCUT; t++) {
       specificity[t] = (Ncount[0][t] + Pcount[0][t]) ? 1.0 * Ncount[0][t] / (Ncount[0][t] + Pcount[0][t]) : 1;
       sensitivity[t] = (Ncount[1][t] + Pcount[1][t]) ? 1.0 * Pcount[1][t] / (Ncount[1][t] + Pcount[1][t]) : 1;
   }    // End of t loop
   double AUC = (1 - specificity[1])*sensitivity[0];
   for (int i = 1; i < TOTALCUT - 1; i++)
       AUC += (specificity[i - 1] - specificity[i + 1])*sensitivity[i];
   AUC += (specificity[TOTALCUT-2] - specificity[TOTALCUT-1])*sensitivity[TOTALCUT-1];
   AUC *= 0.5;
   return 1-AUC;
}


  /*
  double aTrain, aTest, sTrain, sTest;
  for(int r = 0; r < rsCount; r++)
    for(int a = 0; a < affectionCount; a++){
      aTrain = aTest = sTrain = sTest = 0.0;
      for(int run = 0; run < repeat; run++){
        aTrain += AUC[r][run*affectionCount*2+a*2];
        sTrain += AUC[r][run*affectionCount*2+a*2]*AUC[r][run*affectionCount*2+a*2];
        aTest += AUC[r][run*affectionCount*2+a*2+1];
        sTest += AUC[r][run*affectionCount*2+a*2+1]*AUC[r][run*affectionCount*2+a*2+1];
      } // End of run loop
      aucTrain[r][a] = aTrain / repeat;
      aucTest[r][a] = aTest / repeat;
      sTrain = sTrain / repeat - aucTrain[r][a] * aucTrain[r][a];
      sdTrain[r][a] = sTrain < 0? 0: sqrt(sTrain);
      sTest = sTest / repeat - aucTest[r][a] * aucTest[r][a];
      sdTest[r][a] = sTest < 0? 0: sqrt(sTest);
    } // End of a loop
    */
//  Vector meanAUC(affectionCount), sdAUC(affectionCount), freqValidated(affectionCount);
//  IntArray validatedIndex(affectionCount);

void WriteAUCsCV(FILE *fp, StringArray &alltraitNames, IntArray &affections,
  StringArray &rsNames, StringArray &annotation, IntArray &invalidPGS,
  Matrix &aucTrain, Matrix &sdTrain, Matrix &aucTest, Matrix &sdTest)
{
    int affectionCount = affections.Length();
    int rsCount = rsNames.Length();
    int annotCount = annotation.Length();
    char buffer[0x20000];
    int pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "PGS");
    for (int a = 0; a < affectionCount; a++)
        pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)alltraitNames[affections[a]]);
    if (annotCount)
        pbuffer += sprintf(&buffer[pbuffer], "\tAnnotation");
    for (int r = 0; r < rsCount; r++) {
        if (invalidPGS.Find(r) > -1) continue;
        pbuffer += sprintf(&buffer[pbuffer], "\n%s", (const char*)rsNames[r]);
        for (int d = 0; d < affectionCount; d++){
            pbuffer += sprintf(&buffer[pbuffer], "\t%.5lf", aucTrain[r][d]);
            pbuffer += sprintf(&buffer[pbuffer], "\t%.5lf", sdTrain[r][d]);
            pbuffer += sprintf(&buffer[pbuffer], "\t%.5lf", aucTest[r][d]);
            pbuffer += sprintf(&buffer[pbuffer], "\t%.5lf", sdTest[r][d]);
        } // End of d loop
        if (annotCount)
            pbuffer += sprintf(&buffer[pbuffer], "\t%s", (const char*)annotation[r]);
        if (pbuffer > 0xFFFF) {  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
        }
    }   // End of r loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
}


