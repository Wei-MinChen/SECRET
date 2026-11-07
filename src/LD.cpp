////////////////////////////////////////////////////////////// ////////
// LD.cpp
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
// Oct 13, 2026

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "Constant.h"
#include "Error.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

inline unsigned char popcount(unsigned long long int word)
{
   word = word - ((word>>1)&0x5555555555555555);
   word = (word&0x3333333333333333) + ((word>>2)&0x3333333333333333);
   word = (word+(word>>4)) & 0x0F0F0F0F0F0F0F0F;
   word = (word+(word>>8)) & 0x00FF00FF00FF00FF;
   word = (word+(word>>16)) & 0x0000FFFF0000FFFF;
   return((word+(word>>32)) & 0xFF);
}

void ScreenLD1(Matrix &risks, double ldT, IntArray pairList[], int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
  IntArray count1(rsCount);
  int iwordCount = (idCount-1)/64+1;
  unsigned long long int **brisks = new unsigned long long int*[rsCount];
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int r = 0; r < rsCount; r++){
    double mean = risks[r].Average();
    brisks[r] = new unsigned long long int[iwordCount];
    for(int i = 0; i < idCount; i+=64){
      int iMax = i < idCount-64? i+64: idCount;
      int w = (i >> 6);
      unsigned long long int word = 0;
      unsigned long long int word0 = 0x1;
      for(int ii = i; ii < iMax; ii++, word0 <<= 1)
        if(risks[r][ii] > mean)
          word |= word0;
      brisks[r][w] = word;
    } // End of i loop
    int count = 0;
    for(int w = 0; w < iwordCount; w++)
      count += popcount(brisks[r][w]);
    count1[r] = count;
  } // End of OMP r loop
  const int RBLOCKSIZE=64;
  const int RSHIFTSIZE=6;
  const int IBLOCKSIZE=128; // 2^(6+7+3+1) = 128KB, N=2^(6+7)=8K
  int wordCount0 = (rsCount-1)/RBLOCKSIZE+1;
  int wordCount = (rsCount-1)/64+1;
  IntArray loopIndex(0);
  for(int w1 = 0; w1 < wordCount0; w1++)
    for(int w2 = (w1>>(6-RSHIFTSIZE)); w2 < wordCount; w2++)
      loopIndex.Push((w1<<16) | w2);
  int loopCount = loopIndex.Length();
  for(int i = 0; i < coreCount; i++)
    pairList[i].Dimension(0);
  unsigned long long int word, wordN;
  double means[64], vars[64];
  int count11[RBLOCKSIZE][64]; // 2^(6+6+4)=64KB
  int thread = 0;
#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount) private(thread, word,wordN,means,vars,count11)
{
  thread = omp_get_thread_num();
  #pragma omp for
#endif
  for(int l = 0; l < loopCount; l++){
    int r = (loopIndex[l]>>(16-RSHIFTSIZE));
    int s = ((loopIndex[l]&0xFFFF)<<6);
    int rMax = r < rsCount-RBLOCKSIZE? r+RBLOCKSIZE: rsCount;
    int sMax = s < rsCount-64? s+64: rsCount;
    memset(&count11[0][0], 0, RBLOCKSIZE*64*sizeof(int));
    for(int w = 0; w < iwordCount; w+=IBLOCKSIZE){
      int wMax = w < iwordCount-IBLOCKSIZE? w+IBLOCKSIZE: iwordCount;
      for(int rr = r; rr < rMax; rr++)
        for(int ss = s; ss < sMax; ss++){ // CACHE BLOCK: 128KB + 64KB
          wordN = 0;
          for (int ww = w; ww < wMax; ww++) {
            word = brisks[rr][ww] & brisks[ss][ww];
            word = word - ((word >> 1) & 0x5555555555555555);
            word = (word & 0x3333333333333333) + ((word >> 2) & 0x3333333333333333);
            word = (word + (word >> 4)) & 0x0F0F0F0F0F0F0F0F;
            wordN += (word + (word >> 8)) & 0x00FF00FF00FF00FF;
          } // End of ww loop
          wordN = (wordN + (wordN >> 16)) & 0x0000FFFF0000FFFF;
          count11[rr-r][ss-s] += ((wordN + (wordN >> 32)) & 0xFFFFFFFF);
        } // End of ss loop CACHE BLOCK
    } // End of w loop
    for(int ss = s; ss < sMax; ss++){
      means[ss-s] = count1[ss]*1.0/idCount;
      vars[ss-s] = means[ss-s] * (1-means[ss-s]);
    } // End of ss loop
    for(int rr = 0; rr < rMax-r; rr++){
      double meanr = count1[r+rr]*1.0/idCount;
      double varr = meanr * (1-meanr);
      for(int ss = 0; ss < sMax-s; ss++){
        double cor=(count11[rr][ss]*1.0/idCount - meanr * means[ss])/sqrt(varr*vars[ss]);
        if(cor*cor>ldT*0.25 && r+rr < s+ss)
          pairList[thread].Push(((r+rr)<<16)|(s+ss));
      } // End of ss loop
    } // End of rr loop
  } // End of OMP l loop
#ifdef _OPENMP
}
#endif
}

void ScreenLD2(Matrix &risks, double ldT, IntArray pairList[], Vector pairLD[], Vector &rsSD, int coreCount)
{
  ScreenLD1(risks, ldT, pairList, coreCount);
  int rsCount = risks.rows;
  int idCount = risks.cols;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int r = 0; r < rsCount; r++){
    double mean = risks[r].Average();
    //double sd = sqrt(risks[r].Var());
    rsSD[r] = sqrt(risks[r].Var());
    if(rsSD[r] < 1E-100) error("Cannot deal with PGS with no variation");
    for(int i = 0; i < idCount; i++)
      risks[r][i] = (risks[r][i] - mean) / rsSD[r];
  } // End of OMP r loop
  const int BLOCKSIZE=16;
  const int IBLOCKSIZE=512; // 2^(4+9+3+1)=128KB
  double sums[BLOCKSIZE];
  double thres = sqrt(ldT) * idCount;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(sums)
#endif
  for(int c = 0; c < coreCount; c++){
    int pairCount = pairList[c].Length();
    pairLD[c].Dimension(pairCount);
    int newpairCount = 0;
    for(int p = 0; p < pairCount; p+=BLOCKSIZE){
      int pMax = p < pairCount - BLOCKSIZE? p+BLOCKSIZE: pairCount;
      memset(&sums[0], 0.0, BLOCKSIZE*sizeof(double));
      for(int i = 0; i < idCount; i+=IBLOCKSIZE){
        int iMax = i < idCount-IBLOCKSIZE? i+IBLOCKSIZE: idCount;
        for(int pp = p; pp < pMax; pp++){ // CACHE block
          int r = (pairList[c][pp]>>16);
          int s = (pairList[c][pp]&0xFFFF);
          for(int ii = i; ii < iMax; ii++)
            sums[pp-p] += risks[r][ii] * risks[s][ii];
        } // End of pp loop CACHE block
      } // End of i loop
      for(int pp = p; pp < pMax; pp++)
        if((sums[pp-p] > thres || sums[pp-p] < -thres)){
          pairList[c][newpairCount] = pairList[c][pp];
          pairLD[c][newpairCount++] = sums[pp-p] / idCount;
        } // End of pp loop
    } // End of p loop
    pairList[c].Dimension(newpairCount);
    pairLD[c].Dimension(newpairCount);
  } // End of OMP c loop
}

void MakeBinaryLD(int rsCount, IntArray pairList[], unsigned long long int **LD, int coreCount)
{
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int r = 0; r < rsCount; r++)
    LD[r][r>>6] |= ((unsigned long long int)0x1<<(r&0x3F));
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int c = 0; c < coreCount; c++){
    int pairCount = pairList[c].Length();
    for(int p = 0; p < pairCount; p++){
      int r = (pairList[c][p]>>16);
      int s = (pairList[c][p]&0xFFFF);
      LD[r][s>>6] |= ((unsigned long long int)0x1<<(s&0x3F));
    } // End of p loop
  } // End of OMP c loop
  unsigned long long int words[64];
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(words)
#endif
  for(int s = 0; s < rsCount; s += 64){
    int sMax = (s < rsCount - 64? s+64: rsCount);
    int sCol = (s>>6);
    for(int r = 0; r <= s; r += 64){
      int rMax = (r < rsCount - 64? r+64: rsCount);
      int rCol = (r>>6);
      memset(&words[0], 0, 64*sizeof(unsigned long long int));
      for(int rr = r; rr < rMax; rr++){
        unsigned long long int word = LD[rr][sCol];
        if(word){
          unsigned long long int word0 = ((unsigned long long int)0x1 << (rr&0x3F));
          for(int b = 0; b < 64; b++, word >>= 1)
            if( (word & 0x1) && (rr < (s|b)))
              words[b] |= word0;
        }
      } // End of rr loop
      for(int ss = s; ss < sMax; ss++)
        LD[ss][rCol] |= words[ss-s];
    } // End of r loop
  } // End of OMP s loop
}

void MakeBinaryLD(Matrix &risks, double ldT, unsigned long long int **LD, Vector &rsSD, int coreCount)
{
  IntArray *pairList = new IntArray[coreCount];
  Vector *pairLD = new Vector[coreCount];
  ScreenLD1(risks, ldT, pairList, coreCount);
  ScreenLD2(risks, ldT, pairList, pairLD, rsSD, coreCount);
  MakeBinaryLD(risks.rows, pairList, LD, coreCount);
  delete []pairList;
  delete []pairLD;
}

void MakePrunedOrder(unsigned long long int **LD, IntArray &order)
{
  int rsCount = order.Length();
  IntArray flagged(rsCount); flagged.Zero();
  IntArray prunedList(rsCount);
  int prunedCount = 0;
  for (int r = 0; r < rsCount; r++) {
    for (; r < rsCount && flagged[r]; r++);
    if (r == rsCount) break;
    prunedList[prunedCount++] = r;
    int orderR = order[r];
    for (int s = r+1; s < rsCount; s++){
      int orderS = order[s];
      if ((flagged[s] == 0) && ((LD[orderR][orderS>>6] >> (orderS&0x3F))&0x1))
        flagged[s] = 1;
    } // End of s loop
  } // End of r loop
  for (int r = 0; r < prunedCount; r++)
    order[r] = prunedList[r];
  //    order[r] = order[prunedList[r]];
  order.Dimension(prunedCount);
}

int popcounts(unsigned long long int *words, int wordCount)
{
  unsigned long long int wSum = 0;
  int count = 0;
  for(int w = 0; w < wordCount; w += 31){
    int wMax = w < wordCount - 31? w+31: wordCount;
    wSum = 0;
    for(int ww = w; ww < wMax; ww++){
      unsigned long long int word = words[ww];
      for(int b = 0; b < 8; b++, word >>= 1)// intensive
        wSum += (word & 0x0101010101010101);
    } // End of ww loop
    for(int b = 0; b < 8; b++, wSum >>= 8)
      count += (wSum & 0xFF);
  } // End of w loop
  return count;
}

void LDSummary(unsigned long long int **LD, int rsCount, int wordCount, int coreCount)
{
  IntArray correlatedCount(rsCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for(int r = 0; r < rsCount; r++)
      correlatedCount[r] = popcounts(LD[r], wordCount) - 1;
  double average = 0.0;
  for(int r = 0; r < rsCount; r++)
    average += correlatedCount[r];
  average /= rsCount;
  QuickIndex idx;
  idx.Index(correlatedCount);
  printf("  Distribution of connections for each PGS:\n");
  printf("  %5s %5s %5s %5s %5s %5s\n", "MIN", "25\%", "50\%", "75\%", "MAX", "AVE");
  printf("  %5d %5d %5d %5d %5d %5.1lf\n", correlatedCount[idx[0]],
    correlatedCount[idx[rsCount/4]], correlatedCount[idx[rsCount/2]],
    correlatedCount[idx[rsCount*3/4]], correlatedCount[idx[rsCount-1]],
    average);
}

void LDSummaryByAffection(unsigned long long int ***LD, int rsCount, int *wordCount, StringArray &PGSnames, int coreCount)
{
  QuickIndex idx[2];
  IntArray correlatedCount[2];
  double average[2];
  for(int a = 0; a < 2; a++){
    correlatedCount[a].Dimension(rsCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for(int r = 0; r < rsCount; r++)
      correlatedCount[a][r] = popcounts(LD[a][r], wordCount[a]) - 1;
    average[0] = average[1] = 0.0;
    for(int r = 0; r < rsCount; r++)
      average[a] += correlatedCount[a][r];
    average[a] /= rsCount;
    idx[a].Index(correlatedCount[a]);
  } // End of a loop
  printf("  Distribution of connections for each PGS:\n");
  printf("  %10s %5s %5s %5s %5s %5s %5s\n", "Status", "MIN", "25\%", "50\%", "75\%", "MAX", "AVE");
  printf("  %10s %5d %5d %5d %5d %5d %5.1lf\n", "Control", correlatedCount[0][idx[0][0]],
    correlatedCount[0][idx[0][rsCount/4]], correlatedCount[0][idx[0][rsCount/2]],
    correlatedCount[0][idx[0][rsCount*3/4]], correlatedCount[0][idx[0][rsCount-1]],
    average[0]);
  printf("  %10s %5d %5d %5d %5d %5d %5.1lf\n", "Case", correlatedCount[1][idx[1][0]],
    correlatedCount[1][idx[1][rsCount/4]], correlatedCount[1][idx[1][rsCount/2]],
    correlatedCount[1][idx[1][rsCount*3/4]], correlatedCount[1][idx[1][rsCount-1]],
    average[1]);
  printf("  Top 10 hub PGS (with the most connections):\n");
  printf("  %5s %15s %15s\n", "Rank", "Control", "Case");
  for(int i = 0; i < 10; i++)
    printf("  %5d %15s %15s\n", i+1,
      (const char*)PGSnames[correlatedCount[0][idx[0][rsCount-i-1]]],
      (const char*)PGSnames[correlatedCount[1][idx[1][rsCount-i-1]]]);
}





/*
void MakeBinaryLD(Matrix &risks, double ldT, unsigned long long int **LD, int coreCount)
{
  int rsCount = risks.rows;
  int idCount = risks.cols;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int r = 0; r < rsCount; r++){
    double mean = risks[r].Average();
    double sd = sqrt(risks[r].Var());
    if(sd < 1E-100) error("Cannot deal with PGS with no variation");
    for(int i = 0; i < idCount; i++)
      risks[r][i] = (risks[r][i] - mean) / sd;
  } // End of OMP r loop
  const int RBLOCKSIZE=16;
  const int RSHIFTSIZE=4;
  const int IBLOCKSIZE=128;
  double dotprod[RBLOCKSIZE][64];
  int wordCount0 = (rsCount-1)/RBLOCKSIZE+1;
  int wordCount = (rsCount-1)/64+1;
  int memSize = wordCount*sizeof(unsigned long long int);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int r = 0; r < rsCount; r++)
    memset(&LD[r][0], 0, memSize);
  int screenCount = idCount > 256? 256: idCount;
  IntArray loopIndex(0);
  for(int w1 = 0; w1 < wordCount0; w1++)
    for(int w2 = (w1>>(6-RSHIFTSIZE)); w2 < wordCount; w2++)
      loopIndex.Push((w1<<16) | w2);
  int loopCount = loopIndex.Length();
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(dotprod)
#endif
  for(int l = 0; l < loopCount; l++){
    int r = (loopIndex[l]>>(16-RSHIFTSIZE));
    int s = ((loopIndex[l]&0xFFFF)<<6);
    int rMax = r < rsCount-RBLOCKSIZE? RBLOCKSIZE: rsCount-r;
    int sMax = s < rsCount-64? 64: rsCount-s;
    memset(&dotprod[0][0], 0, RBLOCKSIZE*64*sizeof(double));
    for(int i = 0; i < screenCount; i += IBLOCKSIZE){ // Stage 1: every pair & subset id
      int iMax = ((i < screenCount-IBLOCKSIZE)? i+IBLOCKSIZE: screenCount);
      for(int rr = 0; rr < rMax; rr++)
        for(int ss = 0; ss < sMax; ss++){
          double sum = 0.0;
          for(int ii = i; ii < iMax; ii++)
            sum += risks[r+rr][ii] * risks[s+ss][ii];
          dotprod[rr][ss] += sum;
        }  // End of ss loop
    }  // End of i loop
    double thres0 = sqrt(ldT*0.1) * screenCount;
    double thres = sqrt(ldT) * idCount;
    for(int rr = 0; rr < rMax; rr++){ // Stage 2: subset pair & every id
      unsigned long long int word = 0;
      unsigned long long int word1 = 0x1;
      for(int ss = 0; ss < sMax; ss++, word1 <<= 1)
        if(dotprod[rr][ss] > thres0 || dotprod[rr][ss] < -thres0){ // for r2 > T * T
          double sum = dotprod[rr][ss];
          for(int ii = screenCount; ii < idCount && sum >= -thres && sum <= thres; ii++)// selective & sparse
            sum += risks[r+rr][ii] * risks[s+ss][ii];
          if(sum > thres || sum < -thres)
            word |= word1;
        } // End of ss loop
      LD[r+rr][s>>6] = word;
    } // End of rr loop
  } // End of OMP l loop
  unsigned long long int words[64];
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(words)
#endif
  for(int s = 0; s < rsCount; s += 64){
    int sMax = (s < rsCount - 64? s+64: rsCount);
    int sCol = (s>>6);
    for(int r = 0; r < s; r += 64){
      int rMax = (r < rsCount - 64? r+64: rsCount);
      int rCol = (r>>6);
      for(int b = 0; b < 64; b++)
        words[b] = 0;
      for(int rr = r; rr < rMax; rr++){
        unsigned long long int word = LD[rr][sCol];
        if(word){
          unsigned long long int word0 = ((unsigned long long int)0x1 << (rr&0x3F));
          for(int b = 0; b < 64; b++, word >>= 1)
            if(word & 0x1)
              words[b] |= word0;
        }
      } // End of rr loop
      for(int ss = s; ss < sMax; ss++)
        LD[ss][rCol] = words[ss-s];
    } // End of r loop
    int r = s;
    int rMax = (r < rsCount - 64? r+64: rsCount);
    int rCol = (r>>6);
    for(int b = 0; b < 64; b++)
      words[b] = 0;
    for(int rr = r; rr < rMax; rr++){
      unsigned long long int word = LD[rr][sCol];
      if(word){
        unsigned long long int word0 = ((unsigned long long int)0x1 << (rr&0x3F));
        for(int b = 0; b < 64; b++, word >>= 1)
          if( (word & 0x1) && (rr < (s|b)))
            words[b] |= word0;
      }
    } // End of rr loop
    for(int ss = s; ss < sMax; ss++)
      LD[ss][rCol] |= words[ss-s];
  } // End of s loop
}
*/

/*
void PruneSmallLD(Matrix &LD, double ldT, IntArray &prunedList)
{
    int rsCount = LD.rows;
    prunedList.Dimension(rsCount);
    int prunedCount = 0;
    IntArray flagged(rsCount); flagged.Zero();
    for (int r = 0; r < rsCount; r++) {
        for (; r < rsCount && flagged[r]; r++);
        if (r == rsCount) break;
        prunedList[prunedCount++] = r;
        for (int s = r+1; s < rsCount; s++)
            if ((flagged[s] == 0) && (LD[r][s] >= ldT))
                flagged[s] = 1;
    }   // End of r loop
    prunedList.Dimension(prunedCount);
}
*/
/*
void MakeSmallLD(Matrix &risks, IntArray &order, Matrix &LD)
{
  int rsCount = order.Length();
  LD.Dimension(rsCount, rsCount);
  Vector mean(rsCount), var(rsCount);
  int idCount = risks.cols;
  for (int r = 0; r < rsCount; r++) {
        mean[r] = risks[order[r]].Average();
        var[r] = risks[order[r]].Var();
  }   // End of r loop
  for(int r = 0; r < rsCount; r++)
    for (int s = r + 1; s < rsCount; s++){
      double sum = 0.0;
      for (int i = 0; i < idCount; i++)
        sum += risks[order[r]][i] * risks[order[s]][i];
      sum = sum / idCount - mean[r] * mean[s];
      LD[r][s] = sum * sum / (var[r] * var[s]);
    }   // End of s loop
} */
/*
void MakePrunedSmallOrder(Matrix &risks, IntArray &order, double ldT, int &updatedCount)
{
  int rsCount = order.Length();
  Matrix LD(rsCount, rsCount);
  Vector mean(rsCount), var(rsCount);
  for (int r = 0; r < rsCount; r++) {
    mean[r] = risks[order[r]].Average();
    var[r] = risks[order[r]].Var();
  }   // End of r loop
  int idCount = risks.cols;
  for(int r = 0; r < rsCount; r++)
    for (int s = r + 1; s < rsCount; s++){
      double sum = 0.0;
      for (int i = 0; i < idCount; i++)
        sum += risks[order[r]][i] * risks[order[s]][i];
      sum = sum / idCount - mean[r] * mean[s];
      LD[r][s] = sum * sum / (var[r] * var[s]);
    } // End of s loop
  IntArray prunedList(rsCount);
  int prunedCount = 0;
  IntArray flagged(rsCount); flagged.Zero();
  int toprsCount = (rsCount >> 4);
  for (int r = 0; r < rsCount && prunedCount < toprsCount; r++) {
    for (; r < rsCount && flagged[r]; r++);
      if (r == rsCount) break;
      prunedList[prunedCount++] = r;
      for (int s = r+1; s < rsCount; s++)
        if ((flagged[s] == 0) && (LD[r][s] >= ldT))
          flagged[s] = 1;
  }   // End of r loop
  for (int r = 0; r < prunedCount; r++)
    order[r] = order[prunedList[r]];
  updatedCount = prunedCount;
}
*/
/*
  int rsCount = order.Length();
  int toprsCount = (rsCount >> 4);
  IntArray flagged(rsCount); flagged.Zero();
  IntArray prunedList(rsCount);
  int prunedCount = 0;
  for (int r = 0; r < rsCount && prunedCount < toprsCount; r++) {
    for (; r < rsCount && flagged[r]; r++);
    if (r == rsCount) break;
    prunedList[prunedCount++] = r;
    int orderR = order[r];
    for (int s = r+1; s < rsCount; s++){
      int orderS = order[s];
      if ((flagged[s] == 0) && ((LD[orderR][orderS>>6] >> (orderS&0x3F))&0x1))
        flagged[s] = 1;
    } // End of s loop
  } // End of r loop
  for (int r = 0; r < prunedCount; r++)
    order[r] = order[prunedList[r]];
  updatedCount = prunedCount;
*/


// =====================================================================
// Correlations among the top PGSs, overall and by affection status
// =====================================================================

// Pearson correlation among the PGSs rsIndex (rows of risks) over the samples
// sampleIndex (columns of risks).  Each pair of PGSs uses the samples where both
// are observed.  R and N are p x p: N is the number of samples used by the pair,
// R is _NAN_ when that number is below 3 or one of the PGSs does not vary.
void PGSCorrelation(Matrix &risks, IntArray &rsIndex, IntArray &sampleIndex, Matrix &R, Matrix &N, int coreCount)
{
  int p = rsIndex.Length();
  int n = sampleIndex.Length();
  R.Dimension(p, p);
  N.Dimension(p, p);
  R.Set(_NAN_);
  N.Zero();
  if (p == 0 || n == 0) return;
  Matrix X(p, n);
  for(int a = 0; a < p; a++){
    double sum = 0.0;
    int cnt = 0;
    for(int k = 0; k < n; k++){
      double v = risks[rsIndex[a]][sampleIndex[k]];
      X[a][k] = v;
      if(v != _NAN_) { sum += v; cnt++; }
    } // End of k loop
    double mean = cnt ? sum / cnt : 0.0;   // centre first: the sums below stay well conditioned
    for(int k = 0; k < n; k++)
      if(X[a][k] != _NAN_) X[a][k] -= mean;
  } // End of a loop
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(coreCount)
#endif
  for(int a = 0; a < p; a++)
    for(int b = a; b < p; b++){
      double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
      int cnt = 0;
      for(int k = 0; k < n; k++){
        double x = X[a][k], y = X[b][k];
        if(x == _NAN_ || y == _NAN_) continue;
        cnt++;
        sx += x; sy += y;
        sxx += x * x; syy += y * y; sxy += x * y;
      } // End of k loop
      N[a][b] = N[b][a] = cnt;
      if(cnt < 3) continue;
      double vx = sxx - sx * sx / cnt;
      double vy = syy - sy * sy / cnt;
      if(vx < 1E-100 || vy < 1E-100) continue;
      double r = (sxy - sx * sy / cnt) / sqrt(vx * vy);
      if(r > 1.0) r = 1.0;
      if(r < -1.0) r = -1.0;
      R[a][b] = R[b][a] = r;
    } // End of b loop
  for(int a = 0; a < p; a++)
    if(N[a][a] >= 3 && R[a][a] != _NAN_) R[a][a] = 1.0;
}

// Gauss-Jordan inverse of a p x p matrix (row-major); false when singular
static bool InvertMatrix(int p, const double *A, double *inv)
{
  double *M = new double[p * p];
  for(int i = 0; i < p * p; i++) M[i] = A[i];
  for(int i = 0; i < p; i++)
    for(int j = 0; j < p; j++)
      inv[i * p + j] = (i == j) ? 1.0 : 0.0;
  bool ok = true;
  for(int c = 0; c < p && ok; c++){
    int piv = c;
    for(int r = c + 1; r < p; r++)
      if(fabs(M[r * p + c]) > fabs(M[piv * p + c])) piv = r;
    if(fabs(M[piv * p + c]) < 1E-12) { ok = false; break; }
    if(piv != c)
      for(int j = 0; j < p; j++){
        double t = M[c * p + j]; M[c * p + j] = M[piv * p + j]; M[piv * p + j] = t;
        t = inv[c * p + j]; inv[c * p + j] = inv[piv * p + j]; inv[piv * p + j] = t;
      } // End of j loop
    double d = M[c * p + c];
    for(int j = 0; j < p; j++){
      M[c * p + j] /= d;
      inv[c * p + j] /= d;
    } // End of j loop
    for(int r = 0; r < p; r++){
      if(r == c) continue;
      double f = M[r * p + c];
      if(f == 0.0) continue;
      for(int j = 0; j < p; j++){
        M[r * p + j] -= f * M[c * p + j];
        inv[r * p + j] -= f * inv[c * p + j];
      } // End of j loop
    } // End of r loop
  } // End of c loop
  delete []M;
  return ok;
}

// Tests of the difference between the PGS correlations of two independent groups
// (here the controls, R1, and the cases, R2; N1, N2 are the pair sample sizes).
//
//  - Pair by pair: the Fisher z test,
//        z = (atanh(r2) - atanh(r1)) / sqrt(1/(n2-3) + 1/(n1-3)),
//    with Bonferroni and Benjamini-Hochberg (FDR) adjustments over the p(p-1)/2 pairs
//    (Zd, Pd, Padj = BH, Pbon; symmetric, _NAN_ on the diagonal or when untestable).
//  - Whole matrix: Jennrich's (1970, JASA 65:904) chi-square test of H0: R1 = R2,
//        chi2 = tr(Z^2)/2 - diag(Z)' S^-1 diag(Z),
//        Z = sqrt(n1 n2/(n1+n2)) Rbar^-1 (R1 - R2),  S = I + Rbar o Rbar^-1,
//    where Rbar is the sample-size weighted mean of R1 and R2 and o is the
//    elementwise product; chi2 has p(p-1)/2 degrees of freedom.  chi2 is _NAN_
//    when a correlation is missing or Rbar is singular.
void PGSCorrelationTest(Matrix &R1, Matrix &N1, Matrix &R2, Matrix &N2,
  Matrix &Zd, Matrix &Pd, Matrix &Padj, Matrix &Pbon, double &chi2, int &df, double &pJennrich)
{
  int p = R1.rows;
  Zd.Dimension(p, p); Pd.Dimension(p, p); Padj.Dimension(p, p); Pbon.Dimension(p, p);
  Zd.Set(_NAN_); Pd.Set(_NAN_); Padj.Set(_NAN_); Pbon.Set(_NAN_);
  chi2 = _NAN_;
  pJennrich = _NAN_;
  df = p * (p - 1) / 2;
  IntArray pa(0), pb(0);
  Vector pv(0);
  for(int a = 0; a < p; a++)
    for(int b = a + 1; b < p; b++){
      double r1 = R1[a][b], r2 = R2[a][b];
      double n1 = N1[a][b], n2 = N2[a][b];
      if(r1 == _NAN_ || r2 == _NAN_ || n1 <= 3 || n2 <= 3) continue;
      if(r1 > 0.999999) r1 = 0.999999;
      if(r1 < -0.999999) r1 = -0.999999;
      if(r2 > 0.999999) r2 = 0.999999;
      if(r2 < -0.999999) r2 = -0.999999;
      double z1 = 0.5 * log((1.0 + r1) / (1.0 - r1));
      double z2 = 0.5 * log((1.0 + r2) / (1.0 - r2));
      double z = (z2 - z1) / sqrt(1.0 / (n2 - 3.0) + 1.0 / (n1 - 3.0));
      double pval = 2.0 * ndist(fabs(z));
      if(pval > 1.0) pval = 1.0;
      Zd[a][b] = Zd[b][a] = z;
      Pd[a][b] = Pd[b][a] = pval;
      pa.Push(a); pb.Push(b); pv.Push(pval);
    } // End of b loop
  int m = pv.Length();
  if(m > 0){
    QuickIndex idx;
    idx.Index(pv);   // ascending
    double running = 1.0;
    for(int rank = m; rank >= 1; rank--){
      int j = idx[rank - 1];
      double adj = pv[j] * m / rank;
      if(adj < running) running = adj;
      Padj[pa[j]][pb[j]] = Padj[pb[j]][pa[j]] = running;
      double bon = pv[j] * m;
      Pbon[pa[j]][pb[j]] = Pbon[pb[j]][pa[j]] = bon < 1.0 ? bon : 1.0;
    } // End of rank loop
  }
  // Jennrich's test of the whole matrix
  if(p < 2) return;
  double n1 = 1E300, n2 = 1E300;
  for(int a = 0; a < p; a++)
    for(int b = 0; b < p; b++){
      if(R1[a][b] == _NAN_ || R2[a][b] == _NAN_) return;
      if(a != b){
        if(N1[a][b] < n1) n1 = N1[a][b];
        if(N2[a][b] < n2) n2 = N2[a][b];
      }
    } // End of b loop
  if(n1 < 10 || n2 < 10) return;
  double *R = new double[p * p], *Ri = new double[p * p], *Z = new double[p * p];
  double *S = new double[p * p], *Si = new double[p * p], *d = new double[p];
  for(int a = 0; a < p; a++)
    for(int b = 0; b < p; b++)
      R[a * p + b] = (n1 * R1[a][b] + n2 * R2[a][b]) / (n1 + n2);
  if(InvertMatrix(p, R, Ri)){
    double c = sqrt(n1 * n2 / (n1 + n2));
    for(int a = 0; a < p; a++)
      for(int b = 0; b < p; b++){
        double sum = 0.0;
        for(int k = 0; k < p; k++)
          sum += Ri[a * p + k] * (R1[k][b] - R2[k][b]);
        Z[a * p + b] = c * sum;
        S[a * p + b] = (a == b ? 1.0 : 0.0) + R[a * p + b] * Ri[a * p + b];
      } // End of b loop
    if(InvertMatrix(p, S, Si)){
      double trace = 0.0;
      for(int a = 0; a < p; a++){
        d[a] = Z[a * p + a];
        for(int b = 0; b < p; b++)
          trace += Z[a * p + b] * Z[b * p + a];
      } // End of a loop
      double quad = 0.0;
      for(int a = 0; a < p; a++)
        for(int b = 0; b < p; b++)
          quad += d[a] * Si[a * p + b] * d[b];
      chi2 = 0.5 * trace - quad;
      if(chi2 < 0.0) chi2 = 0.0;
      pJennrich = chidist(chi2, (double) df);
    }
  }
  delete []R; delete []Ri; delete []Z; delete []S; delete []Si; delete []d;
}

static void WritePGSMatrix(const char *outfile, Matrix &M, StringArray &names)
{
  FILE *fp = fopen(outfile, "wb");
  if(fp == NULL) error("Cannot open %s to write.", outfile);
  int p = M.rows;
  fprintf(fp, "PGS");
  for(int b = 0; b < p; b++)
    fprintf(fp, "\t%s", (const char*)names[b]);
  for(int a = 0; a < p; a++){
    fprintf(fp, "\n%s", (const char*)names[a]);
    for(int b = 0; b < p; b++)
      if(M[a][b] == _NAN_) fprintf(fp, "\tNA");
      else fprintf(fp, "\t%.4lf", M[a][b]);
  } // End of a loop
  fprintf(fp, "\n");
  fclose(fp);
}

// ---------------------------------------------------------------------
// Clusters among the top PGSs and a recommended --pruner2
// ---------------------------------------------------------------------
// WTSUM prunes the top PGSs in the order of their rank: a PGS is dropped when its
// squared correlation (R2) with a PGS that was kept before it exceeds --pruner2.  The
// PGSs that are highly correlated with each other form the clusters of the dendrogram
// drawn on the PGS correlation plot, so a good --pruner2 is one at which this pruning
// keeps about one PGS per cluster.  It is found in three steps:
//   1. average-linkage clustering on the distance 1 - |r| (the dendrogram of the plot);
//   2. the number of clusters k is the one with the largest mean silhouette width; a
//      mean silhouette below 0.25 means that there is no clear cluster structure
//      (Kaufman and Rousseeuw), and no recommendation is made;
//   3. --pruner2 is the R2 cutoff, on a grid of 0.01, at which the pruning above gives
//      the partition closest to the k clusters (largest adjusted Rand index; a tie goes
//      to the cutoff nearest to the R2 of the dendrogram height half-way between the
//      merges that leave k+1 and k clusters).
// Only the first 200 PGSs (the best ranked) are used, and at least 4 are needed.

double gPrunerR2Rec = _NAN_;   // recommended --pruner2 of the last WritePGSCorrelations call that asked for one (_NAN_: none)
int gPGSClusterCount = 0;      // its number of clusters (0: none)

static double PGSMeanSilhouette(const std::vector<double> &D, int n, const std::vector<int> &label, int k)
{
  std::vector<int> size(k, 0);
  for (int i = 0; i < n; i++) size[label[i]]++;
  std::vector<double> sum(k);
  double total = 0.0;
  for (int i = 0; i < n; i++) {
    if (size[label[i]] < 2) continue;   // a cluster of one PGS has silhouette 0
    std::fill(sum.begin(), sum.end(), 0.0);
    for (int j = 0; j < n; j++)
      if (j != i) sum[label[j]] += D[(size_t)i * n + j];
    double a = sum[label[i]] / (size[label[i]] - 1);
    double b = 1e300;
    for (int c = 0; c < k; c++)
      if (c != label[i] && size[c] > 0 && sum[c] / size[c] < b) b = sum[c] / size[c];
    double mx = a > b ? a : b;
    if (mx > 0.0) total += (b - a) / mx;
  }   // End of i loop
  return total / n;
}

// Average linkage (UPGMA).  Merge s joins the nodes mergeA[s] and mergeB[s] (leaves are
// 0..n-1, merge s is node n+s) at distance height[s]; the heights never decrease.
static void PGSAverageLinkage(const std::vector<double> &D, int n, std::vector<int> &mergeA,
  std::vector<int> &mergeB, std::vector<double> &height)
{
  std::vector<double> d(D);
  std::vector<int> size(n, 1), node(n);
  std::vector<char> active(n, 1);
  for (int i = 0; i < n; i++) node[i] = i;
  mergeA.assign(n - 1, 0); mergeB.assign(n - 1, 0); height.assign(n - 1, 0.0);
  for (int s = 0; s < n - 1; s++) {
    int bi = -1, bj = -1;
    double best = 1e300;
    for (int i = 0; i < n; i++) {
      if (!active[i]) continue;
      for (int j = i + 1; j < n; j++)
        if (active[j] && d[(size_t)i * n + j] < best) { best = d[(size_t)i * n + j]; bi = i; bj = j; }
    }   // End of i loop
    mergeA[s] = node[bi]; mergeB[s] = node[bj]; height[s] = best;
    for (int k = 0; k < n; k++) {
      if (!active[k] || k == bi || k == bj) continue;
      double nd = (size[bi] * d[(size_t)bi * n + k] + size[bj] * d[(size_t)bj * n + k]) / (size[bi] + size[bj]);
      d[(size_t)bi * n + k] = d[(size_t)k * n + bi] = nd;
    }   // End of k loop
    size[bi] += size[bj];
    active[bj] = 0;
    node[bi] = n + s;
  }   // End of s loop
}

// the root of x in the union-find forest `parent`, with path halving (a function, not a lambda:
// C++ Builder 6 has no lambdas)
static int FindRoot(std::vector<int> &parent, int x)
{
  while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
  return x;
}

// Cluster of every PGS when the dendrogram is cut into k clusters (labels 0..k-1)
static void PGSClusterLabels(int n, const std::vector<int> &mergeA, const std::vector<int> &mergeB,
  int k, std::vector<int> &label)
{
  std::vector<int> parent(2 * n - 1);
  for (int i = 0; i < 2 * n - 1; i++) parent[i] = i;
  for (int s = 0; s < n - k; s++) {
    parent[FindRoot(parent, mergeA[s])] = n + s;
    parent[FindRoot(parent, mergeB[s])] = n + s;
  }   // End of s loop
  label.assign(n, -1);
  std::vector<int> rootLabel(2 * n - 1, -1);
  int count = 0;
  for (int i = 0; i < n; i++) {
    int r = FindRoot(parent, i);
    if (rootLabel[r] < 0) rootLabel[r] = count++;
    label[i] = rootLabel[r];
  }   // End of i loop
}

// The pruning of WTSUM at the R2 cutoff c: every PGS is put in the group of the kept PGS
// that dropped it.  Returns the number of kept PGSs (= groups).
static int PGSPruneGroups(Matrix &R, int n, double c, std::vector<int> &group)
{
  group.assign(n, -1);
  int kept = 0;
  for (int a = 0; a < n; a++) {
    if (group[a] >= 0) continue;
    group[a] = kept;
    for (int b = a + 1; b < n; b++)
      if (group[b] < 0 && R[a][b] != _NAN_ && R[a][b] * R[a][b] > c) group[b] = kept;
    kept++;
  }   // End of a loop
  return kept;
}

// number of pairs among m items (a function, not a lambda: C++ Builder 6 has no lambdas)
static double PairCount(double m) { return m * (m - 1.0) / 2.0; }

static double AdjustedRandIndex(const std::vector<int> &x, int kx, const std::vector<int> &y, int ky, int n)
{
  std::vector<double> cell((size_t)kx * ky, 0.0), rowSum(kx, 0.0), colSum(ky, 0.0);
  for (int i = 0; i < n; i++) {
    cell[(size_t)x[i] * ky + y[i]] += 1.0;
    rowSum[x[i]] += 1.0;
    colSum[y[i]] += 1.0;
  }   // End of i loop
  double sij = 0.0, sa = 0.0, sb = 0.0;
  for (size_t t = 0; t < cell.size(); t++) sij += PairCount(cell[t]);
  for (int a = 0; a < kx; a++) sa += PairCount(rowSum[a]);
  for (int b = 0; b < ky; b++) sb += PairCount(colSum[b]);
  double expected = sa * sb / PairCount((double) n), maxIndex = 0.5 * (sa + sb);
  if (maxIndex - expected < 1e-12) return sij == maxIndex ? 1.0 : 0.0;
  return (sij - expected) / (maxIndex - expected);
}

// R: correlations among the top PGSs in rank order.  currentR2: the --pruner2 in use.
// Prints the clusters and the recommended --pruner2, and sets gPrunerR2Rec / gPGSClusterCount.
static void ReportPruneAdvice(Matrix &R, double currentR2)
{
  gPrunerR2Rec = _NAN_;
  gPGSClusterCount = 0;
  int n = R.rows < 200 ? R.rows : 200;
  if (n < 4) return;
  std::vector<double> D((size_t) n * n, 0.0);
  for (int a = 0; a < n; a++)
    for (int b = 0; b < n; b++) {
      if (a == b) continue;
      double x = (R[a][b] == _NAN_ || R[b][a] == _NAN_) ? 0.0 : 0.5 * (fabs(R[a][b]) + fabs(R[b][a]));
      D[(size_t) a * n + b] = 1.0 - x;
    }   // End of b loop
  std::vector<int> mergeA, mergeB, label, group;
  std::vector<double> height;
  PGSAverageLinkage(D, n, mergeA, mergeB, height);
  int bestK = 0;
  double bestSil = -2.0;
  for (int k = 2; k <= n - 1; k++) {
    PGSClusterLabels(n, mergeA, mergeB, k, label);
    double s = PGSMeanSilhouette(D, n, label, k);
    if (s > bestSil + 1e-12) { bestSil = s; bestK = k; }
  }   // End of k loop
  int currentKept = PGSPruneGroups(R, n, currentR2, group);
  if (bestSil < 0.25) {
    printf("  PGS correlation clusters: no clear cluster structure among the top %d PGSs (best mean silhouette %.2f);\n"
           "  no change of --pruner2 %.2f is recommended (it keeps %d of the %d)\n", n, bestSil, currentR2, currentKept, n);
    return;
  }
  PGSClusterLabels(n, mergeA, mergeB, bestK, label);
  double mid = 0.5 * (height[n - bestK - 1] + height[n - bestK]);   // between the merges leaving k+1 and k clusters
  double c0 = (1.0 - mid) * (1.0 - mid);
  if (c0 < 0.01) c0 = 0.01;
  if (c0 > 0.95) c0 = 0.95;
  double bestC = c0, bestARI = -2.0;
  int bestKept = 0, bestStep = 1;
  std::vector<double> ariAt(96, -2.0);
  for (int step = 1; step <= 95; step++) {
    double c = step * 0.01;
    int kept = PGSPruneGroups(R, n, c, group);
    double ari = AdjustedRandIndex(group, kept, label, bestK, n);
    ariAt[step] = ari;
    if (ari > bestARI + 1e-9 || (ari > bestARI - 1e-9 && fabs(c - c0) < fabs(bestC - c0))) {
      bestARI = ari; bestC = c; bestKept = kept; bestStep = step;
    }
  }   // End of step loop
  int lowStep = bestStep, highStep = bestStep;   // the cutoffs around the best one that do as well
  while (lowStep > 1 && ariAt[lowStep - 1] > bestARI - 1e-9) lowStep--;
  while (highStep < 95 && ariAt[highStep + 1] > bestARI - 1e-9) highStep++;
  gPrunerR2Rec = bestC;
  gPGSClusterCount = bestK;
  printf("  PGS correlation clusters: the top %d PGSs form %d clusters (average linkage on 1-|r|, mean silhouette %.2f)\n", n, bestK, bestSil);
  printf("  Recommended --pruner2 %.2f for WTSUM: it keeps %d of the %d PGSs, matching the clusters (adjusted Rand index %.2f;\n"
         "  the same result for --pruner2 %.2f to %.2f); the current --pruner2 %.2f keeps %d\n",
         bestC, bestKept, n, bestARI, lowStep * 0.01, highStep * 0.01, currentR2, currentKept);
}

// Writes the correlations among the top PGSs of one trait.
//
//   outfile0   <prefix>_<trait>; the files are <outfile0>_PGScorr*.txt
//   topIndex   the top PGSs (rows of risks), best first
//   samples    the samples used for the trait (columns of risks)
//   status     empty for a quantitative trait or survival outcome.  For a case-control
//              trait, one entry per sample: 1 = case, 0 = control.
//
// Always written:        _PGScorr.txt           all samples
// Case-control only:     _PGScorr_control.txt   controls only
//                        _PGScorr_case.txt      cases only
//                        _PGScorr_diff.txt      pair-by-pair test of case vs control
//                                               (first line: the Jennrich test of the whole matrix)
// caseControlDone is set when the case/control files were written.
// pruneR2 is the --pruner2 in use by a scan whose WTSUM prunes for correlation (--linear,
// --logistic): the clusters among the top PGSs and a recommended --pruner2 are printed
// (see ReportPruneAdvice).  Pass _NAN_ for the other scans.
void WritePGSCorrelations(const char *outfile0, Matrix &risks, StringArray &rsNames, IntArray &topIndex,
  IntArray &samples, IntArray &status, int coreCount, bool &caseControlDone, double pruneR2)
{
  gPrunerR2Rec = _NAN_;
  gPGSClusterCount = 0;
  caseControlDone = false;
  int p = topIndex.Length();
  if(p < 1 || samples.Length() < 4) return;
  StringArray names(p);
  for(int a = 0; a < p; a++)
    names[a] = rsNames[topIndex[a]];
  Matrix R, N;
  PGSCorrelation(risks, topIndex, samples, R, N, coreCount);
  String file(outfile0);
  file.Add("_PGScorr.txt");
  WritePGSMatrix(file, R, names);
  printf("  Correlations among the top %d PGSs (%d samples) write to %s\n", p, samples.Length(), (const char*)file);
  if(pruneR2 != _NAN_)
    ReportPruneAdvice(R, pruneR2);
  if(status.Length() != samples.Length() || p < 2) return;
  IntArray cas(0), ctl(0);
  for(int k = 0; k < samples.Length(); k++)
    if(status[k] == 1) cas.Push(samples[k]);
    else if(status[k] == 0) ctl.Push(samples[k]);
  if(cas.Length() < 10 || ctl.Length() < 10){
    printf("  Case-only and control-only PGS correlations skipped: %d cases, %d controls (at least 10 each are needed)\n",
      cas.Length(), ctl.Length());
    return;
  }
  Matrix Rt, Nt, Rc, Nc;
  PGSCorrelation(risks, topIndex, ctl, Rt, Nt, coreCount);
  PGSCorrelation(risks, topIndex, cas, Rc, Nc, coreCount);
  String file1(outfile0), file2(outfile0), file3(outfile0);
  file1.Add("_PGScorr_control.txt");
  file2.Add("_PGScorr_case.txt");
  file3.Add("_PGScorr_diff.txt");
  WritePGSMatrix(file1, Rt, names);
  WritePGSMatrix(file2, Rc, names);
  printf("  Control-only (%d) and case-only (%d) PGS correlations write to %s and %s\n",
    ctl.Length(), cas.Length(), (const char*)file1, (const char*)file2);
  Matrix Zd, Pd, Padj, Pbon;
  double chi2, pJ;
  int df;
  PGSCorrelationTest(Rt, Nt, Rc, Nc, Zd, Pd, Padj, Pbon, chi2, df, pJ);
  // pairs, most different first
  IntArray pa(0), pb(0);
  Vector pv(0);
  for(int a = 0; a < p; a++)
    for(int b = a + 1; b < p; b++)
      if(Pd[a][b] != _NAN_) { pa.Push(a); pb.Push(b); pv.Push(Pd[a][b]); }
  int m = pv.Length();
  QuickIndex idx;
  if(m > 0) idx.Index(pv);
  FILE *fp = fopen(file3, "wb");
  if(fp == NULL) error("Cannot open %s to write.", (const char*)file3);
  fprintf(fp, "#Jennrich_chi2=");
  if(chi2 == _NAN_) fprintf(fp, "NA df=%d P=NA", df);
  else fprintf(fp, "%.4lf df=%d P=%.4G", chi2, df, pJ);
  fprintf(fp, " NControl=%d NCase=%d\n", ctl.Length(), cas.Length());
  fprintf(fp, "PGS1\tPGS2\trControl\trCase\tDiff\tZ\tPvalue\tP_BH\tP_Bonferroni\n");
  int sigBH = 0, sigBon = 0;
  for(int j = 0; j < m; j++){
    int k = idx[j];
    int a = pa[k], b = pb[k];
    fprintf(fp, "%s\t%s\t%.4lf\t%.4lf\t%.4lf\t%.3lf\t%.3G\t%.3G\t%.3G\n",
      (const char*)names[a], (const char*)names[b], Rt[a][b], Rc[a][b], Rc[a][b] - Rt[a][b],
      Zd[a][b], Pd[a][b], Padj[a][b], Pbon[a][b]);
    if(Padj[a][b] < 0.05) sigBH++;
    if(Pbon[a][b] < 0.05) sigBon++;
  } // End of j loop
  fclose(fp);
  printf("  Difference in PGS correlations between cases and controls write to %s\n", (const char*)file3);
  if(chi2 == _NAN_)
    printf("  Jennrich test of equal correlation matrices not available (missing correlation or singular matrix)\n");
  else
    printf("  Jennrich test of equal PGS correlation matrices, cases vs controls: chi2=%.2lf (df=%d), P=%.3G\n", chi2, df, pJ);
  printf("  %d of %d PGS pairs differ at FDR < 0.05 (%d at Bonferroni < 0.05)\n", sigBH, m, sigBon);
  caseControlDone = true;
}
