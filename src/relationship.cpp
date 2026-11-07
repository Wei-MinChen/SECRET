////////////////////////////////////////////////////////////// ////////
// relationship.cpp
// (c) 2023-2024 Wei-Min Chen
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
// May 13, 2024

#include <time.h>
#include <math.h>
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

void WriteKinships(StringArray &sampleID, IntArray *relatives, Vector *kinships, const char *prefix, int coreCount)
{
    char buffer[0x20000];
    String outfile(prefix);
    outfile.Add(".kin");
    FILE *fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    int pbuffer = 0;
    pbuffer += sprintf(&buffer[pbuffer], "ID1\tID2\tKinship\n");
    long long int totalpairCount = 0;
    for (int c = 0; c < coreCount; c++) {
        int count = kinships[c].Length();
        for (int i = 0; i < count; i++) {
            pbuffer += sprintf(&buffer[pbuffer], "%s\t%s\t%.4lf\n",
                (const char*)sampleID[relatives[c][i * 2]], 
                (const char*)sampleID[relatives[c][i * 2 + 1]], 
                kinships[c][i]);
            if (pbuffer > 0xFFFF) {
                fwrite(buffer, 1, pbuffer, fp);
                pbuffer = 0;
            }
        }   // End of i loop
        totalpairCount += count;
    }   // End of c loop
    if (pbuffer) fwrite(buffer, 1, pbuffer, fp);
    fclose(fp);
    printf("Kinship data saved in file %s\n", (const char*)outfile);
}

void RelationshipInference(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *covariateList, IntArray *allrelatives, Vector *allkinships, int coreCount)
{
    printf("\nOptions in effect:\n");
    printf("\t--kinship\n");
    if(phefile[0])
        printf("\t--pheno %s\n", phefile);
    if (gCovarFile[0])
        printf("\t--covarfile %s\n", gCovarFile);
    if (covariateList[0])
        printf("\t--covariate %s\n", covariateList);
    printf("\t--score %s\n", rsfile);
    if (coreCount)
        printf("\t--cpus %d\n", coreCount);
    if (prefix)
        printf("\t--prefix %s\n", prefix);
    printf("\n");

    StringArray sampleID, alltraitNames, rsNames, annotation;
    IntArray covariates, affections, quantraits;
    Matrix alltraits, risks;
    SECRET_Read(phefile, dictfile, rsfile, NULL, prefix, covariateList,
        sampleID, alltraitNames, alltraits, covariates, affections, quantraits,
        rsNames, risks, annotation, coreCount);
    int idCount = sampleID.Length();
    int rsCount = rsNames.Length();
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int r = 0; r < rsCount; r++)
        risks[r].Multiply(1 / sqrt(risks[r].Var()));
    Matrix trisks(idCount, rsCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int i = 0; i < idCount; i++)
        for (int r = 0; r < rsCount; r++)
            trisks[i][r] = risks[r][i];
    double minKinship = 0.1767767;
    double T = (2 - minKinship * 4) * rsCount;
    IntArray *relatives = new IntArray[coreCount];
    Vector *kinships = new Vector[coreCount];
    for (int c = 0; c < coreCount; c++) {
        relatives[c].Dimension(0);
        kinships[c].Dimension(0);
    }
    int BLOCKSIZE = 32;
    char SHIFTSIZE = 5; // takes up to 2.1 million samples
    if (idCount > (1 << 20)) {  // sample size over a million
        int count = (idCount >> 20);
        int F = 0;
        for (; count; F++) count >>= 1;
        BLOCKSIZE <<= F;
        SHIFTSIZE += F;
    }
    unsigned int blockCount = (idCount - 1) / BLOCKSIZE + 1;
    unsigned int loopIndexLength = 
        ((blockCount & 1) ? blockCount * ((blockCount + 1) / 2) : (blockCount / 2) * (blockCount + 1));
    unsigned int *loopIndex = new unsigned int[loopIndexLength];
    int index = 0;
    for (unsigned int i = 0; i < idCount; i += BLOCKSIZE) {
        unsigned int iShifted = (i << (16 - SHIFTSIZE));
        for (int j = i; j < idCount; j += BLOCKSIZE)
            loopIndex[index++] = iShifted | (j >> SHIFTSIZE);
    }   // End of i loop
    int firstCount = loopIndexLength / coreCount;
    int onepercent = firstCount ? (firstCount + 99) / 100 : 1;
    time_t timer;
    struct tm *tblock;
    timer = time(NULL);
    tblock = localtime(&timer);
    printf("Relationship inference starts at %s", asctime(tblock));
    int thread = 0;
#ifdef _OPENMP
    printf("%d CPU cores are used.\n", coreCount);
#pragma omp parallel num_threads(coreCount) private(thread)
{
    thread = omp_get_thread_num();
    kinships[thread].Dimension(0);
    relatives[thread].Dimension(0);
#pragma omp for
#endif
    for (long int k = 0; k < loopIndexLength; k++) {
        if (k < firstCount && (k % onepercent) == 0) {
            printf("%d%%\r", k / onepercent);
            fflush(stdout);
        }
        int i = (loopIndex[k] >> 16) << SHIFTSIZE;
        int iMax = i < idCount - BLOCKSIZE ? i + BLOCKSIZE : idCount;
        int j = (loopIndex[k] & 0xFFFF) << SHIFTSIZE;
        int jMax = j < idCount - BLOCKSIZE ? j + BLOCKSIZE : idCount;
        int jMin = j;
        for (int ii = i; ii < iMax; ii++) {
            if (i == j) jMin = ii + 1;
            for (int jj = jMin; jj < jMax; jj++) {
                double dist = 0.0;
                for (int r = 0; r < rsCount && dist < T; r++) {
                    double diff = trisks[ii][r] - trisks[jj][r];
                    dist += diff * diff;
                }   // End of r loop
                double kinship = 0.5 - dist * 0.25 / rsCount;
                if (kinship >= minKinship) {
                    relatives[thread].Push(ii);
                    relatives[thread].Push(jj);
                    kinships[thread].Push(kinship);
                }
            }   // End of i2 loop
        }   // End of i1 loop
    }   // End of k loop
#ifdef _OPENMP
}
#endif
    int relativeCount = kinships[0].Length();
    for (int c = 1; c < coreCount; c++)
        relativeCount += kinships[c].Length();
    IntArray relative1(0);
    if (relativeCount) {
        printf("%d pairs are potentially close relatives, e.g.,\n", relativeCount);
        printf("  ID1\tID2\tKinship\n");
        int totalCount = 0;
        for (int c = 0; c < coreCount; c++) {
            int count = kinships[c].Length();
            for (int i = 0; i < count && totalCount < 10; i++)
                if (kinships[c][i] > 0.25) {
                    printf("  %s\t%s\t%.3lf\n",
                        (const char*)sampleID[relatives[c][i * 2]],
                        (const char*)sampleID[relatives[c][i * 2 + 1]],
                        kinships[c][i]);
                    totalCount++;
                }
        }   // End of c loop
        printf(" A more stringent kinship cutoff may be needed, such as 0.25.\n");
        WriteKinships(sampleID, relatives, kinships, prefix, coreCount);
    }
    else printf("No close relatives are inferred.\n");
    timer = time(NULL);
    tblock = localtime(&timer);
    printf("Relationship inference ends at %s", asctime(tblock));
}
