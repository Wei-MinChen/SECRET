////////////////////////////////////////////////////////////// ////////
// structure.cpp
// (c) 2023-2025 Wei-Min Chen
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
// Nov 6, 2025

#include <time.h>
#include <math.h>
#include <stdint.h>
#include <memory.h>
#include "MathSVD.h"
#include "StringArray.h"
#include "StringHash.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "DataInput.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef WITH_LAPACK
extern "C" void dsyevr_(char*, char*, char*, int*, double*, int*, double*, double*, int*, int*, double*, int*, double*, double*, int*, int*, double*, int*, int*, int*, int*);
extern "C" double dlamch_(char*);
#endif

#ifdef WITH_OPENBLAS
extern "C" void openblas_set_num_threads(int);
#endif

extern int plotRSQ(const char *input, const char *input2, const char *output, const char *rpath);
extern void plotPopStructure(const char *prefix, int projectFlag, const char *rpath);

void TopEigenForPCA(Matrix & D, int nPC, Vector & Values, Matrix & Vectors, int coreCount)
{
    uintptr_t dimN = D.cols;
    int dimPC = dimN > nPC ? nPC : dimN;
    Values.Dimension(dimPC);
    Vectors.Dimension(dimPC, dimN);
#ifdef WITH_LAPACK
    printf("  LAPACK is being used...\n"); fflush(stdout);
#ifdef WITH_OPENBLAS
    openblas_set_num_threads(coreCount);
#endif
    char JOBZ = 'V';
    char RANGE = 'I';
    char UPLO = 'U';
    int N = dimN;
    int info;
    double *A = new double[dimN*dimN];
    int LDA = N;
    for (uintptr_t i = 0; i < dimN; i++)    // transpose is actually used for symmetry
        memcpy(&(A[i*dimN]), &(D[i][0]), dimN * sizeof(double));//assign ith row to ith col
    double VL, VU;
    int M = dimPC;
    int IL = N - dimPC + 1;
    int IU = N;
    char CMACH = 'S';
    double ABSTOL = dlamch_(&CMACH);
    double *W = new double[N];
    double *Z = new double[M*N]; // N rows x M columns
    int LDZ = N;
    int *ISUPPZ = new int[M * 2];
    int INFO;
    double optim_lwork;
    int optim_liwork;
    int LWORK = -1;
    int LIWORK = -1;
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, &optim_lwork, &LWORK, &optim_liwork, &LIWORK, &INFO);
    LWORK = (int)optim_lwork;
    LIWORK = optim_liwork;
    double *WORK = new double[LWORK];
    int *IWORK = new int[LIWORK];
    dsyevr_(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, ISUPPZ, WORK, &LWORK, IWORK, &LIWORK, &INFO);
    if (INFO != 0) error("SVD failed with INFO=%d.", INFO);
    delete[]ISUPPZ;
    delete[]IWORK;
    delete[]WORK;
    delete[]A;
    for (int j = 0; j < dimPC; j++) {
        Values[j] = W[dimPC - 1 - j];   //Next: assign (M-j)th col of Z to jth row of Vectors
        memcpy(&(Vectors[j][0]), &(Z[dimN*(dimPC - 1 - j)]), dimN * sizeof(double));
    }// for(int i = 0; i < N; i++) for(int j = 0; j < dimPC; j++) Vectors[i][j] = Z[(dimPC - 1 - j)*N + i];
    delete[]W;
    delete[]Z;
#else
    printf("  Please re-compile KING with LAPACK library.\n");
    SVD svd;
    int N = dimN;
    svd.Decompose(D);
    printf("done\n");
    if (svd.n == 0) return;
    QuickIndex idx;
    idx.Index(svd.w);
    for (int j = 0; j < dimPC; j++)
        Values[j] = svd.w[idx[N - 1 - j]];
    printf("\n");
    for (int i = 0; i < N; i++)
        for (int j = 0; j < dimPC; j++)
            Vectors[j][i] = svd.v[i][idx[N - 1 - j]];
#endif
}

void localPCA(Matrix &risks, int dimPC, Vector &Values, Matrix &PCs, int coreCount)
{
    int N = risks.rows;
    int rsCount = risks.cols;
    Vector mean(rsCount), sd(rsCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int r = 0; r < rsCount; r++){ 
        double sum = 0.0;
        double ss = 0.0;
        for (int i = 0; i < N; i++) {
            sum += risks[i][r];
            ss += risks[i][r] * risks[i][r];
        }   // End of i loop
        mean[r] = sum / N;
        sd[r] = sqrt((ss - sum * sum / N) / N);
    }   // End of OMP r loop
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int i = 0; i < N; i++)
        for (int r = 0; r < rsCount; r++)
            risks[i][r] = (risks[i][r] - mean[r]) / sd[r];
    Matrix IP(N, N);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for(int i = 0; i < N; i++)
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int r = 0; r < rsCount; r++)// computation O(N*N*R)
                sum += risks[i][r] * risks[j][r];
            IP[i][j] = sum;
        }   // End of OMP ij loop
    TopEigenForPCA(IP, dimPC, Values, PCs, coreCount);
}

void runPCA(const char *rsfile, const char *prefix, bool rplotFlag, int coreCount)
{
    printf("\nOptions in effect:\n");
    printf("\t--popstruct\n");
    printf("\t--score %s\n", rsfile);
    if (rplotFlag)
        printf("\t--rplot\n");
    if (coreCount)
        printf("\t--cpus %d\n", coreCount);
    if (prefix)
        printf("\t--prefix %s\n", prefix);
    printf("\n");
    time_t timer;
    struct tm *tblock;
    timer = time(NULL);
    tblock = localtime(&timer);
    printf("Population structure inference starts at %s", asctime(tblock));
    StringArray sampleID;
    Matrix risks;
    StringArray rsNames;
    ReadScoreMatrix(rsfile, rsNames, sampleID, risks, coreCount);
    int rsCount = rsNames.Length();
    if (rsCount < 20) error("Too few risk scores (M < 20).");
    int idCount = sampleID.Length();
    if (idCount < 10) error("Too few samples (N = %d, less than 10)", idCount);
    printf("  %d risk scores for %d samples loaded.\n", rsCount, idCount);
    Vector vars(rsCount), temprisks(idCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(temprisks)
#endif
    for (int r = 0; r < rsCount; r++) {
        temprisks.Dimension(idCount);
        for (int i = 0; i < idCount; i++)
            temprisks[i] = risks[i][r];
        vars[r] = temprisks.Var();
    }   // End of OMP r loop
    IntArray invalidIndex(0), validIndex(0);
    for (int r = 0; r < rsCount; r++)
        if (vars[r] < 1E-10 || vars[r] > 10000) invalidIndex.Push(r);
        else validIndex.Push(r);
    int validrsCount = validIndex.Length();
    if (validrsCount < rsCount) 
        printf("  %d risk scores removed for too little or too much variation.\n", rsCount - validrsCount);
    printf("  PCA on %d samples x %d risk scores...\n", idCount, validrsCount); fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int i = 0; i < idCount; i++)
        for (int r = 0; r < validrsCount; r++)
            risks[i][r] = (risks[i][validIndex[r]] == _NAN_ ? _NAN_ : risks[i][validIndex[r]] / sqrt(vars[validIndex[r]]));
    int dimPC = 4;
    Vector Values;
    Matrix PCs;
    localPCA(risks, dimPC, Values, PCs, coreCount);

    printf("  Largest %d eigenvalues:", dimPC);
    for (int k = 0; k < dimPC; k++) printf(" %.2lf", sqrt(Values[k]));
    printf("\n");
    String outfile(prefix);
    outfile.Add("_pc.txt");
    FILE *fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    int pbuffer = 0;
    char buffer[0x20000];
    pbuffer += sprintf(&buffer[pbuffer], "ID");
    for (int r = 0; r < dimPC; r++)
        pbuffer += sprintf(&buffer[pbuffer], "\tPC%d", r + 1);
    for (int i = 0; i < idCount; i++) {
        pbuffer += sprintf(&buffer[pbuffer], "\n%s", (const char*)sampleID[i]);
        for (int r = 0; r < dimPC; r++)
            pbuffer += sprintf(&buffer[pbuffer], "\t%.4lf", PCs[r][i]);
        if (pbuffer > 0xFFFF) {  // buffer big enough for writing
            fwrite(buffer, 1, pbuffer, fp);
            pbuffer = 0;
        }
    }   // End of i loop
    pbuffer += sprintf(&buffer[pbuffer], "\n");
    fwrite(buffer, 1, pbuffer, fp);
    fclose(fp);
    printf("  TOP %d PCs of population structure writes to %s\n", dimPC, (const char*)outfile);
    if (rplotFlag) plotPopStructure(prefix, false, "");
    timer = time(NULL);
    tblock = localtime(&timer);
    printf("Population structure inference ends at %s", asctime(tblock));
}
