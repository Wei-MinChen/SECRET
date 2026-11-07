////////////////////////////////////////////////////////////////////////
// Logistic.cpp
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
// Oct 28, 2026
//
// --logistic: logistic regression of a binary trait on one polygenic
// score at a time, with the --covariate columns carried in the model.
// See Logistic.h for the coding conventions.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Logistic.h"
#include "Constant.h"
#include "Error.h"
#include "DataInput.h"
#include "StringArray.h"
#include "StringBasics.h"
#include "QuickIndex.h"
#include "MathStats.h"
#include "PenalizedLogit.h"
#include "Random.h"
#include <stdlib.h>

extern void MakePhenos(Matrix & risks, Matrix & alltraits,
   StringArray & alltraitNames, IntArray & affections, IntArray & quantraits,
   IntArray & covariates, Matrix & allY, Matrix & allX,
   unsigned long long int ** masks, int coreCount);
extern void MakeBinaryLD(Matrix & risks, double ldT,
   unsigned long long int ** LD, Vector & rsSD, int coreCount);
extern void MakePrunedOrder(unsigned long long int ** LD, IntArray & order);
extern void LDSummary(unsigned long long int ** LD, int rsCount,
   int wordCount, int coreCount);
extern void WriteTable(String & outfile, Matrix & table, StringArray & header,
   StringArray & col0);
extern int plotLogistic(const char *input, const char *traitName, const char *rpath, String &msg);
extern int plotCalibrationLogistic(const char *input, const char *traitName, const char *rpath, String &msg);
extern double GenomicControlNumber(Vector &Z, IntArray &sortedrsIndex);
extern void WritePGSCorrelations(const char *outfile0, Matrix &risks, StringArray &rsNames, IntArray &topIndex,
  IntArray &samples, IntArray &status, int coreCount, bool &caseControlDone, double pruneR2);
extern int plotPGSCorr(const char *input, const char *statfile, const char *traitName, const char *rpath, String &msg);

static double sigmoid(double z) {
    if (z >= 0.0) return 1.0 / (1.0 + exp(-z));
    double e = exp(z);
    return e / (1.0 + e);
}

/* In-place Cholesky factorization of SPD matrix A (m x m, row-major).
 * Reads the lower triangle; leaves L in the lower triangle.
 * Returns 0 on success, -1 if A is not (numerically) positive definite:
 * a pivot is rejected when it falls below 1e-10 of its original diagonal. */
static int cholesky_factor(double *A, int m) {
    for (int i = 0; i < m; i++) {
        const double d0 = A[i * m + i];      /* original diagonal (not yet overwritten) */
        for (int j = 0; j <= i; j++) {
            double s = A[i * m + j];
            for (int k = 0; k < j; k++) s -= A[i * m + k] * A[j * m + k];
            if (i == j) {
                /* relative pivot test: scale-free, also rejects NaN */
                if (!(s > 1e-10 * d0)) return -1;
                A[i * m + i] = sqrt(s);
            } else {
                A[i * m + j] = s / A[j * m + j];
            }
        }
    }
    return 0;
}

/* Solve (L L^T) x = b given the factor L from cholesky_factor. */
static void cholesky_solve_factored(const double *L, const double *b,
                                    double *x, int m) {
    for (int i = 0; i < m; i++) {            /* forward: L z = b */
        double s = b[i];
        for (int k = 0; k < i; k++) s -= L[i * m + k] * x[k];
        x[i] = s / L[i * m + i];
    }
    for (int i = m - 1; i >= 0; i--) {       /* backward: L^T x = z */
        double s = x[i];
        for (int k = i + 1; k < m; k++) s -= L[k * m + i] * x[k];
        x[i] = s / L[i * m + i];
    }
}

/* Accumulate gradient g and (lower-triangle) Hessian H at beta.
 * H and g must have length m*m and m, where m = d + 1. */
static void grad_hess(const double *X, const int *y, int n, int d,
                      double l2, const double *beta,
                      double *g, double *H, double *xa) {
    const int m = d + 1;
    for (int j = 0; j < m; j++) g[j] = 0.0;
    for (int j = 0; j < m * m; j++) H[j] = 0.0;
    xa[0] = 1.0;

    for (int i = 0; i < n; i++) {
        const double *xi = X + (size_t)i * d;
        for (int j = 0; j < d; j++) xa[j + 1] = xi[j];

        double z = 0.0;
        for (int j = 0; j < m; j++) z += beta[j] * xa[j];
        double p = sigmoid(z);
        double w = p * (1.0 - p);
        if (w < 1e-10) w = 1e-10;

        double r = y[i] - p;
        for (int j = 0; j < m; j++) {
            g[j] += r * xa[j];
            for (int k = 0; k <= j; k++)
                H[j * m + k] += w * xa[j] * xa[k];
        }
    }
    for (int j = 1; j < m; j++) {            /* L2 on weights, not intercept */
        g[j] -= l2 * beta[j];
        H[j * m + j] += l2;
    }
}
/* Quasi-/complete-separation screen.
 *
 * When a predictor separates the outcome the MLE does not exist and the
 * fitted slope runs off to infinity.  Two cheap indicators are used: an
 * implausibly large slope, or fitted probabilities pinned at 0/1 for a
 * material fraction of the sample.  The intercept is exempt, because it is
 * legitimately large when the outcome is rare.
 * Returns 1 when the fit looks separated, 0 otherwise. */
static int separated_check(const double *X, int n, int d,
                           const double *beta) {
    const int m = d + 1;
    /* A slope is judged per SD of its predictor: raw slopes are not comparable
     * across covariates of different scales (e.g. PCs with SD ~0.002). */
    for (int j = 1; j < m; j++) {
        double sum = 0.0, ss = 0.0;
        for (int i = 0; i < n; i++) {
            double x = X[(size_t)i * d + (j - 1)];
            sum += x;
            ss += x * x;
        }
        double mean = sum / n;
        double var = ss / n - mean * mean;
        double sd = var > 0.0 ? sqrt(var) : 0.0;
        if (fabs(beta[j]) * sd > 15.0) return 1;
    }

    int extreme = 0;
    for (int i = 0; i < n; i++) {
        double z = beta[0];
        for (int j = 0; j < d; j++) z += beta[j + 1] * X[(size_t)i * d + j];
        double p = sigmoid(z);
        if (p < 1e-8 || p > 1.0 - 1e-8) extreme++;
    }
    return (extreme > 0 && extreme * 100 > n) ? 1 : 0;
}

/*
 * Newton / IRLS logistic regression with standard errors.
 *
 *   beta  output, length d + 1 (beta[0] = intercept)
 *   se    output, length d + 1, SE of each coefficient (pass NULL to skip)
 *   cov   output, (d+1) x (d+1) row-major covariance matrix
 *         (pass NULL to skip; requires se != NULL is NOT needed)
 *
 * Returns the number of iterations when the fit converged to a usable
 * estimate, or a negative code otherwise:
 *   -1  allocation failure or non-positive-definite Hessian
 *   -2  hit max_iter without converging
 *   -3  quasi-/complete separation detected (the estimate does not exist)
 * SEs/cov are only filled when the return value is >= 0.
 */
int logreg_irls(const double *X, const int *y, int n, int d,
                double l2, int max_iter, double tol,
                double *beta, double *se, double *cov) {
    const int m = d + 1;
    double *H     = new double[m*m];
    double *g     = new double[m];
    double *delta = new double[m];
    double *xa    = new double[m];
    double *e     = new double[m];
    double *col   = new double[m];
    int iters = -1;

    if (!H || !g || !delta || !xa || !e || !col) goto done;

    for (int j = 0; j < m; j++) beta[j] = 0.0;

    /* ---- Newton iterations ---- */
    for (int it = 1; it <= max_iter; it++) {
        grad_hess(X, y, n, d, l2, beta, g, H, xa);

        if (cholesky_factor(H, m) != 0) { iters = -1; goto done; }
        cholesky_solve_factored(H, g, delta, m);

        double max_step = 0.0;
        for (int j = 0; j < m; j++) {
            beta[j] += delta[j];
            if (fabs(delta[j]) > max_step) max_step = fabs(delta[j]);
        }
        if (max_step < tol) { iters = it; break; }
        if (it == max_iter) iters = -2;      /* hit the cap without converging */
    }

    /* ---- separation screen: never report an estimate that does not exist ---- */
    if (iters >= 0 && separated_check(X, n, d, beta)) {
        iters = -3;
        goto done;
    }

    /* ---- Standard errors from the inverse Hessian at the final beta ---- */
    if (iters >= 0 && (se || cov)) {
        grad_hess(X, y, n, d, l2, beta, g, H, xa);
        if (cholesky_factor(H, m) != 0) { iters = -1; goto done; }

        for (int k = 0; k < m; k++) {
            for (int j = 0; j < m; j++) e[j] = (j == k) ? 1.0 : 0.0;
            cholesky_solve_factored(H, e, col, m);   /* column k of H^-1 */
            if (se) se[k] = sqrt(col[k]);
            if (cov)
                for (int j = 0; j < m; j++) cov[j * m + k] = col[j];
        }
    }

done:
    delete[]H; delete[]g; delete[]delta; delete[]xa; delete[]e; delete[]col;
    return iters;
}

// ---------------------------------------------------------------------
// MakeVariable
// ---------------------------------------------------------------------
// Prepare, for every binary trait a, the sample list, outcome and covariate
// matrix that Logistic_BuildModel / Logistic_ReadModel work on.  Sample i is
// used for trait a when bit a of masks[i] is set AND none of its covariates is
// missing.  Restricting the analysis to a subset of samples (e.g., the train or
// the test set of a cross-validation round) is done by clearing the bits of the
// other samples in masks.
//
//   allX     idCount x covariateCount covariates (_NAN_ = missing)
//   allY     idCount x traitCount phenotypes (SECRET's 1/2 coding)
//   masks    idCount x ceil(traitCount/64) words, bit a = sample has trait a
//   pheno    output, pheno[a] = indices of the samples (columns of risks) used
//   Xs       output, Xs[a] is row-major N x (covariateCount+1): the covariates
//            of the k-th sample are Xs[a][k*(covariateCount+1) + c]; the last
//            column is scratch space for the PGS and is set to 0 here
//   Ys       output, Ys[a][k] = 1 if affected (2), 0 if unaffected (1)
//
// pheno, Xs and Ys must each have traitCount elements.
void MakeVariable(Matrix &allX, Matrix &allY, unsigned long long int **masks,
  IntArray *pheno, Vector *Xs, IntArray *Ys)
{
  const int idCount = allY.rows;
  const int traitCount = allY.cols;
  const int covariateCount = allX.cols;
  const int d = covariateCount + 1;
  for (int a = 0; a < traitCount; a++){
    pheno[a].Dimension(0);
    for (int i = 0; i < idCount; i++){
      if (!((masks[i][a >> 6] >> (a & 0x3F)) & 0x1))
        continue;
      bool ok = true;
      for (int c = 0; c < covariateCount; c++)
        if (allX[i][c] == _NAN_){
          ok = false;
          break;
        }
      if (ok) pheno[a].Push(i);
    } // End of i loop
    const int N = pheno[a].Length();
    Ys[a].Dimension(N);
    Xs[a].Dimension(N * d);
    for (int k = 0; k < N; k++){
      const int i = pheno[a][k];
      Ys[a][k] = (allY[i][a] > 1.5) ? 1 : 0;
      for (int c = 0; c < covariateCount; c++)
        Xs[a][k * d + c] = allX[i][c];
      Xs[a][k * d + covariateCount] = 0.0;
    } // End of k loop
  } // End of a loop
}

// ---------------------------------------------------------------------
// Logistic_ReadModel
// ---------------------------------------------------------------------
// Evaluate a set of PGS models (the columns of coef, as written by --logistic
// or returned by Logistic_BuildModel) on ONE binary trait.  For every model
// the combined score is formed, its AUC is taken ("raw"), and the trait is
// then regressed on the covariates plus the combined score, giving the odds
// per unit of the score (beta, se) and the covariate-adjusted AUC.
//
//   pgsName, coef   one name per row of coef; coef is modelCount x methodCount,
//                   a _NAN_ coefficient is skipped
//   rsNames, risks  the PGSs available; model terms not found in rsNames are
//                   skipped
//   pheno, Y, X     the samples of this trait, as made by MakeVariable (the last
//                   column of X is overwritten with the score)
//   rawAUC, adjAUC, beta, se   output, length methodCount, _NAN_ when the model
//                   could not be evaluated (fewer than 10 samples, no cases or
//                   no controls, no usable term, or a failed fit)
//
// Returns the number of model terms that were not found in rsNames.
int Logistic_ReadModel(StringArray &pgsName, Matrix &coef, StringArray &rsNames,
  Matrix &risks, IntArray &pheno, IntArray &Y, Vector &X, int covariateCount,
  Vector &rawAUC, Vector &adjAUC, Vector &beta, Vector &se)
{
  const int N = pheno.Length();
  const int modelrsCount = coef.rows;
  const int methodCount = coef.cols;
  const int d = covariateCount + 1;     // covariates + combined score
  rawAUC.Dimension(methodCount);
  adjAUC.Dimension(methodCount);
  beta.Dimension(methodCount);
  se.Dimension(methodCount);
  rawAUC.Set(_NAN_);
  adjAUC.Set(_NAN_);
  beta.Set(_NAN_);
  se.Set(_NAN_);
  // locate the model terms once
  IntArray rsIndex(modelrsCount);
  int missingTerms = 0;
  for (int r = 0; r < modelrsCount; r++){
    rsIndex[r] = rsNames.Find(pgsName[r]);
    if (rsIndex[r] == -1)
      missingTerms++;
  } // End of r loop
  if (N < 10) return missingTerms;
  int cases = 0;
  Vector yv(N);
  for (int k = 0; k < N; k++){
    if (Y[k] == 1) cases++;
    yv[k] = Y[k];
  }
  if (cases < 1 || cases == N) return missingTerms;
  int * y = &Y[0];
  double * Xm = &X[0];
  double * bm = new double [d+1];
  double * sm = new double [d+1];
  Vector score(N);
  for (int m = 0; m < methodCount; m++){
    score.Zero();
    int used = 0;
    for (int r = 0; r < modelrsCount; r++){
      int rr = rsIndex[r];
      double bj = coef[r][m];
      if (rr == -1 || bj == _NAN_) continue;
      used++;
      for (int k = 0; k < N; k++)
        score[k] += risks[rr][pheno[k]] * bj;
    } // End of r loop
    if (used == 0) continue;
    rawAUC[m] = LogisticAUC(score, yv);
    for (int k = 0; k < N; k++)
      Xm[k * d + covariateCount] = score[k];
    int it = logreg_irls(Xm, y, N, d, 1e-6, 50, 1e-8, bm, sm, NULL);
    if (it >= 0 && sm[d] > 0.0){
      beta[m] = bm[d];
      se[m] = sm[d];
      for (int k = 0; k < N; k++){
        double e = bm[0];
        for (int c = 0; c < covariateCount; c++)
          e += bm[c + 1] * Xm[k * d + c];
        e += bm[d] * Xm[k * d + covariateCount];
        score[k] = e;
      } // End of k loop
      adjAUC[m] = LogisticAUC(score, yv);
    }
  } // End of m loop
  delete []bm;
  delete []sm;
  return missingTerms;
}

// ---------------------------------------------------------------------
// Logistic_BuildModel
// ---------------------------------------------------------------------
// Scan every PGS against one binary trait, then combine the top PGSs into
// BEST / WTSUM / RIDGE / NET / LASSO models.  See Logistic.h for the interface.
int Logistic_BuildModel(Matrix &risks, IntArray &pheno, IntArray &Y, Vector &X,
   int covariateCount, unsigned long long int ** LD, int nTop, int coreCount,
   Matrix &weights, IntArray &modelIndex, LogisticModelStats *stats)
{
  LogisticModelStats localStats;
  if (stats == NULL) stats = &localStats;
  const int rsCount = risks.rows;
  const int N = pheno.Length();
  const int d = covariateCount + 1;     // covariates + the PGS being tested
  weights.Dimension(0, 5);
  modelIndex.Dimension(0);
  stats->nCase = stats->nControl = 0;
  stats->bestIndex = -1;
  stats->aucBest = stats->aucBestRaw = _NAN_;
  stats->wsumStatus = 0;
  stats->wsBeta = stats->wsZ = stats->wsP = _NAN_;
  stats->aucWsum = stats->aucWsumRaw = _NAN_;
  stats->penalizedSkipped = false;
  stats->penalizedTerms = 0;
  stats->wsumTerms = 0;
  for (int m = 0; m < 3; m++)
    stats->netAuc[m] = stats->netAucRaw[m] = _NAN_;
  if (N < 10)
    return LOGISTIC_TOO_FEW_SAMPLES;
  if (N < 5 * (d+1))
    return LOGISTIC_TOO_FEW_FOR_TERMS;
  int * y = &Y[0];
  Vector yvec(N);
  int cases = 0;
  for (int k = 0; k < N; k++){
    if (y[k] == 1) cases++;
    yvec[k] = y[k];
  }
  int controls = N - cases;
  stats->nCase = cases;
  stats->nControl = controls;
  if (cases < 5 || controls < 5)
    return LOGISTIC_TOO_FEW_CASES_CONTROLS;
  // ---- one-PGS-at-a-time scan, adjusted for the covariates ----
  Vector &beta = stats->beta, &se = stats->se, &zval = stats->zval;
  beta.Dimension(rsCount);
  se.Dimension(rsCount);
  zval.Dimension(rsCount);
  beta.Set(_NAN_);
  se.Zero();
  zval.Set(_NAN_);

#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
{
#endif
  double *Xt = new double[N*d];
  memcpy(&Xt[0], &X[0], d*N*sizeof(double));
  double *b = new double [d+1];
  double *s = new double [d+1];
#ifdef _OPENMP
#pragma omp for
#endif
  for (int r = 0; r < rsCount; r++) {
    for (int k = 0; k < N; k++)
      Xt[k*d+covariateCount] = risks[r][pheno[k]];
    int iters = logreg_irls(Xt, y, N, d, 1e-6, 50, 1e-8, b, s, NULL);
    if(iters>=0 && s[d]>0){
      beta[r] = b[d];
      se[r] = s[d];
      zval[r] = b[d] / s[d];
    }
  } // End of r loop OMP
  delete []b;
  delete []s;
  delete []Xt;
#ifdef _OPENMP
} // End omp parallel
#endif
  IntArray &sortedIndex = stats->sortedIndex;
  GenomicControlNumber(zval, sortedIndex);
  // best PGS
  if (rsCount > 0 && beta[sortedIndex[0]] != _NAN_)
    stats->bestIndex = sortedIndex[0];
  const int bestIndex = stats->bestIndex;
  // ---- the top nTop PGSs are the terms of every model; only WTSUM uses the subset that is
  // left after pruning for correlation (--pruner2), as in --linear ----
  IntArray wsumFlag(0);                  // 1 = the PGS is used by WTSUM
  if (nTop > 0 && rsCount > 0){
    int top = nTop < rsCount ? nTop : rsCount;
    for (int j = 0; j < top; j++){
      int rs = sortedIndex[j];
      if (beta[rs] != _NAN_)
        modelIndex.Push(rs);
    } // End of j loop
    IntArray kept = modelIndex;
    MakePrunedOrder(LD, kept);           // positions in modelIndex of the PGSs that survive
    wsumFlag.Dimension(modelIndex.Length());
    wsumFlag.Zero();
    for (int j = 0; j < kept.Length(); j++)
      wsumFlag[kept[j]] = 1;
    stats->wsumTerms = kept.Length();
  }
  const int modelCount = modelIndex.Length();
  Vector wsumW(modelCount);          // WTSUM weights
  wsumW.Zero();
  Matrix netCoefAll(modelCount, 3);  // ridge / net / lasso coefficients
  netCoefAll.Set(_NAN_);
  // AUC of the best single PGS, taken from the FITTED model: the linear
  // predictor carries the covariates, so this is the discrimination of
  // the adjusted model rather than of the raw score.  With no covariates
  // it reduces to ranking on the score alone (same AUC).
  if(bestIndex < 0) return LOGISTIC_OK;

  double *Xb = new double [N*d];
  memcpy(&Xb[0], &X[0], d*N*sizeof(double));
  double *bb = new double [d+1];
  double *sb = new double [d+1];
  for (int k = 0; k < N; k++)
    Xb[k*d+covariateCount] = risks[bestIndex][pheno[k]];
  int itB = logreg_irls(Xb, y, N, d, 1e-6, 50, 1e-8, bb, sb, NULL);
  Vector pred(N);
  if (itB >= 0)
    for (int k = 0; k < N; k++){
      double z = bb[0];
      for (int c = 0; c < covariateCount; c++)
        z += bb[c+1] * X[k*d+c];
      z += bb[covariateCount+1] * risks[bestIndex][pheno[k]];
      pred[k] = z;
    } // End of k loop
  else // adjusted model could not be fitted: fall back to the raw score
    for (int k = 0; k < N; k++)
      pred[k] = risks[bestIndex][pheno[k]];
  stats->aucBest = LogisticAUC(pred, yvec);
  Vector rawPred(N);
  for (int k = 0; k < N; k++)
    rawPred[k] = risks[bestIndex][pheno[k]];
  stats->aucBestRaw = LogisticAUC(rawPred, yvec);
  // ---- WTSUM: weighted sum of the selected PGSs, refitted ----
  Vector wsum(N);
  wsum.Zero();
  double norm = 0.0;
  for (int j = 0; j < modelCount; j++){
    if (!wsumFlag[j]) continue;          // pruned out of WTSUM
    int rs = modelIndex[j];
    for (int k = 0; k < N; k++)
      wsum[k] += beta[rs] * risks[rs][pheno[k]];
    norm += beta[rs] * beta[rs];
  } // End of j loop
  if (norm > 0.0){
    double sc = sqrt(norm);
    for (int k = 0; k < N; k++)
      wsum[k] /= sc;
    for (int j = 0; j < modelCount; j++)
      if (wsumFlag[j])
        wsumW[j] = beta[modelIndex[j]] / sc;
    for(int k = 0; k < N; k++)
      Xb[k*d+covariateCount] = wsum[k];
    int iters = logreg_irls(Xb, y, N, d, 1e-6, 50, 1e-8, bb, sb, NULL);
    if (iters >= 0 && sb[d] > 0.0){
      stats->wsumStatus = 1;
      stats->wsBeta = bb[d];
      stats->wsZ = bb[d] / sb[d];
      stats->wsP = 2.0 * ndist(fabs(stats->wsZ));
      // AUC from the fitted model: covariates + weighted sum
      Vector predW(N);
      for (int k = 0; k < N; k++){
        double z = bb[0];
        for (int c = 0; c < covariateCount; c++)
          z += bb[c+1] * X[k*d+c];
        z += bb[d] * wsum[k];
        predW[k] = z;
      }
      stats->aucWsumRaw = LogisticAUC(wsum, yvec);
      stats->aucWsum = LogisticAUC(predW, yvec);
    }else
      stats->wsumStatus = 2;
  }
  delete []Xb;
  delete []bb;
  delete []sb;
  // ---- penalized models: ridge / elastic net / lasso, always fitted ----
  int total = modelCount + covariateCount;
  stats->penalizedTerms = total;
  if (total >= N)
    stats->penalizedSkipped = true;
  else{
    Matrix Xnet(N, total);
    for (int i = 0; i < N; i++){
      for (int j = 0; j < modelCount; j++)
        Xnet[i][j] = risks[modelIndex[j]][pheno[i]];
      for (int c = 0; c < covariateCount; c++)
        Xnet[i][modelCount + c] = X[i*d+c];
    }
    double alphas[3] = {0.0, 0.5, 1.0};
    for (int a2 = 0; a2 < 3; a2++){
      PenalizedLogit pl;
      pl.SetData(Xnet, yvec);
      pl.SetAlpha(alphas[a2]);
      pl.SetFolds(10);
      pl.SetSeed(20261020UL + 1000UL * (unsigned long) a2);
      double lambda = 0.0;
      double cvA = _NAN_;      // inner-CV AUC: only used to pick lambda, not reported (optimistic)
      pl.CrossValidate(lambda, cvA);
      Vector bet;
      double icpt = 0.0;
      pl.Fit(lambda, bet, icpt);
      for (int j = 0; j < modelCount; j++)
        netCoefAll[j][a2] = bet[j];
      Vector rawNet(N);
      for (int k = 0; k < N; k++){
        double sum = 0.0;
        for (int j = 0; j < modelCount; j++)
          sum += bet[j] * Xnet[k][j];
        rawNet[k] = sum;
      }
      stats->netAucRaw[a2] = LogisticAUC(rawNet, yvec);
      stats->netAuc[a2] = pl.AUC(bet, icpt);
    } // End of a2 loop
  }
  // ---- the model weights ----
  weights.Dimension(modelCount, 5);
  weights.Zero();
  for (int j = 0; j < modelCount; j++){
    if (modelIndex[j] == bestIndex)
      weights[j][0] = 1.0;
    weights[j][1] = wsumW[j];
    for (int a2 = 0; a2 < 3; a2++)
      weights[j][2 + a2] = netCoefAll[j][a2];
  } // End of j loop
  return LOGISTIC_OK;
}






// ---------------------------------------------------------------------
// Cross-validation
// ---------------------------------------------------------------------
// One round: the samples flagged in testFlag form the test set, all the others
// the train set.  For every binary trait the models (BEST, WTSUM, RIDGE, NET,
// LASSO) are built in the train set by Logistic_BuildModel and evaluated in the
// test set by Logistic_ReadModel.
//
//   rawAUC, adjAUC   output, affectionCount x 5 (the columns of the weights
//                    matrix of Logistic_BuildModel), _NAN_ when a model could
//                    not be built or evaluated in this round
static void ScanCV1Time(Matrix &risks, StringArray &rsNames, Matrix &allX,
  Matrix &allY, unsigned long long int **masks, unsigned long long int **LD,
  int nTop, unsigned long long int *testFlag, Matrix &rawAUC, Matrix &adjAUC,
  int coreCount)
{
  const int idCount = allY.rows;
  const int affectionCount = allY.cols;
  const int covariateCount = allX.cols;
  const int twordCount = (affectionCount - 1) / 64 + 1;
  rawAUC.Dimension(affectionCount, 5);
  adjAUC.Dimension(affectionCount, 5);
  rawAUC.Set(_NAN_);
  adjAUC.Set(_NAN_);
  // train / test split: the masks of the samples that are not in the set are cleared
  unsigned long long int **trainMasks = new unsigned long long int * [idCount];
  unsigned long long int **testMasks = new unsigned long long int * [idCount];
  for (int i = 0; i < idCount; i++){
    trainMasks[i] = new unsigned long long int [twordCount];
    testMasks[i] = new unsigned long long int [twordCount];
    bool inTest = (testFlag[i >> 6] >> (i & 0x3F)) & 0x1;
    for (int w = 0; w < twordCount; w++){
      trainMasks[i][w] = inTest ? 0 : masks[i][w];
      testMasks[i][w] = inTest ? masks[i][w] : 0;
    } // End of w loop
  } // End of i loop
  IntArray *trainPheno = new IntArray[affectionCount];
  IntArray *trainYs = new IntArray[affectionCount];
  Vector *trainXs = new Vector[affectionCount];
  IntArray *testPheno = new IntArray[affectionCount];
  IntArray *testYs = new IntArray[affectionCount];
  Vector *testXs = new Vector[affectionCount];
  MakeVariable(allX, allY, trainMasks, trainPheno, trainXs, trainYs);
  MakeVariable(allX, allY, testMasks, testPheno, testXs, testYs);
  for (int a = 0; a < affectionCount; a++){
    Matrix weights;
    IntArray modelIndex;
    int status = Logistic_BuildModel(risks, trainPheno[a], trainYs[a], trainXs[a],
      covariateCount, LD, nTop, coreCount, weights, modelIndex, NULL);
    int modelCount = modelIndex.Length();
    if (status != LOGISTIC_OK || modelCount == 0 || weights.rows != modelCount)
      continue;
    StringArray pgsName(modelCount);
    for (int j = 0; j < modelCount; j++)
      pgsName[j] = rsNames[modelIndex[j]];
    Vector raw, adj, beta, se;
    Logistic_ReadModel(pgsName, weights, rsNames, risks, testPheno[a], testYs[a],
      testXs[a], covariateCount, raw, adj, beta, se);
    for (int m = 0; m < 5; m++){
      rawAUC[a][m] = raw[m];
      adjAUC[a][m] = adj[m];
    } // End of m loop
  } // End of a loop
  for (int i = 0; i < idCount; i++){
    delete []trainMasks[i];
    delete []testMasks[i];
  }
  delete []trainMasks;
  delete []testMasks;
  delete []trainPheno;
  delete []trainYs;
  delete []trainXs;
  delete []testPheno;
  delete []testYs;
  delete []testXs;
}

// Repeated random-split cross-validation: every run draws idCount/fold samples
// as the test set (testFlags[run], one bit per sample, allocated by the caller
// with ceil(idCount/64) words) and calls ScanCV1Time.  rawAUC[run] and
// adjAUC[run] are affectionCount x 5, see ScanCV1Time.
static void ScanCV(Matrix &risks, StringArray &rsNames, Matrix &allX, Matrix &allY,
  unsigned long long int **masks, unsigned long long int **LD, int fold, int repeat,
  int nTop, Matrix *rawAUC, Matrix *adjAUC, unsigned long long int **testFlags,
  int coreCount)
{
  const int idCount = allY.rows;
  Random ran;
  IntArray testFlag(idCount);
  int imemsize = ((idCount - 1) / 64 + 1) * sizeof(unsigned long long int);
  int testsetCount = idCount * 1.0 / fold;
  for (int run = 0; run < repeat; run++){
    ran.Choose(testFlag, idCount, testsetCount);
    memset(&testFlags[run][0], 0, imemsize);
    for (int i = 0; i < idCount; i += 64){
      int iMax = i < idCount - 64 ? i + 64 : idCount;
      int w = (i >> 6);
      unsigned long long int word0 = 0x1;
      for (int ii = i; ii < iMax; ii++, word0 <<= 1)
        if (testFlag[ii]) testFlags[run][w] |= word0;
    } // End of i loop
  } // End of run loop
  for (int run = 0; run < repeat; run++){
    printf(".");
    fflush(stdout);
    ScanCV1Time(risks, rsNames, allX, allY, masks, LD, nTop, testFlags[run],
      rawAUC[run], adjAUC[run], coreCount);
  } // End of run loop
}

// ---------------------------------------------------------------------
// Logistic_Scan
// ---------------------------------------------------------------------
void Logistic_Scan(const char * phefile, const char * dictfile,
   const char * rsfile, const char * prefix, const char * covariateList,
   int nTop, double ldT, const char * modelfile, bool rplotFlag,
   bool cvFlag, int fold, int repeat, int coreCount){
  printf("\nOptions in effect:\n");
  printf("\t--logistic\n");
  if (modelfile[0])
    printf("\t--readmodel %s\n", modelfile);
  printf("\t--pheno %s\n", phefile);
  if (gCovarFile[0]) printf("\t--covarfile %s\n", gCovarFile);
  if (gTraitList[0]) printf("\t--trait %s\n", gTraitList);
  if (covariateList[0])
    printf("\t--covariate %s\n", covariateList);
  printf("\t--score %s\n", rsfile);
  if (dictfile[0])
    printf("\t--diction %s\n", dictfile);
  if (nTop)
    printf("\t--ntop %d\n", nTop);
  if (ldT != _NAN_ && ldT < 1 && ldT >= 0)
    printf("\t--pruner2 %.3lf\n", ldT);
  else
    ldT = _NAN_;
  if (rplotFlag) printf("\t--rplot\n");
  if (cvFlag) printf("\t--cv\n");
  if (fold) printf("\t--fold %d\n", fold);
  if (repeat) printf("\t--repeat %d\n", repeat);
  if (coreCount)
    printf("\t--cpus %d\n", coreCount);
  printf("\t--prefix %s\n", prefix);
  printf("\n");
  if (cvFlag){
    if (fold == 0) fold = 10;
    if (repeat == 0) repeat = 10;
    if (fold < 2) error("--fold must be at least 2.");
    if (repeat < 1) error("--repeat must be at least 1.");
  }

  StringArray sampleID, alltraitNames, rsNames, annotation;
  IntArray covariates, affections, quantraits;
  Matrix alltraits, risks;
  SECRET_Read(phefile, dictfile, rsfile, "", prefix, covariateList,
    sampleID, alltraitNames, alltraits, covariates, affections, quantraits,
    rsNames, risks, annotation, coreCount);
  int rsCount = rsNames.Length();
  int idCount = sampleID.Length();
  if (idCount == 0) return;
  int quantraitCount = quantraits.Length();
  int affectionCount = affections.Length();
  if (affectionCount == 0){
    printf("No dichotomous traits were found. --logistic needs a trait "
         "coded 0/1 or 1/2.\n");
    return;
  }
  // phenotypes: binary traits only, so allY columns are the affections
  IntArray noQuantraits(0);
  Matrix allY(idCount, affectionCount);
  int twordCount = (affectionCount - 1) / 64 + 1;
  unsigned long long int ** masks = new unsigned long long int * [idCount];
  for (int i = 0; i < idCount; i++)
    masks[i] = new unsigned long long int [twordCount];
  Matrix allX(idCount, covariates.Length());
  MakePhenos(risks, alltraits, alltraitNames, affections, noQuantraits,
    covariates, allY, allX, masks, coreCount);
  int covariateCount = covariates.Length();
  int pcols = 1 + covariateCount + 1;      // intercept, covariates, PGS
  IntArray *pheno = new IntArray[affectionCount];
  IntArray *Ys = new IntArray[affectionCount];
  Vector *Xs = new Vector[affectionCount];
  MakeVariable(allX, allY, masks, pheno, Xs, Ys);
  String msg(256);
  if(modelfile[0]){
    Matrix coef;
    StringArray methodName, pgsName;
    ReadMatrix(modelfile, methodName, pgsName, coef, coreCount);
    int modelrsCount = coef.rows;
    int methodCount = coef.cols;
    if (modelrsCount < 1 || methodCount < 1)
      error("Model file %s does not contain any coefficients.", modelfile);
    Matrix rawAUC(affectionCount, methodCount);
    Matrix adjAUC(affectionCount, methodCount);
    Matrix beta(affectionCount, methodCount);
    Matrix se(affectionCount, methodCount);
    if (cvFlag)
      printf("--cv is ignored when --readmodel is given\n");
    for (int a = 0; a < affectionCount; a++){
      Vector raw, adj, b, sd;
      int missingTerms = Logistic_ReadModel(pgsName, coef, rsNames, risks,
        pheno[a], Ys[a], Xs[a], covariateCount, raw, adj, b, sd);
      if (a == 0 && missingTerms)
        printf("  %d model term(s) not present in the score file were skipped\n",
          missingTerms);
      for (int m = 0; m < methodCount; m++){
        rawAUC[a][m] = raw[m];
        adjAUC[a][m] = adj[m];
        beta[a][m] = b[m];
        se[a][m] = sd[m];
      } // End of m loop
    } // End of a loop
    printf("\nSummary of prediction from %s\n", modelfile);
    printf("%20s%12s%8s%10s%8s%10s%8s%8s\n",
      "TRAIT", "MODEL", "N", "OR", "Z", "P", "AUC", "AUC_adj");
    for (int a = 0; a < affections.Length(); a++){
      int trait = affections[a];
      for(int m = 0; m < methodCount; m++){
        if (beta[a][m] == _NAN_ || se[a][m] == _NAN_){
          printf("%20s%12s%8d%10s%8s%10s%8.3lf%8s\n",
            (const char*)alltraitNames[trait], (const char*)methodName[m],
            pheno[a].Length(), "NA", "NA", "NA",
            rawAUC[a][m] == _NAN_ ? 0.0 : rawAUC[a][m], "NA");
          continue;
        }
        double z = beta[a][m] / se[a][m];
        double p = 2.0 * ndist(fabs(z));
        printf("%20s%12s%8d%10.3lf%8.2lf%10.2G%8.3lf%8.3lf\n",
          (const char*)alltraitNames[trait], (const char*)methodName[m],
          pheno[a].Length(), exp(beta[a][m]), z, p, rawAUC[a][m], adjAUC[a][m]);
      } // End of m loop
    } // End of a loop
    if (rplotFlag){
      // combined score of every model for every trait, then the calibration plots
      IntArray rsIndex(modelrsCount);
      for (int r = 0; r < modelrsCount; r++)
        rsIndex[r] = rsNames.Find(pgsName[r]);
      for (int a = 0; a < affectionCount; a++){
        int trait = affections[a];
        int N = pheno[a].Length();
        int cases = 0;
        for (int k = 0; k < N; k++)
          if (Ys[a][k] == 1) cases++;
        if (N < 10 || cases < 1 || cases == N){
          printf("  Calibration plots for trait %s skipped: %d samples, %d cases\n",
            (const char*)alltraitNames[trait], N, cases);
          continue;
        }
        Matrix scoreTable(N, methodCount + 1);
        for (int k = 0; k < N; k++)
          scoreTable[k][0] = Ys[a][k];
        for (int m = 0; m < methodCount; m++){
          Vector score(N);
          score.Zero();
          bool used = false;
          for (int r = 0; r < modelrsCount; r++){
            int rr = rsIndex[r];
            double bj = coef[r][m];
            if (rr == -1 || bj == _NAN_) continue;
            used = true;
            for (int k = 0; k < N; k++)
              score[k] += risks[rr][pheno[a][k]] * bj;
          } // End of r loop
          for (int k = 0; k < N; k++)
            scoreTable[k][m + 1] = used ? score[k] : _NAN_;
        } // End of m loop
        StringArray sheader(methodCount + 2), scol0(N);
        sheader[0] = "ID";
        sheader[1] = alltraitNames[trait];
        for (int m = 0; m < methodCount; m++)
          sheader[m + 2] = methodName[m];
        for (int k = 0; k < N; k++)
          scol0[k] = sampleID[pheno[a][k]];
        String outfile0(prefix);
        outfile0.Add("_");
        outfile0.Add(alltraitNames[trait]);
        String pgsfile = outfile0;
        pgsfile.Add("_combinedPGS.txt");
        WriteTable(pgsfile, scoreTable, sheader, scol0);
        printf("  %d polygenic scores for trait %s write to %s\n", methodCount,
          (const char*)alltraitNames[trait], (const char*)pgsfile);
        plotCalibrationLogistic(outfile0, alltraitNames[trait], "", msg);
        printf("  %s", (const char*)msg);
      } // End of a loop
    }
    for (int i = 0; i < idCount; i++)
      delete []masks[i];
    delete []masks;
    delete []pheno;
    delete []Ys;
    delete []Xs;
    return;
  }
  if(ldT == _NAN_){
    ldT = 0.2;
    printf("R2 to prune correlated PGSs for WTSUM is missing. "
      "--pruner2=%.1lf is set\n", ldT);
  }
  if (!nTop){
    nTop = rsCount < 30 ? rsCount : 30;
    printf("The number of PGSs for a combined PGS (ntop) is missing. "
      "--ntop=%d is set\n", nTop);
  }
  int annotationCount = annotation.Length();
  int wordCount = (rsCount - 1) / 64 + 1;
  unsigned long long int ** LD = new unsigned long long int * [rsCount];
  for (int r = 0; r < rsCount; r++){
    LD[r] = new unsigned long long int [wordCount];
    memset(&LD[r][0], 0, wordCount * sizeof(unsigned long long int));
  }
  Vector rsSD(rsCount);
  time_t timer;
  struct tm * tblock;
  timer = time(NULL);
  tblock = localtime(&timer);
  printf("Computing PGS correlations starts at %s", asctime(tblock));
  fflush(stdout);
  MakeBinaryLD(risks, ldT, LD, rsSD, coreCount);
  LDSummary(LD, rsCount, wordCount, coreCount);
  timer = time(NULL);
  tblock = localtime(&timer);
  printf("Computing PGS correlations ends at %s", asctime(tblock));
  fflush(stdout);
  printf("\nLogistic regression of %d dichotomous trait(s) on %d PGS(s)",
    affectionCount, rsCount);
  if(covariateCount)
    printf(", adjusted for %d covariate(s)", covariateCount);
  printf("\n");
  Vector bestBeta(affectionCount), bestZ(affectionCount), bestP(affectionCount);
  Vector wsBeta(affectionCount), wsZ(affectionCount), wsP(affectionCount);
  IntArray bestIndex(affectionCount), nCase(affectionCount), nControl(affectionCount);
  Vector aucBest(affectionCount), aucWsum(affectionCount);
  Vector aucBestRaw(affectionCount), aucWsumRaw(affectionCount);
  Matrix netAuc(affectionCount, 3);
  Matrix netAucRaw(affectionCount, 3);
  for (int a = 0; a < affectionCount; a++){
    bestBeta[a] = bestZ[a] = bestP[a] = _NAN_;
    wsBeta[a] = wsZ[a] = wsP[a] = _NAN_;
    aucBest[a] = aucWsum[a] = _NAN_;
    aucBestRaw[a] = aucWsumRaw[a] = _NAN_;
    bestIndex[a] = -1;
    nCase[a] = nControl[a] = 0;
    for (int m = 0; m < 3; m++){
      netAuc[a][m] = _NAN_;
      netAucRaw[a][m] = _NAN_;
    } // End of m loop
  } // End of a loop
  for (int a = 0; a < affectionCount; a++){
    int trait = affections[a];
    printf("\nTrait %s...\n", (const char*)alltraitNames[trait]);
    int N = pheno[a].Length();
    Matrix weights;
    IntArray modelIndex;
    LogisticModelStats st;
    int status = Logistic_BuildModel(risks, pheno[a], Ys[a], Xs[a],
      covariateCount, LD, nTop, coreCount, weights, modelIndex, &st);
    nCase[a] = st.nCase;
    nControl[a] = st.nControl;
    if (status == LOGISTIC_TOO_FEW_SAMPLES){
      printf("  too few samples with phenotype and covariates (N=%d), skipped\n", N);
      continue;
    }
    if (status == LOGISTIC_TOO_FEW_FOR_TERMS){
      printf("  too few samples for %d model terms (N=%d), skipped\n", pcols, N);
      continue;
    }
    if (status == LOGISTIC_TOO_FEW_CASES_CONTROLS){
      printf("  too few cases (%d) or controls (%d), skipped\n", st.nCase, st.nControl);
      continue;
    }
    int cases = st.nCase, controls = st.nControl;
    printf("  %d samples: %d cases, %d controls\n", N, cases, controls);
    Vector &beta = st.beta, &se = st.se, &zval = st.zval;
    IntArray &sortedIndex = st.sortedIndex;
    int sortedCount = sortedIndex.Length();
    String outfile0(prefix);
    outfile0.Add("_");
    outfile0.Add(alltraitNames[trait]);
    String outfile = outfile0;
    outfile.Add(".logistic");
    FILE * fp = fopen(outfile, "wb");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)outfile);
    fprintf(fp, "PGS\tN\tNCase\tNControl\tBeta\tSE\tOR\tZ\tPvalue");
    WriteCoverageHeader(fp);
    if (annotationCount)
      fprintf(fp, "\tAnnotation");
    for (int r = 0; r < sortedCount; r++){
      int rs = sortedIndex[r];
      fprintf(fp, "\n%s\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.3lf\t%.2G",
        (const char*)rsNames[rs], N, cases, controls,
        beta[rs], se[rs], exp(beta[rs]), zval[rs], 2.0 * ndist(fabs(zval[rs])));
      WriteCoverageFields(fp, rs);
      if (annotationCount)
        fprintf(fp, "\t%s", (const char*)annotation[rs]);
    } // End of r loop
    fprintf(fp, "\n");
    fclose(fp);
    printf("  PGS scan statistics write to %s\n", (const char*)outfile);
    {
      // correlations among the top PGSs: all samples, cases only, controls only, and the
      // test of the difference between cases and controls
      int topCount = nTop < sortedCount ? nTop : sortedCount;
      IntArray topIndex(topCount), corrStatus(N);
      for (int r = 0; r < topCount; r++)
        topIndex[r] = sortedIndex[r];
      for (int k = 0; k < N; k++)
        corrStatus[k] = (Ys[a][k] == 1) ? 1 : 0;
      bool caseControlDone = false;
      WritePGSCorrelations(outfile0, risks, rsNames, topIndex, pheno[a], corrStatus, coreCount, caseControlDone, ldT);
      if (rplotFlag){
        plotPGSCorr(outfile0, outfile, alltraitNames[trait], "", msg);
        printf("  %s", (const char*)msg);
      }
    }
    if(rplotFlag){
      // combined score of every model for the plots: BEST, WTSUM, RIDGE, NET, LASSO
      // (the columns of weights), expressed against the same risks the models were
      // built on; the outcome is coded 1 = affected, 0 = unaffected.  A model that
      // could not be fitted (all weights _NAN_) is written as NA.
      int wCount = modelIndex.Length();
      if (wCount > 0 && weights.rows == wCount){
        Matrix scoreTable(N, 6);
        for (int k = 0; k < N; k++)
          scoreTable[k][0] = Ys[a][k];
        for (int m = 0; m < 5; m++){
          bool fitted = false;
          for (int j = 0; j < wCount; j++)
            if (weights[j][m] != _NAN_) fitted = true;
          for (int k = 0; k < N; k++){
            double sum = 0.0;
            for (int j = 0; j < wCount; j++)
              if (weights[j][m] != _NAN_)
                sum += weights[j][m] * risks[modelIndex[j]][pheno[a][k]];
            scoreTable[k][m + 1] = fitted ? sum : _NAN_;
          } // End of k loop
        } // End of m loop
        StringArray sheader(7), scol0(N);
        sheader[0] = "ID";
        sheader[1] = alltraitNames[trait];
        sheader[2] = "BEST";
        sheader[3] = "WTSUM";
        sheader[4] = "RIDGE";
        sheader[5] = "NET";
        sheader[6] = "LASSO";
        for (int k = 0; k < N; k++)
          scol0[k] = sampleID[pheno[a][k]];
        String pgsfile = outfile0;
        pgsfile.Add("_PGS.txt");
        WriteTable(pgsfile, scoreTable, sheader, scol0);
        printf("  Combined PGSs of %d model(s) for trait %s write to %s\n", 5,
          (const char*)alltraitNames[trait], (const char*)pgsfile);
      }
      plotLogistic(outfile0, alltraitNames[trait], "", msg);
      printf("  %s", (const char*)msg);
    }
    // best PGS
    int best = st.bestIndex;
    if (best >= 0){
      bestIndex[a] = best;
      bestBeta[a] = beta[best];
      bestZ[a] = zval[best];
      bestP[a] = 2.0 * ndist(fabs(zval[best]));
      printf("  Top PGS %s: OR=%.3lf Z=%.3f P=%.2G\n",
        (const char*)rsNames[best], exp(beta[best]), zval[best], bestP[a]);
    }
    int modelCount = modelIndex.Length();
    aucBest[a] = st.aucBest;
    aucBestRaw[a] = st.aucBestRaw;
    if (st.wsumStatus == 1){
      wsBeta[a] = st.wsBeta;
      wsZ[a] = st.wsZ;
      wsP[a] = st.wsP;
      aucWsum[a] = st.aucWsum;
      aucWsumRaw[a] = st.aucWsumRaw;
      printf("  WTSUM of %d PGSs: OR=%.3lf Z=%.3f P=%.2G AUC=%.3f\n",
        st.wsumTerms, exp(wsBeta[a]), wsZ[a], wsP[a], aucWsum[a]);
    }else if (st.wsumStatus == 2)
      printf("  WTSUM of %d PGSs: model did not converge\n", st.wsumTerms);
    if (modelCount > 0)
      printf("  %d PGSs in the models (--ntop %d); WTSUM uses %d of them after pruning at R2 > %.2f (--pruner2)\n",
        modelCount, nTop, st.wsumTerms, ldT);
    if (modelCount > 0){
      if (st.penalizedSkipped)
        printf("  penalized models skipped: %d model terms for %d samples\n",
          st.penalizedTerms, N);
      else{
        const char * netName[3] = {"ridge", "elastic net", "lasso"};
        for (int a2 = 0; a2 < 3; a2++){
          netAucRaw[a][a2] = st.netAucRaw[a2];
          netAuc[a][a2] = st.netAuc[a2];
          printf("  %-12s top-%d model: AUC=%.3f\n",
                  netName[a2], modelCount, netAuc[a][a2]);
        } // End of a2 loop
      }
    }
    // ---- write the model file consumed by --readmodel ----
    if (modelCount > 0){
      Matrix fileWeights(modelCount, 5);
      StringArray mheader(6);
      mheader[0] = "PGS_NAME";
      mheader[1] = "BEST";
      mheader[2] = "WTSUM";
      mheader[3] = "RIDGE";
      mheader[4] = "NET";
      mheader[5] = "LASSO";
      StringArray mcol0(modelCount);
      for (int j = 0; j < modelCount; j++){
        mcol0[j] = rsNames[modelIndex[j]];
        // The scan re-standardises risks in place by rsSD (ScreenLD2), but
        // --readmodel returns before that step.  Dividing the fitted
        // coefficients by rsSD expresses them against the risks matrix
        // exactly as --readmodel will see it, so the reconstruction is
        // identical.  BEST stays 1.0, matching --linear.
        double sdj = rsSD[modelIndex[j]];
        if (sdj < 1E-100)
          sdj = 1.0;
        fileWeights[j][0] = weights[j][0];
        for (int c = 1; c < 5; c++)
          fileWeights[j][c] = (weights[j][c] == _NAN_) ? _NAN_ : weights[j][c] / sdj;
      } // End of j loop
      String mfile(prefix);
      mfile.Add("_");
      mfile.Add(alltraitNames[trait]);
      mfile.Add("_PGS.mod");
      WriteTable(mfile, fileWeights, mheader, mcol0);
      printf("  PGS model coefficients for trait %s write to %s\n",
        (const char*)alltraitNames[trait], (const char*)mfile);
    }
  } // End of trait loop
  printf("\nSummary of logistic regression for %d trait(s)\n", affectionCount);
  printf("%20s%10s%10s%10s%24s%10s",
    "TRAIT", "N", "NCase", "NControl", "BEST_PGS", "Best_P");
  if (annotationCount) printf("%11s", "Annotation");
  printf("\n");
  for (int a = 0; a < affectionCount; a++){
    int trait = affections[a];
    const char * bestName = bestIndex[a] >= 0 ?
         (const char*)rsNames[bestIndex[a]] : "NA";
    printf("%20s%10d%10d%10d%24s%10.2G\n",
         (const char*)alltraitNames[trait], nCase[a] + nControl[a],
         nCase[a], nControl[a], bestName,
         bestP[a] == _NAN_ ? 1.0 : bestP[a]);
    if (annotationCount)
      printf("%11s", (const char*)annotation[bestIndex[a]]);
  } // End of a loop

  // AUC summaries (0.000 means the model was not fitted): the unadjusted
  // table ranks on the score alone, the adjusted one adds the covariates
  printf("\nSummary of %d Trait Prediction AUC (unadjusted)\n", affectionCount);
  printf("%20s%8s%10s%10s%10s%10s%10s\n",
     "TRAIT", "N", "AUC_Best", "AUC_WTSUM", "AUC_Ridge", "AUC_Net", "AUC_Lasso");
  for (int a = 0; a < affectionCount; a++){
    int trait = affections[a];
    printf("%20s%8d%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
      (const char*)alltraitNames[trait], nCase[a] + nControl[a],
      aucBestRaw[a] == _NAN_ ? 0.0 : aucBestRaw[a],
      aucWsumRaw[a] == _NAN_ ? 0.0 : aucWsumRaw[a],
      netAucRaw[a][0] == _NAN_ ? 0.0 : netAucRaw[a][0],
      netAucRaw[a][1] == _NAN_ ? 0.0 : netAucRaw[a][1],
      netAucRaw[a][2] == _NAN_ ? 0.0 : netAucRaw[a][2]);
  } // End of a loop
  printf("\nSummary of %d Trait Prediction AUC (adjusted for covariates)\n", affectionCount);
  printf("%20s%8s%10s%10s%10s%10s%10s\n",
    "TRAIT", "N", "AUC_Best", "AUC_WTSUM", "AUC_Ridge", "AUC_Net", "AUC_Lasso");
  for (int a = 0; a < affectionCount; a++){
    int trait = affections[a];
    printf("%20s%8d%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n",
      (const char*)alltraitNames[trait], nCase[a] + nControl[a],
      aucBest[a] == _NAN_ ? 0.0 : aucBest[a],
      aucWsum[a] == _NAN_ ? 0.0 : aucWsum[a],
      netAuc[a][0] == _NAN_ ? 0.0 : netAuc[a][0],
      netAuc[a][1] == _NAN_ ? 0.0 : netAuc[a][1],
      netAuc[a][2] == _NAN_ ? 0.0 : netAuc[a][2]);
  } // End of a loop
  // ---- cross-validation: build the models in a train set, test them in a test set ----
  if (cvFlag){
    int iwordCount = (idCount - 1) / 64 + 1;
    unsigned long long int ** testFlags = new unsigned long long int * [repeat];
    for (int run = 0; run < repeat; run++)
      testFlags[run] = new unsigned long long int [iwordCount];
    Matrix * cvRaw = new Matrix[repeat];
    Matrix * cvAdj = new Matrix[repeat];
    printf("\nCross-validation (test set = 1/%d of the samples, %d round(s)) start",
      fold, repeat);
    ScanCV(risks, rsNames, allX, allY, masks, LD, fold, repeat, nTop,
      cvRaw, cvAdj, testFlags, coreCount);
    printf("done.\n");
    for (int pass = 0; pass < 2; pass++){
      Matrix * res = pass ? cvAdj : cvRaw;
      printf("\nSummary of %d Trait Cross-validation AUC (%s)\n", affectionCount,
        pass ? "adjusted for covariates" : "unadjusted");
      printf("%20s%8s%8s%10s%10s%10s%10s%10s\n", "TRAIT", "N", "Rounds",
        "AUC_Best", "AUC_WTSUM", "AUC_Ridge", "AUC_Net", "AUC_Lasso");
      for (int a = 0; a < affectionCount; a++){
        int trait = affections[a];
        double sum[5] = {0, 0, 0, 0, 0};
        int cnt[5] = {0, 0, 0, 0, 0};
        for (int run = 0; run < repeat; run++)
          for (int m = 0; m < 5; m++)
            if (res[run][a][m] != _NAN_){
              sum[m] += res[run][a][m];
              cnt[m]++;
            }
        printf("%20s%8d%8d", (const char*)alltraitNames[trait],
          nCase[a] + nControl[a], cnt[0]);
        for (int m = 0; m < 5; m++){
          if (cnt[m])
            printf("%10.3lf", sum[m] / cnt[m]);
          else
            printf("%10s", "NA");
        } // End of m loop
        printf("\n");
      } // End of a loop
    } // End of pass loop
    for (int run = 0; run < repeat; run++)
      delete []testFlags[run];
    delete []testFlags;
    delete []cvRaw;
    delete []cvAdj;
  }
  for (int i = 0; i < idCount; i++)
    delete []masks[i];
  delete []masks;
  delete []pheno;
  delete []Ys;
  delete []Xs;
  for (int r = 0; r < rsCount; r++)
    delete []LD[r];
  delete []LD;
}



      /*
    for (int r = 0; r < rsCount; r++){
      int rs = sortedIndex[r];
      if (beta[rs] == _NAN_)
        fprintf(fp, "\n%s\t%d\t%d\t%d\tNA\tNA\tNA\tNA\tNA",
          (const char*)rsNames[rs], N, cases, controls);
      else
        fprintf(fp, "\n%s\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.3lf\t%.2G",
          (const char*)rsNames[rs], N, cases, controls,
          beta[rs], se[rs], exp(beta[rs]), zval[rs], 2.0 * ndist(fabs(zval[rs])));
      if (annotationCount)
        fprintf(fp, "\t%s", (const char*)annotation[rs]);
    } // End of r loop
*/

