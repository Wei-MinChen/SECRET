////////////////////////////////////////////////////////////////////////
// survival.cpp
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
// Cox proportional hazards scan (--coxph --surv time,event).
//
// One Cox model is fitted per polygenic score:
//     Surv(time, event) ~ PGS + covariates
// and the log hazard ratio, its standard error, the hazard ratio,
// Wald z and two-sided p are reported, exactly the columns the previous
// R/coxph script produced, now in <prefix>_<time>2<event>.coxph (like <prefix>_<trait>.linear).
//
// The fit is Newton-Raphson on the Efron approximation to the tied
// partial likelihood (R's coxph default), so ties in the survival times
// are handled the same way. R is no longer needed for the statistics;
// with --rplot it is still used to draw the hazard-ratio plot from the
// table this file writes.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "MathSVD.h"
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#include "Random.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

extern int CheckRout(const char *scriptfile);
extern double GenomicControlNumber(Vector &Z, IntArray &sortedrsIndex);
extern bool gRecodeBinary01;
extern void WriteTopHeatmapR(FILE *fp, const char *corrfile, const char *statfile, int kind, const char *traitName);
extern void WritePGSLabelR(FILE *fp);
extern void MakeBinaryLD(Matrix & risks, double ldT, unsigned long long int ** LD, Vector & rsSD, int coreCount);
extern void MakePrunedOrder(unsigned long long int ** LD, IntArray & order);
extern void LDSummary(unsigned long long int ** LD, int rsCount, int wordCount, int coreCount);
extern int plotCalibrationCox(const char *input, const char *timeName, const char *eventName, const char *rpath, String &msg);
extern int plotPGSCorr(const char *input, const char *statfile, const char *traitName, const char *rpath, String &msg);
extern String RCommandLine(const char *rpath, const char *scriptfile);
extern const char *MakePdf(const char *stem);
extern void WritePGSCorrelations(const char *outfile0, Matrix &risks, StringArray &rsNames, IntArray &topIndex,
   IntArray &samples, IntArray &status, int coreCount, bool &caseControlDone, double pruneR2);
extern void WriteTable(String & outfile, Matrix & table, StringArray & header,
   StringArray & col0);

// ---------------------------------------------------------------------
// dense symmetric positive definite solve (in-place Cholesky)
// ---------------------------------------------------------------------
static bool coxCholesky(double ** A, double * b, int n)
   {
   for (int i = 0; i < n; i++)
      {
      for (int j = 0; j <= i; j++)
         {
         double sum = A[i][j];
         for (int k = 0; k < j; k++)
            sum -= A[i][k] * A[j][k];
         if (i == j)
            {
            if (sum <= 0.0)
               return false;
            A[i][i] = sqrt(sum);
            }
         else
            A[i][j] = sum / A[j][j];
         }
      }
   for (int i = 0; i < n; i++)
      {
      double sum = b[i];
      for (int k = 0; k < i; k++)
         sum -= A[i][k] * b[k];
      b[i] = sum / A[i][i];
      }
   for (int i = n - 1; i >= 0; i--)
      {
      double sum = b[i];
      for (int k = i + 1; k < n; k++)
         sum -= A[k][i] * b[k];
      b[i] = sum / A[i][i];
      }
   return true;
   }

// ---------------------------------------------------------------------
// Efron partial likelihood and its first two derivatives at beta.
//
// X rows are sorted by DECREASING survival time, so the risk set is a
// prefix of the data.  ts and ev are the time and event columns in that
// same order.  When score/info are NULL only the log-likelihood is
// accumulated (used for step halving).
//
// sid (NULL = no strata) is the stratum of each row.  Rows must be sorted by
// stratum first and by decreasing time within a stratum (see coxOrder).  Each
// stratum has its own baseline hazard: its risk sets are formed from its own
// rows only, and the log-likelihood, score and information are summed over
// the strata, exactly as R's coxph(... + strata(s)).
// ---------------------------------------------------------------------
static double coxEfron(int N, int p, double ** X, double * ts, int * ev, const int * sid,
   const double * beta, double * score, double ** info)
   {
   double * eta = new double [N];
   double * wt = new double [N];
   double emax = -1E300;
   for (int i = 0; i < N; i++)
      {
      double e = 0.0;
      for (int j = 0; j < p; j++)
         e += X[i][j] * beta[j];
      eta[i] = e;
      if (e > emax)
         emax = e;
      }
   // shifting eta by a constant leaves ratios and derivatives unchanged
   for (int i = 0; i < N; i++)
      wt[i] = exp(eta[i] - emax);

   double R0 = 0.0;
   double * R1 = new double [p];
   double ** R2 = new double * [p];
   for (int a = 0; a < p; a++)
      {
      R1[a] = 0.0;
      R2[a] = new double [p];
      for (int b = 0; b < p; b++)
         R2[a][b] = 0.0;
      }
   double * D1 = new double [p];
   double ** D2 = new double * [p];
   for (int a = 0; a < p; a++)
      {
      D1[a] = 0.0;
      D2[a] = new double [p];
      for (int b = 0; b < p; b++)
         D2[a][b] = 0.0;
      }
   double * den1 = new double [p];

   if (score)
      for (int a = 0; a < p; a++)
         score[a] = 0.0;
   if (info)
      for (int a = 0; a < p; a++)
         for (int b = 0; b < p; b++)
            info[a][b] = 0.0;

   double ll = 0.0;
   int i = 0;
   while (i < N)
      {
      if (sid && i > 0 && sid[i] != sid[i - 1])
         {
         // a new stratum: its risk sets start from nothing
         R0 = 0.0;
         for (int a = 0; a < p; a++)
            {
            R1[a] = 0.0;
            for (int b = 0; b < p; b++)
               R2[a][b] = 0.0;
            }
         }
      int j = i;
      while (j < N && ts[j] == ts[i] && (!sid || sid[j] == sid[i]))
         j++;

      // add the whole tie group to the risk set
      for (int k = i; k < j; k++)
         {
         R0 += wt[k];
         for (int a = 0; a < p; a++)
            {
            R1[a] += wt[k] * X[k][a];
            for (int b = 0; b <= a; b++)
               R2[a][b] += wt[k] * X[k][a] * X[k][b];
            }
         }

      int m = 0;
      for (int k = i; k < j; k++)
         if (ev[k])
            m++;
      if (m > 0)
         {
         double D0 = 0.0;
         for (int a = 0; a < p; a++)
            {
            D1[a] = 0.0;
            for (int b = 0; b < p; b++)
               D2[a][b] = 0.0;
            }
         for (int k = i; k < j; k++)
            if (ev[k])
               {
               D0 += wt[k];
               for (int a = 0; a < p; a++)
                  {
                  D1[a] += wt[k] * X[k][a];
                  for (int b = 0; b <= a; b++)
                     D2[a][b] += wt[k] * X[k][a] * X[k][b];
                  }
               }

         for (int l = 0; l < m; l++)
            {
            double frac = (double) l / (double) m;
            double den0 = R0 - frac * D0;
            if (den0 <= 0.0)
               den0 = 1E-300;
            for (int a = 0; a < p; a++)
               den1[a] = (R1[a] - frac * D1[a]) / den0;
            ll -= log(den0);
            if (score)
               for (int a = 0; a < p; a++)
                  score[a] -= den1[a];
            if (info)
               for (int a = 0; a < p; a++)
                  {
                  for (int b = 0; b <= a; b++)
                     info[a][b] += (R2[a][b] - frac * D2[a][b]) / den0
                        - den1[a] * den1[b];
                  }
            } // End of l loop

         for (int k = i; k < j; k++)
            if (ev[k])
               {
               ll += eta[k] - emax;
               if (score)
                  for (int a = 0; a < p; a++)
                     score[a] += X[k][a];
               }
         }
      i = j;
      } // End of while

   if (info)
      for (int a = 0; a < p; a++)
         for (int b = 0; b < a; b++)
            info[b][a] = info[a][b];

   for (int a = 0; a < p; a++)
      {
      delete []R2[a];
      delete []D2[a];
      }
   delete []R2;
   delete []D2;
   delete []R1;
   delete []D1;
   delete []den1;
   delete []eta;
   delete []wt;
   return ll;
   }

// ---------------------------------------------------------------------
// Newton-Raphson Cox fit.  X must already be sorted by decreasing time.
// Returns the number of iterations, or -1 if the information matrix
// could not be inverted.
// ---------------------------------------------------------------------
static int CoxFit(int N, int p, double ** X, double * ts, int * ev, const int * sid,
   double * beta, double * se, double tol, int maxIter)
   {
   for (int j = 0; j < p; j++)
      {
      beta[j] = 0.0;
      se[j] = 0.0;
      }

   double ** info = new double * [p];
   double ** A = new double * [p];
   for (int a = 0; a < p; a++)
      {
      info[a] = new double [p];
      A[a] = new double [p];
      }
   double * score = new double [p];
   double * step = new double [p];
   double * cand = new double [p];

   int iter = 0;
   double ll = coxEfron(N, p, X, ts, ev, sid, beta, score, info);
   for (iter = 1; iter <= maxIter; iter++)
      {
      if (iter > 1)
         ll = coxEfron(N, p, X, ts, ev, sid, beta, score, info);

      for (int a = 0; a < p; a++)
         {
         for (int b = 0; b < p; b++)
            A[a][b] = info[a][b];
         step[a] = score[a];
         }
      if (!coxCholesky(A, step, p))
         {
         for (int a = 0; a < p; a++)
            {
            for (int b = 0; b < p; b++)
               A[a][b] = info[a][b];
            A[a][a] += 1E-8 * (1.0 + info[a][a]);
            step[a] = score[a];
            }
         if (!coxCholesky(A, step, p))
            {
            for (int a = 0; a < p; a++)
               {
               delete []info[a];
               delete []A[a];
               }
            delete []info;
            delete []A;
            delete []score;
            delete []step;
            delete []cand;
            return -1;
            }
         }

      double alpha = 1.0;
      double llNew = ll;
      bool accepted = false;
      for (int half = 0; half < 30; half++)
         {
         for (int a = 0; a < p; a++)
            cand[a] = beta[a] + alpha * step[a];
         llNew = coxEfron(N, p, X, ts, ev, sid, cand, NULL, NULL);
         if (llNew >= ll - 1E-12)
            {
            accepted = true;
            break;
            }
         alpha *= 0.5;
         }
      if (!accepted)
         break;

      double maxChange = 0.0;
      for (int a = 0; a < p; a++)
         {
         double d = alpha * step[a];
         beta[a] += d;
         if (fabs(d) > maxChange)
            maxChange = fabs(d);
         }
      if (maxChange < tol)
         break;
      }

   // standard errors from the observed information at the estimate
   coxEfron(N, p, X, ts, ev, sid, beta, score, info);
   for (int a = 0; a < p; a++)
      {
      // coxCholesky overwrites A, so restore the whole matrix (not just row a)
      for (int c = 0; c < p; c++)
         for (int b = 0; b < p; b++)
            A[c][b] = info[c][b];
      for (int b = 0; b < p; b++)
         step[b] = (b == a) ? 1.0 : 0.0;
      if (coxCholesky(A, step, p) && step[a] > 0.0)
         se[a] = sqrt(step[a]);
      else
         se[a] = 0.0;
      }

   for (int a = 0; a < p; a++)
      {
      delete []info[a];
      delete []A[a];
      }
   delete []info;
   delete []A;
   delete []score;
   delete []step;
   delete []cand;
   return iter;
   }

// ---------------------------------------------------------------------
// Penalized Cox models (--coxph)
//
// Maximises the penalized Efron partial log-likelihood
//     l(beta) - lambda [ alpha ||beta||_1 + (1-alpha)/2 ||beta||_2^2 ]
// by Newton's method on the quadratic approximation
//     1/2 b'Hb - c'b,   c = H beta + U
// where U and H are the gradient and observed information from coxEfron.
// The quadratic subproblem is solved by Cholesky (ridge) or coordinate
// descent (elastic net / lasso).
//
// `pen[j]` is a penalty mask: a zero entry (used for the covariates) leaves
// that coefficient unpenalized, so the scores are shrunk while the covariate
// adjustment stays free.
// ---------------------------------------------------------------------
static double coxPenalty(int p, const double * beta, double alpha, double lambda,
   const int * pen)
   {
   double total = 0.0;
   for (int j = 0; j < p; j++)
      if (pen[j])
         total += lambda * alpha * fabs(beta[j])
            + 0.5 * lambda * (1.0 - alpha) * beta[j] * beta[j];
   return total;
   }

// minimise 1/2 b'Hb - c'b + lambda[alpha||b||1 + (1-alpha)/2||b||2^2]
static void coxSolveQuadratic(int p, double ** H, double * c, double alpha,
   double lambda, const int * pen, double * b)
   {
   if (alpha == 0.0)
      {
      double ** A = new double * [p];
      for (int i = 0; i < p; i++)
         {
         A[i] = new double [p];
         for (int j = 0; j < p; j++)
            A[i][j] = H[i][j];
         A[i][i] += lambda * pen[i];
         b[i] = c[i];
         }
      if (!coxCholesky(A, b, p))
         {
         for (int i = 0; i < p; i++)
            {
            for (int j = 0; j < p; j++)
               A[i][j] = H[i][j];
            A[i][i] += lambda * pen[i] + 1E-8 * (1.0 + fabs(H[i][i]));
            b[i] = c[i];
            }
         if (!coxCholesky(A, b, p))
            for (int j = 0; j < p; j++)
               b[j] = 0.0;
         }
      for (int i = 0; i < p; i++)
         delete []A[i];
      delete []A;
      }
   else
      {
      double l1 = lambda * alpha;
      double qp = lambda * (1.0 - alpha);
      for (int sweep = 0; sweep < 200; sweep++)
         {
         double maxChange = 0.0;
         for (int j = 0; j < p; j++)
            {
            double rho = c[j];
            for (int k = 0; k < p; k++)
               if (k != j)
                  rho -= H[j][k] * b[k];
            double nb = 0.0;
            double hjj = H[j][j];
            if (hjj < 1E-12)
               nb = 0.0;
            else if (!pen[j])
               nb = rho / hjj;              // covariates: unpenalized
            else if (rho > l1)
               nb = (rho - l1) / (hjj + qp);
            else if (rho < -l1)
               nb = (rho + l1) / (hjj + qp);
            if (fabs(nb - b[j]) > maxChange)
               maxChange = fabs(nb - b[j]);
            b[j] = nb;
            }
         double scale = 1.0;
         for (int j = 0; j < p; j++)
            if (fabs(b[j]) > scale)
               scale = fabs(b[j]);
         if (maxChange < 1E-9 * scale)
            break;
         }
      }
   }

// Newton fit of the penalized Cox model; beta is warm-started.
static int penalizedCoxFit(int N, int p, double ** X, double * ts, int * ev, const int * sid,
   double alpha, double lambda, const int * pen, double * beta, double tol,
   int maxIter)
   {
   double ** H = new double * [p];
   for (int a = 0; a < p; a++)
      H[a] = new double [p];
   double * U = new double [p];
   double * c = new double [p];
   double * cand = new double [p];

   double ll = coxEfron(N, p, X, ts, ev, sid, beta, U, H);
   double objOld = -ll + coxPenalty(p, beta, alpha, lambda, pen);

   int iter = 0;
   for (iter = 1; iter <= maxIter; iter++)
      {
      ll = coxEfron(N, p, X, ts, ev, sid, beta, U, H);
      for (int j = 0; j < p; j++)
         {
         double s = U[j];
         for (int k = 0; k < p; k++)
            s += H[j][k] * beta[k];
         c[j] = s;
         cand[j] = beta[j];
         }
      coxSolveQuadratic(p, H, c, alpha, lambda, pen, cand);

      // backtrack on the penalized objective; a tight comparison would
      // reject valid Newton steps, so allow a small relative increase
      double llNew = coxEfron(N, p, X, ts, ev, sid, cand, NULL, NULL);
      double objNew = -llNew + coxPenalty(p, cand, alpha, lambda, pen);
      double slack = 1E-8 * (1.0 + fabs(objOld));
      int tries = 0;
      while (objNew > objOld + slack && tries < 30)
         {
         for (int j = 0; j < p; j++)
            cand[j] = 0.5 * (cand[j] + beta[j]);
         llNew = coxEfron(N, p, X, ts, ev, sid, cand, NULL, NULL);
         objNew = -llNew + coxPenalty(p, cand, alpha, lambda, pen);
         tries++;
         }

      double maxChange = 0.0;
      for (int j = 0; j < p; j++)
         {
         double d = cand[j] - beta[j];
         if (fabs(d) > maxChange)
            maxChange = fabs(d);
         beta[j] = cand[j];
         }
      objOld = objNew;
      if (!(objNew == objNew))
         break;
      if (maxChange < tol)
         break;
      }

   for (int a = 0; a < p; a++)
      delete []H[a];
   delete []H;
   delete []U;
   delete []c;
   delete []cand;
   return iter;
   }

// Harrell's concordance index for a risk score
static double coxCIndex(int N, double * risk, double * t, int * ev, const int * sid)
   {
   double concordant = 0.0;
   double comparable = 0.0;
   for (int i = 0; i < N; i++)
      {
      if (!ev[i])
         continue;
      for (int j = 0; j < N; j++)
         {
         if (i == j || t[j] <= t[i])
            continue;
         if (sid && sid[i] != sid[j])      // stratified: compare within a stratum only
            continue;
         comparable += 1.0;
         if (risk[i] > risk[j])
            concordant += 1.0;
         else if (risk[i] == risk[j])
            concordant += 0.5;
         }
      }
   return comparable > 0.0 ? concordant / comparable : 0.5;
   }

// lambda path from the null-model gradient
static void coxLambdaPath(int N, int p, double ** X, double * ts, int * ev, const int * sid,
   double alpha, const int * pen, Vector & path, int nlambda, double minRatio)
   {
   double * beta = new double [p];
   double ** H = new double * [p];
   for (int a = 0; a < p; a++)
      {
      beta[a] = 0.0;
      H[a] = new double [p];
      }
   double * U = new double [p];
   coxEfron(N, p, X, ts, ev, sid, beta, U, H);
   double cmax = 0.0;
   for (int j = 0; j < p; j++)
      if (pen[j] && fabs(U[j]) > cmax)
         cmax = fabs(U[j]);
   if (cmax < 1E-12)
      cmax = 1E-12;
   double lmax = alpha > 0.0 ? cmax / alpha : cmax;
   double lmin = lmax * minRatio;
   if (lmin <= 0.0)
      lmin = lmax * 1E-4;
   path.Dimension(nlambda);
   if (nlambda == 1)
      path[0] = lmax;
   else
      {
      double lr = log(lmin / lmax);
      for (int i = 0; i < nlambda; i++)
         path[i] = exp(log(lmax) + lr * i / (nlambda - 1));
      }
   for (int a = 0; a < p; a++)
      delete []H[a];
   delete []H;
   delete []U;
   delete []beta;
   }

// The order the Cox likelihood needs: by stratum (ascending), and within a
// stratum by DECREASING time.  order[k] is the index, into t (and sid), of the
// k-th sample.  sid may be NULL (one stratum), which gives plain decreasing time.
static void coxOrder(int n, const double * t, const int * sid, IntArray & order)
   {
   order.Dimension(n);
   if (n == 0)
      return;
   Vector tv(n);
   for (int i = 0; i < n; i++)
      tv[i] = t[i];
   QuickIndex qi;
   qi.Index(tv);                       // ascending time
   if (!sid)
      {
      for (int k = 0; k < n; k++)
         order[k] = qi[n - 1 - k];
      return;
      }
   // dense time rank (tied times share a rank, so they stay together), then sort on
   // stratum * (n + 1) + (n - rank): stratum ascending, time descending.  All the
   // keys are integers well inside the exact range of a double.
   Vector key(n);
   int rank = 0;
   for (int k = 0; k < n; k++)
      {
      if (k > 0 && tv[qi[k]] != tv[qi[k - 1]])
         rank++;
      key[qi[k]] = (double) sid[qi[k]] * (double) (n + 1) + (double) (n - rank);
      }
   QuickIndex qk;
   qk.Index(key);
   for (int k = 0; k < n; k++)
      order[k] = qk[k];
   }

// Build the sorted arrays (by stratum, then decreasing time) for a subset of rows.
// `subset` holds row indices and is not modified.  sid (NULL = no strata) is the
// stratum of every row; the strata of the sorted rows are returned in sids.
static void coxSortSubset(int n, int p, double ** X, double * t, int * ev, const int * sid,
   const int * subset, double ** Xs, double * ts, int * evs, int * sids)
   {
   Vector tv(n);
   IntArray sv(n);
   for (int i = 0; i < n; i++)
      {
      tv[i] = t[subset[i]];
      sv[i] = sid ? sid[subset[i]] : 0;
      }
   IntArray order;
   coxOrder(n, &tv[0], sid ? &sv[0] : NULL, order);
   for (int i = 0; i < n; i++)
      {
      int k = subset[order[i]];
      ts[i] = t[k];
      evs[i] = ev[k];
      sids[i] = sid ? sid[k] : 0;
      for (int c = 0; c < p; c++)
         Xs[i][c] = X[k][c];
      }
   }

// Adjusted discrimination of a score: fit  Surv ~ score + covariates  and
// return the C-index of the fitted linear predictor.  `scoreSorted` is in the
// same (decreasing time) order as ts/ev, and Xs holds the covariates of those
// samples, row-major N x (covCount+1) (the last column is not used here).
// Returns _NAN_ when there are no covariates or the fit fails, so the caller
// falls back to the unadjusted C.
static double coxAdjustedC(int N, const double * scoreSorted, int covCount,
   Vector & Xs, double * ts, int * ev, const int * sid)
   {
   if (covCount <= 0)
      return _NAN_;
   int p = covCount + 1;
   double ** X = new double * [N];
   for (int i = 0; i < N; i++)
      X[i] = new double [p];
   for (int k = 0; k < N; k++)
      {
      for (int c = 0; c < covCount; c++)
         X[k][c] = Xs[k * p + c];
      X[k][covCount] = scoreSorted[k];
      }
   double * beta = new double [p];
   double * se = new double [p];
   int it = CoxFit(N, p, X, ts, ev, sid, beta, se, 1E-8, 25);
   double c = _NAN_;
   if (it >= 0 && se[covCount] > 0.0 && beta[covCount] == beta[covCount])
      {
      double * risk = new double [N];
      for (int k = 0; k < N; k++)
         {
         double s = 0.0;
         for (int j = 0; j < p; j++)
            s += X[k][j] * beta[j];
         risk[k] = s;
         }
      c = coxCIndex(N, risk, ts, ev, sid);
      delete []risk;
      }
   for (int i = 0; i < N; i++)
      delete []X[i];
   delete []X;
   delete []beta;
   delete []se;
   return c;
   }

// ---------------------------------------------------------------------
// MakeSurvVariable
// ---------------------------------------------------------------------
// The survival counterpart of MakeVariable in Logistic.cpp: prepare the samples,
// time, event and covariate matrix that Cox_BuildModel / Cox_ReadModel work on.
// Sample i is used when use[i] is non-zero (the flags of the complete cases,
// with the samples of a train or test set switched off in the other set).
//
//   allX     idCount x covariateCount covariates
//   allT     survival time of each sample
//   allEv    event code of each sample (0 = censored, otherwise an event)
//   allStrat stratum of each sample (0, 1, ...), or an empty array for no strata
//   use      idCount flags
//   pheno    output, the used samples (columns of risks) by stratum and, within a
//            stratum, by DECREASING time: the order the Cox likelihood needs
//   ts, ev   output, time and 0/1 event in that order
//   sid      output, the stratum in that order (empty when there are no strata)
//   Xs       output, row-major N x (covariateCount+1): covariates of the k-th
//            sample are Xs[k*(covariateCount+1) + c]; the last column is
//            scratch space for the PGS and is set to 0 here
static void MakeSurvVariable(Matrix & allX, Vector & allT, IntArray & allEv, IntArray & allStrat,
   IntArray & use, IntArray & pheno, Vector & ts, IntArray & ev, IntArray & sid, Vector & Xs)
   {
   const int idCount = allT.Length();
   const int covariateCount = allX.cols;
   const int d = covariateCount + 1;
   IntArray keep(0);
   for (int i = 0; i < idCount; i++)
      if (use[i])
         keep.Push(i);
   const int N = keep.Length();
   const bool stratified = (allStrat.Length() == idCount);
   pheno.Dimension(N);
   ts.Dimension(N);
   ev.Dimension(N);
   sid.Dimension(stratified ? N : 0);
   Xs.Dimension(N * d);
   if (N == 0)
      return;
   Vector timeVec(N);
   IntArray strataVec(N);
   for (int k = 0; k < N; k++)
      {
      timeVec[k] = allT[keep[k]];
      strataVec[k] = stratified ? allStrat[keep[k]] : 0;
      }
   IntArray timeOrder;
   coxOrder(N, &timeVec[0], stratified ? &strataVec[0] : NULL, timeOrder);
   for (int k = 0; k < N; k++)
      {
      int ii = keep[timeOrder[k]];               // by stratum, then descending time
      pheno[k] = ii;
      ts[k] = allT[ii];
      ev[k] = (allEv[ii] != 0) ? 1 : 0;
      if (stratified)
         sid[k] = allStrat[ii];
      for (int c = 0; c < covariateCount; c++)
         Xs[k * d + c] = allX[ii][c];
      Xs[k * d + covariateCount] = 0.0;
      }
   }

// ---------------------------------------------------------------------
// Cox_BuildModel
// ---------------------------------------------------------------------
// Status codes returned by Cox_BuildModel.
enum {
   COX_OK = 0,
   COX_TOO_FEW_SAMPLES = 1,         // N < 10
   COX_TOO_FEW_EVENTS = 2,          // fewer than 5 events
   COX_NO_VALID_PGS = 3             // the Cox fit failed for every PGS
};

// By-products of Cox_BuildModel: everything the scan reports in-sample.
struct CoxModelStats {
   int modelCount;                  // number of top PGSs in the combined models
   // one-PGS-at-a-time Cox scan, length rsCount; a failed fit has coef == _NAN_,
   // se == 0, zval == _NAN_.  The hazard ratio and the p-value are not kept
   // (they are not needed to build or cross-validate the models): exp(coef) and
   // 2 * ndist(|zval|) give them when a table is printed.
   Vector coef, se, zval;
   IntArray sortedIndex;            // PGSs with a fit, by decreasing Z^2 (see GenomicControlNumber)
   double cBest, cBestAdj;          // best single PGS: C / covariate-adjusted C
   int wsumStatus;                  // 0 no usable coefficients, 1 fitted, 2 did not converge
   int wsumTerms;                   // number of the top PGSs WTSUM uses, after pruning for correlation (--pruner2)
   double hrWsum, zWsum, pWsum, cWsum, cWsumAdj;
   double lambda[3], cIn[3], cInAdj[3];   // ridge, elastic net, lasso
   int nonzero[3];
};

// Build the prediction models for ONE survival outcome; nothing is printed.
//
//   risks       rsCount x (all samples) PGS matrix
//   pheno       the samples used by this model (columns of risks), by decreasing
//               survival time; tsV and evA are their time and 0/1 event, and Xs
//               their covariates, as made by MakeSurvVariable
//   nTop        number of top PGSs to combine (<= 0: 30), taken from the PGSs
//               ranked by the square of the Wald Z (not by p-value: it
//               underflows to 0 for |Z| above ~20 and the ranking would tie)
//   LD          rsCount x rsCount bit matrix of the PGS pairs with R2 > --pruner2
//               (MakeBinaryLD).  The top nTop PGSs are the terms of every model; only
//               WTSUM uses the subset left after pruning for correlation, in the order of
//               the ranking (MakePrunedOrder), as --linear and --logistic do.
//
//   weights     output, modelCount x 5: columns BEST, WTSUM, RIDGE, NET, LASSO.
//               Row j belongs to PGS modelIndex[j] (the j-th best by Wald p-value),
//               and the coefficients are expressed against the risks matrix that
//               was passed in.
//   modelIndex  output, the PGS (row of risks) of each weights row
//
// Returns COX_OK, or COX_TOO_FEW_* / COX_NO_VALID_PGS, in which case weights
// has no rows.
// The ridge / elastic net / lasso penalties are chosen by an inner 10-fold
// cross-validation of the C-index; that C is not reported (it is optimistic,
// because the same folds pick lambda), use --cv for an honest estimate.
static int Cox_BuildModel(Matrix & risks, IntArray & pheno, Vector & tsV,
   IntArray & evA, IntArray & sidA, Vector & Xs, int covariateCount, int nTop, int coreCount,
   unsigned long long int ** LD, Matrix & weights, IntArray & modelIndex, CoxModelStats & st)
   {
   const int rsCount = risks.rows;
   const int N = pheno.Length();
   const int p = covariateCount + 1;     // covariates + the PGS being tested
   weights.Dimension(0, 5);
   modelIndex.Dimension(0);
   st.modelCount = 0;
   st.cBest = 0.5;
   st.cBestAdj = _NAN_;
   st.wsumStatus = 0;
   st.wsumTerms = 0;
   st.hrWsum = st.zWsum = st.pWsum = st.cWsum = st.cWsumAdj = _NAN_;
   for (int a2 = 0; a2 < 3; a2++)
      {
      st.lambda[a2] = st.cIn[a2] = st.cInAdj[a2] = _NAN_;
      st.nonzero[a2] = 0;
      }
   if (N < 10)
      return COX_TOO_FEW_SAMPLES;
   int nEvents = 0;
   for (int k = 0; k < N; k++)
      if (evA[k])
         nEvents++;
   if (nEvents < 5)
      return COX_TOO_FEW_EVENTS;
   double * ts = &tsV[0];
   int * ev = &evA[0];
   const int * sid = (sidA.Length() == N) ? &sidA[0] : NULL;   // NULL: no strata

   // ---- one-PGS-at-a-time scan, adjusted for the covariates ----
   Vector & coef = st.coef, & se = st.se, & zval = st.zval;
   coef.Dimension(rsCount);
   se.Dimension(rsCount);
   zval.Dimension(rsCount);
   coef.Set(_NAN_);
   se.Zero();
   zval.Set(_NAN_);

#ifdef _OPENMP
#pragma omp parallel num_threads(coreCount)
   {
#endif
   double ** X = new double * [N];
   for (int i = 0; i < N; i++)
      X[i] = new double [p];
   for (int k = 0; k < N; k++)
      for (int c = 0; c < covariateCount; c++)
         X[k][c] = Xs[k * p + c];
   double * beta = new double [p];
   double * betaSe = new double [p];

#ifdef _OPENMP
#pragma omp for
#endif
   for (int r = 0; r < rsCount; r++)
      {
      for (int k = 0; k < N; k++)
         X[k][covariateCount] = risks[r][pheno[k]];
      int iter = CoxFit(N, p, X, ts, ev, sid, beta, betaSe, 1E-8, 25);
      if (iter >= 0 && betaSe[covariateCount] > 0.0
         && beta[covariateCount] == beta[covariateCount])
         {
         coef[r] = beta[covariateCount];
         se[r] = betaSe[covariateCount];
         zval[r] = beta[covariateCount] / betaSe[covariateCount];
         }
      }

   delete []beta;
   delete []betaSe;
   for (int i = 0; i < N; i++)
      delete []X[i];
   delete []X;
#ifdef _OPENMP
   } // end omp parallel
#endif

   // rank by the square of Z, as --logistic does; the PGSs whose fit failed are left out
   IntArray & sortedIndex = st.sortedIndex;
   sortedIndex.Dimension(0);
   GenomicControlNumber(zval, sortedIndex);
   const int validCount = sortedIndex.Length();
   if (validCount == 0)
      return COX_NO_VALID_PGS;

   // ---- penalized Cox on the top PGSs (always fitted) ------------------
   int m = nTop > 0 ? nTop : 30;
   if (m > validCount)
      m = validCount;
   if (m < 1)
      error("The penalized Cox model needs at least one PGS.");
   st.modelCount = m;
   modelIndex.Dimension(m);
   for (int j = 0; j < m; j++)
      modelIndex[j] = sortedIndex[j];
   // the PGSs WTSUM uses: those left after pruning for correlation (--pruner2)
   IntArray wsumFlag(m);                  // 1 = the PGS is used by WTSUM
   wsumFlag.Zero();
   {
   IntArray kept = modelIndex;
   MakePrunedOrder(LD, kept);             // positions in modelIndex of the PGSs that survive
   for (int j = 0; j < kept.Length(); j++)
      wsumFlag[kept[j]] = 1;
   st.wsumTerms = kept.Length();
   }

   // design: the selected PGS (penalized) followed by the covariates
   // (unpenalized), so the joint model is adjusted like the per-PGS scan
   int total = m + covariateCount;
   Matrix Xnet(N, total);
   for (int i = 0; i < N; i++)
      {
      for (int j = 0; j < m; j++)
         Xnet[i][j] = risks[sortedIndex[j]][pheno[i]];
      for (int c = 0; c < covariateCount; c++)
         Xnet[i][m + c] = Xs[i * p + c];
      }
   int * pen = new int [total];
   for (int j = 0; j < total; j++)
      pen[j] = (j < m) ? 1 : 0;
   double ** Xr = new double * [N];
   for (int i = 0; i < N; i++)
      Xr[i] = &Xnet[i][0];

   double alphas[3] = {0.0, 0.5, 1.0};
   Matrix netCoefAll(m, 3);     // ridge / net / lasso coefficients
   netCoefAll.Set(_NAN_);

   int nlambda = 20;
   int K = 10;
   if (K > N)
      K = N;

   for (int a2 = 0; a2 < 3; a2++)
      {
      double alpha = alphas[a2];
      Vector path;
      coxLambdaPath(N, total, Xr, ts, ev, sid, alpha, pen, path, nlambda, 1E-4);
      int nl = path.Length();
      Vector cvC(nl);
      cvC.Zero();

      // deterministic fold assignment
      IntArray idx(N);
      for (int i = 0; i < N; i++)
         idx[i] = i;
      unsigned long long st64 = 20261020ULL + 1000ULL * (unsigned long long) a2;
      for (int i = N - 1; i > 0; i--)
         {
         st64 ^= st64 >> 12;
         st64 ^= st64 << 25;
         st64 ^= st64 >> 27;
         int k2 = (int) ((st64 * 2685821657736338717ULL)
            % (unsigned long long) (i + 1));
         int tmp = idx[i];
         idx[i] = idx[k2];
         idx[k2] = tmp;
         }

      for (int k = 0; k < K; k++)
         {
         int lo = (int) ((double) k * N / K);
         int hi = (int) ((double) (k + 1) * N / K);
         int ntest = hi - lo;
         int ntrain = N - ntest;
         if (ntest < 5 || ntrain < 20)
            continue;

         int * tr = new int [ntrain];
         int * te = new int [ntest];
         int nt = 0, nv = 0;
         for (int i = 0; i < N; i++)
            {
            if (i >= lo && i < hi)
               te[nv++] = idx[i];
            else
               tr[nt++] = idx[i];
            }

         double ** Xt = new double * [ntrain];
         for (int i = 0; i < ntrain; i++)
            Xt[i] = new double [total];
         double * tt = new double [ntrain];
         int * et = new int [ntrain];
         int * sidT = new int [ntrain];
         coxSortSubset(ntrain, total, Xr, ts, ev, sid, tr, Xt, tt, et, sidT);

         double ** Xv = new double * [ntest];
         for (int i = 0; i < ntest; i++)
            Xv[i] = Xr[te[i]];
         double * tv = new double [ntest];
         int * evv = new int [ntest];
         int * sidV = new int [ntest];
         for (int i = 0; i < ntest; i++)
            {
            tv[i] = ts[te[i]];
            evv[i] = ev[te[i]];
            sidV[i] = sid ? sid[te[i]] : 0;
            }

         double * beta = new double [total];
         for (int j = 0; j < total; j++)
            beta[j] = 0.0;
         double * risk = new double [ntest];
         for (int li = 0; li < nl; li++)
            {
            penalizedCoxFit(ntrain, total, Xt, tt, et, sid ? sidT : NULL, alpha, path[li],
               pen, beta, 1E-7, 50);
            for (int i = 0; i < ntest; i++)
               {
               double sum = 0.0;
               for (int j = 0; j < total; j++)
                  sum += Xv[i][j] * beta[j];
               risk[i] = sum;
               }
            cvC[li] += coxCIndex(ntest, risk, tv, evv, sid ? sidV : NULL) / (double) K;
            }

         delete []risk;
         delete []beta;
         for (int i = 0; i < ntrain; i++)
            delete []Xt[i];
         delete []Xt;
         delete []tt;
         delete []et;
         delete []sidT;
         delete []Xv;
         delete []tv;
         delete []evv;
         delete []sidV;
         delete []tr;
         delete []te;
         } // End of k loop

      int best = 0;
      for (int li = 1; li < nl; li++)
         if (cvC[li] > cvC[best])
            best = li;

      double * beta = new double [total];
      for (int j = 0; j < total; j++)
         beta[j] = 0.0;
      penalizedCoxFit(N, total, Xr, ts, ev, sid, alpha, path[best], pen, beta,
         1E-8, 100);
      for (int j = 0; j < m; j++)
         netCoefAll[j][a2] = beta[j];
      // unadjusted C ranks on the PGS part only; adjusted C uses the whole
      // linear predictor, covariates included
      double * riskP = new double [N];
      double * riskF = new double [N];
      for (int i = 0; i < N; i++)
         {
         double sp = 0.0, sf = 0.0;
         for (int j = 0; j < m; j++)
            sp += Xr[i][j] * beta[j];
         for (int j = m; j < total; j++)
            sf += Xr[i][j] * beta[j];
         riskP[i] = sp;
         riskF[i] = sp + sf;
         }
      st.cIn[a2] = coxCIndex(N, riskP, ts, ev, sid);
      st.cInAdj[a2] = coxCIndex(N, riskF, ts, ev, sid);
      st.lambda[a2] = path[best];
      int nnz = 0;
      for (int j = 0; j < m; j++)
         if (fabs(beta[j]) > 1E-10)
            nnz++;
      st.nonzero[a2] = nnz;
      delete []beta;
      delete []riskP;
      delete []riskF;
      } // End of a2 loop
   delete []Xr;
   delete []pen;

   // ---- WTSUM: weighted sum of the top PGSs, refitted ----
   // Weights are the univariate log hazard ratios from the scan, scaled to
   // unit norm, so the combined score is on a comparable scale to a single
   // PGS.  Unlike the penalized models above, the refit is on the score
   // alone; cWsumAdj adds the covariates back for the adjusted C.
   Vector wsum(N);
   Vector wsumW(m);
   wsum.Zero();
   wsumW.Zero();
   double wnorm = 0.0;
   for (int j = 0; j < m; j++)
      {
      int rs = sortedIndex[j];
      if (coef[rs] == _NAN_ || !wsumFlag[j])   // a PGS pruned out of WTSUM is skipped
         continue;
      for (int i = 0; i < N; i++)
         wsum[i] += coef[rs] * risks[rs][pheno[i]];
      wsumW[j] = coef[rs];
      wnorm += coef[rs] * coef[rs];
      }
   if (wnorm > 0.0)
      {
      double wsc = sqrt(wnorm);
      for (int i = 0; i < N; i++)
         wsum[i] /= wsc;
      for (int j = 0; j < m; j++)
         wsumW[j] /= wsc;
      double ** Xw = new double * [N];
      for (int i = 0; i < N; i++)
         {
         Xw[i] = new double [1];
         Xw[i][0] = wsum[i];
         }
      double betaW = 0.0, seW = 0.0;
      int itW = CoxFit(N, 1, Xw, ts, ev, sid, &betaW, &seW, 1E-8, 25);
      if (itW >= 0 && seW > 0.0 && betaW == betaW)
         {
         st.wsumStatus = 1;
         st.hrWsum = exp(betaW);
         st.zWsum = betaW / seW;
         st.pWsum = 2.0 * ndist(fabs(st.zWsum));
         st.cWsum = coxCIndex(N, &wsum[0], ts, ev, sid);
         st.cWsumAdj = coxAdjustedC(N, &wsum[0], covariateCount, Xs, ts, ev, sid);
         if (st.cWsumAdj == _NAN_)
            st.cWsumAdj = st.cWsum;
         }
      else
         st.wsumStatus = 2;
      for (int i = 0; i < N; i++)
         delete []Xw[i];
      delete []Xw;
      }

   // concordance of the best single PGS, for reference
   if (coef[sortedIndex[0]] != _NAN_)
      {
      double * risk = new double [N];
      for (int i = 0; i < N; i++)
         risk[i] = risks[sortedIndex[0]][pheno[i]];
      st.cBest = coxCIndex(N, risk, ts, ev, sid);
      st.cBestAdj = coxAdjustedC(N, risk, covariateCount, Xs, ts, ev, sid);
      delete []risk;
      }
   if (st.cBestAdj == _NAN_)
      st.cBestAdj = st.cBest;

   // ---- the model weights ----
   weights.Dimension(m, 5);
   weights.Zero();
   for (int j = 0; j < m; j++)
      {
      if (j == 0)
         weights[j][0] = 1.0;          // best single PGS
      // The weights are against the risks matrix as the scan sees it, standardised in
      // place by MakeBinaryLD; survival_Scan divides them by rsSD for the model file.
      weights[j][1] = wsumW[j];
      for (int a2 = 0; a2 < 3; a2++)
         weights[j][2 + a2] = netCoefAll[j][a2];
      }
   return COX_OK;
   }

// ---------------------------------------------------------------------
// Cox_ReadModel
// ---------------------------------------------------------------------
// Evaluate a set of PGS models (the columns of coef, as written by --coxph or
// returned by Cox_BuildModel) on ONE survival outcome: the C-index of every
// combined score (cRaw), and the C-index of the linear predictor of
//     Surv ~ combined score + covariates      (cAdj; = cRaw without covariates).
//
//   pgsName, coef   one name per row of coef (modelCount x methodCount); a _NAN_
//                   coefficient is skipped
//   rsNames, risks  the PGSs available; model terms not found in rsNames are skipped
//   pheno, tsV, evA, Xs   the samples, as made by MakeSurvVariable
//   cRaw, cAdj      output, length methodCount, _NAN_ when the model cannot be
//                   evaluated (fewer than 10 samples or fewer than 2 events)
//
// Returns the number of model terms that were not found in rsNames.
static int Cox_ReadModel(StringArray & pgsName, Matrix & coef, StringArray & rsNames,
   Matrix & risks, IntArray & pheno, Vector & tsV, IntArray & evA, IntArray & sidA, Vector & Xs,
   int covariateCount, Vector & cRaw, Vector & cAdj)
   {
   const int N = pheno.Length();
   const int modelrsCount = coef.rows;
   const int methodCount = coef.cols;
   cRaw.Dimension(methodCount);
   cAdj.Dimension(methodCount);
   cRaw.Set(_NAN_);
   cAdj.Set(_NAN_);
   IntArray rsIndex(modelrsCount);
   int missingTerms = 0;
   for (int r = 0; r < modelrsCount; r++)
      {
      rsIndex[r] = rsNames.Find(pgsName[r]);
      if (rsIndex[r] == -1)
         missingTerms++;
      }
   int nEvents = 0;
   for (int k = 0; k < N; k++)
      if (evA[k])
         nEvents++;
   if (N < 10 || nEvents < 2)
      return missingTerms;
   double * ts = &tsV[0];
   int * ev = &evA[0];
   const int * sid = (sidA.Length() == N) ? &sidA[0] : NULL;   // NULL: no strata
   Vector score(N);
   for (int mm = 0; mm < methodCount; mm++)
      {
      score.Zero();
      for (int r = 0; r < modelrsCount; r++)
         {
         int rr = rsIndex[r];
         double bj = coef[r][mm];
         if (rr == -1 || bj == _NAN_)
            continue;
         for (int k = 0; k < N; k++)
            score[k] += risks[rr][pheno[k]] * bj;
         }
      cRaw[mm] = coxCIndex(N, &score[0], ts, ev, sid);
      cAdj[mm] = coxAdjustedC(N, &score[0], covariateCount, Xs, ts, ev, sid);
      if (cAdj[mm] == _NAN_)
         cAdj[mm] = cRaw[mm];
      }
   return missingTerms;
   }

// ---------------------------------------------------------------------
// Cross-validation
// ---------------------------------------------------------------------
// One round: the samples flagged in testFlag form the test set, all the others
// the train set.  The models (BEST, WTSUM, RIDGE, NET, LASSO) are built in the
// train set by Cox_BuildModel and evaluated in the test set by Cox_ReadModel.
//
//   valid            flags of the complete cases (see survival_Scan)
//   cRaw, cAdj       output, length 5 (the columns of the weights matrix), _NAN_
//                    when a model could not be built or evaluated in this round
static void ScanCV1Time(Matrix & risks, StringArray & rsNames, Matrix & allX,
   Vector & allT, IntArray & allEv, IntArray & allStrat, IntArray & valid, int nTop,
   unsigned long long int ** LD, unsigned long long int * testFlag, Vector & cRaw, Vector & cAdj, int coreCount)
   {
   const int idCount = allT.Length();
   const int covariateCount = allX.cols;
   cRaw.Dimension(5);
   cAdj.Dimension(5);
   cRaw.Set(_NAN_);
   cAdj.Set(_NAN_);
   IntArray trainUse(idCount), testUse(idCount);
   for (int i = 0; i < idCount; i++)
      {
      bool inTest = (testFlag[i >> 6] >> (i & 0x3F)) & 0x1;
      trainUse[i] = (valid[i] && !inTest) ? 1 : 0;
      testUse[i] = (valid[i] && inTest) ? 1 : 0;
      }
   IntArray trainPheno, trainEv, trainSid, testPheno, testEv, testSid;
   Vector trainTs, trainXs, testTs, testXs;
   MakeSurvVariable(allX, allT, allEv, allStrat, trainUse, trainPheno, trainTs, trainEv, trainSid, trainXs);
   MakeSurvVariable(allX, allT, allEv, allStrat, testUse, testPheno, testTs, testEv, testSid, testXs);
   Matrix weights;
   IntArray modelIndex;
   CoxModelStats st;
   int status = Cox_BuildModel(risks, trainPheno, trainTs, trainEv, trainSid, trainXs,
      covariateCount, nTop, coreCount, LD, weights, modelIndex, st);
   int modelCount = modelIndex.Length();
   if (status != COX_OK || modelCount == 0 || weights.rows != modelCount)
      return;
   StringArray pgsName(modelCount);
   for (int j = 0; j < modelCount; j++)
      pgsName[j] = rsNames[modelIndex[j]];
   Vector raw, adj;
   Cox_ReadModel(pgsName, weights, rsNames, risks, testPheno, testTs, testEv,
      testSid, testXs, covariateCount, raw, adj);
   for (int mm = 0; mm < 5; mm++)
      {
      cRaw[mm] = raw[mm];
      cAdj[mm] = adj[mm];
      }
   }

// Repeated random-split cross-validation: every run draws idCount/fold samples
// as the test set (testFlags[run], one bit per sample, allocated by the caller
// with ceil(idCount/64) words) and calls ScanCV1Time.  cvRaw and cvAdj are
// repeat x 5.
static void ScanCV(Matrix & risks, StringArray & rsNames, Matrix & allX,
   Vector & allT, IntArray & allEv, IntArray & allStrat, IntArray & valid, int fold, int repeat,
   int nTop, unsigned long long int ** LD, Matrix & cvRaw, Matrix & cvAdj,
   unsigned long long int ** testFlags, int coreCount)
   {
   const int idCount = allT.Length();
   Random ran;
   IntArray testFlag(idCount);
   int imemsize = ((idCount - 1) / 64 + 1) * sizeof(unsigned long long int);
   int testsetCount = idCount * 1.0 / fold;
   for (int run = 0; run < repeat; run++)
      {
      ran.Choose(testFlag, idCount, testsetCount);
      memset(&testFlags[run][0], 0, imemsize);
      for (int i = 0; i < idCount; i += 64)
         {
         int iMax = i < idCount - 64 ? i + 64 : idCount;
         int w = (i >> 6);
         unsigned long long int word0 = 0x1;
         for (int ii = i; ii < iMax; ii++, word0 <<= 1)
            if (testFlag[ii]) testFlags[run][w] |= word0;
         } // End of i loop
      } // End of run loop
   cvRaw.Dimension(repeat, 5);
   cvAdj.Dimension(repeat, 5);
   for (int run = 0; run < repeat; run++)
      {
      printf(".");
      fflush(stdout);
      Vector raw, adj;
      ScanCV1Time(risks, rsNames, allX, allT, allEv, allStrat, valid, nTop, LD, testFlags[run],
         raw, adj, coreCount);
      for (int mm = 0; mm < 5; mm++)
         {
         cvRaw[run][mm] = raw[mm];
         cvAdj[run][mm] = adj[mm];
         }
      } // End of run loop
   }

// ---------------------------------------------------------------------
// plot-only R script, driven by the table computed above
// ---------------------------------------------------------------------
static int plotCoxHazard(const char * statfile, const char * pdffile,
   const char * title, const char * pgsfile, const char * corrfile, const char * rpath)
   {
   String scriptfile(pdffile);
   scriptfile.Add(".R");
   FILE * fp = fopen(scriptfile, "wt");
   if (fp == NULL)
      error("Cannot open %s to write.", (const char*)scriptfile);
   fprintf(fp, "## %s for SECRET Cox hazard-ratio plot, by Wei-Min Chen\n",
      (const char*)scriptfile);
   fprintf(fp, "rm(list=ls(all=TRUE))\n");
   WritePGSLabelR(fp);
   fprintf(fp, "table <- read.table(\"%s\", header = TRUE, sep = \"\\t\", quote = \"\\\"\", comment.char = \"\")\n",
      statfile);
   fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n",
      pdffile);
   fprintf(fp, "table <- table[order(-table$z * table$z), ]   # rank the PGSs by -Z*Z, as --logistic does\n");
   fprintf(fp, "n <- min(30, nrow(table))   # top 30; fewer when fewer PGSs were fitted\n");
   fprintf(fp, "chisq <- table$z[1:n] * table$z[1:n]\n");
   fprintf(fp, "plot(chisq, 1:n, pch=19, col=c(\"blue\", \"red\")[(table$z[1:n]<0)+1],\n");
   fprintf(fp, "xlab = \"Z-Squared between PGS and %s\", ylab = \"PGS Rank (by Z-Squared)\", main = paste0(\"Top \", n, \" PGS Scores That Predict %s\"))\n", title, title);
   fprintf(fp, "mid <- (min(chisq) + max(chisq))/2\n");
   fprintf(fp, "labs <- pgslabels(table[1:n, , drop=FALSE], table)   # PGS:SNPs:Coverage:trait, as far as known\n");
   fprintf(fp, "text(chisq, 1:n, labels = labs, pos = ifelse(chisq < mid, 4, 2))\n");
   WriteTopHeatmapR(fp, corrfile, statfile, TOPHEAT_HR, title);
   fprintf(fp, "# ---- combined PGS of the models (written by --coxph --rplot): survival ~ score ----\n");
   fprintf(fp, "pgsfile <- \"%s\"\n", pgsfile);
   fprintf(fp, "if(file.exists(pgsfile)){\n");
   fprintf(fp, "  data <- read.table(pgsfile, sep=\"\\t\", header=TRUE, row.names=1)\n");
   fprintf(fp, "  tm <- data[,1]; ev <- data[,2]\n");
   fprintf(fp, "  cindex <- function(risk, tm, ev){\n");
   fprintf(fp, "    num <- 0; den <- 0\n");
   fprintf(fp, "    for(i in which(ev==1)){\n");
   fprintf(fp, "      j <- tm > tm[i]\n");
   fprintf(fp, "      if(any(j)){ den <- den + sum(j); num <- num + sum(risk[i] > risk[j]) + 0.5*sum(risk[i]==risk[j]) }\n");
   fprintf(fp, "    }\n");
   fprintf(fp, "    if(den > 0) num/den else NA\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "  km <- function(tm, ev){\n");
   fprintf(fp, "    o <- order(tm); tm <- tm[o]; ev <- ev[o]\n");
   fprintf(fp, "    ut <- unique(tm[ev==1]); s <- 1; surv <- numeric(length(ut))\n");
   fprintf(fp, "    for(i in seq_along(ut)){ n <- sum(tm >= ut[i]); d <- sum(tm==ut[i] & ev==1); s <- s*(1 - d/n); surv[i] <- s }\n");
   fprintf(fp, "    list(t=c(0, ut), s=c(1, surv))\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "  logrank <- function(tm, ev, g){\n");
   fprintf(fp, "    lev <- levels(g); k <- length(lev); O <- numeric(k); E <- numeric(k); V <- matrix(0, k, k)\n");
   fprintf(fp, "    for(t in sort(unique(tm[ev==1]))){\n");
   fprintf(fp, "      at <- tm >= t; n <- sum(at); d <- sum(tm==t & ev==1)\n");
   fprintf(fp, "      nj <- sapply(lev, function(l) sum(at & g==l)); dj <- sapply(lev, function(l) sum(tm==t & ev==1 & g==l))\n");
   fprintf(fp, "      O <- O + dj; E <- E + d*nj/n\n");
   fprintf(fp, "      if(n > 1) V <- V + d*(n-d)/(n-1) * (diag(nj)/n - outer(nj, nj)/n^2)\n");
   fprintf(fp, "    }\n");
   fprintf(fp, "    z <- (O - E)[-k]\n");
   fprintf(fp, "    pchisq(as.numeric(t(z) %%*%% solve(V[-k,-k]) %%*%% z), k-1, lower.tail=FALSE)\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "  usable <- function(score){\n");
   fprintf(fp, "    ok <- !is.na(score)\n");
   fprintf(fp, "    sum(ok) >= 20 && sum(ev[ok]==1) >= 2 && sd(score[ok]) > 0\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "  combined <- c(\"BEST\", \"WTSUM\", \"RIDGE\", \"NET\", \"LASSO\")\n");
   fprintf(fp, "  labels <- c(paste0(\"Best (\", as.character(table$PGS[1]), \")\"), \"WTSUM\", \"Ridge\", \"Elastic net\", \"Lasso\")   # the best PGS is first in the scan table\n");
   fprintf(fp, "  for(k in 1:5){\n");
   fprintf(fp, "    if(!(combined[k] %%in%% names(data))) next\n");
   fprintf(fp, "    score <- data[[combined[k]]]\n");
   fprintf(fp, "    if(!usable(score)){\n");
   fprintf(fp, "      par(mfrow=c(1,1))\n");
   fprintf(fp, "      plot.new()\n");
   fprintf(fp, "      text(0.5, 0.5, paste0(labels[k], \" PGS on %s is not available or is constant\\n(no PGS left in the model)\"))\n", title);
   fprintf(fp, "      next\n");
   fprintf(fp, "    }\n");
   fprintf(fp, "    ok <- !is.na(score); s <- score[ok]; t2 <- tm[ok]; e2 <- ev[ok]\n");
   fprintf(fp, "    brk <- unique(quantile(s, probs=c(0, 1/3, 2/3, 1)))\n");
   fprintf(fp, "    if(length(brk) < 4){\n");
   fprintf(fp, "      par(mfrow=c(1,1))\n");
   fprintf(fp, "      plot.new()\n");
   fprintf(fp, "      text(0.5, 0.5, paste0(labels[k], \" PGS on %s has too few distinct values to form tertiles\"))\n", title);
   fprintf(fp, "      next\n");
   fprintf(fp, "    }\n");
   fprintf(fp, "    g <- cut(s, breaks=brk, include.lowest=TRUE, labels=c(\"Low\", \"Middle\", \"High\"))\n");
   fprintf(fp, "    lr <- tryCatch(logrank(t2, e2, g), error=function(e) NA)\n");
   fprintf(fp, "    par(mfrow=c(1,2))\n");
   fprintf(fp, "    colr <- c(\"blue\", \"black\", \"red\")\n");
   fprintf(fp, "    plot(c(0, max(t2)), c(0, 1), type=\"n\", xlab=\"Follow-up time\", ylab=\"Survival probability\", cex.main=0.9,\n");
   fprintf(fp, "      main=paste0(labels[k], \" PGS tertiles on %s\\n(C=\", sprintf(\"%%.3f\", cindex(s, t2, e2)), \", log-rank P=\", format.pval(lr, digits=2), \")\"))\n", title);
   fprintf(fp, "    for(i in 1:3){\n");
   fprintf(fp, "      idx <- g == levels(g)[i]; r <- km(t2[idx], e2[idx])\n");
   fprintf(fp, "      lines(c(r$t, max(t2[idx])), c(r$s, tail(r$s, 1)), type=\"s\", col=colr[i], lwd=2)\n");
   fprintf(fp, "    }\n");
   fprintf(fp, "    legend(\"bottomleft\", legend=paste0(levels(g), \" (n=\", as.vector(table(g)), \", events=\", as.vector(tapply(e2, g, sum)), \")\"),\n");
   fprintf(fp, "      col=colr, lty=1, lwd=2, bty=\"n\", cex=0.8)\n");
   fprintf(fp, "    plot(s, t2, pch=ifelse(e2==1, 19, 1), col=ifelse(e2==1, \"red\", \"grey50\"), xlab=paste(labels[k], \"PGS\"),\n");
   fprintf(fp, "      ylab=\"Follow-up time\", cex.main=0.9, main=paste0(labels[k], \" PGS and follow-up time\"))\n");
   fprintf(fp, "    legend(\"topright\", legend=c(\"Event\", \"Censored\"), pch=c(19, 1), col=c(\"red\", \"grey50\"), bty=\"n\")\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "  par(mfrow=c(1,1))\n");
   fprintf(fp, "  cols <- c(\"BEST\", \"WTSUM\", \"RIDGE\", \"NET\", \"LASSO\")\n");
   fprintf(fp, "  lab <- c(\"Best PGS\", \"WTSUM\", \"Ridge\", \"Elastic net\", \"Lasso\")\n");
   fprintf(fp, "  cv <- rep(NA, 5)\n");
   fprintf(fp, "  for(k in 1:5){\n");
   fprintf(fp, "    if(!(cols[k] %%in%% names(data))) next\n");
   fprintf(fp, "    sc <- data[[cols[k]]]\n");
   fprintf(fp, "    if(!usable(sc)) next\n");
   fprintf(fp, "    ok <- !is.na(sc)\n");
   fprintf(fp, "    cv[k] <- cindex(sc[ok], tm[ok], ev[ok])\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "  if(any(!is.na(cv))){\n");
   fprintf(fp, "    dotchart(rev(cv), labels=rev(lab), pch=19, xlim=c(min(0.5, cv, na.rm=TRUE) - 0.05, max(cv, na.rm=TRUE) + 0.05),\n");
   fprintf(fp, "      xlab=\"Harrell's C-index (score only)\", main=\"Concordance of PGSs predicting %s\")\n", title);
   fprintf(fp, "    abline(v=0.5, lty=3)\n");
   fprintf(fp, "  }\n");
   fprintf(fp, "}\n");
   fprintf(fp, "dev.off()\n");
   fprintf(fp, "rm(list=ls(all=TRUE))\n");
   fclose(fp);

   String command = RCommandLine(rpath, scriptfile);   // quotes the paths, uses --rpath
   system((const char*)command);
   int error = CheckRout(scriptfile);
   if (error == 0)
      {
      printf("Hazard ratio plots saved in %s%s\n", pdffile, MakePdf(pdffile));
      }
   return error;
   }

// ---------------------------------------------------------------------
// survival_Scan
// ---------------------------------------------------------------------
void survival_Scan(const char * phefile, const char * dictfile,
   const char * rsfile, const char * prefix, const char * survTime,
   const char * survEvent, const char * strat, const char * covariateList, bool rplotFlag,
   int nTop, double ldT, const char * modelfile, bool cvFlag, int fold, int repeat,
   int coreCount)
   {
   printf("\nOptions in effect:\n");
   printf("\t--coxph\n");
   if (modelfile[0])
      printf("\t--readmodel %s\n", modelfile);
   printf("\t--surv %s,%s\n", survTime, survEvent);
   if (strat[0])
      printf("\t--strat %s\n", strat);
   printf("\t--pheno %s\n", phefile);
   if (gCovarFile[0])
      printf("\t--covarfile %s\n", gCovarFile);
   if (covariateList[0])
      printf("\t--covariate %s\n", covariateList);
   printf("\t--score %s\n", rsfile);
   if (dictfile[0])
      printf("\t--diction %s\n", dictfile);
   if (rplotFlag)
      printf("\t--rplot\n");
   if (cvFlag)
      printf("\t--cv\n");
   if (fold)
      printf("\t--fold %d\n", fold);
   if (repeat)
      printf("\t--repeat %d\n", repeat);
   if (nTop)
      printf("\t--ntop %d\n", nTop);
   if (ldT != _NAN_ && ldT < 1 && ldT >= 0)
      printf("\t--pruner2 %.3lf\n", ldT);
   else
      ldT = _NAN_;
   if (coreCount)
      printf("\t--cpus %d\n", coreCount);
   if (prefix)
      printf("\t--prefix %s\n", prefix);
   printf("\n");
   if (cvFlag)
      {
      if (fold == 0)
         fold = 10;
      if (repeat == 0)
         repeat = 10;
      if (fold < 2)
         error("--fold must be at least 2.");
      if (repeat < 1)
         error("--repeat must be at least 1.");
      }

   StringArray sampleID, alltraitNames, rsNames, annotation;
   IntArray covariates, affections, quantraits;
   Matrix alltraits, risks;
   if (gTraitList[0]) {
      printf("Note: --trait does not apply to --coxph; the outcome is set by --surv.\n");
      gTraitList = "";
   }
   gRecodeBinary01 = false;   // the event column is 0/1 by design, do not recode it to 1/2
   SECRET_Read(phefile, dictfile, rsfile, "", prefix, covariateList,
      sampleID, alltraitNames, alltraits, covariates, affections, quantraits,
      rsNames, risks, annotation, coreCount);

   int idCount = sampleID.Length();
   if (idCount == 0)
      return;
   int rsCount = rsNames.Length();

   int time = alltraitNames.Find(survTime);
   if (time == -1)
      error("Survival time %s does not exist in the phenotype data.\n", survTime);
   int event = alltraitNames.Find(survEvent);
   if (event == -1)
      error("Survival event %s does not exist in the phenotype data.\n", survEvent);

   // ---- --strat: stratified Cox model, like strata() in R's coxph ----
   // A numeric variable of the phenotype (or covariate) file, like --covariate; several
   // variables separated by commas give one stratum per combination, as strata(a, b).
   IntArray stratCols(0);
   String stratLabel;
   if (strat[0])
      {
      StringArray stratNames;
      stratNames.AddTokens(strat, ',');
      for (int s = 0; s < stratNames.Length(); s++)
         {
         int k = alltraitNames.Find(stratNames[s]);
         if (k == -1)
            error("Stratum variable %s given by --strat does not exist in the phenotype data.\n",
               (const char*)stratNames[s]);
         if (k == time || k == event)
            error("--strat cannot be the survival time or the event (%s).\n", (const char*)stratNames[s]);
         stratCols.Push(k);
         if (s)
            stratLabel.Add(",");
         stratLabel.Add(stratNames[s]);
         int cpos = covariates.Find(k);
         if (cpos > -1)
            {
            // a stratum variable is absorbed by the stratified baseline hazard
            printf("WARNING: %s is both a covariate and a stratum variable; it is used as a stratum only.\n",
               (const char*)stratNames[s]);
            covariates.Delete(cpos);
            }
         }
      }
   int covariateCount = covariates.Length();

   printf("Cox proportional hazard model:\n  Surv(%s,%s) ~ PGS", survTime, survEvent);
   for (int c = 0; c < covariateCount; c++)
      printf("+%s", (const char*)alltraitNames[covariates[c]]);
   if (stratCols.Length())
      printf("+strata(%s)", (const char*)stratLabel);
   printf("\n  ties handled by the Efron approximation\n");

   // ---- complete cases -------------------------------------------------
   IntArray valid(0);
   int skippedTime = 0, skippedCov = 0, skippedEvent = 0, skippedStrat = 0;
   for (int i = 0; i < idCount; i++)
      {
      double t = alltraits[i][time];
      if (t == _NAN_ || t <= 0.0)
         {
         skippedTime++;
         continue;
         }
      int e = (int) (alltraits[i][event] + 0.5);
      if (alltraits[i][event] == _NAN_ || (e != 0 && e != 1 && e != 2))
         {
         skippedEvent++;
         continue;
         }
      bool ok = true;
      for (int c = 0; c < covariateCount; c++)
         if (alltraits[i][covariates[c]] == _NAN_)
            {
            ok = false;
            break;
            }
      if (!ok)
         {
         skippedCov++;
         continue;
         }
      for (int s = 0; s < stratCols.Length(); s++)
         if (alltraits[i][stratCols[s]] == _NAN_)
            {
            ok = false;
            break;
            }
      if (!ok)
         {
         skippedStrat++;
         continue;
         }
      valid.Push(i);
      }
   int N = valid.Length();
   if (skippedTime)
      printf("  %d samples dropped for a missing or non-positive survival time\n",
         skippedTime);
   if (skippedEvent)
      printf("  %d samples dropped for a missing or invalid event code\n",
         skippedEvent);
   if (skippedCov)
      printf("  %d samples dropped for missing covariates\n", skippedCov);
   if (skippedStrat)
      printf("  %d samples dropped for a missing stratum (%s)\n", skippedStrat, (const char*)stratLabel);
   if (N < 10)
      error("Too few complete samples for the Cox model (N=%d).", N);

   int nEvents = 0;
   for (int k = 0; k < N; k++)
      if ((int) (alltraits[valid[k]][event] + 0.5) != 0)
         nEvents++;
   printf("  %d samples with follow-up time, %d events\n", N, nEvents);
   if (nEvents < 5)
      error("Too few events for the Cox model (n=%d).", nEvents);

   // ---- the samples: complete cases, sorted by decreasing time (fixed across PGS) ----
   Matrix allX(idCount, covariateCount);
   Vector allT(idCount);
   IntArray allEv(idCount);
   IntArray validFlag(idCount);
   validFlag.Zero();
   for (int i = 0; i < idCount; i++)
      {
      allT[i] = alltraits[i][time];
      allEv[i] = (alltraits[i][event] == _NAN_) ? 0 : (int) (alltraits[i][event] + 0.5);
      for (int c = 0; c < covariateCount; c++)
         allX[i][c] = alltraits[i][covariates[c]];
      }
   for (int k = 0; k < N; k++)
      validFlag[valid[k]] = 1;
   // strata: one number per distinct value (or combination of values) among the complete cases
   IntArray allStrat(0);
   if (stratCols.Length())
      {
      allStrat.Dimension(idCount);
      allStrat.Zero();
      StringIntHash stratLookup;
      stratLookup.Clear();
      int stratCount = 0;
      for (int k = 0; k < N; k++)
         {
         int i = valid[k];
         String key;
         for (int s = 0; s < stratCols.Length(); s++)
            {
            char buffer[64];
            sprintf(buffer, "%.10g|", alltraits[i][stratCols[s]]);
            key.Add(buffer);
            }
         int id = stratLookup.Integer((const char*) key);
         if (id == -1)
            {
            id = stratCount++;
            stratLookup.SetInteger((const char*) key, id);
            }
         allStrat[i] = id;
         }
      IntArray stratSize(stratCount), stratEvents(stratCount);
      stratSize.Zero();
      stratEvents.Zero();
      for (int k = 0; k < N; k++)
         {
         stratSize[allStrat[valid[k]]]++;
         if (allEv[valid[k]] != 0)
            stratEvents[allStrat[valid[k]]]++;
         }
      int minSize = N, maxSize = 0, noEvent = 0;
      for (int s = 0; s < stratCount; s++)
         {
         if (stratSize[s] < minSize) minSize = stratSize[s];
         if (stratSize[s] > maxSize) maxSize = stratSize[s];
         if (stratEvents[s] == 0) noEvent++;
         }
      printf("  Stratified by %s: %d strata with %d to %d samples each\n", (const char*)stratLabel,
         stratCount, minSize, maxSize);
      if (noEvent)
         printf("  note: %d strata have no events and add nothing to the likelihood\n", noEvent);
      if (stratCount > N / 5)
         printf("WARNING: %d strata for %d samples; --strat is meant for a variable with a few distinct values\n",
            stratCount, N);
      }
   IntArray pheno, evA, sidA;
   Vector tsV, Xs;
   MakeSurvVariable(allX, allT, allEv, allStrat, validFlag, pheno, tsV, evA, sidA, Xs);

   // ---- --readmodel: evaluate pre-specified models instead of scanning ----
   if (modelfile[0])
      {
      Matrix coefM;
      StringArray methodName, pgsName;
      ReadMatrix(modelfile, methodName, pgsName, coefM, coreCount);
      int modelrsCount = coefM.rows;
      int methodCount = coefM.cols;
      if (modelrsCount < 1 || methodCount < 1)
         error("Model file %s does not contain any coefficients.", modelfile);
      if (cvFlag)
         printf("--cv is ignored when --readmodel is given\n");

      Vector cR, cA;
      int missingTerms = Cox_ReadModel(pgsName, coefM, rsNames, risks, pheno, tsV,
         evA, sidA, Xs, covariateCount, cR, cA);
      if (missingTerms)
         printf("  %d model term(s) not present in the score file were skipped\n",
            missingTerms);
      if (covariateList[0])
         printf("  note: a .mod file carries score coefficients only, so "
            "covariates are not applied\n");

      printf("\nSummary of prediction from %s\n", modelfile);
      printf("%20s%12s%8s%8s%10s%10s\n", "TIME", "MODEL", "N", "Events",
         "C", "C_adj");
      for (int mm = 0; mm < methodCount; mm++)
         printf("%20s%12s%8d%8d%10.3lf%10.3lf\n", survTime,
            (const char*)methodName[mm], N, nEvents, cR[mm], cA[mm]);
      printf("\n");
      if (rplotFlag)
         {
         // combined score of every model, in decreasing order of time as the samples are; with
         // --strat the stratum is a column too, so that the baseline hazard is estimated per stratum
         IntArray rsIndex(modelrsCount);
         for (int r = 0; r < modelrsCount; r++)
            rsIndex[r] = rsNames.Find(pgsName[r]);
         bool hasStrata = (sidA.Length() == N);
         int lead = hasStrata ? 3 : 2;
         Matrix scoreTable(N, methodCount + lead);
         for (int k = 0; k < N; k++)
            {
            scoreTable[k][0] = tsV[k];
            scoreTable[k][1] = evA[k];
            if (hasStrata)
               scoreTable[k][2] = sidA[k];
            }
         for (int mm = 0; mm < methodCount; mm++)
            {
            Vector score(N);
            score.Zero();
            bool used = false;
            for (int r = 0; r < modelrsCount; r++)
               {
               int rr = rsIndex[r];
               double bj = coefM[r][mm];
               if (rr == -1 || bj == _NAN_)
                  continue;
               used = true;
               for (int k = 0; k < N; k++)
                  score[k] += risks[rr][pheno[k]] * bj;
               }
            for (int k = 0; k < N; k++)
               scoreTable[k][mm + lead] = used ? score[k] : _NAN_;
            }
         StringArray sheader(methodCount + lead + 1), scol0(N);
         sheader[0] = "ID";
         sheader[1] = survTime;
         sheader[2] = survEvent;
         if (hasStrata)
            sheader[3] = "STRATUM";
         for (int mm = 0; mm < methodCount; mm++)
            sheader[mm + lead + 1] = methodName[mm];
         for (int k = 0; k < N; k++)
            scol0[k] = sampleID[pheno[k]];
         String corrbase(prefix);
         corrbase.Add("_");
         corrbase.Add(survTime);
         corrbase.Add("2");
         corrbase.Add(survEvent);
         String pgsfile = corrbase;
         pgsfile.Add("_combinedPGS.txt");
         WriteTable(pgsfile, scoreTable, sheader, scol0);
         printf("  %d polygenic scores write to %s\n", methodCount, (const char*)pgsfile);
         String msg(256);
         plotCalibrationCox(corrbase, survTime, survEvent, "", msg);
         printf("  %s", (const char*)msg);
         }
      return;
      }

   // ---- PGS correlations, to prune the PGSs of WTSUM (--pruner2), as in --linear / --logistic ----
   if (ldT == _NAN_)
      {
      ldT = 0.2;
      printf("R2 to prune correlated PGSs for WTSUM is missing. --pruner2=%.1lf is set\n", ldT);
      }
   int wordCount = (rsCount - 1) / 64 + 1;
   unsigned long long int ** LD = new unsigned long long int * [rsCount];
   for (int r = 0; r < rsCount; r++)
      {
      LD[r] = new unsigned long long int [wordCount];
      memset(&LD[r][0], 0, wordCount * sizeof(unsigned long long int));
      }
   Vector rsSD(rsCount);
   {
   time_t timer = ::time(NULL);
   printf("Computing PGS correlations starts at %s", asctime(localtime(&timer)));
   fflush(stdout);
   MakeBinaryLD(risks, ldT, LD, rsSD, coreCount);   // standardises risks in place
   LDSummary(LD, rsCount, wordCount, coreCount);
   timer = ::time(NULL);
   printf("Computing PGS correlations ends at %s", asctime(localtime(&timer)));
   fflush(stdout);
   }

   printf("\nFitting %d Cox models (one PGS each)...\n", rsCount);

   Matrix weights;
   IntArray modelIndex;
   CoxModelStats st;
   int status = Cox_BuildModel(risks, pheno, tsV, evA, sidA, Xs, covariateCount, nTop,
      coreCount, LD, weights, modelIndex, st);
   if (status == COX_NO_VALID_PGS)
      error("The Cox model could not be fitted for any of the %d PGSs.", rsCount);
   if (status != COX_OK)
      error("The Cox models could not be built (status %d).", status);
   Vector & coef = st.coef, & se = st.se, & zval = st.zval;
   IntArray & sortedIndex = st.sortedIndex;
   const int sortedCount = sortedIndex.Length();
   int m = st.modelCount;

   // ---- rank and write -------------------------------------------------
   String outfile(prefix);
   outfile.Add("_");
   outfile.Add(survTime);
   outfile.Add("2");
   outfile.Add(survEvent);
   outfile.Add(".coxph");
   FILE * fp = fopen(outfile, "wb");
   if (fp == NULL)
      error("Cannot open %s to write.", (const char*)outfile);
   int annotationCount = annotation.Length();
   fprintf(fp, "PGS\tcoef\texp(coef)\tse(coef)\tz\tPr(>|z|)");
   WriteCoverageHeader(fp);
   if (annotationCount)
      fprintf(fp, "\tAnnotation");
   for (int r = 0; r < sortedCount; r++)
      {
      int rs = sortedIndex[r];
      fprintf(fp, "\n%s\t%.6lf\t%.6lf\t%.6lf\t%.4lf\t%.4G",
         (const char*)rsNames[rs], coef[rs], exp(coef[rs]), se[rs], zval[rs],
         2.0 * ndist(fabs(zval[rs])));
      WriteCoverageFields(fp, rs);
      if (annotationCount)
         fprintf(fp, "\t%s", (const char*)annotation[rs]);
      }
   fprintf(fp, "\n");
   fclose(fp);
   printf("Cox scan summary statistics saved in %s\n", (const char*)outfile);
   String corrbase(prefix);
   corrbase.Add("_");
   corrbase.Add(survTime);
   corrbase.Add("2");
   corrbase.Add(survEvent);
   String corrfile(corrbase);
   corrfile.Add("_PGScorr.txt");
   {
      int topCount = (nTop > 0 ? nTop : 30) < sortedCount ? (nTop > 0 ? nTop : 30) : sortedCount;
      IntArray topIndex(topCount), noStatus(0);
      for (int r = 0; r < topCount; r++)
         topIndex[r] = sortedIndex[r];
      bool caseControlDone = false;
      WritePGSCorrelations(corrbase, risks, rsNames, topIndex, pheno, noStatus, coreCount, caseControlDone, ldT);
      if (rplotFlag)
         {
         String corrMsg(256), corrTitle(survTime);
         corrTitle.Add(",");
         corrTitle.Add(survEvent);
         plotPGSCorr(corrbase, outfile, corrTitle, "", corrMsg);
         printf("  %s", (const char*)corrMsg);
         }
   }

   // brief summary of the strongest scores
   printf("\nTop PGS by Cox Wald Z^2\n");
   printf("%28s%12s%10s%12s", "PGS", "HR", "Z", "P");
   if (annotationCount)
      printf("  %s", "Annotation");
   printf("\n");
   for (int r = 0; r < sortedCount && r < 10; r++)
      {
      int rs = sortedIndex[r];
      printf("%28s%12.4lf%10.3lf%12.3G", (const char*)rsNames[rs], exp(coef[rs]),
         zval[rs], 2.0 * ndist(fabs(zval[rs])));
      if (annotationCount)
         printf("  %s", (const char*)annotation[rs]);
      printf("\n");
      }

   // ---- penalized Cox on the top PGSs, and the WTSUM combination ----
   const char * netName[3] = {"ridge", "elastic net", "lasso"};
   printf("\nPenalized Cox models on the top %d PGSs\n", m);
   for (int a2 = 0; a2 < 3; a2++)
      printf("  %-12s lambda=%.4g  C=%.3lf (adj %.3lf)  non-zero=%d/%d\n",
         netName[a2], st.lambda[a2], st.cIn[a2], st.cInAdj[a2], st.nonzero[a2], m);
   if (st.wsumStatus == 1)
      printf("  WTSUM of %d PGSs: HR=%.3lf Z=%.3f P=%.2G C=%.3lf (adj %.3lf)\n",
         st.wsumTerms, st.hrWsum, st.zWsum, st.pWsum, st.cWsum, st.cWsumAdj);
   else if (st.wsumStatus == 2)
      printf("  WTSUM of %d PGSs: model did not converge\n", st.wsumTerms);
   else
      printf("  WTSUM skipped: no usable PGS coefficients\n");
   printf("  %d PGSs in the models (--ntop %d); WTSUM uses %d of them after pruning at R2 > %.2f (--pruner2)\n",
      m, nTop, st.wsumTerms, ldT);

   // ---- write the model file consumed by --readmodel ----
   StringArray mcol0(m);
      {
      StringArray mheader(6);
      mheader[0] = "PGS_NAME";
      mheader[1] = "BEST";
      mheader[2] = "WTSUM";
      mheader[3] = "RIDGE";
      mheader[4] = "NET";
      mheader[5] = "LASSO";
      // The scan standardised risks in place by rsSD (MakeBinaryLD), but --readmodel
      // reads the risks as they were loaded.  Dividing the coefficients by rsSD expresses
      // them against those risks, so that the reconstruction is identical.  BEST stays 1.0.
      Matrix fileWeights(weights.rows, weights.cols);
      for (int j = 0; j < m; j++)
         {
         mcol0[j] = rsNames[modelIndex[j]];
         double sdj = rsSD[modelIndex[j]];
         if (sdj < 1E-100)
            sdj = 1.0;
         fileWeights[j][0] = weights[j][0];
         for (int c = 1; c < weights.cols; c++)
            fileWeights[j][c] = (weights[j][c] == _NAN_) ? _NAN_ : weights[j][c] / sdj;
         }
      String mfile(prefix);
      mfile.Add("_");
      mfile.Add(survTime);
      mfile.Add("2");
      mfile.Add(survEvent);
      mfile.Add(".mod");
      WriteTable(mfile, fileWeights, mheader, mcol0);
      printf("  PGS model coefficients write to %s\n", (const char*)mfile);
      }

   printf("\nSummary of Cox prediction\n");
   printf("%20s%12s%12s\n", "", "C", "C_adj");
   printf("%20s%12.3lf%12.3lf\n", "best single PGS", st.cBest, st.cBestAdj);
   printf("%20s%12.3lf%12.3lf\n", "WTSUM",
      st.cWsum == _NAN_ ? 0.0 : st.cWsum,
      st.cWsumAdj == _NAN_ ? 0.0 : st.cWsumAdj);
   for (int a2 = 0; a2 < 3; a2++)
      printf("%20s%12.3lf%12.3lf\n", netName[a2], st.cIn[a2], st.cInAdj[a2]);

   if (rplotFlag)
      {
      // combined score of every model for the plots: BEST, WTSUM, RIDGE, NET, LASSO
      // (the columns of weights), against the same risks the models were built on;
      // the event is coded 1 / 0.  Samples are in decreasing order of time.
      String pgsfile(prefix);
      pgsfile.Add("_");
      pgsfile.Add(survTime);
      pgsfile.Add("2");
      pgsfile.Add(survEvent);
      pgsfile.Add("_PGS.txt");
      if (m > 0 && weights.rows == m)
         {
         Matrix scoreTable(N, 7);
         for (int k = 0; k < N; k++)
            {
            scoreTable[k][0] = tsV[k];
            scoreTable[k][1] = evA[k];
            }
         for (int mm = 0; mm < 5; mm++)
            for (int k = 0; k < N; k++)
               {
               double sum = 0.0;
               for (int j = 0; j < m; j++)
                  if (weights[j][mm] != _NAN_)
                     sum += weights[j][mm] * risks[modelIndex[j]][pheno[k]];
               scoreTable[k][mm + 2] = sum;
               }
         StringArray sheader(8), scol0(N);
         sheader[0] = "ID";
         sheader[1] = survTime;
         sheader[2] = survEvent;
         sheader[3] = "BEST";
         sheader[4] = "WTSUM";
         sheader[5] = "RIDGE";
         sheader[6] = "NET";
         sheader[7] = "LASSO";
         for (int k = 0; k < N; k++)
            scol0[k] = sampleID[pheno[k]];
         WriteTable(pgsfile, scoreTable, sheader, scol0);
         printf("  Combined PGSs of %d model(s) write to %s\n", 5, (const char*)pgsfile);
         }
      plotCoxHazard(outfile, corrbase, survTime, pgsfile, corrfile, "");
      }

   // ---- cross-validation: build the models in a train set, test them in a test set ----
   if (cvFlag)
      {
      int iwordCount = (idCount - 1) / 64 + 1;
      unsigned long long int ** testFlags = new unsigned long long int * [repeat];
      for (int run = 0; run < repeat; run++)
         testFlags[run] = new unsigned long long int [iwordCount];
      Matrix cvRaw, cvAdj;
      printf("\nCross-validation (test set = 1/%d of the samples, %d round(s)) start",
         fold, repeat);
      ScanCV(risks, rsNames, allX, allT, allEv, allStrat, validFlag, fold, repeat, nTop, LD,
         cvRaw, cvAdj, testFlags, coreCount);
      printf("done.\n");
      for (int pass = 0; pass < 2; pass++)
         {
         Matrix & res = pass ? cvAdj : cvRaw;
         printf("\nSummary of Cox Cross-validation C-index (%s)\n",
            pass ? "adjusted for covariates" : "unadjusted");
         printf("%20s%8s%8s%8s%10s%10s%10s%10s%10s\n", "TIME", "N", "Events",
            "Rounds", "C_Best", "C_WTSUM", "C_Ridge", "C_Net", "C_Lasso");
         double sum[5] = {0, 0, 0, 0, 0};
         int cnt[5] = {0, 0, 0, 0, 0};
         for (int run = 0; run < repeat; run++)
            for (int mm = 0; mm < 5; mm++)
               if (res[run][mm] != _NAN_)
                  {
                  sum[mm] += res[run][mm];
                  cnt[mm]++;
                  }
         printf("%20s%8d%8d%8d", survTime, N, nEvents, cnt[0]);
         for (int mm = 0; mm < 5; mm++)
            {
            if (cnt[mm])
               printf("%10.3lf", sum[mm] / cnt[mm]);
            else
               printf("%10s", "NA");
            }
         printf("\n");
         }
      for (int run = 0; run < repeat; run++)
         delete []testFlags[run];
      delete []testFlags;
      }
   for (int r = 0; r < rsCount; r++)
      delete []LD[r];
   delete []LD;
   }

