////////////////////////////////////////////////////////////////////////
// Penalized.cpp
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
// Oct 22, 2026
//
// In-process ridge / elastic net / lasso solver, replacing the R
// glmnet() calls used by the penalized models.  See Penalized.h.
//
// Algorithm
//   * X is standardized once per fit (column means removed, unit
//     population SD), giving a gram matrix G with unit diagonal.
//   * alpha > 0: cyclic coordinate descent with soft-thresholding,
//     warm-started along the lambda path.
//   * alpha = 0: (G + lambda I) b = c solved by Cholesky.
//   * lambda sequence: log-spaced, anchored at max|c|/alpha for
//     alpha > 0 (glmnet's convention).
//   * K-fold CV over that sequence picks lambda.min; the model is then
//     refit on all samples at lambda.min.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Penalized.h"
#include "Constant.h"
#include "Error.h"
#include "DataInput.h"
#include "StringArray.h"

// ---------------------------------------------------------------------
// deterministic xorshift64* generator, so fold assignment and the
// repeated holdouts are reproducible from a seed
// ---------------------------------------------------------------------
static unsigned long long pnState = 88172645463325252ULL;

static void pnSeed(unsigned long long s)
   {
   pnState = s ? s : 88172645463325252ULL;
   }

static unsigned long long pnNext()
   {
   pnState ^= pnState >> 12;
   pnState ^= pnState << 25;
   pnState ^= pnState >> 27;
   return pnState * 2685821657736338717ULL;
   }

static void shuffleIndex(IntArray & idx, unsigned long seed)
   {
   pnSeed((unsigned long long) seed);
   for (int i = idx.Length() - 1; i > 0; i--)
      {
      int k = (int) (pnNext() % (unsigned long long) (i + 1));
      int tmp = idx[i];
      idx[i] = idx[k];
      idx[k] = tmp;
      }
   }

static double roundTo(double x, int digits)
   {
   double m = pow(10.0, digits);
   if (x >= 0)
      return floor(x * m + 0.5) / m;
   return -floor(-x * m + 0.5) / m;
   }

// R's var(): denominator n-1
static double sampleVariance(Vector & v)
   {
   int n = v.Length();
   if (n < 2)
      return 0.0;
   double mean = 0.0;
   for (int i = 0; i < n; i++)
      mean += v[i];
   mean /= n;
   double ss = 0.0;
   for (int i = 0; i < n; i++)
      {
      double d = v[i] - mean;
      ss += d * d;
      }
   return ss / (n - 1);
   }

// squared correlation, i.e. the R^2 of lm(y ~ yhat)
static double correlationSquare(Vector & a, Vector & b)
   {
   int n = a.Length();
   if (n != b.Length() || n < 2)
      return 0.0;
   double ma = 0.0, mb = 0.0;
   for (int i = 0; i < n; i++)
      {
      ma += a[i];
      mb += b[i];
      }
   ma /= n;
   mb /= n;
   double sa = 0.0, sb = 0.0, sab = 0.0;
   for (int i = 0; i < n; i++)
      {
      double da = a[i] - ma;
      double db = b[i] - mb;
      sa += da * da;
      sb += db * db;
      sab += da * db;
      }
   if (sa <= 0 || sb <= 0)
      return 0.0;
   return sab * sab / (sa * sb);
   }

// ---------------------------------------------------------------------
// symmetric positive definite solve: in-place Cholesky, A is destroyed
// ---------------------------------------------------------------------
static bool choleskySolve(double ** A, double * b, int n)
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
   for (int i = 0; i < n; i++)          // forward substitution
      {
      double sum = b[i];
      for (int k = 0; k < i; k++)
         sum -= A[i][k] * b[k];
      b[i] = sum / A[i][i];
      }
   for (int i = n - 1; i >= 0; i--)     // back substitution
      {
      double sum = b[i];
      for (int k = i + 1; k < n; k++)
         sum -= A[k][i] * b[k];
      b[i] = sum / A[i][i];
      }
   return true;
   }

// G = Xs'Xs/N and c = Xs'(y - ybar)/N
static void buildGram(Matrix & Xs, Vector & yc, Matrix & G, Vector & c)
   {
   int N = Xs.rows;
   int p = Xs.cols;
   G.Dimension(p, p);
   for (int j = 0; j < p; j++)
      {
      for (int k = 0; k <= j; k++)
         {
         double sum = 0.0;
         for (int i = 0; i < N; i++)
            sum += Xs[i][j] * Xs[i][k];
         sum /= N;
         G[j][k] = sum;
         G[k][j] = sum;
         }
      }
   c.Dimension(p);
   for (int j = 0; j < p; j++)
      {
      double sum = 0.0;
      for (int i = 0; i < N; i++)
         sum += Xs[i][j] * yc[i];
      c[j] = sum / N;
      }
   }

// standardize columns of Xin into Xout; also returns ybar for y
static void standardizeDesign(Matrix & Xin, Vector & y, Matrix & Xout,
   Vector & xbar, Vector & xsd, IntArray & constant, double & ybar)
   {
   int N = Xin.rows;
   int p = Xin.cols;
   xbar.Dimension(p);
   xsd.Dimension(p);
   constant.Dimension(p);
   Xout.Dimension(N, p);
   for (int j = 0; j < p; j++)
      {
      double mean = 0.0;
      for (int i = 0; i < N; i++)
         mean += Xin[i][j];
      mean /= N;
      double ss = 0.0;
      for (int i = 0; i < N; i++)
         {
         double d = Xin[i][j] - mean;
         ss += d * d;
         }
      double sd = sqrt(ss / N);            // population sd, as glmnet uses
      xbar[j] = mean;
      if (sd < 1E-14)                      // constant column
         {
         constant[j] = 1;
         xsd[j] = 1.0;
         }
      else
         {
         constant[j] = 0;
         xsd[j] = sd;
         }
      for (int i = 0; i < N; i++)
         Xout[i][j] = (Xin[i][j] - mean) / xsd[j];
      }
   ybar = 0.0;
   for (int i = 0; i < N; i++)
      ybar += y[i];
   ybar /= N;
   }

static void unstandardize(Vector & b, Vector & xbar, Vector & xsd,
   IntArray & constant, double ybar, Vector & beta, double & intercept)
   {
   int p = b.Length();
   beta.Dimension(p);
   double shift = 0.0;
   for (int j = 0; j < p; j++)
      {
      beta[j] = constant[j] ? 0.0 : b[j] / xsd[j];
      shift += xbar[j] * beta[j];
      }
   intercept = ybar - shift;
   }

// soft-threshold  S(z, g) = sign(z) max(|z| - g, 0)
static double softThreshold(double z, double g)
   {
   if (z > g)
      return z - g;
   if (z < -g)
      return z + g;
   return 0.0;
   }

// cyclic coordinate descent; b is both warm start and result
static int coordinateDescent(Matrix & G, Vector & c, IntArray & constant,
   double lambda, double alpha, double tol, int maxSweep, Vector & b)
   {
   int p = c.Length();
   double l1 = lambda * alpha;
   double quad = 1.0 + lambda * (1.0 - alpha);
   int sweep = 0;
   for (sweep = 1; sweep <= maxSweep; sweep++)
      {
      double maxChange = 0.0;
      for (int j = 0; j < p; j++)
         {
         if (constant[j])
            {
            b[j] = 0.0;
            continue;
            }
         double rho = c[j];
         Vector & Gj = G[j];
         for (int k = 0; k < p; k++)
            if (k != j)
               rho -= Gj[k] * b[k];
         double next = softThreshold(rho, l1) / quad;
         double delta = fabs(next - b[j]);
         if (delta > maxChange)
            maxChange = delta;
         b[j] = next;
         }
      double scale = 1.0;
      for (int j = 0; j < p; j++)
         if (fabs(b[j]) > scale)
            scale = fabs(b[j]);
      if (maxChange < tol * scale)
         break;
      }
   return sweep;
   }

// (G + lambda I) b = c
static void ridgeSolve(Matrix & G, Vector & c, IntArray & constant,
   double lambda, Vector & b)
   {
   int p = c.Length();
   double ** A = new double * [p];
   for (int i = 0; i < p; i++)
      A[i] = new double [p];
   double * rhs = new double [p];
   for (int i = 0; i < p; i++)
      {
      for (int j = 0; j < p; j++)
         A[i][j] = G[i][j];
      A[i][i] += lambda;
      rhs[i] = c[i];
      }
   bool ok = choleskySolve(A, rhs, p);
   if (!ok)                              // not positive definite: jitter
      {
      for (int i = 0; i < p; i++)
         {
         for (int j = 0; j < p; j++)
            A[i][j] = G[i][j];
         A[i][i] += lambda + 1E-10 * (fabs(G[i][i]) + 1.0);
         rhs[i] = c[i];
         }
      ok = choleskySolve(A, rhs, p);
      }
   if (!ok)
      error("ElasticNet: ridge system not positive definite (lambda=%g).", lambda);
   for (int i = 0; i < p; i++)
      b[i] = constant[i] ? 0.0 : rhs[i];
   for (int i = 0; i < p; i++)
      delete []A[i];
   delete []A;
   delete []rhs;
   }

// ---------------------------------------------------------------------
// ElasticNet
// ---------------------------------------------------------------------
ElasticNet::ElasticNet()
   {
   N = p = 0;
   ybar = 0.0;
   alpha = 1.0;
   nlambda = 100;
   minRatio = 1E-4;
   tol = 1E-7;
   maxSweep = 10000;
   nfold = 10;
   seed = 20261020UL;
   }

ElasticNet::~ElasticNet()
   {
   }

void ElasticNet::SetAlpha(double a)
   {
   if (a < 0.0 || a > 1.0)
      error("ElasticNet: alpha must be between 0 and 1 (got %g).", a);
   alpha = a;
   }

void ElasticNet::SetLambdas(int n, double ratio)
   {
   if (n < 1)
      error("ElasticNet: nlambda must be positive.");
   nlambda = n;
   minRatio = ratio;
   }

void ElasticNet::SetThresh(double t, int maxIter)
   {
   tol = t;
   if (maxIter > 0)
      maxSweep = maxIter;
   }

void ElasticNet::SetFolds(int k)
   {
   nfold = k < 2 ? 2 : k;
   }

void ElasticNet::SetSeed(unsigned long s)
   {
   seed = s;
   }

void ElasticNet::SetData(Matrix & Xin, Vector & yin)
   {
   N = Xin.rows;
   p = Xin.cols;
   if (N != yin.Length())
      error("ElasticNet: design has %d rows but response has %d entries.",
         N, yin.Length());
   if (N < 3 || p < 1)
      error("ElasticNet: need at least 3 samples and 1 predictor (N=%d, p=%d).",
         N, p);
   if (p > 8000)
      error("ElasticNet: %d predictors requested - the gram matrix alone would "
         "need %d MB. Reduce --ntop.", p, (int) (8.0 * p * (double) p / 1E6));

   X.Dimension(N, p);
   for (int i = 0; i < N; i++)
      for (int j = 0; j < p; j++)
         X[i][j] = Xin[i][j];
   y.Dimension(N);
   for (int i = 0; i < N; i++)
      y[i] = yin[i];

   Standardize();
   BuildGram();
   }

void ElasticNet::Standardize()
   {
   standardizeDesign(X, y, Xs, xbar, xsd, constant, ybar);
   }

void ElasticNet::BuildGram()
   {
   Vector yc(N);
   for (int i = 0; i < N; i++)
      yc[i] = y[i] - ybar;
   buildGram(Xs, yc, G, c);
   }

void ElasticNet::BuildPath()
   {
   double cmax = 0.0;
   for (int j = 0; j < p; j++)
      if (fabs(c[j]) > cmax)
         cmax = fabs(c[j]);
   if (cmax < 1E-12)
      cmax = 1E-12;

   // glmnet anchors the path at max|c|/alpha; for alpha = 0 that
   // expression diverges, so the alpha = 1 anchor is used instead
   double lmax = alpha > 0.0 ? cmax / alpha : cmax;
   double lmin = lmax * minRatio;
   if (lmin <= 0.0)
      lmin = lmax * 1E-4;

   path.Dimension(nlambda);
   if (nlambda == 1)
      {
      path[0] = lmax;
      return;
      }
   double lr = log(lmin / lmax);
   for (int i = 0; i < nlambda; i++)
      path[i] = exp(log(lmax) + lr * i / (nlambda - 1));
   }

bool ElasticNet::Fit(double lambda, Vector & beta, double & intercept)
   {
   Vector b(p);
   b.Zero();
   if (alpha == 0.0)
      ridgeSolve(G, c, constant, lambda, b);
   else
      coordinateDescent(G, c, constant, lambda, alpha, tol, maxSweep, b);
   unstandardize(b, xbar, xsd, constant, ybar, beta, intercept);
   return true;
   }

double ElasticNet::RSquare(Vector & beta, double intercept)
   {
   Vector yhat(N);
   for (int i = 0; i < N; i++)
      {
      double sum = intercept;
      for (int j = 0; j < p; j++)
         sum += X[i][j] * beta[j];
      yhat[i] = sum;
      }
   return correlationSquare(y, yhat);
   }

double ElasticNet::PredictionRSquare(Vector & ynew, Vector & yhat)
   {
   return correlationSquare(ynew, yhat);
   }

void ElasticNet::CrossValidate(double & lambdaMin, double & cvmMin)
   {
   BuildPath();
   int nl = path.Length();
   cvm.Dimension(nl);
   cvm.Zero();

   int K = nfold;
   if (K > N)
      K = N;
   if (K < 2)
      K = 2;

   IntArray idx(N);
   for (int i = 0; i < N; i++)
      idx[i] = i;
   shuffleIndex(idx, seed);

   for (int k = 0; k < K; k++)
      {
      int lo = (int) ((double) k * N / K);
      int hi = (int) ((double) (k + 1) * N / K);
      int ntest = hi - lo;
      if (ntest <= 0)
         continue;
      int ntrain = N - ntest;
      if (ntrain < 2)
         continue;

      Matrix Xt(ntrain, p);
      Vector yt(ntrain);
      int row = 0;
      for (int i = 0; i < N; i++)
         {
         if (i >= lo && i < hi)
            continue;
         int ii = idx[i];
         for (int j = 0; j < p; j++)
            Xt[row][j] = X[ii][j];
         yt[row] = y[ii];
         row++;
         }

      Matrix Xs_t;
      Vector xbar_t, xsd_t;
      IntArray constant_t;
      double ybar_t;
      standardizeDesign(Xt, yt, Xs_t, xbar_t, xsd_t, constant_t, ybar_t);

      Vector yc(ntrain);
      for (int i = 0; i < ntrain; i++)
         yc[i] = yt[i] - ybar_t;
      Matrix Gt;
      Vector ct;
      buildGram(Xs_t, yc, Gt, ct);

      Vector b(p);
      b.Zero();
      Vector beta(p);
      double intercept = 0.0;
      for (int li = 0; li < nl; li++)
         {
         double lambda = path[li];
         if (alpha == 0.0)
            ridgeSolve(Gt, ct, constant_t, lambda, b);
         else
            coordinateDescent(Gt, ct, constant_t, lambda, alpha, tol, maxSweep, b);
         unstandardize(b, xbar_t, xsd_t, constant_t, ybar_t, beta, intercept);

         double sse = 0.0;
         for (int i = lo; i < hi; i++)
            {
            int ii = idx[i];
            double pred = intercept;
            for (int j = 0; j < p; j++)
               pred += X[ii][j] * beta[j];
            double resid = y[ii] - pred;
            sse += resid * resid;
            }
         cvm[li] += sse / (double) ntest / (double) K;
         }
      }

   int best = 0;
   for (int li = 1; li < nl; li++)
      if (cvm[li] < cvm[best])
         best = li;
   lambdaMin = path[best];
   cvmMin = cvm[best];
   }

// ---------------------------------------------------------------------
// high-level drivers
// ---------------------------------------------------------------------

// `repeat` random 90/10 holdouts: CV-select lambda on the training part,
// then the squared correlation between prediction and response on the
// held-out part, averaged over repeats.
static double repeatedHoldout(Matrix & X, Vector & y, double alpha,
   int repeat, unsigned long seed, int nfold)
   {
   int N = X.rows;
   int p = X.cols;
   int ntrain = (int) (0.9 * N);
   if (ntrain < 3 || N - ntrain < 2)
      return 0.0;

   double total = 0.0;
   int used = 0;
   for (int r = 0; r < repeat; r++)
      {
      IntArray idx(N);
      for (int i = 0; i < N; i++)
         idx[i] = i;
      shuffleIndex(idx, seed + 7919UL * (unsigned long) (r + 1));

      Matrix Xt(ntrain, p);
      Vector yt(ntrain);
      for (int i = 0; i < ntrain; i++)
         {
         int ii = idx[i];
         for (int j = 0; j < p; j++)
            Xt[i][j] = X[ii][j];
         yt[i] = y[ii];
         }
      int ntest = N - ntrain;
      Matrix Xv(ntest, p);
      Vector yv(ntest);
      for (int i = 0; i < ntest; i++)
         {
         int ii = idx[ntrain + i];
         for (int j = 0; j < p; j++)
            Xv[i][j] = X[ii][j];
         yv[i] = y[ii];
         }

      ElasticNet en;
      en.SetData(Xt, yt);
      en.SetAlpha(alpha);
      en.SetFolds(nfold);
      en.SetSeed(seed + 104729UL * (unsigned long) (r + 1));
      double lambda, cvmMin;
      en.CrossValidate(lambda, cvmMin);

      Vector beta;
      double intercept = 0.0;
      en.Fit(lambda, beta, intercept);

      Vector yhat(ntest);
      for (int i = 0; i < ntest; i++)
         {
         double sum = intercept;
         for (int j = 0; j < p; j++)
            sum += Xv[i][j] * beta[j];
         yhat[i] = sum;
         }
      total += correlationSquare(yv, yhat);
      used++;
      }
   return used ? total / used : 0.0;
   }

void PenalizedNetR2(Matrix & X, Vector & y, int repeat, unsigned long seed,
   Matrix & rsq, Matrix & coefs)
   {
   double alphas[3] = {0.0, 0.5, 1.0};
   int p = X.cols;

   rsq.Dimension(3, 2);
   rsq.Zero();
   coefs.Dimension(p + 1, 3);
   coefs.Zero();

   double vy = sampleVariance(y);

   for (int a = 0; a < 3; a++)
      {
      ElasticNet en;
      en.SetData(X, y);
      en.SetAlpha(alphas[a]);
      en.SetSeed(seed + 1000UL * (unsigned long) a);

      double lambda, cvmMin;
      en.CrossValidate(lambda, cvmMin);

      Vector beta;
      double intercept = 0.0;
      en.Fit(lambda, beta, intercept);

      // R rounds the in-sample R^2 to 3 digits inside the script
      rsq[a][0] = roundTo(en.RSquare(beta, intercept), 3);

      if (repeat > 0)
         rsq[a][1] = roundTo(repeatedHoldout(X, y, alphas[a], repeat,
            seed + 5000UL * (unsigned long) a, 10), 5);
      else
         rsq[a][1] = roundTo(vy > 0.0 ? 1.0 - cvmMin / vy : 0.0, 5);

      coefs[0][a] = intercept;
      for (int j = 0; j < p; j++)
         coefs[j + 1][a] = beta[j];
      }
   }

int PenalizedRegressionR2(const char * input, const char * traitName,
   int repeat, int coreCount, String & msg)
   {
   String pgsfile(input);
   pgsfile.Add("_PGS.txt");

   StringArray colNames, rowNames;
   Matrix M;
   ReadMatrix(pgsfile, colNames, rowNames, M, coreCount);

   int cols = M.cols;
   if (cols < 3)
      error("PenalizedRegression: %s needs trait, at least one PGS and WTSUM "
         "columns (found %d).", (const char*)pgsfile, cols);
   int p = cols - 2;

   // complete-case: drop samples with a missing sentinel anywhere
   IntArray valid(0);
   for (int i = 0; i < M.rows; i++)
      {
      bool ok = true;
      for (int j = 0; j < cols; j++)
         if (M[i][j] == _NAN_)
            {
            ok = false;
            break;
            }
      if (ok)
         valid.Push(i);
      }
   int N = valid.Length();
   int dropped = M.rows - N;
   if (dropped > 0)
      printf("  %d of %d samples dropped from penalized regression "
         "for missing values.\n", dropped, M.rows);
   if (N < 3 || p < 1)
      error("PenalizedRegression: too few complete samples for %s (N=%d, p=%d).",
         (const char*)pgsfile, N, p);

   Matrix X(N, p);
   Vector y(N);
   for (int i = 0; i < N; i++)
      {
      int ii = valid[i];
      y[i] = M[ii][0];
      for (int j = 0; j < p; j++)
         X[i][j] = M[ii][j + 1];
      }

   Matrix rsq, coefs;
   PenalizedNetR2(X, y, repeat, 20261020UL, rsq, coefs);

   String rsqfile("tmp_");
   rsqfile.Add(input);
   rsqfile.Add("_net_rsq.txt");
   FILE * fp = fopen((const char*)rsqfile, "wt");
   if (fp == NULL)
      error("Cannot open %s to write.", (const char*)rsqfile);
   for (int a = 0; a < 3; a++)
      fprintf(fp, "%.5f\t%.5f\n", rsq[a][0], rsq[a][1]);
   fclose(fp);

   String coeffile("tmp_");
   coeffile.Add(input);
   coeffile.Add("_net_coef.txt");
   fp = fopen((const char*)coeffile, "wt");
   if (fp == NULL)
      error("Cannot open %s to write.", (const char*)coeffile);
   fprintf(fp, "Ridge\tNet\tLasso\n");
   for (int j = 0; j <= p; j++)
      fprintf(fp, "%.6f\t%.6f\t%.6f\n", coefs[j][0], coefs[j][1], coefs[j][2]);
   fclose(fp);

   msg.Copy("Penalized regression (ridge / elastic net / lasso) fitted in C++\n");
   if (traitName != NULL && traitName[0] != '\0')
      printf("  Trait %s: ridge/net/lasso fitted on %d samples x %d scores.\n",
         traitName, N, p);
   return 0;
   }
