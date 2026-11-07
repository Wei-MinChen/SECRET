////////////////////////////////////////////////////////////////////////
// PenalizedLogit.cpp
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
// Penalized logistic regression for --logistic.
// See PenalizedLogit.h for the objective.
////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "PenalizedLogit.h"
#include "Error.h"
#include "QuickIndex.h"

// ---------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------
static double lgSigmoid(double eta)
   {
   if (eta >= 0.0)
      {
      double e = exp(-eta);
      return 1.0 / (1.0 + e);
      }
   double e = exp(eta);
   return e / (1.0 + e);
   }

static double lgSoft(double z, double g)
   {
   if (z > g)
      return z - g;
   if (z < -g)
      return z + g;
   return 0.0;
   }

// stable log(1 + exp(eta))
static double lgLog1pExp(double eta)
   {
   if (eta >= 0.0)
      return eta + log(1.0 + exp(-eta));
   return log(1.0 + exp(eta));
   }

// in-place Cholesky solve; A is destroyed
static bool lgCholesky(double ** A, double * b, int n)
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

static unsigned long long lgState = 88172645463325252ULL;

static void lgSeed(unsigned long long s)
   {
   lgState = s ? s : 88172645463325252ULL;
   }

static unsigned long long lgNext()
   {
   lgState ^= lgState >> 12;
   lgState ^= lgState << 25;
   lgState ^= lgState >> 27;
   return lgState * 2685821657736338717ULL;
   }

static void lgShuffle(IntArray & idx, unsigned long seed)
   {
   lgSeed((unsigned long long) seed);
   for (int i = idx.Length() - 1; i > 0; i--)
      {
      int k = (int) (lgNext() % (unsigned long long) (i + 1));
      int tmp = idx[i];
      idx[i] = idx[k];
      idx[k] = tmp;
      }
   }

// ---------------------------------------------------------------------
// rank-based AUC, ties handled with average ranks
// ---------------------------------------------------------------------
double LogisticAUC(Vector & pred, Vector & y)
   {
   int N = pred.Length();
   if (N != y.Length() || N < 2)
      return 0.5;
   double n1 = 0.0, n0 = 0.0;
   for (int i = 0; i < N; i++)
      if (y[i] > 0.5)
         n1 += 1.0;
      else
         n0 += 1.0;
   if (n1 == 0.0 || n0 == 0.0)
      return 0.5;

   QuickIndex idx;
   idx.Index(pred);

   double sumRanks = 0.0;
   int i = 0;
   while (i < N)
      {
      int j = i;
      while (j + 1 < N && pred[idx[j + 1]] == pred[idx[i]])
         j++;
      double avgRank = (i + 1 + j + 1) / 2.0;   // 1-based average rank
      for (int k = i; k <= j; k++)
         if (y[idx[k]] > 0.5)
            sumRanks += avgRank;
      i = j + 1;
      }
   return (sumRanks - n1 * (n1 + 1.0) / 2.0) / (n1 * n0);
   }

// ---------------------------------------------------------------------
// PenalizedLogit
// ---------------------------------------------------------------------
PenalizedLogit::PenalizedLogit()
   {
   N = p = 0;
   ybar = 0.5;
   alpha = 1.0;
   nlambda = 20;
   minRatio = 1E-4;
   tol = 1E-7;
   maxOuter = 25;
   nfold = 10;
   seed = 20261020UL;
   }

PenalizedLogit::~PenalizedLogit()
   {
   }

void PenalizedLogit::SetAlpha(double a)
   {
   if (a < 0.0 || a > 1.0)
      error("PenalizedLogit: alpha must be between 0 and 1 (got %g).", a);
   alpha = a;
   }

void PenalizedLogit::SetLambdas(int n, double ratio)
   {
   if (n < 1)
      error("PenalizedLogit: nlambda must be positive.");
   nlambda = n;
   minRatio = ratio;
   }

void PenalizedLogit::SetFolds(int k)
   {
   nfold = k < 2 ? 2 : k;
   }

void PenalizedLogit::SetSeed(unsigned long s)
   {
   seed = s;
   }

void PenalizedLogit::SetThresh(double t, int maxIter)
   {
   tol = t;
   if (maxIter > 0)
      maxOuter = maxIter;
   }

void PenalizedLogit::SetData(Matrix & Xin, Vector & yin)
   {
   N = Xin.rows;
   p = Xin.cols;
   if (N != yin.Length())
      error("PenalizedLogit: design has %d rows but response has %d entries.",
         N, yin.Length());
   if (N < 10 || p < 1)
      error("PenalizedLogit: need at least 10 samples and 1 predictor "
         "(N=%d, p=%d).", N, p);

   X.Dimension(N, p);
   for (int i = 0; i < N; i++)
      for (int j = 0; j < p; j++)
         X[i][j] = Xin[i][j];
   y.Dimension(N);
   for (int i = 0; i < N; i++)
      y[i] = yin[i];

   Standardize();
   ResetWarm();
   }

void PenalizedLogit::Standardize()
   {
   xbar.Dimension(p);
   xsd.Dimension(p);
   constant.Dimension(p);
   Xs.Dimension(N, p);
   for (int j = 0; j < p; j++)
      {
      double mean = 0.0;
      for (int i = 0; i < N; i++)
         mean += X[i][j];
      mean /= N;
      double ss = 0.0;
      for (int i = 0; i < N; i++)
         {
         double d = X[i][j] - mean;
         ss += d * d;
         }
      double sd = sqrt(ss / N);
      xbar[j] = mean;
      if (sd < 1E-14)
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
         Xs[i][j] = (X[i][j] - mean) / xsd[j];
      }
   ybar = 0.0;
   for (int i = 0; i < N; i++)
      ybar += y[i];
   ybar /= N;
   }

void PenalizedLogit::ResetWarm()
   {
   bCur.Dimension(p + 1);
   bCur.Zero();
   double m = ybar;
   if (m < 1E-6)
      m = 1E-6;
   if (m > 1.0 - 1E-6)
      m = 1.0 - 1E-6;
   bCur[0] = log(m / (1.0 - m));
   }

void PenalizedLogit::BuildPath()
   {
   double cmax = 0.0;
   for (int j = 0; j < p; j++)
      {
      double s = 0.0;
      for (int i = 0; i < N; i++)
         s += (y[i] - ybar) * Xs[i][j];
      s = fabs(s / N);
      if (s > cmax)
         cmax = s;
      }
   if (cmax < 1E-12)
      cmax = 1E-12;
   double lmax = (alpha > 0.0) ? cmax / alpha : cmax;
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

void PenalizedLogit::BuildWeighted(double * b, double * eta, double * pr,
   double * wt, double * z)
   {
   for (int i = 0; i < N; i++)
      {
      double e = b[0];
      for (int j = 0; j < p; j++)
         e += Xs[i][j] * b[j + 1];
      eta[i] = e;
      pr[i] = lgSigmoid(e);
      wt[i] = pr[i] * (1.0 - pr[i]);
      if (wt[i] < 1E-8)
         wt[i] = 1E-8;
      z[i] = e + (y[i] - pr[i]) / wt[i];
      }
   }

void PenalizedLogit::SolveQuadratic(double lambda, double * wt, double * z,
   double * b)
   {
   int q = p + 1;
   double ** G = new double * [q];
   for (int j = 0; j < q; j++)
      {
      G[j] = new double [q];
      for (int k = 0; k < q; k++)
         G[j][k] = 0.0;
      }
   double * c = new double [q];
   for (int j = 0; j < q; j++)
      c[j] = 0.0;

   for (int i = 0; i < N; i++)
      {
      double wi = wt[i];
      G[0][0] += wi;
      // Fill BOTH triangles as we go.  Filling only one and symmetrising
      // afterwards is wrong here: row 0 (the intercept) is naturally
      // accumulated in the upper triangle while the coefficient block is
      // accumulated in the lower one, so a single symmetrisation pass
      // would overwrite the intercept-covariate cross terms with zeros.
      for (int j = 0; j < p; j++)
         {
         double v = wi * Xs[i][j];
         G[0][j + 1] += v;
         G[j + 1][0] += v;
         }
      for (int j = 0; j < p; j++)
         for (int k = 0; k <= j; k++)
            {
            double v = wi * Xs[i][j] * Xs[i][k];
            G[j + 1][k + 1] += v;
            if (k != j)
               G[k + 1][j + 1] += v;
            }
      c[0] += wi * z[i];
      for (int j = 0; j < p; j++)
         c[j + 1] += wi * Xs[i][j] * z[i];
      }
   for (int j = 0; j < q; j++)
      {
      c[j] /= N;
      for (int k = 0; k < q; k++)
         G[j][k] /= N;
      }

   if (alpha == 0.0)
      {
      // ridge: solve (G + lambda D) b = c, D = diag(0,1,...,1)
      double ** A = new double * [q];
      for (int j = 0; j < q; j++)
         {
         A[j] = new double [q];
         for (int k = 0; k < q; k++)
            A[j][k] = G[j][k];
         if (j > 0)
            A[j][j] += lambda;
         }
      for (int j = 0; j < q; j++)
         b[j] = c[j];
      if (!lgCholesky(A, b, q))
         {
         for (int j = 0; j < q; j++)
            {
            for (int k = 0; k < q; k++)
               A[j][k] = G[j][k];
            if (j > 0)
               A[j][j] += lambda;
            else
               A[j][j] += 1E-8;
            A[j][j] += 1E-10;
            b[j] = c[j];
            }
         if (!lgCholesky(A, b, q))
            error("PenalizedLogit: ridge system not positive definite.");
         }
      for (int j = 0; j < q; j++)
         delete []A[j];
      delete []A;
      }
   else
      {
      double l1 = lambda * alpha;
      double qp = lambda * (1.0 - alpha);
      if (G[0][0] < 1E-12)
         G[0][0] = 1e-12;
      for (int sweep = 0; sweep < 200; sweep++)
         {
         double maxChange = 0.0;
         // intercept, unpenalized
         double num = c[0];
         for (int k = 1; k < q; k++)
            num -= G[0][k] * b[k];
         double nb = num / G[0][0];
         if (fabs(nb - b[0]) > maxChange)
            maxChange = fabs(nb - b[0]);
         b[0] = nb;
         for (int j = 1; j < q; j++)
            {
            if (constant[j - 1])
               {
               b[j] = 0.0;
               continue;
               }
            double rho = c[j];
            for (int k = 0; k < q; k++)
               if (k != j)
                  rho -= G[j][k] * b[k];
            double nbj = lgSoft(rho, l1) / (G[j][j] + qp);
            if (fabs(nbj - b[j]) > maxChange)
               maxChange = fabs(nbj - b[j]);
            b[j] = nbj;
            }
         double scale = 1.0;
         for (int j = 0; j < q; j++)
            if (fabs(b[j]) > scale)
               scale = fabs(b[j]);
         if (maxChange < 1E-9 * scale)
            break;
         }
      }

   for (int j = 0; j < q; j++)
      delete []G[j];
   delete []G;
   delete []c;
   }

double PenalizedLogit::Objective(double * b, double lambda)
   {
   double obj = 0.0;
   for (int i = 0; i < N; i++)
      {
      double e = b[0];
      for (int j = 0; j < p; j++)
         e += Xs[i][j] * b[j + 1];
      obj += lgLog1pExp(e) - y[i] * e;
      }
   obj /= N;
   double l1 = lambda * alpha;
   double qp = 0.5 * lambda * (1.0 - alpha);
   for (int j = 1; j <= p; j++)
      obj += l1 * fabs(b[j]) + qp * b[j] * b[j];
   return obj;
   }

int PenalizedLogit::IRLSFit(double lambda, double tolIn, int maxOuterIn)
   {
   int q = p + 1;
   double * eta = new double [N];
   double * pr = new double [N];
   double * wt = new double [N];
   double * z = new double [N];
   double * b = new double [q];
   double * bOld = new double [q];
   for (int j = 0; j < q; j++)
      b[j] = bCur[j];

   double objOld = Objective(b, lambda);
   int iter = 0;
   for (iter = 1; iter <= maxOuterIn; iter++)
      {
      BuildWeighted(b, eta, pr, wt, z);
      for (int j = 0; j < q; j++)
         bOld[j] = b[j];

      SolveQuadratic(lambda, wt, z, b);

      // Undamped IRLS can raise the objective by a tiny amount on the way,
      // so only damp when the increase is material (a divergence guard);
      // a tight comparison here would halve every step and stall the fit.
      double objNew = Objective(b, lambda);
      double slack = 1E-6 * (1.0 + fabs(objOld));
      int tries = 0;
      while (objNew > objOld + slack && tries < 25)
         {
         for (int j = 0; j < q; j++)
            b[j] = 0.5 * (b[j] + bOld[j]);
         objNew = Objective(b, lambda);
         tries++;
         }

      double delta = 0.0;
      for (int j = 0; j < q; j++)
         if (fabs(b[j] - bOld[j]) > delta)
            delta = fabs(b[j] - bOld[j]);
      objOld = objNew;
      if (delta < tolIn)
         break;
      }

   for (int j = 0; j < q; j++)
      bCur[j] = b[j];

   delete []eta;
   delete []pr;
   delete []wt;
   delete []z;
   delete []b;
   delete []bOld;
   return iter;
   }

void PenalizedLogit::Unstandardize(double * b, Vector & beta, double & intercept)
   {
   beta.Dimension(p);
   double shift = 0.0;
   for (int j = 0; j < p; j++)
      {
      beta[j] = constant[j] ? 0.0 : b[j + 1] / xsd[j];
      shift += xbar[j] * beta[j];
      }
   intercept = b[0] - shift;
   }

bool PenalizedLogit::Fit(double lambda, Vector & beta, double & intercept)
   {
   IRLSFit(lambda, tol, maxOuter);
   int q = p + 1;
   double * b = new double [q];
   for (int j = 0; j < q; j++)
      b[j] = bCur[j];
   Unstandardize(b, beta, intercept);
   delete []b;
   return true;
   }

double PenalizedLogit::AUC(Vector & beta, double intercept)
   {
   Vector pred(N);
   for (int i = 0; i < N; i++)
      {
      double sum = intercept;
      for (int j = 0; j < p; j++)
         sum += X[i][j] * beta[j];
      pred[i] = sum;
      }
   return LogisticAUC(pred, y);
   }

void PenalizedLogit::CrossValidate(double & lambdaMin, double & cvAUCout)
   {
   BuildPath();
   int nl = path.Length();
   Vector cvDev(nl);
   Vector cvAuc(nl);
   cvDev.Zero();
   cvAuc.Zero();

   int K = nfold;
   if (K > N)
      K = N;
   if (K < 2)
      K = 2;

   IntArray idx(N);
   for (int i = 0; i < N; i++)
      idx[i] = i;
   lgShuffle(idx, seed);

   for (int k = 0; k < K; k++)
      {
      int lo = (int) ((double) k * N / K);
      int hi = (int) ((double) (k + 1) * N / K);
      int ntest = hi - lo;
      if (ntest <= 0)
         continue;
      int ntrain = N - ntest;
      if (ntrain < 10)
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

      PenalizedLogit fold;
      fold.SetData(Xt, yt);
      fold.SetAlpha(alpha);
      fold.SetThresh(tol, maxOuter);
      fold.SetFolds(nfold);
      fold.ResetWarm();

      for (int li = 0; li < nl; li++)
         {
         Vector beta;
         double intercept = 0.0;
         fold.Fit(path[li], beta, intercept);

         Vector pred(ntest);
         for (int i = lo; i < hi; i++)
            {
            int ii = idx[i];
            double sum = intercept;
            for (int j = 0; j < p; j++)
               sum += X[ii][j] * beta[j];
            pred[i - lo] = sum;
            }
         Vector yv(ntest);
         for (int i = lo; i < hi; i++)
            yv[i - lo] = y[idx[i]];

         double dev = 0.0;
         for (int i = 0; i < ntest; i++)
            dev += 2.0 * (lgLog1pExp(pred[i]) - yv[i] * pred[i]);
         cvDev[li] += dev / (double) ntest / (double) K;
         cvAuc[li] += LogisticAUC(pred, yv) / (double) K;
         }
      }

   int best = 0;
   for (int li = 1; li < nl; li++)
      if (cvDev[li] < cvDev[best])
         best = li;
   lambdaMin = path[best];
   cvAUCout = cvAuc[best];
   }
