////////////////////////////////////////////////////////////////////////
// Penalized.h
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
// Penalized (ridge / elastic net / lasso) linear regression for SECRET.
// This replaces the R glmnet() calls with an in-process C++ solver, so
// the penalized models no longer require R to compute the RSQ_Ridge /
// RSQ_Net / RSQ_Lasso columns.
//
// Objective, for the Gaussian family (glmnet's parameterization):
//
//   min_{b0,b}  (1/(2N)) || y - b0 - X b ||^2
//             + lambda * [ alpha ||b||_1 + (1-alpha)/2 ||b||_2^2 ]
//
//   alpha = 0   ridge
//   alpha = 0.5 elastic net
//   alpha = 1   lasso
//
// As in glmnet, X is standardized internally (column means removed,
// columns scaled to unit population standard deviation), the penalty is
// applied to the standardized coefficients, and coefficients are mapped
// back to the original scale before they are reported.
////////////////////////////////////////////////////////////////////////

#ifndef __PENALIZED_H__
#define __PENALIZED_H__

#include "MathMatrix.h"
#include "MathVector.h"
#include "IntArray.h"
#include "StringBasics.h"

class ElasticNet
   {
   public:
      ElasticNet();
      ~ElasticNet();

      // X: N rows (samples) x p columns (risk scores); y: N responses.
      // The data are copied, so the caller may reuse X and y.
      void SetData(Matrix & X, Vector & y);

      void SetAlpha(double a);            // 0 ridge, 0.5 net, 1 lasso
      void SetLambdas(int n, double minRatio);
      void SetThresh(double tol, int maxSweep);
      void SetFolds(int nfold);
      void SetSeed(unsigned long seed);

      int  Observations() const { return N; }
      int  Dimension() const { return p; }
      double Alpha() const { return alpha; }

      // Fit a single lambda on the full data.  On input beta may hold a
      // warm start (length p); on output it holds the coefficients on the
      // original scale, and intercept the fitted intercept.
      bool Fit(double lambda, Vector & beta, double & intercept);

      // Build the lambda sequence from the data and run K-fold cross
      // validation over it.  Returns the lambda minimizing mean squared
      // prediction error (lambdaMin) and that error (cvmMin).
      void CrossValidate(double & lambdaMin, double & cvmMin);

      // In-sample R^2 of lm(y ~ intercept + X beta).
      double RSquare(Vector & beta, double intercept);

      // R^2 of the fitted model applied to new data: squared correlation
      // between the prediction and ynew (what R's lm(y ~ yhat) reports).
      static double PredictionRSquare(Vector & ynew, Vector & yhat);

      const Vector & LambdaPath() const { return path; }
      const Vector & CVError() const { return cvm; }

   private:
      void   Standardize();
      void   BuildGram();
      void   BuildPath();

      int    N, p;
      Matrix X;                 // N x p, original scale
      Matrix Xs;                // N x p, standardized
      Vector y;                 // N

      Vector xbar, xsd;         // column means and population sd
      IntArray constant;        // 1 if the column has no variation
      double ybar;

      Matrix G;                 // p x p, G = Xs'Xs/N  (unit diagonal)
      Vector c;                 // p,      c = Xs'(y-ybar)/N

      Vector path;              // lambda sequence, decreasing
      Vector cvm;               // CV mean squared prediction error

      double alpha;
      int    nlambda;
      double minRatio;
      double tol;
      int    maxSweep;
      int    nfold;
      unsigned long seed;
   };

// ---------------------------------------------------------------------
// High-level drivers
// ---------------------------------------------------------------------

// Compute the 3 x 2 table that SECRET reports for a single trait:
//   rsq[row][0] = in-sample R^2, row = 0 ridge, 1 elastic net, 2 lasso
//   rsq[row][1] = 10-fold CV R^2 when repeat == 0, otherwise the mean
//                 test R^2 over `repeat` random 90/10 holdouts
// coefs is (p+1) x 3: row 0 is the intercept, rows 1..p the coefficients.
void PenalizedNetR2(Matrix & X, Vector & y, int repeat, unsigned long seed,
   Matrix & rsq, Matrix & coefs);

// Drop-in replacement for rplot.cpp's glmnet(): reads <input>_PGS.txt
// (written by WriteTable in OLS.cpp) and writes
//   tmp_<input>_net_rsq.txt
//   tmp_<input>_net_coef.txt
// in the same format the existing code consumes.
int PenalizedRegressionR2(const char * input, const char * traitName,
   int repeat, int coreCount, String & msg);

#endif
