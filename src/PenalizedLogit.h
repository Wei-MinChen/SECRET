////////////////////////////////////////////////////////////////////////
// PenalizedLogit.h
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
// Penalized logistic regression, the binary-trait counterpart of the
// Gaussian solver in Penalized.cpp.  It provides the penalized models
// reported by --logistic.
//
// Objective (glmnet's parameterization, intercept unpenalized):
//
//   min_{b0,b} (1/N) sum_i [ log(1+exp(eta_i)) - y_i eta_i ]
//             + lambda * [ alpha ||b||_1 + (1-alpha)/2 ||b||_2^2 ]
//
//   eta_i = b0 + x_i'b      alpha = 0 ridge, 0.5 elastic net, 1 lasso
//
// X is standardized internally (unit population SD) and coefficients are
// mapped back to the original scale before they are returned.
//
// Fitting uses IRLS (Newton) with an inner coordinate-descent solve of the
// weighted elastic-net least-squares problem, warm-started along the lambda
// path and safeguarded by backtracking on the true penalized objective.
////////////////////////////////////////////////////////////////////////

#ifndef __PENALIZEDLOGIT_H__
#define __PENALIZEDLOGIT_H__

#include "MathMatrix.h"
#include "MathVector.h"
#include "IntArray.h"

class PenalizedLogit
   {
   public:
      PenalizedLogit();
      ~PenalizedLogit();

      // X: N x p design, y: N responses in {0,1}.  Data are copied.
      void SetData(Matrix & X, Vector & y);

      void SetAlpha(double a);              // 0 ridge, 0.5 net, 1 lasso
      void SetLambdas(int n, double ratio);
      void SetFolds(int k);
      void SetSeed(unsigned long s);
      void SetThresh(double tol, int maxOuter);

      int Observations() const { return N; }
      int Dimension() const { return p; }

      // Forget the warm start (used when the data set changes).
      void ResetWarm();

      // Fit one lambda.  beta (length p) and intercept are on the original
      // scale.  Repeated calls reuse the previous solution as a warm start.
      bool Fit(double lambda, Vector & beta, double & intercept);

      // K-fold CV over a lambda path built from the data.  Selects the
      // lambda with the smallest mean deviance and reports the AUC there.
      void CrossValidate(double & lambdaMin, double & cvAUC);

      // In-sample AUC of intercept + X beta.
      double AUC(Vector & beta, double intercept);

      const Vector & LambdaPath() const { return path; }

   private:
      void   Standardize();
      void   BuildPath();
      void   BuildWeighted(double * b, double * eta, double * pr, double * wt,
                 double * z);
      void   SolveQuadratic(double lambda, double * wt, double * z, double * b);
      double Objective(double * b, double lambda);
      int    IRLSFit(double lambda, double tol, int maxOuter);
      void   Unstandardize(double * b, Vector & beta, double & intercept);

      int    N, p;
      Matrix X;                 // original scale
      Matrix Xs;                // standardized
      Vector y;
      Vector xbar, xsd;
      IntArray constant;
      double ybar;

      Vector bCur;              // warm start, bCur[0] = intercept
      Vector path;

      double alpha;
      int    nlambda;
      double minRatio;
      double tol;
      int    maxOuter;
      int    nfold;
      unsigned long seed;
   };

// Rank-based AUC with ties handled by average ranks.  y must be {0,1}.
double LogisticAUC(Vector & pred, Vector & y);

#endif
