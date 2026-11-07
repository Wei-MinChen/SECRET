////////////////////////////////////////////////////////////////////////
// Logistic.h
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
// Logistic regression for binary traits, the --logistic counterpart of
// --linear.  Each polygenic score is tested against each dichotomous trait
// by Wald test in a logistic model that also carries the --covariate
// columns, so covariate adjustment is done inside the model (binary traits
// cannot be residualized the way quantitative traits are).
//
// The outcome is coded from SECRET's 1/2 convention: 2 (affected) -> 1,
// 1 (unaffected) -> 0.  A positive beta therefore means the score raises
// the odds of being affected.
//
// Risk scores are already scaled to unit variance by ReadRisks, so the
// reported odds ratio is per standard deviation of the score.
////////////////////////////////////////////////////////////////////////

#ifndef __LOGISTIC_H__
#define __LOGISTIC_H__

#include "MathMatrix.h"
#include "MathVector.h"
#include "IntArray.h"
#include "StringArray.h"

// Status codes returned by Logistic_BuildModel.
enum {
   LOGISTIC_OK = 0,
   LOGISTIC_TOO_FEW_SAMPLES = 1,         // N < 10
   LOGISTIC_TOO_FEW_FOR_TERMS = 2,       // N < 5 * (intercept + covariates + PGS)
   LOGISTIC_TOO_FEW_CASES_CONTROLS = 3   // fewer than 5 cases or 5 controls
};

// Optional by-products of Logistic_BuildModel.  They are only needed to report
// the in-sample results of a scan; callers that just want the weights (e.g.,
// cross-validation folds) can pass NULL.
struct LogisticModelStats {
   int nCase, nControl;
   // one-PGS-at-a-time Wald scan, length rsCount; a failed fit has
   // beta == _NAN_, se == 0, zval == 0, pval == 1, odds == _NAN_
   Vector beta, se, zval;
   IntArray sortedIndex;                 // PGS indices by increasing p-value
   int bestIndex;                        // top PGS, -1 if its fit failed
   double aucBest, aucBestRaw;           // best PGS: covariate-adjusted / score only
   int wsumStatus;                       // 0 not fitted, 1 fitted, 2 did not converge
   double wsBeta, wsZ, wsP, aucWsum, aucWsumRaw;
   int wsumTerms;                        // PGSs in WTSUM: the top PGSs left after pruning for correlation
   bool penalizedSkipped;                // too many terms for the sample size
   int penalizedTerms;
   double netAuc[3], netAucRaw[3];   // ridge, elastic net, lasso
};

// Build the prediction models for ONE binary trait; nothing is printed.
//
//   risks          rsCount x (all samples) PGS matrix, as used by the caller
//   pheno          columns of risks used by this model (the samples)
//   Y              0/1 outcome of each of those samples (1 = affected)
//   X              covariates, row-major N x (covariateCount + 1): the covariates
//                  of sample k are X[k*(covariateCount+1) + c]; the last column is
//                  scratch space for the PGS and is overwritten internally
//   LD             binary LD flags from MakeBinaryLD, used to prune the PGSs of WTSUM only
//   nTop           number of top PGSs in the models (already resolved, > 0).  BEST is
//                  the top PGS, RIDGE / NET / LASSO use all nTop PGSs (a PGS whose
//                  scan fit failed is skipped), and WTSUM uses the ones left after
//                  pruning for correlation (stats->wsumTerms of them), as --linear does
//
//   weights        output, modelCount x 5: columns BEST, WTSUM, RIDGE, NET, LASSO.
//                  Row j belongs to PGS modelIndex[j].  Coefficients are expressed
//                  against the `risks` matrix that was passed in, so they can be
//                  fed straight to Logistic_ReadModel with the same matrix.  A
//                  column that could not be fitted (penalized models skipped) is
//                  _NAN_.
//                  The WTSUM weight of a PGS removed by pruning is 0.
//   modelIndex     output, the PGS (row of risks) of each weights row: the top nTop PGSs
//   stats          optional, see above
//
// Returns LOGISTIC_OK, or a LOGISTIC_TOO_FEW_* code, in which case weights has
// no rows and modelIndex is empty.
int Logistic_BuildModel(Matrix & risks, IntArray & pheno, IntArray & Y, Vector & X,
   int covariateCount, unsigned long long int ** LD, int nTop, int coreCount,
   Matrix & weights, IntArray & modelIndex, LogisticModelStats * stats = NULL);

// Prepare, for every binary trait, the samples, outcome and covariate matrix
// used by Logistic_BuildModel / Logistic_ReadModel.  A sample is used for trait
// a when bit a of masks[i] is set and none of its covariates is missing, so a
// train or test subset is selected by clearing the mask bits of the others.
//   pheno[a]  indices of the samples (columns of risks)
//   Xs[a]     row-major N x (covariateCount+1), last column = scratch for the PGS
//   Ys[a]     0/1 outcome (1 = affected)
// pheno, Xs and Ys need allY.cols elements each.
void MakeVariable(Matrix & allX, Matrix & allY, unsigned long long int ** masks,
   IntArray * pheno, Vector * Xs, IntArray * Ys);

// Evaluate models (the columns of coef, e.g. the weights of Logistic_BuildModel
// with pgsName = the PGS names of its rows) on ONE binary trait: AUC of each
// combined score (rawAUC), and the Wald estimate of the score (beta, se) and the
// AUC (adjAUC) in a logistic model that also carries the covariates.  The four
// outputs have one element per column of coef, _NAN_ when not evaluable.  X is
// as made by MakeVariable (its last column is overwritten).  Returns the number
// of model terms not found in rsNames.
int Logistic_ReadModel(StringArray & pgsName, Matrix & coef, StringArray & rsNames,
   Matrix & risks, IntArray & pheno, IntArray & Y, Vector & X, int covariateCount,
   Vector & rawAUC, Vector & adjAUC, Vector & beta, Vector & se);

// Scan every polygenic score against every binary trait (--logistic).
// Ridge / elastic net / lasso are always fitted on the selected top PGSs,
// and their in-sample AUC is reported.  With cvFlag
// (--cv, --fold, --repeat; 0 = default of 10) every model is also built in a
// random train set and its AUC measured in the held-out test set.
void Logistic_Scan(const char * phefile, const char * dictfile,
   const char * rsfile, const char * prefix, const char * covariateList,
   int nTop, double ldT, const char * modelfile, bool rplotFlag,
   bool cvFlag, int fold, int repeat, int coreCount);

#endif
