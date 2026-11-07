////////////////////////////////////////////////////////////////////////
// Main.cpp
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

#define VERSION "0.9.3"

#include <math.h>
#include <time.h>
#include "IntArray.h"
#include "StringArray.h"
#include "MathStats.h"
#include "Parameters.h"
#include "Error.h"
#include "StringHash.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

extern const char *gTraitList;
extern const char *gCovarFile;
extern const char *gRpath;
extern bool CheckRInstallation(const char *rpath, String &version, String &exe);
extern bool CheckPs2pdf();
extern String CovarFileColumns(const char *covarfile);
extern void predictionAUC(const char *phefile, const char *dictfile, const char *rsfile,
  bool cvFlag, int fold, int repeat, const char *modelfile,
  const char *prefix, const char *strat, const char *covariateList, bool rplotFlag, const char *rpath, int defaultMaxCoreCount);
extern void RelationshipInference(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *covariateList, IntArray *allrelatives, Vector *allkinships, int coreCount);
extern void GDT_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *strat, bool rplotFlag, int coreCount);
extern void Linear_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *modelfile,
  const char *prefix, const char *covariateList, bool rplotFlag, bool cvFlag,
  int nTop, double ldT, int fold, int repeat, bool nowriteFlag, int coreCount);
extern void TDT_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *strat, int nTop, int coreCount);
extern void Logistic_Scan(const char *phefile, const char *dictfile, const char *rsfile,
  const char *prefix, const char *covariateList, int nTop, double ldT, const char *modelfile, bool rplotFlag,
  bool cvFlag, int fold, int repeat, int coreCount);
extern void runPCA(const char *rsfile, const char *prefix, bool rplotFlag, int coreCount);
extern void survival_Scan(const char *phefile, const char *dictfile, const char *rsfile, const char *prefix, const char *survTime, const char *survEvent, const char *strat, const char *covariateList, bool rplotFlag, int nTop, double ldT, const char *modelfile, bool cvFlag, int fold, int repeat, int coreCount);

void ShowBanner()
{
#ifndef VERSION
   printf("SECRET (c) 10/28/2026 Wei-Min Chen");
#else
   printf("SECRET " VERSION " - (c) 2023-2026 Wei-Min Chen");
#endif
   printf("\nScreen and Evaluate Catalogued Risk scores to Enhance Trait predictions\n");
}

int main(int argc, char **argv)
{
  ShowBanner();
  ParameterList pl;
  String phefile, dictfile, rsfile, modelfile(""), rpath, prefix("secret"), covariateList(""), covarFile("");
  String survList(""), survTime(""), survEvent(""), traitList(""), stratVar("");
  bool aucFlag, rplotFlag, kinshipFlag, pcaFlag, gdtFlag, tdtFlag, linearFlag,
        survFlag, cvFlag, nowriteFlag, logisticFlag;
  aucFlag = rplotFlag = kinshipFlag = pcaFlag = gdtFlag = tdtFlag = linearFlag
        = survFlag = cvFlag = nowriteFlag = logisticFlag = false;
  int coreCount=0, nTop = 0, nfold=0, nrepeat=0;
  double ldT = _NAN_;
  BEGIN_LONG_PARAMETERS(longParameters)
    LONG_PARAMETER_GROUP("Input Files")
         LONG_STRINGPARAMETER("score", &rsfile)
         LONG_STRINGPARAMETER("pheno", &phefile)
         LONG_STRINGPARAMETER("covarfile", &covarFile)
         LONG_STRINGPARAMETER("diction", &dictfile)
         LONG_STRINGPARAMETER("readmodel", &modelfile)
    LONG_PARAMETER_GROUP("PGS Scan")
         LONG_PARAMETER("linear", &linearFlag)
         LONG_PARAMETER("logistic", &logisticFlag)
         LONG_PARAMETER("coxph", &survFlag)
         LONG_PARAMETER("auc", &aucFlag)
         LONG_PARAMETER("tdt", &tdtFlag)
#ifndef VERSION
         LONG_PARAMETER("gdt", &gdtFlag)
#endif
    LONG_PARAMETER_GROUP("Scan Parameter")
         LONG_PARAMETER("cv", &cvFlag)
         LONG_INTPARAMETER("fold", &nfold)
         LONG_INTPARAMETER("repeat", &nrepeat)
         LONG_INTPARAMETER("ntop", &nTop)
         LONG_DOUBLEPARAMETER("pruner2", &ldT)
    LONG_PARAMETER_GROUP("Statistical Model")
         LONG_STRINGPARAMETER("trait", &traitList)
         LONG_STRINGPARAMETER("covariate", &covariateList)
         LONG_STRINGPARAMETER("surv", &survList)
         LONG_STRINGPARAMETER("strat", &stratVar)
#ifndef VERSION
    LONG_PARAMETER_GROUP("Genome Inference")
         LONG_PARAMETER("kinship", &kinshipFlag)
         LONG_PARAMETER("popstruct", &pcaFlag)
#endif
    LONG_PARAMETER_GROUP("Computing Parameter")
         LONG_INTPARAMETER("cpus", &coreCount)
    LONG_PARAMETER_GROUP("Output")
         LONG_PARAMETER("rplot", &rplotFlag)
#ifndef VERSION
         LONG_PARAMETER("nowrite", &nowriteFlag)
#endif
    LONG_PARAMETER_GROUP("Output Parameter")
         LONG_STRINGPARAMETER("prefix", &prefix)
         LONG_STRINGPARAMETER("rpath", &rpath)
  END_LONG_PARAMETERS()

  pl.Add(new LongParameters("", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  if (rsfile.IsEmpty())
    error("Polygenic score files are required. e.g.,\n  secret --score scores.txt --pheno pheno.txt --auc\n\nPlease check the SECRET website at kingrelatedness.com/SECRET/");
  if (phefile.IsEmpty() && aucFlag)
    error("Phenotype files are required. e.g.,\n  secret --score scores.txt --pheno pheno.txt --auc\n\nPlease check the SECRET website at kingrelatedness.com/SECRET/");
  // --surv time,event: exactly two trait names separated by a comma
  if (!survList.IsEmpty() && !survFlag)
    printf("Note: --surv only applies to --coxph and is ignored.\n");
  if (survFlag) {
    int commas = 0;
    for (int i = 0; i < survList.Length(); i++)
      if (survList[i] == ',') commas++;
    StringArray survNames;
    survNames.AddTokens(survList, ',');
    if (commas != 1 || survNames.Length() != 2) {
      if (survList.IsEmpty())
        error("--coxph requires --surv time,event. e.g.,\n  secret --score scores.txt --pheno pheno.txt --coxph --surv time,event");
      else
        error("--surv takes exactly two trait names separated by a comma, the survival time and the event. e.g.,\n  --surv time,event");
    }
    survTime = survNames[0];
    survEvent = survNames[1];
  }
  // --strat: the trait that stratifies the results of --auc, --tdt and --gdt (one trait),
  // and the stratum variable(s) of --coxph
  if (!stratVar.IsEmpty() && !aucFlag && !survFlag && !tdtFlag && !gdtFlag)
    printf("Note: --strat only applies to --auc, --coxph and --tdt, and is ignored.\n");
  if (!stratVar.IsEmpty() && (aucFlag || tdtFlag || gdtFlag) && stratVar.Find(',') > -1)
    error("--strat takes a single trait for --auc and --tdt; several variables separated by commas are for --coxph only.");
  if (!traitList.IsEmpty() && aucFlag && !stratVar.IsEmpty()) {   // the --strat trait must be analyzed too
    StringArray listed;
    listed.AddTokens(traitList, ',');
    if (listed.Find(stratVar) == -1) {
      traitList.Add(',');
      traitList.Add(stratVar);
      printf("Note: the --strat trait %s is added to the traits given by --trait.\n", (const char*)stratVar);
    }
  }
  gTraitList = (const char*)traitList;   // read by MakeTraits; traitList is not modified after this point
  if (!covarFile.IsEmpty()) {
    if (phefile.IsEmpty()) {
      printf("Note: --covarfile is ignored because --pheno is not given.\n");
      covarFile = "";
    }
    else {
      // without --covariate, every column of the covariate file is a covariate
      if (covariateList.IsEmpty()) covariateList = CovarFileColumns(covarFile);
      gCovarFile = (const char*)covarFile;   // read by SECRET_Read; covarFile is not modified after this point
    }
  }
  // R for the figures: --rpath, or the R found in the PATH.  Every plot reads gRpath.
  gRpath = (const char*)rpath;   // rpath is not modified after this point
  if (rplotFlag || pcaFlag) {
    String rVersion, rExe;
    if (CheckRInstallation(rpath, rVersion, rExe)) {
      printf("R for the figures: %s (%s)\n", (const char*)rVersion, rpath.IsEmpty() ? "from the PATH" : (const char*)rExe);
      // R writes PostScript; ps2pdf (Ghostscript) makes the PDF
      if (!CheckPs2pdf())
        printf("WARNING: ps2pdf (Ghostscript) is not working, so the figures are kept as .ps files and not made into .pdf.\n"
               "  Install Ghostscript and put ps2pdf in the PATH to get PDF files.\n");
    }
    else if (!rpath.IsEmpty())
      error("--rpath %s is not a working R.\n  Give the R executable (e.g. /opt/R/4.3.1/bin/R) or the directory that contains it.", (const char*)rpath);
    else
      printf("WARNING: R was not found in the PATH, so the figures cannot be made.\n  Install R, or give its location with --rpath (e.g. --rpath /opt/R/4.3.1/bin/R).\n");
  }
  else if (!rpath.IsEmpty())
    printf("Note: --rpath is ignored because --rplot is not given.\n");
  time_t timer;
  struct tm *tblock;
  timer = time(NULL);
  tblock = localtime(&timer);
  printf("SECRET starts at %s", asctime(tblock));

  int defaultMaxCoreCount = 1;
#ifdef _OPENMP
  int availableCoreCount = omp_get_max_threads();
  if (coreCount > 0)
    defaultMaxCoreCount = coreCount > availableCoreCount? availableCoreCount: coreCount;
  else
    defaultMaxCoreCount = availableCoreCount==1? 1: availableCoreCount / 2;
#endif
  if (aucFlag)    // the --rplot figures are made there too, for each trait
    predictionAUC(phefile, dictfile, rsfile, cvFlag, nfold, nrepeat, modelfile, prefix, stratVar, covariateList, rplotFlag, rpath, defaultMaxCoreCount);
  if (gdtFlag)
    GDT_Scan(phefile, dictfile, rsfile, prefix, stratVar, rplotFlag, defaultMaxCoreCount);
  if (tdtFlag)
    TDT_Scan(phefile, dictfile, rsfile, prefix, stratVar, nTop, defaultMaxCoreCount);
  if (linearFlag)
    Linear_Scan(phefile, dictfile, rsfile, modelfile, prefix, covariateList,
    rplotFlag, cvFlag, nTop, ldT, nfold, nrepeat, nowriteFlag, defaultMaxCoreCount);
  if (logisticFlag)
    Logistic_Scan(phefile, dictfile, rsfile, prefix, covariateList, nTop, ldT, modelfile, rplotFlag,
      cvFlag, nfold, nrepeat, defaultMaxCoreCount);
  if (survFlag)
    survival_Scan(phefile, dictfile, rsfile, prefix, survTime, survEvent, stratVar, covariateList, rplotFlag, nTop, ldT, modelfile,
      cvFlag, nfold, nrepeat, defaultMaxCoreCount);
  if (pcaFlag)
    runPCA(rsfile, prefix, rplotFlag, defaultMaxCoreCount);
  if (kinshipFlag)
    RelationshipInference(phefile, dictfile, rsfile, prefix, covariateList, NULL, NULL, defaultMaxCoreCount);

  timer = time(NULL);
  tblock = localtime(&timer);
  printf("SECRET ends at %s\n", asctime(tblock));
}


