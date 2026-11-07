#include <stdio.h>
#include "StringArray.h"
#include "MathStats.h"

extern const char *gTraitList;   // --trait: comma-separated traits to analyze; "" means all traits
extern const char *gCovarFile;   // --covarfile: comma-separated file names, covariates kept apart from the traits; "" means none
extern bool gHasCoverage;        // the dictionary (--diction) file carries SNPused and/or Coverage information
extern Vector gSNPused, gCoverage;   // per PGS, in the order of rsNames; _NAN_ when unknown
void WriteCoverageHeader(FILE *fp);           // writes "\tSNPused\tCoverage" when coverage information is available
void WriteCoverageFields(FILE *fp, int rs);   // the two values of PGS number rs (NA when unknown)
int FormatCoverageHeader(char *buffer);       // the same text as WriteCoverageHeader, put in buffer; returns its length
int FormatCoverageFields(char *buffer, int rs);   // the same text as WriteCoverageFields, put in buffer; returns its length

// What a top-PGS heatmap (WriteTopHeatmapR in rplot.cpp) shows in its cells
#define TOPHEAT_RSQ  0   // --linear:   Rsq
#define TOPHEAT_OR   1   // --logistic: odds ratio
#define TOPHEAT_HR   2   // --coxph:    hazard ratio
String CovarFileColumns(const char *covarfile);   // names of the covariate columns in --covarfile (one or more comma-separated files)
void SplitFileList(const char *list, StringArray &files);   // comma-separated file names -> one name per element
void ReadScoreMatrix(const char *rsfile, StringArray &rsNames, StringArray &rowNames, Matrix &M, int coreCount);   // --score: one or more comma-separated files; rows are samples, columns are PGS

void ReadLines(const char *filename, StringArray &lines, bool printFlag);
void ReadTraits(const char *phefile, StringArray &alltraitNames, StringArray &familyID, StringArray &sampleID, Matrix &alltraits, int coreCount);
void ReadTraits4(const char *phefile, StringArray &alltraitNames, StringArray &familyID, StringArray &sampleID, StringArray &FAID, StringArray &MOID, Matrix &alltraits, int coreCount);
void ReadRisks(const char *rsfile, StringArray &rsNames, StringArray &familyID, StringArray &sampleID, Matrix &risks, IntArray &rsKey, int coreCount);
void ReadMatrix(const char *phefile, StringArray &colNames, Matrix &M, int coreCount);
void ReadMatrix(const char *phefile, StringArray &colNames, StringArray &rowNames, Matrix &M, int coreCount);
void ReadKeys(const char *filename, StringArray &rsNames, StringArray &annotation);
void MakeCovariates(const char *covariateList, StringArray &alltraitNames, StringArray &covList, IntArray &covariates);
void MakeTraits(Matrix &alltraits, StringArray &alltraitNames, IntArray &covariates, IntArray &affections, IntArray &quantraits);
void AdjustCovariate(Matrix &risks, Matrix &alltraits, StringArray &alltraitNames, IntArray &quantraits, IntArray &covariates, Matrix &allY, unsigned long long int **masks, int coreCount);
void AdjustNoCovariate(Matrix &risks, Matrix &alltraits, StringArray &alltraitNames, IntArray &quantraits, Matrix &allY, unsigned long long int **masks, int coreCount);
void SECRET_Read(const char *phefile, const char *dictfile, const char *rsfile, const char *pcarsfile, const char *prefix, const char *covariateList,
    StringArray &sampleID, StringArray &alltraitNames, Matrix &alltraits, IntArray &covariates, IntArray &affections,
    IntArray &quantraits, StringArray &rsNames, Matrix &risks, StringArray &annotation, int coreCount);
