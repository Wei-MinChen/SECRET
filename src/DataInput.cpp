////////////////////////////////////////////////////////////// ////////
// DataInput.cpp
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
// Oct 6, 2026

#include <time.h>
#include <math.h>
#include <string.h>
#include "MathSVD.h"
#include "StringArray.h"
#include "MathStats.h"
#include "QuickIndex.h"
#include "StringHash.h"
#include "DataInput.h"
#ifdef _OPENMP
  #include <omp.h>
#endif

//extern void AdjustPC1(Matrix &risks, int coreCount);
const char *gTraitList = "";     // set from --trait in Main.cpp
const char *gCovarFile = "";     // set from --covarfile in Main.cpp
bool gHasCoverage = false;       // set by ReadKeys when the dictionary has SNPused / Coverage columns
Vector gSNPused, gCoverage;      // per PGS (same order as rsNames), _NAN_ when unknown
bool gRecodeBinary01 = true;   // survival scans turn this off: their event column is 0/1 by design
extern void localPCA(Matrix &risks, int dimPC, Vector &Values, Matrix &PCs, int coreCount);

void ReadLines(const char *filename, StringArray &lines, bool printFlag=true)
{
  FILE *fp = fopen(filename, "rb");
  if (fp == NULL) error("File %s cannot open to read.", filename);
  if(printFlag) {printf("Reading %s...", filename); fflush(stdout);}
  lines.Dimension(0);
  const size_t CHUNK = 0x100000;  // 1M bytes; lines of any length are allowed
  char *buffer = new char[CHUNK + 1];
  String current;   // the line being assembled; may span several chunks
  size_t bytes;
  while ((bytes = fread(buffer, 1, CHUNK, fp)) > 0) {
    buffer[bytes] = '\0';
    size_t pos = 0;
    while (pos < bytes) {
      char *nl = (char *) memchr(buffer + pos, '\n', bytes - pos);
      if (nl == NULL) {   // no line end in this chunk, carry the rest over
        current.Add(buffer + pos);
        break;
      }
      *nl = '\0';
      current.Add(buffer + pos);
      if (current.Length() > 0) lines.Push(current);  // blank lines are skipped
      current.Clear();
      pos = (nl - buffer) + 1;
    }
  } // End of while
  if (current.Length() > 0) lines.Push(current);      // last line without a newline
  delete [] buffer;
  fclose(fp);
  if(printFlag){printf("done.\n"); fflush(stdout);}
}

void LineToMatrix(StringArray &lines, StringArray &colNames, StringArray &rowNames, Matrix &M, char &sep, int coreCount)
{
  int lineCount = lines.Length();
  if (lines[0].Find('\r') > -1 || lines[0].Find('"') > -1)
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for (int i = 0; i < lineCount; i++) {
    lines[i].ExcludeCharacters('\r');
    lines[i].ExcludeCharacters('"');
  }   // End of i loop
  sep = '\t';
  String line(lines[0]);
  StringArray tokens(0);
  tokens.AddTokens(line, sep);
  int colCount = tokens.Length() - 1;
  if (colCount < 1) {
   sep = ' ';
   tokens.Clear();
   tokens.AddTokens(line, sep);
   colCount = tokens.Length() - 1;
   if (colCount < 1) {
     sep = ',';
     tokens.Clear();
     tokens.AddTokens(line, sep);
     colCount = tokens.Length() - 1;
     if (colCount < 1) error("Please use tab, space, or commas-delimited format for input file.");
   }
  }
  colNames.Dimension(colCount);
  for (int r = 0; r < colCount; r++)
    colNames[r] = tokens[r + 1];
  int idCount = lineCount - 1;
  rowNames.Dimension(idCount);
  if (idCount == 0) return;
  M.Dimension(idCount, colCount);
  M.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(tokens)
#endif
  for(int i = 0; i < idCount; i++) {
    tokens.Clear();
    tokens.AddTokens(lines[i+1], sep);
    if(tokens.Length() != colCount + 1) continue;
    rowNames[i] = tokens[0];
    for(int r = 0; r < colCount; r++)
      M[i][r] = tokens[r + 1].IsNumber() ? tokens[r + 1].AsDouble() : _NAN_;
  } // End of i loop
  IntArray invalid(0);
  for(int i = 0; i < idCount; i++)
    if (rowNames[i].IsEmpty()) invalid.Push(i);
  int invalidCount = invalid.Length();
  if (invalidCount)
    printf("WARNING: %d line(s) were skipped because the number of fields differs from the header. Empty cells must be written as NA.\n", invalidCount);
  for(int i = invalidCount - 1; i >= 0; i--) {
    rowNames.Delete(invalid[i]);
    M.DeleteRow(invalid[i]);
  } // End of i loop
}

void LineToMatrix(StringArray &lines, StringArray &colNames, Matrix &M, char &sep, int coreCount)
{
  int lineCount = lines.Length();
  if (lines[0].Find('\r') > -1 || lines[0].Find('"') > -1)
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for (int i = 0; i < lineCount; i++) {
    lines[i].ExcludeCharacters('\r');
    lines[i].ExcludeCharacters('"');
  }   // End of i loop
  sep = '\t';
  String line(lines[0]);
  StringArray tokens(0);
  tokens.AddTokens(line, sep);
  int colCount = tokens.Length();
  if (colCount < 2) {
   sep = ' ';
   tokens.Clear();
   tokens.AddTokens(line, sep);
   colCount = tokens.Length();
   if (colCount < 2) {
     sep = ',';
     tokens.Clear();
     tokens.AddTokens(line, sep);
     colCount = tokens.Length();
     if (colCount < 2) error("Please use tab, space, or commas-delimited format for input file.");
   }
  }
  colNames.Dimension(colCount);
  for (int r = 0; r < colCount; r++)
    colNames[r] = tokens[r];
  int idCount = lineCount-1;
  if (idCount == 0) return;
  M.Dimension(idCount, colCount);
  M.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(tokens)
#endif
  for(int i = 0; i < idCount; i++) {
    tokens.Clear();
    tokens.AddTokens(lines[i+1], sep);
    if(tokens.Length() != colCount) continue;
    for(int r = 0; r < colCount; r++)
      M[i][r] = tokens[r].IsNumber() ? tokens[r].AsDouble() : _NAN_;
  } // End of i loop
}

/*
  IntArray invalid(0);
  for(int i = 0; i < idCount; i++)
    if (rowNames[i].IsEmpty()) invalid.Push(i);
  int invalidCount = invalid.Length();
  for(int i = invalidCount - 1; i >= 0; i--) {
    rowNames.Delete(invalid[i]);
    M.DeleteRow(invalid[i]);
  } // End of i loop
  */

void LineToMatrix2(StringArray &lines, StringArray &colNames, StringArray &rowNames1, StringArray &rowNames2, Matrix &M, char &sep, int coreCount)
{
    int lineCount = lines.Length();
    if (lines[0].Find('\r') > -1 || lines[0].Find('"') > -1)
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int i = 0; i < lineCount; i++) {
            lines[i].ExcludeCharacters('\r');
            lines[i].ExcludeCharacters('"');
        }   // End of i loop
    sep = '\t';
    String line(lines[0]);
    StringArray tokens(0);
    tokens.AddTokens(line, sep);
    int colCount = tokens.Length() - 2;
    if (colCount < 1) {
        sep = ' ';
        tokens.Clear();
        tokens.AddTokens(line, sep);
        colCount = tokens.Length() - 2;
        if (colCount < 1) {
            sep = ',';
            tokens.Clear();
            tokens.AddTokens(line, sep);
            colCount = tokens.Length() - 2;
            if (colCount < 1) error("Please use tab, space, or commas-delimited format for input file.");
        }
    }
    colNames.Dimension(colCount);
    for (int r = 0; r < colCount; r++)
        colNames[r] = tokens[r + 2];
    int idCount = lineCount - 1;
    rowNames1.Dimension(idCount);
    rowNames2.Dimension(idCount);
    if (idCount == 0) return;
    M.Dimension(idCount, colCount);
    M.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(tokens)
#endif
    for (int i = 0; i < idCount; i++) {
        tokens.Clear();
        tokens.AddTokens(lines[i + 1], sep);
        if (tokens.Length() != colCount + 2) continue;
        rowNames1[i] = tokens[0];
        rowNames2[i] = tokens[1];
        for (int r = 0; r < colCount; r++)
            M[i][r] = tokens[r + 2].IsNumber() ? tokens[r + 2].AsDouble() : _NAN_;
    }   // End of while
    IntArray invalid(0);
    for (int i = 0; i < idCount; i++)
        if (rowNames1[i].IsEmpty()) invalid.Push(i);
    int invalidCount = invalid.Length();
    if (invalidCount)
        printf("WARNING: %d line(s) were skipped because the number of fields differs from the header. Empty cells must be written as NA.\n", invalidCount);
    for (int i = invalidCount - 1; i >= 0; i--) {
        rowNames1.Delete(invalid[i]);
        rowNames2.Delete(invalid[i]);
        M.DeleteRow(invalid[i]);
    }   // End of i loop
}

void LineToMatrix4(StringArray &lines, StringArray &colNames, StringArray &rowNames1, StringArray &rowNames2, StringArray &rowNames3, StringArray &rowNames4, Matrix &M, char &sep, int coreCount)
{
    int lineCount = lines.Length();
    if (lines[0].Find('\r') > -1 || lines[0].Find('"') > -1)
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int i = 0; i < lineCount; i++) {
            lines[i].ExcludeCharacters('\r');
            lines[i].ExcludeCharacters('"');
        }   // End of i loop
    sep = '\t';
    String line(lines[0]);
    StringArray tokens(0);
    tokens.AddTokens(line, sep);
    int colCount = tokens.Length() - 4;
    if (colCount < 1) {
        sep = ' ';
        tokens.Clear();
        tokens.AddTokens(line, sep);
        colCount = tokens.Length() - 4;
        if (colCount < 1) {
            sep = ',';
            tokens.Clear();
            tokens.AddTokens(line, sep);
            colCount = tokens.Length() - 4;
            if (colCount < 1) error("Please use tab, space, or commas-delimited format for input file.");
        }
    }
    colNames.Dimension(colCount);
    for (int r = 0; r < colCount; r++)
        colNames[r] = tokens[r + 4];
    int idCount = lineCount - 1;
    rowNames1.Dimension(idCount);
    rowNames2.Dimension(idCount);
    rowNames3.Dimension(idCount);
    rowNames4.Dimension(idCount);
    if (idCount == 0) return;
    M.Dimension(idCount, colCount);
    M.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(tokens)
#endif
    for (int i = 0; i < idCount; i++) {
        tokens.Clear();
        tokens.AddTokens(lines[i + 1], sep);
        if (tokens.Length() != colCount + 4) continue;
        rowNames1[i] = tokens[0];
        rowNames2[i] = tokens[1];
        rowNames3[i] = tokens[2];
        rowNames4[i] = tokens[3];
        for (int r = 0; r < colCount; r++)
            M[i][r] = tokens[r + 4].IsNumber() ? tokens[r + 4].AsDouble() : _NAN_;
    }   // End of while
    IntArray invalid(0);
    for (int i = 0; i < idCount; i++)
        if (rowNames1[i].IsEmpty()) invalid.Push(i);
    int invalidCount = invalid.Length();
    if (invalidCount)
        printf("WARNING: %d line(s) were skipped because the number of fields differs from the header. Empty cells must be written as NA.\n", invalidCount);
    for (int i = invalidCount - 1; i >= 0; i--) {
        rowNames1.Delete(invalid[i]);
        rowNames2.Delete(invalid[i]);
        rowNames3.Delete(invalid[i]);
        rowNames4.Delete(invalid[i]);
        M.DeleteRow(invalid[i]);
    }   // End of i loop
}

void ReadMatrix(const char *phefile, StringArray &colNames, StringArray &rowNames, Matrix &M, int coreCount)
{
  StringArray lines;
  ReadLines(phefile, lines);
  char sep;
  LineToMatrix(lines, colNames, rowNames, M, sep, coreCount);
}

void ReadMatrix(const char *phefile, StringArray &colNames, Matrix &M, int coreCount)
{
  StringArray lines;
  ReadLines(phefile, lines, false);
  char sep;
  LineToMatrix(lines, colNames, M, sep, coreCount);
}

// ---------------------------------------------------------------------
// Several input files for one option (--pheno, --covarfile, --score)
// ---------------------------------------------------------------------
// The option value is a comma-separated list of file names.  The first file
// defines the samples.  Every further file only contributes its columns
// (traits, covariates or PGS), matched to the samples of the first file by
// ID: the union is taken over the column names, never over the samples.  A
// sample absent from a further file gets missing values; a column name that
// is already present is skipped, so the earlier file wins.

void SplitFileList(const char *list, StringArray &files)
{
  files.Clear();
  StringArray tokens;
  tokens.AddTokens(String(list), ',');
  for (int t = 0; t < tokens.Length(); t++) {
    String name(tokens[t]);
    name.Trim();
    if (!name.IsEmpty()) files.Push(name);
  }   // End of t loop
}

static void MergeColumnFile(const char *file, const char *unit, StringArray &alltraitNames,
  StringArray &familyID, StringArray &sampleID, Matrix &alltraits, int coreCount);

void ReadTraits(const char *phefile, StringArray &alltraitNames,
  StringArray &familyID, StringArray &sampleID, Matrix &alltraits, int coreCount)
{
  StringArray files;
  SplitFileList(phefile, files);
  if (files.Length() == 0) error("Phenotype file is not specified.");
  StringArray lines;
  ReadLines(files[0], lines);
  char sep;
  if (lines[0].Find("FID") == 0)  //FID IID
    LineToMatrix2(lines, alltraitNames, familyID, sampleID, alltraits, sep, coreCount);
  else
    LineToMatrix(lines, alltraitNames, sampleID, alltraits, sep, coreCount);
  if (alltraitNames[0] == "FA" && alltraitNames[1] == "MO") {
    int idCount = sampleID.Length();
    alltraitNames.Delete(1);
    alltraitNames.Delete(0);
    int alltraitCount = alltraitNames.Length();
    for (int i = 0; i < idCount; i++)
      for (int t = 0; t < alltraitCount; t++)
        alltraits[i][t] = alltraits[i][t+2];
    alltraits.Dimension(idCount, alltraitCount);
  }
  for (int f = 1; f < files.Length(); f++)   // further phenotype files add traits, not samples
    MergeColumnFile(files[f], "trait", alltraitNames, familyID, sampleID, alltraits, coreCount);
}

// --score: one or more files.  Rows are samples, columns are PGS.
void ReadScoreMatrix(const char *rsfile, StringArray &rsNames, StringArray &rowNames, Matrix &M, int coreCount)
{
  StringArray files;
  SplitFileList(rsfile, files);
  if (files.Length() == 0) error("Polygenic score file is not specified.");
  ReadMatrix(files[0], rsNames, rowNames, M, coreCount);
  int rowCount = rowNames.Length();
  for (int f = 1; f < files.Length(); f++) {
    StringArray names2, rowNames2;
    Matrix M2;
    ReadMatrix(files[f], names2, rowNames2, M2, coreCount);
    int rowCount2 = rowNames2.Length();
    if (rowCount2 == 0 || rowCount == 0)
      error("Score file %s has no valid data lines.", (const char*)files[f]);
    StringIntHash byRow2, haveName;
    byRow2.Clear();
    haveName.Clear();
    for (int k = 0; k < rowCount2; k++)
      if (byRow2.Integer(rowNames2[k]) == -1) byRow2.SetInteger((const char*)rowNames2[k], k);
    IntArray row(rowCount);   // row of this file for each sample of the first file (-1: not found)
    int matched = 0;
    for (int i = 0; i < rowCount; i++) {
      row[i] = byRow2.Integer(rowNames[i]);
      if (row[i] > -1) matched++;
    }   // End of i loop
    if (matched == 0) {
      printf("None of the samples in %s have matching IDs in the first score file %s.\n", (const char*)files[f], (const char*)files[0]);
      printf("IDs in %s look like:\n", (const char*)files[f]);
      for (int i = 0; i < 5 && i < rowCount2; i++)
        printf("  %s\n", (const char*)rowNames2[i]);
      printf("IDs in %s look like:\n", (const char*)files[0]);
      for (int i = 0; i < 5 && i < rowCount; i++)
        printf("  %s\n", (const char*)rowNames[i]);
      printf("Please modify IDs in one of the score files and rerun again.\n");
      exit(0);
    }
    int oldCount = rsNames.Length();
    for (int c = 0; c < oldCount; c++)
      haveName.SetInteger((const char*)rsNames[c], c);
    IntArray add(0);
    int skipped = 0;
    for (int c = 0; c < names2.Length(); c++)
      if (haveName.Integer(names2[c]) > -1) skipped++;
      else {
        haveName.SetInteger((const char*)names2[c], oldCount + add.Length());
        add.Push(c);
      }
    if (skipped)
      printf("WARNING: %d PGS in %s are skipped because the names are already in an earlier score file.\n", skipped, (const char*)files[f]);
    int addCount = add.Length();
    if (addCount > 0) {
      Matrix merged(rowCount, oldCount + addCount);
      merged.Set(_NAN_);
      for (int i = 0; i < rowCount; i++) {
        for (int c = 0; c < oldCount; c++)
          merged[i][c] = M[i][c];
        if (row[i] > -1)
          for (int a = 0; a < addCount; a++)
            merged[i][oldCount + a] = M2[row[i]][add[a]];
      }   // End of i loop
      M = merged;
      for (int a = 0; a < addCount; a++)
        rsNames.Push(names2[add[a]]);
    }
    printf("%d PGS for %d of %d samples loaded from %s.\n", addCount, matched, rowCount, (const char*)files[f]);
  }   // End of f loop
}

void ReadRisks(const char *rsfile, StringArray &rsNames, StringArray &familyID, StringArray &sampleID, Matrix &risks, IntArray &rsKey, int coreCount)
{
    int idCount = sampleID.Length();
    StringArray rowNames;
    Matrix M;
    ReadScoreMatrix(rsfile, rsNames, rowNames, M, coreCount);
    int rsCount = rsNames.Length();
    if (rowNames[0].Find("_") > -1 && familyID.Length()==idCount) {
        if (sampleID.Find(rowNames[0]) > -1)
            error("Please remove the _ in IID, or use the sampleID format");
        StringArray IID = sampleID;
        sampleID = familyID;
        for (int i = 0; i < idCount; i++) {
            sampleID[i].Add('_');
            sampleID[i].Add(IID[i]);
        }   // End of i loop
    }
    StringIntHash sampleLookup;
    sampleLookup.Clear();
    for (int i = 0; i < idCount; i++)
        sampleLookup.SetInteger((const char*)sampleID[i], i);
    IntArray index(idCount); index.Set(-1);
    int keyCount = 0;
    int rowCount = rowNames.Length();
    for (int k = 0; k < rowCount; k++) {
        int i = sampleLookup.Integer(rowNames[k]);
        if (i > -1) {
            index[i] = k;
            keyCount++;
        }
    }   // End of k loop
    if (keyCount) {
        risks.Dimension(rsCount, idCount);
        risks.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int r = 0; r < rsCount; r++)
            for (int i = 0; i < idCount; i++)
                if (index[i] > -1)
                    risks[r][i] = M[index[i]][r];
        printf("%d samples have both phenotypes and risk scores.\n", keyCount);
        rsKey.Dimension(keyCount);
        if (idCount > keyCount) {   // redo risks
            int count = 0;
            for (int i = 0; i < idCount; i++) {
                if (index[i] == -1) continue;
                rsKey[count] = i;
                for (int r = 0; r < rsCount; r++)
                    risks[r][count] = risks[r][i];
                count++;
            }   // End of i loop
            risks.Dimension(rsCount, keyCount);
        }
        else {
            for (int i = 0; i < idCount; i++) rsKey[i] = i;
        }
    }
    else {
        printf("None of the risk scores have matching IDs with the phenotype data.\n");
        printf("IDs in the risk score file %s look like:\n", rsfile);
        for (int i = 0; i < 5 && i < rowCount; i++)
            printf("  %s\n", (const char*)rowNames[i]);
        printf("IDs in the phenotype data look like:\n");
        for (int i = 0; i < 5 && i < idCount; i++)
            printf("  %s\n", (const char*)sampleID[i]);
        printf("Please modify IDs in one of the input files and rerun again.\n");
        exit(0);
    }
    IntArray validrsIndex(0);
    for (int r = 0; r < rsCount; r++) {
        // mean and variance over non-missing values only; missing values stay _NAN_
        double sum = 0.0, ss = 0.0;
        int n = 0;
        for (int i = 0; i < keyCount; i++)
            if (risks[r][i] != _NAN_) { sum += risks[r][i]; n++; }
        double mean = n ? sum / n : 0.0;
        for (int i = 0; i < keyCount; i++)
            if (risks[r][i] != _NAN_) ss += (risks[r][i] - mean) * (risks[r][i] - mean);
        double sd0 = n > 1 ? sqrt(ss / (n - 1)) : 0.0;
        if (sd0 > 1E-10){
            for (int i = 0; i < keyCount; i++)
                if (risks[r][i] != _NAN_) risks[r][i] /= sd0;
            validrsIndex.Push(r);
        }
    }   // End of r loop
    int validrsCount = validrsIndex.Length();
    if (validrsCount < rsCount) {
        printf("%d risk scores removed for lacking variation.\n", rsCount-validrsCount);
        for (int r = 0; r < validrsCount; r++) {
            rsNames[r] = rsNames[validrsIndex[r]];
            for (int i = 0; i < keyCount; i++)
                risks[r][i] = risks[validrsIndex[r]][i];
        }   // End of r loop
        rsNames.Dimension(validrsCount);
        risks.Dimension(validrsCount, keyCount);
    }
}

// Split one comma-delimited line into fields.  A field may be wrapped in double
// quotes to contain commas, and "" inside quotes stands for a literal quote.
static void SplitCSVLine(const String &line, StringArray &fields)
{
    fields.Clear();
    String field;
    bool quoted = false;
    int len = line.Length();
    for (int i = 0; i < len; i++) {
        char c = line[i];
        if (quoted) {
            if (c == '\"') {
                if (i + 1 < len && line[i + 1] == '\"') { field.Add('\"'); i++; }
                else quoted = false;
            }
            else field.Add(c);
        }
        else if (c == '\"') quoted = true;
        else if (c == ',') { fields.Push(field); field.Clear(); }
        else if (c != '\r') field.Add(c);
    }   // End of i loop
    fields.Push(field);
}

// A column header reduced to lower-case letters and digits, so that "SNP_used",
// "SNPs used" and "snpused" are the same name.
static String NormalizeHeader(const String &s)
{
    String out;
    int len = s.Length();
    for (int i = 0; i < len; i++) {
        char c = s[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        if ((c >= 'a' && c <= 'z') || (c >= '0' && c <= '9')) out.Add(c);
    }   // End of i loop
    return out;
}

// First column whose normalized name is in the alias list (the list is in order of
// preference), or -1.
static int FindHeaderColumn(StringArray &normalized, const char *const *aliases)
{
    for (int a = 0; aliases[a] != NULL; a++)
        for (int c = 0; c < normalized.Length(); c++)
            if (normalized[c] == aliases[a]) return c;
    return -1;
}

// A number such as 1234, 0.93 or 93.1%; _NAN_ for anything else (NA, empty, text)
static double ParseCoverageNumber(const String &field)
{
    String t(field);
    t.Trim();
    t.ExcludeCharacters("%,");
    if (t.IsEmpty() || !t.IsNumber()) return _NAN_;
    return t.AsDouble();
}

// Comma-delimited dictionary with a header line (e.g., the scores.coverage file
// from PGS Catalog scoring).  The annotation is the column named "trait"; the
// score name is the column named "score", or the first column if there is none.
//
// When the file also tells how many variants of each score were found in the
// data and what fraction of the score that is, these are kept in gSNPused and
// gCoverage and appear as the SNPused and Coverage columns of the scan tables.
// The columns are recognised by name (case, spaces, "_" and "-" are ignored):
//   SNPused   snpused, snpsused, variantsused, nvariantsused, used, matched,
//             nmatched, found, variantsfound, overlap, included, ...
//   Coverage  coverage, coveragepct, percentcoverage, matchrate, percentmatched, ...
//   (total)   variants, nvariants, snps, totalvariants, ...  When a file has the
//             SNPused and the total but no coverage, Coverage = SNPused / total.
// Coverage is reported as given in the file (a fraction or a percentage).
static void ReadKeysInCSV(const char *filename, StringArray &lines, StringArray &rsNames, StringArray &annotation)
{
    int rsCount = rsNames.Length();
    annotation.Dimension(rsCount);
    for (int r = 0; r < rsCount; r++)
        annotation[r] = rsNames[r];
    StringIntHash rsLookup;
    rsLookup.Clear();
    for (int r = 0; r < rsCount; r++)
        rsLookup.SetInteger((const char*)rsNames[r], r);
    StringArray fields;
    SplitCSVLine(lines[0], fields);
    int idCol = 0, traitCol = -1;
    for (int c = 0; c < fields.Length(); c++) {
        String name(fields[c]);
        name.ToLower();
        if (name == "trait") traitCol = c;
        else if (name == "score" && idCol == 0) idCol = c;
    }   // End of c loop
    if (traitCol == -1) {
        String found;
        for (int c = 0; c < fields.Length(); c++) {
            if (c) found.Add(", ");
            found.Add(fields[c]);
        }   // End of c loop
        error("Dictionary file %s is comma-delimited but has no column named \"trait\".\nColumns found: %s",
            filename, (const char*)found);
    }
    // SNPused / Coverage columns, when the file has them
    static const char *const usedAliases[] = {"variantsused", "snpused", "snpsused", "nsnpused", "nsnpsused",
        "variantused", "nvariantsused", "usedvariants", "usedsnps", "nused", "used", "snpsmatched", "nsnpsmatched",
        "variantsmatched", "matchedvariants", "matchedsnps", "nmatched", "matched", "variantsfound", "snpsfound",
        "nfound", "found", "variantsincluded", "snpsincluded", "included", "noverlap", "overlap", NULL};
    static const char *const coverageAliases[] = {"coverage", "coveragepct", "coveragepercent", "percentcoverage",
        "pctcoverage", "coveragerate", "fractioncovered", "covered", "matchrate", "percentmatched", "pctmatched",
        "matchpercent", "matchedpercent", NULL};
    static const char *const totalAliases[] = {"variants", "nvariants", "totalvariants", "variantstotal",
        "numvariants", "snps", "nsnps", "totalsnps", "variantsinscore", "scorevariants", "variantsinscorefile",
        "scoresize", NULL};
    StringArray headerFields = fields;   // fields is reused for the data lines below
    StringArray normalized(fields.Length());
    for (int c = 0; c < fields.Length(); c++)
        normalized[c] = NormalizeHeader(fields[c]);
    int usedCol = FindHeaderColumn(normalized, usedAliases);
    int coverageCol = FindHeaderColumn(normalized, coverageAliases);
    int totalCol = FindHeaderColumn(normalized, totalAliases);
    if (totalCol == usedCol) totalCol = -1;
    bool deriveCoverage = (coverageCol == -1 && usedCol > -1 && totalCol > -1);
    gHasCoverage = (usedCol > -1 || coverageCol > -1 || deriveCoverage);
    if (gHasCoverage) {
        gSNPused.Dimension(rsCount);
        gCoverage.Dimension(rsCount);
        gSNPused.Set(_NAN_);
        gCoverage.Set(_NAN_);
    }
    int needed = (idCol > traitCol ? idCol : traitCol) + 1;
    int lineCount = lines.Length();
    int matched = 0;
    for (int i = 1; i < lineCount; i++) {
        SplitCSVLine(lines[i], fields);
        if (fields.Length() < needed) continue;
        int newIndex = rsLookup.Integer(fields[idCol]);
        if (newIndex > -1 && !fields[traitCol].IsEmpty()) {
            annotation[newIndex] = fields[traitCol];
            matched++;
        }
        if (gHasCoverage && newIndex > -1) {
            if (usedCol > -1 && usedCol < fields.Length())
                gSNPused[newIndex] = ParseCoverageNumber(fields[usedCol]);
            if (coverageCol > -1 && coverageCol < fields.Length())
                gCoverage[newIndex] = ParseCoverageNumber(fields[coverageCol]);
            else if (deriveCoverage && totalCol < fields.Length()) {
                double total = ParseCoverageNumber(fields[totalCol]);
                if (gSNPused[newIndex] != _NAN_ && total != _NAN_ && total > 0.0)
                    gCoverage[newIndex] = gSNPused[newIndex] / total;
            }
        }
    }   // End of i loop
    printf("  %d of %d scores annotated by the \"trait\" column of %s\n", matched, rsCount, filename);
    if (gHasCoverage) {
        printf("  Coverage information found in %s:", filename);
        if (usedCol > -1) printf(" SNPused = column \"%s\"", (const char*)headerFields[usedCol]);
        if (coverageCol > -1) printf("%s Coverage = column \"%s\"", usedCol > -1 ? ";" : "", (const char*)headerFields[coverageCol]);
        else if (deriveCoverage) printf("; Coverage = \"%s\" / \"%s\"", (const char*)headerFields[usedCol], (const char*)headerFields[totalCol]);
        printf("\n");
    }
}

static void ReadKeysInText(StringArray &lines, StringArray &rsNames, StringArray &annotation)
{
    int rsCount = rsNames.Length();
    annotation.Dimension(rsCount);
    for (int r = 0; r < rsCount; r++)
        annotation[r] = rsNames[r];
    StringIntHash rsLookup;
    rsLookup.Clear();
    for (int r = 0; r < rsCount; r++)
        rsLookup.SetInteger((const char*)rsNames[r], r);
    int lineCount = lines.Length();
    StringArray tokens;
    for (int i = 0; i < lineCount; i++) {
        tokens.Clear();
        tokens.AddTokens(lines[i], '\t');
        if (tokens.Length() < 2) continue;
        int newIndex = rsLookup.Integer(tokens[0]);
        if (newIndex > -1) annotation[newIndex] = tokens[1];
    }   // End of i loop
}

void ReadKeysInHTML(const char *filename, StringArray &rsNames, StringArray &annotation)
{
    int rsCount = rsNames.Length();
    annotation.Dimension(rsCount);
    for (int r = 0; r < rsCount; r++)
        annotation[r] = rsNames[r];
    StringIntHash rsLookup;
    rsLookup.Clear();
    for (int r = 0; r < rsCount; r++)
        rsLookup.SetInteger((const char*)rsNames[r], r);
    StringArray lines;
    ReadLines(filename, lines);
    int lineCount = lines.Length();
    String line, token;
    for (int i = 0; i < lineCount; i++) {
        int start = lines[i].Find("data-meta");
        if (start == -1) continue;
        line = lines[i].SubStr(start + 10);
        line.ExcludeCharacters('"');
        int commaPos = line.Find(",");
        if (commaPos == -1) continue;
        token = line.SubStr(0, commaPos);
        int newIndex = rsLookup.Integer(token);
        token = line.SubStr(commaPos + 1);
        if (newIndex > -1) annotation[newIndex] = token;
    }   // End of i loop
}

// SNPused and Coverage columns of the scan tables (.linear, .logistic, .coxph): they are
// written, before the Annotation column, only when the dictionary file has coverage information.
// The same two pieces of text as WriteCoverageHeader / WriteCoverageFields, put in a
// character buffer (at least 64 characters) for the writers that build their lines in memory.
// The return value is the number of characters written (0 without coverage information).
int FormatCoverageHeader(char *buffer)
{
    buffer[0] = '\0';
    if (!gHasCoverage) return 0;
    return sprintf(buffer, "\tSNPused\tCoverage");
}

int FormatCoverageFields(char *buffer, int rs)
{
    buffer[0] = '\0';
    if (!gHasCoverage) return 0;
    double used = (rs < gSNPused.Length()) ? gSNPused[rs] : _NAN_;
    double coverage = (rs < gCoverage.Length()) ? gCoverage[rs] : _NAN_;
    int len = 0;
    if (used == _NAN_) len += sprintf(buffer + len, "\tNA");
    else if (used == floor(used) && fabs(used) < 1E15) len += sprintf(buffer + len, "\t%.0lf", used);
    else len += sprintf(buffer + len, "\t%.6g", used);
    if (coverage == _NAN_) len += sprintf(buffer + len, "\tNA");
    else len += sprintf(buffer + len, "\t%.5g", coverage);
    return len;
}

void WriteCoverageHeader(FILE *fp)
{
    char buffer[64];
    if (FormatCoverageHeader(buffer)) fputs(buffer, fp);
}

void WriteCoverageFields(FILE *fp, int rs)
{
    char buffer[128];
    if (FormatCoverageFields(buffer, rs)) fputs(buffer, fp);
}

void ReadKeys(const char *filename, StringArray &rsNames, StringArray &annotation)
{
    gHasCoverage = false;       // set again by ReadKeysInCSV when the file has coverage information
    gSNPused.Dimension(0);
    gCoverage.Dimension(0);
    String string(filename);
    if (string.Find(".html") > -1)
        ReadKeysInHTML(filename, rsNames, annotation);
    else {
        StringArray lines;
        ReadLines(filename, lines);
        if (lines.Length() > 0) {
            // drop a UTF-8 byte-order mark left by spreadsheet exports
            if (lines[0].Length() >= 3 && (unsigned char) lines[0][0] == 0xEF
                && (unsigned char) lines[0][1] == 0xBB && (unsigned char) lines[0][2] == 0xBF)
                lines[0] = lines[0].SubStr(3);
            // comma-delimited with a header: the tab-delimited format never looks like this
            if (lines[0].FindChar('\t') == -1 && lines[0].FindChar(',') > -1) {
                ReadKeysInCSV(filename, lines, rsNames, annotation);
                return;
            }
        }
        ReadKeysInText(lines, rsNames, annotation);
    }
}

void ReadTraits4(const char *phefile, StringArray &alltraitNames, StringArray &familyID, StringArray &sampleID, StringArray &FAID, StringArray &MOID, Matrix &alltraits, int coreCount)
{
  StringArray files;
  SplitFileList(phefile, files);
  if (files.Length() == 0) error("Phenotype file is not specified.");
  StringArray lines;
  ReadLines(files[0], lines);
  char sep;
  if (lines[0].Find("FID") == 0)  //FID IID
    LineToMatrix4(lines, alltraitNames, familyID, sampleID, FAID, MOID, alltraits, sep, coreCount);
  else
    LineToMatrix(lines, alltraitNames, sampleID, alltraits, sep, coreCount);
  for (int f = 1; f < files.Length(); f++)   // further phenotype files add traits, not samples
    MergeColumnFile(files[f], "trait", alltraitNames, familyID, sampleID, alltraits, coreCount);
}   

void MakeCovariates(const char *covariateList, StringArray &alltraitNames, StringArray &covList, IntArray &covariates)
{
    covList.AddTokens(covariateList, ',');
    covariates.Dimension(0);
    for (int i = 0; i < covList.Length(); i++) {
        int k = alltraitNames.Find(covList[i]);
        if (k > -1) covariates.Push(k);
        else printf("WARNING: covariate %s was not found in the phenotype or covariate file and is ignored.\n", (const char*)covList[i]);
    }   // End of i loop
}

void MakeTraits(Matrix &alltraits, StringArray &alltraitNames, IntArray &covariates, IntArray &affections, IntArray &quantraits)
{
  int idCount = alltraits.rows;
  int alltraitCount = alltraits.cols;
  int covariateCount = covariates.Length();
  // --trait: restrict the analysis to the listed traits (all traits when not given)
  IntArray selected(alltraitCount);
  selected.Set(1);
  if (gTraitList[0]) {
    selected.Zero();
    StringArray wanted;
    wanted.AddTokens(gTraitList, ',');
    String missing;
    for (int w = 0; w < wanted.Length(); w++) {
      int k = alltraitNames.Find(wanted[w]);
      if (k == -1) {
        if (!missing.IsEmpty()) missing.Add(", ");
        missing.Add(wanted[w]);
      }
      else if (covariateCount && covariates.Find(k) != -1)
        printf("WARNING: %s is a covariate and is not analyzed as a trait.\n", (const char*)wanted[w]);
      else selected[k] = 1;
    }   // End of w loop
    if (!missing.IsEmpty()) {
      String available;
      for (int t = 0; t < alltraitCount && t < 30; t++) {
        if (t) available.Add(", ");
        available.Add(alltraitNames[t]);
      }
      if (alltraitCount > 30) available.Add(", ...");
      error("Trait(s) given by --trait not found in the phenotype file: %s\nAvailable columns: %s",
        (const char*)missing, (const char*)available);
    }
    if (selected.Sum() == 0)
      error("None of the traits given by --trait can be analyzed (covariates are excluded).");
    printf("%d trait(s) selected by --trait %s\n", (int)selected.Sum(), gTraitList);
  }
  for (int t = 0; t < alltraitCount; t++) {
    if (covariateCount && covariates.Find(t) != -1) continue;
    if (!selected[t]) continue;
    bool binaryFlag = true;
    bool naFlag = false;
    for (int i = 0; i < idCount; i++)
      if (alltraits[i][t] != 0 && alltraits[i][t] != 1 && alltraits[i][t] != 2 && alltraits[i][t] != _NAN_) {
        binaryFlag = false;  break;
      }
      else if (alltraits[i][t] == _NAN_) naFlag = true;
      if (binaryFlag) {
        // a trait coded 0/1 (no 2s, both 0 and 1 present) is recoded to 1=control, 2=case;
        // otherwise 0 would be read as "missing" and every 1 as a control
        bool has0 = false, has1 = false, has2 = false;
        for (int i = 0; i < idCount; i++) {
          if (alltraits[i][t] == 0) has0 = true;
          else if (alltraits[i][t] == 1) has1 = true;
          else if (alltraits[i][t] == 2) has2 = true;
        }
        if (gRecodeBinary01 && has0 && has1 && !has2) {
          for (int i = 0; i < idCount; i++)
            if (alltraits[i][t] != _NAN_) alltraits[i][t] += 1.0;
          printf("Trait %s is coded 0/1 and was recoded to 1/2 (1=control, 2=case).\n", (const char*)alltraitNames[t]);
        }
        affections.Push(t);
        if (naFlag)
          for (int i = 0; i < idCount; i++)
        if (alltraits[i][t] == _NAN_) alltraits[i][t] = 0;
      }
      else quantraits.Push(t);
  }   // End of t loop
  int affectionCount = affections.Length();
  if (affectionCount) {
    printf("%d dichotomous traits for %d samples loaded. Here is the distribution:\n", affectionCount, idCount);
    printf("%20s%10s%10s%10s\n", "Trait", "N2", "N1", "NMISS");
    int statusCount[2];
    for (int a = 0; a < affectionCount; a++) {
      int trait = affections[a];
      printf("%20s", (const char*)alltraitNames[trait]);
      statusCount[0] = statusCount[1] = 0;
      for (int i = 0; i < idCount; i++)
        if (alltraits[i][trait] == 2) statusCount[1]++;
        else if (alltraits[i][trait] == 1) statusCount[0]++;
      printf("%10d%10d%10d\n", statusCount[1], statusCount[0], idCount - statusCount[0] - statusCount[1]);
    }   // End of a loop
  }
}

void AdjustCovariate(Matrix &risks, Matrix &alltraits, StringArray &alltraitNames,
  IntArray &quantraits, IntArray &covariates, Matrix &allY, unsigned long long int **masks, int coreCount)
{
  int idCount = risks.cols;
  int covariateCount = covariates.Length();
  int quanCount = quantraits.Length();
  Matrix X0(idCount, covariateCount+1); X0.Zero();
  IntArray cValid(idCount); cValid.Zero();
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int i = 0; i < idCount; i++){
    bool validFlag = true;
    for(int c = 0; c < covariateCount; c++){
      X0[i][c+1] = alltraits[i][covariates[c]];
      if(X0[i][c+1]==_NAN_) validFlag = false;
    }
    if(validFlag){
      X0[i][0] = 1;
      cValid[i] = 1;
    }else
      for(int c = 0; c < covariateCount; c++)
        X0[i][c+1] = 0;
  } // End of i loop
  allY.Zero();
  int twordCount = (quanCount-1)/64+1;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int i = 0; i < idCount; i++){
    for(int w = 0; w < twordCount; w++)
      masks[i][w] = 0;
    if(!cValid[i]) continue;
    for (int t = 0; t < quanCount; t+=64) {
      int tMax = t < quanCount-64? t+64: quanCount;
      int tw = (t>>6);
      unsigned long long int word0 = 0x1;
      for (int tt = t; tt < tMax; tt++, word0<<=1) {
        int trait = quantraits[tt];
        if(alltraits[i][trait] != _NAN_ && risks[0][i] != _NAN_){
          allY[i][tt] = alltraits[i][trait];
          masks[i][tw] |= word0;
        }
      } // End of tt loop
    } // End of t loop
  } // End of OMP i loop
  Matrix betas(quanCount, covariateCount+1);
  Vector beta;
  SVD svd;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount) private(svd,beta)
#endif
  for(int t = 0; t < quanCount; t++){
    int tWord = (t>>6);
    int tShift = (t&0x3F);
    Matrix X(idCount, covariateCount+1); X.Zero();
    IntArray mask(idCount);
    for(int i = 0; i < idCount; i++){
      mask[i] = (masks[i][tWord]>>tShift)&0x1;
      if(mask[i]) X[i] = X0[i];
    } // End of i loop
    svd.Decompose(X);
    svd.Edit();
    Vector Y(idCount);
    for(int i = 0; i < idCount; i++)
      Y[i] = allY[i][t];
    svd.BackSubst(Y);
    beta = svd.x;
    for (int c = 0; c < covariateCount; c++)
      if (svd.w[c + 1] == 0.0) beta[c + 1] = 0.0;
    for (int i = 0; i < idCount; i++)
      if(mask[i]){
        double sum = 0.0;
        for (int c = 0; c < covariateCount + 1; c++)
          sum += beta[c] * X[i][c];
        allY[i][t] -= sum;
      }
    betas[t] = beta;
  } // End of OMP t loop
  printf("Distribution of %d trait residuals after regression (N=%d):\n", quanCount, idCount);
  printf("%16s%7s", "Trait", "Mu");
  for(int c = 0; c < covariateCount; c++)
    printf("%7s", (const char*)alltraitNames[covariates[c]]);
  printf("\n");
  for(int t = 0; t < quanCount; t++) {
    int trait = quantraits[t];
    printf("%16s", (const char*)alltraitNames[trait]);
    for(int c = 0; c < covariateCount + 1; c++)
      printf("%7.2lf", betas[t][c]);
      printf("\n");
  }   // End of t loop
}

void AdjustNoCovariate(Matrix &risks, Matrix &alltraits, StringArray &alltraitNames,
  IntArray &quantraits, Matrix &allY, unsigned long long int **masks, int coreCount)
{
  int idCount = risks.cols;
  int quanCount = quantraits.Length();
  allY.Zero();
  int twordCount = (quanCount-1)/64+1;
  for(int i = 0; i < idCount; i++)
    for(int w = 0; w < twordCount; w++)
      masks[i][w] = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
  for(int i = 0; i < idCount; i++){
    for (int t = 0; t < quanCount; t+=64) {
      int tMax = t < quanCount-64? t+64: quanCount;
      int tw = (t>>6);
      unsigned long long int word0 = 0x1;
      for (int tt = t; tt < tMax; tt++, word0<<=1) {
        int trait = quantraits[tt];
        if(alltraits[i][trait] != _NAN_ && risks[0][i] != _NAN_){
          allY[i][tt] = alltraits[i][trait];
          masks[i][tw] |= word0;
        }
      } // End of tt loop
    } // End of t loop
  } // End of OMP i loop
}
/*
    printf("%d quantitative traits for %d samples loaded. Here is the distribution:\n", quanCount, idCount);
    printf("%20s%10s%10s%10s\n", "Trait", "N", "Mean", "SD");
    for (int t = 0; t < quanCount; t++) {
        int trait = quantraits[t];
        printf("%20s", (const char*)alltraitNames[trait]);
        double sum = 0.0;
        double ss = 0.0;
        int N = Ys[t].Length();
        for (int i = 0; i < N; i++)
            if (Ys[t][i] != _NAN_) {
                sum += Ys[t][i];
                ss += Ys[t][i] * Ys[t][i];
            }   // End of i loop
        printf("%10d%10.3lf%10.3lf\n", N, sum / N, sqrt((ss - sum * sum / N) / N));
    }   // End of t loop
*/

void UseResiduals1(Matrix &rs, Vector &PC, int coreCount)
{
    int rsCount = rs.rows;
    int idCount = PC.Length();
    if (rs.cols < idCount) error("Column number %d of risk score matrix less than PC length %d", rs.cols, idCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
    for (int r = 0; r < rsCount; r++) {
        double lxx = 0.0;
        double lxy = 0.0;
        double meanX = 0.0;
        double meanG = 0.0;
        for (int i = 0; i < idCount; i++) {
            lxx += PC[i] * PC[i];
            lxy += PC[i] * rs[r][i];
            meanX += PC[i];
            meanG += rs[r][i];
        }   // End of i loop
        meanX /= idCount;
        meanG /= idCount;
        double beta = (lxy - meanX * meanG * idCount) / (lxx - meanX * meanX * idCount);
        double mu = meanG - meanX * beta;
        for (int i = 0; i < idCount; i++)
            rs[r][i] -= (mu + beta * PC[i]);
    }   // End of r loop
}

void AdjustPC1(Matrix &risks, int coreCount)
{
    int rsCount = risks.rows;
    int idCount = risks.cols;
    IntArray validIndex(0);
    for (int r = 0; r < rsCount; r++)
        if (risks[r].Var() > 1E-10) validIndex.Push(r);
    int validrsCount = validIndex.Length();
    printf("%d scores are used for PCA\n", validrsCount); fflush(stdout);
    Matrix vrisks(idCount, validrsCount);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
for (int i = 0; i < idCount; i++)
    for (int r = 0; r < validrsCount; r++)
        vrisks[i][r] = risks[validIndex[r]][i];
Vector Values;
Matrix PCs;
localPCA(vrisks, 1, Values, PCs, coreCount);
printf("Adjusting PC1..."); fflush(stdout);
UseResiduals1(risks, PCs[0], coreCount);
printf("done\n"); fflush(stdout);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
for (int r = 0; r < rsCount; r++) {
    double v = risks[r].Var();
    if (v > 1E-100) risks[r].Multiply(1 / sqrt(v));
}   // End of r loop
}

void AdjustPC(Matrix &risks, Matrix &alltraits, StringArray &alltraitNames, StringArray &covariateList, int coreCount)
{
    int idCount = risks.cols;
    char temp[256];
    SVD svd;
    Vector beta, PC(idCount);
    int kPC1 = alltraitNames.Find("PC1");
    if (kPC1 == -1) {   // PC1 requested but not available
        AdjustPC1(risks, coreCount);
        printf("PC1 inferred from PGS data and adjusted in PGS.\n");
    }
    else {  // PC1 data available
        int pcCount = 1;
        for (; pcCount < 31; pcCount++) {
            sprintf(temp, "PC%d", pcCount + 1);
            if (covariateList.Find(temp) == -1 || alltraitNames.Find(temp) == -1) break;
        }   // End of pcCount loop
        for (int k = 0; k < pcCount; k++) {
            sprintf(temp, "PC%d", k + 1);
            int trait = alltraitNames.Find(temp);
            bool naFlag = false;
            for (int i = 0; i < idCount; i++) {
                PC[i] = alltraits[i][trait];
                if (PC[i] == _NAN_)
                    naFlag = true;
            }   // End of i loop
            if (naFlag) {
                double mean = 0.0;
                int valididCount = 0;
                for (int i = 0; i < idCount; i++)
                    if (PC[i] != _NAN_) {
                        mean += PC[i]; valididCount++;
                    }
                mean /= valididCount;
                for (int i = 0; i < idCount; i++)
                    if (PC[i] == _NAN_) PC[i] = mean;
            }
            UseResiduals1(risks, PC, coreCount);
        }   // End of k loop
        if (pcCount == 1)
            printf("PC1 adjusted in PGS.\n");
        else
            printf("PC1-%d adjusted in PGS.\n", pcCount);
    }
}

void ReadModels(const char*filename, StringArray &rowNames, Matrix &weights, int coreCount)
{
  StringArray colNames;
  ReadMatrix(filename, colNames, rowNames, weights, coreCount);
/*
  int pcarsCol = colNames.Find("WTSUM");
  if (pcarsCol > -1) {
    int rsCount = rowNames.Length();
    weight.Dimension(rsCount);
    for (int r = 0; r < rsCount; r++)
      weight[r] = M[r][pcarsCol];
  }*/


}

// ---------------------------------------------------------------------
// --covarfile: covariates kept in a file separate from the traits
// ---------------------------------------------------------------------
// The file has the same layout as the phenotype file: a header line and one
// line per sample, with the ID in the first column (or FID and IID in the
// first two columns, when the header starts with FID).  Every other column is
// a covariate.  Tab, space and comma delimiters are all accepted.
// --covarfile takes one or more files separated by commas; see the note above
// ReadTraits on how several files are combined.

// Names of the covariate columns of all files in the list, read from the header
// lines only.  A name that occurs in several files is listed once.
String CovarFileColumns(const char *covarfile)
{
  StringArray files;
  SplitFileList(covarfile, files);
  String names;
  StringIntHash seen;
  seen.Clear();
  for (int f = 0; f < files.Length(); f++) {
    FILE *fp = fopen(files[f], "rb");
    if (fp == NULL) error("File %s cannot open to read.", (const char*)files[f]);
    String line;
    int c;
    while ((c = fgetc(fp)) != EOF && c != '\n')
      if (c != '\r' && c != '"') line.Add((char)c);
    fclose(fp);
    StringArray tokens;
    char seps[3] = {'\t', ' ', ','};
    int skip = (line.Find("FID") == 0) ? 2 : 1;
    for (int s = 0; s < 3; s++) {
      tokens.Clear();
      tokens.AddTokens(line, seps[s]);
      if (tokens.Length() - skip >= 1) break;
    }   // End of s loop
    if (tokens.Length() - skip < 1)
      error("Covariate file %s must have an ID column and at least one covariate column,\n  in tab, space, or commas-delimited format.", (const char*)files[f]);
    for (int t = skip; t < tokens.Length(); t++) {
      if (seen.Integer(tokens[t]) > -1) continue;
      seen.SetInteger((const char*)tokens[t], 1);
      if (names.Length() > 0) names.Add(',');
      names.Add(tokens[t]);
    }   // End of t loop
  }   // End of f loop
  return names;
}

// Append the columns of a further phenotype file or of a covariate file to the
// phenotype matrix, matching samples by ID.  A sample absent from the file gets
// missing values, and so is left out of any analysis that uses those columns.
// unit is "trait" or "covariate", for the messages.
static void MergeColumnFile(const char *covarfile, const char *unit, StringArray &alltraitNames,
  StringArray &familyID, StringArray &sampleID, Matrix &alltraits, int coreCount)
{
  int idCount = sampleID.Length();
  if (idCount == 0) return;
  StringArray lines;
  ReadLines(covarfile, lines);
  if (lines.Length() < 2)
    error("File %s has no data lines.", covarfile);
  bool covHasFID = (lines[0].Find("FID") == 0);
  bool pheHasFID = (familyID.Length() == idCount);
  StringArray covNames, covFID, covIID;
  Matrix C;
  char sep;
  if (covHasFID)
    LineToMatrix2(lines, covNames, covFID, covIID, C, sep, coreCount);
  else
    LineToMatrix(lines, covNames, covIID, C, sep, coreCount);
  int covRows = covIID.Length();
  int covCount = covNames.Length();
  if (covRows == 0)
    error("File %s has no valid data lines.", covarfile);
  StringIntHash byIID, byBoth;
  for (int r = 0; r < covRows; r++) {
    if (byIID.Integer(covIID[r]) == -1) byIID.SetInteger((const char*)covIID[r], r);
    if (covHasFID) {
      String key(covFID[r]);
      key.Add('_');
      key.Add(covIID[r]);
      if (byBoth.Integer(key) == -1) byBoth.SetInteger((const char*)key, r);
    }
  }   // End of r loop
  // row of the covariate file for every phenotype sample (-1: not found)
  IntArray row(idCount);
  row.Set(-1);
  int matched = 0;
  for (int i = 0; i < idCount; i++) {
    int k = -1;
    String both;
    if (pheHasFID) {
      both = familyID[i];
      both.Add('_');
      both.Add(sampleID[i]);
    }
    if (covHasFID && pheHasFID)         // FID and IID in both files: match on the pair
      k = byBoth.Integer(both);
    else {
      k = byIID.Integer(sampleID[i]);   // otherwise IID, then FID_IID
      if (k == -1 && pheHasFID) k = byIID.Integer(both);
    }
    row[i] = k;
    if (k > -1) matched++;
  }   // End of i loop
  if (matched == 0) {
    printf("None of the samples in the phenotype data have matching IDs in %s.\n", covarfile);
    printf("IDs in %s look like:\n", covarfile);
    for (int i = 0; i < 5 && i < covRows; i++)
      printf("  %s\n", (const char*)covIID[i]);
    printf("IDs in the phenotype data look like:\n");
    for (int i = 0; i < 5 && i < idCount; i++)
      printf("  %s\n", (const char*)sampleID[i]);
    printf("Please modify IDs in one of the input files and rerun again.\n");
    exit(0);
  }
  // columns to add; a name that is already present keeps the values of the earlier file
  IntArray add(0);
  for (int j = 0; j < covCount; j++) {
    if (alltraitNames.Find(covNames[j]) > -1)
      printf("WARNING: %s in %s is skipped because the name is already in an earlier file.\n", (const char*)covNames[j], covarfile);
    else
      add.Push(j);
  }   // End of j loop
  int oldCount = alltraitNames.Length();
  int addCount = add.Length();
  if (addCount > 0) {
    Matrix merged(idCount, oldCount + addCount);
    merged.Set(_NAN_);
    for (int i = 0; i < idCount; i++) {
      for (int t = 0; t < oldCount; t++)
        merged[i][t] = alltraits[i][t];
      if (row[i] > -1)
        for (int a = 0; a < addCount; a++)
          merged[i][oldCount + a] = C[row[i]][add[a]];
    }   // End of i loop
    alltraits = merged;
    for (int a = 0; a < addCount; a++)
      alltraitNames.Push(covNames[add[a]]);
  }
  printf("%d %s(s) for %d of %d samples loaded from %s.\n", addCount, unit, matched, idCount, covarfile);
}

void SECRET_Read(const char *phefile, const char *dictfile, const char *rsfile, const char *pcarsfile, const char *prefix, const char *covariateList,
    StringArray &sampleID, StringArray &alltraitNames, Matrix &alltraits, IntArray &covariates, IntArray &affections,
    IntArray &quantraits, StringArray &rsNames, Matrix &risks, StringArray &annotation, int coreCount)
{
    StringArray familyID;   // Step 1: Read in traits
    if (phefile[0])
        ReadTraits(phefile, alltraitNames, familyID, sampleID, alltraits, coreCount);
    if (phefile[0] && gCovarFile[0]) {  // covariates kept in their own file(s) are appended to the traits
        StringArray covarFiles;
        SplitFileList(gCovarFile, covarFiles);
        for (int f = 0; f < covarFiles.Length(); f++)
            MergeColumnFile(covarFiles[f], "covariate", alltraitNames, familyID, sampleID, alltraits, coreCount);
    }
    int alltraitCount = alltraitNames.Length();
    int idCount = sampleID.Length();
    StringArray covList;    // Step 2: Make covariates
    covariates.Dimension(0);
    if (covariateList[0])
        MakeCovariates(covariateList, alltraitNames, covList, covariates);
    if (alltraitCount) {
        MakeTraits(alltraits, alltraitNames, covariates, affections, quantraits);//Step 3
        IntArray rsKey;
        ReadRisks(rsfile, rsNames, familyID, sampleID, risks, rsKey, coreCount);//Step 4
        int keyCount = rsKey.Length();
        if (keyCount < idCount){
            bool hasFamilyID = (familyID.Length() == idCount);   // no FID column otherwise
            for (int i = 0; i < keyCount; i++) {
                int i0 = rsKey[i];
                if (hasFamilyID) familyID[i] = familyID[i0];
                sampleID[i] = sampleID[i0];
                for(int t = 0; t < alltraitCount; t++)
                    alltraits[i][t] = alltraits[i0][t];
            }   // End of i loop
            idCount = keyCount;
            if (hasFamilyID) familyID.Dimension(idCount);
            sampleID.Dimension(idCount);
            alltraits.Dimension(idCount, alltraitCount);
        }
    }
    else {  // no traits
        Matrix M;
        ReadScoreMatrix(rsfile, rsNames, sampleID, M, coreCount);
        int rsCount = M.cols;
        idCount = M.rows;
        risks.Dimension(rsCount, idCount);
        risks.Set(_NAN_);
#ifdef _OPENMP
#pragma omp parallel for num_threads(coreCount)
#endif
        for (int r = 0; r < rsCount; r++)
            for (int i = 0; i < idCount; i++)
                    risks[r][i] = M[i][r];
    }
    int rsCount = rsNames.Length();
    printf("%d risk scores for %d samples loaded.\n", rsCount, idCount);
    annotation.Dimension(0);    // Step 5: Read in risk score annotations
    if (dictfile[0]) {
        ReadKeys(dictfile, rsNames, annotation);
        printf("Data dictionary loaded.\n");
    }
    /*
    if (idCount && pcarsfile && pcarsfile[0]) { // Step 6: Read in pcars file, optional
        StringArray rowNames;
        Vector weight, pcars(idCount);
        pcars.Zero();
        ReadPCARS(pcarsfile, rowNames, weight, coreCount);
        int rowCount = weight.Length();
        for (int r = 0; r < rowCount; r++) {
            int k = rsNames.Find(rowNames[r]);
            if(k>-1)
                for (int i = 0; i < idCount; i++)
                    pcars[i] += risks[k][i] * weight[r];
        }   // End of r loop
        rsCount++;
        risks.Dimension(rsCount, idCount);
        risks[rsCount - 1] = pcars;
        rsNames.Add("PGS_WTSUM");
        if (annotation.Length()) annotation.Add("PGS_WTSUM");
    } */
    if (covList.Find("PC1") > -1)   // Step 7: Adjust PCs
        AdjustPC(risks, alltraits, alltraitNames, covList, coreCount);
}
