//////////////////////////////////////////////////////////////////////
// rplot.cpp
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

#include "MathMatrix.h"
#include "StringArray.h"
#include "IntArray.h"
#include "DataInput.h"
#include "Constant.h"
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>

// ---------------------------------------------------------------------
// PDF figures: R writes <stem>.ps, ps2pdf (Ghostscript) makes <stem>.pdf
// ---------------------------------------------------------------------
// Main.cpp calls CheckPs2pdf once; without a working ps2pdf the figures stay as .ps files and
// every message says so.
bool gPs2pdfOK = true;

static bool FileHasContent(const char *filename)
{
  FILE *fp = fopen(filename, "rb");
  if (fp == NULL) return false;
  fseek(fp, 0, SEEK_END);
  long size = ftell(fp);
  fclose(fp);
  return size > 0;
}

// Converts a tiny PostScript file, which tests ps2pdf and Ghostscript together.
bool CheckPs2pdf()
{
  const char *psfile = "secret_ps2pdf_test.ps", *pdffile = "secret_ps2pdf_test.pdf";
  FILE *fp = fopen(psfile, "wt");
  if (fp == NULL) { gPs2pdfOK = false; return false; }
  fputs("%!PS\n/Helvetica findfont 12 scalefont setfont 72 720 moveto (SECRET) show showpage\n", fp);
  fclose(fp);
  remove(pdffile);
  String command;
  command = "ps2pdf \"";
  command += psfile;
  command += "\" \"";
  command += pdffile;
  command += "\"";
#ifdef _WIN32
  command += " >NUL 2>&1";
#else
  command += " >/dev/null 2>&1";
#endif
  system((const char*)command);
  gPs2pdfOK = FileHasContent(pdffile);
  remove(psfile);
  remove(pdffile);
  return gPs2pdfOK;
}

// Makes <stem>.pdf from <stem>.ps.  Returns the extension of the figure that exists: ".pdf", or
// ".ps" when it could not be converted.  Use it in the message that tells where the figure is.
const char *MakePdf(const char *stem)
{
  if (!gPs2pdfOK) return ".ps";
  String pdffile(stem), command;
  pdffile.Add(".pdf");
  remove((const char*)pdffile);
  command = "ps2pdf -sPAPERSIZE=letter \"";
  command += stem;
  command += ".ps\" \"";
  command += stem;
  command += ".pdf\"";
  system((const char*)command);
  return FileHasContent((const char*)pdffile) ? ".pdf" : ".ps";
}

// ---------------------------------------------------------------------
// Which R runs the plots (--rpath)
// ---------------------------------------------------------------------
// --rpath is for an R that is not the first one in the PATH: the R executable, or the
// directory that holds it (e.g. /opt/R/4.3.1/bin or "C:/Program Files/R/R-4.3.1/bin").
// Main.cpp sets gRpath, which every plot uses when it is not given a path of its own.
const char *gRpath = "";

static bool IsDirectory(const char *path)
{
  struct stat st;
  return stat(path, &st) == 0 && (st.st_mode & S_IFDIR) != 0;
}

// The R executable to run: rpath (or gRpath when rpath is empty), or "R" from the PATH.
static String RExecutable(const char *rpath)
{
  const char *given = (rpath != NULL && rpath[0]) ? rpath : gRpath;
  String exe(given);
  if (exe.IsEmpty()) return String("R");
  while (exe.Length() > 1 && (exe[exe.Length() - 1] == '/' || exe[exe.Length() - 1] == '\\'))
    exe = exe.Left(exe.Length() - 1);                 // a directory may end with a separator
  if (IsDirectory(exe)) {
    exe += "/R";
#ifdef _WIN32
    exe += ".exe";
#endif
  }
  return exe;
}

// The command line that runs an R script: R CMD BATCH <script> <script>out, with every path
// quoted so that a directory name with spaces works.
String RCommandLine(const char *rpath, const char *scriptfile)
{
  String exe = RExecutable(rpath), command;
  if (exe == "R") command = "R";
  else { command = "\""; command += exe; command += "\""; }
  command += " CMD BATCH \"";
  command += scriptfile;
  command += "\" \"";
  command += scriptfile;
  command += "out\"";
#ifdef _WIN32
  command = "\"" + command + "\"";                   // cmd.exe removes the outer pair of quotes
#endif
  return command;
}

// Runs "R --version".  Returns true when R works, and puts its first line (e.g. "R version
// 4.3.1 (2023-06-16)") and the executable that was run in version / exe.
bool CheckRInstallation(const char *rpath, String &version, String &exe)
{
  exe = RExecutable(rpath);
  String command;
  if (exe == "R") command = "R";
  else { command = "\""; command += exe; command += "\""; }
  command += " --version 2>&1";
#ifdef _WIN32
  command = "\"" + command + "\"";
  FILE *fp = _popen((const char*)command, "r");
#else
  FILE *fp = popen((const char*)command, "r");
#endif
  version.Clear();
  if (fp == NULL) return false;
  char line[512];
  line[0] = '\0';
  if (fgets(line, sizeof(line), fp) != NULL) {
    int len = (int)strlen(line);
    while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r')) line[--len] = '\0';
    version = line;
  }
#ifdef _WIN32
  int status = _pclose(fp);
#else
  int status = pclose(fp);
#endif
  return status == 0 && version.Find("R version") > -1;
}

// R stopped because a package is missing: find the R error line (it starts with
// "Error" and says "is not installed", see the requireNamespace() checks before
// library()) and put its text, without the "Error: " prefix, into message.
// Returns false when the .Rout file has no such line.
static bool RMissingPackage(const char *scriptfile, char *message, int size)
{
    String outfile = scriptfile;
    outfile.Add("out");
    message[0] = '\0';
    FILE *fp = fopen((const char*)outfile, "rb");
    if (fp == NULL) return false;
    char line[512];
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "Error", 5) != 0 || strstr(line, "is not installed") == NULL) continue;
        char *start = strstr(line, "R package");
        if (start == NULL) start = line;
        int len = (int)strlen(start);
        while (len > 0 && (start[len-1] == '\n' || start[len-1] == '\r')) start[--len] = '\0';
        snprintf(message, size, "R plots were not made. %s\n", start);
        fclose(fp);
        return true;
    }
    fclose(fp);
    return false;
}

int CheckRout(const char *scriptfile)
{
    String outfile = scriptfile;
    outfile.Add("out");
    int errorFlag = 1;
    char buff[256];   // define the buffer and allocate the length
    FILE *fp = fopen((const char*)outfile, "rb");
    if (fp != NULL) {
        fseek(fp, -17, SEEK_END);// set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
        fread(buff, 16, 1, fp); // read the contents of the file starting from where fseek() positioned us
        buff[16] = '\0';        // close the string
        String lastline = buff;
        if (lastline.Find("Execution halted") == -1) errorFlag = 0;
        else {
            fseek(fp, -100, SEEK_END);  // set pointer to the end of file minus the length you need. Presumably there can be more than one new line caracter
            fread(buff, 80, 1, fp); // read the contents of the file starting from where fseek() positioned us
            buff[80] = '\0';           // close the string
            lastline = buff;
            char missing[256];
            if (lastline.Find("Error in library") > -1 || RMissingPackage(scriptfile, missing, sizeof(missing))) errorFlag = 3;
        }
        fclose(fp);                   // close the file
    }
    else errorFlag = 2;             // unable to open .Rout file
    return errorFlag;
}

// ---------------------------------------------------------------------
// R code shared by the plots of the top PGSs
// ---------------------------------------------------------------------

// Defines pgslabels(tab, ref) in the R script: the label of every row of a scan table
//     PGS000023:5(71%):trait      (SNPs used, coverage in parentheses; large counts as 100K, 1.2M)
// The SNPused / Coverage columns (written when --diction is a scores.coverage file) and the
// Annotation column are used when they exist, and a part that is NA for a PGS is left out.
// ref is the whole table, so that a Coverage column given as a fraction (all values <= 1) is
// recognised as such even when tab is only the top rows: it is then shown as a percentage.
void WritePGSLabelR(FILE *fp)
{
  fprintf(fp, "# ---- PGS label: PGS:<SNPs used>(<coverage>%%):<trait>, a part only when it is known ----\n");
  fputs("pgslabels <- function(tab, ref=tab){\n", fp);
  fputs("  n <- nrow(tab)\n", fp);
  fputs("  num <- function(x) suppressWarnings(as.numeric(as.character(x)))\n", fp);
  fputs("  cs <- 1\n", fp);
  fputs("  if(\"Coverage\" %in% names(ref)){ cr <- num(ref$Coverage); if(any(!is.na(cr)) && max(cr, na.rm=TRUE) <= 1) cs <- 100 }\n", fp);
  fputs("  used <- if(\"SNPused\" %in% names(tab)) num(tab$SNPused) else rep(NA_real_, n)\n", fp);
  fputs("  cov <- if(\"Coverage\" %in% names(tab)) num(tab$Coverage) * cs else rep(NA_real_, n)\n", fp);
  fputs("  fmtn <- function(x){\n", fp);
  fputs("    s <- signif(x, 3)\n", fp);
  fputs("    if(s >= 1e6) paste0(format(signif(x/1e6, 3), scientific=FALSE, trim=TRUE), \"M\")\n", fp);
  fputs("    else if(s >= 1e3) paste0(format(signif(x/1e3, 3), scientific=FALSE, trim=TRUE), \"K\")\n", fp);
  fputs("    else format(round(x), scientific=FALSE, trim=TRUE)\n", fp);
  fputs("  }\n", fp);
  fputs("  anno <- if(\"Annotation\" %in% names(tab)) as.character(tab$Annotation) else rep(NA_character_, n)\n", fp);
  fputs("  id <- as.character(tab$PGS)\n", fp);
  fputs("  anno[!is.na(anno) & anno == id] <- NA     # no annotation was found for this PGS: it is only its own name\n", fp);
  fputs("  vapply(seq_len(n), function(i){\n", fp);
  fputs("    st <- paste0(if(!is.na(used[i])) fmtn(used[i]) else \"\", if(!is.na(cov[i])) paste0(\"(\", sprintf(\"%.0f\", cov[i]), \"%)\") else \"\")\n", fp);
  fputs("    paste(c(id[i], if(nzchar(st)) st, if(!is.na(anno[i])) anno[i]), collapse=\":\")\n", fp);
  fputs("  }, \"\")\n", fp);
  fputs("}\n", fp);
}

// Defines heatcore(corrfile, ids, val, lab, kind, ttl, dirn) in the R script: a one-column
// heatmap of the top (at most 30) PGSs.  The rows are clustered by the correlation among the
// PGSs (<corrfile>, dendrogram on the left), lab is the row label on the right, and val is both
// the colour and the number:
//   kind "rsq"       Rsq, white to red from 0; dirn (Positive / Negative) is an annotation bar
//   kind "or", "hr"  odds / hazard ratio, coloured on the log scale: blue < 1 < red, white at 1
//   kind "auc"       AUC, white to red from 0.5
// Needs the pheatmap package; without it, or without the correlation file, nothing is drawn.
void WriteHeatmapCoreR(FILE *fp)
{
  fprintf(fp, "# ---- heatmap of the top PGSs: correlation cluster on the left, PGS labels on the right ----\n");
  fputs("heatcore <- function(corrfile, ids, val, lab, kind, ttl, dirn=NULL){\n", fp);
  fputs("  if(!file.exists(corrfile)) return(invisible(NULL))\n", fp);
  fputs("  if(!requireNamespace(\"pheatmap\", quietly=TRUE)){ message(\"R package pheatmap is not installed; the top PGS heatmap is skipped\"); return(invisible(NULL)) }\n", fp);
  fputs("  C <- as.matrix(read.table(corrfile, sep=\"\\t\", header=TRUE, row.names=1, check.names=FALSE, quote=\"\\\"\", comment.char=\"\"))\n", fp);
  fputs("  keep <- head(which(ids %in% rownames(C)), 30)\n", fp);
  fputs("  ids <- ids[keep]; val <- val[keep]; lab <- make.unique(lab[keep]); if(!is.null(dirn)) dirn <- dirn[keep]\n", fp);
  fputs("  n <- length(ids)\n", fp);
  fputs("  if(n < 2) return(invisible(NULL))\n", fp);
  fputs("  C <- C[ids, ids, drop=FALSE]\n", fp);
  fputs("  D <- 1 - abs(C); D[is.na(D)] <- 1; D <- (D + t(D))/2; diag(D) <- 0\n", fp);
  fputs("  hc <- hclust(as.dist(D), method=\"average\")\n", fp);
  fputs("  vname <- c(rsq=\"Rsq\", or=\"OR\", hr=\"HR\", auc=\"AUC\")[[kind]]\n", fp);
  fputs("  txt <- matrix(ifelse(is.finite(val), sprintf(\"%.3f\", val), \"\"), ncol=1)\n", fp);
  fputs("  M <- matrix(val, ncol=1, dimnames=list(lab, vname))\n", fp);
  fputs("  args <- list(M, cluster_rows=hc, cluster_cols=FALSE, display_numbers=txt, number_color=\"black\",\n", fp);
  fputs("    cellwidth=70, cellheight=12, fontsize_row=7, treeheight_row=70, border_color=\"grey70\")\n", fp);
  fputs("  if(kind %in% c(\"or\", \"hr\")){\n", fp);
  fputs("    lg <- log(val); lim <- max(abs(lg[is.finite(lg)]), 0)\n", fp);
  fputs("    if(!is.finite(lim) || lim <= 0) lim <- 0.1\n", fp);
  fputs("    M[, 1] <- lg; args[[1]] <- M\n", fp);
  fputs("    tk <- seq(-lim, lim, length.out=5)\n", fp);
  fputs("    args$breaks <- seq(-lim, lim, length.out=101)\n", fp);
  fputs("    args$color <- colorRampPalette(c(\"#2166AC\", \"white\", \"#B2182B\"))(100)\n", fp);
  fputs("    args$legend_breaks <- tk; args$legend_labels <- sprintf(\"%.2f\", exp(tk))\n", fp);
  fputs("  } else {\n", fp);
  fputs("    lo <- if(kind == \"auc\") min(0.5, val, na.rm=TRUE) else 0\n", fp);
  fputs("    top <- max(val, na.rm=TRUE); if(!is.finite(top) || top <= lo) top <- lo + 1\n", fp);
  fputs("    args$breaks <- seq(lo, top, length.out=101)\n", fp);
  fputs("    args$color <- colorRampPalette(c(\"white\", \"#FEE090\", \"#D73027\"))(100)\n", fp);
  fputs("  }\n", fp);
  fputs("  if(!is.null(dirn)){\n", fp);
  fputs("    args$annotation_row <- data.frame(Direction=factor(dirn, levels=c(\"Positive\", \"Negative\")), row.names=lab)\n", fp);
  fputs("    args$annotation_colors <- list(Direction=c(Positive=\"blue\", Negative=\"red\")); args$annotation_legend <- TRUE\n", fp);
  fputs("  }\n", fp);
  fputs("  args$main <- paste0(\"Top \", n, \" PGSs\", if(kind == \"auc\") \" by AUC\" else \"\", \" for \", ttl, \"\\n(rows clustered by PGS correlation)\")\n", fp);
  fputs("  do.call(pheatmap::pheatmap, args)\n", fp);
  fputs("  invisible(NULL)\n", fp);
  fputs("}\n", fp);
}

// ---------------------------------------------------------------------
// Heatmap of the top PGSs of a scan (one column): rows are clustered by the correlation
// among the PGSs (dendrogram on the left), the PGS label (see WritePGSLabelR) is on the
// right, and the effect is both the colour and the number.  Writes R code into an open
// script; needs the _PGScorr.txt file and the scan table (.linear / .logistic / .coxph).
// kind: TOPHEAT_RSQ  Rsq of --linear
//       TOPHEAT_OR   odds ratio of --logistic (colour on the log scale, centred at 1)
//       TOPHEAT_HR   hazard ratio of --coxph  (colour on the log scale, centred at 1)
// ---------------------------------------------------------------------
void WriteTopHeatmapR(FILE *fp, const char *corrfile, const char *statfile, int kind, const char *traitName)
{
  const char *kindName = (kind == TOPHEAT_OR) ? "or" : (kind == TOPHEAT_HR) ? "hr" : "rsq";
  WritePGSLabelR(fp);
  WriteHeatmapCoreR(fp);
  fprintf(fp, "topheat <- function(corrfile, statfile, kind, ttl){\n");
  fputs("  if(!file.exists(corrfile)) return(invisible(NULL))\n", fp);
  fputs("  tab <- read.table(statfile, sep=\"\\t\", header=TRUE, quote=\"\\\"\", comment.char=\"\", stringsAsFactors=FALSE)\n", fp);
  fputs("  zc <- if(\"Z\" %in% names(tab)) \"Z\" else \"z\"\n", fp);
  fputs("  val <- switch(kind,\n", fp);
  fputs("    rsq = tab$Rsq,\n", fp);
  fputs("    or  = if(\"OR\" %in% names(tab)) tab$OR else exp(tab$Beta),\n", fp);
  fputs("    hr  = if(\"exp.coef.\" %in% names(tab)) tab[[\"exp.coef.\"]] else exp(tab$coef))\n", fp);
  fputs("  dirn <- if(kind == \"rsq\") ifelse(tab[[zc]] < 0, \"Negative\", \"Positive\") else NULL\n", fp);
  fputs("  heatcore(corrfile, as.character(tab$PGS), val, pgslabels(tab), kind, ttl, dirn)\n", fp);
  fputs("  invisible(NULL)\n", fp);
  fputs("}\n", fp);
  fprintf(fp, "tryCatch(topheat(\"%s\", \"%s\", \"%s\", \"%s\"), error=function(e) message(\"top PGS heatmap skipped: \", conditionMessage(e)))\n", corrfile, statfile, kindName, traitName);
}

// ---------------------------------------------------------------------
// PGS correlation plots of the top 30 PGSs of a scan (<input>_PGScorrplot.pdf)
// ---------------------------------------------------------------------
// Every scan (--linear, --logistic, --coxph, --auc) has
//   page 1  correlations among the top 30 PGSs of all samples, rows and columns clustered
//           (average linkage on 1 - |r|).  When a --pruner2 was recommended (--linear,
//           --logistic), the dendrogram is cut into the clusters it was based on.
// --logistic adds, from the case and control correlations,
//   page 2  controls in the lower left triangle, cases in the upper right, in the order of
//           page 1 (no dendrogram: a clustering of all samples does not describe cases or
//           controls)
//   page 3  page 1 - page 2: the correlation in all samples minus the correlation in controls
//           (lower left) and minus the correlation in cases (upper right)
//   page 4  controls and cases, each with its own correlation cluster
//   page 5  case - control difference, in the order of page 1; * marks FDR < 0.05
// Reads <input>_PGScorr.txt (and, for the case-control pages, _PGScorr_control.txt,
// _PGScorr_case.txt, _PGScorr_diff.txt) and the scan table statfile (for the trait
// annotation of the labels).

extern double gPrunerR2Rec;
extern int gPGSClusterCount;

// R function pgscorr_pages(input, statfile, ttl, recR2, recK): the pages of one trait
static void WritePGSCorrPagesR(FILE *fp)
{
  fputs("rd <- function(f) as.matrix(read.table(f, sep=\"\\t\", header=TRUE, row.names=1, check.names=FALSE, quote=\"\\\"\", comment.char=\"\"))\n", fp);
  fputs("cl <- function(M){ D <- 1 - abs(M); D[is.na(D)] <- 1; D <- (D + t(D))/2; diag(D) <- 0; hclust(as.dist(D), method=\"average\") }\n", fp);
  fputs("cols <- colorRampPalette(c(\"navy\", \"white\", \"firebrick\"))(101); brk <- seq(-1, 1, length.out=102)\n", fp);
  fputs("nums <- function(M, star=NULL){ T <- matrix(ifelse(is.na(M), \"\", sprintf(\"%.2f\", M)), nrow(M)); if(!is.null(star)) T <- matrix(paste0(T, star), nrow(M)); T }\n", fp);
  fputs("tri <- function(Lo, Up, ord, lab, ids){   # lower left from Lo, upper right from Up, both in the order ord\n", fp);
  fputs("  A <- Lo[ord, ord]; B <- Up[ord, ord]\n", fp);
  fputs("  Tm <- matrix(NA_real_, length(ord), length(ord))\n", fp);
  fputs("  Tm[lower.tri(Tm)] <- A[lower.tri(A)]; Tm[upper.tri(Tm)] <- B[upper.tri(B)]\n", fp);
  fputs("  dimnames(Tm) <- list(lab[ord], ids[ord]); Tm }\n", fp);
  fputs("pgscorr_pages <- function(input, statfile, ttl, recR2=NA, recK=0){\n", fp);
  fputs("  corrfile <- paste0(input, \"_PGScorr.txt\")\n", fp);
  fputs("  if(!file.exists(corrfile)){ message(\"no PGS correlations for \", ttl); return(invisible(NULL)) }\n", fp);
  fputs("  Call <- rd(corrfile)\n", fp);
  fputs("  nAll <- nrow(Call)\n", fp);
  fputs("  ids <- head(rownames(Call), 30); n <- length(ids)          # the top 30 PGSs\n", fp);
  fputs("  if(n < 2) return(invisible(NULL))\n", fp);
  fputs("  Call <- Call[ids, ids, drop=FALSE]\n", fp);
  fputs("  tab <- if(file.exists(statfile)) read.table(statfile, sep=\"\\t\", header=TRUE, quote=\"\\\"\", comment.char=\"\", stringsAsFactors=FALSE, check.names=FALSE) else NULL\n", fp);
  fputs("  k <- if(!is.null(tab) && \"PGS\" %in% names(tab)) match(ids, as.character(tab$PGS)) else rep(NA_integer_, n)\n", fp);
  fputs("  an <- if(!is.null(tab) && \"Annotation\" %in% names(tab)) as.character(tab$Annotation)[k] else rep(NA_character_, n)\n", fp);
  fputs("  lab <- make.unique(ifelse(is.na(an) | !nzchar(an) | an == ids, ids, paste0(ids, \": \", an)))\n", fp);
  fputs("  fs <- if(n > 25) 4 else 5\n", fp);
  fputs("  hc <- cl(Call); ord <- hc$order\n", fp);
  fputs("  # page 1: all samples, clustered\n", fp);
  fputs("  Tm <- Call; dimnames(Tm) <- list(lab, ids)\n", fp);
  fputs("  args <- list(Tm, cluster_rows=hc, cluster_cols=hc, color=cols, breaks=brk, display_numbers=nums(Tm), number_color=\"black\",\n", fp);
  fputs("    fontsize_number=fs, fontsize_row=6, fontsize_col=6, na_col=\"grey85\", border_color=\"grey80\")\n", fp);
  fputs("  sub <- \"\"\n", fp);
  fputs("  if(!is.na(recR2) && recK >= 2 && recK < n && nAll == n){   # the clusters were found among these PGSs\n", fp);
  fputs("    args$cutree_rows <- recK; args$cutree_cols <- recK\n", fp);
  fputs("    sub <- paste0(\"\\ngaps: the \", recK, \" clusters behind the recommended --pruner2 \", sprintf(\"%.2f\", recR2))\n", fp);
  fputs("  }\n", fp);
  fputs("  args$main <- paste0(\"Top \", n, \" PGS correlations for \", ttl, \", all samples (rows and columns clustered)\", sub)\n", fp);
  fputs("  do.call(pheatmap, args)\n", fp);
  fputs("  ctlfile <- paste0(input, \"_PGScorr_control.txt\"); casfile <- paste0(input, \"_PGScorr_case.txt\"); difffile <- paste0(input, \"_PGScorr_diff.txt\")\n", fp);
  fputs("  if(!(file.exists(ctlfile) && file.exists(casfile) && file.exists(difffile))) return(invisible(NULL))\n", fp);
  fputs("  Ctl <- rd(ctlfile)[ids, ids, drop=FALSE]; Cas <- rd(casfile)[ids, ids, drop=FALSE]\n", fp);
  fputs("  dif <- read.table(difffile, sep=\"\\t\", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE, skip=1)\n", fp);
  fputs("  hdr <- readLines(difffile, n=1)\n", fp);
  fputs("  getv <- function(key) sub(paste0(\".*\", key, \"=([^ ]+).*\"), \"\\\\1\", hdr)\n", fp);
  fputs("  pJ <- getv(\"P\"); nCtl <- getv(\"NControl\"); nCas <- getv(\"NCase\")\n", fp);
  fputs("  glab <- paste0(\"Jennrich test of equal correlation matrices: P=\", pJ)\n", fp);
  fputs("  # page 2: lower left controls, upper right cases, in the order of page 1\n", fp);
  fputs("  pheatmap(tri(Ctl, Cas, ord, lab, ids), cluster_rows=FALSE, cluster_cols=FALSE, color=cols, breaks=brk,\n", fp);
  fputs("    display_numbers=nums(tri(Ctl, Cas, ord, lab, ids)), number_color=\"black\", fontsize_number=fs, fontsize_row=6, fontsize_col=6,\n", fp);
  fputs("    na_col=\"grey85\", border_color=\"grey80\",\n", fp);
  fputs("    main=paste0(\"Top \", n, \" PGS correlations for \", ttl, \"\\nlower left: controls (N=\", nCtl, \"); upper right: cases (N=\", nCas, \")\\nPGSs in the order of the clustering of all samples (page 1)\\n\", glab))\n", fp);
  fputs("  # page 3: page 1 - page 2\n", fp);
  fputs("  Dm <- tri(Call - Ctl, Call - Cas, ord, lab, ids)\n", fp);
  fputs("  lim <- max(abs(Dm), na.rm=TRUE); if(!is.finite(lim) || lim <= 0) lim <- 1\n", fp);
  fputs("  pheatmap(Dm, cluster_rows=FALSE, cluster_cols=FALSE, color=cols, breaks=seq(-lim, lim, length.out=102),\n", fp);
  fputs("    display_numbers=nums(Dm), number_color=\"black\", fontsize_number=fs, fontsize_row=6, fontsize_col=6,\n", fp);
  fputs("    na_col=\"grey85\", border_color=\"grey80\",\n", fp);
  fputs("    main=paste0(\"Correlation in all samples minus correlation in controls (lower left, N=\", nCtl, \") and in cases (upper right, N=\", nCas, \")\\nfor \", ttl, \", top \", n, \" PGSs in the order of page 1\"))\n", fp);
  fputs("  # page 4: controls and cases, each with its own cluster\n", fp);
  fputs("  mk <- function(M, nm, N){ rownames(M) <- colnames(M) <- ids\n", fp);
  fputs("    pheatmap(M, cluster_rows=cl(M), cluster_cols=cl(M), color=cols, breaks=brk, display_numbers=nums(M), number_color=\"black\",\n", fp);
  fputs("      fontsize_number=fs-1, fontsize_row=6, fontsize_col=6, na_col=\"grey85\", border_color=\"grey80\", legend=FALSE,\n", fp);
  fputs("      main=paste0(nm, \" (N=\", N, \")\"), silent=TRUE)$gtable }\n", fp);
  fputs("  g1 <- mk(Ctl, \"Controls\", nCtl); g2 <- mk(Cas, \"Cases\", nCas)\n", fp);
  fputs("  grid.newpage()\n", fp);
  fputs("  pushViewport(viewport(layout=grid.layout(2, 2, heights=unit(c(1, 12), \"null\"))))\n", fp);
  fputs("  grid.text(paste0(\"PGS correlations among the top \", n, \" PGSs for \", ttl, \", each panel clustered separately\\n\", glab),\n", fp);
  fputs("    vp=viewport(layout.pos.row=1, layout.pos.col=1:2), gp=gpar(fontsize=11, fontface=\"bold\"))\n", fp);
  fputs("  pushViewport(viewport(layout.pos.row=2, layout.pos.col=1)); grid.draw(g1); popViewport()\n", fp);
  fputs("  pushViewport(viewport(layout.pos.row=2, layout.pos.col=2)); grid.draw(g2); popViewport()\n", fp);
  fputs("  popViewport()\n", fp);
  fputs("  # page 5: case - control difference (order of page 1); * = FDR < 0.05\n", fp);
  fputs("  Dd <- (Cas - Ctl)[ord, ord]\n", fp);
  fputs("  St <- matrix(\"\", n, n, dimnames=list(ids[ord], ids[ord]))\n", fp);
  fputs("  sig <- dif[!is.na(dif$P_BH) & dif$P_BH < 0.05, , drop=FALSE]\n", fp);
  fputs("  for(i in seq_len(nrow(sig))){ a <- match(sig$PGS1[i], ids[ord]); b <- match(sig$PGS2[i], ids[ord]); if(!is.na(a) && !is.na(b)) St[a, b] <- St[b, a] <- \"*\" }\n", fp);
  fputs("  lim <- max(abs(Dd), na.rm=TRUE); if(!is.finite(lim) || lim <= 0) lim <- 1\n", fp);
  fputs("  dimnames(Dd) <- list(lab[ord], ids[ord])\n", fp);
  fputs("  pheatmap(Dd, cluster_rows=FALSE, cluster_cols=FALSE, color=cols, breaks=seq(-lim, lim, length.out=102),\n", fp);
  fputs("    display_numbers=nums(Dd, St), number_color=\"black\", fontsize_number=fs, fontsize_row=6, fontsize_col=6,\n", fp);
  fputs("    na_col=\"grey85\", border_color=\"grey80\",\n", fp);
  fputs("    main=paste0(\"Case - control PGS correlation for \", ttl, \" (* FDR < 0.05)\\n\", glab))\n", fp);
  fputs("  invisible(NULL)\n", fp);
  fputs("}\n", fp);
}

// Writes the R script that draws the pages of the given (input, title) pairs into
// outfile.ps, runs R and makes outfile.pdf.
static int RunPGSCorrPlots(const char *outfile, const char *statfile, StringArray &inputs, StringArray &titles,
  bool recommended, const char *rpath, String &msg)
{
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET PGS correlation plots, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "statfile <- \"%s\"\n", statfile);
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", outfile);
  fputs("if(!requireNamespace(\"pheatmap\", quietly=TRUE)) stop(\"R package pheatmap is not installed. Please install it in R with: install.packages('pheatmap')\", call.=FALSE)\n", fp);
  fputs("library(pheatmap); library(grid)\n", fp);
  WritePGSCorrPagesR(fp);
  for (int i = 0; i < inputs.Length(); i++) {
    if (recommended && gPrunerR2Rec != _NAN_ && gPGSClusterCount >= 2)
      fprintf(fp, "tryCatch(pgscorr_pages(\"%s\", statfile, \"%s\", %.2f, %d), error=function(e) message(\"PGS correlation plots skipped: \", conditionMessage(e)))\n",
        (const char*)inputs[i], (const char*)titles[i], gPrunerR2Rec, gPGSClusterCount);
    else
      fprintf(fp, "tryCatch(pgscorr_pages(\"%s\", statfile, \"%s\"), error=function(e) message(\"PGS correlation plots skipped: \", conditionMessage(e)))\n",
        (const char*)inputs[i], (const char*)titles[i]);
  }   // End of i loop
  fputs("dev.off()\n", fp);
  fputs("rm(list=ls())\n", fp);
  fputs("rm(list=ls(all.names=TRUE))\n", fp);
  fputs("q(save=\"no\")\n", fp);
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  message[0]='\0';
  int error = CheckRout(scriptfile);
  switch (error) {
    case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      break;
    case 1:  sprintf(message, "R code %s failed.\n  Please check %sout for details.\n", (const char*)scriptfile, (const char*)scriptfile);
      break;
    case 3:  if (!RMissingPackage(scriptfile, message, sizeof(message)))
          sprintf(message, "R plots were not made: an R package is not installed. Please check %sout for details.\n", (const char*)scriptfile);
          break;
    case 0:
      sprintf(message, "PGS correlation plots saved in %s%s\n", outfile, MakePdf(outfile));
      break;
    default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}

// One trait of --linear, --logistic or --coxph: <input>_PGScorrplot.pdf.
// statfile is the scan table (.linear / .logistic / .coxph).
int plotPGSCorr(const char *input, const char *statfile, const char *traitName, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_PGScorrplot");
  StringArray inputs, titles;
  inputs.Push(String(input));
  titles.Push(String(traitName));
  return RunPGSCorrPlots(outfile, statfile, inputs, titles, true, rpath, msg);
}

// --auc: one page per column of the AUC table, in <outprefix>_PGScorrplot.pdf.  inputs are
// the <prefix>_<column>_AUC base names of the correlation files.
int plotPGSCorrList(const char *outprefix, const char *statfile, StringArray &inputs, StringArray &titles,
  const char *rpath, String &msg)
{
  String outfile = outprefix;
  outfile.Add("_PGScorrplot");
  return RunPGSCorrPlots(outfile, statfile, inputs, titles, false, rpath, msg);
}

int plotBestHeatmap(const char *input, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("plot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Best PGS heapmap plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "data <- read.table(\"%s.txt\", sep=\"\\t\", header = TRUE, row.names=1)\n", input);
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  fprintf(fp, "if(!requireNamespace(\"pheatmap\", quietly=TRUE)) stop(\"R package pheatmap is not installed. Please install it in R with: install.packages('pheatmap')\", call.=FALSE)\n");
  fprintf(fp, "library(pheatmap)\n");
  fprintf(fp, "pheatmap(data, main=\"Heatmap Plot for R-Squared of the Best PGSs\")\n");
  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  int error = CheckRout(scriptfile);
  switch (error) {
  case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
    break;
  case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
    sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
    break;
  case 3:  if (!RMissingPackage(scriptfile, message, sizeof(message)))
        sprintf(message, "R plots were not made: an R package is not installed. Please check %sout for details.\n", (const char*)scriptfile);
        break;
  case 0:
    sprintf(message, "Rsq plots saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
    break;
  default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}

int plotCalibration(const char *input, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_calibrationplot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Calibration plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  fprintf(fp, "data <- read.table(\"%s_combinedPGS.txt\", sep=\"\\t\", header = TRUE, row.names=1)\n", input);
  fprintf(fp, "N.r <- dim(data)[2]-1\n");
  fprintf(fp, "x <- as.matrix(data[,-1])\n");
  fprintf(fp, "y <- data[,1]\n");
  fprintf(fp, "for(i in 1:N.r){\n");
  fprintf(fp, "model <- lm(y ~ x[,i], na.action=\"na.exclude\")\n");
  fprintf(fp, "y.fit <- fitted(model)\n");
  fprintf(fp, "corr <- cor(y, y.fit, use=\"complete.obs\")\n");
  fprintf(fp, "rsq <- round(corr * corr, digit=3)\n");
  fprintf(fp, "model2 <- lm(y ~ y.fit)\n");
  fprintf(fp, "plot(y.fit, y, xlab=paste(\"Predicted\", names(data)[1]), ylab=paste(\"Observed\", names(data)[1]), main=paste0(names(data)[i+1], \" (R-squared=\", rsq, \")\"))\n");
  fprintf(fp, "abline(a=model2$coefficient[1], b=model2$coefficient[2])\n");
  fprintf(fp, "}\n");
  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  message[0]='\0';
  int error = CheckRout(scriptfile);
  switch (error) {
    case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      break;
    case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
      break;
    case 0:
      sprintf(message, "Calibration plots saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
      break;
    default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}

int plotLogistic(const char *input, const char *traitName, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_logisticplot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Rsq plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  WritePGSLabelR(fp);
  fprintf(fp, "full <- read.table(\"%s.logistic\", sep=\"\\t\", header = TRUE, quote=\"\\\"\", comment.char=\"\")\n", input);
  fprintf(fp, "data <- head(full, 30)   # the top 30 PGSs; the whole table decides how Coverage is shown\n");
  fprintf(fp, "chisq <- data$Z * data$Z\n");
  fprintf(fp, "bestname <- as.character(data$PGS[1])   # the best PGS: first in the scan table, column BEST of the _PGS.txt file\n");
  fprintf(fp, "n <- nrow(data)   # at most 30 PGSs; fewer when fewer were fitted\n");
  fprintf(fp, "plot(chisq, 1:n, pch=19, col=c(\"blue\", \"red\")[(data$Z<0)+1],\n");
  fprintf(fp, "xlab = \"Z-Squared between PGS and %s\", ylab = \"PGS Rank (by Z-Squared)\", main = paste0(\"Top \", n, \" PGS Scores That Predict %s\"))\n", traitName, traitName);
  fprintf(fp, "mid <- (min(chisq) + max(chisq))/2\n");
  fprintf(fp, "labs <- pgslabels(data, full)   # PGS:SNPs:Coverage:trait, as far as known\n");
  fprintf(fp, "text(chisq, 1:n, labels = labs, pos = ifelse(chisq < mid, 4, 2))\n");
  { String corrfile = input; corrfile.Add("_PGScorr.txt"); String statfile = input; statfile.Add(".logistic");
    WriteTopHeatmapR(fp, corrfile, statfile, TOPHEAT_OR, traitName); }
  
  fprintf(fp, "# ---- combined PGS of the models (written by --logistic --rplot): trait ~ score ----\n");
  fprintf(fp, "pgsfile <- \"%s_PGS.txt\"\n", input);
  fprintf(fp, "if(file.exists(pgsfile)){\n");
  fprintf(fp, "  data <- read.table(pgsfile, sep=\"\\t\", header=TRUE, row.names=1)\n");
  fprintf(fp, "  y <- data[,1]\n");
  fprintf(fp, "  wilson <- function(k, n){\n");
  fprintf(fp, "    p <- k/n; z <- 1.96; d <- 1 + z^2/n\n");
  fprintf(fp, "    h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d\n");
  fprintf(fp, "    cbind((p + z^2/(2*n))/d - h, (p + z^2/(2*n))/d + h)\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  aucfun <- function(score, y){\n");
  fprintf(fp, "    r <- rank(score); n1 <- sum(y==1); n0 <- sum(y==0)\n");
  fprintf(fp, "    (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0)\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  rocfun <- function(score, y){\n");
  fprintf(fp, "    o <- order(score, decreasing=TRUE); y <- y[o]\n");
  fprintf(fp, "    list(fpr=c(0, cumsum(1-y)/sum(1-y)), tpr=c(0, cumsum(y)/sum(y)))\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  usable <- function(score, y){\n");
  fprintf(fp, "    ok <- !is.na(score)\n");
  fprintf(fp, "    sum(ok) >= 20 && sum(y[ok]==1) >= 2 && sum(y[ok]==0) >= 2 && sd(score[ok]) > 0\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  combined <- c(\"BEST\", \"WTSUM\", \"RIDGE\", \"NET\", \"LASSO\")\n");
  fprintf(fp, "  labels <- c(\"Best\", \"WTSUM\", \"Ridge\", \"Elastic net\", \"Lasso\")\n");
  fprintf(fp, "  pgsname <- function(k) if(k==1) paste0(\"Best PGS (\", bestname, \")\") else paste0(labels[k], \" PGS\")\n");
  fprintf(fp, "  for(k in 1:5){\n");
  fprintf(fp, "    if(!(combined[k] %%in%% names(data))) next\n");
  fprintf(fp, "    score <- data[[combined[k]]]\n");
  fprintf(fp, "    if(!usable(score, y)){\n");
  fprintf(fp, "      par(mfrow=c(1,1))\n");
  fprintf(fp, "      plot.new()\n");
  fprintf(fp, "      text(0.5, 0.5, paste0(pgsname(k), \" on %s is not available or is constant\\n(no PGS left in the model)\"))\n", traitName);
  fprintf(fp, "      next\n");
  fprintf(fp, "    }\n");
  fprintf(fp, "    ok <- !is.na(score); s <- score[ok]; yy <- y[ok]\n");
  fprintf(fp, "    nb <- if(length(s) >= 500) 10 else 5\n");
  fprintf(fp, "    brk <- unique(quantile(s, probs=seq(0, 1, length.out=nb+1)))\n");
  fprintf(fp, "    bin <- cut(s, breaks=brk, include.lowest=TRUE)\n");
  fprintf(fp, "    nk <- as.vector(tapply(yy, bin, length)); kk <- as.vector(tapply(yy, bin, sum))\n");
  fprintf(fp, "    px <- as.vector(tapply(s, bin, mean))\n");
  fprintf(fp, "    ci <- wilson(kk, nk)\n");
  fprintf(fp, "    fit <- tryCatch(suppressWarnings(glm(yy ~ s, family=binomial)), error=function(e) NULL)\n");
  fprintf(fp, "    orsd <- if(is.null(fit)) NA else exp(coef(fit)[2]*sd(s))\n");
  fprintf(fp, "    par(mfrow=c(1,2))\n");
  fprintf(fp, "    plot(px, kk/nk, ylim=c(0, max(ci[,2])), pch=19, xlab=pgsname(k),\n");
  fprintf(fp, "      ylab=\"Proportion affected (95%% CI)\", cex.main=0.9,\n");
  fprintf(fp, "      main=paste0(pgsname(k), \" on %s\\n(AUC=\", sprintf(\"%%.3f\", aucfun(s, yy)),\n", traitName);
  fprintf(fp, "      \", OR per SD=\", sprintf(\"%%.2f\", orsd), \")\"))\n");
  fprintf(fp, "    segments(px, ci[,1], px, ci[,2])\n");
  fprintf(fp, "    if(!is.null(fit)){\n");
  fprintf(fp, "      xs <- seq(min(s), max(s), length.out=100)\n");
  fprintf(fp, "      lines(xs, predict(fit, data.frame(s=xs), type=\"response\"))\n");
  fprintf(fp, "    }\n");
  fprintf(fp, "    abline(h=mean(yy), lty=3)\n");
  fprintf(fp, "    dc <- density(s[yy==0]); da <- density(s[yy==1])\n");
  fprintf(fp, "    plot(range(dc$x, da$x), c(0, max(dc$y, da$y)), type=\"n\", xlab=pgsname(k),\n");
  fprintf(fp, "      ylab=\"Density\", cex.main=0.9, main=paste0(pgsname(k), \" in affected and unaffected\"))\n");
  fprintf(fp, "    lines(dc, col=\"blue\"); lines(da, col=\"red\")\n");
  fprintf(fp, "    legend(\"topleft\", legend=c(paste0(\"Unaffected (n=\", sum(yy==0), \")\"),\n");
  fprintf(fp, "      paste0(\"Affected (n=\", sum(yy==1), \")\")), col=c(\"blue\", \"red\"), lty=1, bty=\"n\")\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  par(mfrow=c(1,1))\n");
  fprintf(fp, "  cols <- c(\"BEST\", \"WTSUM\", \"RIDGE\", \"NET\", \"LASSO\")\n");
  fprintf(fp, "  lab <- c(\"Best PGS\", \"WTSUM\", \"Ridge\", \"Elastic net\", \"Lasso\")\n");
  fprintf(fp, "  colr <- c(\"grey40\", \"black\", \"blue\", \"darkgreen\", \"red\")\n");
  fprintf(fp, "  lty <- c(2, 1, 1, 1, 1)\n");
  fprintf(fp, "  plot(c(0,1), c(0,1), type=\"n\", xlab=\"False positive rate (1 - specificity)\",\n");
  fprintf(fp, "    ylab=\"True positive rate (sensitivity)\", main=\"ROC curves of PGSs predicting %s (score only)\")\n", traitName);
  fprintf(fp, "  abline(0, 1, lty=3, col=\"grey\")\n");
  fprintf(fp, "  leg <- c(); keep <- c()\n");
  fprintf(fp, "  for(k in 1:5){\n");
  fprintf(fp, "    if(!(cols[k] %%in%% names(data))) next\n");
  fprintf(fp, "    sc <- data[[cols[k]]]\n");
  fprintf(fp, "    if(!usable(sc, y)) next\n");
  fprintf(fp, "    ok <- !is.na(sc)\n");
  fprintf(fp, "    r <- rocfun(sc[ok], y[ok])\n");
  fprintf(fp, "    lines(r$fpr, r$tpr, col=colr[k], lty=lty[k], lwd=2)\n");
  fprintf(fp, "    leg <- c(leg, paste0(lab[k], \" (AUC=\", sprintf(\"%%.3f\", aucfun(sc[ok], y[ok])), \")\"))\n");
  fprintf(fp, "    keep <- c(keep, k)\n");
  fprintf(fp, "  }\n");
  fprintf(fp, "  if(length(keep)) legend(\"bottomright\", legend=leg, col=colr[keep], lty=lty[keep], lwd=2, bty=\"n\")\n");
  fprintf(fp, "}\n");

  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  message[0]='\0';
  int error = CheckRout(scriptfile);
  switch (error) {
    case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      break;
    case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
      break;
    case 0:
      sprintf(message, "PGS scan plots saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
      break;
    default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}

int plotRSQ(const char *input, const char *traitName, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_rsqplot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Rsq plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  WritePGSLabelR(fp);
  fprintf(fp, "full <- read.table(\"%s.linear\", sep=\"\\t\", header = TRUE, quote=\"\\\"\", comment.char=\"\")\n", input);
  fprintf(fp, "data <- head(full, 30)   # the top 30 PGSs; the whole table decides how Coverage is shown\n");
  fprintf(fp, "n <- nrow(data)   # at most 30 PGSs; fewer when fewer were scanned\n");
  fprintf(fp, "plot(data$Rsq, 1:n, pch=19, col=c(\"blue\", \"red\")[(data$Z<0)+1],\n");
  fprintf(fp, "xlab = \"R-Squared of Association between PGS and %s\", ylab = \"PGS Rank (by R-Squared)\", main = paste0(\"Top \", n, \" PGS Scores That Predict %s\"))\n", traitName, traitName);
  fprintf(fp, "mid <- (min(data$Rsq) + max(data$Rsq))/2\n");
  fprintf(fp, "labs <- pgslabels(data, full)   # PGS:SNPs:Coverage:trait, as far as known\n");
  fprintf(fp, "text(data$Rsq, 1:n, labels = labs, pos = ifelse(data$Rsq < mid, 4, 2))\n");
  { String corrfile = input; corrfile.Add("_PGScorr.txt"); String statfile = input; statfile.Add(".linear");
    WriteTopHeatmapR(fp, corrfile, statfile, TOPHEAT_RSQ, traitName); }

  fprintf(fp, "data <- read.table(\"%s_PGS.txt\", sep=\"\\t\", header = TRUE, row.names=1)\n", input);
  fprintf(fp, "# columns: trait, top PGSs, WTSUM, then RIDGE / NET / LASSO combined scores (found by name)\n");
  fprintf(fp, "w <- which(names(data)==\"WTSUM\")[1]\n");
  fprintf(fp, "N.r <- w-2\n");
  fprintf(fp, "x <- as.matrix(data[,2:(w-1), drop=FALSE])\n");
  fprintf(fp, "y <- data[,1]\n");
  fprintf(fp, "if(!requireNamespace(\"pheatmap\", quietly=TRUE)) stop(\"R package pheatmap is not installed. Please install it in R with: install.packages('pheatmap')\", call.=FALSE)\n");
  fprintf(fp, "library(pheatmap)\n");
  fprintf(fp, "order.y <- order(y)\n");
  fprintf(fp, "index.ordered <- order.y[(1:50)*floor(length(y)/50)]\n");
  fprintf(fp, "x.ordered <- t(x[index.ordered, 1:min(30,N.r), drop=FALSE])\n");
  fprintf(fp, "colnames(x.ordered) <- paste0(\"Q\", (1:50)*2)\n");
  fprintf(fp, "if(N.r >= 2) pheatmap(x.ordered, cluster_cols=FALSE, main=paste0(\"Heatmap Plot for Top \", min(30,N.r), \" Polygenic Scores Across 50 %s-Sorted Individuals\"))   # a single PGS cannot be clustered\n", traitName);
  fprintf(fp, "model <- lm(y ~ data[,2])\n");
  fprintf(fp, "rsq <- round(summary(model)$r.squared, digits=3)\n");
  fprintf(fp, "plot(data[,2], y, xlab=names(data)[2], ylab=\"%s\", main=paste0(names(data)[2], \" on %s (R-squared=\", rsq, \")\"))\n", traitName, traitName);
  fprintf(fp, "abline(a=model$coefficients[1], b=model$coefficients[2])\n");
  fprintf(fp, "# trait ~ combined PGS for the weighted sum and the 3 penalized regressions\n");
  fprintf(fp, "combined <- c(\"WTSUM\", \"RIDGE\", \"NET\", \"LASSO\")\n");
  fprintf(fp, "labels <- c(\"WTSUM\", \"Ridge\", \"Elastic net\", \"Lasso\")\n");
  fprintf(fp, "for(k in 1:4){\n");
  fprintf(fp, "  if(!(combined[k] %%in%% names(data))) next\n");
  fprintf(fp, "  score <- data[[combined[k]]]\n");
  fprintf(fp, "  flat <- !(sd(score, na.rm=TRUE) > 0)   # e.g., lasso that shrinks every PGS to zero\n");
  fprintf(fp, "  rsq <- if(flat) 0 else round(summary(lm(y ~ score))$r.squared, digits=3)\n");
  fprintf(fp, "  plot(score, y, xlab=paste(labels[k], \"PGS\"), ylab=\"%s\", main=paste0(labels[k], \" PGS on %s (R-squared=\", rsq, \")\"))\n", traitName, traitName);
  fprintf(fp, "  if(!flat){model <- lm(y ~ score); abline(a=model$coefficients[1], b=model$coefficients[2])}\n");
  fprintf(fp, "}\n");

  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  message[0]='\0';
  int error = CheckRout(scriptfile);
  switch (error) {
    case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      break;
    case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
      break;
    case 3:  if (!RMissingPackage(scriptfile, message, sizeof(message)))
          sprintf(message, "R plots were not made: an R package is not installed. Please check %sout for details.\n", (const char*)scriptfile);
          break;
    case 0:
      sprintf(message, "Rsq plots saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
      break;
    default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}

// The columns of an AUC table (.auc): PGS, one column per trait, then SNPused and Coverage
// (only with a coverage dictionary) and Annotation (only with a dictionary).  Defines
// auctraits(data), the numbers of the trait columns, in the R script.
static void WriteAUCTraitColumnsR(FILE *fp)
{
    fputs("auctraits <- function(data){\n", fp);
    fputs("  nm <- names(data)\n", fp);
    fputs("  extra <- nm %in% c(\"PGS\", \"Annotation\")\n", fp);
    fputs("  if(all(c(\"SNPused\", \"Coverage\") %in% nm)) extra <- extra | nm %in% c(\"SNPused\", \"Coverage\")\n", fp);
    fputs("  which(!extra)\n", fp);
    fputs("}\n", fp);
}

int plotAUCs(const char *input, const char *output, const char *rpath)
{
    String scriptfile = output;
    scriptfile.Add("_aucplot.R");
    FILE *fp = fopen(scriptfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
    fprintf(fp, "## %s for SECRET AUC plot, by Wei-Min Chen\n", (const char*)scriptfile);
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fprintf(fp, "data <- read.table(\"%s\", sep=\"\\t\", header = TRUE, quote=\"\\\"\", comment.char=\"\", check.names=FALSE)\n", input);
    fprintf(fp, "outprefix <- \"%s\"\n", output);
    fprintf(fp, "pgsCount <- dim(data)[1]\n");
    fprintf(fp, "postscript(\"%s_aucplot.ps\", paper = \"letter\", horizontal = T)\n", output);
    WritePGSLabelR(fp);
    WriteHeatmapCoreR(fp);
    WriteAUCTraitColumnsR(fp);
    // input is the .auc file of one trait (<prefix>_<trait>.auc, the largest AUC first) and output is
    // <prefix>_<trait>.  The figures: the top 1-30 PGSs as a dot plot, then the heatmap of the top 30
    // (rows clustered by the PGS correlations of <output>_AUC_PGScorr.txt)
    fputs("for(trait in auctraits(data)){\n", fp);
    fputs("  x <- data[, trait]\n", fp);
    fputs("  O <- order(-x)\n", fp);
    fputs("  x.ordered <- x[O]\n", fp);
    fputs("  name <- names(data)[trait]\n", fp);
    fputs("  nOK <- sum(!is.na(x))\n", fp);
    fputs("  if(nOK == 0) next\n", fp);
    fputs("  range <- 1:min(30, nOK)   # the top 30 PGSs\n", fp);
    fputs("  plot(x.ordered[range], range, pch=19, xlab = paste(\"AUC of Predicting\", name),\n", fp);
    fputs("   ylab = \"PGS Rank (by AUC)\", main = paste0(\"Top 1-\", max(range), \" PGS Scores That Predict \", name))\n", fp);
    fputs("  mid <- (min(x.ordered[range]) + max(x.ordered[range]))/2\n", fp);
    fputs("  labs <- pgslabels(data[O[range], , drop=FALSE], data)   # PGS:SNPs:Coverage:trait, as far as known\n", fp);
    fputs("  text(x.ordered[range], range, labels = labs, pos = ifelse(x.ordered[range] < mid, 4, 2))\n", fp);
    fputs("  top <- head(O[!is.na(x[O])], 30)\n", fp);
    fputs("  tryCatch(heatcore(paste0(outprefix, \"_AUC_PGScorr.txt\"), as.character(data$PGS[top]), x[top],\n", fp);
    fputs("    pgslabels(data[top, , drop=FALSE], data), \"auc\", name),\n", fp);
    fputs("    error=function(e) message(\"top PGS heatmap skipped: \", conditionMessage(e)))\n", fp);
    fputs("}\n", fp);
    fputs("dev.off()\n", fp);
    fputs("rm(list=ls(all=TRUE))\n", fp);
    fclose(fp);
    char command[4096];
    snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
    system(command);
    int error = CheckRout(scriptfile);
    switch (error) {
    case 2:  printf("R code %s failed.\n", (const char*)scriptfile);
        break;
    case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
        printf("  Please check %sout for details.\n", (const char*)scriptfile);
        break;
    case 0:
        { String stem(output); stem.Add("_aucplot");
          printf("AUC plots saved in %s%s\n", (const char*)stem, MakePdf(stem)); }
        break;
    default: printf("Unexpected error. Please contact SECRET authors.\n");
    }
    return(error);
}

// The AUCs of one trait in the two strata of --strat, from <prefix>_<strat>_1.auc and
// <prefix>_<strat>_2.auc, where prefix is <prefix>_<trait>; output <prefix>_<strat>_aucplot.
int plotTuneTraitAUCs(const char *prefix, const char *tuneTrait, const char *rpath)
{
    String output(prefix);
    output.Add("_");
    output.Add(tuneTrait);
    String inputs[2];
    inputs[0] = output;
    inputs[0].Add("_1.auc");
    inputs[1] = output;
    inputs[1].Add("_2.auc");
    String scriptfile = output;
    scriptfile.Add("_aucplot.R");
    FILE *fp = fopen(scriptfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
    fprintf(fp, "## %s for SECRET AUC plot, by Wei-Min Chen\n", (const char*)scriptfile);
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fprintf(fp, "data1 <- read.table(\"%s\", sep=\"\\t\", header = TRUE, quote=\"\\\"\", comment.char=\"\", check.names=FALSE)\n", (const char*)inputs[0]);
    fprintf(fp, "data2 <- read.table(\"%s\", sep=\"\\t\", header = TRUE, quote=\"\\\"\", comment.char=\"\", check.names=FALSE)\n", (const char*)inputs[1]);
    fprintf(fp, "pgsCount <- dim(data1)[1]\n");
    fprintf(fp, "postscript(\"%s_aucplot.ps\", paper = \"letter\", horizontal = T)\n", (const char*)output);
    WritePGSLabelR(fp);
    WriteAUCTraitColumnsR(fp);
    fputs("dict <- \"Annotation\" %in% names(data1)\n", fp);
    fputs("for(trait in auctraits(data1)){\n", fp);
    fputs("  x <- data1[, trait]\n", fp);
    fputs("  y <- data2[match(data1$PGS, data2$PGS), trait]   # the two files are sorted by their own AUCs\n", fp);
    fputs("  cand <- function(v){ z <- (x > v) | (y > v); z[is.na(z)] <- FALSE; z }\n", fp);
    fputs("  vT <- 0.8\n", fp);
    fputs("  candidate <- cand(vT)\n", fp);
    fputs("  while(sum(candidate) < 20 & vT > 0.56){vT <- vT - 0.01\n", fp);
    fputs("  candidate <- cand(vT)}\n", fp);
    fputs("  name <- names(data1)[trait]\n", fp);
    fputs("  if(sum(candidate) == 0){   # no PGS stands out: say so, a plot of nothing cannot be drawn\n", fp);
    fputs("    plot.new()\n", fp);
    fputs("    title(main = paste(\"Top PGS Scores That Predict\", name))\n", fp);
    fputs("    text(0.5, 0.5, paste(\"No PGS has an AUC above\", round(vT, 2), \"for\", name, \"in either stratum\"))\n", fp);
    fputs("    next\n", fp);
    fputs("  }\n", fp);
    fprintf(fp, "  plot(x[candidate], y[candidate], pch = 19, xlab = paste(\"AUC of Predicting\", name, \"for %s=1\"),\n", tuneTrait);
    fprintf(fp, "  ylab = paste(\"AUC of Predicting\", name, \"for %s=2\"),\n", tuneTrait);
    fputs("  main = paste(\"Top PGS Scores That Predict\", name))\n", fp);
    fputs("  middle <- (min(x[candidate]) + max(x[candidate]))/2\n", fp);
    fputs("  if(sum(candidate)>0){abline(a=0,b=1)\n", fp);
    fputs("  labs <- pgslabels(data1[candidate, , drop=FALSE], data1)   # PGS:SNPs:Coverage:trait, as far as known\n", fp);
    fputs("  text(x = x[candidate], y = y[candidate], pos = (x[candidate]<middle)*2+2, labels = labs)}\n", fp);
    fputs("}\n", fp);
    fputs("dev.off()\n", fp);
    fputs("rm(list=ls(all=TRUE))\n", fp);
    fclose(fp);
    char command[4096];
    snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
    system(command);
    int error = CheckRout(scriptfile);
    switch (error) {
    case 2:  printf("R code %s failed.\n", (const char*)scriptfile);
        break;
    case 1:  printf("R code %s failed.\n", (const char*)scriptfile);
        printf("  Please check %sout for details.\n", (const char*)scriptfile);
        break;
    case 0:
        { String stem(output); stem.Add("_aucplot");
          printf("AUC plots saved in %s%s\n", (const char*)stem, MakePdf(stem)); }
        break;
    default: printf("Unexpected error. Please contact SECRET authors.\n");
    }
    return(error);
}

// ---------------------------------------------------------------------
// Calibration plots of --logistic --readmodel and --coxph --readmodel
// ---------------------------------------------------------------------

// Runs the R script, makes the pdf and sets msg; returns the CheckRout code.
static int RunCalibrationScript(const char *scriptfile, const char *outfile, const char *rpath, const char *what, String &msg)
{
  String command = RCommandLine(rpath, scriptfile);
  system((const char*)command);
  char message[2048];
  message[0] = '\0';
  int error = CheckRout(scriptfile);
  switch (error) {
    case 0: {
      const char *figureExt = MakePdf(outfile);
      snprintf(message, sizeof(message), "%s saved in %s%s\n", what, outfile, figureExt);
      break; }
    case 3:  if (!RMissingPackage(scriptfile, message, sizeof(message)))
               snprintf(message, sizeof(message), "R plots were not made: an R package is not installed. Please check %sout for details.\n", scriptfile);
             break;
    case 2:  snprintf(message, sizeof(message), "R code %s failed.\n", scriptfile); break;
    default: snprintf(message, sizeof(message), "R code %s failed.\n  Please check %sout for details.\n", scriptfile, scriptfile);
  }
  msg.Copy(message);
  return error;
}

// --logistic --readmodel --rplot.  Reads <input>_combinedPGS.txt (ID, outcome 1/0, one column per
// model: the combined score).  The model file has no intercept, so the predicted risk of a model is
//     expit(a0 + score),   a0 fitted in the test sample with the slope fixed at 1,
// i.e. the average predicted risk equals the observed proportion.  Calibration is then the logistic
// recalibration  logit P(affected) = a + b * logit(predicted):  b = calibration slope (ideal 1), a =
// calibration intercept (ideal 0).  Per model: observed vs predicted risk by decile (Wilson 95% CI),
// the fitted calibration curve, and the distribution of the predicted risk in cases and controls;
// then a page with the slopes of all models.  Output <input>_calibrationplot.pdf and <input>_calibration.txt.
int plotCalibrationLogistic(const char *input, const char *traitName, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_calibrationplot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET logistic calibration plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "input <- \"%s\"; trait <- \"%s\"\n", input, traitName);
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  fputs("data <- read.table(paste0(input, \"_combinedPGS.txt\"), sep=\"\\t\", header=TRUE, row.names=1, check.names=FALSE, quote=\"\\\"\", comment.char=\"\", stringsAsFactors=FALSE)\n", fp);
  fputs("y <- as.numeric(data[, 1]); models <- names(data)[-1]\n", fp);
  fputs("expit <- function(x) 1/(1 + exp(-x)); logit <- function(p) log(p/(1 - p))\n", fp);
  fputs("wilson <- function(k, n){ p <- k/n; z <- 1.96; d <- 1 + z^2/n; h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d\n", fp);
  fputs("  cbind((p + z^2/(2*n))/d - h, (p + z^2/(2*n))/d + h) }\n", fp);
  fputs("aucfun <- function(score, y){ r <- rank(score); n1 <- sum(y==1); n0 <- sum(y==0); (sum(r[y==1]) - n1*(n1+1)/2)/(n1*n0) }\n", fp);
  fputs("usable <- function(score, y){ ok <- !is.na(score); sum(ok) >= 20 && sum(y[ok]==1) >= 2 && sum(y[ok]==0) >= 2 && sd(score[ok]) > 0 }\n", fp);
  fputs("note <- function(txt){ par(mfrow=c(1,1)); plot.new(); text(0.5, 0.5, txt) }\n", fp);
  fputs("foot <- \"Predicted risk = expit(a0 + score); the model file has no intercept, so a0 is fitted in this sample (mean predicted risk = observed proportion)\"\n", fp);
  fputs("res <- NULL\n", fp);
  fputs("for(m in models){\n", fp);
  fputs("  s <- as.numeric(data[[m]])\n", fp);
  fputs("  if(!usable(s, y)){ note(paste0(m, \" on \", trait, \" is not available or is constant\\n(no PGS left in the model)\")); next }\n", fp);
  fputs("  ok <- !is.na(s); s <- s[ok]; yy <- y[ok]; n <- length(s)\n", fp);
  fputs("  a0 <- tryCatch(unname(coef(suppressWarnings(glm(yy ~ 1 + offset(s), family=binomial)))[1]), error=function(e) NA)\n", fp);
  fputs("  if(is.na(a0)){ note(paste0(m, \" on \", trait, \": the model could not be fitted\")); next }\n", fp);
  fputs("  lp <- a0 + s; pr <- expit(lp)\n", fp);
  fputs("  fit <- tryCatch(suppressWarnings(glm(yy ~ lp, family=binomial)), error=function(e) NULL)\n", fp);
  fputs("  cf <- if(is.null(fit)) rep(NA_real_, 2) else unname(coef(fit))\n", fp);
  fputs("  ci <- if(is.null(fit)) matrix(NA_real_, 2, 2) else suppressWarnings(confint.default(fit))\n", fp);
  fputs("  auc <- aucfun(s, yy)\n", fp);
  fputs("  res <- rbind(res, data.frame(Model=m, N=n, NCase=sum(yy==1), Slope=cf[2], SlopeLCL=ci[2,1], SlopeUCL=ci[2,2],\n", fp);
  fputs("    Intercept=cf[1], InterceptLCL=ci[1,1], InterceptUCL=ci[1,2], AUC=auc, stringsAsFactors=FALSE))\n", fp);
  fputs("  nb <- if(n >= 500) 10 else 5\n", fp);
  fputs("  brk <- unique(quantile(pr, probs=seq(0, 1, length.out=nb + 1)))\n", fp);
  fputs("  if(length(brk) < 3){ note(paste0(m, \" on \", trait, \" has too few distinct predicted risks\")); next }\n", fp);
  fputs("  bin <- cut(pr, breaks=brk, include.lowest=TRUE)\n", fp);
  fputs("  nk <- as.vector(tapply(yy, bin, length)); kk <- as.vector(tapply(yy, bin, sum)); px <- as.vector(tapply(pr, bin, mean))\n", fp);
  fputs("  keep <- !is.na(nk) & nk > 0; nk <- nk[keep]; kk <- kk[keep]; px <- px[keep]\n", fp);
  fputs("  ob <- kk/nk; wc <- wilson(kk, nk)\n", fp);
  fputs("  par(mfrow=c(1,2), oma=c(2, 0, 0, 0))\n", fp);
  fputs("  lim <- c(0, max(c(px, ob, wc[,2]), na.rm=TRUE)*1.05)\n", fp);
  fputs("  plot(px, ob, xlim=lim, ylim=lim, pch=19, xlab=\"Predicted risk\", ylab=\"Observed proportion affected (95% CI)\", cex.main=0.9,\n", fp);
  fputs("    main=paste0(m, \" on \", trait, \"\\nCalibration (N=\", n, \", cases=\", sum(yy==1), \")\"))\n", fp);
  fputs("  segments(px, wc[,1], px, wc[,2])\n", fp);
  fputs("  abline(0, 1, lty=2, col=\"grey40\")\n", fp);
  fputs("  if(!any(is.na(cf))){\n", fp);
  fputs("    xs <- seq(max(min(pr), 1e-6), min(max(pr), 1 - 1e-6), length.out=100)\n", fp);
  fputs("    lines(xs, expit(cf[1] + cf[2]*logit(xs)), col=\"red\", lwd=2)\n", fp);
  fputs("  }\n", fp);
  fputs("  legend(\"topleft\", bty=\"n\", cex=0.8, legend=c(\n", fp);
  fputs("    sprintf(\"Calibration slope: %.2f (%.2f, %.2f)\", cf[2], ci[2,1], ci[2,2]),\n", fp);
  fputs("    sprintf(\"Calibration intercept: %.2f (%.2f, %.2f)\", cf[1], ci[1,1], ci[1,2]),\n", fp);
  fputs("    sprintf(\"AUC: %.3f\", auc)))\n", fp);
  fputs("  legend(\"bottomright\", bty=\"n\", cex=0.8, legend=c(\"Ideal\", \"Fitted calibration curve\"), lty=c(2, 1), col=c(\"grey40\", \"red\"), lwd=c(1, 2))\n", fp);
  fputs("  dc <- density(pr[yy==0]); da <- density(pr[yy==1])\n", fp);
  fputs("  plot(range(dc$x, da$x), c(0, max(dc$y, da$y)), type=\"n\", xlab=\"Predicted risk\", ylab=\"Density\", cex.main=0.9,\n", fp);
  fputs("    main=paste0(m, \": predicted risk in unaffected and affected\"))\n", fp);
  fputs("  lines(dc, col=\"blue\"); lines(da, col=\"red\")\n", fp);
  fputs("  legend(\"topright\", bty=\"n\", legend=c(paste0(\"Unaffected (n=\", sum(yy==0), \")\"), paste0(\"Affected (n=\", sum(yy==1), \")\")), col=c(\"blue\", \"red\"), lty=1)\n", fp);
  fputs("  mtext(foot, side=1, outer=TRUE, line=0.5, cex=0.6)\n", fp);
  fputs("  par(oma=c(0, 0, 0, 0))\n", fp);
  fputs("}\n", fp);
  fputs("par(mfrow=c(1,1), oma=c(0, 0, 0, 0))\n", fp);
  fputs("if(!is.null(res)){\n", fp);
  fputs("  write.table(res, paste0(input, \"_calibration.txt\"), sep=\"\\t\", quote=FALSE, row.names=FALSE)\n", fp);
  fputs("  r <- res[!is.na(res$Slope), , drop=FALSE]\n", fp);
  fputs("  if(nrow(r) > 0){\n", fp);
  fputs("    k <- nrow(r)\n", fp);
  fputs("    dotchart(rev(r$Slope), labels=rev(r$Model), pch=19, xlim=range(c(1, r$SlopeLCL, r$SlopeUCL), na.rm=TRUE),\n", fp);
  fputs("      xlab=\"Calibration slope (95% CI); 1 = calibrated\", main=paste0(\"Calibration slope of the PGS models for \", trait))\n", fp);
  fputs("    segments(rev(r$SlopeLCL), 1:k, rev(r$SlopeUCL), 1:k)\n", fp);
  fputs("    abline(v=1, lty=3)\n", fp);
  fputs("  }\n", fp);
  fputs("}\n", fp);
  fputs("dev.off()\n", fp);
  fputs("rm(list=ls(all=TRUE))\n", fp);
  fclose(fp);
  return RunCalibrationScript(scriptfile, outfile, rpath, "Calibration plots", msg);
}

// --coxph --readmodel --rplot.  Reads <input>_combinedPGS.txt (ID, time, event, optional STRATUM, one
// column per model: the combined score = log hazard ratio scale).  The model file has no baseline
// hazard, so the baseline is the Breslow estimate with the score as an offset (slope fixed at 1) in
// the test sample, per stratum.  The predicted risk at the horizon t (the 75th percentile of the event
// times) is 1 - exp(-H0(t) exp(score)).  Per model: Kaplan-Meier observed risk at t vs mean predicted
// risk by risk group (Greenwood log-log 95% CI), and the calibration line on the complementary
// log-log scale, cloglog(observed) = a + b cloglog(predicted), fitted across the groups weighted by
// the inverse variance: b = calibration slope (ideal 1), a = calibration intercept (ideal 0).  The
// Cox slope of the score (the usual prognostic-index slope) is shown too.
// Output <input>_calibrationplot.pdf and <input>_calibration.txt.
int plotCalibrationCox(const char *input, const char *timeName, const char *eventName, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_calibrationplot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Cox calibration plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "input <- \"%s\"; trait <- \"%s to %s\"\n", input, timeName, eventName);
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  fputs("data <- read.table(paste0(input, \"_combinedPGS.txt\"), sep=\"\\t\", header=TRUE, row.names=1, check.names=FALSE, quote=\"\\\"\", comment.char=\"\", stringsAsFactors=FALSE)\n", fp);
  fputs("tm <- as.numeric(data[, 1]); ev <- as.numeric(data[, 2])\n", fp);
  fputs("hasS <- ncol(data) >= 3 && names(data)[3] == \"STRATUM\"\n", fp);
  fputs("st <- if(hasS) data[[3]] else rep(1, nrow(data))\n", fp);
  fputs("models <- names(data)[-(1:(if(hasS) 3 else 2))]\n", fp);
  fputs("note <- function(txt){ par(mfrow=c(1,1)); plot.new(); text(0.5, 0.5, txt) }\n", fp);
  fputs("# Kaplan-Meier risk set sums at the event times <= tstar: survival and the Greenwood sum\n", fp);
  fputs("km_at <- function(t, e, tstar){\n", fp);
  fputs("  o <- order(t); t <- t[o]; e <- e[o]\n", fp);
  fputs("  ut <- sort(unique(t[e==1 & t <= tstar])); if(!length(ut)) return(c(1, 0))\n", fp);
  fputs("  nr <- length(t) - findInterval(ut, t, left.open=TRUE)\n", fp);
  fputs("  d <- tabulate(match(t[e==1 & t <= tstar], ut), nbins=length(ut))\n", fp);
  fputs("  c(prod(1 - d/nr), sum(ifelse(nr > d, d/(nr*(nr - d)), 0)))\n", fp);
  fputs("}\n", fp);
  fputs("# Breslow baseline cumulative hazard at tstar with the score as an offset, per stratum\n", fp);
  fputs("base_haz <- function(t, e, s, g, tstar){\n", fp);
  fputs("  out <- c()\n", fp);
  fputs("  for(k in unique(g)){\n", fp);
  fputs("    i <- g == k; tt <- t[i]; ee <- e[i]; w <- exp(s[i])\n", fp);
  fputs("    ut <- sort(unique(tt[ee==1 & tt <= tstar]))\n", fp);
  fputs("    if(!length(ut)){ out[as.character(k)] <- 0; next }\n", fp);
  fputs("    o <- order(tt); ts <- tt[o]; cs <- rev(cumsum(rev(w[o])))\n", fp);
  fputs("    den <- cs[findInterval(ut, ts, left.open=TRUE) + 1]\n", fp);
  fputs("    d <- tabulate(match(tt[ee==1 & tt <= tstar], ut), nbins=length(ut))\n", fp);
  fputs("    out[as.character(k)] <- sum(d/den)\n", fp);
  fputs("  }\n", fp);
  fputs("  out\n", fp);
  fputs("}\n", fp);
  fputs("# Cox slope of one score (Breslow ties, strata): Newton-Raphson\n", fp);
  fputs("cox1 <- function(t, e, s, g){\n", fp);
  fputs("  parts <- lapply(unique(g), function(k){ i <- which(g == k); o <- order(t[i]); list(t=t[i][o], e=e[i][o], s=s[i][o]) })\n", fp);
  fputs("  b <- 0; I <- NA\n", fp);
  fputs("  for(it in 1:30){\n", fp);
  fputs("    U <- 0; I <- 0\n", fp);
  fputs("    for(p in parts){\n", fp);
  fputs("      w <- exp(b*p$s); c0 <- rev(cumsum(rev(w))); c1 <- rev(cumsum(rev(w*p$s))); c2 <- rev(cumsum(rev(w*p$s^2)))\n", fp);
  fputs("      k <- findInterval(p$t, p$t, left.open=TRUE) + 1\n", fp);
  fputs("      j <- which(p$e == 1)\n", fp);
  fputs("      m1 <- c1[k[j]]/c0[k[j]]; m2 <- c2[k[j]]/c0[k[j]]\n", fp);
  fputs("      U <- U + sum(p$s[j] - m1); I <- I + sum(m2 - m1^2)\n", fp);
  fputs("    }\n", fp);
  fputs("    if(!is.finite(I) || I <= 0) return(c(NA, NA))\n", fp);
  fputs("    step <- max(-5, min(5, U/I)); b <- b + step\n", fp);
  fputs("    if(abs(step) < 1e-8) break\n", fp);
  fputs("  }\n", fp);
  fputs("  c(b, 1/sqrt(I))\n", fp);
  fputs("}\n", fp);
  fputs("res <- NULL\n", fp);
  fputs("foot <- \"Predicted risk = 1 - exp(-H0(t) exp(score)); the model file has no baseline hazard, so H0 is estimated in this sample with the score as an offset\"\n", fp);
  fputs("for(m in models){\n", fp);
  fputs("  s <- as.numeric(data[[m]]); ok <- !is.na(s) & !is.na(tm) & !is.na(ev)\n", fp);
  fputs("  if(sum(ok) < 20 || sum(ev[ok]==1) < 5 || !(sd(s[ok]) > 0)){ note(paste0(m, \" on \", trait, \" is not available or is constant\\n(no PGS left in the model)\")); next }\n", fp);
  fputs("  t2 <- tm[ok]; e2 <- ev[ok]; g2 <- st[ok]; s2 <- s[ok] - mean(s[ok]); n <- length(s2)\n", fp);
  fputs("  tstar <- as.numeric(quantile(t2[e2==1], 0.75))\n", fp);
  fputs("  pr <- 1 - exp(-base_haz(t2, e2, s2, g2, tstar)[as.character(g2)] * exp(s2))\n", fp);
  fputs("  cx <- cox1(t2, e2, s2, g2)\n", fp);
  fputs("  G <- if(sum(e2) >= 100 && n >= 500) 10 else 5\n", fp);
  fputs("  brk <- unique(quantile(pr, probs=seq(0, 1, length.out=G + 1)))\n", fp);
  fputs("  if(length(brk) < 3){ note(paste0(m, \" on \", trait, \" has too few distinct predicted risks\")); next }\n", fp);
  fputs("  bin <- cut(pr, breaks=brk, include.lowest=TRUE)\n", fp);
  fputs("  grp <- do.call(rbind, lapply(levels(bin), function(l){\n", fp);
  fputs("    i <- which(bin == l)\n", fp);
  fputs("    if(length(i) < 2 || max(t2[i]) < tstar) return(NULL)\n", fp);
  fputs("    k <- km_at(t2[i], e2[i], tstar); S <- k[1]; g <- k[2]\n", fp);
  fputs("    lo <- hi <- NA_real_\n", fp);
  fputs("    if(S > 0 && S < 1 && g > 0){ se <- sqrt(g)/abs(log(S)); lo <- 1 - S^exp(-1.96*se); hi <- 1 - S^exp(1.96*se) }\n", fp);
  fputs("    data.frame(pred=mean(pr[i]), obs=1 - S, lo=lo, hi=hi, S=S, g=g, n=length(i))\n", fp);
  fputs("  }))\n", fp);
  fputs("  a <- b <- NA_real_; ci <- matrix(NA_real_, 2, 2)\n", fp);
  fputs("  if(!is.null(grp)){\n", fp);
  fputs("    v <- grp$S > 0 & grp$S < 1 & grp$g > 0 & grp$pred > 0\n", fp);
  fputs("    if(sum(v) >= 3){\n", fp);
  fputs("      x <- log(-log(1 - grp$pred[v])); yv <- log(-log(grp$S[v])); w <- log(grp$S[v])^2/grp$g[v]\n", fp);
  fputs("      f <- tryCatch(lm(yv ~ x, weights=w), error=function(e) NULL)\n", fp);
  fputs("      if(!is.null(f) && all(is.finite(coef(f)))){ a <- unname(coef(f)[1]); b <- unname(coef(f)[2])\n", fp);
  fputs("        ci <- tryCatch(confint(f), error=function(e) matrix(NA_real_, 2, 2)) }\n", fp);
  fputs("    }\n", fp);
  fputs("  }\n", fp);
  fputs("  res <- rbind(res, data.frame(Model=m, N=n, Events=sum(e2), Horizon=tstar, Slope=b, SlopeLCL=ci[2,1], SlopeUCL=ci[2,2],\n", fp);
  fputs("    Intercept=a, InterceptLCL=ci[1,1], InterceptUCL=ci[1,2], CoxSlope=cx[1], CoxSlopeSE=cx[2], stringsAsFactors=FALSE))\n", fp);
  fputs("  if(is.null(grp) || nrow(grp) < 2){ note(paste0(m, \" on \", trait, \": too few risk groups followed up to the horizon\")); next }\n", fp);
  fputs("  par(mfrow=c(1,2), oma=c(2, 0, 0, 0))\n", fp);
  fputs("  lim <- c(0, max(c(grp$pred, grp$obs, grp$hi), na.rm=TRUE)*1.05)\n", fp);
  fputs("  plot(grp$pred, grp$obs, xlim=lim, ylim=lim, pch=19, xlab=paste0(\"Predicted risk by time \", signif(tstar, 3)),\n", fp);
  fputs("    ylab=\"Observed risk, Kaplan-Meier (95% CI)\", cex.main=0.9,\n", fp);
  fputs("    main=paste0(m, \" on \", trait, \"\\nCalibration (N=\", n, \", events=\", sum(e2), \")\"))\n", fp);
  fputs("  ok2 <- !is.na(grp$lo); segments(grp$pred[ok2], grp$lo[ok2], grp$pred[ok2], grp$hi[ok2])\n", fp);
  fputs("  abline(0, 1, lty=2, col=\"grey40\")\n", fp);
  fputs("  if(!is.na(a)){\n", fp);
  fputs("    xs <- seq(max(min(grp$pred), 1e-6), min(max(grp$pred), 1 - 1e-6), length.out=100)\n", fp);
  fputs("    lines(xs, 1 - exp(-exp(a)*(-log(1 - xs))^b), col=\"red\", lwd=2)\n", fp);
  fputs("  }\n", fp);
  fputs("  legend(\"topleft\", bty=\"n\", cex=0.8, legend=c(\n", fp);
  fputs("    sprintf(\"Calibration slope: %.2f (%.2f, %.2f)\", b, ci[2,1], ci[2,2]),\n", fp);
  fputs("    sprintf(\"Calibration intercept: %.2f (%.2f, %.2f)\", a, ci[1,1], ci[1,2]),\n", fp);
  fputs("    sprintf(\"Cox slope of the score: %.2f (%.2f, %.2f)\", cx[1], cx[1] - 1.96*cx[2], cx[1] + 1.96*cx[2])))\n", fp);
  fputs("  legend(\"bottomright\", bty=\"n\", cex=0.8, legend=c(\"Ideal\", \"Fitted calibration line\"), lty=c(2, 1), col=c(\"grey40\", \"red\"), lwd=c(1, 2))\n", fp);
  fputs("  hist(pr, breaks=30, col=\"grey80\", border=\"white\", xlab=\"Predicted risk\", main=paste0(m, \": distribution of predicted risk\"), cex.main=0.9)\n", fp);
  fputs("  mtext(foot, side=1, outer=TRUE, line=0.5, cex=0.6)\n", fp);
  fputs("  par(oma=c(0, 0, 0, 0))\n", fp);
  fputs("}\n", fp);
  fputs("par(mfrow=c(1,1), oma=c(0, 0, 0, 0))\n", fp);
  fputs("if(!is.null(res)){\n", fp);
  fputs("  write.table(res, paste0(input, \"_calibration.txt\"), sep=\"\\t\", quote=FALSE, row.names=FALSE)\n", fp);
  fputs("  r <- res[!is.na(res$Slope), , drop=FALSE]\n", fp);
  fputs("  if(nrow(r) > 0){\n", fp);
  fputs("    k <- nrow(r)\n", fp);
  fputs("    dotchart(rev(r$Slope), labels=rev(r$Model), pch=19, xlim=range(c(1, r$SlopeLCL, r$SlopeUCL), na.rm=TRUE),\n", fp);
  fputs("      xlab=\"Calibration slope (95% CI); 1 = calibrated\", main=paste0(\"Calibration slope of the PGS models for \", trait))\n", fp);
  fputs("    segments(rev(r$SlopeLCL), 1:k, rev(r$SlopeUCL), 1:k)\n", fp);
  fputs("    abline(v=1, lty=3)\n", fp);
  fputs("  }\n", fp);
  fputs("}\n", fp);
  fputs("dev.off()\n", fp);
  fputs("rm(list=ls(all=TRUE))\n", fp);
  fclose(fp);
  return RunCalibrationScript(scriptfile, outfile, rpath, "Calibration plots", msg);
}

void plotPopStructure(const char *prefix, int projectFlag, const char *rpath)
{
    String scriptfile = prefix;
    scriptfile.Add("_pcplot.R");
    FILE *fp = fopen(scriptfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
    fprintf(fp, "## %s for secret --pca, by Wei-Min Chen\n", (const char*)scriptfile);
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fprintf(fp, "postscript(\"%s_pcplot.ps\", paper=\"letter\", horizontal=T)\n", (const char*)prefix);
    fprintf(fp, "data <- read.table(file=\"%s_pc.txt\", header=T)\n", (const char*)prefix);
    fprintf(fp, "plot(data$PC1, data$PC2, type=\"p\", xlab=\"PC1\", ylab=\"PC2\", main = \"Population Structure in %s\")\n", (const char*)prefix);
    if (projectFlag) {
        fprintf(fp, "points(data[data[,6]==2,7], data[data[,6]==2,8], col = \"red\")\n");
        fprintf(fp, "legend(\"topright\", c(\"Reference Sample to Generate PCs\", \"Study Sample Projected to Reference PC Space\"),\n");
        fprintf(fp, "col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = 19, cex = 1.2)\n");
    }
    fprintf(fp, "plot(data$PC3, data$PC4, type=\"p\", xlab=\"PC3\", ylab=\"PC4\", main = \"Population Structure in %s\")\n", (const char*)prefix);
    fprintf(fp, "dev.off()\n");
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fclose(fp);
    char command[4096];
    snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
    system(command);
    int error = CheckRout(scriptfile);
    if (error == 2)
        printf("--pca is done but R code %s failed.\n\n", (const char*)scriptfile);
    else if (error) {
        printf("--pca is done but R code %s failed.\n", (const char*)scriptfile);
        printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
    }
    else {
        { String stem(prefix); stem.Add("_pcplot");
          printf("Population structure plot saved in %s%s\n", (const char*)stem, MakePdf(stem)); }
    }
}









/*
void plotPGS(const char *input, const char *output, const char *rpath)
{
    StringArray tokens;
    tokens.Clear();
    tokens.AddTokens(input, '_');
    String outfile = tokens[0];
    outfile.Add("_");
    outfile.Add(output);
    outfile.Add("_grsplot");
    String scriptfile(outfile);
    scriptfile.Add(".R");
    FILE *fp = fopen(scriptfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
    fprintf(fp, "## %s for SECRET PGS plot, by Wei-Min Chen\n", (const char*)scriptfile);
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fprintf(fp, "postscript(\"%s.ps\", paper=\"letter\", horizontal=T)\n", (const char*)outfile);
    fprintf(fp, "data <- read.table(file=\"%s\", header=T)\n", input);
// time-consuming for large sample size
//    fprintf(fp, "pc <- cmdscale(dist(data[, 3:12]), k = 2)\n");
//    fprintf(fp, "plot(pc[, 1], pc[, 2], xlab = \"PC1\", ylab = \"PC2\", main = \"MDS PCs of Top 10 PGS\")\n");
//    fprintf(fp, "points(pc[data[, 2]==2, 1], pc[data[, 2]==2, 2], col = \"red\", pch = 19)\n");
//    fprintf(fp, "rsq <- round(summary(lm(pc[,1]~data[,2]))$r.squared,digit=3)\n");
//    fprintf(fp, "boxplot(pc[,1]~data[,2], xlab=\"%s\", ylab=\"PC1\", main=paste0(\"PC1 of Top 10 PGS (R-squared=\",rsq,\")\"))\n", output);
    fprintf(fp, "IsBinary <- sum(data[,2]==2 | data[,2]==1 | data[,2]==0)==dim(data)[1]\n");
    fprintf(fp, "if(IsBinary){plot(data[,3], data[,4], xlab=names(data)[3], ylab=names(data)[4], main = \"Top 1 vs Top 2 for %s\")\n", output);
    fprintf(fp, "points(data[data[,2]==2,3], data[data[,2]==2,4], col = \"red\", pch=19)\n");
    fprintf(fp, "legend(\"topright\", c(\"Controls\", \"Cases\"),col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = c(1,19), cex = 1.2)\n");
    fprintf(fp, "plot(data[,5], data[,6], xlab=names(data)[5], ylab=names(data)[6], main = \"Top 3 vs Top 4 for %s\")\n", output);
    fprintf(fp, "points(data[data[,2]==2,5], data[data[,2]==2,6], col = \"red\", pch=19)\n");
    fprintf(fp, "legend(\"topright\", c(\"Controls\", \"Cases\"),col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = c(1,19), cex = 1.2)\n");
    fprintf(fp, "for(i in 1:10){boxplot(data[,i+2]~data[,2], xlab=\"%s\", ylab=names(data)[i+2], main=names(data)[i+2])}\n", output);
    fprintf(fp, "}else{for(i in 1:10){\n");
    fprintf(fp, "model <- lm(data[,2] ~ data[,i+2])\n");
    fprintf(fp, "rsq <- round(summary(model)$r.squared, digit=3)\n");
    fprintf(fp, "plot(data[,i+2], data[,2], xlab=names(data)[i+2], ylab=\"%s\", main=paste0(names(data)[i+2], \" on %s (R-squared=\", rsq, \")\"))\n", output, output);
    fprintf(fp, "abline(a=model$coefficient[1], b=model$coefficient[2])}}\n");
    fprintf(fp, "dev.off()\n");
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fclose(fp);
    char command[4096];
    snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
    system(command);
    int error = CheckRout(scriptfile);
    if (error == 2)
        printf("SECRET analysis is done but R code %s failed.\n\n", (const char*)scriptfile);
    else if (error) {
        printf("SECRET analysis is done but R code %s failed.\n", (const char*)scriptfile);
        printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
    }
    else {
        printf("Top PGS plot saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
    }
}

void plotPGS2(const char *input, const char *output, const char *rpath)
{
    StringArray tokens;
    tokens.Clear();
    tokens.AddTokens(input, '_');
    String outfile = tokens[0];
    outfile.Add("_");
    outfile.Add(output);
    outfile.Add("_grsplot");
    String scriptfile(outfile);
    scriptfile.Add(".R");
    FILE *fp = fopen(scriptfile, "wt");
    if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
    fprintf(fp, "## %s for SECRET PGS plot, by Wei-Min Chen\n", (const char*)scriptfile);
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fprintf(fp, "postscript(\"%s.ps\", paper=\"letter\", horizontal=T)\n", (const char*)outfile);
    fprintf(fp, "data <- read.table(file=\"%s\", header=T)\n", input);
    fprintf(fp, "for(i in 1:10){\n");
    fprintf(fp, "plot(data[,i*2+1]-data[,i*2+2], data[,i*2+2], xlab=paste(\"Family Average of\",names(data)[i*2+1]), ylab=\"Adjusted PGS Score\", main = paste(\"Top\", i, \"PGS Scores for %s\"))\n", output);
    fprintf(fp, "points(data[data[,2]==2,i*2+1]-data[data[,2]==2,i*2+2], data[data[,2]==2,i*2+2], col = \"red\", pch=19)\n");
    fprintf(fp, "abline(h=0)\n");
    fprintf(fp, "legend(\"topright\", c(\"Controls\", \"Cases\"),col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = c(1,19), cex = 1.2)}\n");
    fprintf(fp, "dev.off()\n");
    fprintf(fp, "rm(list=ls(all=TRUE))\n");
    fclose(fp);
    char command[4096];
    snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
    system(command);
    int error = CheckRout(scriptfile);
    if (error == 2)
        printf("SECRET analysis is done but R code %s failed.\n\n", (const char*)scriptfile);
    else if (error) {
        printf("SECRET analysis is done but R code %s failed.\n", (const char*)scriptfile);
        printf("  Please check %sout for details.\n\n", (const char*)scriptfile);
    }
    else {
        printf("Top PGS plot is generated in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
    }
}

*/

//    fprintf(fp, "plot(data[,3], data[,4], xlab=names(data)[3], ylab=names(data)[4], main = \"Top 1 vs Top 2 for %s\")\n", traitName);
//    fprintf(fp, "points(data[data[,2]==2,3], data[data[,2]==2,4], col = \"red\", pch=19)\n");
//    fprintf(fp, "legend(\"topright\", c(\"Controls\", \"Cases\"),col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = c(1,19), cex = 1.2)\n");
//    fprintf(fp, "plot(data[,5], data[,6], xlab=names(data)[5], ylab=names(data)[6], main = \"Top 3 vs Top 4 for %s\")\n", traitName);
//    fprintf(fp, "points(data[data[,2]==2,5], data[data[,2]==2,6], col = \"red\", pch=19)\n");
//    fprintf(fp, "legend(\"topright\", c(\"Controls\", \"Cases\"),col=c(\"black\",\"red\"),text.col = c(\"black\", \"red\"), pch = c(1,19), cex = 1.2)\n");
//    fprintf(fp, "for(i in 1:4){boxplot(data[,i+2]~data[,2], xlab=\"%s\", ylab=names(data)[i+2], main=names(data)[i+2])}\n", traitName);

/*
int plotPGSHeatmap(const char *input, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_PGSin50plot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET PGS heapmap plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "data <- read.table(\"%s_PGSin50.txt\", sep=\"\\t\", header = TRUE, row.names=1)\n", input);
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);
  fprintf(fp, "library(pheatmap)\n");
  fprintf(fp, "pheatmap(data, cluster_cols=FALSE, main=\"Heatmap Plot for Top 30 PGS Across 50 Individuals\")\n");
  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  int error = CheckRout(scriptfile);
  switch (error) {
  case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
    break;
  case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
    sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
    break;
  case 0:
    sprintf(message, "Rsq plots saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
    break;
  default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}
*/

    /*
    fprintf(fp, "data <- read.table(file=\"%s.pcars\", header=T, sep=\"\\t\", quote=\"\\\"\")\n", input);
    fprintf(fp, "nPC <- dim(data)[2] - 1\n");
    fprintf(fp, "data <- data[-dim(data)[1], ]\n");
    fprintf(fp, "IsAnno <- colnames(data)[nPC + 1] == \"Annotation\"\n");
    fprintf(fp, "if(IsAnno){nPC <- nPC-1}\n");
    fprintf(fp, "data.heatmap <- as.matrix(data[, 2:(nPC + 1)])\n");
    fprintf(fp, "name <- colnames(data)[2:(nPC+1)]\n");
    fprintf(fp, "if(IsAnno){rownames(data.heatmap) <- paste0(data$RiskScore, \":\", data$Annotation)\n");
    fprintf(fp, "}else{rownames(data.heatmap) <- data$RiskScore}\n");
    fprintf(fp, "heatmap(data.heatmap, scale = \"none\", Colv = NA, Rowv = NA, labCol=name, margins=c(5,15),\n");
    fprintf(fp, "col = colorRampPalette(c(\"red\", \"white\", \"blue\"))(256))\n");
    */
/*
int glmnet(const char *input, const char *traitName, int repeat, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_glmnet");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Rsq plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");

  fprintf(fp, "data <- read.table(\"%s_PGS.txt\", sep=\"\\t\", header = TRUE, row.names=1)\n", input);
  fprintf(fp, "N.r <- dim(data)[2]-2\n");
  fprintf(fp, "x <- as.matrix(data[,c(-1,-N.r-2)])\n");
  fprintf(fp, "y <- data[,1]\n");

  fprintf(fp, "library(glmnet)\n");
  fprintf(fp, "R2 <- R2.CV <- c()\n");
  fprintf(fp, "a.all <- c(0,0.5,1)\n");
  fprintf(fp, "methods <- c(\"Ridge PGS\", \"Elastic Net PGS\", \"Lasso PGS\")\n");
  fprintf(fp, "coefs <- matrix(0, nrow=N.r+1, ncol=0)\n");
  fprintf(fp, "for(i in 1:3){\n");
  fprintf(fp, "  cv.fit <- cv.glmnet(x,y,alpha=a.all[i])\n");
  fprintf(fp, "  R2.cv <- 1 - cv.fit$cvm[cv.fit$lambda==cv.fit$lambda.min] / var(y)\n");
  fprintf(fp, "  R2.CV <- c(R2.CV, R2.cv)\n");
  fprintf(fp, "  model.best <- glmnet(x, y, a=a.all[i], lambda=cv.fit$lambda.min)\n");
  fprintf(fp, "  y.predict <- predict(model.best, newx=x, s=\"lambda.min\")\n");
  fprintf(fp, "  model <- lm(y ~ y.predict)\n");
  fprintf(fp, "  rsq <- round(summary(model)$r.squared, digit=3)\n");
  fprintf(fp, "  plot(y.predict, y, xlab=methods[i], ylab=\"%s\", main=paste0(methods[i], \" on %s (R-squared=\", rsq, \")\"))\n", traitName,traitName);
  fprintf(fp, "  abline(a=model$coefficient[1], b=model$coefficient[2])\n");
  fprintf(fp, "  R2 <- c(R2, rsq)\n");
  fprintf(fp, "  coefs <- cbind(coefs, as.vector(coef(model.best)))\n");
  fprintf(fp, "}\n");
  fprintf(fp, "  colnames(coefs) <- c(\"Ridge\", \"Net\", \"Lasso\"\n");

  if(repeat){
    fprintf(fp, "N <- nrow(x)\n");
    fprintf(fp, "R2.cvs <- matrix(0, nrow=3, ncol=0)\n");
    fprintf(fp, "set.seed(42)\n");
    fprintf(fp, "for(r in 1:%d){\n", repeat);
    fprintf(fp, "  R2.cv <- c()\n");
    fprintf(fp, "  train.idx <- sample(1:N, 0.9*N)\n");
    fprintf(fp, "  x.train <- x[train.idx,]\n");
    fprintf(fp, "  y.train <- y[train.idx]\n");
    fprintf(fp, "  x.test <- x[-train.idx,]\n");
    fprintf(fp, "  y.test <- y[-train.idx]\n");
    fprintf(fp, "  for(a in c(0,0.5,1)){\n");
    fprintf(fp, "    model.cv <- cv.glmnet(x.train, y.train, alpha=a, nfolds=10)\n");
    fprintf(fp, "    corr <- cor(predict(model.cv, newx=x.test, s=\"lambda.min\"), y.test)\n");
    fprintf(fp, "    R2.cv <- c(R2.cv, corr*corr)\n");
    fprintf(fp, "   }\n");
    fprintf(fp, "   R2.cvs <- cbind(R2.cvs, R2.cv)\n");
    fprintf(fp, "}\n");
    fprintf(fp, "R2.table <- round(cbind(R2, apply(R2.cvs, 1, mean)), digit=5)\n");
  }else
    fprintf(fp, "R2.table <- round(cbind(R2, R2.CV), digit=5)\n");
  String outfile2 = "tmp_"; outfile2.Add(input);
  outfile2.Add("_net_rsq.txt");
  fprintf(fp, "write.table(R2.table, file=\"%s\", row.names=FALSE, col.names=FALSE, sep=\"\\t\")\n", (const char*)outfile2);
  outfile2 = "tmp_"; outfile2.Add(input);
  outfile2.Add("_net_coef.txt");
  fprintf(fp, "write.table(coefs, file=\"%s\")\n", (const char*)outfile2);

//  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  message[0]='\0';
  int error = CheckRout(scriptfile);
  switch (error) {
    case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      break;
    case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
      break;
    case 0:
//      sprintf(command, "ps2pdf -sPAPERSIZE=letter %s.ps %s.pdf", (const char*)outfile, (const char*)outfile);
//      system(command);
//      sprintf(message, "Rsq plots saved in %s.pdf\n", (const char*)outfile);
      break;
    default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}
*/
/*
int plotRSQ(const char *input, const char *traitName, const char *rpath, String &msg)
{
  String outfile = input;
  outfile.Add("_rsqplot");
  String scriptfile(outfile);
  scriptfile.Add(".R");
  FILE *fp = fopen(scriptfile, "wt");
  if (fp == NULL) error("Cannot open %s to write.", (const char*)scriptfile);
  fprintf(fp, "## %s for SECRET Rsq plot, by Wei-Min Chen\n", (const char*)scriptfile);
  fprintf(fp, "rm(list=ls(all=TRUE))\n");
  fprintf(fp, "postscript(\"%s.ps\", paper = \"letter\", horizontal = T)\n", (const char*)outfile);

  fprintf(fp, "library(pheatmap)\n");
  fprintf(fp, "data <- read.table(\"%s_PGSin50.txt\", sep=\"\\t\", header = TRUE, row.names=1)\n", input);
  fprintf(fp, "pheatmap(data, cluster_cols=FALSE, main=\"Heatmap Plot for Top 30 Polygenic Scores Across 50 %s-Sorted Individuals\")\n", traitName);

  fprintf(fp, "data <- read.table(\"%s.linear\", sep=\"\\t\", header = TRUE, nrows=30, quote=\"\\\"\")\n", input);
  fprintf(fp, "plot(data$Rsq[1:30], 1:30, pch=19, col=c(\"blue\", \"red\")[(data$Z<0)+1],\n");
  fprintf(fp, "xlab = \"R-Squared of Association between PGS and %s\", ylab = \"PGS Rank (by R-Squared)\", main = \"Top 30 PGS Scores That Predict %s\")\n", traitName, traitName);
  fprintf(fp, "left <- min((1:30)[data$Rsq[1:30] < (min(data$Rsq[1:30]) + max(data$Rsq[1:30]))/2])\n");
  fprintf(fp, "if(names(data)[dim(data)[2]]==\"Annotation\"){\n");
  fprintf(fp, "text(data$Rsq[1:(left-1)], 1:(left-1), labels = paste0(data$PGS[1:(left-1)], \":\", data$Annotation[1:(left-1)]), pos = 2)\n");
  fprintf(fp, "text(data$Rsq[left:30], left:30, labels = paste0(data$PGS[left:30], \":\", data$Annotation[left:30]), pos = 4)\n");
  fprintf(fp, "}else{text(data$Rsq[1:(left-1)], 1:(left-1), labels = data$PGS[1:(left-1)], pos = 2)\n");
  fprintf(fp, "text(data$Rsq[left:30], left:30, labels = data$PGS[left:30], pos = 4)}\n");

  fprintf(fp, "data <- read.table(file=\"%s_wtsum.score\", header=T)\n", input);
  fprintf(fp, "data <- data[-1,]\n");
  fprintf(fp, "N.rs <- dim(data)[2] - 2\n");
  fprintf(fp, "IsBinary <- sum(data[,2]==2 | data[,2]==1 | data[,2]==0)==dim(data)[1]\n");
  fprintf(fp, "if(IsBinary){\n");
  fprintf(fp, "for(i in 1:N.rs){Q <- quantile(data[,i+2], probs=seq(0,1,0.2))\n");
  fprintf(fp, "M <- matrix(0,5,2)\n");
  fprintf(fp, "for(k in 1:5){IsQ <- data[,i+2]>Q[k] & data[,i+2]<=Q[k+1]\n");
  fprintf(fp, "M[k,1] <- sum(data[IsQ,2]==1)\n");
  fprintf(fp, "M[k,2] <- sum(data[IsQ,2]==2)}\n");
  fprintf(fp, "OR <- c(1)\n");
  fprintf(fp, "SE <- c(0)\n");
  fprintf(fp, "for(k in 2:5){OR <- c(OR, (M[k,2]*M[1,1])/(M[k,1]*M[1,2]))\n");
  fprintf(fp, "SE <- c(SE, sqrt(1/M[k,2]+1/M[1,1]+1/M[k,1]+1/M[1,2]))}\n");
  fprintf(fp, "if(sum(OR==Inf)==0&&sum(SE==Inf)==0){plot(1:5, OR, type=\"l\", xlab=paste(names(data)[i+2],\"Quintiles\"),ylab=\"OR %s (95%% CI)\", \n", traitName);
  fprintf(fp, "ylim=c(0, max(OR*exp(1.96*SE))),main=paste(names(data)[i+2],\"on %s\"))\n", traitName);
  fprintf(fp, "for(k in 2:5){lines(c(k,k), OR[k] * exp(1.96*SE[k]*c(-1,1)))}}}\n");
  fprintf(fp, "}else{for(i in 1:N.rs){\n");
  fprintf(fp, "model <- lm(data[,2] ~ data[,i+2])\n");
  fprintf(fp, "rsq <- round(summary(model)$r.squared, digit=3)\n");
  fprintf(fp, "plot(data[,i+2], data[,2], xlab=names(data)[i+2], ylab=\"%s\", main=paste0(names(data)[i+2], \" on %s (R-squared=\", rsq, \")\"))\n", traitName, traitName);
  fprintf(fp, "abline(a=model$coefficient[1], b=model$coefficient[2])}}\n");

  fprintf(fp, "dev.off()\n");
  fprintf(fp, "rm(list=ls())\n");
  fprintf(fp, "rm(list=ls(all.names=TRUE))\n");
  fprintf(fp, "q(save=\"no\")\n");
  fclose(fp);
  char command[4096];
  snprintf(command, sizeof(command), "%s", (const char*)RCommandLine(rpath, scriptfile));
  system(command);
  char message[256];
  int error = CheckRout(scriptfile);
  switch (error) {
    case 2:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      break;
    case 1:  sprintf(message, "R code %s failed.\n", (const char*)scriptfile);
      sprintf(message, "  Please check %sout for details.\n", (const char*)scriptfile);
      break;
    case 0:
      sprintf(message, "Rsq plots saved in %s%s\n", (const char*)outfile, MakePdf((const char*)outfile));
      break;
    default: sprintf(message, "Unexpected error. Please contact SECRET authors.\n");
  }
  msg.Copy(message);
  return(error);
}
*/

