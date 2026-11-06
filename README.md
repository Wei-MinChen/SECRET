# SECRET

**S**creen and **E**valuate **C**atalogued **R**isk scores to **E**nhance **T**rait predictions.

Version 0.9.3 — (c) 2023-2026 Wei-Min Chen. See the license notice at the top of each source file.
More about SECRET: <https://kingrelatedness.com/SECRET/>

SECRET screens many polygenic scores (PGS), for example scores from the
[PGS Catalog](https://www.pgscatalog.org/), against the traits in your data set. For every trait it
ranks the scores and builds combined predictions from the top ones. The scores themselves can be
calculated from PLINK genotypes with [secretscores](https://github.com/Wei-MinChen/SECRET-scores), which
writes the score file that SECRET reads.

## Build

You need a C++ compiler and zlib (`zlib1g-dev` on Debian/Ubuntu). OpenMP is optional but gives
multithreading.

```bash
make
```

This creates the `secret` executable. Without OpenMP, use `make CXXFLAGS=-O2`. The code avoids C++11
features, so older compilers such as C++ Builder 6 can build it too.

For the figures (`--rplot`) you also need [R](https://www.r-project.org/) with the package `pheatmap`
(`install.packages("pheatmap")`). The figures are written as `.ps` files and turned into PDF by
`ps2pdf`, so install Ghostscript to get PDF files.

## Quick start

```bash
./secret --score scores.txt --pheno pheno.txt --auc
./secret --score scores.txt --pheno pheno.txt --logistic --trait T1 --covariate age,sex --rplot
./secret --score scores.txt --pheno pheno.txt --linear --trait height
./secret --score scores.txt --pheno pheno.txt --coxph --surv time,event
```

Add `--diction` to name the scores in the output and `--prefix` to choose the output file names
(default `secret`).

## Input files

| Option | File |
|---|---|
| `--score` | The polygenic scores: a table with a header line, the sample ID (IID) in the first column and one column per score. This is the `<prefix>_scores.txt` file of secretscores. If IIDs repeat across families, the IDs are `FID_IID`. |
| `--pheno` | The traits: a header line, then one row per sample. The first column is the IID, or the header starts with `FID IID` when family IDs are given. Binary traits are coded 1/2 (1 = unaffected, 2 = affected); a trait coded 0/1 is recoded automatically. Tab, space or comma delimited. |
| `--covarfile` | Optional: covariates in a separate file. Without `--covariate`, every column of it is a covariate. |
| `--diction` | Optional: a dictionary of the scores. Either two tab-delimited columns (score name, description), or a comma-delimited file with the columns `score` and `trait`. The `<prefix>_scores.coverage` file of secretscores can be used as it is; its SNPs used and coverage then appear in the output too. The description goes to the `Annotation` column. |
| `--readmodel` | Optional: a model file (`.mod`) written by an earlier scan, to evaluate it in another data set. |

For `--tdt`, the phenotype file needs `FID IID FA MO` as its first four columns (family, individual,
father and mother), followed by the traits.

## Scans

| Option | What it does |
|---|---|
| `--auc` | AUC of every score for every binary trait. Writes one table per trait, `<prefix>_<trait>.auc`, with the best score first, and prints the cross-validated AUC of the best score and how consistently the same score is chosen (`--fold`, `--repeat`). With `--strat`, the AUCs are also computed in each stratum. |
| `--logistic` | Logistic regression of each binary trait on each score: OR, Z and P, written to `<prefix>_<trait>.logistic`. Reports the AUC of the best score and of combinations of the top `--ntop` scores: a weighted sum (after pruning correlated scores at R² > `--pruner2`), ridge, elastic net and lasso. `--cv` evaluates them by cross-validation. |
| `--linear` | The same for quantitative traits: `<prefix>_<trait>.linear` with the R² of each score. |
| `--coxph` | Cox regression of each score for a time-to-event outcome given by `--surv time,event`: hazard ratio, Z and P, written to `<prefix>_<time>2<event>.coxph`, and penalized models of the top scores. |
| `--tdt` | Transmission test of the scores in parent–offspring trios: `<prefix>_<trait>.tdt`. |

## Options

| Option | Meaning |
|---|---|
| `--trait <list>` | Only analyze these traits (comma-separated). Default: all |
| `--covariate <list>` | Covariates (comma-separated columns of the phenotype or covariate file) |
| `--surv <time,event>` | Survival time and event columns for `--coxph` |
| `--strat <trait>` | Stratify `--auc` and `--tdt` by a binary trait; for `--coxph`, the stratum variable(s) |
| `--cv`, `--fold <n>`, `--repeat <n>` | Cross-validation, its folds and repeats |
| `--ntop <n>` | Number of top scores in the combined models (default 30) |
| `--pruner2 <r2>` | R² above which a correlated score is dropped from the weighted sum (default 0.2) |
| `--cpus <n>` | Threads (default: half of the available cores) |
| `--rplot` | Make the figures with R |
| `--rpath <path>` | The R executable, if R is not in the PATH |
| `--prefix <name>` | Output file prefix (default `secret`) |

## Output

Every scan writes its table to `<prefix>_<trait>.<scan>` (see the table above); in the tables of `--auc`,
`--logistic`, `--linear` and `--coxph` the strongest scores come first. The correlations among the top scores
are written to `*_PGScorr.txt`, and the model files that `--readmodel` reads to `*.mod`. With `--rplot`, each
trait gets its own figures, for example `<prefix>_<trait>_aucplot.pdf` for the top 30 scores of `--auc`, with a
heatmap of their correlations.
