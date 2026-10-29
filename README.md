# How much data do you need to learn a Bayesian network?

A simulation-based sample-size / power analysis for Bayesian network
structure learning with [`bnlearn`](https://www.bnlearn.com/), built around a
worked example: predicting response to oral ketamine from clinical,
neuropsychological, and actigraphy variables with n = 60.

The repo has two parts:

1. **Structure recovery (scripts `R/01`-`R/05`)**: can a BN (`hc()`) recover the arcs of a
   known 16- or 19-node network (15 candidate predictors in scenario 1; plus 3 early-change
   variables in scenario 2) at a given n? (20-predictor variant in `R/01_define_networks.R`)
2. **Prediction performance (scripts `R/06`-`R/09`)**: if we train on n participants, how
   well do we predict response in *new* participants? Reports out-of-sample AUC, calibration 
   and decision metrics for the 11 pre-specified predictors, compares BN pipeline to lasso,
   to a clinical-only model and to an oracle, and last estimates the chance that an 
   independent n = 60 validation set confirms useful discrimination.

**[Read the tutorial →](https://MatthieuVignes.github.io/ketamine-bn-tutorial/)**

## What's here

| Path | What it is |
|---|---|
| `index.qmd` | The Quarto tutorial - narrative + code + figures |
| `R/` | The underlying R scripts, sourced by the tutorial |
| `data/` | Simulation outputs (CSV/RDS) - pre-computed, checked in |
| `figures/` | Diagrams and plots - pre-computed, checked in |
| `_freeze/` | Quarto's cached execution results (see below) |
| `install.R` | Installs required packages |
| `run_all.R` | Runs the full pipeline standalone (without Quarto) |

## Quick start

```r
# From the repo root in R:
source("install.R")
source("run_all.R") # regenerates everything in data/ and figures/
```

To render the tutorial itself:

```sh
quarto render index.qmd
```

## Publishing to GitHub Pages

This repo is already set up for the simplest Quarto → GitHub Pages workflow.

## Requirements to actually run the R code

- R ≥ 4.4 (bnlearn's current CRAN requirement)
- Packages: `bnlearn`, `igraph`, `ggplot2`, `glmnet`, `MASS` (see `install.R`)
- [Quarto](https://quarto.org/docs/get-started/) ≥ 1.4, if rendering the tutorial

## License

Licensed under [CeCILL-B](https://cecill.info/licences/Licence_CeCILL-B_V1-en.html) a permissive, BSD/MIT-equivalent license requiring attribution, chosen so the code is freely reusable for your own project.

