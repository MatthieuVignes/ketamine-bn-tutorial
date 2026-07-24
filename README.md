# How much data do you need to learn a Bayesian network?

A simulation-based sample-size / power analysis for Bayesian network
structure learning with [`bnlearn`](https://www.bnlearn.com/), built around a
worked example: predicting treatment response from 15 (or 20 in a different version) clinical,
neuropsychological, and actigraphy variables with n = 60.

**[Read the tutorial →](https://MatthieuVignes.github.io/ketamine-bn-tutorial/index.html)**

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
source("run_all.R")   # regenerates everything in data/ and figures/
```

To render the tutorial itself:

```sh
quarto render index.qmd
```

## Publishing to GitHub Pages

This repo is already set up for the simplest Quarto → GitHub Pages workflow.

**About `_freeze/`:** heavy chunks in `index.qmd` (the main simulation,
~1 min) use Quarto's [freeze](https://quarto.org/docs/projects/code-execution.html#freeze)
feature, which caches computed output in `_freeze/`. That folder is
committed to the repo on purpose — it means GitHub Pages (or anyone cloning
the repo) gets the tutorial with all its output already rendered, without
needing R, `bnlearn`, or several minutes of compute. If you change the code
in a cached chunk, re-render locally (`quarto render index.qmd`) to
refresh `_freeze/` before pushing.

## Requirements to actually run the R code

- R ≥ 4.4 (bnlearn's current CRAN requirement)
- Packages: `bnlearn`, `igraph`, `ggplot2` (see `install.R`)
- [Quarto](https://quarto.org/docs/get-started/) ≥ 1.4, if rendering the tutorial

## License

Licensed under [CeCILL-B](https://cecill.info/licences/Licence_CeCILL-B_V1-en.html) a permissive, BSD/MIT-equivalent license requiring attribution, chosen so the code is freely reusable for your own project.

