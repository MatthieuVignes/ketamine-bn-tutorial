# Run from the repo root
#  setwd("~/ketamine-bn-tutorial-main")
#  source("run_all.R")

source("install.R") # installs bnlearn/ggplot2/igraph if missing

## Estimating network structure recovery ##
# Version 1 - 20+1(response) = 21 nodes and 20+3(repeated measures from the 20)+1(response) = 24 nodes
# source("R/01_define_networks.R") # builds the two true DAGs + diagrams
# Version 2 - 15+1(response) = 16 nodes and 15+3(repeated measures from the 20)+1(response) = 19 nodes
source("R/01_define_networks(smaller).R") # builds the two true DAGs + diagrams

source("R/02_fit_parameters.R") # attaches linear-Gaussian parameters, sanity-checks R^2
source("R/03_metrics.R") # computes performance of prediction vs true networks
# This can take a few minutes
source("R/04_run_simulation.R") # Monte Carlo (200 reps) simulation/hc() reconstruction power curve for n=20..500
source("R/05_summarise_and_plot.R") # tables (data/) + plots (figures/)

## Estimating predictive performance with varied noise in network relationships ##
Sys.setenv(N_REPS = 100, N_REPS_BOOT = 40) # quick test:  Sys.setenv(N_REPS = 5, N_REPS_BOOT = 2)
source("R/08_run_prediction_simulation.R") # trains on n=60/120/200, tests on n=60 validation + large test set
source("R/09_summarise_prediction.R")      # table2 + plots 6-9
# Full run ~20-40 min on 2 cores