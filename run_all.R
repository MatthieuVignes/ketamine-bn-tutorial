##################################################
## Run the entire pipeline end-to-end, in order.
## Run from the repo root, e.g.:
##   setwd("~/ketamine-bn-power-analysis")
##   source("run_all.R")
##################################################

source("install.R")                    # installs bnlearn/ggplot2/igraph if missing
source("R/01_define_networks.R")       # builds the two true DAGs + diagrams
source("R/02_fit_parameters.R")        # attaches linear-Gaussian parameters, sanity-checks R^2
source("R/04_run_simulation.R")        # PRIMARY: single hc() power curve, n=20..500, 200 MC reps
source("R/05_summarize_and_plot.R")    # tables (data/) + plots (figures/)
