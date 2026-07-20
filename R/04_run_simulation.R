suppressPackageStartupMessages({
  library(bnlearn)
  library(parallel)
})
source("R/03_metrics.R")

dags   <- readRDS("data/true_dags.rds")
fitted <- readRDS("data/fitted_true_bns.rds")

n_grid   <- c(20, 30, 40, 50, 60, 75, 100, 150, 200, 300, 500) # tested sample sizes
n_reps   <- 200 # Monte Carlo repeats per (network, n) cell
n_cores  <- max(1, detectCores() - 1)

networks <- list(A_baseline_only = list(dag = dags$dagA, fitted = fitted$fittedA),
  B_with_time     = list(dag = dags$dagB, fitted = fitted$fittedB))

run_cell <- function(net_label, net, n, rep_id) {
  set.seed(net_label |> utils::URLencode() |> nchar() |> (\(x) x * 100000 + n * 1000 + rep_id)())
  d <- rbn(net$fitted, n = n)
  learned <- tryCatch(hc(d, score = "bic-g"), error = function(e) NULL)
  if (is.null(learned)) return(NULL)
  res <- evaluate_reconstruction(net$dag, learned)
  res$network <- net_label
  res$n <- n
  res$rep <- rep_id
  res
}

jobs <- expand.grid(net_label = names(networks),
  n = n_grid,
  rep_id = seq_len(n_reps),
  stringsAsFactors = FALSE)
cat(sprintf("Total simulation cells to run: %d (networks x n-values x reps)\n", nrow(jobs)))

t0 <- Sys.time()
results_list <- mclapply(seq_len(nrow(jobs)), function(i) {
  j <- jobs[i, ]
  run_cell(j$net_label, networks[[j$net_label]], j$n, j$rep_id)
}, mc.cores = n_cores)
t1 <- Sys.time()
cat("Simulation wall time:", round(as.numeric(t1 - t0, units = "secs"), 1), "s\n")

results_list <- results_list[!vapply(results_list, is.null, logical(1))]
results <- do.call(rbind, results_list)

saveRDS(results, "data/simulation_results.rds")
write.csv(results, "data/simulation_results.csv", row.names = FALSE)
cat("Saved", nrow(results), "rows.\n")

##################################################
## Edge-frequency data for ROC / precision-recall curves (see
## 05_summarize_and_plot.R). At a handful of focal sample sizes, record
## how often EACH possible node pair was learned as an edge across the
## n_reps independent simulated datasets. That per-pair frequency
## (0-1) is a natural continuous "score" for whether an edge is real,
## which is what a threshold-sweep ROC/PR curve needs. A single hc()
## fit only gives a present/absent decision, not a score, so this
## reuses the same repeated-simulation design already run above rather
## than requiring a separate bootstrap procedure.
##################################################
focal_n <- c(20, 60, 100, 200)
focal_n <- focal_n[focal_n %in% n_grid]

compute_edge_frequencies <- function(net_label, net, n, n_reps) {
  universe <- skeleton_universe(net$dag)
  hits <- integer(nrow(universe))
  names(hits) <- universe$pair
  net_id <- match(net_label, names(networks))

  for (rep_id in seq_len(n_reps)) {
    set.seed(net_id * 1000000 + n * 1000 + rep_id)
    d <- rbn(net$fitted, n = n)
    learned <- tryCatch(hc(d, score = "bic-g"), error = function(e) NULL)
    if (is.null(learned)) next
    present <- skeleton_pairs(learned)
    hits[present] <- hits[present] + 1L
  }
  universe$frequency <- hits[universe$pair] / n_reps
  universe$network <- net_label
  universe$n <- n
  universe
}

cat("Computing edge frequencies for ROC/PR curves at n =", paste(focal_n, collapse = ", "), "...\n")
edge_freq_list <- list()
for (net_label in names(networks)) {
  for (n in focal_n) {
    edge_freq_list[[paste(net_label, n)]] <- compute_edge_frequencies(
      net_label, networks[[net_label]], n, n_reps
    )
  }
}
edge_freq <- do.call(rbind, edge_freq_list)
write.csv(edge_freq, "data/edge_frequency.csv", row.names = FALSE)

