suppressPackageStartupMessages({
  library(bnlearn)
  library(parallel)
})

source("R/06_prediction_generators.R")
source("R/07_prediction_metrics.R")

N_REPS      <- as.integer(Sys.getenv("N_REPS", 200)) # replicates per cell
N_REPS_BOOT <- as.integer(Sys.getenv("N_REPS_BOOT", 60)) # BN_hc_boot is ~20x slower: fewer reps
MC_CORES    <- as.integer(Sys.getenv("MC_CORES", max(1, detectCores() - 1)))
N_VAL       <- 60 # size of the independent validation set (K-BAT)
N_BIG       <- 3000 # size of the large "truth" test set to derive perf metrics
PT          <- 0.40 # decision threshold: treat if P(response) >= 0.40

R2_GRID     <- c(0.10, 0.20, 0.30, 0.45) # signal strength of Outcome
NTRAIN_GRID <- c(60, 120, 200) # new study / new + K-BAT pooled / ample (yeah?)

# If we have it, we can load the baseline-predictor correlation matrix from K-BAT 
# (no outcome information), so that simulated predictors have realistic collinearity:
#   pcor <- cor(kbat[, PRIMARY_PREDICTORS], use = "pairwise.complete.obs")
pcor <- NULL
scenarios <- scenario_library(predictor_cor = pcor)

jobs <- expand.grid(scenario = names(scenarios), r2 = R2_GRID, n_train = NTRAIN_GRID,
                    rep = seq_len(N_REPS), stringsAsFactors = FALSE)
cat(sprintf("Cells: %d  (scenarios x r2 x n_train x reps).  Cores: %d\n", nrow(jobs), MC_CORES))

run_job <- function(i) {
  j  <- jobs[i, ]
  sc <- scenarios[[j$scenario]]
  set.seed(1e6 * match(j$scenario, names(scenarios)) + 1e4 * which(R2_GRID == j$r2) +
             10 * j$n_train + j$rep)
  tr  <- gen_data(sc, j$n_train, j$r2)
  val <- gen_data(sc, N_VAL, j$r2)
  big <- gen_data(sc, N_BIG, j$r2)
  n1 <- sum(val$Responder == 1); n0 <- sum(val$Responder == 0)

  learners <- LEARNERS
  if (j$rep > N_REPS_BOOT) learners[["BN_hc_boot"]] <- NULL
  # the early marker is only a legitimate predictor in the D_early scenario
  out <- lapply(names(learners), function(nm) {
    f <- tryCatch(learners[[nm]](tr), error = function(e) NULL)
    if (is.null(f)) return(NULL)
    pv <- f(val); pb <- f(big)
    p_val <- to_prob(pv, sc$thr); p_big <- to_prob(pb, sc$thr)
    a_val <- auc_rank(p_val, val$Responder)
    ci <- auc_hm_ci(a_val, n1, n0)
    dec <- decision_summary(p_big, big$Responder, PT)
    data.frame(scenario = j$scenario, r2 = j$r2, n_train = j$n_train, rep = j$rep,
               learner = nm,
               auc_val = a_val, auc_val_lb = ci[1],
               auc_true = auc_rank(p_big, big$Responder),
               cal_slope = calibration_slope(p_big, big$Responder),
               brier = mean((p_big - big$Responder)^2),
               nb_gain = dec[["nb_gain"]], spared = dec[["spared"]], denied = dec[["denied"]],
               n_blanket = if (!is.null(pb$blanket)) length(pb$blanket) else NA_integer_,
               shrink = if (!is.null(pb$shrink)) pb$shrink else NA_real_)
  })
  do.call(rbind, out)
}

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(jobs)), run_job, mc.cores = MC_CORES)
cat("Wall time:", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])

dir.create("data", showWarnings = FALSE)
write.csv(res, "data/prediction_results.csv", row.names = FALSE)
cat("Saved", nrow(res), "rows to data/prediction_results.csv\n")
