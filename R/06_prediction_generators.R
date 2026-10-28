suppressPackageStartupMessages({
  library(bnlearn)
  library(MASS)
})

PRIMARY_PREDICTORS <- c("BaselineSeverity", "Age", "Gender", "NPriorTrials",
                        "PsychomotorSpeed", "ExecFunction",
                        "RA", "IV", "IS", "SE", "M10")
CLINICAL_PREDICTORS <- c("BaselineSeverity", "Age", "Gender", "NPriorTrials")

# Assumed true direct effects on Outcome (relative sizes/signs only; the overall 
# scale is set by r2). 6 of the 11 predictors truly matter, 5 are null (arbitrary)
# Can be changed (or future: run a sensitivity analysis) if we have lit-based prior
TRUE_EFFECTS <- c(BaselineSeverity = -0.30, NPriorTrials = -0.25,
                  PsychomotorSpeed = 0.35, RA = 0.30, M10 = 0.20, SE = 0.15)

# 14 extra, correlated, NULL candidates for the "wide" scenario
# (other actigraphy metrics, individual cognitive subtests, etc.)
EXTRA_PREDICTORS <- c("L5", "L5Time", "M10Time", "SleepOnsetLatency", "TotalSleepTime",
                      "WASO", "TrailsA", "TrailsB", "DigitSymbol", "DigitSpan",
                      "Fluency", "Chronicity", "Anxiety", "ChildhoodTrauma")

PREDICTOR_DAG_STRING <- paste0(
  "[Age][Gender][NPriorTrials|Age]",
  "[BaselineSeverity|NPriorTrials]",
  "[PsychomotorSpeed|Age:BaselineSeverity]",
  "[ExecFunction|Age:PsychomotorSpeed]",
  "[RA|BaselineSeverity:Age]",
  "[M10|RA:PsychomotorSpeed]",
  "[IS|RA][IV|RA]",
  "[SE|BaselineSeverity:RA]"
)

EXTRA_DAG_STRING <- paste0(
  "[L5|RA][L5Time|RA:Age][M10Time|M10][SleepOnsetLatency|SE]",
  "[TotalSleepTime|SE][WASO|SE:IV]",
  "[TrailsA|PsychomotorSpeed][TrailsB|ExecFunction:PsychomotorSpeed]",
  "[DigitSymbol|PsychomotorSpeed][DigitSpan|ExecFunction][Fluency|ExecFunction]",
  "[Chronicity|Age:NPriorTrials][Anxiety|BaselineSeverity][ChildhoodTrauma|Anxiety]"
)

# Arc coefficients among predictors: |beta| in [0.20, 0.45], random sign, fixed 
# seed. Residual sd chosen so that each node's marginal variance about 1 
# (exactly standardised afterwards in make_predictor_reference())
assign_predictor_coefs <- function(dag, seed) {
  set.seed(seed)
  dist <- list()
  for (nd in nodes(dag)) {
    pars <- parents(dag, nd)
    if (length(pars) == 0) {
      dist[[nd]] <- list(coef = c("(Intercept)" = 0), sd = 1)
    } else {
      b <- runif(length(pars), 0.20, 0.45) * sample(c(-1, 1), length(pars), replace = TRUE)
      names(b) <- pars
      dist[[nd]] <- list(coef = c("(Intercept)" = 0, b),
                         sd = sqrt(max(0.3, 1 - sum(b^2))))
    }
  }
  dist
}

# Function to simulate 1 scenario
#   wide: TRUE adds the 14 NULL extra candidates
#   early_rho: NULL (no early marker) or correlation of the early marker 
#             `EarlyChange` with the final Outcome
#   predictor_cor: NULL, or a correlation matrix (dimnames = predictor names; 
#                  must cover the 11 primary predictors) from real, outcome-blind 
#                  baseline data
make_scenario <- function(wide = FALSE, early_rho = NULL, predictor_cor = NULL,
                          response_rate = 0.60, seed = 2026) {
  predictors <- if (wide) c(PRIMARY_PREDICTORS, EXTRA_PREDICTORS) else PRIMARY_PREDICTORS
  dag <- model2network(paste0(PREDICTOR_DAG_STRING, if (wide) EXTRA_DAG_STRING else ""))
  fit <- custom.fit(dag, dist = assign_predictor_coefs(dag, seed))

  if (!is.null(predictor_cor)) {
    stopifnot(all(PRIMARY_PREDICTORS %in% colnames(predictor_cor)))
    if (wide) stop("predictor_cor currently supports the 11 primary predictors only (wide = FALSE).")
    predictor_cor <- predictor_cor[PRIMARY_PREDICTORS, PRIMARY_PREDICTORS]
  }

  # reference sample to standardise predictors and get Var(X beta0)
  set.seed(seed + 1)
  ref <- if (is.null(predictor_cor)) rbn(fit, 50000) else
    as.data.frame(mvrnorm(50000, rep(0, ncol(predictor_cor)), predictor_cor))
  ref <- ref[, predictors, drop = FALSE]
  mu <- colMeans(ref); sdv <- apply(ref, 2, sd)
  refz <- scale(ref, center = mu, scale = sdv)

  beta0 <- setNames(rep(0, length(predictors)), predictors)
  beta0[names(TRUE_EFFECTS)] <- TRUE_EFFECTS
  z_sd <- sd(as.vector(refz %*% beta0))

  list(fit = fit, dag = dag, predictors = predictors, mu = mu, sdv = sdv,
       beta0 = beta0, z_sd = z_sd, early_rho = early_rho,
       predictor_cor = predictor_cor,
       thr = qnorm(1 - response_rate), # Outcome > thr <=> responder
       response_rate = response_rate)
}

# Simulate n participants from a scenario at signal strength r2
#   Outcome: continuous, standardised improvement (var 1)
#   Responder: 0/1, fixed population threshold (60% responders)
gen_data <- function(sc, n, r2) {
  X <- if (is.null(sc$predictor_cor)) {
    rbn(sc$fit, n)[, sc$predictors, drop = FALSE]
  } else {
    d <- as.data.frame(mvrnorm(n, rep(0, ncol(sc$predictor_cor)), sc$predictor_cor))
    d[, sc$predictors, drop = FALSE]
  }
  X <- as.data.frame(scale(X, center = sc$mu, scale = sc$sdv))
  lp <- as.vector(as.matrix(X) %*% sc$beta0) / sc$z_sd      # unit-variance linear predictor
  Outcome <- sqrt(r2) * lp + sqrt(1 - r2) * rnorm(n)
  if (!is.null(sc$early_rho)) {
    X$EarlyChange <- sc$early_rho * Outcome + sqrt(1 - sc$early_rho^2) * rnorm(n)
  }
  X$Outcome   <- Outcome
  X$Responder <- as.integer(Outcome > sc$thr)
  X
}

# The 3 scenarios used in the application's power section.
#   "C_primary": the 11 pre-specified predictors (baseline model)
#   "C_wide": 11 + 14 correlated null candidates (cost of selection)
#   "D_early": 11 + early-change marker (baseline-then-update model)
scenario_library <- function(predictor_cor = NULL) {
  list(
    C_primary = make_scenario(wide = FALSE, predictor_cor = predictor_cor),
    C_wide    = make_scenario(wide = TRUE),
    D_early   = make_scenario(wide = FALSE, early_rho = 0.5, predictor_cor = predictor_cor)
  )
}
