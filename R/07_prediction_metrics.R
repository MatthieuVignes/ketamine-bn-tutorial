suppressPackageStartupMessages({
  library(bnlearn)
  library(glmnet)
})

##################################################
# Metrics

# AUC via the rank (Mann-Whitney) formula
auc_rank <- function(p, y) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Hanley-McNeil standard error of an AUC (approximate 95% CI)
auc_hm_ci <- function(auc, n1, n0) {
  if (is.na(auc) || n1 < 2 || n0 < 2) return(c(NA_real_, NA_real_))
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc * (1 - auc) + (n1 - 1) * (q1 - auc^2) + (n0 - 1) * (q2 - auc^2)) / (n1 * n0))
  c(auc - 1.96 * se, auc + 1.96 * se)
}

# Calibration slope: coef of logit(p) in a logistic regression of the obs resp 
# 1 = perfect; <1 = over-confident (small samples?)
calibration_slope <- function(p, y) {
  p <- pmin(pmax(p, 1e-4), 1 - 1e-4)
  fit <- tryCatch(suppressWarnings(glm(y ~ qlogis(p), family = binomial)), error = function(e) NULL)
  if (is.null(fit)) NA_real_ else unname(coef(fit)[2])
}

# Decision-oriented summary at a "treat if P(response) >= pt" rule
#   nb_gain: net benefit of the model minus net benefit of treating everyone
##  spared: proportion of NON-responders the rule would spare the trial
#   denied: proportion of RESPONDERS the rule would wrongly withhold treatment from
decision_summary <- function(p, y, pt = 0.40) {
  treat <- p >= pt
  N <- length(y); prev <- mean(y)
  tp <- sum(treat & y == 1); fp <- sum(treat & y == 0)
  nb_model <- tp / N - fp / N * pt / (1 - pt)
  nb_all   <- prev - (1 - prev) * pt / (1 - pt)
  c(nb_gain = nb_model - nb_all,
    spared  = if (sum(y == 0) > 0) sum(!treat & y == 0) / sum(y == 0) else NA_real_,
    denied  = if (sum(y == 1) > 0) sum(!treat & y == 1) / sum(y == 1) else NA_real_)
}

###################################################
# Learners which all return a function f(newdata) -> list(mu, s)
# mu: predicted Outcome and s: predictive SD.

outcome_predictors <- function(d) setdiff(names(d), c("Outcome", "Responder"))

# BN pipeline
# 1. learn structure by score-based hill climbing (BIC-Gaussian, restarts)
# 2. take the Markov blanket of Outcome (in a linear-Gaussian BN this is
#    exactly the set of variables the prediction of Outcome depends on)
# 3. fit the Outcome given its MB by least squares
# Note: `hc_args` are the same arguments as the structure-recovery simulation.
fit_bn_core <- function(train, restart = 5, perturb = 5) {
  vars <- c(outcome_predictors(train), "Outcome")
  g <- hc(train[, vars, drop = FALSE], score = "bic-g", restart = restart, perturb = perturb)
  blanket <- mb(g, "Outcome")
  fml <- if (length(blanket) == 0) Outcome ~ 1 else
    as.formula(paste("Outcome ~", paste(sprintf("`%s`", blanket), collapse = " + ")))
  m <- lm(fml, data = train)
  list(model = m, blanket = blanket, graph = g,
       s = sqrt(sum(residuals(m)^2) / max(1, m$df.residual)))
}

learner_bn <- function(train, restart = 5, perturb = 5) {
  core <- fit_bn_core(train, restart, perturb)
  function(newdata) list(mu = as.vector(predict(core$model, newdata)), s = core$s,
                         blanket = core$blanket)
}

# BN + bootstrap optimism correction ("stable network" step) estimate: 
# (a) a global shrinkage slope and (b) an honest predictive SD from out-of-bag 
# predictions of B bootstrap re-fits of the entire pipeline
learner_bn_boot <- function(train, B = 20, restart = 3, perturb = 5) {
  core <- fit_bn_core(train, restart, perturb)
  n <- nrow(train); ybar <- mean(train$Outcome)
  slopes <- numeric(0); oob_pred <- numeric(0); oob_y <- numeric(0)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    oob <- setdiff(seq_len(n), idx)
    cb <- tryCatch(fit_bn_core(train[idx, , drop = FALSE], restart, perturb), error = function(e) NULL)
    if (is.null(cb) || length(oob) < 3) next
    mu_all <- as.vector(predict(cb$model, train))
    sl <- tryCatch(unname(coef(lm(train$Outcome ~ mu_all))[2]), error = function(e) NA_real_)
    if (is.finite(sl)) slopes <- c(slopes, sl)
    oob_pred <- c(oob_pred, mu_all[oob]); oob_y <- c(oob_y, train$Outcome[oob])
  }
  shrink <- if (length(slopes) > 0) min(1, max(0, mean(slopes))) else 1
  s_oob <- if (length(oob_y) > 5) sqrt(mean((oob_y - (ybar + shrink * (oob_pred - ybar)))^2)) else core$s
  function(newdata) {
    mu <- as.vector(predict(core$model, newdata))
    list(mu = ybar + shrink * (mu - ybar), s = max(core$s, s_oob), blanket = core$blanket,
         shrink = shrink)
  }
}

# Competitors
learner_lasso_cont <- function(train) {
  vars <- outcome_predictors(train)
  cv <- cv.glmnet(as.matrix(train[, vars]), train$Outcome, nfolds = 5)
  s <- sqrt(min(cv$cvm))                           # cross-validated RMSE
  function(newdata) list(mu = as.vector(predict(cv, as.matrix(newdata[, vars]), s = "lambda.min")), s = s)
}

# Dichotomised-outcome lasso: returns probabilities directly
learner_lasso_logit <- function(train) {
  vars <- outcome_predictors(train)
  y <- train$Responder
  if (min(table(factor(y, levels = 0:1))) < 6) {
    p0 <- mean(y)
    return(function(newdata) list(p = rep(p0, nrow(newdata))))
  }
  cv <- cv.glmnet(as.matrix(train[, vars]), y, family = "binomial", nfolds = 5)
  function(newdata) list(p = as.vector(predict(cv, as.matrix(newdata[, vars]), s = "lambda.min", type = "response")))
}

# Clinical-only linear model ("does anything beyond cheap clinical variables 
# help?" benchmark)
learner_clinical <- function(train, clinical = CLINICAL_PREDICTORS) {
  m <- lm(as.formula(paste("Outcome ~", paste(clinical, collapse = " + "))), data = train)
  s <- sqrt(sum(residuals(m)^2) / m$df.residual)
  function(newdata) list(mu = as.vector(predict(m, newdata)), s = s)
}

# Oracle knows which predictors matter (not their coefficients). Upper
# bound on what any learner can achieve (at this n)?
learner_oracle <- function(train, true_vars = names(TRUE_EFFECTS)) {
  m <- lm(as.formula(paste("Outcome ~", paste(true_vars, collapse = " + "))), data = train)
  s <- sqrt(sum(residuals(m)^2) / m$df.residual)
  function(newdata) list(mu = as.vector(predict(m, newdata)), s = s)
}

# Convert learner's output to P(response)
to_prob <- function(pred, thr) {
  if (!is.null(pred$p)) return(pmin(pmax(pred$p, 1e-4), 1 - 1e-4))
  pmin(pmax(pnorm((pred$mu - thr) / pred$s), 1e-4), 1 - 1e-4)
}

# Which learners are run in the main simulation
LEARNERS <- list(
  BN_hc        = function(tr) learner_bn(tr),
  BN_hc_boot   = function(tr) learner_bn_boot(tr),
  Lasso_cont   = function(tr) learner_lasso_cont(tr),
  Lasso_logit  = function(tr) learner_lasso_logit(tr),
  Clinical_OLS = function(tr) learner_clinical(tr),
  Oracle       = function(tr) learner_oracle(tr)
)
