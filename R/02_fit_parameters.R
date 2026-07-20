suppressPackageStartupMessages(library(bnlearn))
set.seed(2026)

dags <- readRDS("data/true_dags.rds")

## Assign each arc a "true" standardised path coefficient. Coefficients are
## drawn once (fixed seed) from +/-[0.30, 0.55]: moderate, realistic effect
## sizes for clinical/cognitive/actigraphy predictors (comparable to the
## small-moderate correlations reported in refs 10-15 of the EOI), not the
## unrealistically strong signals (e.g. |beta| > 0.8) that make structure
## learning artificially easy.
assign_coefs <- function(dag, sd_root = 1, sd_child = 1, seed = 1) {
  set.seed(seed)
  node_names <- nodes(dag)
  dist <- vector("list", length(node_names))
  names(dist) <- node_names

  for (nd in node_names) {
    pars <- parents(dag, nd)
    if (length(pars) == 0) {
      dist[[nd]] <- list(coef = c("(Intercept)" = 0), sd = sd_root)
    } else {
      betas <- runif(length(pars), 0.30, 0.55) * sample(c(-1, 1), length(pars), replace = TRUE)
      names(betas) <- pars
      coefs <- c("(Intercept)" = 0, betas)
      dist[[nd]] <- list(coef = coefs, sd = sd_child)
    }
  }
  dist
}

distA <- assign_coefs(dags$dagA, seed = 2024)
distB <- assign_coefs(dags$dagB, seed = 2025)

fittedA <- custom.fit(dags$dagA, dist = distA)
fittedB <- custom.fit(dags$dagB, dist = distB)

## Sanity check: simulate a large sample and look at R^2 for the Response
## node under each network, so the "signal" being asked of the learning
## algorithm is inspectable and plausible (not near-zero, not near-perfect).
chk <- function(fitted, label) {
  d <- rbn(fitted, n = 20000)
  parents_resp <- parents(fitted, "Response")
  fmla <- as.formula(paste("Response ~", paste(parents_resp, collapse = " + ")))
  r2 <- summary(lm(fmla, data = d))$r.squared
  cat(sprintf("[%s] Response ~ %s ; R^2 = %.3f\n", label, paste(parents_resp, collapse=", "), r2))
}
chk(fittedA, "Network A")
chk(fittedB, "Network B")

saveRDS(list(fittedA = fittedA, fittedB = fittedB), "data/fitted_true_bns.rds")

