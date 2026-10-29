suppressPackageStartupMessages(library(ggplot2))

res <- read.csv("data/prediction_results.csv")

scen_lab <- c(C_primary = "11 pre-specified predictors",
              C_wide    = "25 candidates (14 extra null)",
              D_early   = "11 predictors + early-change marker")
res$scenario_label <- factor(scen_lab[res$scenario], levels = scen_lab)
res$learner <- factor(res$learner,
  levels = c("BN_hc", "BN_hc_boot", "Lasso_cont", "Lasso_logit", "Clinical_OLS", "Oracle"))

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
summ <- do.call(rbind, lapply(
  split(res, list(res$scenario, res$r2, res$n_train, res$learner), drop = TRUE),
  function(d) data.frame(
    scenario = d$scenario[1], r2 = d$r2[1], n_train = d$n_train[1],
    learner = as.character(d$learner[1]), n_sim = nrow(d),
    auc_true = mean(d$auc_true, na.rm = TRUE),
    auc_true_p10 = q(d$auc_true, .10), auc_true_p90 = q(d$auc_true, .90),
    cal_slope_median = median(d$cal_slope, na.rm = TRUE),
    brier = mean(d$brier, na.rm = TRUE),
    nb_gain = mean(d$nb_gain, na.rm = TRUE),
    spared = mean(d$spared, na.rm = TRUE), denied = mean(d$denied, na.rm = TRUE),
    p_auc_ge_070 = mean(d$auc_val >= 0.70, na.rm = TRUE),
    p_lb_gt_060  = mean(d$auc_val_lb > 0.60, na.rm = TRUE),
    mean_blanket = mean(d$n_blanket, na.rm = TRUE))))
summ <- summ[order(summ$scenario, summ$n_train, summ$r2, summ$learner), ]
num <- vapply(summ, is.numeric, logical(1)); summ[num] <- lapply(summ[num], round, 3)
write.csv(summ, "data/table2_prediction_summary.csv", row.names = FALSE)

cat("\n=== True out-of-sample AUC / calibration slope / assurance (Pr[lower 95% CI of validation AUC > 0.60]) ===\n")
print(summ[summ$learner %in% c("BN_hc", "BN_hc_boot", "Lasso_cont", "Clinical_OLS", "Oracle"),
           c("scenario", "n_train", "r2", "learner", "n_sim", "auc_true", "auc_true_p10", "auc_true_p90",
             "cal_slope_median", "p_auc_ge_070", "p_lb_gt_060")], row.names = FALSE)

dir.create("figures", showWarnings = FALSE)
cols <- c(BN_hc = "#1b6ca8", BN_hc_boot = "#0b3d63", Lasso_cont = "#d98a1c",
          Lasso_logit = "#a0a0a0", Clinical_OLS = "#6a994e", Oracle = "black")

# Plot 6: true out-of-sample AUC vs signal strength
p6 <- ggplot(summ, aes(r2, auc_true, colour = learner, group = learner)) +
  geom_hline(yintercept = 0.70, linetype = "dashed", colour = "grey50") +
  geom_ribbon(data = subset(summ, learner == "BN_hc_boot" | learner == "BN_hc"),
              aes(ymin = auc_true_p10, ymax = auc_true_p90, fill = learner),
              alpha = 0.10, colour = NA, show.legend = FALSE) +
  geom_line() + geom_point(size = 1.4) +
  facet_grid(n_train ~ scenario, labeller = labeller(
    n_train = function(x) paste("train n =", x), scenario = scen_lab)) +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(x = "True signal strength (R² of improvement score)", y = "True out-of-sample AUC",
       colour = NULL, title = "How well can we predict response in new patients?",
       subtitle = "Mean over simulated studies; band = 10-90% range for the BN pipeline; dashed = AUC 0.70. (Oracle ignores the early marker.)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/plot6_prediction_auc.png", p6, width = 12, height = 7.5, dpi = 150)

# Plot 7: assurance - chance the n = 60 validation unambiguously shows AUC > 0.60
p7 <- ggplot(subset(summ, learner %in% c("BN_hc", "BN_hc_boot", "Lasso_cont", "Oracle")),
             aes(r2, p_lb_gt_060, colour = learner, group = learner)) +
  geom_line() + geom_point(size = 1.4) + ylim(0, 1) +
  facet_grid(n_train ~ scenario, labeller = labeller(
    n_train = function(x) paste("train n =", x), scenario = scen_lab)) +
  scale_colour_manual(values = cols) +
  labs(x = "True signal strength (R²)", y = "Pr(lower 95% bound of validation AUC > 0.60)",
       colour = NULL, title = "Chance the n = 60 validation set confirms useful discrimination",
       subtitle = "Observed AUC in a simulated n = 60 independent validation set (Hanley-McNeil 95% CI)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/plot7_validation_assurance.png", p7, width = 12, height = 7.5, dpi = 150)

# Plot 8: calibration slope (1 = ideal; <1 = over-confident predicted probabilities)
p8 <- ggplot(subset(summ, learner %in% c("BN_hc", "BN_hc_boot", "Lasso_cont", "Lasso_logit")),
             aes(r2, cal_slope_median, colour = learner, group = learner)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line() + geom_point(size = 1.4) +
  facet_grid(n_train ~ scenario, labeller = labeller(
    n_train = function(x) paste("train n =", x), scenario = scen_lab)) +
  scale_colour_manual(values = cols) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(x = "True signal strength (R²)", y = "Median calibration slope (large test set)",
       colour = NULL, title = "Are predicted probabilities trustworthy?",
       subtitle = "Slope < 1: predictions too extreme. Bootstrap shrinkage (BN_hc_boot) is the correction we will apply.") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/plot8_calibration.png", p8, width = 12, height = 7.5, dpi = 150)

# Plot 9: how much does an early-change marker add? (BN only, same r2 grid)
d9 <- subset(summ, learner == "BN_hc_boot" & scenario %in% c("C_primary", "D_early"))
p9 <- ggplot(d9, aes(r2, auc_true, colour = scenario, group = scenario)) +
  geom_line() + geom_point() + facet_wrap(~n_train, labeller = labeller(n_train = function(x) paste("train n =", x))) +
  scale_colour_manual(values = c(C_primary = "#1b6ca8", D_early = "#c1440e"), labels = scen_lab[c("C_primary", "D_early")]) +
  labs(x = "True signal strength (R²)", y = "True out-of-sample AUC", colour = NULL,
       title = "Baseline model vs baseline + early-change update (BN pipeline)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("figures/plot9_early_update.png", p9, width = 10, height = 4.5, dpi = 150)
