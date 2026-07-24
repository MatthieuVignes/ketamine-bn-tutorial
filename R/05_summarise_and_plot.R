suppressPackageStartupMessages(library(ggplot2))

res <- read.csv("data/simulation_results.csv")

## Derive network labels (with correct node counts) directly from the true
## DAGs, rather than hardcoding node counts as text -- avoids the labels
## going stale if the network definitions in 01_define_networks.R change.
dags <- readRDS("data/true_dags.rds")
net_labels <- c(
  A_baseline_only = sprintf("A: Baseline-only (%d nodes)", length(nodes(dags$dagA))),
  B_with_time     = sprintf("B: Baseline + early-change (%d nodes)", length(nodes(dags$dagB)))
)
res$network_label <- net_labels[res$network]

## Table: primary single-hc power curve
summ <- aggregate(
  cbind(precision, sensitivity, specificity, f1, shd,
        response_precision, response_recall) ~ network_label + n,
  data = res, FUN = function(x) mean(x, na.rm = TRUE)
)
summ <- summ[order(summ$network_label, summ$n), ]
write.csv(summ, "data/table1_primary_summary.csv", row.names = FALSE)
cat("=== Mean reconstruction metrics (single hc, 200 MC reps) ===\n")
print(round(summ[, -1], 3))

## Plot 1: skeleton sensitivity & precision vs n, by network
p1 <- ggplot(summ, aes(x = n)) +
  geom_line(aes(y = sensitivity, colour = "Sensitivity (recall)")) +
  geom_point(aes(y = sensitivity, colour = "Sensitivity (recall)")) +
  geom_line(aes(y = precision, colour = "Precision")) +
  geom_point(aes(y = precision, colour = "Precision")) +
  geom_line(aes(y = specificity, colour = "Specificity")) +
  geom_point(aes(y = specificity, colour = "Specificity")) +
  facet_wrap(~network_label) +
  scale_x_log10(breaks = c(20,30,40,50,60,75,100,150,200,300,500)) +
  labs(x = "Sample size (n, log scale)", y = "Mean value across 200 simulations",
       colour = NULL,
       title = "Skeleton reconstruction accuracy vs sample size",
       subtitle = "Score-based hill-climbing (BIC-Gaussian), whole-network skeleton") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 60, linetype = "dashed", colour = "grey40")
ggsave("figures/plot1_whole_network_accuracy.png", p1, width = 9, height = 5, dpi = 150)

## Plot 2: Response-node-specific precision/recall (the actual predictive model)
p2 <- ggplot(summ, aes(x = n)) +
  geom_line(aes(y = response_recall, colour = "Recall (sensitivity)")) +
  geom_point(aes(y = response_recall, colour = "Recall (sensitivity)")) +
  geom_line(aes(y = response_precision, colour = "Precision")) +
  geom_point(aes(y = response_precision, colour = "Precision")) +
  facet_wrap(~network_label) +
  scale_x_log10(breaks = c(20,30,40,50,60,75,100,150,200,300,500)) +
  ylim(0, 1) +
  labs(x = "Sample size (n, log scale)", y = "Mean value across 200 simulations",
       colour = NULL,
       title = "Recovery of Response's true predictors specifically",
       subtitle = "i.e. accuracy of the clinically-relevant part of the model") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 60, linetype = "dashed", colour = "grey40")
ggsave("figures/plot2_response_predictors_accuracy.png", p2, width = 9, height = 5, dpi = 150)

## Plot 3: SHD vs sample size
p3 <- ggplot(summ, aes(x = n, y = shd, colour = network_label)) +
  geom_line() + geom_point() +
  scale_x_log10(breaks = c(20,30,40,50,60,75,100,150,200,300,500)) +
  labs(x = "Sample size (n, log scale)", y = "Structural Hamming Distance (lower = better)",
       colour = NULL,
       title = "Structural Hamming Distance vs sample size",
       subtitle = "Score-based hill-climbing (BIC-Gaussian), 200 MC reps per point") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom") +
  geom_vline(xintercept = 60, linetype = "dashed", colour = "grey40")
ggsave("figures/plot3_shd_vs_n.png", p3, width = 8, height = 5, dpi = 150)

## ROC and Precision-Recall curves
## Score per possible edge = frequency with which it was learned across
## 200 independently simulated datasets at a given n (see edge_frequency.csv,
## produced in 04_run_simulation.R). Sweeping a threshold on that frequency
## gives a standard ROC / PR curve for "is this a real edge", at several
## sample sizes so the curves themselves show how much information n=60
## buys you relative to larger n.
edge_freq <- read.csv("data/edge_frequency.csv")
edge_freq$network_label <- net_labels[edge_freq$network]

roc_pr_from_scores <- function(score, label, thresholds = seq(0, 1, by = 0.01)) {
  out <- lapply(thresholds, function(t) {
    pred <- as.integer(score >= t)
    tp <- sum(pred == 1 & label == 1)
    fp <- sum(pred == 1 & label == 0)
    fn <- sum(pred == 0 & label == 1)
    tn <- sum(pred == 0 & label == 0)
    data.frame(
      threshold = t,
      tpr = if ((tp + fn) > 0) tp / (tp + fn) else NA,   # sensitivity/recall
      fpr = if ((fp + tn) > 0) fp / (fp + tn) else NA,   # 1 - specificity
      precision = if ((tp + fp) > 0) tp / (tp + fp) else NA
    )
  })
  do.call(rbind, out)
}

curve_list <- list()
for (net_label in unique(edge_freq$network)) {
  for (n_val in unique(edge_freq$n[edge_freq$network == net_label])) {
    sub <- edge_freq[edge_freq$network == net_label & edge_freq$n == n_val, ]
    curve <- roc_pr_from_scores(sub$frequency, sub$true_label)
    curve$network <- net_label
    curve$network_label <- net_labels[net_label]
    curve$n <- n_val
    curve_list[[paste(net_label, n_val)]] <- curve
  }
}
curves <- do.call(rbind, curve_list)
curves$n_f <- factor(curves$n, levels = sort(unique(curves$n)))
write.csv(curves, "data/roc_pr_curve_points.csv", row.names = FALSE)

## AUC (trapezoidal) per network x n, for a compact summary alongside the plots
auc_trap <- function(x, y) {
  o <- order(x)
  x <- x[o]; y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}
auc_summary <- do.call(rbind, lapply(split(curves, list(curves$network, curves$n)), function(sub) {
  sub <- sub[order(sub$fpr), ]
  data.frame(network_label = sub$network_label[1], n = sub$n[1],
             roc_auc = auc_trap(sub$fpr, sub$tpr))
}))
auc_summary <- auc_summary[order(auc_summary$network_label, auc_summary$n), ]
write.csv(auc_summary, "data/roc_auc_summary.csv", row.names = FALSE)
cat("\n=== ROC AUC by network and sample size ===\n")
print(auc_summary, row.names = FALSE)

## Plot 4: ROC curves
p4 <- ggplot(curves, aes(x = fpr, y = tpr, colour = n_f)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60") +
  geom_line(linewidth = 0.8) +
  facet_wrap(~network_label) +
  labs(x = "False positive rate (1 - specificity)", y = "True positive rate (sensitivity)",
       colour = "n",
       title = "ROC curves for edge recovery",
       subtitle = "Score = frequency of edge across 200 independent simulated datasets") +
  coord_equal() +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("figures/plot4_roc_curves.png", p4, width = 9, height = 5.5, dpi = 150)

## Plot 5: Precision-Recall curves
p5 <- ggplot(curves, aes(x = tpr, y = precision, colour = n_f)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~network_label) +
  ylim(0, 1) +
  labs(x = "Recall (sensitivity)", y = "Precision",
       colour = "n",
       title = "Precision-recall curves for edge recovery",
       subtitle = "Score = frequency of edge across 200 independent simulated datasets") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("figures/plot5_precision_recall_curves.png", p5, width = 9, height = 5.5, dpi = 150)

