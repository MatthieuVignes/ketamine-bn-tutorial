suppressPackageStartupMessages(library(bnlearn))

## Undirected skeleton (set of unordered node pairs with an edge) of a BN
skeleton_pairs <- function(bn_obj) {
  a <- arcs(bn_obj)
  if (nrow(a) == 0) return(character(0))
  pairs <- apply(a, 1, function(r) paste(sort(r), collapse = "--"))
  unique(pairs)
}

all_possible_pairs <- function(node_names) {
  cmb <- combn(sort(node_names), 2)
  apply(cmb, 2, function(r) paste(r, collapse = "--"))
}

## Data frame of every possible undirected node pair for a DAG's node set,
## with a true_label (1 = arc exists between them in the true DAG, 0 = not).
## Used as the labelled "universe" for ROC / precision-recall curves, where
## the score per pair is how often that edge was learned across repeated
## simulated datasets (see 04_run_simulation.R's edge-frequency output).
skeleton_universe <- function(true_dag) {
  node_names <- nodes(true_dag)
  cmb <- combn(sort(node_names), 2)
  pair    <- apply(cmb, 2, function(r) paste(r, collapse = "--"))
  node1   <- cmb[1, ]
  node2   <- cmb[2, ]
  true_edges <- skeleton_pairs(true_dag)
  data.frame(node1 = node1, node2 = node2, pair = pair,
             true_label = as.integer(pair %in% true_edges),
             stringsAsFactors = FALSE)
}

## Skeleton-level precision / recall / specificity / F1, plus SHD (on CPDAGs,
## via bnlearn::shd) and a targeted metric for recovery of of the edges
## touching Response. More specifically its parent set, since Response is a sink node
## (no children) in both networks here. Note this is NOT the general
## definition of a Markov blanket (parents + children + spouses/co-parents
## of children) -- it only coincides with it because Response has no
## children. If Response ever gains a child (e.g. a downstream outcome),
## this metric would need to include that child's other parents too.
## That later is probably what a clinician-facing predictive model actually needs.
evaluate_reconstruction <- function(true_dag, learned_dag) {
  nodes_all   <- nodes(true_dag)
  universe    <- all_possible_pairs(nodes_all)
  true_edges  <- skeleton_pairs(true_dag)
  learned_edges <- skeleton_pairs(learned_dag)

  tp <- length(intersect(true_edges, learned_edges))
  fp <- length(setdiff(learned_edges, true_edges))
  fn <- length(setdiff(true_edges, learned_edges))
  tn <- length(universe) - tp - fp - fn

  precision   <- if ((tp + fp) > 0) tp / (tp + fp) else NA
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  f1          <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * precision * sensitivity / (precision + sensitivity)
  } else NA

  shd_val <- shd(learned_dag, true_dag)

  ## Response-specific recovery: edges (undirected) touching "Response"
  resp_true    <- true_edges[grepl("(^|--)Response(--|$)", true_edges)]
  resp_learned <- learned_edges[grepl("(^|--)Response(--|$)", learned_edges)]
  resp_universe <- universe[grepl("(^|--)Response(--|$)", universe)]
  r_tp <- length(intersect(resp_true, resp_learned))
  r_fp <- length(setdiff(resp_learned, resp_true))
  r_fn <- length(setdiff(resp_true, resp_learned))
  r_precision <- if ((r_tp + r_fp) > 0) r_tp / (r_tp + r_fp) else NA
  r_recall    <- if ((r_tp + r_fn) > 0) r_tp / (r_tp + r_fn) else NA

  data.frame(tp = tp, fp = fp, fn = fn, tn = tn,
    precision = precision, sensitivity = sensitivity,
    specificity = specificity, f1 = f1, shd = shd_val,
    response_precision = r_precision, response_recall = r_recall,
    response_true_edges = length(resp_true)
  )
}

