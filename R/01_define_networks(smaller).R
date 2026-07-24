suppressPackageStartupMessages(library(bnlearn))
set.seed(2026)

## ------------------------------------------------------------------
## NETWORK A: baseline-only predictors (15 predictors + Response = 16 nodes)
##
## Reduced from an original 20-variable candidate list by:
##  - collapsing 4 neuropsych domain scores down to 2: PsychomotorSpeed is
##    kept as its own node (specific ketamine-response evidence in the
##    cited literature), while ExecutiveFunction/LearningMemory/AttentionWM
##    are merged into a single HigherCognition composite
##  - collapsing 3 actigraphy sleep metrics (SleepEfficiency,
##    SleepOnsetLatency, TotalSleepTime) into a single SleepQuality
##    composite
##  - dropping FamilyHistory (weakest/most removable predictor)
## This keeps the domains that have the most direct ketamine-specific
## evidence as standalone nodes, and folds generic/correlated predictors
## into composites rather than discarding them outright.
## ------------------------------------------------------------------

modelstringA <- paste0(
  "[Age]",
  "[ChildhoodTrauma]",
  "[Chronicity|Age]",
  "[ComorbidAnxiety|ChildhoodTrauma]",
  "[TxResistance|Chronicity]",
  "[PriorECT|TxResistance]",
  "[MedicationLoad|TxResistance]",
  "[BaselineSeverity|Chronicity:ComorbidAnxiety]",
  "[SuicidalIdeation|BaselineSeverity:ChildhoodTrauma]",
  "[FunctionalImpairment|BaselineSeverity]",
  "[PsychomotorSpeed|Age:BaselineSeverity]",
  "[HigherCognition|Age:PsychomotorSpeed]",
  "[IS|BaselineSeverity]",
  "[RA|IS]",
  "[SleepQuality|BaselineSeverity:RA:ComorbidAnxiety]",
  "[Response|TxResistance:BaselineSeverity:PsychomotorSpeed:IS:ComorbidAnxiety]"
)

dagA <- model2network(modelstringA)
cat("Network A -- nodes:", length(nodes(dagA)), " arcs:", nrow(arcs(dagA)), "\n")

## ------------------------------------------------------------------
## NETWORK B: baseline predictors + 3 early within-treatment "time"
## variables (early MADRS change, early actigraphy change, early sleep
## change), which in turn feed Response alongside baseline predictors.
## 15 + 3 = 18, + Response = 19 nodes.
## ------------------------------------------------------------------

modelstringB <- paste0(
  "[Age]",
  "[ChildhoodTrauma]",
  "[Chronicity|Age]",
  "[ComorbidAnxiety|ChildhoodTrauma]",
  "[TxResistance|Chronicity]",
  "[PriorECT|TxResistance]",
  "[MedicationLoad|TxResistance]",
  "[BaselineSeverity|Chronicity:ComorbidAnxiety]",
  "[SuicidalIdeation|BaselineSeverity:ChildhoodTrauma]",
  "[FunctionalImpairment|BaselineSeverity]",
  "[PsychomotorSpeed|Age:BaselineSeverity]",
  "[HigherCognition|Age:PsychomotorSpeed]",
  "[IS|BaselineSeverity]",
  "[RA|IS]",
  "[SleepQuality|BaselineSeverity:RA:ComorbidAnxiety]",
  "[EarlyMADRSChange|BaselineSeverity:TxResistance]",
  "[EarlyActigraphyChange|RA:ComorbidAnxiety]",
  "[EarlySleepChange|SleepQuality]",
  "[Response|TxResistance:BaselineSeverity:PsychomotorSpeed:IS:ComorbidAnxiety:EarlyMADRSChange:EarlyActigraphyChange]"
)

dagB <- model2network(modelstringB)
cat("Network B -- nodes:", length(nodes(dagB)), " arcs:", nrow(arcs(dagB)), "\n")

saveRDS(list(dagA = dagA, dagB = dagB), "data/true_dags.rds")

## ------------------------------------------------------------------
## Visualise the two true DAGs (layered/Sugiyama layout so edge
## direction reads top-to-bottom; Response highlighted in red as the
## outcome node).
## ------------------------------------------------------------------
suppressPackageStartupMessages(library(igraph))

plot_bn <- function(dag, title, file, highlight = "Response") {
  el <- arcs(dag)
  g  <- graph_from_edgelist(el, directed = TRUE)
  missing_nodes <- setdiff(nodes(dag), V(g)$name)
  if (length(missing_nodes) > 0) g <- add_vertices(g, length(missing_nodes), name = missing_nodes)

  lay <- layout_with_sugiyama(g, hgap = 1, vgap = 1.4)$layout

  vcol    <- ifelse(V(g)$name == highlight, "#e31a1c", "#a6cee3")
  vsize   <- ifelse(V(g)$name == highlight, 20, 14)
  lab_col <- ifelse(V(g)$name == highlight, "white", "black")

  png(file, width = 1600, height = 1100, res = 150)
  par(mar = c(1, 1, 3, 1))
  plot(g, layout = lay, vertex.size = vsize, vertex.color = vcol,
       vertex.label.color = lab_col, vertex.label.cex = 0.68,
       vertex.label.family = "sans", edge.arrow.size = 0.35,
       edge.color = "grey50", vertex.frame.color = "grey30",
       main = title)
  dev.off()
}

plot_bn(dagA, sprintf("Network A: baseline-only (%d nodes, %d arcs)",
                       length(nodes(dagA)), nrow(arcs(dagA))),
        "figures/networkA_diagram.png")
plot_bn(dagB, sprintf("Network B: baseline + early-change (%d nodes, %d arcs)",
                       length(nodes(dagB)), nrow(arcs(dagB))),
        "figures/networkB_diagram.png")
cat("Saved network diagrams: networkA_diagram.png, networkB_diagram.png\n")
