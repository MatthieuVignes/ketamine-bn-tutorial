suppressPackageStartupMessages(library(bnlearn))
set.seed(2026)

##################################################
## NETWORK A: baseline-only predictors (20 predictors + Response = 21 nodes)
## Mirrors the 20-variable list discussed with colleagues: demographics/illness course,
## comorbidity/history, treatment history, suicidality/functioning,
## neuropsychological domains, actigraphy RAR parameters -> Response
##################################################

modelstringA <- paste0(
  "[Age]",
  "[ChildhoodTrauma]",
  "[FamilyHistory]",
  "[Chronicity|Age]",
  "[ComorbidAnxiety|ChildhoodTrauma]",
  "[TxResistance|Chronicity]",
  "[PriorECT|TxResistance]",
  "[MedicationLoad|TxResistance]",
  "[BaselineSeverity|Chronicity:ComorbidAnxiety]",
  "[SuicidalIdeation|BaselineSeverity:ChildhoodTrauma]",
  "[FunctionalImpairment|BaselineSeverity]",
  "[PsychomotorSpeed|Age:BaselineSeverity]",
  "[ExecutiveFunction|Age:PsychomotorSpeed]",
  "[LearningMemory|Age:ExecutiveFunction]",
  "[AttentionWM|ExecutiveFunction]",
  "[IS|BaselineSeverity]",
  "[RA|IS]",
  "[SleepEfficiency|BaselineSeverity:RA]",
  "[SleepOnsetLatency|ComorbidAnxiety:SleepEfficiency]",
  "[TotalSleepTime|SleepEfficiency:SleepOnsetLatency]",
  "[Response|FamilyHistory:TxResistance:BaselineSeverity:PsychomotorSpeed:IS:ComorbidAnxiety]"
)

dagA <- model2network(modelstringA)
cat("Network A -- nodes:", length(nodes(dagA)), " arcs:", nrow(arcs(dagA)), "\n")

##################################################
## NETWORK B: baseline predictors + 3 early within-treatment "time"
## variables (early MADRS change, early actigraphy change, early sleep
## change), which in turn feed Response alongside baseline predictors.
## 21 + 3 = 24 nodes.
##################################################

modelstringB <- paste0(
  "[Age]",
  "[ChildhoodTrauma]",
  "[FamilyHistory]",
  "[Chronicity|Age]",
  "[ComorbidAnxiety|ChildhoodTrauma]",
  "[TxResistance|Chronicity]",
  "[PriorECT|TxResistance]",
  "[MedicationLoad|TxResistance]",
  "[BaselineSeverity|Chronicity:ComorbidAnxiety]",
  "[SuicidalIdeation|BaselineSeverity:ChildhoodTrauma]",
  "[FunctionalImpairment|BaselineSeverity]",
  "[PsychomotorSpeed|Age:BaselineSeverity]",
  "[ExecutiveFunction|Age:PsychomotorSpeed]",
  "[LearningMemory|Age:ExecutiveFunction]",
  "[AttentionWM|ExecutiveFunction]",
  "[IS|BaselineSeverity]",
  "[RA|IS]",
  "[SleepEfficiency|BaselineSeverity:RA]",
  "[SleepOnsetLatency|ComorbidAnxiety:SleepEfficiency]",
  "[TotalSleepTime|SleepEfficiency:SleepOnsetLatency]",
  "[EarlyMADRSChange|BaselineSeverity:TxResistance]",
  "[EarlyActigraphyChange|RA:ComorbidAnxiety]",
  "[EarlySleepChange|SleepEfficiency]",
  "[Response|FamilyHistory:TxResistance:BaselineSeverity:PsychomotorSpeed:IS:ComorbidAnxiety:EarlyMADRSChange:EarlyActigraphyChange]"
)

dagB <- model2network(modelstringB)
cat("Network B -- nodes:", length(nodes(dagB)), " arcs:", nrow(arcs(dagB)), "\n")

saveRDS(list(dagA = dagA, dagB = dagB), "data/true_dags.rds")

##################################################
## Visualise the two true DAGs (layered/Sugiyama layout so edge
## direction reads top-to-bottom; Response highlighted in red as the
## outcome node).
##################################################
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

