##################################################
## Install required packages.
## Run this once before anything else.
##################################################

required <- c("bnlearn", "ggplot2", "igraph", "parallel")

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
