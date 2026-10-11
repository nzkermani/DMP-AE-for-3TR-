## Run once in RStudio:  source("setup/install_packages.R")
## Installs everything the skills, scripts and MCP server use.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran <- c("jsonlite", "ggplot2", "dplyr", "tidyr", "tibble", "data.table",
          "lme4", "mgcv", "emmeans", "flexmix", "cluster", "patchwork",
          "msigdbr", "ellmer", "mcptools",
          # statistics, survival, ML, dimensionality reduction
          "survival", "survminer", "car", "MASS", "glmnet", "ranger", "e1071",
          "pls", "pROC", "mclust", "Boruta", "Rtsne", "uwot", "pheatmap")
bioc <- c("SummarizedExperiment", "MultiAssayExperiment", "S4Vectors",
          "BiocParallel", "limma", "edgeR", "variancePartition", "GSVA",
          "fgsea", "MOFA2", "mixOmics", "maSigPro", "ConsensusClusterPlus")

need_cran <- setdiff(cran, rownames(installed.packages()))
if (length(need_cran)) install.packages(need_cran)
need_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(need_bioc)) BiocManager::install(need_bioc, update = FALSE, ask = FALSE)

## LongOmicsTraj (declares R >= 4.6.0)
if (getRversion() < "4.6.0") {
  message("LongOmicsTraj requires R >= 4.6.0; you have ", getRversion(), ". Upgrade R first.")
} else if (!requireNamespace("LongOmicsTraj", quietly = TRUE)) {
  remotes::install_github("nzkermani/LongOmicsTraj")
}

all_pk <- c(cran, bioc, "LongOmicsTraj")
ok <- vapply(all_pk, requireNamespace, logical(1), quietly = TRUE)
print(data.frame(package = all_pk, installed = ok), row.names = FALSE)
if (!all(ok)) message("Missing: ", paste(all_pk[!ok], collapse = ", "))
