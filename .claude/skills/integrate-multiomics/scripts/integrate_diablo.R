#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## integrate_diablo.R — supervised multi-omics integration (mixOmics DIABLO)
## discriminating VISITS within subjects, i.e. which coordinated
## multi-omic signature changes over time.
##
## Usage (from the project root):
##   Rscript .claude/skills/integrate-multiomics/scripts/integrate_diablo.R \
##       [config] [group|all] [ncomp] [keepX] [top_var_features_per_layer]
##   group  restrict to one arm (recommended: run each arm separately, then
##          compare signatures); "all" pools arms
##   defaults: all, ncomp = 3, keepX = 20 per block per component, top_var = 2000
##
## Repeated measures: each block is first passed through
## mixOmics::withinVariation() with subject as the design, removing
## between-subject variation (multilevel decomposition). This requires
## samples present in EVERY layer; subjects with only one such sample are
## dropped because their within-subject variation is undefined.
##
## Outputs (results/integration/diablo[_<group>]/):
##   selected_features.csv, sample_scores.csv, block_correlations.csv,
##   diablo_plots.pdf, perf.rds (if run_perf = TRUE)
##
## Status: written against the mixOmics API (withinVariation, block.splsda,
## selectVar, plotIndiv, circosPlot, perf); NOT executed in the environment
## this kit was built in (mixOmics unavailable there). keepX is NOT tuned
## here — use mixOmics::tune.block.splsda before reporting results.
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(mixOmics)
})
source("R/lomics_utils.R")

args <- commandArgs(trailingOnly = TRUE)
a <- function(i, d) if (length(args) >= i && nzchar(args[i])) args[i] else d
cfg     <- lomics_read_config(a(1, "study_config.json"))
grp     <- a(2, "all")
ncomp   <- as.integer(a(3, "3"))
keepX_n <- as.integer(a(4, "20"))
top_var <- as.integer(a(5, "2000"))
run_perf <- FALSE
layers  <- lomics_layer_names(cfg)

ses <- lapply(stats::setNames(layers, layers), function(l) lomics_load_layer(cfg, l, "processed"))
common <- Reduce(intersect, lapply(ses, colnames))
meta <- as.data.frame(colData(ses[[1]]))[common, ]
if (grp != "all") meta <- meta[meta$group == grp, ]
## keep subjects with >= 2 visits in the complete-case set
n_per <- table(meta$subject_id)
meta <- meta[meta$subject_id %in% names(n_per)[n_per >= 2], ]
meta <- meta[order(meta$subject_id, meta$visit), ]
message(sprintf("DIABLO on %d samples from %d subjects (all layers present)",
                nrow(meta), length(unique(meta$subject_id))))

X <- lapply(ses, function(se) {
  x <- assay(se, "values")[, meta$sample_id, drop = FALSE]
  v <- apply(x, 1, stats::var)
  x <- x[order(v, decreasing = TRUE)[seq_len(min(top_var, nrow(x)))], , drop = FALSE]
  t(x)                                                 # samples x features
})
X <- lapply(X, function(b) mixOmics::withinVariation(X = b, design = data.frame(meta$subject_id)))
Y <- factor(meta$visit, levels = cfg$visits)

design <- matrix(0.1, length(X), length(X), dimnames = list(names(X), names(X)))
diag(design) <- 0
keepX <- lapply(X, function(b) rep(min(keepX_n, ncol(b)), ncomp))
fit <- mixOmics::block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design)

outdir <- dirname(lomics_out(cfg, "integration",
                             if (grp == "all") "diablo" else paste0("diablo_", make.names(grp)), "x"))
sel <- do.call(rbind, lapply(seq_len(ncomp), function(cmp) {
  do.call(rbind, lapply(names(X), function(b) {
    s <- mixOmics::selectVar(fit, block = b, comp = cmp)[[b]]$value
    if (is.null(s) || !nrow(s)) return(NULL)
    data.frame(component = cmp, layer = b, feature = rownames(s), loading = s$value.var)
  }))
}))
utils::write.csv(sel, file.path(outdir, "selected_features.csv"), row.names = FALSE)

scores <- do.call(cbind, lapply(names(fit$variates)[names(fit$variates) != "Y"], function(b) {
  v <- fit$variates[[b]]; colnames(v) <- paste(b, colnames(v), sep = "_"); v
}))
utils::write.csv(data.frame(meta[, c("sample_id", "subject_id", "visit", "group")], scores),
                 file.path(outdir, "sample_scores.csv"), row.names = FALSE)

## Correlation between block components (component 1)
c1 <- sapply(names(X), function(b) fit$variates[[b]][, 1])
utils::write.csv(stats::cor(c1), file.path(outdir, "block_correlations.csv"))

pdf(file.path(outdir, "diablo_plots.pdf"), width = 9, height = 7)
try(mixOmics::plotIndiv(fit, ind.names = FALSE, legend = TRUE, title = "DIABLO (within-subject)"))
try(mixOmics::plotDiablo(fit, ncomp = 1))
try(mixOmics::circosPlot(fit, cutoff = 0.7))
invisible(dev.off())

if (run_perf) {
  pf <- mixOmics::perf(fit, validation = "Mfold", folds = 5, nrepeat = 10)
  saveRDS(pf, file.path(outdir, "perf.rds"))
}
lomics_session_log(cfg, "diablo")
