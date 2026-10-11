#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## longitudinal_de.R — repeated-measures differential analysis over 4 visits.
##
## Usage (from the project root):
##   Rscript .claude/skills/stats-longitudinal-de/scripts/longitudinal_de.R \
##       [config] [layer|all] [method]
##   method = "limma" (default; duplicateCorrelation on subject)
##          | "dream" (variancePartition; subject random intercept; better for
##                     unbalanced follow-up or when you want extra random terms)
##
## Model: ~ 0 + cell + covariates, where cell = group x visit, with the
## within-subject correlation handled by duplicateCorrelation (limma) or a
## (1 | subject_id) random effect (dream).
##
## Contrasts written per layer (results/de/<layer>/):
##   change_<group>_<Vk>vsV1      within-group change from baseline
##   step_<group>_<Vk>vs<Vk-1>    consecutive-visit change (transition)
##   interaction_<G2>vs<G1>_<Vk>  difference in change-from-baseline between
##                                groups (only when >= 2 groups)
##   overall_time_<group>         F-test across all visits within group
## plus de_summary.csv with counts at FDR < 0.05.
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
})
source("R/lomics_utils.R")

args   <- commandArgs(trailingOnly = TRUE)
cfg    <- lomics_read_config(if (length(args) >= 1) args[1] else "study_config.json")
layers <- if (length(args) >= 2 && args[2] != "all") args[2] else lomics_layer_names(cfg)
method <- if (length(args) >= 3) args[3] else "limma"
fdr    <- 0.05

safe <- function(x) make.names(x)

build_design <- function(meta, covars) {
  meta$cell <- factor(paste(safe(meta$group), safe(as.character(meta$visit)), sep = "."))
  covars <- intersect(covars, colnames(meta))
  ## drop covariates that are constant or fully missing
  covars <- covars[vapply(covars, function(cv) length(unique(stats::na.omit(meta[[cv]]))) > 1, logical(1))]
  f <- stats::as.formula(paste("~ 0 + cell", if (length(covars)) paste("+", paste(covars, collapse = " + "))))
  list(meta = meta, formula = f, covars = covars)
}

make_contrasts <- function(groups, visits, design_cols) {
  g <- safe(groups); v <- safe(visits)
  cn <- function(gg, vv) paste0("cell", gg, ".", vv)
  ok <- function(...) all(c(...) %in% design_cols)
  cons <- list()
  for (gg in g) {
    for (k in seq_along(v)[-1]) {
      if (ok(cn(gg, v[k]), cn(gg, v[1])))
        cons[[sprintf("change_%s_%svs%s", gg, v[k], v[1])]] <- sprintf("%s - %s", cn(gg, v[k]), cn(gg, v[1]))
      if (k > 2 && ok(cn(gg, v[k]), cn(gg, v[k - 1])))
        cons[[sprintf("step_%s_%svs%s", gg, v[k], v[k - 1])]] <- sprintf("%s - %s", cn(gg, v[k]), cn(gg, v[k - 1]))
    }
  }
  if (length(g) >= 2) {
    ref <- g[1]
    for (gg in g[-1]) for (k in seq_along(v)[-1]) {
      if (ok(cn(gg, v[k]), cn(gg, v[1]), cn(ref, v[k]), cn(ref, v[1])))
        cons[[sprintf("interaction_%svs%s_%s", gg, ref, v[k])]] <-
          sprintf("(%s - %s) - (%s - %s)", cn(gg, v[k]), cn(gg, v[1]), cn(ref, v[k]), cn(ref, v[1]))
    }
  }
  cons
}

summary_rows <- list()
for (layer in layers) {
  se <- lomics_load_layer(cfg, layer, "processed")
  meta <- as.data.frame(colData(se))
  L <- cfg$layers[cfg$layers$name == layer, ]
  scale_l <- if (nrow(L) == 1) L$scale else "scores"   # derived layers, e.g. pathways_*
  is_counts <- identical(scale_l, "counts") && "counts" %in% assayNames(se)
  bd <- build_design(meta, cfg$covariates)
  meta <- bd$meta
  cc <- stats::complete.cases(meta[, c("cell", bd$covars), drop = FALSE])
  if (any(!cc)) message(layer, ": dropping ", sum(!cc), " samples with missing covariates")
  meta <- meta[cc, , drop = FALSE]; se <- se[, cc]
  design <- stats::model.matrix(bd$formula, data = meta)
  colnames(design) <- make.names(colnames(design))
  cons <- make_contrasts(unique(meta$group), cfg$visits, colnames(design))
  cmat <- limma::makeContrasts(contrasts = unlist(cons), levels = design)
  colnames(cmat) <- names(cons)
  outdir <- lomics_out(cfg, "de", layer, "x"); outdir <- dirname(outdir)

  if (method == "limma") {
    if (is_counts) {
      dge <- edgeR::DGEList(assay(se, "counts"))
      dge$samples$norm.factors <- S4Vectors::metadata(se)$norm_factors[cc]
      v <- limma::voom(dge, design)
      dc <- limma::duplicateCorrelation(v, design, block = meta$subject_id)
      v <- limma::voom(dge, design, block = meta$subject_id, correlation = dc$consensus.correlation)
      dc <- limma::duplicateCorrelation(v, design, block = meta$subject_id)
      y <- v
    } else {
      y <- assay(se, "values")
      dc <- limma::duplicateCorrelation(y, design, block = meta$subject_id)
    }
    message(sprintf("%s: consensus within-subject correlation = %.3f", layer, dc$consensus.correlation))
    fit <- limma::lmFit(y, design, block = meta$subject_id, correlation = dc$consensus.correlation)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cmat), trend = !is_counts, robust = TRUE)
    get_tab <- function(coef) limma::topTable(fit2, coef = coef, number = Inf, sort.by = "P")
  } else if (method == "dream") {
    if (!requireNamespace("variancePartition", quietly = TRUE)) stop("Install variancePartition for method = 'dream'.")
    form <- stats::as.formula(paste(paste(deparse(bd$formula), collapse = ""), "+ (1 | subject_id)"))
    if (is_counts) {
      dge <- edgeR::DGEList(assay(se, "counts"))
      dge$samples$norm.factors <- S4Vectors::metadata(se)$norm_factors[cc]
      y <- variancePartition::voomWithDreamWeights(dge, form, meta)
    } else {
      y <- assay(se, "values")
    }
    L_mat <- variancePartition::makeContrastsDream(form, meta, contrasts = unlist(cons))
    colnames(L_mat) <- names(cons)
    fit <- variancePartition::dream(y, form, meta, L = L_mat)
    fit2 <- variancePartition::eBayes(fit)
    get_tab <- function(coef) variancePartition::topTable(fit2, coef = coef, number = Inf, sort.by = "P")
  } else stop("Unknown method: ", method)

  for (nm in names(cons)) {
    tt <- get_tab(nm)
    tt <- data.frame(feature = rownames(tt), tt, row.names = NULL, check.names = FALSE)
    utils::write.csv(tt, file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      layer = layer, contrast = nm, n_features = nrow(tt),
      n_up = sum(tt$adj.P.Val < fdr & tt$logFC > 0),
      n_down = sum(tt$adj.P.Val < fdr & tt$logFC < 0))
  }
  ## Overall time effect within each group (F-test over change-from-baseline)
  if (method == "limma") {
    for (gg in safe(unique(meta$group))) {
      cols <- grep(sprintf("^change_%s_", gg), colnames(cmat), value = TRUE)
      if (length(cols) >= 2) {
        tt <- limma::topTable(fit2, coef = cols, number = Inf, sort.by = "F")
        tt <- data.frame(feature = rownames(tt), tt, row.names = NULL, check.names = FALSE)
        utils::write.csv(tt, file.path(outdir, sprintf("overall_time_%s.csv", gg)), row.names = FALSE)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          layer = layer, contrast = sprintf("overall_time_%s", gg), n_features = nrow(tt),
          n_up = sum(tt$adj.P.Val < fdr), n_down = NA)
      }
    }
  }
  writeLines(c(sprintf("method: %s", method), sprintf("formula: %s", deparse(bd$formula)),
               if (method == "limma") sprintf("within-subject correlation: %.3f", dc$consensus.correlation),
               "contrasts:", paste(" ", names(cons), "=", unlist(cons))),
             file.path(outdir, "model_spec.txt"))
}
smry <- do.call(rbind, summary_rows)
sfile <- lomics_out(cfg, "de", "de_summary.csv")
if (file.exists(sfile)) {                     # keep layers from earlier runs
  prev <- utils::read.csv(sfile)
  smry_all <- rbind(prev[!prev$layer %in% layers, , drop = FALSE], smry)
} else smry_all <- smry
utils::write.csv(smry_all, sfile, row.names = FALSE)
print(smry, row.names = FALSE)
lomics_session_log(cfg, paste0("de_", method))
