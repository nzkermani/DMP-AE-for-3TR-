#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## qc_normalise.R — per-layer QC, filtering, imputation and normalisation.
##
## Usage (from the project root):
##   Rscript .claude/skills/prep-qc-normalisation/scripts/qc_normalise.R [config] [layer]
##   (omit `layer` to process every layer in the config)
##
## Behaviour is chosen by the layer's `scale` in study_config.json and can be
## overridden by optional per-layer fields:
##   missing_max  (default 0.3)   drop features missing in > this fraction of
##                                samples in EVERY visit (visit-aware filter)
##   impute       "halfmin" | "minprob" | "none"   (default: halfmin; counts: none)
##   normalise    "tmm" (counts) | "median" | "pqn" | "none"
##                defaults: counts -> tmm, log2 -> none, raw_intensity -> pqn
##   qc_rsd_max   (default 0.3)   only used when config$qc_sample_column marks
##                                pooled-QC injections ("QC") in the metadata
##
## Outputs:
##   results/data/<layer>_processed_se.rds   assay "values" = log2 scale
##                                           (counts layers also keep "counts")
##   results/qc/<layer>_qc_summary.txt
##   results/qc/<layer>_pca.pdf              PCA coloured by visit, group, batch
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
})
source("R/lomics_utils.R")

args  <- commandArgs(trailingOnly = TRUE)
cfg   <- lomics_read_config(if (length(args) >= 1) args[1] else "study_config.json")
which_layers <- if (length(args) >= 2) args[2] else lomics_layer_names(cfg)

opt <- function(L, field, default) {
  v <- if (field %in% names(L)) L[[field]] else NULL
  if (is.null(v) || length(v) == 0 || is.na(v)) default else v
}

## Visit-aware missingness filter: keep a feature if, in at least one visit,
## its missing fraction is <= missing_max. Avoids discarding features that
## switch on/off over time (e.g. treatment-induced metabolites).
visit_aware_keep <- function(x, visit, missing_max) {
  per_visit <- sapply(split(seq_len(ncol(x)), visit), function(j) {
    rowMeans(is.na(x[, j, drop = FALSE]))
  })
  per_visit <- matrix(per_visit, nrow = nrow(x))
  apply(per_visit, 1, min) <= missing_max
}

impute_halfmin <- function(x, logged) {
  for (i in which(rowSums(is.na(x)) > 0)) {
    mn <- min(x[i, ], na.rm = TRUE)
    x[i, is.na(x[i, ])] <- if (logged) mn - 1 else mn / 2
  }
  x
}

## MinProb-style left-censored imputation on the log scale: draws from a
## normal centred on the 1st percentile of each sample, sd = 0.3 * median
## feature sd. Seeded for reproducibility.
impute_minprob <- function(x, q = 0.01, tune_sigma = 0.3, seed = 1) {
  set.seed(seed)
  sdf <- stats::median(apply(x, 1, stats::sd, na.rm = TRUE), na.rm = TRUE)
  for (j in seq_len(ncol(x))) {
    na <- is.na(x[, j])
    if (any(na)) {
      mu <- stats::quantile(x[, j], q, na.rm = TRUE)
      x[na, j] <- stats::rnorm(sum(na), mu, tune_sigma * sdf)
    }
  }
  x
}

pqn <- function(x) {
  ref <- apply(x, 1, stats::median, na.rm = TRUE)
  q <- apply(x / ref, 2, stats::median, na.rm = TRUE)
  sweep(x, 2, q, "/")
}

pca_plots <- function(v, meta, layer, batch_col) {
  v <- v[apply(v, 1, stats::var) > 0, , drop = FALSE]
  pc <- stats::prcomp(t(v), scale. = TRUE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  d <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], meta)
  cols <- intersect(c("visit", "group", batch_col), colnames(d))
  lapply(cols, function(cl) {
    ggplot(d, aes(PC1, PC2, colour = .data[[cl]])) +
      geom_point(size = 2, alpha = .8) +
      labs(title = sprintf("%s — PCA coloured by %s", layer, cl),
           x = sprintf("PC1 (%s%%)", ve[1]), y = sprintf("PC2 (%s%%)", ve[2])) +
      theme_minimal()
  })
}

for (layer in which_layers) {
  L <- cfg$layers[cfg$layers$name == layer, , drop = FALSE]
  if (nrow(L) != 1) stop("Layer not in config: ", layer, call. = FALSE)
  L <- as.list(L)
  se <- lomics_load_layer(cfg, layer, "raw")
  x <- assay(se, "raw"); meta <- as.data.frame(colData(se))
  log <- c(sprintf("Layer: %s  (scale: %s)", layer, L$scale),
           sprintf("Input: %d features x %d samples", nrow(x), ncol(x)))

  missing_max <- as.numeric(opt(L, "missing_max", 0.3))
  default_norm <- switch(L$scale, counts = "tmm", log2 = "none",
                         raw_intensity = "pqn", "none")
  normalise <- opt(L, "normalise", default_norm)
  impute <- opt(L, "impute", if (L$scale == "counts") "none" else "halfmin")

  ## Pooled-QC RSD filter (metabolomics / MS proteomics), if QC samples exist
  qc_col <- cfg$qc_sample_column
  if (!is.null(qc_col) && qc_col %in% colnames(meta)) {
    is_qc <- meta[[qc_col]] == "QC"
    if (sum(is_qc) >= 3) {
      raw_lin <- if (L$scale == "log2") 2^x else x
      rsd <- apply(raw_lin[, is_qc, drop = FALSE], 1,
                   function(r) stats::sd(r, na.rm = TRUE) / mean(r, na.rm = TRUE))
      keep_rsd <- !is.na(rsd) & rsd <= as.numeric(opt(L, "qc_rsd_max", 0.3))
      log <- c(log, sprintf("QC-RSD filter: kept %d / %d features", sum(keep_rsd), length(keep_rsd)))
      x <- x[keep_rsd, !is_qc, drop = FALSE]; meta <- meta[!is_qc, , drop = FALSE]
    }
  }

  if (L$scale == "counts") {
    if (!requireNamespace("edgeR", quietly = TRUE)) stop("Install edgeR.")
    grp <- interaction(meta$group, meta$visit, drop = TRUE)
    dge <- edgeR::DGEList(counts = round(x))
    keep <- edgeR::filterByExpr(dge, group = grp)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    dge <- edgeR::calcNormFactors(dge, method = if (normalise == "tmm") "TMM" else "none")
    values <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
    log <- c(log, sprintf("filterByExpr (group x visit): kept %d / %d genes", sum(keep), length(keep)),
             sprintf("Normalisation: %s; values = log2-CPM (prior.count = 2)", normalise))
    out <- SummarizedExperiment(
      assays = list(values = values, counts = dge$counts),
      colData = S4Vectors::DataFrame(meta, row.names = colnames(values)),
      rowData = rowData(se)[rownames(values), , drop = FALSE])
    S4Vectors::metadata(out)$norm_factors <- dge$samples$norm.factors
  } else {
    keep <- visit_aware_keep(x, meta$visit, missing_max)
    x <- x[keep, , drop = FALSE]
    log <- c(log, sprintf("Visit-aware missingness filter (<= %.0f%% missing in >= 1 visit): kept %d / %d",
                          100 * missing_max, sum(keep), length(keep)))
    if (L$scale == "raw_intensity") {
      if (any(x <= 0, na.rm = TRUE)) x[!is.na(x) & x <= 0] <- NA
      if (impute == "halfmin") x <- impute_halfmin(x, logged = FALSE)
      if (normalise == "pqn") x <- pqn(x)
      if (normalise == "median") x <- sweep(x, 2, apply(x, 2, stats::median, na.rm = TRUE), "/") *
                                       stats::median(x, na.rm = TRUE)
      values <- log2(x)
      if (impute == "minprob") values <- impute_minprob(values)
    } else {                                   # already log2
      values <- x
      if (normalise == "median") {
        values <- sweep(values, 2, apply(values, 2, stats::median, na.rm = TRUE)) +
          stats::median(values, na.rm = TRUE)
      }
      if (impute == "halfmin") values <- impute_halfmin(values, logged = TRUE)
      if (impute == "minprob") values <- impute_minprob(values)
    }
    log <- c(log, sprintf("Imputation: %s; normalisation: %s; values on log2 scale", impute, normalise),
             sprintf("Remaining NA: %d", sum(is.na(values))))
    out <- SummarizedExperiment(
      assays = list(values = values),
      colData = S4Vectors::DataFrame(meta, row.names = colnames(values)),
      rowData = rowData(se)[rownames(values), , drop = FALSE])
  }
  S4Vectors::metadata(out)$qc_log <- log
  saveRDS(out, lomics_out(cfg, "data", paste0(layer, "_processed_se.rds")))

  ## Diagnostics
  pdf(lomics_out(cfg, "qc", paste0(layer, "_pca.pdf")), width = 7, height = 5)
  for (p in pca_plots(assay(out, "values"), meta, layer, cfg$batch_column)) print(p)
  invisible(dev.off())
  pc <- stats::prcomp(t(assay(out, "values")[apply(assay(out, "values"), 1, stats::var) > 0, ]), scale. = TRUE)
  assoc <- sapply(intersect(c("visit", "group", cfg$batch_column, cfg$covariates), colnames(meta)),
                  function(cl) sapply(1:min(5, ncol(pc$x)), function(k) {
                    f <- meta[[cl]]
                    p <- tryCatch(stats::anova(stats::lm(pc$x[, k] ~ f))[["Pr(>F)"]][1],
                                  error = function(e) NA)
                    signif(p, 2)
                  }))
  assoc <- matrix(assoc, nrow = min(5, ncol(pc$x)),
                  dimnames = list(paste0("PC", 1:min(5, ncol(pc$x))),
                                  intersect(c("visit", "group", cfg$batch_column, cfg$covariates), colnames(meta))))
  log <- c(log, "", "PC association p-values (lm PCk ~ variable):",
           utils::capture.output(print(assoc)))
  writeLines(log, lomics_out(cfg, "qc", paste0(layer, "_qc_summary.txt")))
  cat(log, sep = "\n"); cat("\n")
}
lomics_session_log(cfg, "qc")
