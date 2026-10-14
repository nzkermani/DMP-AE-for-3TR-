#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## breath_preprocess.R — exhaled-breath GC-MS (breathomics) feature table ->
## processed layer: contaminant removal, optional zero filter, log2, batch
## correction with limma::removeBatchEffect (visit preserved), PCA before/after.
##
## Input is the aligned feature table (e.g. XCMS peak table: fragments x
## samples, raw intensities) already loaded by prep-study-setup
## (results/data/<layer>_raw_se.rds). Raw chromatogram processing (XCMS
## denoising, peak picking, alignment) is NOT done here.
##
## Usage (from the project root):
##   Rscript .claude/skills/prep-breathomics-gcms/scripts/breath_preprocess.R \
##     --config=study_config.json --layer=breath_voc \
##     [--annotation=path.csv --contaminant_col=contaminant] \
##     [--contaminant_regex="siloxane|phthalate"] [--zero_max=1] \
##     [--batch=gcms_batch|none] [--preserve=visit] [--force_batch=FALSE]
##
## Outputs (results/ = output_dir in the config):
##   data/<layer>_processed_se.rds   assay "values" = log2 (batch-corrected),
##                                   assay "log2_uncorrected"; rowData has
##                                   the annotation
##   qc/<layer>_breath_qc.txt        counts at each step, batch x visit table,
##                                   PC1-5 vs batch / visit p-values before/after
##   qc/<layer>_breath_pca.pdf/.png  PCA before / after batch correction
## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", layer = NULL, annotation = NULL,
                      contaminant_col = "contaminant", contaminant_regex = NULL,
                      compound_col = "compound", zero_max = "1", batch = NULL,
                      preserve = "visit", force_batch = "FALSE"))
cfg <- lomics_read_config(a$config)
layer <- if (is.null(a$layer)) cfg$layers$name[1] else a$layer
Lcfg <- cfg$layers[cfg$layers$name == layer, , drop = FALSE]
if (is.null(a$annotation) && "annotation" %in% names(Lcfg)) a$annotation <- Lcfg$annotation
batch_col <- if (is.null(a$batch)) cfg$batch_column else a$batch
if (identical(batch_col, "none")) batch_col <- NULL

se <- lomics_load_layer(cfg, layer, "raw")
x <- assay(se, "raw"); meta <- as.data.frame(colData(se))
log_lines <- c(sprintf("Layer: %s | samples: %d | fragments in: %d", layer, ncol(x), nrow(x)))

## 1. annotation + contaminant removal ---------------------------------------
annot <- NULL
if (!is.null(a$annotation)) {
  annot <- utils::read.csv(a$annotation, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(annot) <- annot[[1]]
  miss <- setdiff(rownames(x), rownames(annot))
  if (length(miss)) log_lines <- c(log_lines, sprintf("%d fragments have no annotation row", length(miss)))
}
is_contam <- rep(FALSE, nrow(x)); names(is_contam) <- rownames(x)
if (!is.null(annot) && a$contaminant_col %in% colnames(annot)) {
  v <- annot[rownames(x), a$contaminant_col]
  is_contam <- is_contam | (!is.na(v) & as.logical(v))
}
if (!is.null(a$contaminant_regex) && !is.null(annot) && a$compound_col %in% colnames(annot)) {
  is_contam <- is_contam | grepl(a$contaminant_regex, annot[rownames(x), a$compound_col], ignore.case = TRUE)
}
x <- x[!is_contam, , drop = FALSE]
log_lines <- c(log_lines, sprintf("Known contaminants removed: %d -> %d fragments remain",
                                  sum(is_contam), nrow(x)))

## 2. optional zero / missing filter (visit-aware: kept if OK in any visit) ---
zero_max <- as.numeric(a$zero_max)
if (zero_max < 1) {
  z <- is.na(x) | x <= 0
  ok_any <- Reduce(`|`, lapply(split(seq_len(ncol(x)), meta$visit), function(ix)
    rowMeans(z[, ix, drop = FALSE]) <= zero_max))
  log_lines <- c(log_lines, sprintf("Zero/missing filter (<= %.0f%% in at least one visit): %d removed",
                                    100 * zero_max, sum(!ok_any)))
  x <- x[ok_any, , drop = FALSE]
} else log_lines <- c(log_lines, "Zero/missing filter: off (--zero_max=1)")

## 3. log2 with per-fragment half-minimum for zeros ---------------------------
n_zero <- sum(x <= 0 | is.na(x))
lx <- t(apply(x, 1, function(r) {
  pos <- r[!is.na(r) & r > 0]
  r[is.na(r) | r <= 0] <- if (length(pos)) min(pos) / 2 else 1
  log2(r)
}))
dimnames(lx) <- dimnames(x)
log_lines <- c(log_lines, sprintf("log2 transform; %d zero/NA values (%.2f%%) set to half the fragment minimum",
                                  n_zero, 100 * n_zero / length(x)))

## 4. PCA association helper ------------------------------------------------
pc_assoc <- function(m, vars) {
  pc <- stats::prcomp(t(m), scale. = TRUE)$x[, 1:5]
  out <- sapply(vars, function(v) apply(pc, 2, function(s) {
    g <- meta[[v]]
    if (is.numeric(g)) stats::cor.test(s, g)$p.value else stats::kruskal.test(s ~ factor(g))$p.value
  }))
  signif(out, 2)
}
assoc_vars <- unique(c(batch_col, a$preserve, "group", if ("site" %in% colnames(meta)) "site"))
assoc_vars <- intersect(assoc_vars, colnames(meta))

## 5. batch correction ---------------------------------------------------------
corrected <- lx
if (!is.null(batch_col)) {
  if (!batch_col %in% colnames(meta)) stop("Batch column not in metadata: ", batch_col, call. = FALSE)
  b <- factor(meta[[batch_col]]); pv <- factor(meta[[a$preserve]])
  tab <- table(b, pv)
  log_lines <- c(log_lines, "", sprintf("Batch (%s) x %s:", batch_col, a$preserve),
                 utils::capture.output(print(tab)))
  nested <- all(rowSums(tab > 0) == 1) || all(colSums(tab > 0) == 1)
  cv <- suppressWarnings(sqrt(stats::chisq.test(tab)$statistic / (sum(tab) * (min(dim(tab)) - 1))))
  log_lines <- c(log_lines, sprintf("Cramer's V (batch vs %s) = %.2f", a$preserve, cv))
  if (nested && !as.logical(a$force_batch)) {
    stop("Batch is confounded with ", a$preserve, " (each batch holds a single level or vice versa). ",
         "Batch correction would remove the effect of interest. Model batch as a covariate instead, ",
         "or pass --force_batch=TRUE if you are sure.", call. = FALSE)
  }
  if (cv > 0.5) log_lines <- c(log_lines, "WARNING: batch strongly associated with the preserved factor; interpret with care.")
  design <- stats::model.matrix(~ pv)
  corrected <- limma::removeBatchEffect(lx, batch = b, design = design)
  log_lines <- c(log_lines, sprintf("limma::removeBatchEffect(batch = %s, design = ~ %s)", batch_col, a$preserve))
  if ("site" %in% colnames(meta)) {
    st <- table(b, meta$site)
    if (all(rowSums(st > 0) == 1)) log_lines <- c(log_lines,
      "NOTE: each batch belongs to one site, so batch correction also removes between-site differences.")
  }
} else log_lines <- c(log_lines, "Batch correction: none")

log_lines <- c(log_lines, "", "PC1-5 association p-values BEFORE batch correction:",
               utils::capture.output(print(pc_assoc(lx, assoc_vars))),
               "", "PC1-5 association p-values AFTER batch correction:",
               utils::capture.output(print(pc_assoc(corrected, assoc_vars))))

## 6. plots ----------------------------------------------------------------------
col_by <- if (!is.null(batch_col)) batch_col else "visit"
draw <- function() {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  for (m in list(list(lx, "Before batch correction"), list(corrected, "After batch correction"))) {
    pc <- stats::prcomp(t(m[[1]]), scale. = TRUE)
    ve <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2), 1)
    g <- factor(meta[[col_by]]); cols <- grDevices::hcl.colors(nlevels(g), "Dark 3")
    graphics::plot(pc$x[, 1], pc$x[, 2], col = cols[g], pch = c(16, 17, 15)[as.integer(factor(meta$visit))],
                   xlab = sprintf("PC1 (%.1f%%)", ve[1]), ylab = sprintf("PC2 (%.1f%%)", ve[2]), main = m[[2]])
    graphics::legend("topright", legend = levels(g), col = cols, pch = 16, cex = 0.7, title = col_by, bty = "n")
  }
  graphics::par(op)
}
grDevices::pdf(lomics_out(cfg, "qc", paste0(layer, "_breath_pca.pdf")), width = 10, height = 5); draw(); invisible(grDevices::dev.off())
grDevices::png(lomics_out(cfg, "qc", paste0(layer, "_breath_pca.png")), width = 1500, height = 750, res = 150); draw(); invisible(grDevices::dev.off())

## 7. save --------------------------------------------------------------------------
rd <- S4Vectors::DataFrame(feature_id = rownames(corrected), row.names = rownames(corrected))
if (!is.null(annot)) {
  an <- annot[rownames(corrected), -1, drop = FALSE]
  for (cc in colnames(an)) rd[[cc]] <- an[[cc]]
}
out <- SummarizedExperiment(assays = list(values = corrected, log2_uncorrected = lx),
                            colData = colData(se)[colnames(corrected), ], rowData = rd)
S4Vectors::metadata(out)$processing <- log_lines
saveRDS(out, lomics_out(cfg, "data", paste0(layer, "_processed_se.rds")))
writeLines(log_lines, lomics_out(cfg, "qc", paste0(layer, "_breath_qc.txt")))
lomics_session_log(cfg, paste0("breath_preprocess_", layer))
cat(log_lines, sep = "\n")
