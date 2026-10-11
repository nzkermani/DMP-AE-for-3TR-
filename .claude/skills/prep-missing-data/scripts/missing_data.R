#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## missing_data.R — diagnose, benchmark and handle missing data in omic layers
## and clinical variables of a longitudinal study.
##
## Usage (from the project root):
##   S=.claude/skills/prep-missing-data/scripts/missing_data.R
##   Rscript $S --mode=diagnose  --layer=proteomics [--zero_as_na=FALSE] [--clinical=feno,ige,...]
##   Rscript $S --mode=benchmark --layer=proteomics [--methods=halfmin,minprob,knn,mean[,bpca]] \
##        [--frac=0.05] [--reps=3] [--max_missing=0.5]
##   Rscript $S --mode=impute    --layer=proteomics --method=knn [--max_missing=0.5] [--normalise=none|median]
##        -> new layer "<layer>_imp_<method>" usable by every skill (--layer=...)
##   Rscript $S --mode=mice --visit=BL --formula="fev1_pp ~ age + sex + feno + ige" [--m=20] [--log=feno,ige]
##        -> multiple imputation (mice, PMM) + Rubin-pooled lm vs complete-case
##
## Omic values are taken from results/data/<layer>_raw_se.rds (prep-study-setup);
## raw_intensity / gcms_intensity layers are log2-transformed, counts layers
## are refused (zeros in counts are not missing). --zero_as_na treats zeros
## as missing (below detection), e.g. for GC-MS / untargeted MS tables.
##
## Outputs: <output_dir>/missing/<mode>_<tag>/ (+ data/<layer>_imp_<method>_processed_se.rds)
## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lomics_utils.R"); source("R/lomics_missing.R")
a <- lomics_args(list(config = "study_config.json", mode = "diagnose", layer = NULL, zero_as_na = "FALSE",
                      clinical = NULL, methods = "halfmin,minprob,knn,mean", frac = "0.05", reps = "3",
                      max_missing = "0.5", method = "halfmin", normalise = "none", visit = NULL, formula = NULL,
                      m = "20", log = NULL, tag = NULL, seed = NULL))
cfg <- lomics_read_config(a$config)
seed <- if (is.null(a$seed)) (if (is.null(cfg$seed)) 1 else cfg$seed) else as.integer(a$seed)
split_arg <- function(x) if (is.null(x) || !nzchar(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]
layer <- if (is.null(a$layer)) cfg$layers$name[1] else a$layer
tag <- if (is.null(a$tag)) (if (a$mode == "mice") (if (is.null(a$visit)) cfg$default_visit else a$visit) else layer) else a$tag
out_dir <- dirname(lomics_out(cfg, "missing", paste0(a$mode, "_", tag), "x"))
unlink(list.files(out_dir, full.names = TRUE))
lines <- character()
save_plot <- function(fun, name, w = 10, h = 4) {
  grDevices::pdf(file.path(out_dir, paste0(name, ".pdf")), width = w, height = h); fun(); invisible(grDevices::dev.off())
  grDevices::png(file.path(out_dir, paste0(name, ".png")), width = w * 150, height = h * 150, res = 150); fun(); invisible(grDevices::dev.off())
}

load_log_matrix <- function() {
  se <- lomics_load_layer(cfg, layer, "raw")
  x <- assay(se, "raw"); sc <- S4Vectors::metadata(se)$scale
  if (identical(sc, "counts")) stop(layer, " is a counts layer: zeros are observations, not missing values.", call. = FALSE)
  if (as.logical(a$zero_as_na)) x[!is.na(x) & x <= 0] <- NA
  if (sc %in% c("raw_intensity", "gcms_intensity")) { x[!is.na(x) & x <= 0] <- NA; x <- log2(x) }
  list(x = x, meta = as.data.frame(colData(se)), scale = sc, se = se)
}

## ================================================================ diagnose
if (a$mode == "diagnose") {
  L <- load_log_matrix(); x <- L$x; meta <- L$meta; meta$visit <- as.character(meta$visit)
  M <- is.na(x)
  fmiss <- rowMeans(M); smiss <- colMeans(M)
  by_visit <- sapply(split(seq_len(ncol(x)), factor(meta$visit, cfg$visits)), function(ix) if (length(ix)) rowMeans(M[, ix, drop = FALSE]) else rep(NA, nrow(x)))
  feat <- data.frame(feature = rownames(x), pct_missing = round(100 * fmiss, 1),
                     mean_observed = rowMeans(x, na.rm = TRUE), round(100 * by_visit, 1), check.names = FALSE)
  colnames(feat)[-(1:3)] <- paste0("pct_missing_", colnames(by_visit))
  utils::write.csv(feat[order(-feat$pct_missing), ], file.path(out_dir, "features.csv"), row.names = FALSE)
  samp <- data.frame(sample_id = colnames(x), subject_id = meta$subject_id, visit = meta$visit, group = meta$group,
                     pct_missing = round(100 * smiss, 1))
  utils::write.csv(samp, file.path(out_dir, "samples.csv"), row.names = FALSE)
  ## MNAR signature: features with lower intensity have more missing values
  sub <- fmiss > 0 & fmiss < 1
  mnar_rho <- if (sum(sub) > 10) suppressWarnings(stats::cor(feat$mean_observed[sub], fmiss[sub], method = "spearman")) else NA
  ## does sample missingness depend on design variables?
  dv <- intersect(c("visit", "group", cfg$batch_column, "site"), colnames(meta))
  dep <- sapply(dv, function(v) stats::kruskal.test(smiss ~ factor(meta[[v]]))$p.value)
  ## subject follow-up pattern (all visits in the metadata)
  md <- lomics_read_metadata(cfg); md$visit <- as.character(md$visit)
  has <- sapply(cfg$visits, function(v) tapply(md$visit == v, md$subject_id, any))
  pattern <- apply(has, 1, function(h) {
    if (all(h)) return("complete")
    last <- max(which(h)); if (all(h[seq_len(last)])) "monotone dropout" else "intermittent" })
  pat_tab <- table(pattern)
  lines <- c(sprintf("Missing-data diagnosis: layer %s (%s scale%s)", layer, L$scale, if (as.logical(a$zero_as_na)) ", zeros treated as missing" else ""),
    sprintf("Overall missing: %.2f%% of %d values | features with any missing: %d of %d | > 50%% missing: %d",
            100 * mean(M), length(M), sum(fmiss > 0), nrow(x), sum(fmiss > 0.5)),
    sprintf("Samples: median %.1f%% missing (max %.1f%%, %s)", 100 * stats::median(smiss), 100 * max(smiss), names(which.max(smiss))),
    sprintf("Missing by visit (mean %% per feature): %s", paste(sprintf("%s %.1f", colnames(by_visit), 100 * colMeans(by_visit, na.rm = TRUE)), collapse = "; ")),
    sprintf("Intensity vs missingness (Spearman, features with partial missingness): %.2f", mnar_rho),
    if (!is.na(mnar_rho) && mnar_rho < -0.3) "  -> strong negative: low-abundance values go missing = left-censored (MNAR). halfmin/minprob underestimate, kNN/BPCA (MAR) overestimate; compare error AND bias with --mode=benchmark." else
      "  -> weak: missingness not driven by intensity; MAR methods (knn, bpca) are reasonable. Check with --mode=benchmark.",
    "Sample missingness vs design (Kruskal-Wallis p):", paste(sprintf("  %s: %.3g", names(dep), dep), collapse = "\n"),
    if (any(dep < 0.05)) "  -> missingness differs by a design variable: report it, check that imputed values do not create differences between those levels, and keep that variable in the models." else "",
    "Subject follow-up pattern (whole visits missing):", utils::capture.output(print(pat_tab)),
    "  Whole missing visits are handled by mixed models (stats-mixed-models) under MAR; do not impute whole visits.")
  if (!is.null(a$clinical)) {
    cv <- split_arg(a$clinical)
    ct <- sapply(cv, function(v) tapply(is.na(md[[v]]), factor(md$visit, cfg$visits), sum))
    lines <- c(lines, "", "Clinical variables: number missing by visit", utils::capture.output(print(ct)))
    utils::write.csv(data.frame(visit = rownames(ct), ct, check.names = FALSE), file.path(out_dir, "clinical_missing.csv"), row.names = FALSE)
  }
  save_plot(function() {
    op <- graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    graphics::plot(feat$mean_observed, 100 * fmiss, pch = 16, cex = 0.4, col = grDevices::adjustcolor("grey20", 0.4),
                   xlab = "mean observed log2 value", ylab = "% missing", main = sprintf("Intensity vs missing (rho %.2f)", mnar_rho))
    graphics::hist(100 * fmiss[fmiss > 0], breaks = 30, main = "Features with missing values", xlab = "% missing", col = "grey80")
    graphics::boxplot(100 * smiss ~ factor(meta$visit, cfg$visits), main = "Sample missingness by visit", xlab = "", ylab = "% missing", col = "grey92")
    graphics::par(op)
  }, "diagnosis", 12, 4)
}

## ================================================================ benchmark
if (a$mode == "benchmark") {
  L <- load_log_matrix(); x <- L$x
  x <- x[rowMeans(is.na(x)) <= as.numeric(a$max_missing), , drop = FALSE]
  meths <- split_arg(a$methods)
  avail <- c(halfmin = TRUE, minprob = TRUE, mean = TRUE, knn = requireNamespace("impute", quietly = TRUE),
             bpca = requireNamespace("pcaMethods", quietly = TRUE))
  skipped <- meths[!avail[meths] %in% TRUE]; meths <- meths[avail[meths] %in% TRUE]
  b <- rbind(lomics_mask_benchmark(x, meths, as.numeric(a$frac), "MCAR", reps = as.integer(a$reps), seed = seed),
             lomics_mask_benchmark(x, meths, as.numeric(a$frac), "MNAR", reps = as.integer(a$reps), seed = seed))
  utils::write.csv(b, file.path(out_dir, "benchmark_reps.csv"), row.names = FALSE)
  s <- stats::aggregate(cbind(rmse, nrmse, bias, spearman) ~ mechanism + method, b, mean)
  s <- s[order(s$mechanism, s$nrmse), ]
  utils::write.csv(s, file.path(out_dir, "benchmark_summary.csv"), row.names = FALSE)
  best <- sapply(split(s, s$mechanism), function(d) d$method[1])
  lines <- c(sprintf("Imputation benchmark: layer %s, %d features (<= %s missing), %s%% of observed values masked, %s reps, seed %d",
                     layer, nrow(x), a$max_missing, 100 * as.numeric(a$frac), a$reps, seed),
             "MCAR = masked at random; MNAR = masked among each feature's lowest 20% (left-censoring).",
             if (length(skipped)) paste("Skipped (package missing):", paste(skipped, collapse = ", ")) else "",
             utils::capture.output(print(s, row.names = FALSE, digits = 3)),
             sprintf("Lowest NRMSE: MCAR -> %s; MNAR -> %s. Choose by the mechanism --mode=diagnose suggests, and weigh the bias sign (negative = imputed too low).", best["MCAR"], best["MNAR"]))
  save_plot(function() {
    op <- graphics::par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
    for (mm in c("MCAR", "MNAR")) {
      d <- b[b$mechanism == mm, ]; d$method <- factor(d$method, levels = s$method[s$mechanism == mm])
      graphics::boxplot(nrmse ~ method, d, las = 2, xlab = "", ylab = "NRMSE = RMSE / SD of all values (lower is better)", main = mm, col = "grey90")
    }
    graphics::par(op)
  }, "benchmark", 9, 4)
}

## ================================================================ impute
if (a$mode == "impute") {
  L <- load_log_matrix(); x <- L$x
  keep <- rowMeans(is.na(x)) <= as.numeric(a$max_missing)
  mask <- is.na(x[keep, , drop = FALSE])
  xi <- lomics_impute(x[keep, , drop = FALSE], a$method, seed = seed)
  if (a$normalise == "median") xi <- sweep(xi, 2, apply(xi, 2, stats::median) - stats::median(xi))
  new_layer <- paste0(layer, "_imp_", a$method)
  out <- SummarizedExperiment(assays = list(values = xi, missing_mask = mask * 1),
                              colData = colData(L$se), rowData = rowData(L$se)[keep, , drop = FALSE])
  S4Vectors::metadata(out)$imputation <- list(method = a$method, seed = seed, max_missing = a$max_missing,
                                                normalise = a$normalise, source_layer = layer)
  saveRDS(out, lomics_out(cfg, "data", paste0(new_layer, "_processed_se.rds")))
  lines <- c(sprintf("Imputed layer written: %s (method %s, seed %d)", new_layer, a$method, seed),
             sprintf("Features kept (<= %s missing): %d of %d; values imputed: %d (%.2f%%)", a$max_missing,
                     sum(keep), length(keep), sum(mask), 100 * mean(mask)),
             sprintf("Normalisation after imputation: %s. Use --layer=%s in any skill; assay 'missing_mask' marks imputed cells.",
                     a$normalise, new_layer))
}

## ================================================================ mice
if (a$mode == "mice") {
  if (!requireNamespace("mice", quietly = TRUE)) stop("Package 'mice' is required.", call. = FALSE)
  if (is.null(a$formula)) stop("--formula is required for --mode=mice", call. = FALSE)
  md <- lomics_read_metadata(cfg); md$visit <- as.character(md$visit)
  visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
  d <- md[md$visit == visit, , drop = FALSE]
  fml <- stats::as.formula(a$formula); vars <- all.vars(fml)
  dd <- d[, vars, drop = FALSE]
  for (v in intersect(split_arg(a$log), vars)) dd[[v]] <- log(dd[[v]])
  for (v in vars) if (is.character(dd[[v]])) dd[[v]] <- factor(dd[[v]])
  miss <- colSums(is.na(dd))
  imp <- mice::mice(dd, m = as.integer(a$m), method = "pmm", seed = seed, printFlag = FALSE)
  fit <- mice::as.mira(lapply(seq_len(imp$m), function(i) stats::lm(fml, data = mice::complete(imp, i))))
  pooled <- summary(mice::pool(fit), conf.int = TRUE)
  cc <- stats::lm(fml, data = dd)
  ccs <- data.frame(term = names(stats::coef(cc)), estimate = stats::coef(cc), stats::confint(cc),
                    p.value = summary(cc)$coefficients[, 4], row.names = NULL)
  names(ccs)[3:4] <- c("2.5 %", "97.5 %")
  utils::write.csv(pooled, file.path(out_dir, "pooled_mice.csv"), row.names = FALSE)
  utils::write.csv(ccs, file.path(out_dir, "complete_case.csv"), row.names = FALSE)
  saveRDS(imp, file.path(out_dir, "mice_object.rds"))
  lines <- c(sprintf("Multiple imputation (mice PMM, m = %s, seed %d) at visit %s: %s", a$m, seed, visit, deparse(fml)),
             sprintf("n = %d; complete cases = %d; missing per variable: %s", nrow(dd), sum(stats::complete.cases(dd)),
                     paste(sprintf("%s %d", names(miss), miss), collapse = ", ")),
             if (!is.null(a$log)) paste("log-transformed:", a$log) else "",
             "", "Pooled (Rubin's rules):", utils::capture.output(print(pooled[, c("term", "estimate", "std.error", "2.5 %", "97.5 %", "p.value")], digits = 3, row.names = FALSE)),
             "", "Complete-case:", utils::capture.output(print(ccs, digits = 3, row.names = FALSE)),
             "", "Imputation model = the analysis variables only; add auxiliary variables that predict missingness for a stronger MAR assumption.")
}

if (!length(lines)) stop("Unknown --mode: ", a$mode, " (diagnose | benchmark | impute | mice)", call. = FALSE)
writeLines(lines, file.path(out_dir, "summary.txt"))
lomics_session_log(cfg, paste0("missing_", a$mode, "_", tag))
cat(lines, sep = "\n")
