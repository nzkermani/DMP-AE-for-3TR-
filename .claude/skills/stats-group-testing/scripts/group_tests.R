#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## group_tests.R — per-feature hypothesis tests with FDR, plus a clinical
## characteristics table.
##
## Usage (from the project root):
##   Rscript .claude/skills/stats-group-testing/scripts/group_tests.R --design=<design> [options]
##
##   --design=between   compare groups at one visit (independent samples)
##       --by=group|categorical outcome column   --visit=V1 [--change_from=V1]
##       2 groups : Welch t-test + Wilcoxon rank-sum; effect = mean diff, Cohen's d
##       >2 groups: one-way ANOVA (Welch) + Kruskal-Wallis; effect = eta^2
##   --design=paired    compare two visits within subjects
##       --visits=V1,V3  [--subset_col=group --subset=Treatment]
##       paired t-test + Wilcoxon signed-rank; effect = mean paired diff, dz
##   --design=repeated  all visits within subjects (complete cases)
##       [--subset_col=group --subset=Treatment]
##       Friedman test + repeated-measures mixed model LRT (lme4)
##   --design=clinical  table of metadata variables by --by at --visit
##       (chi-square / Fisher for categorical; t / Wilcoxon / ANOVA / Kruskal for numeric)
##
##   common: --layers=proteomics|all  --fdr=0.05
##
## Outputs: results/tests/<design>_<tag>/<layer>_tests.csv, summary.txt
## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", design = "between", layers = "all",
                      by = "group", visit = NULL, change_from = NULL, visits = NULL,
                      subset_col = NULL, subset = NULL, fdr = "0.05"))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
fdr <- as.numeric(a$fdr)
by <- lomics_outcome(cfg, a$by)

subset_meta <- function(meta) {
  if (!is.null(a$subset_col)) meta[as.character(meta[[a$subset_col]]) == a$subset, , drop = FALSE] else meta
}
p_safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
bh <- function(p) stats::p.adjust(p, "BH")

## ---------------------------------------------------------------- clinical
if (a$design == "clinical") {
  se <- lomics_load_layer(cfg, layers[1], "processed")
  meta <- as.data.frame(colData(se)); meta <- meta[as.character(meta$visit) == visit, ]
  g <- factor(meta[[by]])
  vars <- setdiff(colnames(meta), c("sample_id", "subject_id", "visit", "time", by))
  rows <- lapply(vars, function(v) {
    x <- meta[[v]]
    if (is.numeric(x) && length(unique(x)) > 2) {
      summ <- tapply(x, g, function(z) sprintf("%.2f (%.2f)", mean(z, na.rm = TRUE), stats::sd(z, na.rm = TRUE)))
      p_par <- p_safe(if (nlevels(g) == 2) stats::t.test(x ~ g)$p.value else stats::oneway.test(x ~ g)$p.value)
      p_np  <- p_safe(stats::kruskal.test(x ~ g)$p.value)
      data.frame(variable = v, type = "numeric: mean (SD)", t(summ), p_parametric = p_par, p_nonparametric = p_np, check.names = FALSE)
    } else {
      tb <- table(x, g)
      summ <- apply(tb, 2, function(cl) paste(sprintf("%s: %d (%.0f%%)", rownames(tb), cl, 100 * cl / sum(cl)), collapse = "; "))
      p_chi <- p_safe(suppressWarnings(stats::chisq.test(tb)$p.value))
      p_fis <- p_safe(stats::fisher.test(tb, simulate.p.value = nrow(tb) * ncol(tb) > 4)$p.value)
      data.frame(variable = v, type = "categorical: n (%)", t(summ), p_parametric = p_chi, p_nonparametric = p_fis, check.names = FALSE)
    }
  })
  tab <- do.call(rbind, rows)
  outdir <- dirname(lomics_out(cfg, "tests", paste0("clinical_by_", by, "_", visit), "x"))
  utils::write.csv(tab, file.path(outdir, "clinical_table.csv"), row.names = FALSE)
  print(tab, row.names = FALSE)
  message("p_parametric = t/ANOVA or chi-square; p_nonparametric = Kruskal-Wallis or Fisher. Not FDR-adjusted.")
  quit(save = "no")
}

## ---------------------------------------------------------------- features
res_all <- list()
tag <- switch(a$design,
  between  = lomics_tag(by, visit, a$change_from),
  paired   = paste0(gsub(",", "vs", a$visits), if (!is.null(a$subset)) paste0("_", a$subset)),
  repeated = paste0("allvisits", if (!is.null(a$subset)) paste0("_", a$subset)),
  stop("Unknown --design: ", a$design))
outdir <- dirname(lomics_out(cfg, "tests", paste0(a$design, "_", tag), "x"))

for (layer in layers) {
  if (a$design == "between") {
    d <- lomics_xy(cfg, layer, visit, a$change_from)
    X <- d$X; g <- factor(d$meta[[by]]); ok <- !is.na(g); X <- X[ok, , drop = FALSE]; g <- droplevels(g[ok])
    if (nlevels(g) < 2) stop("Need >= 2 levels of ", by)
    if (nlevels(g) == 2) {
      l1 <- levels(g)[1]; l2 <- levels(g)[2]
      res <- do.call(rbind, lapply(colnames(X), function(f) {
        x1 <- X[g == l1, f]; x2 <- X[g == l2, f]
        sp <- sqrt(((length(x1) - 1) * stats::var(x1) + (length(x2) - 1) * stats::var(x2)) / (length(x1) + length(x2) - 2))
        data.frame(feature = f,
                   mean_diff = mean(x2) - mean(x1),
                   cohens_d = (mean(x2) - mean(x1)) / sp,
                   p_t = p_safe(stats::t.test(x2, x1)$p.value),
                   p_wilcox = p_safe(stats::wilcox.test(x2, x1, exact = FALSE)$p.value))
      }))
      res$comparison <- sprintf("%s - %s", l2, l1)
    } else {
      res <- do.call(rbind, lapply(colnames(X), function(f) {
        x <- X[, f]; ss_b <- sum(tapply(x, g, function(z) length(z) * (mean(z) - mean(x))^2))
        data.frame(feature = f, eta_sq = ss_b / sum((x - mean(x))^2),
                   p_anova = p_safe(stats::oneway.test(x ~ g)$p.value),
                   p_kruskal = p_safe(stats::kruskal.test(x ~ g)$p.value))
      }))
    }
    n_info <- paste(names(table(g)), table(g), sep = "=", collapse = ", ")
  } else {
    se <- lomics_load_layer(cfg, layer, "processed")
    x <- assay(se, "values"); meta <- subset_meta(as.data.frame(colData(se)))
    meta$visit <- as.character(meta$visit)
    if (a$design == "paired") {
      vv <- strsplit(a$visits, ",")[[1]]; stopifnot(length(vv) == 2)
      m1 <- meta[meta$visit == vv[1], ]; m2 <- meta[meta$visit == vv[2], ]
      s <- intersect(m1$subject_id, m2$subject_id)
      x1 <- x[, m1$sample_id[match(s, m1$subject_id)], drop = FALSE]
      x2 <- x[, m2$sample_id[match(s, m2$subject_id)], drop = FALSE]
      dd <- x2 - x1
      res <- data.frame(feature = rownames(x), mean_paired_diff = rowMeans(dd),
                        cohens_dz = rowMeans(dd) / apply(dd, 1, stats::sd),
                        p_paired_t = apply(dd, 1, function(z) p_safe(stats::t.test(z)$p.value)),
                        p_signed_rank = apply(dd, 1, function(z) p_safe(stats::wilcox.test(z, exact = FALSE)$p.value)))
      res$comparison <- sprintf("%s - %s", vv[2], vv[1])
      n_info <- sprintf("%d subjects with both visits", length(s))
    } else {
      tabv <- table(meta$subject_id)
      complete <- names(tabv)[tabv == length(cfg$visits)]
      mc <- meta[meta$subject_id %in% complete, ]
      mc <- mc[order(mc$subject_id, factor(mc$visit, cfg$visits)), ]
      xm <- x[, mc$sample_id, drop = FALSE]
      res <- do.call(rbind, lapply(rownames(xm), function(f) {
        y <- xm[f, ]
        mat <- matrix(y, ncol = length(cfg$visits), byrow = TRUE)
        df <- data.frame(y = y, visit = factor(mc$visit, cfg$visits), subject_id = mc$subject_id)
        p_lrt <- p_safe({
          f1 <- lme4::lmer(y ~ visit + (1 | subject_id), df, REML = FALSE)
          f0 <- lme4::lmer(y ~ 1 + (1 | subject_id), df, REML = FALSE)
          stats::anova(f0, f1)$`Pr(>Chisq)`[2]
        })
        data.frame(feature = f, max_visit_mean_diff = diff(range(colMeans(mat))),
                   p_friedman = p_safe(stats::friedman.test(mat)$p.value), p_mixed_lrt = p_lrt)
      }))
      n_info <- sprintf("%d subjects with all %d visits", length(complete), length(cfg$visits))
    }
  }
  pcols <- grep("^p_", names(res), value = TRUE)
  for (p in pcols) res[[sub("^p_", "fdr_", p)]] <- bh(res[[p]])
  res <- res[order(res[[pcols[1]]]), ]
  res$layer <- layer
  utils::write.csv(res, file.path(outdir, paste0(layer, "_tests.csv")), row.names = FALSE)
  res_all[[layer]] <- data.frame(layer = layer, n = n_info, features = nrow(res),
                                 t(sapply(sub("^p_", "fdr_", pcols), function(fc) sum(res[[fc]] < fdr, na.rm = TRUE))),
                                 check.names = FALSE)
}
smry <- do.call(rbind, lapply(res_all, function(z) { names(z) <- sub("^fdr_", paste0("n_fdr<", fdr, "_"), names(z)); z }))
writeLines(c(sprintf("design: %s | %s", a$design, tag), utils::capture.output(print(smry, row.names = FALSE))),
           file.path(outdir, "summary.txt"))
print(smry, row.names = FALSE)
message("Outputs in ", outdir)
lomics_session_log(cfg, paste0("tests_", a$design))
