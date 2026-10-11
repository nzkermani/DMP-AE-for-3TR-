#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## clinical_association.R — clinical endpoints and feature–clinical links
## for small longitudinal cohorts (one script, four modes).
##
## Usage (from the project root):
##   S=.claude/skills/stats-clinical-association/scripts/clinical_association.R
##   Rscript $S --mode=table1 --visit=BL [--vars=age,sex,...] [--by=]
##       Cohort description. Numeric: mean (SD) if Shapiro-Wilk p > 0.05,
##       else median [IQR]; categorical: n (%). Optional --by gives columns
##       per group (no tests).
##   Rscript $S --mode=change --vars=blood_eos,feno,... --from=BL --to=M4 [--show_visits=BL,M1,M4]
##       Within-subject change, complete pairs. Shapiro-Wilk on the paired
##       differences decides: paired t-test (mean difference, 95% CI) or
##       Wilcoxon signed-rank (Hodges-Lehmann median difference, 95% CI).
##       Per-visit summaries and a boxplot panel per variable.
##   Rscript $S --mode=correlation --layer=breath_voc --features=f1,f2|--panel=path/panel.csv \
##       --vars=acq5,blood_eos,... --visit=BL [--boot=2000]
##       Kendall's tau (normal approximation: z and p as in stats::cor.test
##       with exact = FALSE), percentile bootstrap 95% CI for tau.
##   Rscript $S --mode=exposure --layer=breath_voc --features=...|--panel=... \
##       --exposures=ocs_maintenance,lama,ltra[,ics_dose] --visit=BL
##       Feature levels by medication exposure: Mann-Whitney U (2 levels) or
##       Kruskal-Wallis (> 2). Numeric exposures are split at the median.
##
## p-values are NOT adjusted unless --adjust=BH (default none: the
## exploratory convention of the source paper; say so when reporting).
## Outputs: <output_dir>/clinical/<mode>_<tag>/  (csv + summary.txt + plots)
## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", mode = "table1", visit = NULL, vars = NULL,
                      by = NULL, from = NULL, to = NULL, show_visits = NULL, layer = NULL,
                      features = NULL, panel = NULL, exposures = NULL, boot = "2000",
                      adjust = "none", tag = NULL, seed = NULL))
cfg <- lomics_read_config(a$config)
meta_all <- lomics_read_metadata(cfg); meta_all$visit <- as.character(meta_all$visit)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
split_arg <- function(x) if (is.null(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]
seed <- if (is.null(a$seed)) (if (is.null(cfg$seed)) 1 else cfg$seed) else as.integer(a$seed)
id_cols <- c("sample_id", "subject_id", "visit", "time", "group", cfg$columns$sample,
             cfg$columns$subject, cfg$columns$visit)
tag <- if (is.null(a$tag)) switch(a$mode, change = paste0(a$from, "to", a$to), visit) else a$tag
out_dir <- dirname(lomics_out(cfg, "clinical", paste0(a$mode, "_", tag), "x"))
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
is_normal <- function(x) { x <- stats::na.omit(x); length(x) >= 3 && length(unique(x)) > 2 && stats::shapiro.test(x)$p.value > 0.05 }
lines <- character()

features_from_args <- function() {
  f <- split_arg(a$features)
  if (!is.null(a$panel)) f <- c(f, utils::read.csv(a$panel, stringsAsFactors = FALSE)$feature)
  if (!length(f)) stop("Give --features=... or --panel=<panel.csv>", call. = FALSE)
  unique(f)
}
feature_values <- function(layer, feats, v) {
  se <- lomics_load_layer(cfg, layer, "processed")
  x <- assay(se, "values"); m <- as.data.frame(colData(se)); m$visit <- as.character(m$visit)
  miss <- setdiff(feats, rownames(x)); if (length(miss)) stop("Features not found: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- m$visit == v
  rd <- as.data.frame(rowData(se))
  lab <- if ("compound" %in% colnames(rd)) ifelse(rd[feats, "compound"] %in% c(NA, "unidentified"), feats, rd[feats, "compound"]) else feats
  list(X = t(x[feats, keep, drop = FALSE]), meta = m[keep, ], label = stats::setNames(lab, feats))
}

## ---------------------------------------------------------------- table1
if (a$mode == "table1") {
  m <- meta_all[meta_all$visit == visit, ]
  vars <- if (is.null(a$vars)) setdiff(colnames(m), id_cols) else split_arg(a$vars)
  groups <- if (is.null(a$by)) list(all = m) else split(m, m[[a$by]])
  rows <- list()
  for (v in vars) {
    x <- m[[v]]
    if (is.numeric(x) && length(unique(stats::na.omit(x))) > 2) {
      nm <- is_normal(x)
      cells <- vapply(groups, function(g) { z <- g[[v]]
        if (nm) sprintf("%s (%s)", fmt(mean(z, na.rm = TRUE)), fmt(stats::sd(z, na.rm = TRUE)))
        else { q <- stats::quantile(z, c(.5, .25, .75), na.rm = TRUE); sprintf("%s [%s, %s]", fmt(q[1]), fmt(q[2]), fmt(q[3])) } }, "")
      rows[[length(rows) + 1]] <- data.frame(variable = v, level = "", summary = if (nm) "mean (SD)" else "median [IQR]",
                                             t(cells), n_missing = sum(is.na(x)), check.names = FALSE)
    } else {
      for (lv in sort(unique(stats::na.omit(as.character(x))))) {
        cells <- vapply(groups, function(g) { z <- as.character(g[[v]]); k <- sum(z == lv, na.rm = TRUE)
          sprintf("%d (%.0f%%)", k, 100 * k / sum(!is.na(z))) }, "")
        rows[[length(rows) + 1]] <- data.frame(variable = v, level = lv, summary = "n (%)", t(cells),
                                               n_missing = sum(is.na(x)), check.names = FALSE)
      }
    }
  }
  tab <- do.call(rbind, rows)
  colnames(tab)[4:(3 + length(groups))] <- sprintf("%s (n = %d)", names(groups), vapply(groups, nrow, 1L))
  utils::write.csv(tab, file.path(out_dir, "table1.csv"), row.names = FALSE)
  lines <- c(sprintf("Table 1 at visit %s (n = %d)%s; normality by Shapiro-Wilk (p > 0.05).", visit, nrow(m),
                     if (is.null(a$by)) "" else paste(" by", a$by)),
             utils::capture.output(print(tab, row.names = FALSE)))
}

## ---------------------------------------------------------------- change
if (a$mode == "change") {
  vars <- split_arg(a$vars); from <- a$from; to <- a$to
  sv <- if (is.null(a$show_visits)) unique(c(from, to)) else split_arg(a$show_visits)
  res <- list(); perv <- list()
  for (v in vars) {
    w <- stats::reshape(meta_all[meta_all$visit %in% sv, c("subject_id", "visit", v)],
                        idvar = "subject_id", timevar = "visit", direction = "wide")
    cf <- paste0(v, ".", from); ct <- paste0(v, ".", to)
    d <- w[[ct]] - w[[cf]]; ok <- !is.na(d); d <- d[ok]
    nm <- is_normal(d)
    if (nm) {
      tt <- stats::t.test(w[[ct]][ok], w[[cf]][ok], paired = TRUE)
      est <- unname(tt$estimate); ci <- tt$conf.int; p <- tt$p.value; test <- "paired t-test"; lab <- "mean difference"
    } else {
      wt <- suppressWarnings(stats::wilcox.test(w[[ct]][ok], w[[cf]][ok], paired = TRUE, conf.int = TRUE, exact = FALSE))
      est <- unname(wt$estimate); ci <- wt$conf.int; p <- wt$p.value; test <- "Wilcoxon signed-rank"; lab <- "median difference (Hodges-Lehmann)"
    }
    summ <- vapply(sv, function(s) { z <- w[[paste0(v, ".", s)]][ok]
      if (all(is.na(z))) return("NA")
      if (nm) sprintf("%s (%s)", fmt(mean(z, na.rm = TRUE)), fmt(stats::sd(z, na.rm = TRUE)))
      else { q <- stats::quantile(z, c(.5, .25, .75), na.rm = TRUE); sprintf("%s (%s-%s)", fmt(q[1]), fmt(q[2]), fmt(q[3])) } }, "")
    res[[v]] <- data.frame(variable = v, n_pairs = length(d), shapiro_normal_diff = nm, test = test,
                           t(summ), estimate_type = lab, estimate = est, ci_low = ci[1], ci_high = ci[2],
                           median_of_differences = stats::median(d), p_value = p, check.names = FALSE)
    perv[[v]] <- w
  }
  tab <- do.call(rbind, res)
  colnames(tab)[5:(4 + length(sv))] <- paste0("summary_", sv)
  if (a$adjust != "none") tab$p_adj <- stats::p.adjust(tab$p_value, a$adjust)
  utils::write.csv(tab, file.path(out_dir, "change.csv"), row.names = FALSE)
  draw <- function() {
    k <- length(vars); op <- graphics::par(mfrow = c(ceiling(k / 3), min(3, k)), mar = c(3, 4, 3, 1))
    for (v in vars) {
      mm <- meta_all[meta_all$visit %in% sv, ]; mm$visit <- factor(mm$visit, levels = sv)
      graphics::boxplot(mm[[v]] ~ mm$visit, main = sprintf("%s (n = %d pairs)", v, tab[v, "n_pairs"]),
                        xlab = "", ylab = v, col = "grey92", outline = FALSE, ylim = range(mm[[v]], na.rm = TRUE))
      graphics::points(jitter(as.integer(mm$visit), 0.4), mm[[v]], pch = 16, cex = 0.5, col = "grey35")
    }
    graphics::par(op)
  }
  grDevices::pdf(file.path(out_dir, "clinical_boxplots.pdf"), width = 10, height = 3.3 * ceiling(length(vars) / 3)); draw(); invisible(grDevices::dev.off())
  grDevices::png(file.path(out_dir, "clinical_boxplots.png"), width = 1500, height = 500 * ceiling(length(vars) / 3), res = 150); draw(); invisible(grDevices::dev.off())
  lines <- c(sprintf("Within-subject change %s -> %s (complete pairs; normality of differences by Shapiro-Wilk).", from, to),
             utils::capture.output(print(tab, row.names = FALSE, digits = 3)))
}

## ---------------------------------------------------------------- correlation
if (a$mode == "correlation") {
  layer <- if (is.null(a$layer)) cfg$layers$name[1] else a$layer
  fv <- feature_values(layer, features_from_args(), visit); vars <- split_arg(a$vars)
  B <- as.integer(a$boot); set.seed(seed)
  res <- list()
  for (f in colnames(fv$X)) for (v in vars) {
    x <- fv$X[, f]; z <- fv$meta[[v]]; ok <- !is.na(x) & !is.na(z)
    ct <- suppressWarnings(stats::cor.test(x[ok], z[ok], method = "kendall", exact = FALSE))
    bt <- if (B > 0) replicate(B, { i <- sample(which(ok), replace = TRUE); suppressWarnings(stats::cor(x[i], z[i], method = "kendall")) }) else NA
    res[[length(res) + 1]] <- data.frame(feature = f, label = fv$label[[f]], clinical = v, n = sum(ok),
      tau = unname(ct$estimate), z = unname(ct$statistic), p_value = ct$p.value,
      tau_ci_low = if (B > 0) stats::quantile(bt, 0.025, na.rm = TRUE) else NA,
      tau_ci_high = if (B > 0) stats::quantile(bt, 0.975, na.rm = TRUE) else NA, row.names = NULL)
  }
  tab <- do.call(rbind, res)
  if (a$adjust != "none") tab$p_adj <- stats::p.adjust(tab$p_value, a$adjust)
  utils::write.csv(tab, file.path(out_dir, "kendall.csv"), row.names = FALSE)
  wide <- stats::reshape(data.frame(label = tab$label, clinical = tab$clinical,
                                    cell = sprintf("tau %.2f, z = %.2f, p = %.3g", tab$tau, tab$z, tab$p_value)),
                         idvar = "clinical", timevar = "label", direction = "wide")
  colnames(wide) <- sub("^cell\\.", "", colnames(wide))
  utils::write.csv(wide, file.path(out_dir, "kendall_table.csv"), row.names = FALSE)
  lines <- c(sprintf("Kendall rank correlation at %s (layer %s); z = normal-approximation statistic; tau CI = percentile bootstrap (B = %d, seed %d); p %s.",
                     visit, layer, B, seed, if (a$adjust == "none") "unadjusted" else paste("adjusted:", a$adjust)),
             utils::capture.output(print(tab, row.names = FALSE, digits = 3)))
}

## ---------------------------------------------------------------- exposure
if (a$mode == "exposure") {
  layer <- if (is.null(a$layer)) cfg$layers$name[1] else a$layer
  fv <- feature_values(layer, features_from_args(), visit); ex <- split_arg(a$exposures)
  res <- list()
  for (e in ex) {
    g <- fv$meta[[e]]
    if (is.numeric(g) && length(unique(stats::na.omit(g))) > 2) {
      g <- ifelse(g > stats::median(g, na.rm = TRUE), paste0("> median (", stats::median(g, na.rm = TRUE), ")"), "<= median")
    }
    g <- factor(g)
    for (f in colnames(fv$X)) {
      x <- fv$X[, f]; ok <- !is.na(g)
      if (nlevels(droplevels(g[ok])) < 2) next
      p <- if (nlevels(droplevels(g[ok])) == 2) suppressWarnings(stats::wilcox.test(x[ok] ~ droplevels(g[ok]), exact = FALSE)$p.value)
           else stats::kruskal.test(x[ok] ~ droplevels(g[ok]))$p.value
      med <- tapply(x[ok], droplevels(g[ok]), stats::median)
      res[[length(res) + 1]] <- data.frame(exposure = e, feature = f, label = fv$label[[f]],
        groups = paste(sprintf("%s: n=%d, median %.2f", names(med), as.vector(table(droplevels(g[ok]))), med), collapse = "; "),
        test = if (nlevels(droplevels(g[ok])) == 2) "Mann-Whitney U" else "Kruskal-Wallis", p_value = p)
    }
  }
  tab <- do.call(rbind, res)
  if (a$adjust != "none") tab$p_adj <- stats::p.adjust(tab$p_value, a$adjust)
  utils::write.csv(tab, file.path(out_dir, "exposure_tests.csv"), row.names = FALSE)
  lines <- c(sprintf("Feature levels by exposure at %s (layer %s).", visit, layer),
             utils::capture.output(print(tab, row.names = FALSE, digits = 3)))
}

if (!length(lines)) stop("Unknown --mode: ", a$mode, " (table1 | change | correlation | exposure)", call. = FALSE)
writeLines(lines, file.path(out_dir, "summary.txt"))
lomics_session_log(cfg, paste0("clinical_", a$mode, "_", tag))
cat(lines, sep = "\n")
