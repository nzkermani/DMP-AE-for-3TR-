#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## mixed_models.R — linear mixed-effects models for longitudinal data, per omic
## feature or per clinical variable, using ALL available visits (subjects with
## missing visits are kept).
##
##   y ~ time [* group] [+ covariates] + (1 | subject)      (--random=intercept)
##   y ~ time [* group] [+ covariates] + (1 + t | subject)  (--random=slope)
##
## Usage (from the project root):
##   S=.claude/skills/stats-mixed-models/scripts/mixed_models.R
##   Rscript $S --mode=features --layer=breath_voc [--features=f1,f2|--panel=panel.csv] \
##        [--time=factor|continuous] [--by=group|none] [--covariates=age,sex] \
##        [--random=intercept|slope] [--visits=BL,M1,M4] [--report_top=20]
##   Rscript $S --mode=clinical --vars=blood_eos,feno,ige,fev1_pp [--log=blood_eos,feno,ige] ...
##
## Inference
##   * time effect and time x group interaction: likelihood-ratio tests of
##     nested models fitted by ML (lme4), BH-FDR across features/variables;
##   * group-specific time effect: LRT of time within each group level;
##   * reported estimates (REML): estimated marginal means per visit (x group)
##     and contrasts vs the first visit with 95% CI (emmeans, Kenward-Roger df
##     when pbkrtest is installed), for the variables in --report_top or the
##     given feature list;
##   * ICC (between-subject share of variance) and singular-fit flags.
##
## Outputs: <output_dir>/mixed/<mode>_<tag>/
##   lmm_results.csv  one row per feature/variable (n, LRT stats, p, FDR, ICC)
##   emmeans.csv      marginal means by visit (x group)
##   contrasts.csv    change vs first visit, per group, with CI and p
##   trajectories.(pdf|png)  observed subject lines + model means with CI
##   summary.txt
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(SummarizedExperiment); library(lme4) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", mode = "features", layer = NULL, features = NULL,
                      panel = NULL, vars = NULL, log = NULL, time = "factor", by = "none",
                      covariates = NULL, random = "intercept", visits = NULL, report_top = "20",
                      fdr = "0.05", tag = NULL))
cfg <- lomics_read_config(a$config)
split_arg <- function(x) if (is.null(x) || identical(x, "none") || !nzchar(x)) NULL else strsplit(x, ",", fixed = TRUE)[[1]]
vis <- if (is.null(a$visits)) cfg$visits else split_arg(a$visits)
by <- if (identical(a$by, "none")) NULL else lomics_outcome(cfg, a$by)
covs <- split_arg(a$covariates)
tag <- if (is.null(a$tag)) paste(c(if (a$mode == "features") (if (is.null(a$layer)) cfg$layers$name[1] else a$layer),
                                    a$time, if (!is.null(by)) paste0("by", by), a$random), collapse = "_") else a$tag
out_dir <- dirname(lomics_out(cfg, "mixed", paste0(a$mode, "_", tag), "x"))
unlink(list.files(out_dir, full.names = TRUE))

## ---------------------------------------------------------------- data
if (a$mode == "features") {
  layer <- if (is.null(a$layer)) cfg$layers$name[1] else a$layer
  se <- lomics_load_layer(cfg, layer, "processed")
  Y <- t(assay(se, "values")); meta <- as.data.frame(colData(se))
  feats <- split_arg(a$features)
  if (!is.null(a$panel)) feats <- unique(c(feats, utils::read.csv(a$panel, stringsAsFactors = FALSE)$feature))
  if (!is.null(feats)) Y <- Y[, feats, drop = FALSE]
  rd <- as.data.frame(rowData(se))
  label <- if ("compound" %in% colnames(rd)) stats::setNames(ifelse(rd$compound %in% c(NA, "unidentified"), rownames(rd), rd$compound), rownames(rd)) else stats::setNames(colnames(Y), colnames(Y))
} else {
  meta <- lomics_read_metadata(cfg)
  vars <- split_arg(a$vars); if (is.null(vars)) stop("--vars is required for --mode=clinical", call. = FALSE)
  Y <- as.matrix(meta[, vars, drop = FALSE])
  for (v in intersect(split_arg(a$log), vars)) Y[, v] <- log(pmax(Y[, v], 0) + if (any(Y[, v] <= 0, na.rm = TRUE)) 0.01 else 0)
  label <- stats::setNames(ifelse(vars %in% split_arg(a$log), paste0("log ", vars), vars), vars)
}
meta$visit <- as.character(meta$visit)
keep <- meta$visit %in% vis
Y <- Y[keep, , drop = FALSE]; meta <- meta[keep, , drop = FALSE]
meta$visit <- factor(meta$visit, levels = vis)
meta$t <- unname(cfg$time_values[as.character(meta$visit)])
meta$subject_id <- factor(meta$subject_id)
tterm <- if (a$time == "factor") "visit" else "t"
if (!is.null(by)) meta$.grp <- factor(meta[[by]])
rhs_cov <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""
re <- if (a$random == "slope") "(1 + t | subject_id)" else "(1 | subject_id)"
f_full <- stats::as.formula(paste("y ~", if (is.null(by)) tterm else paste(tterm, "* .grp"), rhs_cov, "+", re))
f_noint <- if (!is.null(by)) stats::as.formula(paste("y ~", tterm, "+ .grp", rhs_cov, "+", re)) else NULL
f_notime <- stats::as.formula(paste("y ~", if (is.null(by)) "1" else ".grp", rhs_cov, "+", re))
ctrl <- lme4::lmerControl(calc.derivs = FALSE, check.conv.singular = "ignore")

## ---------------------------------------------------------------- fits
lrt <- function(m1, m0) { s <- stats::anova(m0, m1); c(chisq = s$Chisq[2], df = s$Df[2], p = s$`Pr(>Chisq)`[2]) }
fit_one <- function(j) {
  d <- meta; d$y <- Y[, j]; d <- d[!is.na(d$y), , drop = FALSE]
  out <- data.frame(variable = colnames(Y)[j], label = label[[colnames(Y)[j]]], n_obs = nrow(d),
                    n_subjects = length(unique(d$subject_id)),
                    n_by_visit = paste(sprintf("%s=%d", vis, as.vector(table(d$visit))), collapse = ";"))
  r <- tryCatch({
    m1 <- suppressMessages(lme4::lmer(f_full, data = d, REML = FALSE, control = ctrl))
    m0 <- suppressMessages(lme4::lmer(f_notime, data = d, REML = FALSE, control = ctrl))
    tt <- lrt(m1, m0)
    vc <- as.data.frame(lme4::VarCorr(m1)); icc <- vc$vcov[1] / sum(vc$vcov[is.na(vc$var2)])
    res <- data.frame(time_chisq = tt[["chisq"]], time_df = tt[["df"]], p_time = tt[["p"]], icc = icc,
                      singular = lme4::isSingular(m1))
    if (!is.null(by)) {
      mi <- suppressMessages(lme4::lmer(f_noint, data = d, REML = FALSE, control = ctrl))
      ti <- lrt(m1, mi); res$interaction_chisq <- ti[["chisq"]]; res$interaction_df <- ti[["df"]]; res$p_interaction <- ti[["p"]]
      for (g in levels(d$.grp)) {      # time effect within each group (separate model)
        dg <- d[d$.grp == g, , drop = FALSE]
        fg1 <- stats::as.formula(paste("y ~", tterm, rhs_cov, "+", re)); fg0 <- stats::as.formula(paste("y ~ 1", rhs_cov, "+", re))
        pg <- tryCatch(lrt(suppressMessages(lme4::lmer(fg1, dg, REML = FALSE, control = ctrl)),
                           suppressMessages(lme4::lmer(fg0, dg, REML = FALSE, control = ctrl)))[["p"]], error = function(e) NA)
        res[[paste0("p_time_", g)]] <- pg
      }
    }
    res
  }, error = function(e) data.frame(error = conditionMessage(e)))
  cbind(out, r)
}
res <- do.call(rbind, lapply(seq_len(ncol(Y)), function(j) {
  r <- fit_one(j); for (cc in setdiff(c("p_time", "error"), names(r))) r[[cc]] <- NA; r
}))
res$fdr_time <- stats::p.adjust(res$p_time, "BH")
if ("p_interaction" %in% names(res)) res$fdr_interaction <- stats::p.adjust(res$p_interaction, "BH")
res <- res[order(res$p_time), ]
utils::write.csv(res, file.path(out_dir, "lmm_results.csv"), row.names = FALSE)

## ---------------------------------------------------------------- estimates
rep_vars <- if (a$mode == "clinical" || !is.null(split_arg(a$features)) || !is.null(a$panel)) res$variable else
  utils::head(res$variable, as.integer(a$report_top))
emm_rows <- list(); con_rows <- list(); fits <- list()
has_emm <- requireNamespace("emmeans", quietly = TRUE)
for (v in rep_vars) {
  d <- meta; d$y <- Y[, v]; d <- d[!is.na(d$y), , drop = FALSE]
  m <- tryCatch(suppressMessages(lme4::lmer(f_full, data = d, REML = TRUE, control = ctrl)), error = function(e) NULL)
  if (is.null(m) || !has_emm) next
  fits[[v]] <- list(d = d, m = m)
  spec <- if (a$time == "factor") { if (is.null(by)) ~ visit else ~ visit | .grp } else NULL
  if (is.null(spec)) {                                  # continuous time: slope per group
    tr <- emmeans::emtrends(m, if (is.null(by)) ~ 1 else ~ .grp, var = "t")
    s <- as.data.frame(summary(tr, infer = TRUE))
    con_rows[[v]] <- data.frame(variable = v, label = label[[v]], s, check.names = FALSE)
    next
  }
  em <- suppressMessages(emmeans::emmeans(m, spec))
  emm_rows[[v]] <- data.frame(variable = v, label = label[[v]], as.data.frame(summary(em)), check.names = FALSE)
  ct <- emmeans::contrast(em, method = "trt.vs.ctrl", ref = 1, adjust = "none")
  con_rows[[v]] <- data.frame(variable = v, label = label[[v]], as.data.frame(summary(ct, infer = TRUE)), check.names = FALSE)
}
if (length(emm_rows)) utils::write.csv(do.call(rbind, emm_rows), file.path(out_dir, "emmeans.csv"), row.names = FALSE)
if (length(con_rows)) utils::write.csv(do.call(rbind, con_rows), file.path(out_dir, "contrasts.csv"), row.names = FALSE)

## ---------------------------------------------------------------- plot
plot_vars <- utils::head(names(fits), 12)
draw <- function() {
  k <- length(plot_vars); if (!k) return(invisible())
  op <- graphics::par(mfrow = c(ceiling(k / 3), min(3, k)), mar = c(3, 4, 3, 1))
  for (v in plot_vars) {
    d <- fits[[v]]$d; x <- if (a$time == "factor") as.integer(d$visit) else d$t
    g <- if (is.null(by)) factor(rep("all", nrow(d))) else d$.grp
    cols <- grDevices::hcl.colors(max(2, nlevels(g)), "Dark 3")
    graphics::plot(x, d$y, type = "n", xaxt = "n", xlab = "", ylab = label[[v]],
                   main = sprintf("%s\nLRT time p = %.2g", label[[v]], res$p_time[res$variable == v]), cex.main = 0.85)
    graphics::axis(1, at = if (a$time == "factor") seq_along(vis) else unique(d$t), labels = if (a$time == "factor") vis else unique(d$t))
    for (s in unique(d$subject_id)) { ix <- which(d$subject_id == s); ix <- ix[order(x[ix])]
      graphics::lines(x[ix], d$y[ix], col = grDevices::adjustcolor(cols[g[ix[1]]], 0.35)) }
    if (!is.null(emm_rows[[v]])) {
      e <- emm_rows[[v]]; eg <- if (is.null(by)) factor(rep("all", nrow(e))) else factor(e$.grp, levels = levels(g))
      for (lv in levels(eg)) { ee <- e[eg == lv, ]; xx <- as.integer(factor(ee$visit, levels = vis)) + (as.integer(factor(lv, levels(eg))) - 1) * 0.08
        graphics::arrows(xx, ee$lower.CL, xx, ee$upper.CL, angle = 90, code = 3, length = 0.03, col = cols[match(lv, levels(eg))], lwd = 2)
        graphics::lines(xx, ee$emmean, col = cols[match(lv, levels(eg))], lwd = 3) }
    }
    if (v == plot_vars[1] && !is.null(by)) graphics::legend("topright", levels(g), col = cols, lwd = 3, bty = "n", cex = 0.7)
  }
  graphics::par(op)
}
h <- 3.3 * max(1, ceiling(length(plot_vars) / 3))
grDevices::pdf(file.path(out_dir, "trajectories.pdf"), width = 10, height = h); draw(); invisible(grDevices::dev.off())
grDevices::png(file.path(out_dir, "trajectories.png"), width = 1500, height = h * 150, res = 150); draw(); invisible(grDevices::dev.off())

## ---------------------------------------------------------------- summary
fdr <- as.numeric(a$fdr)
lines <- c(sprintf("Linear mixed models (%s): %s", a$mode, deparse(f_full)),
           sprintf("visits: %s | time: %s | random: %s | covariates: %s | by: %s",
                   paste(vis, collapse = ","), a$time, a$random, if (length(covs)) paste(covs, collapse = ",") else "none",
                   if (is.null(by)) "none" else by),
           sprintf("%d variables; time effect FDR < %.2f: %d%s; singular fits: %d; failed: %d",
                   nrow(res), fdr, sum(res$fdr_time < fdr, na.rm = TRUE),
                   if ("fdr_interaction" %in% names(res)) sprintf("; time x %s FDR < %.2f: %d", by, fdr, sum(res$fdr_interaction < fdr, na.rm = TRUE)) else "",
                   sum(res$singular %in% TRUE), sum(!is.na(res$error))),
           "Tests: ML likelihood-ratio tests of nested models; estimates: REML emmeans (vs first visit, unadjusted within variable).",
           "", utils::capture.output(print(utils::head(res[, intersect(c("variable", "label", "n_obs", "n_subjects", "p_time", "fdr_time",
                                                                        "p_interaction", "fdr_interaction", grep("^p_time_", names(res), value = TRUE), "icc", "singular"), names(res))], 15),
                                           row.names = FALSE, digits = 3)))
writeLines(lines, file.path(out_dir, "summary.txt"))
lomics_session_log(cfg, paste0("mixed_", a$mode, "_", tag))
cat(lines, sep = "\n")
