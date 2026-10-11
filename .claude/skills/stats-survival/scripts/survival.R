#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## survival.R — time-to-event analysis (Kaplan-Meier, Cox, penalised Cox).
##
## Usage (from the project root):
##   Rscript .claude/skills/stats-survival/scripts/survival.R --mode=<mode> [options]
##
##   --mode=km        Kaplan-Meier + log-rank by a grouping
##       --strata=group|<metadata column>   or
##       --strata_file=results/clustering/.../assignments.csv --strata_col=cluster
##   --mode=cox_features   univariable-per-feature Cox models adjusted for
##       covariates: Surv(time, event) ~ feature(z) + covariates; HR per 1 SD, BH FDR
##   --mode=cox_penalised  elastic-net Cox (glmnet) on all features of the
##       chosen layers; repeated CV C-index, selected features, risk score;
##       KM of risk-score tertiles from cross-validated (out-of-fold) scores
##       [--alpha=0.5] [--folds=5] [--repeats=10]
##   --mode=cox_formula    one multivariable Cox model
##       --formula="age + sex + proteomics__PROT21"
##
##   common: --time=<col> --event=<col> (default from config outcomes)
##           --layers=all --visit=V1 [--change_from=V1] [--covariates=age,sex]
##
## Baseline (V1) measurements are the default: survival outcomes are
## subject-level, so predictors measured after follow-up starts create
## immortal-time bias. Use later visits only with a landmark analysis
## (time measured from that visit) — the script does NOT do this for you.
##
## Outputs: results/survival/<mode>_<tag>/
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(survival); library(SummarizedExperiment) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", mode = "cox_features", layers = "all",
                      visit = NULL, change_from = NULL, time = NULL, event = NULL,
                      strata = "group", strata_file = NULL, strata_col = "cluster",
                      covariates = NULL, no_covariates = "FALSE", formula = NULL,
                      alpha = "0.5", folds = "5", repeats = "10", fdr = "0.05", seed = NULL))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
if (identical(visit, "all")) stop("Survival analysis needs one row per subject: choose a single --visit.")
tcol <- if (is.null(a$time)) cfg$outcomes$survival_time else a$time
ecol <- if (is.null(a$event)) cfg$outcomes$survival_event else a$event
seed <- as.integer(if (is.null(a$seed)) cfg$seed else a$seed)
if (visit != cfg$visits[1]) message("NOTE: predictors from ", visit, " — check this is a landmark analysis (see header).")

d <- lomics_xy(cfg, layers, visit, a$change_from)
meta <- d$meta; X <- d$X
ok <- !is.na(meta[[tcol]]) & !is.na(meta[[ecol]]) & meta[[tcol]] > 0
meta <- meta[ok, , drop = FALSE]; X <- X[ok, , drop = FALSE]
S <- Surv(meta[[tcol]], as.integer(meta[[ecol]]))
covs <- if (isTRUE(as.logical(a$no_covariates))) character() else
  if (!is.null(a$covariates)) strsplit(a$covariates, ",")[[1]] else lomics_covariates(cfg, meta)
tag <- lomics_tag(layers, visit, a$change_from)
outdir <- dirname(lomics_out(cfg, "survival", paste0(a$mode, "_", tag), "x"))
message(sprintf("%s | n = %d, events = %d", d$description, nrow(meta), sum(S[, 2])))

km_plot <- function(fit, title, file) {
  pdf(file, width = 7, height = 5)
  if (requireNamespace("survminer", quietly = TRUE)) {
    print(survminer::ggsurvplot(fit, data = attr(fit, "data"), risk.table = TRUE, pval = TRUE, title = title))
  } else {
    plot(fit, col = seq_along(fit$strata), xlab = "Time", ylab = "Survival probability", main = title)
    legend("bottomleft", names(fit$strata), col = seq_along(fit$strata), lty = 1, bty = "n")
  }
  invisible(dev.off())
}

if (a$mode == "km") {
  if (!is.null(a$strata_file)) {
    sf <- utils::read.csv(a$strata_file)
    key <- intersect(c("observation", "subject_id", "sample_id"), names(sf))[1]
    meta$strata <- factor(sf[[a$strata_col]][match(rownames(meta), sf[[key]])])
    sname <- a$strata_col
  } else { meta$strata <- factor(meta[[a$strata]]); sname <- a$strata }
  df <- data.frame(time = meta[[tcol]], event = as.integer(meta[[ecol]]), strata = meta$strata)
  df <- df[!is.na(df$strata), ]
  fit <- survfit(Surv(time, event) ~ strata, data = df); attr(fit, "data") <- df
  lr <- survdiff(Surv(time, event) ~ strata, data = df)
  p <- stats::pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  cx <- coxph(Surv(time, event) ~ strata, data = df)
  out <- c(sprintf("Kaplan-Meier by %s (n = %d, events = %d)", sname, nrow(df), sum(df$event)),
           utils::capture.output(print(fit)), sprintf("Log-rank p = %.3g", p),
           "Cox HR vs reference level:", utils::capture.output(print(summary(cx)$conf.int)))
  km_plot(fit, sprintf("KM by %s (log-rank p = %.3g)", sname, p), file.path(outdir, "km.pdf"))
} else if (a$mode == "cox_features") {
  Xz <- scale(X)
  res <- do.call(rbind, lapply(colnames(Xz), function(f) {
    df <- data.frame(feat = Xz[, f], meta[, covs, drop = FALSE])
    fit <- tryCatch(coxph(stats::as.formula(paste("S ~ feat", if (length(covs)) paste("+", paste(covs, collapse = "+")))), data = df),
                    error = function(e) NULL, warning = function(w) NULL)
    if (is.null(fit)) return(data.frame(feature = f, HR_per_SD = NA, HR_low = NA, HR_high = NA, p = NA, ph_test_p = NA))
    ci <- summary(fit)$conf.int["feat", ]
    ph <- tryCatch(cox.zph(fit)$table["feat", "p"], error = function(e) NA)
    data.frame(feature = f, HR_per_SD = ci[1], HR_low = ci[3], HR_high = ci[4],
               p = summary(fit)$coefficients["feat", "Pr(>|z|)"], ph_test_p = ph)
  }))
  res$fdr <- stats::p.adjust(res$p, "BH"); res$layer <- d$feature_layer[res$feature]
  res <- res[order(res$p), ]
  utils::write.csv(res, file.path(outdir, "cox_per_feature.csv"), row.names = FALSE)
  out <- c(sprintf("Per-feature Cox: Surv(%s, %s) ~ feature(z)%s", tcol, ecol, if (length(covs)) paste(" +", paste(covs, collapse = " + ")) else ""),
           sprintf("n = %d, events = %d, features = %d, FDR < %s: %d", nrow(meta), sum(S[, 2]), nrow(res), a$fdr, sum(res$fdr < as.numeric(a$fdr), na.rm = TRUE)),
           sprintf("PH assumption violated (cox.zph p < 0.05) for %d features — check before reporting.", sum(res$ph_test_p < 0.05, na.rm = TRUE)),
           "Top 15:", utils::capture.output(print(utils::head(res, 15), row.names = FALSE)))
  if (sum(S[, 2]) < 10 * (1 + length(covs))) out <- c(out, "WARNING: fewer than ~10 events per model parameter; estimates are unstable.")
} else if (a$mode == "cox_penalised") {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install glmnet.")
  alpha <- as.numeric(a$alpha); K <- as.integer(a$folds); R <- as.integer(a$repeats)
  Xz <- scale(X); y <- cbind(time = meta[[tcol]], status = as.integer(meta[[ecol]]))
  set.seed(seed)
  cidx <- c(); oof_all <- matrix(NA, nrow(Xz), R)
  for (r in seq_len(R)) {
    fold <- sample(rep(seq_len(K), length.out = nrow(Xz)))
    for (k in seq_len(K)) {
      tr <- fold != k
      cv <- glmnet::cv.glmnet(Xz[tr, ], y[tr, ], family = "cox", alpha = alpha, nfolds = min(5, sum(tr) %/% 3))
      oof_all[!tr, r] <- as.numeric(stats::predict(cv, Xz[!tr, , drop = FALSE], s = "lambda.min", type = "link"))
    }
    cidx[r] <- concordance(Surv(y[, 1], y[, 2]) ~ oof_all[, r], reverse = TRUE)$concordance
  }
  final <- glmnet::cv.glmnet(Xz, y, family = "cox", alpha = alpha)
  b <- as.matrix(stats::coef(final, s = "lambda.min"))[, 1]; b <- b[b != 0]
  sel <- data.frame(feature = names(b), coef = b, HR_per_SD = exp(b), layer = d$feature_layer[names(b)])
  sel <- sel[order(-abs(sel$coef)), ]
  utils::write.csv(sel, file.path(outdir, "selected_features.csv"), row.names = FALSE)
  oof <- rowMeans(oof_all)
  grp <- cut(oof, stats::quantile(oof, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("low", "mid", "high"))
  df <- data.frame(time = y[, 1], event = y[, 2], strata = grp)
  utils::write.csv(data.frame(observation = rownames(Xz), oof_risk = oof, risk_group = grp),
                   file.path(outdir, "cv_risk_scores.csv"), row.names = FALSE)
  fit <- survfit(Surv(time, event) ~ strata, data = df); attr(fit, "data") <- df
  km_plot(fit, "Out-of-fold risk tertiles", file.path(outdir, "km_risk_tertiles.pdf"))
  out <- c(sprintf("Elastic-net Cox (alpha = %s), %d-fold CV x %d repeats, features z-scored", alpha, K, R),
           sprintf("n = %d, events = %d, candidate features = %d", nrow(Xz), sum(y[, 2]), ncol(Xz)),
           sprintf("Out-of-fold C-index: mean %.3f (SD %.3f; range %.3f-%.3f)", mean(cidx), stats::sd(cidx), min(cidx), max(cidx)),
           sprintf("Features selected in the final full-data model: %d (selection is unstable at small n; see ml-feature-selection for stability selection)", nrow(sel)),
           utils::capture.output(print(utils::head(sel, 20), row.names = FALSE)))
} else if (a$mode == "cox_formula") {
  if (is.null(a$formula)) stop("--formula is required")
  df <- data.frame(meta, X, check.names = FALSE)
  fit <- coxph(stats::as.formula(paste("Surv(", tcol, ",", ecol, ") ~", a$formula)), data = df)
  zph <- cox.zph(fit)
  out <- c(utils::capture.output(print(summary(fit))), "Proportional hazards test (cox.zph):", utils::capture.output(print(zph)))
  pdf(file.path(outdir, "ph_diagnostics.pdf")); plot(zph); invisible(dev.off())
} else stop("Unknown --mode")

writeLines(out, file.path(outdir, "summary.txt")); cat(out, sep = "\n")
message("Outputs in ", outdir)
lomics_session_log(cfg, paste0("survival_", a$mode))
