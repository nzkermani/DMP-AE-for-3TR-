#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## regression.R — linear / logistic regression for omics and clinical data.
##
## Usage (from the project root):
##
## (1) --mode=features   omics-wide association with an outcome (one model per
##     feature, limma moderated statistics):  feature ~ outcome + covariates
##       --outcome=continuous|binary|<column>  --layers=all  --visit=V1
##       [--change_from=V1] [--covariates=age,sex] [--no_covariates]
##     For visit=all the subject is a blocking factor (duplicateCorrelation).
##
## (2) --mode=outcome    per-feature regression with the clinical outcome as
##     response:  outcome ~ feature + covariates
##       linear (gaussian) for continuous outcomes, logistic for binary
##       (odds ratio per 1 SD of feature)
##
## (3) --mode=formula    one multivariable model you specify
##       --formula="fev1 ~ age + sex + proteomics__PROT21"
##       [--family=gaussian|binomial|poisson]  --layers=...  --visit=V1
##     Feature columns are named <layer>__<feature>. Reports coefficients with
##     95% CI, VIF (car), and diagnostic plots.
##
## Outputs: results/regression/<mode>_<outcome>_<tag>/...
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(SummarizedExperiment); library(limma) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", mode = "features", outcome = "continuous",
                      layers = "all", visit = NULL, change_from = NULL, covariates = NULL,
                      no_covariates = "FALSE", formula = NULL, family = "gaussian", fdr = "0.05"))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
fdr <- as.numeric(a$fdr)

## ---------------------------------------------------------------- formula
if (a$mode == "formula") {
  if (is.null(a$formula)) stop("--formula is required for mode=formula")
  fo <- stats::as.formula(a$formula)
  d <- lomics_xy(cfg, layers, visit, a$change_from)
  df <- data.frame(d$meta, d$X, check.names = FALSE)
  missing_vars <- setdiff(all.vars(fo), names(df))
  if (length(missing_vars)) stop("Variables not found: ", paste(missing_vars, collapse = ", "),
                                 ". Feature columns are named <layer>__<feature> when >1 layer is loaded.")
  fam <- switch(a$family, gaussian = stats::gaussian(), binomial = stats::binomial(), poisson = stats::poisson())
  y <- df[[all.vars(fo)[1]]]
  if (a$family == "binomial" && !is.numeric(y)) df[[all.vars(fo)[1]]] <- factor(y)
  fit <- stats::glm(fo, data = df, family = fam)
  co <- summary(fit)$coefficients
  ci <- tryCatch(suppressMessages(stats::confint(fit)), error = function(e) stats::confint.default(fit))
  tab <- data.frame(term = rownames(co), estimate = co[, 1], se = co[, 2], ci_low = ci[, 1], ci_high = ci[, 2],
                    statistic = co[, 3], p = co[, 4], row.names = NULL)
  if (a$family == "binomial") { tab$odds_ratio <- exp(tab$estimate); tab$or_low <- exp(tab$ci_low); tab$or_high <- exp(tab$ci_high) }
  outdir <- dirname(lomics_out(cfg, "regression", paste0("formula_", lomics_tag(layers, visit, a$change_from)), "x"))
  utils::write.csv(tab, file.path(outdir, "coefficients.csv"), row.names = FALSE)
  out <- c(sprintf("Model: %s  (family %s)", a$formula, a$family), sprintf("n = %d", stats::nobs(fit)),
           utils::capture.output(print(tab, row.names = FALSE)))
  if (a$family == "gaussian") {
    s <- summary(stats::lm(fo, data = df))
    out <- c(out, sprintf("R2 = %.3f, adjusted R2 = %.3f", s$r.squared, s$adj.r.squared))
  } else out <- c(out, sprintf("AIC = %.1f; deviance %.1f on %d df", stats::AIC(fit), fit$deviance, fit$df.residual))
  if (length(attr(stats::terms(fo), "term.labels")) > 1 && requireNamespace("car", quietly = TRUE)) {
    v <- tryCatch(car::vif(fit), error = function(e) NULL)
    if (!is.null(v)) out <- c(out, "VIF (values > 5 suggest collinearity):", utils::capture.output(print(round(v, 2))))
  }
  pdf(file.path(outdir, "diagnostics.pdf")); graphics::par(mfrow = c(2, 2)); plot(fit); invisible(dev.off())
  writeLines(out, file.path(outdir, "model_summary.txt")); cat(out, sep = "\n")
  lomics_session_log(cfg, "regression_formula"); quit(save = "no")
}

## ---------------------------------------------------------------- features / outcome
outcome <- lomics_outcome(cfg, a$outcome)
outdir <- dirname(lomics_out(cfg, "regression", paste0(a$mode, "_", outcome, "_", lomics_tag(layers, visit, a$change_from)), "x"))
smry <- list()
for (layer in layers) {
  d <- lomics_xy(cfg, layer, visit, a$change_from)
  meta <- d$meta; X <- d$X
  covs <- if (isTRUE(as.logical(a$no_covariates))) character() else
    if (!is.null(a$covariates)) strsplit(a$covariates, ",")[[1]] else lomics_covariates(cfg, meta)
  covs <- setdiff(covs, outcome)
  ok <- stats::complete.cases(meta[, c(outcome, covs), drop = FALSE])
  meta <- meta[ok, , drop = FALSE]; X <- X[ok, , drop = FALSE]
  y <- meta[[outcome]]
  is_binary <- length(unique(y)) == 2
  if (a$mode == "features") {
    meta$.y <- if (is.numeric(y) && !is_binary) y else factor(y)
    design <- stats::model.matrix(stats::as.formula(paste("~ .y", if (length(covs)) paste("+", paste(covs, collapse = "+")))), meta)
    coef_name <- grep("^\\.y", colnames(design), value = TRUE)
    if (identical(visit, "all")) {
      dc <- duplicateCorrelation(t(X), design, block = meta$subject_id)
      fit <- lmFit(t(X), design, block = meta$subject_id, correlation = dc$consensus.correlation)
    } else fit <- lmFit(t(X), design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = if (length(coef_name) > 1) "F" else "P")
    tab <- data.frame(feature = rownames(tt), tt, row.names = NULL, check.names = FALSE)
    note <- if (length(coef_name) == 1) sprintf("logFC = change in feature per unit of %s%s", outcome,
                                                 if (is_binary) sprintf(" (%s vs %s)", levels(meta$.y)[2], levels(meta$.y)[1]) else "")
            else "F-test across outcome levels"
  } else if (a$mode == "outcome") {
    fam <- if (is_binary) stats::binomial() else stats::gaussian()
    yy <- if (is_binary) as.integer(factor(y)) - 1 else y
    Xz <- scale(X)
    tab <- do.call(rbind, lapply(colnames(Xz), function(f) {
      df <- data.frame(y = yy, feat = Xz[, f], meta[, covs, drop = FALSE])
      fit <- tryCatch(stats::glm(stats::as.formula(paste("y ~ feat", if (length(covs)) paste("+", paste(covs, collapse = "+")))),
                                 data = df, family = fam), error = function(e) NULL)
      if (is.null(fit)) return(data.frame(feature = f, estimate = NA, se = NA, p = NA))
      co <- summary(fit)$coefficients["feat", ]
      data.frame(feature = f, estimate = co[1], se = co[2], p = co[4])
    }))
    if (is_binary) { tab$odds_ratio_per_SD <- exp(tab$estimate); tab$or_low <- exp(tab$estimate - 1.96 * tab$se); tab$or_high <- exp(tab$estimate + 1.96 * tab$se) }
    tab$adj.P.Val <- stats::p.adjust(tab$p, "BH"); tab <- tab[order(tab$p), ]
    note <- if (is_binary) sprintf("logistic: P(%s = %s) per 1 SD of feature", outcome, levels(factor(y))[2])
            else sprintf("linear: change in %s per 1 SD of feature", outcome)
    if (identical(visit, "all")) warning("mode=outcome with visit=all treats repeated samples as independent; use one visit.")
  } else stop("Unknown --mode")
  tab$layer <- layer
  utils::write.csv(tab, file.path(outdir, paste0(layer, "_regression.csv")), row.names = FALSE)
  smry[[layer]] <- data.frame(layer = layer, n = nrow(X), covariates = paste(covs, collapse = "+"),
                              n_fdr = sum(tab$adj.P.Val < fdr, na.rm = TRUE), interpretation = note)
}
smry <- do.call(rbind, smry); names(smry)[4] <- paste0("n_fdr<", fdr)
utils::write.csv(smry, file.path(outdir, "summary.csv"), row.names = FALSE)
print(smry, row.names = FALSE); message("Outputs in ", outdir)
lomics_session_log(cfg, paste0("regression_", a$mode))
