#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## feature_selection.R — rank and select features for an outcome with several
## complementary methods, then combine into a consensus ranking.
##
## Usage (from the project root):
##   Rscript .claude/skills/ml-feature-selection/scripts/feature_selection.R \
##     --outcome=binary|continuous|survival|<column> --layers=all --visit=V1 \
##     [--change_from=V1] [--methods=variance,univariate,stability_enet,rf,boruta] \
##     [--alpha=0.5] [--n_boot=100] [--pi_thr=0.6] [--top_n=50] [--seed=42]
##
## Methods
##   variance        unsupervised MAD ranking (does not use the outcome)
##   univariate      limma moderated t / F (classes or continuous) or
##                   per-feature Cox (survival); BH FDR
##   stability_enet  stability selection: elastic net (glmnet) on n_boot
##                   half-subsamples, lambda chosen by CV on each; selection
##                   frequency per feature; selected if freq >= pi_thr
##   rf              random forest (ranger) permutation importance; for
##                   survival a random survival forest
##   boruta          all-relevant selection (Boruta; skipped if not installed)
##
## IMPORTANT: features selected here on all subjects must NOT then be used to
## estimate classification accuracy on the same subjects (selection bias).
## ml-classification repeats selection inside each CV fold for that reason.
##
## Outputs (results/feature_selection/<outcome>_<tag>/):
##   ranking_<method>.csv, consensus_ranking.csv, selected_features.txt, plots.pdf
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(SummarizedExperiment); library(limma) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", outcome = "binary", layers = "all", visit = NULL,
                      change_from = NULL, methods = "variance,univariate,stability_enet,rf,boruta",
                      alpha = "0.5", n_boot = "100", pi_thr = "0.6", top_n = "50", seed = NULL))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
seed <- as.integer(if (is.null(a$seed)) cfg$seed else a$seed)
methods <- strsplit(a$methods, ",")[[1]]; top_n <- as.integer(a$top_n)
if (identical(visit, "all")) stop("Use a single --visit (one row per subject); repeated samples would leak between subsamples.")

d <- lomics_xy(cfg, layers, visit, a$change_from)
meta <- d$meta; X <- d$X[, apply(d$X, 2, stats::var) > 0, drop = FALSE]
is_surv <- a$outcome == "survival"
if (is_surv) {
  tcol <- cfg$outcomes$survival_time; ecol <- cfg$outcomes$survival_event
  ok <- !is.na(meta[[tcol]]) & !is.na(meta[[ecol]]) & meta[[tcol]] > 0
  y <- survival::Surv(meta[[tcol]][ok], meta[[ecol]][ok]); type <- "survival"; oname <- "survival"
} else {
  oname <- lomics_outcome(cfg, a$outcome); yy <- meta[[oname]]; ok <- !is.na(yy)
  type <- if (is.numeric(yy) && length(unique(yy)) > 2) "continuous" else "class"
  y <- if (type == "class") factor(yy[ok]) else yy[ok]
}
X <- X[ok, , drop = FALSE]; Xz <- scale(X)
outdir <- dirname(lomics_out(cfg, "feature_selection", paste0(oname, "_", lomics_tag(layers, visit, a$change_from)), "x"))
message(d$description, " | outcome: ", oname, " (", type, ")", if (type == "class") paste0(" ", paste(names(table(y)), table(y), sep = "=", collapse = ", ")))
rk <- list()
save_rank <- function(df, m) { df$rank <- seq_len(nrow(df)); df$layer <- d$feature_layer[df$feature]
  utils::write.csv(df, file.path(outdir, paste0("ranking_", m, ".csv")), row.names = FALSE); rk[[m]] <<- df }

if ("variance" %in% methods) {
  v <- apply(X, 2, stats::mad); o <- order(v, decreasing = TRUE)
  save_rank(data.frame(feature = colnames(X)[o], score = v[o]), "variance")
}
if ("univariate" %in% methods) {
  if (is_surv) {
    p <- apply(Xz, 2, function(x) tryCatch(summary(survival::coxph(y ~ x))$coefficients[1, 5], error = function(e) NA))
    df <- data.frame(feature = colnames(X), score = -log10(p), p = p)
  } else {
    des <- if (type == "class") stats::model.matrix(~ y) else stats::model.matrix(~ y)
    fit <- eBayes(lmFit(t(X), des), trend = TRUE)
    tt <- topTable(fit, coef = 2:ncol(des), number = Inf, sort.by = "none")
    df <- data.frame(feature = rownames(tt), score = -log10(tt$P.Value), p = tt$P.Value)
  }
  df$fdr <- stats::p.adjust(df$p, "BH"); df <- df[order(df$p), ]
  save_rank(df, "univariate")
}
if ("stability_enet" %in% methods) {
  fam <- if (is_surv) "cox" else if (type == "continuous") "gaussian" else if (nlevels(y) == 2) "binomial" else "multinomial"
  set.seed(seed); B <- as.integer(a$n_boot); sel <- matrix(0, B, ncol(Xz), dimnames = list(NULL, colnames(Xz)))
  for (b in seq_len(B)) {
    idx <- if (type == "class") unlist(lapply(split(seq_along(y), y), function(i) sample(i, ceiling(length(i) / 2))))
           else sample(nrow(Xz), floor(nrow(Xz) / 2))
    yb <- if (is_surv) y[idx] else y[idx]
    cv <- tryCatch(glmnet::cv.glmnet(Xz[idx, ], yb, family = fam, alpha = as.numeric(a$alpha), nfolds = 5),
                   error = function(e) NULL)
    if (is.null(cv)) next
    co <- stats::coef(cv, s = "lambda.min")
    nz <- if (is.list(co)) unique(unlist(lapply(co, function(m) rownames(m)[as.vector(m) != 0])))
          else rownames(co)[as.vector(co) != 0]
    sel[b, intersect(nz, colnames(Xz))] <- 1
  }
  fr <- colMeans(sel); o <- order(fr, decreasing = TRUE)
  save_rank(data.frame(feature = names(fr)[o], score = fr[o], selected = fr[o] >= as.numeric(a$pi_thr)), "stability_enet")
}
if ("rf" %in% methods) {
  if (!requireNamespace("ranger", quietly = TRUE)) message("ranger not installed — skipping rf") else {
    df <- data.frame(Xz, check.names = FALSE); colnames(df) <- make.names(colnames(Xz), unique = TRUE)
    map <- stats::setNames(colnames(Xz), colnames(df))
    rf <- if (is_surv) ranger::ranger(x = df, y = y, num.trees = 1000, importance = "permutation", seed = seed)
          else ranger::ranger(x = df, y = y, num.trees = 1000, importance = "permutation", seed = seed,
                              classification = type == "class")
    imp <- ranger::importance(rf); o <- order(imp, decreasing = TRUE)
    save_rank(data.frame(feature = map[names(imp)[o]], score = imp[o]), "rf")
    message(sprintf("RF out-of-bag %s: %.3f", if (type == "class") "error" else if (is_surv) "1 - C-index" else "MSE", rf$prediction.error))
  }
}
if ("boruta" %in% methods) {
  if (is_surv || !requireNamespace("Boruta", quietly = TRUE)) message("Boruta skipped (not installed or survival outcome)") else {
    set.seed(seed); bo <- Boruta::Boruta(x = data.frame(Xz, check.names = FALSE), y = y, maxRuns = 200)
    st <- Boruta::attStats(bo); st <- st[order(-st$meanImp), ]
    save_rank(data.frame(feature = rownames(st), score = st$meanImp, decision = st$decision), "boruta")
  }
}

## Consensus: mean rank across outcome-based methods (variance excluded)
sup <- setdiff(names(rk), "variance")
if (length(sup)) {
  ranks <- sapply(sup, function(m) rk[[m]]$rank[match(colnames(X), rk[[m]]$feature)])
  rownames(ranks) <- colnames(X)
  cons <- data.frame(feature = colnames(X), layer = d$feature_layer[colnames(X)], mean_rank = rowMeans(ranks, na.rm = TRUE),
                     ranks, check.names = FALSE)
  cons$n_methods_top_n <- rowSums(ranks <= top_n, na.rm = TRUE)
  cons <- cons[order(cons$mean_rank), ]
  utils::write.csv(cons, file.path(outdir, "consensus_ranking.csv"), row.names = FALSE)
  chosen <- if ("stability_enet" %in% names(rk) && any(rk$stability_enet$selected))
    rk$stability_enet$feature[rk$stability_enet$selected] else utils::head(cons$feature, top_n)
  writeLines(chosen, file.path(outdir, "selected_features.txt"))
  cat("Top 20 consensus features:\n"); print(utils::head(cons, 20), row.names = FALSE)
  cat(sprintf("\nselected_features.txt: %d features (%s)\n", length(chosen),
              if ("stability_enet" %in% names(rk) && any(rk$stability_enet$selected)) sprintf("stability selection freq >= %s", a$pi_thr) else sprintf("top %d by consensus rank", top_n)))
  if (length(sup) > 1) {
    cat("\nSpearman correlation between method rankings:\n")
    print(round(stats::cor(ranks, method = "spearman", use = "pairwise.complete.obs"), 2))
  }
  pdf(file.path(outdir, "plots.pdf"), width = 7, height = 6)
  if (length(sup) > 1) graphics::pairs(ranks, pch = 16, cex = .5, main = "Rank agreement between methods")
  if ("stability_enet" %in% names(rk)) { s <- rk$stability_enet
    plot(s$score, type = "h", xlab = "Feature (ordered)", ylab = "Selection frequency", main = "Stability selection")
    graphics::abline(h = as.numeric(a$pi_thr), lty = 2, col = 2) }
  invisible(dev.off())
}
message("Outputs in ", outdir)
lomics_session_log(cfg, "feature_selection")
