#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## classification.R — compare classifiers with repeated, nested cross-validation.
##
## Usage (from the project root):
##   Rscript .claude/skills/ml-classification/scripts/classification.R \
##     --outcome=binary|<column> --layers=all --visit=V1 [--change_from=V1] \
##     [--models=plsda,splsda,lda,enet,rf,svm] [--folds=5] [--repeats=10] \
##     [--prefilter=200] [--features_file=path/to/list.txt] [--seed=42]
##
## Models (all tuned only on the training part of each outer fold)
##   plsda   PLS-DA (pls::plsr on class indicator matrix; mixOmics::plsda if
##           installed); ncomp 1-5 tuned by inner CV. VIP scores reported.
##   splsda  sparse PLS-DA (mixOmics only; keepX tuned by inner CV)
##   lda     PCA-LDA: PCA on the training fold, LDA (MASS) on the top PCs
##           (number tuned); plain LDA is singular when p > n
##   enet    elastic-net logistic / multinomial regression (glmnet, alpha 0.5,
##           lambda by inner CV)
##   rf      random forest (ranger, 1000 trees, probability forest)
##   svm     linear SVM (e1071), cost tuned by inner CV
##
## Leakage control: scaling, the univariate prefilter (top --prefilter
## features by limma F on the training fold), PCA and tuning are all done
## inside each outer training fold. With --visit=all, folds are grouped by
## subject so no subject is in both training and test.
## --features_file restricts candidate features to a pre-specified list —
## only legitimate if that list was NOT derived from these same subjects.
##
## Outputs (results/classification/<outcome>_<tag>/):
## --permute_seed=N shuffles the outcome between subjects before CV; run with
## N = 1..20 (e.g. --repeats=2) to build a null distribution of AUC and
## report an empirical p-value for the real result.
##
##   cv_metrics.csv (per model x repeat), cv_summary.csv (mean, SD),
##   oof_predictions.csv, final_importance_<model>.csv, roc.pdf
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(SummarizedExperiment); library(limma) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", outcome = "binary", layers = "all", visit = NULL,
                      change_from = NULL, models = "plsda,splsda,lda,enet,rf,svm", folds = "5",
                      repeats = "10", prefilter = "200", features_file = NULL, seed = NULL,
                      permute_seed = NULL))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
seed <- as.integer(if (is.null(a$seed)) cfg$seed else a$seed)
models <- strsplit(a$models, ",")[[1]]
K <- as.integer(a$folds); R <- as.integer(a$repeats); pref <- as.integer(a$prefilter)
has_mix <- requireNamespace("mixOmics", quietly = TRUE)
if ("splsda" %in% models && !has_mix) { message("mixOmics not installed — dropping splsda"); models <- setdiff(models, "splsda") }

oname <- lomics_outcome(cfg, a$outcome)
d <- lomics_xy(cfg, layers, visit, a$change_from)
meta <- d$meta; X <- d$X
ok <- !is.na(meta[[oname]]); meta <- meta[ok, , drop = FALSE]; X <- X[ok, , drop = FALSE]
if (!is.null(a$features_file)) {
  keep <- intersect(readLines(a$features_file), colnames(X))
  if (!length(keep)) stop("No features from --features_file match (multi-layer names are <layer>__<feature>).")
  X <- X[, keep, drop = FALSE]
}
X <- X[, apply(X, 2, stats::var) > 0, drop = FALSE]
y <- factor(meta[[oname]]); lev <- levels(y); binary <- nlevels(y) == 2
if (!is.null(a$permute_seed)) {               # null run: shuffle labels between subjects
  set.seed(as.integer(a$permute_seed)); us <- unique(meta$subject_id)
  lab <- stats::setNames(as.character(y)[match(us, meta$subject_id)], us)
  lab <- stats::setNames(sample(lab), us); y <- factor(lab[meta$subject_id], levels = lev)
  message("PERMUTED LABELS (seed ", a$permute_seed, ") — this run estimates the null distribution")
}
if (min(table(y)) < K) stop("Smallest class has fewer subjects than folds.")
groups <- meta$subject_id
outdir <- dirname(lomics_out(cfg, "classification", paste0(oname, "_", lomics_tag(layers, visit, a$change_from,
                                if (!is.null(a$permute_seed)) paste0("perm", a$permute_seed))), "x"))
message(d$description, " | ", oname, ": ", paste(names(table(y)), table(y), sep = "=", collapse = ", "),
        " | models: ", paste(models, collapse = ", "))

## ---- helpers -------------------------------------------------------------
strat_folds <- function(y, groups, k) {             # stratified, grouped by subject
  g <- unique(data.frame(groups, y = as.character(y)))
  g <- g[!duplicated(g$groups), ]
  f <- integer(nrow(g))
  for (cl in unique(g$y)) { i <- which(g$y == cl); f[i] <- sample(rep(seq_len(k), length.out = length(i))) }
  f[match(groups, g$groups)]
}
prefilter <- function(Xtr, ytr, n) {
  if (ncol(Xtr) <= n) return(colnames(Xtr))
  fit <- eBayes(lmFit(t(Xtr), stats::model.matrix(~ ytr)))
  tt <- topTable(fit, coef = 2:nlevels(ytr), number = Inf, sort.by = "none")
  rownames(tt)[order(tt$P.Value)[seq_len(n)]]
}
inner_cv <- function(Xtr, ytr, fit_fun, grid, k = 5) {   # returns best grid value (accuracy)
  f <- strat_folds(ytr, seq_along(ytr), min(k, min(table(ytr))))
  acc <- sapply(grid, function(g) mean(unlist(lapply(unique(f), function(j) {
    p <- tryCatch(fit_fun(Xtr[f != j, , drop = FALSE], ytr[f != j], Xtr[f == j, , drop = FALSE], g), error = function(e) NULL)
    if (is.null(p)) return(NA)
    lev[max.col(p, ties.method = "first")] == as.character(ytr[f == j])
  })), na.rm = TRUE))
  grid[which.max(acc)]
}
softmax_clip <- function(P) { P <- pmax(P, 1e-6); P / rowSums(P) }

## each fitter: (Xtr, ytr, Xte, tuning) -> probability matrix (n_te x levels)
fitters <- list(
  plsda = function(Xtr, ytr, Xte, nc) {
    if (has_mix) { m <- mixOmics::plsda(Xtr, ytr, ncomp = nc, scale = TRUE)
      P <- stats::predict(m, Xte)$predict[, , nc] } else {
      Y <- stats::model.matrix(~ ytr - 1); colnames(Y) <- levels(ytr)
      m <- pls::plsr(Y ~ Xtr, ncomp = nc, scale = TRUE, method = "oscorespls")
      P <- stats::predict(m, newdata = Xte, ncomp = nc)[, , 1] }
    P <- matrix(P, nrow = nrow(Xte)); colnames(P) <- levels(ytr); softmax_clip(P) },
  splsda = function(Xtr, ytr, Xte, kx) {
    m <- mixOmics::splsda(Xtr, ytr, ncomp = 2, keepX = rep(kx, 2))
    P <- stats::predict(m, Xte)$predict[, , 2]; P <- matrix(P, nrow = nrow(Xte)); colnames(P) <- levels(ytr); softmax_clip(P) },
  lda = function(Xtr, ytr, Xte, npc) {
    pc <- stats::prcomp(Xtr, center = TRUE, scale. = TRUE); npc <- min(npc, ncol(pc$x))
    m <- MASS::lda(pc$x[, 1:npc, drop = FALSE], ytr)
    stats::predict(m, stats::predict(pc, Xte)[, 1:npc, drop = FALSE])$posterior },
  enet = function(Xtr, ytr, Xte, g) {
    m <- glmnet::cv.glmnet(Xtr, ytr, family = if (nlevels(ytr) == 2) "binomial" else "multinomial", alpha = 0.5,
                           nfolds = min(5, min(table(ytr))), type.measure = "deviance")
    p <- stats::predict(m, Xte, s = "lambda.min", type = "response")
    P <- if (nlevels(ytr) == 2) { p <- as.numeric(p); cbind(1 - p, p) } else p[, , 1]
    colnames(P) <- levels(ytr); P },
  rf = function(Xtr, ytr, Xte, g) {
    colnames(Xtr) <- colnames(Xte) <- make.names(colnames(Xtr), unique = TRUE)
    m <- ranger::ranger(x = data.frame(Xtr), y = ytr, probability = TRUE, num.trees = 1000, seed = seed)
    stats::predict(m, data.frame(Xte))$predictions },
  svm = function(Xtr, ytr, Xte, cost) {
    m <- e1071::svm(Xtr, ytr, kernel = "linear", cost = cost, probability = TRUE, scale = TRUE)
    attr(stats::predict(m, Xte, probability = TRUE), "probabilities")[, levels(ytr), drop = FALSE] })
grids <- list(plsda = 1:5, splsda = c(5, 10, 20, 50), lda = c(2, 3, 5, 8, 10), enet = NA, rf = NA, svm = c(0.001, 0.01, 0.1, 1))

metrics <- function(P, ytrue) {
  pred <- factor(lev[max.col(P, ties.method = "first")], lev)
  sens <- sapply(lev, function(l) mean(pred[ytrue == l] == l))
  auc <- tryCatch(if (binary) as.numeric(pROC::auc(pROC::roc(ytrue, P[, 2], levels = lev, direction = "<", quiet = TRUE)))
                  else as.numeric(pROC::multiclass.roc(ytrue, P, quiet = TRUE)$auc), error = function(e) NA)
  c(AUC = auc, accuracy = mean(pred == ytrue), balanced_accuracy = mean(sens),
    if (binary) c(sensitivity = sens[[2]], specificity = sens[[1]]),
    brier = mean(rowSums((P - stats::model.matrix(~ ytrue - 1))^2)))
}

## ---- repeated nested CV ----------------------------------------------------
set.seed(seed)
all_m <- list(); oof <- list()
for (r in seq_len(R)) {
  fold <- strat_folds(y, groups, K)
  Pm <- lapply(stats::setNames(models, models), function(m) matrix(NA, nrow(X), nlevels(y), dimnames = list(NULL, lev)))
  for (k in seq_len(K)) {
    tr <- fold != k; te <- !tr
    ytr <- droplevels(y[tr]); if (nlevels(ytr) < nlevels(y)) next
    mu <- colMeans(X[tr, ]); sdv <- apply(X[tr, ], 2, stats::sd); sdv[sdv == 0] <- 1
    Xtr <- sweep(sweep(X[tr, , drop = FALSE], 2, mu), 2, sdv, "/")
    Xte <- sweep(sweep(X[te, , drop = FALSE], 2, mu), 2, sdv, "/")
    fs <- prefilter(Xtr, ytr, pref); Xtr <- Xtr[, fs, drop = FALSE]; Xte <- Xte[, fs, drop = FALSE]
    for (m in models) {
      g <- grids[[m]]
      best <- if (length(g) > 1) inner_cv(Xtr, ytr, fitters[[m]], g) else g
      P <- tryCatch(fitters[[m]](Xtr, ytr, Xte, best), error = function(e) { message(m, " fold ", k, ": ", conditionMessage(e)); NULL })
      if (!is.null(P)) Pm[[m]][te, ] <- P[, lev, drop = FALSE]
    }
  }
  for (m in models) {
    okp <- stats::complete.cases(Pm[[m]])
    all_m[[length(all_m) + 1]] <- data.frame(model = m, repeat_id = r, t(metrics(Pm[[m]][okp, , drop = FALSE], y[okp])))
    oof[[length(oof) + 1]] <- data.frame(observation = rownames(X), repeat_id = r, model = m, truth = y, Pm[[m]], check.names = FALSE)
  }
  message(sprintf("repeat %d/%d done", r, R))
}
cvm <- do.call(rbind, all_m); oofd <- do.call(rbind, oof)
utils::write.csv(cvm, file.path(outdir, "cv_metrics.csv"), row.names = FALSE)
utils::write.csv(oofd, file.path(outdir, "oof_predictions.csv"), row.names = FALSE)
mcols <- setdiff(names(cvm), c("model", "repeat_id"))
summ <- do.call(rbind, lapply(split(cvm, cvm$model), function(z)
  data.frame(model = z$model[1], t(sapply(mcols, function(cn) sprintf("%.3f (%.3f)", mean(z[[cn]], na.rm = TRUE), stats::sd(z[[cn]], na.rm = TRUE)))), check.names = FALSE)))
summ <- summ[order(-sapply(split(cvm$AUC, cvm$model), mean, na.rm = TRUE)[summ$model]), ]
utils::write.csv(summ, file.path(outdir, "cv_summary.csv"), row.names = FALSE)
chance <- max(table(y)) / length(y)

## ---- final models on all data (for interpretation only) ------------------
Xa <- scale(X); fs <- prefilter(Xa, y, pref); Xa <- Xa[, fs, drop = FALSE]
imp_out <- function(df, m) { df$layer <- d$feature_layer[df$feature]; utils::write.csv(df, file.path(outdir, paste0("final_importance_", m, ".csv")), row.names = FALSE) }
if ("plsda" %in% models) {
  nc <- inner_cv(Xa, y, fitters$plsda, 1:5)
  Y <- stats::model.matrix(~ y - 1)
  pm <- pls::plsr(Y ~ Xa, ncomp = nc, scale = FALSE, method = "oscorespls")
  W <- pm$loading.weights[, 1:nc, drop = FALSE]; Tt <- pm$scores[, 1:nc, drop = FALSE]
  Q <- pm$Yloadings[, 1:nc, drop = FALSE]
  ss <- colSums(Q^2) * colSums(Tt^2)                    # Y-variance explained per component
  vip <- sqrt(nrow(W) * rowSums(sweep(W^2 / matrix(colSums(W^2), nrow(W), nc, byrow = TRUE), 2, ss, "*")) / sum(ss))
  imp_out(data.frame(feature = names(vip), VIP = vip)[order(-vip), ], "plsda")
  pdf(file.path(outdir, "plsda_scores.pdf"), width = 6, height = 5)
  plot(Tt[, 1], if (nc > 1) Tt[, 2] else rep(0, nrow(Tt)), col = as.integer(y), pch = 16,
       xlab = "PLS comp 1", ylab = if (nc > 1) "PLS comp 2" else "", main = sprintf("PLS-DA scores (ncomp = %d, full data — apparent fit)", nc))
  legend("topright", lev, col = seq_along(lev), pch = 16, bty = "n"); invisible(dev.off())
}
if ("enet" %in% models) {
  m <- glmnet::cv.glmnet(Xa, y, family = if (binary) "binomial" else "multinomial", alpha = 0.5)
  co <- stats::coef(m, s = "lambda.min")
  co <- if (is.list(co)) do.call(cbind, lapply(co, as.matrix)) else as.matrix(co)
  co <- co[-1, , drop = FALSE]; keep <- rowSums(co != 0) > 0
  imp_out(data.frame(feature = rownames(co)[keep], co[keep, , drop = FALSE], check.names = FALSE), "enet")
}
if ("rf" %in% models) {
  dfa <- data.frame(Xa, check.names = FALSE); nm <- colnames(dfa); colnames(dfa) <- make.names(nm, unique = TRUE)
  rfm <- ranger::ranger(x = dfa, y = y, num.trees = 1000, importance = "permutation", seed = seed)
  im <- ranger::importance(rfm); imp_out(data.frame(feature = nm, importance = im)[order(-im), ], "rf")
}

## ---- ROC plot (pooled out-of-fold predictions, first repeat) -------------
if (binary) {
  pdf(file.path(outdir, "roc.pdf"), width = 6, height = 6)
  first <- oofd[oofd$repeat_id == 1, ]
  for (i in seq_along(models)) {
    z <- first[first$model == models[i] & !is.na(first[[lev[2]]]), ]
    rc <- pROC::roc(z$truth, z[[lev[2]]], levels = lev, direction = "<", quiet = TRUE)
    if (i == 1) plot(rc, col = i, main = "Out-of-fold ROC (repeat 1)") else plot(rc, col = i, add = TRUE)
  }
  legend("bottomright", sprintf("%s (AUC %.2f)", models, sapply(models, function(m) mean(cvm$AUC[cvm$model == m], na.rm = TRUE))),
         col = seq_along(models), lwd = 2, bty = "n")
  invisible(dev.off())
}
cat(sprintf("\n%d-fold CV x %d repeats; prefilter top %d features inside folds; chance accuracy = %.3f\n", K, R, pref, chance))
cat("Mean (SD) over repeats:\n"); print(summ, row.names = FALSE)
message("Outputs in ", outdir)
lomics_session_log(cfg, "classification")
