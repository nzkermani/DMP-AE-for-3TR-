## =============================================================================
## lomics_spls.R — multilevel sparse PLS-DA helpers (shared by the
## ml-multilevel-splsda skill and the use-case notebooks).
##
## Engine
##   * "mixOmics" — mixOmics::withinVariation() + mixOmics::splsda() + predict().
##     Used automatically when mixOmics is installed.
##   * "internal" — a dependency-free re-implementation of the same algorithm
##     (PLS2 / NIPALS, regression-mode deflation, L1 soft-thresholding that
##     keeps `keepX` variables per component, as in Le Cao et al. 2011,
##     BMC Bioinformatics 12:253). With keepX = all variables it reproduces
##     pls::plsr(method = "oscorespls") on the scaled class-indicator matrix
##     (checked in tests/test_spls_internal.R). Used when mixOmics is missing,
##     or on request. Results can differ slightly from mixOmics (convergence
##     tolerance, tie handling); report which engine was used.
##
## Multilevel decomposition (Westerhuis et al. 2010; Liquet et al. 2012):
## each sample minus its subject's mean, so only within-subject variation is
## modelled. For two visits per subject this is the paired design.
## =============================================================================

lomics_within_variation <- function(X, subject) {
  subject <- as.character(subject)
  if (requireNamespace("mixOmics", quietly = TRUE)) {
    Xw <- mixOmics::withinVariation(X = X, design = data.frame(subject))
    dimnames(Xw) <- dimnames(X)
    return(Xw)
  }
  means <- rowsum(X, subject) / as.vector(table(subject)[unique(sort(subject))])
  X - means[subject, , drop = FALSE]
}

.lomics_soft_keep <- function(a, keep) {
  p <- length(a)
  if (keep >= p) return(a)
  thr <- sort(abs(a), decreasing = TRUE)[keep + 1]
  sign(a) * pmax(abs(a) - thr, 0)
}

## Internal sPLS-DA fit. X: n x p numeric; y: factor (2+ levels).
.lomics_spls_internal <- function(X, y, ncomp = 2, keepX = rep(ncol(X), ncomp),
                                  scale = TRUE, tol = 1e-06, max_iter = 500) {
  y <- droplevels(as.factor(y))
  Y <- stats::model.matrix(~ y - 1); colnames(Y) <- levels(y)
  n <- nrow(X)
  xm <- colMeans(X); Xc <- X - rep(xm, each = n)
  xs <- if (scale) sqrt(colSums(Xc^2) / (n - 1)) else rep(1, ncol(X))
  xs[!is.finite(xs) | xs == 0] <- 1
  ym <- colMeans(Y); Yc <- Y - rep(ym, each = n)
  ys <- sqrt(colSums(Yc^2) / (n - 1)); ys[ys == 0] <- 1
  Xh <- Xc / rep(xs, each = n)
  Yh <- Yc / rep(ys, each = n)
  p <- ncol(X); q <- ncol(Y)
  A <- matrix(0, p, ncomp, dimnames = list(colnames(X), paste0("comp", 1:ncomp)))
  Pm <- A; Cm <- matrix(0, q, ncomp); Tm <- matrix(0, nrow(X), ncomp)
  keepX <- rep_len(pmin(keepX, p), ncomp)
  for (h in seq_len(ncomp)) {
    M <- crossprod(Xh, Yh)
    sv <- svd(M, nu = 1, nv = 1)
    a <- sv$u[, 1]; b <- sv$v[, 1]
    for (it in seq_len(max_iter)) {
      a_new <- .lomics_soft_keep(drop(M %*% b), keepX[h])
      if (all(a_new == 0)) break
      a_new <- a_new / sqrt(sum(a_new^2))
      b_new <- drop(crossprod(M, a_new)); b_new <- b_new / sqrt(sum(b_new^2))
      conv <- sum((a_new - a)^2) < tol
      a <- a_new; b <- b_new
      if (conv) break
    }
    if (sum(a) == 0 && all(a == 0)) stop("sPLS: empty loading vector on component ", h)
    ## sign convention: largest absolute loading positive
    if (a[which.max(abs(a))] < 0) a <- -a
    t <- drop(Xh %*% a); tt <- sum(t^2)
    pl <- drop(crossprod(Xh, t)) / tt
    cl <- drop(crossprod(Yh, t)) / tt
    Xh <- Xh - outer(t, pl)
    Yh <- Yh - outer(t, cl)                       # regression-mode deflation
    A[, h] <- a; Pm[, h] <- pl; Cm[, h] <- cl; Tm[, h] <- t
  }
  Wstar <- A %*% solve(crossprod(Pm, A))
  structure(list(engine = "internal", loadings = A, variates = Tm, Wstar = Wstar, C = Cm,
                 xm = xm, xs = xs, ym = ym, ys = ys, levels = levels(y), ncomp = ncomp,
                 keepX = keepX, features = colnames(X)),
            class = "lomics_spls")
}

## Fit sPLS-DA with either engine. Returns an object usable by
## lomics_spls_predict() and lomics_spls_selected().
lomics_spls_fit <- function(X, y, ncomp = 2, keepX = rep(ncol(X), ncomp),
                            engine = c("auto", "mixOmics", "internal")) {
  engine <- match.arg(engine)
  if (engine == "auto") engine <- if (requireNamespace("mixOmics", quietly = TRUE)) "mixOmics" else "internal"
  ncomp <- min(ncomp, ncol(X)); keepX <- rep_len(pmin(keepX, ncol(X)), ncomp)
  if (engine == "internal") return(.lomics_spls_internal(X, y, ncomp, keepX))
  m <- mixOmics::splsda(X, as.factor(y), ncomp = ncomp, keepX = keepX, scale = TRUE)
  A <- m$loadings$X; colnames(A) <- paste0("comp", seq_len(ncol(A)))
  structure(list(engine = "mixOmics", fit = m, loadings = A, variates = m$variates$X,
                 levels = levels(droplevels(as.factor(y))), ncomp = ncomp, keepX = keepX,
                 features = colnames(X)), class = "lomics_spls")
}

## Predicted class-indicator values (n_new x n_classes) using all components.
lomics_spls_predict <- function(fit, Xnew) {
  Xnew <- Xnew[, fit$features, drop = FALSE]
  if (fit$engine == "mixOmics") {
    pr <- stats::predict(fit$fit, newdata = Xnew)
    out <- pr$predict[, , fit$ncomp, drop = TRUE]
    if (is.null(dim(out))) out <- matrix(out, nrow = 1, dimnames = list(rownames(Xnew), fit$levels))
    colnames(out) <- fit$levels
    return(out)
  }
  m <- nrow(Xnew)
  Xs <- (Xnew - rep(fit$xm, each = m)) / rep(fit$xs, each = m)
  Yhat <- Xs %*% (fit$Wstar %*% t(fit$C))
  Yhat <- Yhat * rep(fit$ys, each = m) + rep(fit$ym, each = m)
  colnames(Yhat) <- fit$levels; rownames(Yhat) <- rownames(Xnew)
  Yhat
}

## Score for the second class level (used for AUC; higher = more like level 2)
lomics_spls_score <- function(fit, Xnew) {
  P <- lomics_spls_predict(fit, Xnew)
  P[, fit$levels[2]] - P[, fit$levels[1]]
}

## Selected variables with their loadings (all components)
lomics_spls_selected <- function(fit) {
  A <- fit$loadings
  sel <- rownames(A)[rowSums(A != 0) > 0]
  d <- data.frame(feature = sel, A[sel, , drop = FALSE], check.names = FALSE)
  d$importance <- apply(abs(A[sel, , drop = FALSE]), 1, max)
  d[order(-d$importance), , drop = FALSE]
}

## AUC (pROC if present, else Mann-Whitney identity)
lomics_auc <- function(y, score, positive) {
  yy <- as.integer(as.character(y) == positive)
  if (length(unique(yy)) < 2) return(NA_real_)
  if (requireNamespace("pROC", quietly = TRUE)) {
    return(as.numeric(pROC::auc(pROC::roc(yy, score, levels = c(0, 1), direction = "<", quiet = TRUE))))
  }
  r <- rank(score); n1 <- sum(yy); n0 <- sum(1 - yy)
  (sum(r[yy == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

lomics_auc_ci <- function(y, score, positive) {
  yy <- as.integer(as.character(y) == positive)
  if (!requireNamespace("pROC", quietly = TRUE) || length(unique(yy)) < 2) return(c(NA, NA, NA))
  r <- pROC::roc(yy, score, levels = c(0, 1), direction = "<", quiet = TRUE)
  as.numeric(pROC::ci.auc(r, method = "delong"))      # lower, auc, upper
}

## Paired accuracy: fraction of subjects whose later-visit sample scores
## higher than their earlier-visit sample (calibration-free, so it is valid
## when scores come from different fold models), with exact binomial 95% CI.
lomics_paired_accuracy <- function(y, score, subject, positive) {
  d <- vapply(split(seq_along(y), subject), function(ix) {
    if (length(ix) != 2) return(NA_real_)
    score[ix][as.character(y[ix]) == positive] - score[ix][as.character(y[ix]) != positive]
  }, numeric(1))
  d <- d[!is.na(d)]; k <- sum(d > 0); n <- length(d)
  ci <- stats::binom.test(k, n)$conf.int
  c(acc = k / n, ci_low = ci[1], ci_high = ci[2], n = n)
}

## Leave-one-subject-out predictions for a FIXED feature set. The multilevel
## decomposition is per subject and uses no labels, so it is applied to the
## full matrix once; scaling and the model are refit without the held-out
## subject.
lomics_spls_loso <- function(Xw, y, subject, ncomp = 2, keepX = NULL, engine = "auto") {
  subject <- as.character(subject); score <- rep(NA_real_, nrow(Xw))
  for (s in unique(subject)) {
    tr <- subject != s
    kx <- if (is.null(keepX)) rep(ncol(Xw), ncomp) else keepX
    f <- lomics_spls_fit(Xw[tr, , drop = FALSE], y[tr], ncomp, kx, engine)
    score[!tr] <- lomics_spls_score(f, Xw[!tr, , drop = FALSE])
  }
  score
}

## Sequential keepX tuning (component by component, as mixOmics::tune.splsda)
## using leave-one-subject-out AUC.
lomics_spls_tune <- function(Xw, y, subject, ncomp = 2, grid = c(5, 10, 20, 30, 50, 100),
                             engine = "auto", positive = levels(as.factor(y))[2]) {
  grid <- sort(unique(pmin(grid, ncol(Xw)))); chosen <- integer(); res <- list()
  for (h in seq_len(ncomp)) {
    aucs <- vapply(grid, function(k) {
      sc <- lomics_spls_loso(Xw, y, subject, ncomp = h, keepX = c(chosen, k), engine = engine)
      lomics_auc(y, sc, positive)
    }, numeric(1))
    res[[h]] <- data.frame(comp = h, keepX = grid, loso_auc = aucs)
    chosen <- c(chosen, grid[which.max(aucs)])   # ties -> smallest keepX
  }
  list(keepX = chosen, table = do.call(rbind, res))
}

## Backward elimination of a fixed panel (Vassileva et al. 2026, Methods 2.5):
## remove one feature at a time, least important first; a feature is retained
## only if its removal lowers the AUC by more than `drop_max` (relative by
## default) or would take the AUC below `target_auc`.
##   auc_on = "train" (as in the paper: resubstitution AUC of the full-data
##            model) or "loso" (leave-one-subject-out AUC; slower, less optimistic)
lomics_spls_backward <- function(Xw, y, subject, panel, ncomp = 2, drop_max = 0.05,
                                 drop_type = c("relative", "absolute"), target_auc = 0.9,
                                 auc_on = c("train", "loso"), engine = "auto",
                                 positive = levels(as.factor(y))[2]) {
  drop_type <- match.arg(drop_type); auc_on <- match.arg(auc_on)
  eval_auc <- function(f) {
    Xs <- Xw[, f, drop = FALSE]; nc <- min(ncomp, length(f))
    sc <- if (auc_on == "train") {
      lomics_spls_score(lomics_spls_fit(Xs, y, nc, rep(length(f), nc), engine), Xs)
    } else lomics_spls_loso(Xs, y, subject, nc, NULL, engine)
    lomics_auc(y, sc, positive)
  }
  imp <- function(f) {
    fit <- lomics_spls_fit(Xw[, f, drop = FALSE], y, min(ncomp, length(f)), rep(length(f), ncomp), engine)
    a <- abs(fit$loadings[, 1]); f[order(a)]
  }
  cur <- panel; base <- eval_auc(cur)
  log <- data.frame(step = 0, removed = NA, n_features = length(cur), auc = base, decision = "start")
  repeat {
    if (length(cur) <= 1) break
    removed <- FALSE
    for (f in imp(cur)) {
      new <- eval_auc(setdiff(cur, f))
      drop <- if (drop_type == "relative") (base - new) / base else base - new
      ok <- drop <= drop_max && (is.na(target_auc) || new >= target_auc)
      if (ok) {
        cur <- setdiff(cur, f)
        log <- rbind(log, data.frame(step = nrow(log), removed = f, n_features = length(cur),
                                     auc = new, decision = sprintf("removed (AUC change %+.3f)", new - base)))
        base <- new; removed <- TRUE; break
      }
    }
    if (!removed) break
  }
  list(panel = cur, auc = base, log = log)
}

## Full selection pipeline on one data set: tune keepX -> sPLS-DA ->
## backward elimination. Returns the fixed panel.
lomics_spls_select <- function(Xw, y, subject, ncomp = 2, grid = c(5, 10, 20, 30, 50, 100),
                               drop_max = 0.05, drop_type = "relative", target_auc = 0.9,
                               auc_on = "train", engine = "auto") {
  tn <- lomics_spls_tune(Xw, y, subject, ncomp, grid, engine)
  fit <- lomics_spls_fit(Xw, y, ncomp, tn$keepX, engine)
  start <- lomics_spls_selected(fit)$feature
  be <- lomics_spls_backward(Xw, y, subject, start, ncomp, drop_max, drop_type,
                             target_auc, auc_on, engine)
  list(tuning = tn, start_panel = start, panel = be$panel, train_auc = be$auc,
       elimination = be$log, engine = fit$engine)
}
