## =============================================================================
## lomics_missing.R — missing-value imputation and benchmarking helpers
## (shared by the prep-missing-data skill).
## All functions take a log-scale matrix x (features in rows, samples in
## columns) with NA for missing values and return the same shape.
## =============================================================================

## Half-minimum on the log2 scale: per-feature minimum - 1 (= half the
## minimum on the linear scale). Left-censoring (MNAR) assumption.
lomics_impute_halfmin <- function(x) {
  t(apply(x, 1, function(r) { if (all(is.na(r))) return(r); r[is.na(r)] <- min(r, na.rm = TRUE) - 1; r }))
}

## MinProb as described for imputeLCMD::impute.MinProb (Lazar et al. 2016):
## per sample, draw from N(q-th quantile of that sample, sigma), with sigma =
## median of the feature SDs x tune_sigma. Random: seed it and report it.
lomics_impute_minprob <- function(x, q = 0.01, tune_sigma = 1, seed = 1) {
  set.seed(seed)
  sig <- stats::median(apply(x, 1, stats::sd, na.rm = TRUE), na.rm = TRUE) * tune_sigma
  for (j in seq_len(ncol(x))) {
    m <- is.na(x[, j]); if (!any(m)) next
    x[m, j] <- stats::rnorm(sum(m), stats::quantile(x[, j], q, na.rm = TRUE), sig)
  }
  x
}

## Feature mean (MCAR baseline; shrinks variance — for comparison only)
lomics_impute_mean <- function(x) {
  t(apply(x, 1, function(r) { r[is.na(r)] <- mean(r, na.rm = TRUE); r }))
}

## k-nearest neighbours over features (impute::impute.knn; MAR assumption)
lomics_impute_knn <- function(x, k = 10) {
  if (!requireNamespace("impute", quietly = TRUE)) stop("Package 'impute' (Bioconductor) is required for knn.", call. = FALSE)
  out <- utils::capture.output(r <- impute::impute.knn(x, k = k, rowmax = 0.8, colmax = 0.8)$data)
  r
}

## Bayesian PCA (pcaMethods::pca(method = "bpca"); MAR assumption)
lomics_impute_bpca <- function(x, npcs = 5) {
  if (!requireNamespace("pcaMethods", quietly = TRUE)) stop("Package 'pcaMethods' (Bioconductor) is required for bpca.", call. = FALSE)
  keep <- rowSums(!is.na(x)) > 0
  p <- pcaMethods::pca(t(x[keep, , drop = FALSE]), method = "bpca", nPcs = min(npcs, ncol(x) - 1), verbose = FALSE)
  x[keep, ] <- t(pcaMethods::completeObs(p))
  x
}

lomics_impute <- function(x, method = c("halfmin", "minprob", "knn", "bpca", "mean", "none"), seed = 1) {
  method <- match.arg(method)
  switch(method, halfmin = lomics_impute_halfmin(x), minprob = lomics_impute_minprob(x, seed = seed),
         knn = lomics_impute_knn(x), bpca = lomics_impute_bpca(x), mean = lomics_impute_mean(x), none = x)
}

## Masking benchmark: hide a fraction of OBSERVED values under a given
## mechanism, impute, and score against the hidden truth.
##   mechanism "MCAR": uniformly at random
##   mechanism "MNAR": left-censoring — values below the feature's
##             `lod_q` quantile are masked with high probability
## Returns RMSE, NRMSE (RMSE / SD of all observed values), bias (mean imputed
## - true) and Spearman correlation per method.
lomics_mask_benchmark <- function(x, methods, frac = 0.05, mechanism = c("MCAR", "MNAR"),
                                  lod_q = 0.2, reps = 3, seed = 1) {
  mechanism <- match.arg(mechanism)
  res <- list()
  for (r in seq_len(reps)) {
    set.seed(seed + r)
    obs <- which(!is.na(x))
    n_mask <- round(frac * length(obs))
    if (mechanism == "MCAR") {
      idx <- sample(obs, n_mask)
    } else {
      thr <- apply(x, 1, stats::quantile, lod_q, na.rm = TRUE)
      low <- obs[x[obs] <= thr[(obs - 1) %% nrow(x) + 1]]
      idx <- if (length(low) > n_mask) sample(low, n_mask) else low
    }
    truth <- x[idx]; xm <- x; xm[idx] <- NA; sd_all <- stats::sd(x, na.rm = TRUE)
    for (m in methods) {
      imp <- tryCatch(lomics_impute(xm, m, seed = seed + r), error = function(e) NULL)
      if (is.null(imp)) next
      est <- imp[idx]; ok <- !is.na(est)
      res[[length(res) + 1]] <- data.frame(mechanism = mechanism, rep = r, method = m, n_masked = length(idx),
        rmse = sqrt(mean((est[ok] - truth[ok])^2)),
        nrmse = sqrt(mean((est[ok] - truth[ok])^2)) / sd_all,
        bias = mean(est[ok] - truth[ok]),
        spearman = suppressWarnings(stats::cor(est[ok], truth[ok], method = "spearman")))
    }
  }
  do.call(rbind, res)
}
