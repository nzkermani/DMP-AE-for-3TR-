#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## cluster_samples.R — unsupervised clustering of subjects/samples
## (endotype / phenotype discovery) with k selection, stability and clinical
## characterisation.
##
## Usage (from the project root):
##   Rscript .claude/skills/cluster-samples/scripts/cluster_samples.R \
##     --layers=all --visit=V1 [--change_from=V1] \
##     [--methods=kmeans,pam,hclust,gmm,consensus] [--k=2:6] [--top_var=1000] \
##     [--n_pcs=10] [--boot=50] [--seed=42]
##
## Pipeline: autoscale (+ block scaling for >1 layer) -> top PCs ->
##   k-means (25 starts) | PAM | Ward hierarchical | Gaussian mixture (mclust, on PCs)
##   | ConsensusClusterPlus (if installed)
## k chosen per method by mean silhouette (GMM: BIC). Stability: bootstrap
## Jaccard similarity of each cluster (clusterboot-style, resampling subjects).
## Agreement between methods: adjusted Rand index.
##
## Outputs (results/clustering/<tag>/):
##   k_selection.csv, assignments.csv (one column per method + "cluster" =
##   recommended solution), stability.csv, method_agreement_ARI.csv,
##   clinical_association.csv, top_features_by_cluster.csv, plots.pdf
## Use assignments.csv with stats-survival --mode=km --strata_file=...
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(SummarizedExperiment); library(cluster); library(ggplot2) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", layers = "all", visit = NULL, change_from = NULL,
                      methods = "kmeans,pam,hclust,gmm,consensus", k = "2:6", top_var = "1000",
                      n_pcs = "10", boot = "50", seed = NULL))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
seed <- as.integer(if (is.null(a$seed)) cfg$seed else a$seed)
ks <- eval(parse(text = a$k)); methods <- strsplit(a$methods, ",")[[1]]
if (identical(visit, "all")) message("NOTE: visit=all clusters samples, not subjects; the same subject at different visits will tend to co-cluster.")

d <- lomics_xy(cfg, layers, visit, a$change_from, top_var = as.integer(a$top_var))
X <- d$X[, apply(d$X, 2, stats::var) > 0, drop = FALSE]; meta <- d$meta
Xs <- scale(X)
if (length(layers) > 1) { fl <- d$feature_layer[colnames(Xs)]; for (l in unique(fl)) Xs[, fl == l] <- Xs[, fl == l] / sqrt(sum(fl == l)) }
npc <- min(as.integer(a$n_pcs), nrow(Xs) - 1, ncol(Xs))
pc <- stats::prcomp(Xs, center = FALSE); Z <- pc$x[, 1:npc, drop = FALSE]
D <- stats::dist(Z)
tag <- lomics_tag(layers, visit, a$change_from)
outdir <- dirname(lomics_out(cfg, "clustering", tag, "x"))
message(d$description, sprintf(" | clustering on %d PCs (%.0f%% variance)", npc, 100 * sum(pc$sdev[1:npc]^2) / sum(pc$sdev^2)))

cl_fun <- list(
  kmeans = function(Z, k) { set.seed(seed); stats::kmeans(Z, k, nstart = 25)$cluster },
  pam    = function(Z, k) cluster::pam(Z, k)$clustering,
  hclust = function(Z, k) stats::cutree(stats::hclust(stats::dist(Z), "ward.D2"), k),
  gmm    = function(Z, k) mclust::Mclust(Z, G = k, verbose = FALSE)$classification)
## Mclust() calls mclustBIC() unqualified, so mclust must be attached
if ("gmm" %in% methods && requireNamespace("mclust", quietly = TRUE)) suppressPackageStartupMessages(library(mclust))
sil <- function(cl) if (length(unique(cl)) < 2) NA else mean(cluster::silhouette(cl, D)[, 3])

ksel <- list(); best <- list()
for (m in intersect(methods, names(cl_fun))) {
  if (m == "gmm" && !requireNamespace("mclust", quietly = TRUE)) { message("mclust not installed — skipping gmm"); next }
  sols <- lapply(ks, function(k) tryCatch(cl_fun[[m]](Z, k), error = function(e) NULL))
  s <- vapply(sols, function(cl) if (is.null(cl)) NA else sil(cl), numeric(1))
  bic <- if (m == "gmm") vapply(ks, function(k) tryCatch(mclust::Mclust(Z, G = k, verbose = FALSE)$bic, error = function(e) NA), numeric(1)) else NA
  ksel[[m]] <- data.frame(method = m, k = ks, mean_silhouette = s, BIC = bic)
  crit <- if (m == "gmm") bic else s
  if (all(is.na(crit))) { message(m, ": no valid solution — skipped"); next }
  best[[m]] <- sols[[which.max(crit)]]
}
if ("consensus" %in% methods) {
  if (requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    cc <- ConsensusClusterPlus::ConsensusClusterPlus(t(Z), maxK = max(ks), reps = 200, pItem = 0.8, pFeature = 1,
                                                     clusterAlg = "km", distance = "euclidean", seed = seed,
                                                     title = file.path(outdir, "consensus"), plot = "pdf", verbose = FALSE)
    ## PAC (proportion of ambiguous clustering), lower = more stable
    pac <- vapply(ks, function(k) { cm <- cc[[k]]$consensusMatrix; v <- cm[lower.tri(cm)]; mean(v > 0.1 & v < 0.9) }, numeric(1))
    s <- vapply(ks, function(k) sil(cc[[k]]$consensusClass), numeric(1))
    ksel$consensus <- data.frame(method = "consensus", k = ks, mean_silhouette = s, BIC = NA, PAC = pac)
    best$consensus <- cc[[ks[which.min(pac)]]]$consensusClass
  } else message("ConsensusClusterPlus not installed — skipping consensus")
}
ksel <- do.call(rbind, lapply(ksel, function(z) { if (!"PAC" %in% names(z)) z$PAC <- NA; z }))
utils::write.csv(ksel, file.path(outdir, "k_selection.csv"), row.names = FALSE)

## Bootstrap stability (Jaccard) of each method's chosen solution
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
set.seed(seed); B <- as.integer(a$boot)
stab <- do.call(rbind, lapply(names(best), function(m) {
  if (m == "consensus") return(NULL)
  cl0 <- best[[m]]; k <- length(unique(cl0))
  J <- matrix(NA, B, k)
  for (b in seq_len(B)) {
    idx <- sort(unique(sample(nrow(Z), replace = TRUE)))
    clb <- tryCatch(cl_fun[[m]](Z[idx, , drop = FALSE], k), error = function(e) NULL)
    if (is.null(clb)) next
    for (c0 in seq_len(k)) {
      orig <- intersect(which(cl0 == c0), idx)
      if (!length(orig)) next
      J[b, c0] <- max(vapply(unique(clb), function(cb) jacc(orig, idx[clb == cb]), numeric(1)))
    }
  }
  data.frame(method = m, cluster = seq_len(k), size = as.vector(table(factor(cl0, seq_len(k)))),
             mean_jaccard = colMeans(J, na.rm = TRUE))
}))
utils::write.csv(stab, file.path(outdir, "stability.csv"), row.names = FALSE)

## Agreement between methods
ari <- function(x, y) {  # adjusted Rand index
  tb <- table(x, y); n <- sum(tb); s <- function(v) sum(choose(v, 2))
  e <- s(rowSums(tb)) * s(colSums(tb)) / choose(n, 2)
  (s(tb) - e) / ((s(rowSums(tb)) + s(colSums(tb))) / 2 - e)
}
M <- outer(names(best), names(best), Vectorize(function(i, j) ari(best[[i]], best[[j]])))
dimnames(M) <- list(names(best), names(best))
utils::write.csv(round(M, 3), file.path(outdir, "method_agreement_ARI.csv"))

## Recommended solution: most stable method (mean Jaccard), ties -> silhouette
rec <- if (!is.null(stab) && nrow(stab)) {
  ms <- stats::aggregate(mean_jaccard ~ method, stab, mean); ms$method[which.max(ms$mean_jaccard)]
} else names(best)[1]
assign_df <- data.frame(observation = rownames(Z), meta[, intersect(c("subject_id", "sample_id", "visit"), names(meta)), drop = FALSE],
                        do.call(cbind, lapply(best, as.integer)), cluster = paste0("C", best[[rec]]), check.names = FALSE)
utils::write.csv(assign_df, file.path(outdir, "assignments.csv"), row.names = FALSE)

## Clinical association of the recommended clusters
cl <- factor(assign_df$cluster)
vars <- setdiff(colnames(meta), c("sample_id", "subject_id", "visit", "time"))
assoc <- do.call(rbind, lapply(vars, function(v) {
  x <- meta[[v]]; if (length(unique(stats::na.omit(x))) < 2) return(NULL)
  if (is.numeric(x) && length(unique(x)) > 2) {
    data.frame(variable = v, test = "Kruskal-Wallis", p = stats::kruskal.test(x ~ cl)$p.value,
               summary = paste(sprintf("%s: %.2f", levels(cl), tapply(x, cl, stats::median, na.rm = TRUE)), collapse = "; "))
  } else {
    tb <- table(x, cl)
    data.frame(variable = v, test = "Fisher", p = stats::fisher.test(tb, simulate.p.value = TRUE, B = 5000)$p.value,
               summary = paste(apply(tb, 2, function(z) paste(z, collapse = "/")), collapse = "; "))
  }
}))
if (!is.null(cfg$outcomes$survival_time) && all(c(cfg$outcomes$survival_time, cfg$outcomes$survival_event) %in% names(meta))) {
  sv <- survival::survdiff(survival::Surv(meta[[cfg$outcomes$survival_time]], meta[[cfg$outcomes$survival_event]]) ~ cl)
  assoc <- rbind(assoc, data.frame(variable = "survival", test = "log-rank",
                                   p = stats::pchisq(sv$chisq, length(sv$n) - 1, lower.tail = FALSE), summary = ""))
}
assoc$fdr <- stats::p.adjust(assoc$p, "BH")
utils::write.csv(assoc, file.path(outdir, "clinical_association.csv"), row.names = FALSE)

## Features driving clusters (one-vs-rest Wilcoxon on original features)
topf <- do.call(rbind, lapply(levels(cl), function(c0) {
  p <- apply(X, 2, function(x) stats::wilcox.test(x[cl == c0], x[cl != c0], exact = FALSE)$p.value)
  dif <- colMeans(X[cl == c0, , drop = FALSE]) - colMeans(X[cl != c0, , drop = FALSE])
  r <- data.frame(cluster = c0, feature = colnames(X), layer = d$feature_layer[colnames(X)], diff_vs_rest = dif, p = p, fdr = stats::p.adjust(p, "BH"))
  utils::head(r[order(r$p), ], 30)
}))
utils::write.csv(topf, file.path(outdir, "top_features_by_cluster.csv"), row.names = FALSE)

pdf(file.path(outdir, "plots.pdf"), width = 8, height = 6)
print(ggplot(ksel, aes(k, mean_silhouette, colour = method)) + geom_line() + geom_point() +
        labs(title = "Mean silhouette by k", y = "Mean silhouette") + theme_minimal())
pdd <- data.frame(PC1 = Z[, 1], PC2 = Z[, 2], cluster = cl)
print(ggplot(pdd, aes(PC1, PC2, colour = cluster)) + geom_point(size = 2) +
        labs(title = sprintf("Recommended: %s, k = %d", rec, nlevels(cl))) + theme_minimal())
plot(cluster::silhouette(as.integer(cl), D), main = "Silhouette (recommended solution)", border = NA)
top_feat <- unique(topf$feature[topf$fdr < 0.05])[1:min(50, sum(topf$fdr < 0.05))]
if (length(top_feat) >= 2 && requireNamespace("pheatmap", quietly = TRUE)) {
  o <- order(cl)
  pheatmap::pheatmap(t(scale(X[o, stats::na.omit(top_feat), drop = FALSE])), cluster_cols = FALSE, show_colnames = FALSE,
                     annotation_col = data.frame(cluster = cl[o], row.names = rownames(X)[o]),
                     breaks = seq(-3, 3, length.out = 101), main = "Top cluster-discriminating features")
}
invisible(dev.off())

cat("k selection:\n"); print(ksel, row.names = FALSE)
cat("\nStability (bootstrap Jaccard; >0.75 stable, <0.6 dissolved):\n"); print(stab, row.names = FALSE)
cat("\nMethod agreement (ARI):\n"); print(round(M, 2))
cat(sprintf("\nRecommended: %s (k = %d)\n", rec, nlevels(cl)))
cat("\nClinical associations:\n"); print(assoc[order(assoc$p), ], row.names = FALSE)
message("Outputs in ", outdir)
lomics_session_log(cfg, "clustering")
