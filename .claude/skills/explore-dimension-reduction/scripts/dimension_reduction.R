#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## dimension_reduction.R — PCA, t-SNE and UMAP for one or more omic layers.
##
## Usage (from the project root):
##   Rscript .claude/skills/explore-dimension-reduction/scripts/dimension_reduction.R \
##     --layers=proteomics[,metabolomics|all] --visit=V1|all [--change_from=V1] \
##     [--methods=pca,tsne,umap] [--colour_by=group,visit,batch,responder] \
##     [--top_var=2000] [--perplexity=auto] [--n_neighbors=15] [--seed=42]
##
## Multi-layer input: each layer is autoscaled and then divided by
## sqrt(n_features) (block scaling) so every layer has equal total weight.
## visit=all plots every sample and joins each subject's visits with a path.
##
## Outputs (results/dimred/<tag>/):
##   pca_scores.csv, pca_variance.csv, pca_loadings.csv (top 50 per PC),
##   pc_associations.csv (PC1-10 vs metadata: lm / ANOVA p-values),
##   tsne_coords.csv, umap_coords.csv, plots.pdf
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({ library(SummarizedExperiment); library(ggplot2) })
source("R/lomics_utils.R")
a <- lomics_args(list(config = "study_config.json", layers = "all", visit = NULL,
                      change_from = NULL, methods = "pca,tsne,umap", colour_by = NULL,
                      top_var = "2000", perplexity = "auto", n_neighbors = "15", seed = NULL))
cfg <- lomics_read_config(a$config)
layers <- lomics_layers_arg(cfg, a$layers)
visit <- if (is.null(a$visit)) cfg$default_visit else a$visit
seed <- as.integer(if (is.null(a$seed)) cfg$seed else a$seed)
methods <- strsplit(a$methods, ",")[[1]]
d <- lomics_xy(cfg, layers, visit, a$change_from, top_var = as.integer(a$top_var))
X <- d$X; meta <- d$meta
X <- X[, apply(X, 2, stats::var, na.rm = TRUE) > 0, drop = FALSE]
if (anyNA(X)) stop("NA values in the matrix; impute in prep-qc-normalisation first.")

## autoscale + block scaling
Xs <- scale(X)
if (length(layers) > 1) {
  fl <- d$feature_layer[colnames(Xs)]
  for (l in unique(fl)) Xs[, fl == l] <- Xs[, fl == l] / sqrt(sum(fl == l))
}

colour_by <- if (!is.null(a$colour_by)) strsplit(a$colour_by, ",")[[1]] else
  unique(c("group", if (identical(visit, "all")) "visit", cfg$batch_column,
           unlist(cfg$outcomes[c("binary", "categorical")])))
colour_by <- intersect(colour_by, colnames(meta))
outdir <- dirname(lomics_out(cfg, "dimred", lomics_tag(layers, visit, a$change_from), "x"))
message(d$description)

plot_embed <- function(coords, xlab, ylab, title) {
  df <- data.frame(coords[, 1:2], meta, check.names = FALSE)
  names(df)[1:2] <- c("D1", "D2")
  lapply(colour_by, function(cl) {
    p <- ggplot(df, aes(D1, D2))
    if (identical(visit, "all")) {
      df <- df[order(df$subject_id, df$visit), ]
      p <- p + geom_path(data = df, aes(group = subject_id), colour = "grey80", linewidth = .3)
    }
    p + geom_point(aes(colour = .data[[cl]]), size = 2, alpha = .85) +
      labs(x = xlab, y = ylab, colour = cl, title = sprintf("%s — %s", title, d$description)) +
      theme_minimal(base_size = 10)
  })
}
plots <- list()

if ("pca" %in% methods) {
  pc <- stats::prcomp(Xs, center = FALSE, scale. = FALSE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  npc <- min(10, ncol(pc$x))
  utils::write.csv(data.frame(observation = rownames(pc$x), pc$x[, 1:npc]),
                   file.path(outdir, "pca_scores.csv"), row.names = FALSE)
  utils::write.csv(data.frame(PC = seq_along(ve), variance = ve, cumulative = cumsum(ve))[1:npc, ],
                   file.path(outdir, "pca_variance.csv"), row.names = FALSE)
  load <- do.call(rbind, lapply(1:min(5, npc), function(k) {
    l <- pc$rotation[, k]; o <- order(abs(l), decreasing = TRUE)[1:min(50, length(l))]
    data.frame(PC = k, feature = names(l)[o], loading = l[o], layer = d$feature_layer[names(l)[o]])
  }))
  utils::write.csv(load, file.path(outdir, "pca_loadings.csv"), row.names = FALSE)
  vars <- setdiff(colnames(meta), c("sample_id", "subject_id"))
  vars <- vars[vapply(vars, function(v) { u <- length(unique(stats::na.omit(meta[[v]]))); u > 1 && u < nrow(meta) }, logical(1))]
  assoc <- do.call(rbind, lapply(vars, function(v) {
    sapply(1:npc, function(k) tryCatch(stats::anova(stats::lm(pc$x[, k] ~ meta[[v]]))[["Pr(>F)"]][1],
                                       error = function(e) NA))
  }))
  dimnames(assoc) <- list(vars, paste0("PC", 1:npc))
  utils::write.csv(signif(assoc, 3), file.path(outdir, "pc_associations.csv"))
  plots <- c(plots,
             list(ggplot(data.frame(PC = 1:npc, v = ve[1:npc]), aes(factor(PC), v)) +
                    geom_col(fill = "grey40") + labs(x = "PC", y = "Variance explained", title = "Scree") +
                    theme_minimal()),
             plot_embed(pc$x, sprintf("PC1 (%.1f%%)", 100 * ve[1]), sprintf("PC2 (%.1f%%)", 100 * ve[2]), "PCA"))
  message("PCA: PC1-PC3 explain ", paste(round(100 * ve[1:3], 1), collapse = "/"), "%")
}

## t-SNE / UMAP on the top PCs (standard denoising step)
pc_in <- stats::prcomp(Xs, center = FALSE)$x[, 1:min(30, nrow(Xs) - 1, ncol(Xs)), drop = FALSE]

if ("tsne" %in% methods) {
  if (!requireNamespace("Rtsne", quietly = TRUE)) message("Rtsne not installed — skipping t-SNE.") else {
    perp <- if (a$perplexity == "auto") max(2, min(30, floor((nrow(pc_in) - 1) / 3) - 1)) else as.numeric(a$perplexity)
    set.seed(seed)
    ts <- Rtsne::Rtsne(pc_in, dims = 2, perplexity = perp, pca = FALSE, check_duplicates = FALSE)$Y
    rownames(ts) <- rownames(Xs)
    utils::write.csv(data.frame(observation = rownames(ts), tSNE1 = ts[, 1], tSNE2 = ts[, 2]),
                     file.path(outdir, "tsne_coords.csv"), row.names = FALSE)
    plots <- c(plots, plot_embed(ts, "tSNE1", "tSNE2", sprintf("t-SNE (perplexity %s)", perp)))
  }
}
if ("umap" %in% methods) {
  if (!requireNamespace("uwot", quietly = TRUE)) message("uwot not installed — skipping UMAP.") else {
    nn <- min(as.integer(a$n_neighbors), nrow(pc_in) - 1)
    set.seed(seed)
    um <- uwot::umap(pc_in, n_neighbors = nn, min_dist = 0.1, n_threads = 1, seed = seed)
    rownames(um) <- rownames(Xs)
    utils::write.csv(data.frame(observation = rownames(um), UMAP1 = um[, 1], UMAP2 = um[, 2]),
                     file.path(outdir, "umap_coords.csv"), row.names = FALSE)
    plots <- c(plots, plot_embed(um, "UMAP1", "UMAP2", sprintf("UMAP (n_neighbors %d)", nn)))
  }
}
pdf(file.path(outdir, "plots.pdf"), width = 8, height = 6)
for (p in plots) print(p)
invisible(dev.off())
if (exists("assoc")) { cat("\nPC associations (p-values):\n"); print(signif(assoc[, 1:min(5, ncol(assoc))], 2)) }
message("Outputs in ", outdir)
lomics_session_log(cfg, "dimred")
