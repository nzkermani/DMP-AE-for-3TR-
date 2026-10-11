#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## run_longomicstraj.R — Omics Trajectory Signatures (OTS) for each omic layer
## with LongOmicsTraj (github.com/nzkermani/LongOmicsTraj).
##
## Usage (from the project root):
##   Rscript .claude/skills/cluster-trajectories-longomicstraj/scripts/run_longomicstraj.R \
##       [config] [layer|all] [estimators] [threshold_method] [delta_min] [k]
##   estimators       comma list from empirical,lmm,gam,masigpro
##                    (default "empirical,lmm,gam"; masigpro needs maSigPro)
##   threshold_method fixed | global_mad | gene_mad | hybrid | quantile (default hybrid)
##   delta_min        minimum |delta| on the log2 scale (default 0.15)
##   k                MAD multiplier (default 0.5)
##
## IMPORTANT design choice: one LongOmicsTraj object per omic layer.
## In LongOmicsTraj 0.99.0, lot_from_*() stops when a VMM already exists unless
## overwrite = TRUE, and overwrite = TRUE replaces the VMM for ALL assays, so
## several assays cannot be accumulated in one object. Layers are therefore
## processed separately and their OTS tables combined afterwards.
##
## Outputs (results/trajectories/):
##   <layer>/lot_<estimator>.rds, <layer>/ots_all_estimators.csv,
##   <layer>/consensus.csv, <layer>/summary_<estimator>.txt,
##   <layer>/ots_distribution.pdf
##   ots_all_layers.csv              consensus OTS for every layer (long format)
##   ots_layer_group_proportions.csv topology proportions per layer x group
##   discordant_between_groups.csv   features whose consensus OTS differs by arm
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(LongOmicsTraj)
  library(SummarizedExperiment)
  library(ggplot2)
})
source("R/lomics_utils.R")

args <- commandArgs(trailingOnly = TRUE)
a <- function(i, d) if (length(args) >= i && nzchar(args[i])) args[i] else d
cfg        <- lomics_read_config(a(1, "study_config.json"))
layers     <- if (a(2, "all") == "all") lomics_layer_names(cfg) else a(2, "all")
estimators <- strsplit(a(3, "empirical,lmm,gam"), ",")[[1]]
th_method  <- a(4, "hybrid")
delta_min  <- as.numeric(a(5, "0.15"))
k_mad      <- as.numeric(a(6, "0.5"))
visits     <- cfg$visits

if ("masigpro" %in% estimators && !requireNamespace("maSigPro", quietly = TRUE)) {
  message("maSigPro not installed — dropping the masigpro estimator.")
  estimators <- setdiff(estimators, "masigpro")
}

finish <- function(obj, layer) {
  obj <- lot_compute_deltas(obj, assay = layer, visits = visits, overwrite = TRUE)
  obj <- lot_compute_thresholds(obj, assay = layer, method = th_method,
                                delta_min = delta_min, k = k_mad,
                                visits = visits, overwrite = TRUE, verbose = FALSE)
  lot_map_to_ots(obj, visits = visits, overwrite = TRUE)
}

all_layers <- list(); discord <- list()
for (layer in layers) {
  se <- lomics_load_layer(cfg, layer, "processed")
  x <- assay(se, "values"); meta <- as.data.frame(colData(se))
  if (anyNA(x)) stop(layer, ": NA values remain; impute in the QC step first.")
  L <- cfg$layers[cfg$layers$name == layer, ]
  ftype <- if (nrow(L) == 1) L$feature_type else "pathway"   # derived layers
  base <- lot_create(expression = x, metadata = meta, assay = layer,
                     visit_levels = visits, time_values = cfg$time_values,
                     feature_type = ftype)
  outdir <- dirname(lomics_out(cfg, "trajectories", layer, "x"))
  fits <- list()
  for (est in estimators) {
    message(sprintf("[%s] estimator: %s", layer, est))
    obj <- tryCatch(switch(est,
      empirical = lot_from_empirical(base, assay = layer, visits = visits, overwrite = TRUE, verbose = FALSE),
      lmm = lot_from_lmm(base, assay = layer,
                         formula_fixed = y ~ visit * group + (1 | subject_id),
                         visits = visits, overwrite = TRUE, verbose = FALSE),
      gam = lot_from_gam(base, assay = layer,
                         formula_fixed = y ~ group +
                           mgcv::s(time, by = group, k = min(3, length(visits) - 1), bs = "cr") +
                           mgcv::s(subject_id, bs = "re"),
                         visits = visits, overwrite = TRUE, verbose = FALSE),
      masigpro = lot_from_masigpro(base, assay = layer, visits = visits, overwrite = TRUE, verbose = FALSE),
      stop("Unknown estimator ", est)),
      error = function(e) { message("  failed: ", conditionMessage(e)); NULL })
    if (is.null(obj)) next
    obj <- finish(obj, layer)
    fits[[est]] <- obj
    saveRDS(obj, file.path(outdir, paste0("lot_", est, ".rds")))
    writeLines(utils::capture.output(print(lot_summary(obj, top_n = 8L))),
               file.path(outdir, paste0("summary_", est, ".txt")))
  }
  if (!length(fits)) stop(layer, ": every estimator failed.")

  ots_long <- do.call(rbind, lapply(names(fits), function(est) {
    o <- as.data.frame(lot_ots(fits[[est]]))
    data.frame(layer = layer, feature = o$object_id, group = o$group,
               estimator = est, topology = o$topology_label)
  }))
  utils::write.csv(ots_long, file.path(outdir, "ots_all_estimators.csv"), row.names = FALSE)

  if (length(fits) >= 2) {
    cmp <- do.call(lot_compare_estimators, c(unname(fits), list(.names = names(fits))))
    ## lot_consensus() (v0.99.0) votes across estimator x group columns, i.e.
    ## it pools treatment arms. Run it within each group so arms stay separate.
    cons_tab <- do.call(rbind, lapply(unique(as.character(cmp$group)), function(g) {
      cg <- as.data.frame(lot_consensus(cmp[as.character(cmp$group) == g, ],
                                        min_agreement = 2L, mode = "majority"))
      data.frame(layer = layer, feature = cg$object_id, group = g,
                 topology = cg$consensus_topology,
                 agreement_count = cg$agreement_count, n_estimators = cg$n_available)
    }))
  } else {
    cons_tab <- transform(ots_long, estimator = NULL)
  }
  utils::write.csv(cons_tab, file.path(outdir, "consensus.csv"), row.names = FALSE)
  all_layers[[layer]] <- cons_tab[, c("layer", "feature", "group", "topology")]

  ## Features whose consensus topology differs between arms
  grp <- unique(stats::na.omit(cons_tab$group))
  if (length(grp) >= 2) {
    w <- stats::reshape(cons_tab[, c("feature", "group", "topology")], idvar = "feature",
                        timevar = "group", direction = "wide")
    tcols <- grep("^topology\\.", names(w), value = TRUE)
    w$n_distinct <- apply(w[, tcols, drop = FALSE], 1, function(r) length(unique(stats::na.omit(r))))
    discord[[layer]] <- data.frame(layer = layer, w[w$n_distinct > 1, , drop = FALSE], check.names = FALSE)
  }

  pdf(file.path(outdir, "ots_distribution.pdf"), width = 10, height = 6)
  for (est in names(fits)) {
    p <- tryCatch(lot_plot_ots_distribution(fits[[est]], type = "proportion") +
                    ggtitle(sprintf("%s — %s", layer, est)), error = function(e) NULL)
    if (!is.null(p)) print(p)
  }
  invisible(dev.off())
}

ots_all <- do.call(rbind, all_layers)

## Significance gating: OTS thresholds describe the SHAPE of a trajectory but
## not whether the change exceeds sampling noise. If the longitudinal-
## differential step has been run, attach the within-group overall-time FDR so
## shapes of non-changing features can be set aside.
ots_all$time_fdr <- NA_real_
for (layer in unique(ots_all$layer)) for (g in unique(ots_all$group)) {
  f <- file.path(cfg$output_dir, "de", layer, sprintf("overall_time_%s.csv", make.names(g)))
  if (file.exists(f)) {
    de <- utils::read.csv(f)
    i <- ots_all$layer == layer & ots_all$group == g
    ots_all$time_fdr[i] <- de$adj.P.Val[match(ots_all$feature[i], de$feature)]
  }
}
ots_all$topology_gated <- ifelse(is.na(ots_all$time_fdr), NA,
                                 ifelse(ots_all$time_fdr < 0.05, ots_all$topology, "no_time_effect"))
## Merge with layers processed in earlier runs (re-running one layer must not
## wipe the others)
prev_f <- lomics_out(cfg, "trajectories", "ots_all_layers.csv")
if (file.exists(prev_f)) {
  prev <- utils::read.csv(prev_f)
  prev <- prev[!prev$layer %in% unique(ots_all$layer), intersect(names(prev), names(ots_all)), drop = FALSE]
  for (cl in setdiff(names(ots_all), names(prev))) prev[[cl]] <- NA
  ots_all <- rbind(prev[, names(ots_all)], ots_all)
}
utils::write.csv(ots_all, prev_f, row.names = FALSE)
props <- stats::aggregate(feature ~ layer + group + topology, data = ots_all, FUN = length)
names(props)[4] <- "n"
props$proportion <- stats::ave(props$n, props$layer, props$group, FUN = function(n) n / sum(n))
props <- props[order(props$layer, props$group, -props$n), ]
utils::write.csv(props, lomics_out(cfg, "trajectories", "ots_layer_group_proportions.csv"), row.names = FALSE)
if (length(discord)) {
  d <- do.call(rbind, lapply(discord, function(z) { names(z) <- sub("^topology\\.", "", names(z)); z }))
  dfile <- lomics_out(cfg, "trajectories", "discordant_between_groups.csv")
  if (file.exists(dfile)) {
    prev <- utils::read.csv(dfile, check.names = FALSE)
    prev <- prev[!prev$layer %in% unique(d$layer), , drop = FALSE]
    if (identical(names(prev), names(d))) d <- rbind(prev, d)
  }
  utils::write.csv(d, dfile, row.names = FALSE)
}

p <- ggplot(props, aes(x = stats::reorder(topology, -proportion), y = proportion, fill = layer)) +
  geom_col(position = "dodge") + facet_wrap(~ group, ncol = 1) +
  labs(x = "Consensus OTS", y = "Proportion of features",
       title = "Trajectory topology composition across omic layers") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
ggsave(lomics_out(cfg, "trajectories", "ots_layers_comparison.pdf"), p, width = 12, height = 7)
print(utils::head(props, 20), row.names = FALSE)
lomics_session_log(cfg, "longomicstraj")
