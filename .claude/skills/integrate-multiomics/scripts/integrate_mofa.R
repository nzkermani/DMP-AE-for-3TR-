#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## integrate_mofa.R — unsupervised multi-omics integration across 4 visits with
## MOFA2, followed by repeated-measures testing of each latent factor.
##
## Usage (from the project root):
##   Rscript .claude/skills/integrate-multiomics/scripts/integrate_mofa.R \
##       [config] [n_factors] [top_var_features_per_layer] [layers]
##   defaults: 15 factors, 2000 most variable features per layer, all layers
##
## Design choices
##   * All samples (subject x visit) enter one MOFA model; MOFA tolerates
##     samples missing a whole layer, so incomplete follow-up is kept.
##   * Each layer is reduced to its most variable features so large layers
##     (RNA-seq) do not dominate; views are scaled to unit variance.
##   * Factors are then tested with  Factor ~ visit * group + covariates +
##     (1 | subject_id)  (lme4, likelihood-ratio tests), giving time, group and
##     time x group effects per factor.
##
## Outputs (results/integration/mofa/):
##   mofa_model.hdf5, factors.csv, weights_<layer>.csv,
##   variance_explained.csv, factor_tests.csv, factor_trajectories.pdf
##
## Status: written against the MOFA2 API (create_mofa / prepare_mofa /
## run_mofa); NOT executed in the environment this kit was built in, because
## MOFA2 was unavailable there. Run on a small subset first.
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(MOFA2)
  library(lme4)
  library(ggplot2)
})
source("R/lomics_utils.R")

args <- commandArgs(trailingOnly = TRUE)
a <- function(i, d) if (length(args) >= i && nzchar(args[i])) args[i] else d
cfg     <- lomics_read_config(a(1, "study_config.json"))
n_fac   <- as.integer(a(2, "15"))
top_var <- as.integer(a(3, "2000"))
layers  <- if (a(4, "all") == "all") lomics_layer_names(cfg) else strsplit(a(4, ""), ",")[[1]]
outdir  <- dirname(lomics_out(cfg, "integration", "mofa", "x"))

long <- list(); meta_all <- NULL
for (layer in layers) {
  se <- lomics_load_layer(cfg, layer, "processed")
  x <- assay(se, "values")
  v <- apply(x, 1, stats::var, na.rm = TRUE)
  x <- x[order(v, decreasing = TRUE)[seq_len(min(top_var, nrow(x)))], , drop = FALSE]
  long[[layer]] <- data.frame(
    sample  = rep(colnames(x), each = nrow(x)),
    feature = paste(layer, rownames(x), sep = "__"),   # unique across views
    view    = layer,
    value   = as.vector(x))
  m <- as.data.frame(colData(se))
  meta_all <- rbind(meta_all, m[, intersect(c("sample_id", "subject_id", "visit", "group",
                                              cfg$covariates), names(m))])
}
long <- do.call(rbind, long)
meta_all <- unique(meta_all)
meta_all <- meta_all[!duplicated(meta_all$sample_id), ]

mofa <- create_mofa(long)
d_opts <- get_default_data_options(mofa); d_opts$scale_views <- TRUE
m_opts <- get_default_model_options(mofa); m_opts$num_factors <- n_fac
t_opts <- get_default_training_options(mofa)
t_opts$seed <- 42; t_opts$convergence_mode <- "medium"; t_opts$drop_factor_threshold <- 0.01
mofa <- prepare_mofa(mofa, data_options = d_opts, model_options = m_opts,
                     training_options = t_opts)
mofa <- run_mofa(mofa, outfile = file.path(outdir, "mofa_model.hdf5"), use_basilisk = TRUE)

Z <- get_factors(mofa, factors = "all")[[1]]
fac <- data.frame(sample_id = rownames(Z), Z, check.names = FALSE)
fac <- merge(meta_all, fac, by = "sample_id")
fac$visit <- factor(fac$visit, levels = cfg$visits)
utils::write.csv(fac, file.path(outdir, "factors.csv"), row.names = FALSE)

W <- get_weights(mofa, views = "all", factors = "all")
for (vw in names(W)) {
  w <- data.frame(feature = sub("^.*?__", "", rownames(W[[vw]])), W[[vw]], check.names = FALSE)
  utils::write.csv(w, file.path(outdir, paste0("weights_", vw, ".csv")), row.names = FALSE)
}
r2 <- get_variance_explained(mofa)$r2_per_factor[[1]]
utils::write.csv(data.frame(factor = rownames(r2), r2, check.names = FALSE),
                 file.path(outdir, "variance_explained.csv"), row.names = FALSE)

## Repeated-measures tests per factor (ML fits + likelihood-ratio tests)
covs <- intersect(cfg$covariates, names(fac))
covs <- covs[vapply(covs, function(cv) length(unique(stats::na.omit(fac[[cv]]))) > 1, logical(1))]
cv_term <- if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""
multi_group <- length(unique(fac$group)) > 1
tests <- do.call(rbind, lapply(colnames(Z), function(f) {
  d <- fac[stats::complete.cases(fac[, c(f, "visit", "group", covs)]), ]
  d$y <- d[[f]]
  fit <- function(rhs) lme4::lmer(stats::as.formula(paste("y ~", rhs, cv_term, "+ (1 | subject_id)")),
                                  data = d, REML = FALSE)
  if (multi_group) {
    full <- fit("visit * group"); add <- fit("visit + group")
    no_t <- fit("group"); no_g <- fit("visit")
    data.frame(factor = f,
               p_visit = stats::anova(no_t, add)$`Pr(>Chisq)`[2],
               p_group = stats::anova(no_g, add)$`Pr(>Chisq)`[2],
               p_visit_x_group = stats::anova(add, full)$`Pr(>Chisq)`[2])
  } else {
    full <- fit("visit"); null <- fit("1")
    data.frame(factor = f, p_visit = stats::anova(null, full)$`Pr(>Chisq)`[2],
               p_group = NA, p_visit_x_group = NA)
  }
}))
for (p in grep("^p_", names(tests), value = TRUE)) tests[[sub("^p_", "fdr_", p)]] <- stats::p.adjust(tests[[p]], "BH")
utils::write.csv(tests, file.path(outdir, "factor_tests.csv"), row.names = FALSE)
print(tests, row.names = FALSE)

## Factor trajectories (mean +/- SE per group and visit)
fl <- stats::reshape(fac[, c("sample_id", "subject_id", "visit", "group", colnames(Z))],
                     direction = "long", varying = colnames(Z), v.names = "score",
                     timevar = "factor", times = colnames(Z), idvar = "sample_id")
agg <- stats::aggregate(score ~ factor + group + visit, fl, function(s) c(m = mean(s), se = stats::sd(s) / sqrt(length(s))))
agg <- do.call(data.frame, agg)
p <- ggplot(agg, aes(visit, score.m, colour = group, group = group)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = score.m - score.se, ymax = score.m + score.se), width = .15) +
  facet_wrap(~ factor, scales = "free_y") +
  labs(y = "Factor score (mean ± SE)", title = "MOFA factor trajectories") + theme_minimal()
ggsave(file.path(outdir, "factor_trajectories.pdf"), p, width = 12, height = 8)
lomics_session_log(cfg, "mofa")
