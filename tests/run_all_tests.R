## Run every skill script on the synthetic test data — works on Windows,
## macOS and Linux. In RStudio, with the project open:
##     source("tests/run_all_tests.R")
## or from a terminal at the project root:  Rscript tests/run_all_tests.R
## Writes results/ and tests/test_log.txt. Steps needing packages that are not
## installed are reported as SKIP, not FAIL.

if (!file.exists("study_config.json")) stop("Run this from the project root (open the RStudio project first).")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
log_file <- "tests/test_log.txt"; cat("", file = log_file)
K <- ".claude/skills"
res <- data.frame(step = character(), status = character(), seconds = numeric())

run <- function(name, pkg, script, args = character()) {
  if (!is.na(pkg) && !requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("SKIP  %-32s (needs %s)", name, pkg))
    res[nrow(res) + 1, ] <<- list(name, paste("SKIP: needs", pkg), NA); return(invisible())
  }
  t0 <- Sys.time()
  out <- suppressWarnings(system2(rscript, shQuote(c(file.path(K, script), args)), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); ok <- is.null(st) || st == 0
  cat(sprintf("\n===== %s =====\n", name), out, sep = "\n", file = log_file, append = TRUE)
  secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  message(sprintf("%s  %-32s %5.1fs", if (ok) "OK  " else "FAIL", name, secs))
  if (!ok) message("      last lines: ", paste(utils::tail(out, 3), collapse = " | "))
  res[nrow(res) + 1, ] <<- list(name, if (ok) "OK" else "FAIL", secs)
}

unlink("results", recursive = TRUE)
run("prep: setup", NA, "prep-study-setup/scripts/build_layers.R", "study_config.json")
run("prep: QC + normalisation", NA, "prep-qc-normalisation/scripts/qc_normalise.R", "study_config.json")
run("explore: PCA/tSNE/UMAP V1", NA, "explore-dimension-reduction/scripts/dimension_reduction.R", c("--layers=all", "--visit=V1"))
run("explore: all visits", NA, "explore-dimension-reduction/scripts/dimension_reduction.R", c("--layers=proteomics", "--visit=all", "--methods=pca,umap"))
run("stats: between groups", NA, "stats-group-testing/scripts/group_tests.R", c("--design=between", "--by=endotype"))
run("stats: change between arms", NA, "stats-group-testing/scripts/group_tests.R", c("--design=between", "--by=group", "--visit=V3", "--change_from=V1"))
run("stats: paired visits", NA, "stats-group-testing/scripts/group_tests.R", c("--design=paired", "--visits=V1,V3", "--subset_col=group", "--subset=Treatment"))
run("stats: repeated (Friedman)", NA, "stats-group-testing/scripts/group_tests.R", c("--design=repeated", "--layers=proteomics"))
run("stats: clinical table", NA, "stats-group-testing/scripts/group_tests.R", c("--design=clinical", "--by=group"))
run("regression: omics-wide", NA, "stats-regression/scripts/regression.R", c("--mode=features", "--outcome=continuous"))
run("regression: logistic/feature", NA, "stats-regression/scripts/regression.R", c("--mode=outcome", "--outcome=binary", "--layers=proteomics"))
run("regression: formula", NA, "stats-regression/scripts/regression.R", c("--mode=formula", "--layers=proteomics", "--formula=fev1 ~ age + sex + PROT21"))
run("longitudinal: limma", NA, "stats-longitudinal-de/scripts/longitudinal_de.R", c("study_config.json", "all", "limma"))
run("longitudinal: dream", "variancePartition", "stats-longitudinal-de/scripts/longitudinal_de.R", c("study_config.json", "proteomics", "dream"))
run("survival: KM", NA, "stats-survival/scripts/survival.R", c("--mode=km", "--strata=group"))
run("survival: Cox per feature", NA, "stats-survival/scripts/survival.R", c("--mode=cox_features", "--layers=proteomics"))
run("survival: penalised Cox", "glmnet", "stats-survival/scripts/survival.R", c("--mode=cox_penalised", "--layers=proteomics", "--repeats=3"))
run("survival: Cox formula", NA, "stats-survival/scripts/survival.R", c("--mode=cox_formula", "--layers=proteomics", "--formula=age + PROT21"))
run("cluster: samples", NA, "cluster-samples/scripts/cluster_samples.R", c("--layers=transcriptomics,metabolomics", "--boot=20"))
run("cluster: trajectories (LOT)", "LongOmicsTraj", "cluster-trajectories-longomicstraj/scripts/run_longomicstraj.R", c("study_config.json", "all", "empirical,lmm"))
run("ml: feature selection", "glmnet", "ml-feature-selection/scripts/feature_selection.R", c("--outcome=endotype", "--n_boot=30"))
run("ml: classification", "glmnet", "ml-classification/scripts/classification.R", c("--outcome=endotype", "--repeats=2"))
run("ml: classification null", "glmnet", "ml-classification/scripts/classification.R", c("--outcome=endotype", "--repeats=1", "--permute_seed=1", "--layers=proteomics"))
run("pathways: ssGSEA layer", "GSVA", "interpret-pathway-enrichment/scripts/pathway_scores.R", c("study_config.json", "transcriptomics", "testdata/test_sets.gmt"))
run("integration: MOFA2", "MOFA2", "integrate-multiomics/scripts/integrate_mofa.R", c("study_config.json", "5", "300"))
run("integration: DIABLO", "mixOmics", "integrate-multiomics/scripts/integrate_diablo.R", c("study_config.json", "Treatment", "2", "10", "300"))

## Planted-signal checks (see testdata/simulate_testdata.R)
chk <- tryCatch({
  cl <- utils::read.csv("results/clustering/transcriptomics+metabolomics_V1/clinical_association.csv")
  stopifnot(cl$p[cl$variable == "endotype"] < 1e-6)
  fs <- utils::read.csv(Sys.glob("results/feature_selection/endotype_*/consensus_ranking.csv")[1])
  top <- sub(".*__", "", utils::head(fs$feature, 20))
  stopifnot(mean(top %in% c(paste0("GENE", 100:139), paste0("MET", 30:39))) >= 0.9)
  cs <- utils::read.csv("results/classification/endotype_transcriptomics+proteomics+metabolomics_V1/cv_metrics.csv")
  stopifnot(mean(cs$AUC) > 0.9)
  if (file.exists("results/trajectories/ots_all_layers.csv")) {
    o <- utils::read.csv("results/trajectories/ots_all_layers.csv")
    trt <- o[o$layer == "proteomics" & o$group == "Treatment" & o$feature %in% paste0("PROT", 1:15), ]
    stopifnot(mean(trt$topology_gated != "no_time_effect", na.rm = TRUE) > 0.5)
  }
  TRUE
}, error = function(e) { message("check error: ", conditionMessage(e)); FALSE })
message(sprintf("%s  %-32s", if (chk) "OK  " else "FAIL", "check: planted signals"))
res[nrow(res) + 1, ] <- list("check: planted signals", if (chk) "OK" else "FAIL", NA)

n_ok <- sum(res$status == "OK"); n_fail <- sum(res$status == "FAIL"); n_skip <- sum(startsWith(res$status, "SKIP"))
message(sprintf("\npassed %d, failed %d, skipped %d  (details: %s)", n_ok, n_fail, n_skip, log_file))
utils::write.csv(res, "tests/test_summary.csv", row.names = FALSE)
invisible(res)
