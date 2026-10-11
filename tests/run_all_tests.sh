#!/usr/bin/env bash
# Run every skill script on the synthetic test data, from the project root:
#   bash tests/run_all_tests.sh            (uses study_config.json -> testdata/)
# Writes results/ and tests/test_log.txt; exits non-zero if any step fails.
# Steps needing packages that are not installed (MOFA2, mixOmics, fgsea) are
# reported as SKIP, not FAIL.
set -u
LOG=tests/test_log.txt; : > "$LOG"
pass=0; fail=0; skip=0
run() {  # name, needs-package (or -), command...
  local name=$1 pkg=$2; shift 2
  if [ "$pkg" != "-" ] && ! Rscript -e "quit(status = !requireNamespace('$pkg', quietly = TRUE))" >/dev/null 2>&1; then
    echo "SKIP  $name (needs $pkg)" | tee -a "$LOG"; skip=$((skip+1)); return
  fi
  if "$@" >>"$LOG" 2>&1; then echo "OK    $name" | tee -a "$LOG"; pass=$((pass+1))
  else echo "FAIL  $name" | tee -a "$LOG"; fail=$((fail+1)); fi
}
K=.claude/skills
rm -rf results
run "prep: setup"                 -  Rscript $K/prep-study-setup/scripts/build_layers.R study_config.json
run "prep: QC + normalisation"    -  Rscript $K/prep-qc-normalisation/scripts/qc_normalise.R study_config.json
run "explore: PCA/tSNE/UMAP V1"   -  Rscript $K/explore-dimension-reduction/scripts/dimension_reduction.R --layers=all --visit=V1
run "explore: all visits"         -  Rscript $K/explore-dimension-reduction/scripts/dimension_reduction.R --layers=proteomics --visit=all --methods=pca,umap
run "stats: between groups"       -  Rscript $K/stats-group-testing/scripts/group_tests.R --design=between --by=endotype
run "stats: change between arms"  -  Rscript $K/stats-group-testing/scripts/group_tests.R --design=between --by=group --visit=V3 --change_from=V1
run "stats: paired visits"        -  Rscript $K/stats-group-testing/scripts/group_tests.R --design=paired --visits=V1,V3 --subset_col=group --subset=Treatment
run "stats: repeated (Friedman)"  -  Rscript $K/stats-group-testing/scripts/group_tests.R --design=repeated --layers=proteomics
run "stats: clinical table"       -  Rscript $K/stats-group-testing/scripts/group_tests.R --design=clinical --by=group
run "regression: omics-wide"      -  Rscript $K/stats-regression/scripts/regression.R --mode=features --outcome=continuous
run "regression: logistic/feature" - Rscript $K/stats-regression/scripts/regression.R --mode=outcome --outcome=binary --layers=proteomics
run "regression: formula"         -  Rscript $K/stats-regression/scripts/regression.R --mode=formula --layers=proteomics "--formula=fev1 ~ age + sex + PROT21"
run "longitudinal: limma"         -  Rscript $K/stats-longitudinal-de/scripts/longitudinal_de.R study_config.json all limma
run "longitudinal: dream"         variancePartition Rscript $K/stats-longitudinal-de/scripts/longitudinal_de.R study_config.json proteomics dream
run "survival: KM"                -  Rscript $K/stats-survival/scripts/survival.R --mode=km --strata=group
run "survival: Cox per feature"   -  Rscript $K/stats-survival/scripts/survival.R --mode=cox_features --layers=proteomics
run "survival: penalised Cox"     glmnet Rscript $K/stats-survival/scripts/survival.R --mode=cox_penalised --layers=proteomics --repeats=3
run "survival: Cox formula"       -  Rscript $K/stats-survival/scripts/survival.R --mode=cox_formula --layers=proteomics "--formula=age + PROT21"
run "cluster: samples"            -  Rscript $K/cluster-samples/scripts/cluster_samples.R --layers=transcriptomics,metabolomics --boot=20
run "cluster: trajectories (LOT)" LongOmicsTraj Rscript $K/cluster-trajectories-longomicstraj/scripts/run_longomicstraj.R study_config.json all empirical,lmm
run "ml: feature selection"       glmnet Rscript $K/ml-feature-selection/scripts/feature_selection.R --outcome=endotype --n_boot=30
run "ml: classification"          glmnet Rscript $K/ml-classification/scripts/classification.R --outcome=endotype --repeats=2
run "ml: classification null"     glmnet Rscript $K/ml-classification/scripts/classification.R --outcome=endotype --repeats=1 --permute_seed=1 --layers=proteomics
run "pathways: ssGSEA layer"      GSVA Rscript $K/interpret-pathway-enrichment/scripts/pathway_scores.R study_config.json transcriptomics testdata/test_sets.gmt
run "integration: MOFA2"          MOFA2 Rscript $K/integrate-multiomics/scripts/integrate_mofa.R study_config.json 5 300
run "integration: DIABLO"         mixOmics Rscript $K/integrate-multiomics/scripts/integrate_diablo.R study_config.json Treatment 2 10 300

# Planted-signal checks (see testdata/simulate_testdata.R)
run "check: planted signals" - Rscript -e '
  cl <- read.csv("results/clustering/transcriptomics+metabolomics_V1/clinical_association.csv")
  stopifnot(cl$p[cl$variable == "endotype"] < 1e-6)
  fs <- read.csv(Sys.glob("results/feature_selection/endotype_*/consensus_ranking.csv"))
  top <- sub(".*__", "", head(fs$feature, 20))
  stopifnot(mean(top %in% c(paste0("GENE", 100:139), paste0("MET", 30:39))) >= 0.9)
  cs <- read.csv(Sys.glob("results/classification/endotype_transcriptomics+proteomics+metabolomics_V1/cv_metrics.csv"))
  stopifnot(mean(cs$AUC) > 0.9)
  if (file.exists("results/trajectories/ots_all_layers.csv")) {
    o <- read.csv("results/trajectories/ots_all_layers.csv")
    trt <- o[o$layer == "proteomics" & o$group == "Treatment" & o$feature %in% paste0("PROT", 1:15), ]
    stopifnot(mean(trt$topology_gated != "no_time_effect", na.rm = TRUE) > 0.5)
  }
  cat("planted signals recovered\n")'

echo "passed $pass, failed $fail, skipped $skip  (details: $LOG)"
[ "$fail" -eq 0 ]
