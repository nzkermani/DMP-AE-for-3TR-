# DMP-AE for 3TR — Data Management Platform Analytical Engine

The **Analytical Engine (AE)** of the 3TR **Data Management Platform (DMP)**:
Claude Code skills, an R MCP server and tested R scripts for analysing
transcriptomic, proteomic and metabolomic data collected at four visits in
RStudio. Methods are grouped into six categories; every script reads one
study config and writes to `results/<category>/`.

## Method catalogue

| # | Category | Skill | Methods |
|---|---|---|---|
| 1 | **Data preparation** | `prep-study-setup` | align layers to metadata, design audit, SummarizedExperiment / MultiAssayExperiment |
| | | `prep-qc-normalisation` | visit-aware missingness filter, half-min / MinProb imputation, TMM, PQN, median, QC-RSD, QC PCA |
| 2 | **Exploration & dimensionality reduction** | `explore-dimension-reduction` | PCA (scree, loadings, PC–metadata tests), t-SNE, UMAP; single / multi-layer (block-scaled) |
| 3 | **Statistics** | `stats-group-testing` | Welch t, Wilcoxon, ANOVA, Kruskal-Wallis, paired t, signed-rank, Friedman, mixed-model LRT, chi-square / Fisher Table 1 |
| | | `stats-regression` | limma omics-wide regression, per-feature linear / logistic (OR per SD), multivariable lm / glm with CI, VIF, diagnostics |
| | | `stats-longitudinal-de` | repeated-measures limma (duplicateCorrelation, voom) / dream: change from baseline, steps, arm × visit, F-tests |
| | | `stats-survival` | Kaplan-Meier + log-rank, per-feature Cox (HR, PH test), elastic-net Cox (CV C-index, risk groups), multivariable Cox |
| 4 | **Clustering** | `cluster-samples` | k-means, PAM, Ward hierarchical, Gaussian mixture, consensus; silhouette / BIC / PAC; bootstrap Jaccard; ARI |
| | | `cluster-trajectories-longomicstraj` | [LongOmicsTraj](https://github.com/nzkermani/LongOmicsTraj) trajectory signatures (empirical / LMM / GAM / maSigPro), per-arm consensus, significance gating |
| 5 | **Machine learning** | `ml-feature-selection` | variance filter, limma / Cox ranking, elastic-net stability selection, RF permutation importance, Boruta, consensus ranking |
| | | `ml-classification` | PLS-DA (VIP), sPLS-DA, PCA-LDA, elastic net, random forest, linear SVM; repeated nested CV, grouped folds, permutation null |
| 6 | **Interpretation & integration** | `interpret-pathway-enrichment` | ssGSEA / GSVA pathway layers, fgsea on contrasts |
| | | `integrate-multiomics` | MOFA2 + mixed-model factor tests; DIABLO with multilevel decomposition |
| — | Router | `omics-analysis-guide` | question → skill map, shared conventions |

## Repository layout
```
CLAUDE.md                  rules Claude Code follows in this project
.mcp.json                  registers the r-longomics MCP server (mcptools)
study_config.json          visits, times, arms, covariates, outcomes, layers
R/lomics_utils.R           shared helpers (config, loading, lomics_xy(), args)
mcp/longomics_tools.R      MCP tools: run_r, project_status, list_skills,
                           run_skill_script, query_trajectories,
                           feature_report, read_result
.claude/skills/<skill>/    SKILL.md + scripts/ (+ references/)
setup/                     install_packages.R, Rprofile_snippet.R
testdata/                  synthetic 60-subject, 4-visit, 3-omic study + generator
tests/run_all_tests.sh     runs every script and checks planted signals
```

## Setup
1. Clone into (or open as) an RStudio project.
2. In RStudio: `source("setup/install_packages.R")`. LongOmicsTraj needs R ≥ 4.6.
3. Optional (lets Claude use your live R session): add
   `setup/Rprofile_snippet.R` to the project `.Rprofile` and restart R.
4. Terminal at the project root: `claude`. Approve the `r-longomics` MCP
   server; `/mcp` shows its status.
5. Check everything works: `bash tests/run_all_tests.sh` (about 4 min).
6. Point `study_config.json` at your data. Patient-level data go in
   `data/`, which is git-ignored.

## Shared conventions
- `--visit=V1` gives one row per subject. `--visit=V3 --change_from=V1`
  gives change scores. `--visit=all` gives every sample.
- `--layers=proteomics,metabolomics` joins layers by subject. Features
  are named `<layer>__<feature>`.
- `--outcome=binary|continuous|survival` refers to `outcomes` in the
  config; any metadata column name also works.

## Test status
`tests/run_all_tests.sh` on the synthetic data (R 4.3.3; LongOmicsTraj
0.99.0 installed with its R-version requirement relaxed): **24 passed,
0 failed, 3 skipped**. Planted signals were recovered:

- the baseline endotype by PCA (PC1) and by k-means, with perfect
  agreement and bootstrap Jaccard 1.0
- 30/30 top consensus features by feature selection
- AUC 1.0 by classification
- treatment trajectories by LongOmicsTraj with gating

On six random-label outcomes, mean cross-validated AUC per model was
0.48–0.53, so the classification pipeline does not leak.

**Not executed** (packages unavailable in the build environment): dream,
fgsea, msigdbr, MOFA2, DIABLO / sPLS-DA (mixOmics), maSigPro, Boruta, and
the MCP server through mcptools/ellmer. The MCP tool functions were tested
by calling them directly in R. Expect small API fixes on first run of the
untested parts.

## Notes on LongOmicsTraj 0.99.0 (found while building)
- `lot_from_*()` with `overwrite = TRUE` replaces the VMM for all assays,
  so the scripts use one object per layer.
- `lot_consensus()` votes across estimator × group, merging arms. The
  scripts call it within each arm.
- OTS tables report `object_type = "gene"` for every layer.
- With the default thresholds, null features mostly get non-flat labels,
  hence the FDR gating against `stats-longitudinal-de`.
