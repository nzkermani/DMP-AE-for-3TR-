# Multi-omics analysis project (4 visits) — instructions for Claude

## Study
Design, visits, arms, covariates, outcomes and omic layers are defined in
`study_config.json`. Read it at the start of every session (or call MCP
`project_status`). Never assume visit order, visit timing, group labels,
outcome coding or the scale of a layer. Take them from the config or ask.

## Skills (.claude/skills/), by category
Start with `omics-analysis-guide`: it maps questions to skills.

| category | skills |
|---|---|
| 1 Data preparation | `prep-study-setup` → `prep-qc-normalisation` (always first) |
| 2 Exploration | `explore-dimension-reduction` (PCA, t-SNE, UMAP) |
| 3 Statistics | `stats-group-testing`, `stats-regression`, `stats-longitudinal-de`, `stats-survival` |
| 4 Clustering | `cluster-samples` (subjects), `cluster-trajectories-longomicstraj` (features by trajectory) |
| 5 Machine learning | `ml-feature-selection`, `ml-classification` |
| 6 Interpretation & integration | `interpret-pathway-enrichment`, `integrate-multiomics` |

Run `stats-longitudinal-de` before `cluster-trajectories-longomicstraj`
(it uses the time-effect FDR for gating).

## How to run things
- Analyses: run the skill scripts with `Rscript` from the project root
  (Bash), or with MCP `run_skill_script`. Every run writes to
  `results/<category>/<run>/` plus session info in `results/logs/`.
- Inspection: MCP `run_r`, `read_result`, `project_status`,
  `list_skills`, `query_trajectories`, `feature_report`.
- Don't run long model fits through `run_r` in the user's interactive
  session; it blocks RStudio.
- New or one-off code goes in `analysis/` as numbered R scripts that
  source `R/lomics_utils.R` (use `lomics_xy()` to build analysis
  matrices). Don't edit skill scripts for one-off variants.

## Scientific rules
- Report numbers only from result files read in this session; never from memory.
- Never invent function names or arguments. Check with
  `run_r("args(pkg::fn)")` if unsure.
- Ask before dropping samples, subjects or features beyond the documented filters.
- Keep selection and evaluation separate: no accuracy, AUC, C-index or p-value
  reported on the same data used to select features, clusters or models,
  unless it comes from the nested-CV / out-of-fold outputs.
- Distinguish within-arm change from between-arm difference in change
  (interaction); trajectory shape (OTS) from statistical evidence (FDR);
  exploratory from confirmatory.
- Survival predictors must be measured at or before time zero (baseline by default).
- Flag confounding (batch aligned with visit or arm) explicitly.
- State seeds and any random imputation (MinProb).
