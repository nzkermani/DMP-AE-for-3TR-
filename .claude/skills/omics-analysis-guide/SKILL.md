---
name: omics-analysis-guide
description: Catalogue and decision guide for every analysis in this project (data preparation, missing data and imputation, exploration and dimensionality reduction, statistical testing, regression, longitudinal and mixed-effects models, survival analysis, clustering, feature selection, classification, pathway and multi-omics integration). Use at the start of any analysis request, when the user asks what methods are available or which analysis to run, or when a question could map to several methods.
---

# Analysis catalogue

All skills read `study_config.json`, work on the processed layers in
`results/data/`, and are run from the project root with `Rscript` (or the
MCP tool `run_skill_script`). Options for the newer scripts are
`--key=value`.

## 1. Data preparation — run first, in order
| skill | does |
|---|---|
| `prep-study-setup` | align matrices to metadata, audit visits/arms, build SummarizedExperiments |
| `prep-qc-normalisation` | per-omic filtering, imputation, normalisation (TMM / PQN / median), QC PCA |
| `prep-breathomics-gcms` | breath GC-MS / VOC layers instead of the above: contaminant removal, log2, limma batch correction (visit preserved), PCA before/after |

## 2. Missing data — diagnose before choosing an imputation
| skill | methods |
|---|---|
| `prep-missing-data` | diagnose (MNAR signature, by visit / arm / batch, dropout patterns) → benchmark imputations (half-min, MinProb, kNN, BPCA, mean; MCAR + MNAR masking) → imputed layer; mice multiple imputation + Rubin pooling for clinical covariates |

## 3. Exploration & dimensionality reduction
| skill | methods |
|---|---|
| `explore-dimension-reduction` | PCA (scree, loadings, PC–metadata tests), t-SNE, UMAP; single or multi-layer; visit, change or all visits |

## 4. Statistics (cross-sectional, paired, survival)
| skill | methods |
|---|---|
| `stats-group-testing` | t / Wilcoxon, ANOVA / Kruskal, paired t / signed-rank, Friedman, Table 1 (chi-square / Fisher) |
| `stats-regression` | linear and logistic regression: limma omics-wide association, per-feature glm (OR / beta per SD), multivariable lm/glm with VIF and diagnostics |
| `stats-survival` | Kaplan-Meier + log-rank, per-feature Cox (HR, PH test), elastic-net Cox (CV C-index), multivariable Cox |
| `stats-clinical-association` | Table 1 (normality-based), clinical change between visits (paired t / signed-rank), Kendall feature–clinical correlation, feature levels by medication |

## 5. Longitudinal & mixed-effects models
| skill | methods |
|---|---|
| `stats-mixed-models` | lme4 mixed models per feature or clinical variable: visit / time × arm, random intercept or slope, LRTs + FDR, ICC, emmeans change vs baseline; keeps subjects with missing visits |
| `stats-longitudinal-de` | repeated-measures limma / dream for thousands of features: change from baseline, visit steps, arm × visit interaction, overall time F-test |

Pick `stats-mixed-models` for clinical variables, panels, unbalanced or
missing visits and per-feature random effects; `stats-longitudinal-de` for
fast omics-wide screens (limma with a shared within-subject correlation).

## 6. Clustering
| skill | clusters |
|---|---|
| `cluster-samples` | **subjects** (endotypes): k-means, PAM, Ward, GMM, consensus; silhouette / BIC / PAC; bootstrap stability; ARI |
| `cluster-trajectories-longomicstraj` | **features** by longitudinal shape: LongOmicsTraj Omics Trajectory Signatures, k-means / FlexMix trajectory clusters |

## 7. Machine learning
| skill | methods |
|---|---|
| `ml-feature-selection` | variance filter, limma / Cox ranking, elastic-net stability selection, RF importance, Boruta → consensus |
| `ml-classification` | PLS-DA (VIP), sPLS-DA, PCA-LDA, elastic net, random forest, SVM; nested repeated CV; permutation null |
| `ml-multilevel-splsda` | paired (multilevel) sPLS-DA between two visits: keepX tuning, backward elimination to a panel, training / LOSO / nested CV / permutation |

## 8. Interpretation & integration
| skill | methods |
|---|---|
| `interpret-pathway-enrichment` | ssGSEA / GSVA pathway layers (then run any skill on them), fgsea on contrasts |
| `integrate-multiomics` | MOFA2 factors + mixed models; DIABLO multilevel |

## Question → skill
| the user asks… | use |
|---|---|
| "what does the data look like / any batch effect / outliers?" | explore-dimension-reduction, then prep-qc-normalisation diagnostics |
| "is X different between groups at baseline?" | stats-group-testing `between` (or stats-regression `features` for omics with covariates) |
| "what changes over time / with treatment?" | stats-longitudinal-de (interaction = treatment effect) |
| "which features are associated with FEV1 / response?" | stats-regression |
| "what predicts time to exacerbation?" | stats-survival |
| "are there patient subgroups / endotypes?" | cluster-samples → stats-survival km / stats-group-testing on clusters |
| "which features rise then fall / rebound?" | cluster-trajectories-longomicstraj |
| "can omics predict response, how well?" | ml-classification (+ permutation null) |
| "which features are the best biomarkers?" | ml-feature-selection, then ml-classification for honest accuracy |
| "which pathways?" | interpret-pathway-enrichment |
| "which VOCs / metabolites separate baseline from follow-up within patients?" | ml-multilevel-splsda (paired) |
| "mixed model / random effects / change over time with missing visits / treatment × time" | stats-mixed-models |
| "how much is missing / below detection / which imputation / dropout / MCAR-MAR-MNAR / multiple imputation?" | prep-missing-data |
| "baseline table / did clinical scores change / do the markers correlate with ACQ, eosinophils?" | stats-clinical-association |
| "breath (GC-MS) data: contaminants, batch effect?" | prep-breathomics-gcms |
| "reproduce paper X" / "show me a complete workflow" | see `usecases/` (worked notebooks on synthetic data) |
| "what is shared across omics?" | integrate-multiomics |

## Shared conventions
- `--visit=V1` gives one row per subject. `--visit=V3 --change_from=V1`
  gives subject change scores. `--visit=all` gives every sample; the
  scripts that support it account for repeated subjects where stated.
- `--layers=a,b` joins layers by subject. Feature names become
  `<layer>__<feature>`.
- Outcomes: `--outcome=binary|continuous|survival` refer to
  `study_config.json` → `outcomes`; any metadata column name also works.
- If no applicable skill exists, write a new numbered script in
  `analysis/` that sources `R/lomics_utils.R` and uses `lomics_xy()`.
  Don't improvise inside the skill scripts.

## Before answering with results
Read the output files of the run you just did. State n, the design
(visit / change / all), covariates, multiple-testing method and
validation scheme. Flag anything exploratory as exploratory.
