---
name: ml-multilevel-splsda
description: Paired / repeated-measures biomarker panel discovery between two visits with multilevel sparse PLS-DA (within-subject decomposition), keepX tuning by leave-one-subject-out CV, stepwise backward elimination to a small panel (AUC-drop rule), training AUROC with DeLong CI, LOSO-CV of the fixed panel, nested LOSO and within-subject permutation null, panel boxplots. Use when the user asks which features change between baseline and follow-up within patients, for a paired sPLS-DA / multilevel PLS-DA, a breath VOC or metabolite panel that separates pre- vs post-treatment, or to reproduce such an analysis from a paper.
---

# Multilevel sPLS-DA (paired visits)

Category: **Machine learning**. Needs a processed layer (`prep-*`).

```bash
S=.claude/skills/ml-multilevel-splsda/scripts/multilevel_splsda.R
Rscript $S --config=study_config.json --layer=breath_voc --visits=BL,M4 \
  [--ncomp=2] [--grid=5,10,15,20,30,50,75,100] [--zscore=TRUE] \
  [--drop_max=0.05] [--drop_type=relative|absolute] [--target_auc=0.9|NA] \
  [--elim_auc=train|loso] [--nested=TRUE] [--nested_grid=5,10,20,50] \
  [--nested_folds=5|loso] [--nested_repeats=5] \
  [--n_perm=20] [--plot_visits=BL,M1,M4] [--colour_by=group] \
  [--engine=auto|mixOmics|internal] [--panel=f1,f2,...] [--tag=name]
```
`--panel` evaluates a given feature set with no selection (e.g. a published panel).
Runtime with 3,000 features and 25 pairs: about 3 minutes, plus about 15 s per permutation.

Outputs `<output_dir>/splsda/<v1>vs<v2>_<layer>[_<tag>]/`: `tuning.csv`,
`elimination_log.csv`, `panel.csv` (compound, m/z, RT, loading, direction,
n up/down), `performance.csv`, `predictions_*.csv`, `nested_panels.csv`,
`nested_repeats.csv`, `permutation.csv`, `sample_plot`, `roc`, `panel_boxplots` (.pdf/.png), `summary.txt`.

## Method
1. Subjects with both visits only; optional z-score per feature (paper: caret center/scale).
2. Multilevel decomposition: `mixOmics::withinVariation` (or subject-mean centring).
3. keepX tuned per component by leave-one-subject-out (LOSO) AUC.
4. Backward elimination from the tuned panel: least important feature first. A
   feature is removed if the AUC drops by ≤ `drop_max` (relative by default)
   and stays ≥ `target_auc`. The AUC is the training AUC (`train`, as published) or LOSO (`loso`).
5. Engine: mixOmics if installed, otherwise the internal sPLS implementation in
   `R/lomics_spls.R`. It is checked against `pls::plsr` in `tests/test_spls_internal.R`. Report the engine.

## Four performance estimates — report all of them
| | what | honest? |
|---|---|---|
| A training | panel chosen and evaluated on all subjects | no — optimistic |
| B LOSO fixed panel | refit without each subject, panel fixed beforehand (paper's LOSO-CV) | no — selection leaked |
| C nested CV | tuning + elimination redone inside every outer fold (subject-grouped 5-fold × 5 repeats by default) | yes |
| D permutation | A and B after flipping visit labels within random subjects | shows what A/B give with **no** effect |

- **Paired accuracy** is the share of subjects whose later sample scores higher. Prefer it to pooled AUC for C, because fold models have different score scales.
- Avoid `--nested_folds=loso` when selection is aggressive. It is biased **downwards**: when a non-responder is held out, its feature looks stronger in the remaining subjects and gets selected, so that subject is then misclassified. Use repeated grouped K-fold (the default) and report the range over repeats.
- Multilevel AUC compares samples across subjects after within-subject centring. It runs high even for modest paired effects. With thousands of features and ~25 pairs, A and B are often above 0.9 under the null. Always show D next to them.
- A panel from A/B is hypothesis-generating. State the nested estimate (mean and range over repeats) and how often each panel feature was re-selected across fold fits (`summary.txt`).
- Descriptive per-feature Wilcoxon p-values in `panel.csv` are post-selection. Never present them as tests.
- Direction comes from mean within-subject change in z-scores (`mean_change_z`).
