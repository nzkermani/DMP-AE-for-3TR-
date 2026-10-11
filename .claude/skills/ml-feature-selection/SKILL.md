---
name: ml-feature-selection
description: Feature selection and biomarker ranking for a binary, multi-class, continuous or survival outcome — variance filtering, univariate limma / Cox ranking, elastic-net stability selection, random-forest permutation importance and Boruta, combined into a consensus ranking. Use when the user asks for feature selection, biomarker discovery, the most important or predictive features, or a reduced feature panel.
---

# Feature selection

Category: **Machine learning**.

```bash
Rscript .claude/skills/ml-feature-selection/scripts/feature_selection.R \
  --outcome=binary|continuous|survival|<column> --layers=all --visit=V1 \
  [--change_from=V1] [--methods=variance,univariate,stability_enet,rf,boruta] \
  [--alpha=0.5] [--n_boot=100] [--pi_thr=0.6] [--top_n=50] [--seed=42]
```
Outputs `results/feature_selection/<outcome>_<tag>/`:
`ranking_<method>.csv`, `consensus_ranking.csv`, `selected_features.txt`,
`plots.pdf`.

## Methods by family
| family | method | uses outcome | notes |
|---|---|---|---|
| filter, unsupervised | variance (MAD) | no | safe pre-step before any supervised analysis |
| filter, univariate | limma moderated t/F; Cox for survival | yes | ignores correlation between features |
| embedded | elastic-net stability selection | yes | selection frequency over half-subsamples; controls false selections better than one lasso fit |
| wrapper / ensemble | random forest permutation importance | yes | captures non-linearity and interactions; random survival forest for survival |
| all-relevant | Boruta | yes | optional; skipped if not installed |

`selected_features.txt` holds the stability-selected set (frequency ≥
`pi_thr`), or the top-`top_n` consensus features if none pass.

## Rules
- **Selection bias**: accuracy measured on the same subjects used for
  selection is inflated. For predictive performance, run
  `ml-classification`, which repeats selection inside every CV fold.
  Passing `selected_features.txt` to it via `--features_file` re-introduces
  the bias. Only do that for a list derived from an independent cohort.
- Report method agreement (Spearman between rankings). Features that only
  one method finds are weak.
- Correlated features (e.g. gene and its protein) split importance and
  selection frequency. A low frequency does not mean irrelevant.
- Use one visit (one row per subject). Change scores (`--change_from`)
  select features whose change relates to the outcome.
