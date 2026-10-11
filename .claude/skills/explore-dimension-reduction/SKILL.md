---
name: explore-dimension-reduction
description: Dimensionality reduction and visual exploration of omics data — PCA (scree, loadings, PC-metadata association tests), t-SNE and UMAP for one or several layers at one visit, change scores or all visits with subject paths. Use when the user asks for PCA, t-SNE, UMAP, embeddings, sample structure, outliers or what drives the main variation.
---

# Dimensionality reduction (PCA, t-SNE, UMAP)

Category: **Exploration**. Needs `prep-qc-normalisation` outputs.

```bash
Rscript .claude/skills/explore-dimension-reduction/scripts/dimension_reduction.R \
  --layers=proteomics,metabolomics --visit=V1 [--change_from=V1] \
  [--methods=pca,tsne,umap] [--colour_by=group,endotype] [--top_var=2000] \
  [--perplexity=auto] [--n_neighbors=15] [--seed=42]
```
Outputs `results/dimred/<layers>_<visit>/`: `pca_scores.csv`,
`pca_variance.csv`, `pca_loadings.csv`, `pc_associations.csv`,
`tsne_coords.csv`, `umap_coords.csv`, `plots.pdf`.

## Choices
| Situation | Settings |
|---|---|
| Baseline structure / endotypes | `--visit=V1` |
| What changes over time | `--visit=V3 --change_from=V1` (subject-level change) |
| All samples, trajectories | `--visit=all` (draws each subject's path) |
| Joint view of layers | `--layers=a,b` — autoscaled then block-scaled by 1/sqrt(p) so a large layer cannot dominate |

## Interpretation rules
- **PCA is the analysis; t-SNE/UMAP are pictures.** Distances between
  t-SNE/UMAP clusters and cluster sizes are not meaningful; never report
  "separation" from them without a PCA or statistical test behind it.
  Both run on the top 30 PCs; perplexity auto = min(30, (n−1)/3 − 1).
  Results depend on the seed — re-run with another seed before
  interpreting a pattern.
- `pc_associations.csv`: p-values for PC1–10 vs each metadata variable
  (one-way lm/ANOVA). Use it to say what the main axes reflect (batch,
  sex, visit, endotype). These are exploratory and unadjusted. With
  `visit=all`, samples from the same subject are not independent, so the
  p-values are anticonservative.
- Batch on PC1/PC2 → add it to covariates in `study_config.json` for every
  model, and report it.
- Outliers: list them by ID and ask the user; don't drop them silently.
