---
name: prep-qc-normalisation
description: Quality control, filtering, imputation and normalisation of transcriptomic (RNA-seq counts), proteomic (Olink/SomaScan/MS) and metabolomic layers in a repeated-visit study, producing log2-scale processed SummarizedExperiments plus PCA diagnostics. Use after prep-study-setup and before any differential, trajectory or integration analysis, or when the user asks about outliers, batch effects, missing values or normalisation.
---

# QC and normalisation per omic layer

Run (all layers, or name one):
```bash
Rscript .claude/skills/prep-qc-normalisation/scripts/qc_normalise.R study_config.json [layer]
```
Outputs: `results/data/<layer>_processed_se.rds` (assay `values`, log2 scale;
counts layers also keep `counts` + TMM factors), `results/qc/<layer>_qc_summary.txt`,
`results/qc/<layer>_pca.pdf`.

## Defaults by `scale` (override per layer in study_config.json)
| scale | filter | impute | normalise | output |
|---|---|---|---|---|
| counts | `edgeR::filterByExpr` with group × visit cells | none | TMM | log2-CPM (prior 2) |
| log2 | visit-aware missingness (`missing_max`, 0.3) | half-min (min − 1 on log2) | none (`median` for MS) | log2 |
| raw_intensity | pooled-QC RSD ≤ 0.3 if QC samples, then visit-aware missingness | half-min | PQN | log2 |

Per-layer optional fields: `missing_max`, `impute` (`halfmin` | `minprob` | `none`),
`normalise` (`tmm` | `median` | `pqn` | `none`), `qc_rsd_max`.

**Visit-aware filter**: a feature is kept if it is adequately measured in at
least one visit. This deliberately retains features that appear or disappear
over time; a global missingness filter would remove exactly those.

## Choosing options — ask the user when unclear
- Olink NPX / SomaScan: already normalised by the vendor → `normalise: none`.
  Olink values below LOD: ask whether LOD-flagged values are kept or set to NA.
- MS proteomics (log2 LFQ): `normalise: median`, `impute: minprob`
  (missing-not-at-random). MinProb draws random values — seeded, but note
  this in methods and consider a sensitivity analysis without imputation.
- Metabolomics: PQN is the default for dilution effects (urine, sputum, BAL);
  for plasma the user may prefer `median`. Half-minimum is conventional but
  biases variance downward for heavily-missing features.

## Reading the diagnostics
`*_qc_summary.txt` ends with p-values for PC1–PC5 vs visit, group, batch and
covariates. Report:
- PCs associated with **batch** → add batch to `covariates` in the config so
  it enters the models. Do NOT run `removeBatchEffect`/ComBat on the processed
  data unless batch is fully crossed with visit and group; if batch is
  confounded with visit (e.g. each visit run as a separate plate), say so
  plainly — it cannot be separated statistically.
- Outlier samples on the PCA: identify by `sample_id`, show the user, and let
  them decide. Never drop samples without confirmation; removing one sample in
  a 4-visit design breaks that subject's trajectory.
- Subject effects dominating PC1 are expected in repeated-measures data.
