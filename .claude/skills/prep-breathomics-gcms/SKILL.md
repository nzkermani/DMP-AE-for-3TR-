---
name: prep-breathomics-gcms
description: Preprocess an exhaled-breath GC-MS (breathomics, VOC) feature table — known-contaminant removal (siloxanes, phthalates, sorbent artefacts) from an annotation file, optional zero filter, log2 with half-minimum for zeros, limma::removeBatchEffect for GC-MS analysis batches with the visit effect preserved (refuses when batch is confounded with visit), PCA before/after. Use instead of prep-qc-normalisation for breath VOC / TD-GC-MS / eNose-style layers, or when the user asks about contaminants or batch correction of breath data.
---

# Breathomics (GC-MS) preprocessing

Category: **1 Data preparation** (replaces `prep-qc-normalisation` for breath layers).
Run `prep-study-setup` first so `results/data/<layer>_raw_se.rds` exists.

```bash
Rscript .claude/skills/prep-breathomics-gcms/scripts/breath_preprocess.R \
  --config=study_config.json --layer=breath_voc \
  [--annotation=data/feature_annotation.csv --contaminant_col=contaminant] \
  [--contaminant_regex="siloxane|phthalate"] [--zero_max=1] \
  [--batch=gcms_batch|none] [--preserve=visit] [--force_batch=FALSE]
```
The annotation path can also be set as `"annotation"` on the layer in
`study_config.json`. Batch defaults to the config's `batch_column`.

Outputs `data/<layer>_processed_se.rds` (assay `values` = log2 batch-corrected,
`log2_uncorrected`; rowData = annotation: compound, m/z, RT, class),
`qc/<layer>_breath_qc.txt` (counts per step, batch × visit table, Cramér's V,
PC1–5 association p-values before/after), `qc/<layer>_breath_pca.(pdf|png)`.

## Steps and defaults
| step | default | note |
|---|---|---|
| raw chromatograms (XCMS) | **not done here** | start from the aligned peak table (fragments × samples) |
| contaminants | annotation column `contaminant` (TRUE/FALSE) and/or `--contaminant_regex` on `compound` | report fragments before → after |
| zero filter | off (`--zero_max=1`) | `--zero_max=0.5` keeps a fragment if ≤ 50 % zeros in at least one visit |
| zeros | half the fragment's minimum positive value, then log2 | state this in methods |
| batch | `limma::removeBatchEffect(batch, design = ~ visit)` | stops if batch is nested in visit; warns if Cramér's V > 0.5 |

## Rules
- The script **refuses** to batch-correct when each batch holds one visit (or vice versa); the correction would remove the effect of interest. Then put batch in the model instead, and say plainly that batch and visit cannot be separated.
- If batches are nested in site, correction also removes site differences — the QC file says so; mention it.
- `removeBatchEffect` output is for exploration and for models without batch terms. For linear models, prefer batch as a covariate on `log2_uncorrected`.
- z-scoring is done in the analysis step (`ml-multilevel-splsda --zscore=TRUE`), not here.
- Ask before removing extra samples (failed runs) — show them on the PCA first.
