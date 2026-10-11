---
name: prep-study-setup
description: Set up a longitudinal multi-omics study (transcriptomics, proteomics, metabolomics at repeated visits) for analysis in R — writes study_config.json, aligns each omic matrix to the sample metadata, audits the visit design and saves SummarizedExperiment / MultiAssayExperiment objects. Use first, before QC, differential, trajectory or integration analysis, or whenever new data files or metadata arrive.
---

# Omics study setup (4-visit multi-omics)

## What this step produces
`results/data/<layer>_raw_se.rds`, `results/data/raw_mae.rds`,
`results/setup/design_audit.txt`, `results/setup/sample_overlap.csv`.
Every later skill reads these, so nothing else runs until this succeeds.

## Data contract (check before writing any config)
- **Metadata**: one row per sample (subject × visit). Needs sample, subject,
  visit and group columns; any names are fine, they are mapped in the config.
  One sample per subject per visit — collapse technical replicates first.
- **Each omic layer**: matrix with features in rows, samples in columns, first
  column = feature IDs, header = sample IDs (CSV/TSV), or an `.rds` matrix /
  SummarizedExperiment.
- `scale` per layer is one of `counts` (RNA-seq raw counts), `log2` (Olink NPX,
  SomaScan log2 RFU, MS log2 intensities), `raw_intensity` (untransformed
  metabolomics / MS peak areas). Getting this wrong silently breaks QC — if the
  user has not said, inspect value ranges (`summary(as.vector(x))`) and ASK.

## Steps
1. Ask for (or locate) the metadata file and one file per layer. Confirm:
   visit labels and their chronological order, the real time of each visit
   (e.g. months: `{"V1":0,"V2":3,"V3":6,"V4":12}`; unequal spacing matters for
   the GAM/maSigPro estimators), the group/arm column, covariates, batch column,
   and whether pooled-QC injections exist (set `"qc_sample_column"`).
2. Edit `study_config.json` at the project root (template already present).
3. Run:
   ```bash
   Rscript .claude/skills/prep-study-setup/scripts/build_layers.R study_config.json
   ```
   (or MCP tool `run_skill_script` with skill `prep-study-setup`).
4. Read `results/setup/design_audit.txt` and report to the user:
   samples per group × visit for each layer, number of subjects with all four
   visits, missingness, and how many samples have all three layers (this bounds
   DIABLO / complete-case integration).

## Things to flag rather than fix silently
- Sample IDs that appear in a matrix but not the metadata (dropped, count shown).
- Visit labels not in `visits` (script stops) — never guess the ordering.
- Layers measured on different subsets of visits (e.g. metabolomics only at
  V1 and V4): trajectory signatures need ≥ 3 visits to be informative; tell
  the user which analyses that layer can join.
- Feature IDs that will need mapping for pathway analysis (Ensembl → symbol,
  UniProt/Olink ID → gene symbol, metabolite names → HMDB/KEGG).
