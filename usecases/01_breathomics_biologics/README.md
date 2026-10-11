# Use case 01 — Exhaled VOCs before and after biologics (breathomics, paired sPLS-DA)

Reproduces the workflow of **Vassileva SM *et al.* The impact of benralizumab
and dupilumab on exhaled volatile organic compound profiles of severe asthma
patients – a preliminary analysis. *Respir Med* 2026;262:109067**
(doi:10.1016/j.rmed.2026.109067; 3TR-ABC, NCT06126692) with DMP-AE skills, on
**synthetic data**, because the patient-level data are available only on request.

## Contents
```
README.md                     this file
study_config.json             design: visits BL / M1 / M4, group = biologic,
                              batch = gcms_batch, layer breath_voc (+ annotation)
breathomics_biologics.Rmd     the notebook (Knit in RStudio)
breathomics_biologics.html    rendered notebook
data/simulate_breath_data.R   synthetic cohort generator (seeded)
data/metadata.csv             65 samples: clinical data per visit
data/breath_gcms_raw.csv      4270 fragments x 65 samples (XCMS-style table)
data/feature_annotation.csv   m/z, RT, compound, class, contaminant flag
data/planted_truth.csv        the effects built into the simulation
results/                      everything the skills write (created by the notebook)
```

## Run it
- **RStudio:** open the project, open `breathomics_biologics.Rmd`, then **Knit**.
  A run from scratch takes about 10 minutes, mostly the permutation null; the
  `n_perm` parameter sets the number of permutations. Steps with results on
  disk are reused; set `rerun: true` to recompute them.
- **Claude / MCP:** `run_usecase("01_breathomics_biologics")`, then
  `panel_report("usecases/01_breathomics_biologics/results/splsda/BLvsM4_breath_voc")`.
- **Step by step:** each notebook chunk is one `Rscript .claude/skills/...` call
  with `--config=usecases/01_breathomics_biologics/study_config.json`.

## Paper step → engine skill
| paper | skill | status |
|---|---|---|
| clinical tables, Table 1 / 3, Fig. 2 | `stats-clinical-association` (table1, change) | **new** |
| XCMS table → contaminants → log → batch (limma), Suppl. Fig. 1 | `prep-breathomics-gcms` | **new** |
| multilevel sPLS-DA, stepwise elimination, AUROC, LOSO-CV, Fig. 3–4, Table 2 | `ml-multilevel-splsda` (+ `R/lomics_spls.R`) | **new** |
| Kendall correlations (Table 4) | `stats-clinical-association --mode=correlation` | **new** |
| medication Mann–Whitney | `stats-clinical-association --mode=exposure` | **new** |
| VOC change by biologic | `stats-group-testing --design=between --change_from=BL` | existing |
| raw chromatogram processing (XCMS) | — | out of scope (starts from the peak table) |

New MCP tools: `list_usecases`, `run_usecase`, `panel_report`.

## Added beyond the paper
- **Nested leave-one-subject-out CV:** keepX tuning and elimination are
  repeated inside each fold, giving an honest performance estimate. Results
  are reported as paired accuracy.
- **Within-patient permutation null:** shows the training and fixed-panel LOSO
  AUROC that the same pipeline produces when there is no effect.
- **Re-selection frequency:** how often each panel VOC is selected again
  across the nested folds.

## Caveats
- All numbers describe the simulation. Planted effects and the seed set
  them, and they are not estimates for the 3TR-ABC cohort.
- mixOmics was not available where this was built, so the internal sPLS-DA
  engine was used; it is checked against `pls::plsr` in
  `tests/test_spls_internal.R`. With mixOmics installed, the skill uses it
  automatically, and results can differ slightly.
