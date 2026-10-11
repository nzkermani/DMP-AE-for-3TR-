---
name: stats-longitudinal-de
description: Repeated-measures differential analysis of omic features across four visits with limma (duplicateCorrelation / voom) or variancePartition dream — change from baseline, visit-to-visit steps, group × visit interaction and overall time F-tests per arm. Use when the user asks which genes, proteins, metabolites or pathways change over time, differ between arms over time, or respond to treatment.
---

# Longitudinal differential analysis

```bash
Rscript .claude/skills/stats-longitudinal-de/scripts/longitudinal_de.R study_config.json [layer|all] [limma|dream]
```
Outputs in `results/de/<layer>/`: one CSV per contrast (limma `topTable`
columns: logFC, AveExpr, t, P.Value, adj.P.Val, B), `model_spec.txt`, and
`results/de/de_summary.csv` (counts at FDR < 0.05; re-running one layer keeps
the others).

## Model
`~ 0 + cell + covariates`, `cell` = group × visit (cell-means coding), so
every contrast is explicit:
- `change_<group>_<Vk>vsV1` — within-arm change from baseline
- `step_<group>_<Vk>vs<Vk-1>` — consecutive transition (V3vsV2, V4vsV3; V2vsV1 is the first change contrast)
- `interaction_<G2>vs<G1>_<Vk>` — difference between arms in change from baseline (the treatment-effect contrast for an RCT)
- `overall_time_<group>` — F-test across all change contrasts in that arm (limma only)

Within-subject correlation:
- **limma** (default): `duplicateCorrelation(block = subject_id)`; for counts,
  the standard voom → duplicateCorrelation → voom → duplicateCorrelation
  sequence. One consensus correlation for all features. `trend = TRUE`,
  `robust = TRUE` for non-count layers.
- **dream**: per-feature `(1 | subject_id)` mixed model (Satterthwaite df
  by default). Prefer when follow-up is unbalanced (many missing visits),
  or for extra random terms. Slower. Requires `variancePartition`.

## Interpretation rules
- An interaction contrast, not "significant in arm A but not in arm B", is the
  evidence for a differential treatment response. Point this out if the user
  compares lists of significant features between arms.
- FDR is controlled within each contrast; across many contrasts × layers,
  report this structure in methods.
- Covariates that are constant within subject (sex, baseline age) are
  absorbed largely by the subject correlation; they mostly matter for the
  between-arm comparison at baseline.
- If within-subject correlation is near 0 or negative, check sample labelling
  (swapped visits/subjects) before trusting results.
- The pathway layer `pathways_<layer>` (from interpret-pathway-enrichment) can be run
  through this same script.

## Before reporting numbers
Read `de_summary.csv` and the relevant CSV — never quote counts or feature
names from memory. Use MCP `feature_report` to show a single feature's full
profile across contrasts.
