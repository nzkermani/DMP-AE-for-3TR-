---
name: stats-survival
description: Survival / time-to-event analysis — Kaplan-Meier with log-rank by arm, clinical group or omics-derived clusters, per-feature Cox models with hazard ratios and proportional-hazards checks, elastic-net penalised Cox with cross-validated C-index and risk groups, and multivariable Cox models. Use when the user asks about survival, time to exacerbation / event, hazard ratios, Kaplan-Meier, Cox regression or prognostic biomarkers.
---

# Survival analysis

Category: **Statistics**. Needs `prep-qc-normalisation` outputs and
`survival_time` / `survival_event` in `study_config.json` → `outcomes`
(event coded 1 = event, 0 = censored).

```bash
S=.claude/skills/stats-survival/scripts/survival.R
Rscript $S --mode=km --strata=group
Rscript $S --mode=km --strata_file=results/clustering/<tag>/assignments.csv --strata_col=cluster
Rscript $S --mode=cox_features  --layers=all                     # HR per SD, FDR, PH test
Rscript $S --mode=cox_penalised --layers=proteomics --repeats=10 # CV C-index, risk tertiles
Rscript $S --mode=cox_formula   --layers=proteomics --formula="age + sex + PROT21"
```
Options: `--time=`, `--event=` override the config. `--visit=V1` is the
default. `--covariates=`, `--alpha=0.5`, `--folds=5`.
Outputs `results/survival/<mode>_<tag>/`: `summary.txt`, `km*.pdf`,
`cox_per_feature.csv`, `selected_features.csv`, `cv_risk_scores.csv`,
`ph_diagnostics.pdf`.

## Rules
- **Predictors must be measured at or before time zero.** Use baseline
  (V1) by default. Using V3 omics to predict events counted from V1 is
  immortal-time bias. Do a landmark analysis instead (recompute time from
  V3, exclude earlier events) and ask the user for the event dates.
- Check proportional hazards (`ph_test_p`, `cox.zph`). If violated,
  report it and consider stratification or time-varying effects.
- Events per variable: per-feature Cox with k covariates needs about
  10 × (k + 1) events. The script warns when below that.
- The penalised-Cox C-index is **out-of-fold** (honest). The selected
  feature list comes from one full-data fit and is unstable at small n.
  For a robust list, run `ml-feature-selection --outcome=survival`
  (stability selection).
- Risk tertiles are built from out-of-fold scores. KM curves of in-sample
  risk groups would be optimistic; never present those as validation.
- For KM by clusters, the clusters must not have been derived using the
  outcome.
