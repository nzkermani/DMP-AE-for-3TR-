---
name: stats-regression
description: Linear and logistic regression for omics and clinical data — omics-wide association of every feature with a continuous or categorical outcome adjusted for covariates (limma), per-feature outcome-on-feature linear/logistic models with odds ratios, and single multivariable lm/glm models with confidence intervals, VIF and diagnostics. Use when the user asks for linear regression, logistic regression, association with a clinical variable (e.g. FEV1, response), adjusted models or multivariable models.
---

# Regression (linear / logistic)

Category: **Statistics**. Needs `prep-qc-normalisation` outputs.

```bash
S=.claude/skills/stats-regression/scripts/regression.R
# 1. omics-wide: feature ~ outcome + covariates (limma, moderated t, BH)
Rscript $S --mode=features --outcome=continuous --visit=V1
# 2. outcome ~ feature + covariates, one model per feature (OR or beta per 1 SD)
Rscript $S --mode=outcome  --outcome=binary --layers=proteomics
# 3. one multivariable model
Rscript $S --mode=formula --layers=proteomics --formula="fev1 ~ age + sex + PROT21" [--family=binomial]
```
`--outcome` accepts a config key (`continuous`, `binary`) or any metadata
column. Covariates default to `covariates` in the config. Override with
`--covariates=age,sex,batch`, or use `--no_covariates`.
`--change_from=V1` uses change scores. `--visit=all` (mode features)
blocks on subject via `duplicateCorrelation`.

Outputs `results/regression/<mode>_<outcome>_<tag>/`: `<layer>_regression.csv`,
`summary.csv`; formula mode: `coefficients.csv`, `model_summary.txt`,
`diagnostics.pdf`.

## Which mode
- **features**: the standard omics association scan. The coefficient
  (`logFC`) is the change in feature (log2) per unit of outcome, or the
  difference between outcome levels. Most powerful at small n.
- **outcome**: when the clinical variable is naturally the response. It
  gives odds ratios or betas per 1 SD of each feature, which is useful for
  candidate biomarkers. Use only on one visit (it warns otherwise).
- **formula**: confirmatory multivariable models on a few pre-chosen
  variables. Multi-layer feature names are `<layer>__<feature>`. The
  output includes VIF (>5 means collinear) and residual diagnostics. Check
  that the plots look reasonable before reporting.

## Rules
- Don't build a formula from features selected by p-value on the same
  data and then report its p-values or R² as confirmatory. That is
  double-dipping. Say it is exploratory, or validate in held-out data.
- Logistic models: keep events per variable ≥ 10. With complete or
  quasi-complete separation, glm estimates explode — use penalised models
  (`ml-classification` enet) instead.
- Report adjusted estimates with 95% CI and the covariate set used.
