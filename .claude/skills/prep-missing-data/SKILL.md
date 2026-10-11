---
name: prep-missing-data
description: Missing-data handling for longitudinal omics and clinical data — diagnosis (missingness per feature, sample, visit and group; intensity–missingness correlation as an MNAR / left-censoring signature; missingness vs visit, arm and batch; subject follow-up patterns), masking benchmark of imputation methods (half-minimum, MinProb, kNN, BPCA, mean) under MCAR and MNAR, writing an imputed layer, and multiple imputation with mice plus Rubin-pooled regression vs complete-case. Use when the user asks about missing values, NAs, values below detection, imputation, dropout, incomplete visits, MCAR/MAR/MNAR, or multiple imputation.
---

# Missing data

Category: **Missing data**. Run after `prep-study-setup`. The diagnosis then
decides the imputation in `prep-qc-normalisation` (`impute:` per layer) or here.

```bash
S=.claude/skills/prep-missing-data/scripts/missing_data.R
Rscript $S --mode=diagnose  --layer=proteomics [--zero_as_na=TRUE] [--clinical=feno,ige,acq5]
Rscript $S --mode=benchmark --layer=proteomics [--methods=halfmin,minprob,knn,mean,bpca] [--frac=0.05] [--reps=3]
Rscript $S --mode=impute    --layer=proteomics --method=knn [--max_missing=0.5] [--normalise=none|median]
Rscript $S --mode=mice --visit=V1 --formula="fev1 ~ age + sex + feno" [--log=feno] [--m=20]
```
Common options: `--config`, `--seed`, `--tag`. Outputs go to `<output_dir>/missing/<mode>_<tag>/`.
`impute` writes a new layer `<layer>_imp_<method>`, which any skill can use via `--layer=`.

## Decision guide
| what diagnose shows | missingness type | handle with |
|---|---|---|
| strong negative intensity–missingness correlation; values below LOD / zeros in MS tables | **MNAR** (left-censored) | halfmin or minprob, or a censored-data model. MAR methods overestimate these values |
| no intensity link, random spread | MCAR / MAR | knn or bpca (benchmark); complete-case is unbiased only under MCAR |
| sample missingness differs by visit, arm or batch | informative by design | report it and keep that variable in the model; do not impute across arms differently |
| whole visits missing (dropout, intermittent) | MAR at subject level | **do not impute whole visits**: use `stats-mixed-models` (likelihood-based, valid under MAR). LOCF is biased and is not offered |
| clinical covariates missing for some subjects | MAR | `--mode=mice` (PMM, m ≥ 20, Rubin pooling); report it next to complete-case |

## Rules
- MinProb and mice are random. Report the seed, and for important results rerun with another seed (or more `--m`).
- The benchmark masks *observed* values. Its MNAR scenario hides each feature's lowest 20%, which approximates left-censoring but cannot prove the true mechanism. Report both error (NRMSE) and bias: negative bias means imputed values are too low.
- Single imputation (halfmin, knn, …) understates uncertainty. Downstream p-values treat imputed values as observed, so say how many values were imputed (`missing_mask`).
- Filter before imputing: features missing in > 50% of samples (`--max_missing`) carry little information. The visit-aware filter in `prep-qc-normalisation` keeps features that appear or disappear over time; impute those with halfmin/minprob, not kNN.
- `bpca` (pcaMethods) is accurate when features are correlated, but it is slow above roughly 1,000 features, so it is off by default in the benchmark.
- Counts layers (RNA-seq) are refused, because a zero count is an observation, not a missing value.
