---
name: stats-clinical-association
description: Clinical endpoint statistics for small longitudinal cohorts — Table 1 with normality-based summaries, within-subject change in clinical variables between two visits (Shapiro-Wilk → paired t-test or Wilcoxon signed-rank with Hodges-Lehmann CI), Kendall rank correlation of selected features with clinical variables (tau, z, p, bootstrap CI), and feature levels by medication exposure (Mann-Whitney / Kruskal-Wallis). Use for baseline characteristics tables, "did FeNO / ACQ / FEV1 change after treatment", correlations of a biomarker panel with clinical scores, or checking medication confounding.
---

# Clinical change, correlation and exposure checks

Category: **Statistics**.

```bash
S=.claude/skills/stats-clinical-association/scripts/clinical_association.R
Rscript $S --mode=table1 --visit=BL [--vars=age,sex,...] [--by=biologic]
Rscript $S --mode=change --vars=blood_eos,feno,ige,fev1_pp,acq5,aqlq --from=BL --to=M4 [--show_visits=BL,M1,M4]
Rscript $S --mode=correlation --layer=breath_voc --panel=<run>/panel.csv --vars=acq5,blood_eos,ige,feno,fev1_pp --visit=BL [--boot=2000]
Rscript $S --mode=exposure --layer=breath_voc --panel=<run>/panel.csv --exposures=ocs_maintenance,lama,ltra,ics_dose --visit=BL
```
Common options: `--config`, `--adjust=none|BH`, `--tag`, `--seed`.
Outputs `<output_dir>/clinical/<mode>_<tag>/` (`table1.csv`, `change.csv` +
`clinical_boxplots`, `kendall.csv` + `kendall_table.csv`, `exposure_tests.csv`, `summary.txt`).

## Choices made by the script (report them)
- Normality: Shapiro-Wilk p > 0.05 → mean (SD) / paired t. Otherwise → median [IQR] / Wilcoxon signed-rank. At n ≈ 25 the test has low power; say which test was used for each variable.
- The Wilcoxon estimate is the Hodges-Lehmann pseudo-median of differences, with its CI. `median_of_differences` is the plain median. Papers often report either; say which.
- Change uses complete pairs only. n per variable is in `n_pairs`.
- Kendall: `cor.test(..., exact = FALSE)`, so z is the normal-approximation statistic (as in R's stats). The tau CI is a percentile bootstrap; the source paper does not state its CI method.
- p-values are unadjusted by default (exploratory). With 5 features × 5 clinical variables, about 1 of 25 tests is expected below 0.05 by chance. Say so, or use `--adjust=BH`.
- Exposure: numeric exposures are split at the median. The test is between patients at one visit, so it cannot separate drug effects from disease severity.
