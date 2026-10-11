---
name: stats-group-testing
description: Classical statistical tests on omic features and clinical variables with FDR control — t-test / Wilcoxon, ANOVA / Kruskal-Wallis between groups at a visit or on change scores, paired t / signed-rank between two visits, Friedman and mixed-model LRT across all visits, and a clinical characteristics (Table 1) with chi-square / Fisher. Use when the user asks for statistical testing, group comparisons, p-values, paired comparisons or a baseline characteristics table.
---

# Statistical testing

Category: **Statistics**. Needs `prep-qc-normalisation` outputs.

```bash
S=.claude/skills/stats-group-testing/scripts/group_tests.R
Rscript $S --design=between  --by=group --visit=V1                   # independent groups
Rscript $S --design=between  --by=group --visit=V3 --change_from=V1  # change scores between arms
Rscript $S --design=paired   --visits=V1,V3 --subset_col=group --subset=Treatment
Rscript $S --design=repeated --subset_col=group --subset=Treatment   # all 4 visits
Rscript $S --design=clinical --by=group --visit=V1                   # Table 1
```
Add `--layers=proteomics` to restrict layers. Outputs in
`results/tests/<design>_<tag>/`: `<layer>_tests.csv`, `summary.txt`
(or `clinical_table.csv`).

## What each design runs
| design | parametric | non-parametric | effect size |
|---|---|---|---|
| between, 2 groups | Welch t | Wilcoxon rank-sum | mean difference, Cohen's d |
| between, >2 groups | Welch one-way ANOVA | Kruskal-Wallis | eta² |
| paired (2 visits) | paired t | Wilcoxon signed-rank | mean paired difference, d_z |
| repeated (all visits, complete cases) | mixed-model LRT `visit + (1|subject)` | Friedman | max visit-mean range |
| clinical | t / ANOVA, chi-square | Kruskal / Fisher | mean (SD) or n (%) |

BH FDR is applied within each layer and each test. Values are on the
processed log2 scale, so mean differences are log2 fold changes.

## Choosing and reporting
- For omics, the moderated statistics in `stats-regression` (limma) and
  `stats-longitudinal-de` are more powerful at small n. Use this skill
  when classical tests are explicitly wanted, for non-parametric
  confirmation, or for clinical variables.
- Pick parametric vs non-parametric **before** looking at results; report
  one as primary and the other as sensitivity. Never choose per feature
  by whichever p is smaller.
- `repeated` drops subjects missing any visit (count reported). For
  incomplete follow-up, prefer `stats-longitudinal-de` (all samples).
- Testing "within Treatment" and "within Placebo" separately is not a test
  of a treatment effect. Use `between` on change scores, or the
  interaction contrasts in `stats-longitudinal-de`.
- Clinical table p-values are unadjusted and descriptive. In an RCT,
  baseline p-values between randomised arms are not meaningful (CONSORT);
  say so if asked.
