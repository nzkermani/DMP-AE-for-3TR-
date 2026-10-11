---
name: stats-mixed-models
description: Linear mixed-effects models (lme4) for longitudinal data with unbalanced or missing visits — per omic feature or per clinical variable, y ~ visit (or continuous time) × group + covariates + random intercept (or slope) per subject; likelihood-ratio tests for the time effect, time × group interaction and time within each group, BH-FDR, ICC, estimated marginal means and change-vs-baseline contrasts with CI (emmeans), trajectory plots. Use when the user asks for mixed models, LMM, lmer, random effects, repeated measures with missing visits, change over time in clinical scores or omics, or a treatment × time interaction.
---

# Mixed-effects models for longitudinal data

Category: **Longitudinal & mixed-effects models**. Unlike paired tests and multilevel sPLS-DA, it uses **every
visit of every subject**: subjects with a missing visit stay in the model.

```bash
S=.claude/skills/stats-mixed-models/scripts/mixed_models.R
# omics: every feature (or --features= / --panel=), visit as factor, arm interaction
Rscript $S --config=study_config.json --mode=features --layer=breath_voc --by=group
# clinical variables, log-scale for skewed ones
Rscript $S --mode=clinical --vars=blood_eos,feno,ige,fev1_pp,acq5 --log=blood_eos,feno,ige --by=group
# options: --time=factor|continuous  --random=intercept|slope  --covariates=age,sex
#          --visits=BL,M1,M4  --report_top=20  --fdr=0.05  --tag=name
```
Runtime is about 0.1 s per feature (about 5 min for 3,000 features with `--by`).
Outputs go to `<output_dir>/mixed/<mode>_<tag>/`: `lmm_results.csv`, `emmeans.csv`,
`contrasts.csv`, `trajectories.(pdf|png)` and `summary.txt`.

## Model and tests
| quantity | how |
|---|---|
| overall time effect | ML likelihood-ratio test: full model vs a model without time terms (including the interaction) |
| time × group | LRT: full model vs additive model |
| time within each group | separate LRT fitted in each group |
| estimates | REML refit; `emmeans` per visit (× group) and contrasts vs the first visit, with 95% CI (Kenward-Roger df when pbkrtest is installed) |
| ICC | random-intercept variance / total; high ICC means stable subject differences, which is why paired designs gain power |

The FDR is BH across features or variables. Contrasts are unadjusted within a variable.

## Rules
- Use `--time=factor` for few, irregular visits; it makes no shape assumption. Use `continuous` (time_values from the config) for a linear trend, and add `--random=slope` only with ≥ 3 visits per subject.
- Report the interaction as *the* treatment-difference test. "Significant in arm A and not in arm B" is not evidence of a difference between arms.
- Small or unbalanced groups (e.g. 5 vs 20) make interaction p-values unstable, and ML LRTs are slightly anticonservative at n ≈ 25. Treat interaction hits as exploratory and check the trajectory plot.
- Report singular fits (`singular = TRUE`). They usually mean the random-effect variance is close to 0, so the model reduces to ordinary regression.
- MAR assumption: subjects with missing visits are retained, which is valid when missingness does not depend on unobserved values. Say why visits are missing.
- For thousands of features, `stats-longitudinal-de` (limma with duplicateCorrelation, using one shared correlation) is faster. Use this skill for per-feature random effects and interpretable per-variable models.
