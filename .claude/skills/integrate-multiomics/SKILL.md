---
name: integrate-multiomics
description: Integrate transcriptomics, proteomics and metabolomics measured at repeated visits — MOFA2 latent factors tested with mixed models (time, arm, time × arm) and mixOmics DIABLO with multilevel (within-subject) decomposition to find coordinated multi-omic signatures that change across visits. Use when the user asks for joint/integrated multi-omics analysis, shared factors, cross-omic signatures or networks over time.
---

# Multi-omics integration over visits

Both scripts were written against the documented MOFA2 / mixOmics APIs but
were **not executed** where this kit was built (packages unavailable there);
the downstream factor-testing code was tested. On first use, run on a small
subset (`top_var = 300`, few factors) and fix any API drift before the full run.

## MOFA2 (default, unsupervised; tolerates missing layers)
```bash
Rscript .claude/skills/integrate-multiomics/scripts/integrate_mofa.R study_config.json [n_factors=15] [top_var=2000] [layers=all]
```
- All subject × visit samples in one model; a sample missing a layer is kept.
- Top-variance features per layer, views scaled, seeded training.
- Each factor tested with `lme4` likelihood-ratio tests:
  `Factor ~ visit * group + covariates + (1 | subject_id)` → `p_visit`,
  `p_group`, `p_visit_x_group` + BH FDR (`factor_tests.csv`).
- `variance_explained.csv`: which layers drive each factor — a factor active
  in ≥ 2 layers is the cross-omic signal; single-layer factors are layer-
  specific variation.
- Requires Python via basilisk (`use_basilisk = TRUE`) on first run; this
  downloads a conda environment — warn the user it takes time.

Interpretation: factors associated with visit (and especially visit × group)
are the longitudinal multi-omic axes. Pass a factor's score table to
LongOmicsTraj (treat factors as features of a layer named `mofa_factors`)
to label factor trajectories. Inspect top weights per layer in
`weights_<layer>.csv` and run pathway enrichment on them.

Alternative not scripted here: MEFISTO (`set_covariates(mofa, "time")`) for
smooth temporal factors when visit timing varies across subjects.

## DIABLO (supervised; complete cases only)
```bash
Rscript .claude/skills/integrate-multiomics/scripts/integrate_diablo.R study_config.json [group|all] [ncomp=3] [keepX=20] [top_var=2000]
```
- Y = visit; each block passed through `withinVariation()` (subject design)
  so the model contrasts visits within subjects.
- Uses only samples present in all layers and subjects with ≥ 2 such samples
  — report how many were lost.
- Run per arm and compare selected features; pooling arms mixes treatment
  response with natural time effects.
- keepX is not tuned in the script. Before reporting, tune with
  `mixOmics::tune.block.splsda` and assess with `perf()`. Standard Mfold CV
  does not keep a subject's visits together, so performance estimates are
  optimistic; say this if reporting classification error.
- Design matrix off-diagonal 0.1 (favours discrimination over inter-block
  correlation); justify any change.

## Choosing
MOFA2 when the question is "what are the main coordinated sources of
variation and do they change over time / differ by arm" and data are
incomplete. DIABLO when the question is "which multi-omic signature best
separates visits" and complete-case loss is acceptable. Report both only if
the user wants convergent evidence; do not cherry-pick.
