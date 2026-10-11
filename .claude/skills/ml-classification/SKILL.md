---
name: ml-classification
description: Supervised classification of a binary or multi-class outcome from omics — PLS-DA (with VIP), sparse PLS-DA, PCA-LDA, elastic-net logistic regression, random forest and linear SVM compared with repeated nested cross-validation (AUC, balanced accuracy, sensitivity, specificity, Brier), subject-grouped folds, label-permutation null and final-model importances. Use when the user asks for classification, prediction, PLS-DA, LDA, a classifier, AUC/ROC or which model predicts an outcome best.
---

# Classification

Category: **Machine learning**.

```bash
S=.claude/skills/ml-classification/scripts/classification.R
Rscript $S --outcome=binary --layers=all --visit=V1 \
  [--models=plsda,splsda,lda,enet,rf,svm] [--folds=5] [--repeats=10] [--prefilter=200]
# null distribution: same call with labels shuffled between subjects
for i in $(seq 1 20); do Rscript $S --outcome=binary --repeats=2 --permute_seed=$i; done
```
Outputs `results/classification/<outcome>_<tag>/`: `cv_summary.csv`,
`cv_metrics.csv`, `oof_predictions.csv`, `roc.pdf`,
`final_importance_<model>.csv` (PLS-DA VIP, enet coefficients, RF
importance), `plsda_scores.pdf`.

## Models
| model | type | tuned (inner CV) |
|---|---|---|
| plsda | PLS-DA (pls; mixOmics if installed) | ncomp 1–5 |
| splsda | sparse PLS-DA (mixOmics only) | keepX |
| lda | PCA → LDA (MASS) | number of PCs |
| enet | elastic-net logistic/multinomial (glmnet) | lambda |
| rf | random forest (ranger) | — |
| svm | linear SVM (e1071) | cost |

## Validation design (do not weaken)
- Outer: stratified K-fold × R repeats. Inside each training fold:
  scaling, univariate prefilter (top `--prefilter` by limma), PCA and
  hyperparameter tuning. Test folds never influence the model.
- `--visit=all`: folds are grouped by subject (no subject in both train
  and test).
- Tested: on 6 random-label outcomes, mean AUC was 0.48–0.53 per model.
  Single null runs ranged 0.37–0.63 at n = 58, so **one** AUC around 0.6
  is not evidence. Use the permutation loop and report an empirical
  p-value: (1 + #null ≥ observed) / (1 + #null).

## Reporting rules
- Report mean (SD) over repeats for AUC and balanced accuracy, chance
  accuracy (majority class), n per class, and the permutation p-value.
- Compare models on the same folds (they are). Differences within one SD
  are not meaningful; prefer the simpler model.
- `final_importance_*` and `plsda_scores.pdf` come from a model fit on
  all data. They are for interpretation only. The PLS-DA score plot shows
  apparent separation, which is always optimistic. Never present it as
  validation.
- Performance generalises only to the population sampled. An external
  cohort is needed before any claims about a clinical biomarker.
- Class imbalance: use balanced accuracy and AUC, not accuracy.
