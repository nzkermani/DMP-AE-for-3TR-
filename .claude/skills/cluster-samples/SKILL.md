---
name: cluster-samples
description: Unsupervised clustering of subjects or samples for endotype / phenotype discovery — k-means, PAM, Ward hierarchical, Gaussian mixture (mclust) and consensus clustering on one or several omic layers, with k selection (silhouette, BIC, PAC), bootstrap stability, method agreement (ARI), clinical and survival characterisation and cluster-driving features. Use when the user asks to cluster patients, find endotypes, subgroups, or compare clustering methods.
---

# Sample clustering (endotypes)

Category: **Clustering** (sample-level). For clustering features by their
trajectory shape, use `cluster-trajectories-longomicstraj`.

```bash
Rscript .claude/skills/cluster-samples/scripts/cluster_samples.R \
  --layers=transcriptomics,metabolomics --visit=V1 \
  [--methods=kmeans,pam,hclust,gmm,consensus] [--k=2:6] [--top_var=1000] \
  [--n_pcs=10] [--boot=50] [--seed=42]
```
Outputs `results/clustering/<tag>/`: `k_selection.csv`,
`assignments.csv`, `stability.csv`, `method_agreement_ARI.csv`,
`clinical_association.csv`, `top_features_by_cluster.csv`, `plots.pdf`,
`consensus/`.

## How it works
Data are autoscaled (block-scaled for >1 layer) and reduced to the top PCs.
Each method chooses k by silhouette (GMM by BIC, consensus by lowest PAC).
Stability is mean bootstrap Jaccard per cluster: >0.75 is stable, <0.6
means the cluster dissolves. The **recommended** solution in
`assignments.csv$cluster` is the method with highest mean stability.
Always show the user the agreement table too.

## Reporting rules
- Report k-selection evidence, stability, ARI between methods, and
  cluster sizes. A cluster found by one method only, or with Jaccard < 0.6,
  is not a finding.
- `clinical_association.csv` tests clusters against every metadata
  variable (BH across variables). A cluster that tracks batch, sex or
  visit is a technical/confounding axis, not an endotype. Check this
  first.
- `top_features_by_cluster.csv` (one-vs-rest Wilcoxon) is circular:
  features differ between clusters because they were used to make them.
  Use it to describe clusters, not as evidence that the clusters exist.
- Clusters are defined at one visit. To see how clusters evolve, re-run
  at later visits and cross-tabulate subjects. Don't cluster `visit=all`
  for endotypes, because subjects co-cluster with themselves.
- Downstream: `stats-survival --mode=km --strata_file=.../assignments.csv`,
  `stats-group-testing --by=cluster` (after merging assignments), or
  `ml-classification` to build a classifier that reproduces the clusters.
