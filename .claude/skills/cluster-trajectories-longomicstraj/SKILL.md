---
name: cluster-trajectories-longomicstraj
description: Classify the longitudinal shape of every gene, protein, metabolite or pathway across visits into Omics Trajectory Signatures (e.g. up_flat_down, flat_up_up) with the LongOmicsTraj R package (nzkermani/LongOmicsTraj) — empirical, LMM, GAM and maSigPro estimators, per-arm consensus, cross-omic topology comparison, trajectory queries and cluster-level topology. Use when the user asks about trajectory shapes, rebound / transient / sustained responses, OTS, topology, or which features share a temporal pattern across omic layers or arms.
---

# Trajectory signatures with LongOmicsTraj

## Run
```bash
Rscript .claude/skills/cluster-trajectories-longomicstraj/scripts/run_longomicstraj.R \
    study_config.json [layer|all] [empirical,lmm,gam,masigpro] [hybrid] [0.15] [0.5]
```
Arguments: layer, estimators, threshold method (`fixed`, `global_mad`,
`gene_mad`, `hybrid`, `quantile`), `delta_min` (log2 units), `k` (MAD multiplier).
Run **stats-longitudinal-de first** so results can be significance-gated.

Outputs (`results/trajectories/`): per layer `lot_<estimator>.rds` (full
LongOmicsTraj objects), `ots_all_estimators.csv`, `consensus.csv`,
`summary_<estimator>.txt`, `ots_distribution.pdf`; across layers
`ots_all_layers.csv` (with `time_fdr`, `topology_gated`),
`ots_layer_group_proportions.csv`, `discordant_between_groups.csv`,
`ots_layers_comparison.pdf`.

## Four visits → three transitions
Labels have three parts (`state1_state2_state3`, each up/flat/down): 27
possible signatures. Examples: `up_flat_flat` early sustained response,
`up_flat_down` response then return, `flat_flat_up` late response,
`up_down_flat` transient, `down_up_flat` rebound.
`lot_query(object, topology = "...")` matches exact labels; for wildcard
queries across layers use MCP `query_trajectories` (`"up_*_down"`).

## Package behaviours to respect (LongOmicsTraj 0.99.0, checked in source)
- **One object per omic layer.** `lot_from_*()` errors if a VMM exists and
  `overwrite = TRUE` replaces the VMM for all assays, so assays cannot be
  accumulated in one object. The script builds one object per layer and
  merges OTS tables afterwards. Do not "fix" this by adding experiments to
  one MAE.
- **`lot_consensus()` pools arms.** It votes across estimator × group
  columns, giving one label per feature. The script calls it within each
  group so arm-specific trajectories are kept. If writing new code, do the same.
- Ties in consensus are broken by `which.max` (first label alphabetically).
- `lot_plot_ots_distribution()` uses Brewer `Set3` (max 12 colours); with up
  to 27 labels it warns and recycles colours — harmless but legends repeat.
- `lot_create(feature_type = ...)` — OTS output still reports
  `object_type = "gene"` for all layers; use the `layer` column instead.
- Package declares `R (>= 4.6.0)`.

## Significance gating — important for interpretation
Thresholds define how big a visit-to-visit change must be to count as
up/down, but they are not a test. On noisy group means, null features
receive non-flat labels (in synthetic null data most did). `topology_gated`
replaces the label with `no_time_effect` when the arm's overall-time FDR from
the limma F-test is ≥ 0.05. Report gated labels by default and say so;
show ungated distributions only as descriptive.

## Estimator choice
- `empirical`: raw visit means per arm — transparent, noisiest.
- `lmm`: `y ~ visit * group + (1 | subject_id)` — handles missing visits;
  recommended primary estimator.
- `gam`: smooth in real `time` per arm, `k = min(3, n_visits - 1)` — use when
  visit spacing is unequal; smooths sharp transients.
- `masigpro`: polynomial regression (needs `maSigPro`).
Consensus = majority (≥ 2 estimators agree) within each arm. Report
`agreement_count` / `n_estimators`.

## Further analyses available in the package (call from R directly)
- `lot_from_clusters(obj, assay, visits, k, estimator = "kmeans"|"pam")`,
  `lot_plot_cluster_trajectory_arms()`, `lot_plot_topology_deviation()` —
  topology composition of profile clusters.
- `lot_from_flexmix_clusters(obj, assay, group, visits, k = 2:5, model = "quadratic")`,
  `lot_plot_cluster_selection()`, `lot_plot_cluster_topology()`.
- `lot_query(obj, conditions = list(Placebo = "flat_flat_flat", Treatment = NOT("flat_flat_flat")))`
  — cross-arm queries.
- `lot_plot_trajectory(obj, features = c(...), assay = layer, show_subjects = TRUE)`.
Full signatures: `references/api.md`. Check `?function` in R before using an
argument not listed there.

## Cross-omic questions
`ots_layer_group_proportions.csv` compares how often each topology occurs per
layer and arm. Matching features across layers (gene ↔ protein) needs an ID
map from the user; join on mapped IDs and cross-tabulate topologies. Do not
claim layer concordance from name matching alone.
