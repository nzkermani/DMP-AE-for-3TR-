# LongOmicsTraj API reference (v0.99.0, from man/*.Rd usage sections)

Extracted from the package source. `\S4method{fn}{LongOmicsTraj}(...)` lines give the full argument list for S4 generics. Check `?fn` in R if your installed version differs.

```r
compute_vmm_deltas(vmm, visits)
lot_mae(object)

lot_vmm(object)

lot_deltas(object)

lot_thresholds(object)

lot_ots(object)

lot_relationships(object)

lot_results(object)

\S4method{lot_mae}{LongOmicsTraj}(object)

\S4method{lot_vmm}{LongOmicsTraj}(object)

\S4method{lot_deltas}{LongOmicsTraj}(object)

\S4method{lot_thresholds}{LongOmicsTraj}(object)

\S4method{lot_ots}{LongOmicsTraj}(object)

\S4method{lot_relationships}{LongOmicsTraj}(object)

\S4method{lot_results}{LongOmicsTraj}(object)
lot_cluster_purity(assignments, purity_threshold = 0.9)
lot_cluster_topology_score(
  purity,
  k = NULL,
  bic = NA_real_,
  aic = NA_real_,
  minimum_cluster_size = NA_integer_
)
lot_compare_estimators(object, ..., .names = NULL)

\S4method{lot_compare_estimators}{LongOmicsTraj}(object, ..., .names = NULL)
lot_compute_deltas(object, ...)

\S4method{lot_compute_deltas}{LongOmicsTraj}(object, assay, visits, overwrite = FALSE)
lot_compute_thresholds(object, ...)

\S4method{lot_compute_thresholds}{LongOmicsTraj}(
  object,
  assay,
  method = c("fixed", "global_mad", "gene_mad", "hybrid", "quantile"),
  delta_min = NULL,
  k = 0.25,
  prob = 0.75,
  visits,
  overwrite = FALSE,
  verbose = TRUE
)
lot_consensus(x, min_agreement = 2L, mode = c("strict", "majority"))

\S4method{lot_consensus}{DFrame}(x, min_agreement = 2L, mode = c("strict", "majority"))
lot_create(
  expression,
  metadata,
  assay = "transcriptomics",
  sample_col = "sample_id",
  subject_col = "subject_id",
  visit_col = "visit",
  group_col = "group",
  visit_levels = NULL,
  time_values = NULL,
  time_col = "time",
  feature_type = "feature"
)
lot_from_clusters(object, ...)

\S4method{lot_from_clusters}{LongOmicsTraj}(
  object,
  assay,
  group_col = "group",
  visit_col = "visit",
  visits,
  k = 4L,
  estimator = c("kmeans", "pam"),
  assay_data = NULL,
  scale_features = TRUE,
  standardise_ids = FALSE,
  seed = 1L,
  require_gene_ots = FALSE,
  overwrite = FALSE,
  verbose = TRUE
)
lot_from_flexmix_clusters(
  object,
  assay,
  group,
  visits,
  k = 2:5,
  assay_data = NULL,
  group_col = "group",
  visit_col = "visit",
  sample_col = "sample_id",
  trajectory_source = c("estimated", "mean"),
  trajectory_estimator = NULL,
  model = c("quadratic", "linear"),
  scale_features = TRUE,
  purity_threshold = 0.9,
  required_passing_fraction = 0.8,
  min_cluster_size = 5L,
  nrep = 10L,
  seed = 42L,
  store_fits = FALSE,
  overwrite = FALSE,
  verbose = TRUE
)
lot_map_to_ots(object, visits, overwrite = FALSE, ...)

\S4method{lot_map_to_ots}{LongOmicsTraj}(object, visits, overwrite = FALSE, ...)
lot_plot_cluster_selection(
  object,
  metrics = c("BIC", "weighted_purity", "passing_fraction", "minimum_cluster_size"),
  purity_threshold = 0.9,
  required_passing_fraction = 0.8,
  min_cluster_size = 5L,
  show_selected = TRUE
)
lot_plot_cluster_topology(
  object,
  clusters = NULL,
  group = NULL,
  value_type = c("scaled", "raw"),
  show_features = TRUE,
  max_feature_lines = 50L,
  purity_threshold = 0.9,
  sort_by = c("cluster", "purity", "size", "dominant_topology"),
  seed = 1L
)
lot_plot_cluster_trajectory_arms(
  object,
  clusters = NULL,
  groups = NULL,
  top_n = NULL,
  min_proportion = 0,
  show_counts = TRUE,
  show_percent = TRUE,
  scale_width = TRUE,
  ncol = NULL
)
lot_plot_ots_distribution(
  object,
  estimator = NULL,
  type = c("proportion", "count"),
  palette = "Set3"
)
lot_plot_topology_deviation(
  object,
  clusters = NULL,
  group = NULL,
  min_abs_deviation = 0,
  top_n = NULL,
  sort_topologies = TRUE
)
lot_plot_trajectory(
  object,
  features,
  assay = "transcriptomics",
  assay_data = NULL,
  groups = NULL,
  estimator = NULL,
  show_subjects = FALSE,
  summary = c("mean_se", "mean_ci", "mean"),
  facet_scales = c("free_y", "fixed"),
  ncol = NULL
)
lot_query(object, ...)

NOT(x)

\S4method{lot_query}{LongOmicsTraj}(
  object,
  object_type = NULL,
  assay = NULL,
  group = NULL,
  topology = NULL,
  conditions = NULL,
  estimator = NULL
)
lot_summary(object, estimator = NULL, top_n = 5L)
lot_from_empirical(object, ...)

lot_from_lmm(object, ...)

lot_from_gam(object, ...)

lot_from_masigpro(object, ...)

\S4method{lot_from_empirical}{LongOmicsTraj}(
  object,
  assay,
  visits,
  group_col = "group",
  visit_col = "visit",
  estimator = "empirical",
  overwrite = FALSE,
  verbose = TRUE
)

\S4method{lot_from_gam}{LongOmicsTraj}(
  object,
  assay,
  formula_fixed,
  group_col = "group",
  visit_col = "visit",
  subject_col = "subject_id",
  visits,
  estimator = "gam",
  overwrite = FALSE,
  verbose = TRUE
)

\S4method{lot_from_lmm}{LongOmicsTraj}(
  object,
  assay,
  formula_fixed,
  subject_col = "subject_id",
  group_col = "group",
  visit_col = "visit",
  visits,
  estimator = "lmm",
  BPPARAM = NULL,
  overwrite = FALSE,
  verbose = TRUE
)

\S4method{lot_from_masigpro}{LongOmicsTraj}(
  object,
  assay,
  group_col = "group",
  visit_col = "visit",
  time_col = "time",
  subject_col = "subject_id",
  visits,
  degree = NULL,
  estimator = "masigpro",
  overwrite = FALSE,
  verbose = TRUE
)
plot_cluster_spaghetti(object, assay_name, visits_vector, max_clusters = 3L)
plot_cluster_spaghetti_single(object, cluster_id, assay_name, visits_vector)
\S4method{show}{LongOmicsTraj}(object)
```
