## =============================================================================
## lomics_utils.R
## Shared helpers for 4-visit multi-omics analysis (transcriptomics,
## proteomics, metabolomics). Sourced by every script in .claude/skills/ and by
## mcp/longomics_tools.R.
##
## Data contract
##   * One metadata CSV: one row per sample, with sample / subject / visit /
##     group columns (names set in study_config.json).
##   * One matrix CSV per omic layer: features in rows, samples in columns,
##     first column = feature IDs, column names = sample IDs.
## =============================================================================

lomics_read_config <- function(path = "study_config.json") {
  if (!file.exists(path)) {
    stop("Config file not found: ", path, call. = FALSE)
  }
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  needed <- c("metadata_file", "columns", "visits", "layers")
  miss <- setdiff(needed, names(cfg))
  if (length(miss)) {
    stop("study_config.json is missing: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  if (is.null(cfg$output_dir)) cfg$output_dir <- "results"
  if (!is.data.frame(cfg$layers)) cfg$layers <- as.data.frame(cfg$layers)
  cfg$time_values <- if (!is.null(cfg$time_values)) {
    unlist(cfg$time_values)
  } else {
    stats::setNames(seq_along(cfg$visits) - 1, cfg$visits)
  }
  cfg
}

lomics_out <- function(cfg, ...) {
  p <- file.path(cfg$output_dir, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

lomics_read_metadata <- function(cfg) {
  meta <- utils::read.csv(cfg$metadata_file, check.names = FALSE,
                          stringsAsFactors = FALSE)
  cn <- cfg$columns
  req <- unlist(cn[c("sample", "subject", "visit", "group")])
  miss <- setdiff(req, colnames(meta))
  if (length(miss)) {
    stop("Metadata is missing columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  ## Standardise column names to the LongOmicsTraj defaults
  meta$sample_id  <- as.character(meta[[cn$sample]])
  meta$subject_id <- as.character(meta[[cn$subject]])
  meta$visit      <- as.character(meta[[cn$visit]])
  meta$group      <- as.character(meta[[cn$group]])
  if (anyDuplicated(meta$sample_id)) {
    stop("Duplicated sample IDs in metadata.", call. = FALSE)
  }
  bad_visits <- setdiff(unique(meta$visit), cfg$visits)
  if (length(bad_visits)) {
    stop("Visit labels not listed in config$visits: ",
         paste(bad_visits, collapse = ", "), call. = FALSE)
  }
  meta$visit <- factor(meta$visit, levels = cfg$visits)
  meta$time  <- unname(cfg$time_values[as.character(meta$visit)])
  dup_sv <- duplicated(meta[, c("subject_id", "visit")])
  if (any(dup_sv)) {
    stop("Some subjects have more than one sample at the same visit ",
         "(technical replicates?). Collapse them first.", call. = FALSE)
  }
  rownames(meta) <- meta$sample_id
  meta
}

lomics_read_matrix <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (ext == "rds") {
    x <- readRDS(file)
    if (methods::is(x, "SummarizedExperiment")) {
      x <- SummarizedExperiment::assay(x)
    }
    x <- as.matrix(x)
  } else {
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    df <- utils::read.table(file, sep = sep, header = TRUE, row.names = 1,
                            check.names = FALSE, quote = "\"",
                            comment.char = "")
    x <- as.matrix(df)
  }
  storage.mode(x) <- "numeric"
  if (anyDuplicated(rownames(x))) {
    stop("Duplicated feature IDs in ", file, call. = FALSE)
  }
  x
}

## Align one layer to metadata; keeps only samples present in both.
lomics_align <- function(x, meta, layer = "layer") {
  common <- intersect(colnames(x), meta$sample_id)
  dropped_x <- setdiff(colnames(x), meta$sample_id)
  if (length(common) == 0) {
    stop(layer, ": no sample IDs match the metadata.", call. = FALSE)
  }
  if (length(dropped_x)) {
    message(layer, ": ", length(dropped_x),
            " matrix columns have no metadata and were dropped.")
  }
  m <- meta[common, , drop = FALSE]
  m <- m[order(m$subject_id, m$visit), , drop = FALSE]
  list(x = x[, m$sample_id, drop = FALSE], meta = m)
}

## Per-layer design audit: samples per visit x group, subjects with complete
## follow-up, missingness. Returns a list of data frames.
lomics_design_audit <- function(meta, x = NULL, visits) {
  tab <- as.data.frame.matrix(table(meta$group, factor(meta$visit, visits)))
  per_subj <- tapply(meta$visit, meta$subject_id,
                     function(v) sum(visits %in% as.character(v)))
  complete <- data.frame(
    n_subjects = length(per_subj),
    n_all_visits = sum(per_subj == length(visits)),
    n_one_visit_only = sum(per_subj == 1)
  )
  res <- list(samples_by_group_visit = tab, follow_up = complete)
  if (!is.null(x)) {
    res$missingness <- data.frame(
      features = nrow(x),
      samples = ncol(x),
      pct_missing_overall = round(100 * mean(is.na(x)), 2),
      features_gt20pct_missing = sum(rowMeans(is.na(x)) > 0.2),
      samples_gt20pct_missing = sum(colMeans(is.na(x)) > 0.2)
    )
  }
  res
}

lomics_build_se <- function(x, meta, assay_name = "values",
                            feature_type = "feature") {
  rd <- S4Vectors::DataFrame(feature_id = rownames(x),
                             feature_type = feature_type,
                             row.names = rownames(x))
  SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(list(x), assay_name),
    colData = S4Vectors::DataFrame(meta, row.names = meta$sample_id),
    rowData = rd
  )
}

## Load the processed (QC'd, log-scale) SummarizedExperiment for one layer
lomics_load_layer <- function(cfg, layer, stage = c("processed", "raw")) {
  stage <- match.arg(stage)
  f <- file.path(cfg$output_dir, "data", paste0(layer, "_", stage, "_se.rds"))
  if (!file.exists(f)) {
    stop("Missing ", f, ". Run the prep-study-setup",
         if (stage == "processed") " and prep-qc-normalisation", " step(s).",
         call. = FALSE)
  }
  readRDS(f)
}

lomics_layer_names <- function(cfg) cfg$layers$name

lomics_session_log <- function(cfg, tag) {
  f <- lomics_out(cfg, "logs", paste0(tag, "_sessionInfo.txt"))
  writeLines(utils::capture.output(utils::sessionInfo()), f)
  invisible(f)
}

## =============================================================================
## Cross-sectional / subject-level analysis helpers (used by the explore-,
## stats-, cluster- and ml- skills)
## =============================================================================

## Parse --key=value command-line arguments into a list, filling defaults.
##   a <- lomics_args(list(config = "study_config.json", visit = NULL))
## Values are character; convert where needed. "--flag" alone gives "TRUE".
lomics_args <- function(defaults = list(), args = commandArgs(trailingOnly = TRUE)) {
  out <- defaults
  for (a in args) {
    if (!startsWith(a, "--")) stop("Arguments must look like --key=value, got: ", a, call. = FALSE)
    kv <- strsplit(sub("^--", "", a), "=", fixed = TRUE)[[1]]
    key <- kv[1]; val <- if (length(kv) > 1) paste(kv[-1], collapse = "=") else "TRUE"
    out[[key]] <- val
  }
  out
}

## Split "a,b,c" -> c("a","b","c"); "all" -> all layer names
lomics_layers_arg <- function(cfg, x) {
  if (is.null(x) || identical(x, "all")) return(lomics_layer_names(cfg))
  strsplit(x, ",", fixed = TRUE)[[1]]
}

## Build an analysis matrix: rows = observations, columns = features.
##   layers       character vector of processed layers (joined side by side)
##   visit        a visit label (one row per subject), or "all" (one row per
##                sample; rows of the same subject are NOT independent)
##   change_from  optional reference visit: returns value(visit) - value(ref)
##                per subject (e.g. visit = "V3", change_from = "V1")
##   top_var      keep the N most variable features per layer (NULL = all)
## Returns list(X = matrix, meta = data.frame, feature_layer = named vector).
## Feature names are prefixed "<layer>__" when more than one layer is used.
lomics_xy <- function(cfg, layers, visit = cfg$default_visit, change_from = NULL,
                      top_var = NULL) {
  if (is.null(visit)) visit <- cfg$visits[1]
  blocks <- list(); metas <- list()
  for (layer in layers) {
    se <- lomics_load_layer(cfg, layer, "processed")
    x <- SummarizedExperiment::assay(se, "values")
    m <- as.data.frame(SummarizedExperiment::colData(se))
    m$visit <- as.character(m$visit)
    if (identical(visit, "all")) {
      if (!is.null(change_from)) stop("change_from cannot be used with visit = 'all'.", call. = FALSE)
      xm <- t(x); rownames(xm) <- m$sample_id; key <- m$sample_id
      mm <- m
    } else {
      sel <- m$visit == visit
      if (!any(sel)) stop(layer, ": no samples at visit ", visit, call. = FALSE)
      xv <- t(x[, sel, drop = FALSE]); rownames(xv) <- m$subject_id[sel]
      mm <- m[sel, , drop = FALSE]; rownames(mm) <- mm$subject_id
      if (!is.null(change_from)) {
        selr <- m$visit == change_from
        xr <- t(x[, selr, drop = FALSE]); rownames(xr) <- m$subject_id[selr]
        both <- intersect(rownames(xv), rownames(xr))
        xv <- xv[both, , drop = FALSE] - xr[both, , drop = FALSE]
        mm <- mm[both, , drop = FALSE]
      }
      xm <- xv; key <- rownames(xm)
    }
    if (!is.null(top_var) && ncol(xm) > top_var) {
      v <- apply(xm, 2, stats::var, na.rm = TRUE)
      xm <- xm[, order(v, decreasing = TRUE)[seq_len(top_var)], drop = FALSE]
    }
    if (length(layers) > 1) colnames(xm) <- paste(layer, colnames(xm), sep = "__")
    rownames(xm) <- key
    blocks[[layer]] <- xm; metas[[layer]] <- mm
  }
  keep <- Reduce(intersect, lapply(blocks, rownames))
  if (length(keep) < 5) stop("Fewer than 5 observations have all requested layers.", call. = FALSE)
  dropped <- setdiff(unique(unlist(lapply(blocks, rownames))), keep)
  if (length(dropped)) message(length(dropped), " observations lack at least one layer and were dropped.")
  X <- do.call(cbind, lapply(blocks, function(b) b[keep, , drop = FALSE]))
  meta <- metas[[1]][keep, , drop = FALSE]
  fl <- unlist(lapply(names(blocks), function(l) stats::setNames(rep(l, ncol(blocks[[l]])), colnames(blocks[[l]]))))
  list(X = X, meta = meta, feature_layer = fl,
       description = sprintf("layers=%s; visit=%s%s; n=%d; p=%d",
                             paste(layers, collapse = "+"), visit,
                             if (!is.null(change_from)) paste0(" minus ", change_from) else "",
                             nrow(X), ncol(X)))
}

## Resolve an outcome name: "binary"/"continuous"/... keys from config, or a
## literal metadata column name.
lomics_outcome <- function(cfg, x) {
  if (!is.null(cfg$outcomes) && x %in% names(cfg$outcomes)) return(cfg$outcomes[[x]][1])
  x
}

## Output folder for an analysis, tagged with layers/visit so runs don't clash
lomics_tag <- function(layers, visit, change_from = NULL, extra = NULL) {
  paste(c(paste(layers, collapse = "+"), visit,
          if (!is.null(change_from)) paste0("minus", change_from), extra), collapse = "_")
}

## Covariates from the config that exist and vary in `meta`
lomics_covariates <- function(cfg, meta, extra = NULL) {
  cv <- unique(c(cfg$covariates, extra))
  cv <- intersect(cv, colnames(meta))
  cv[vapply(cv, function(c) length(unique(stats::na.omit(meta[[c]]))) > 1, logical(1))]
}

lomics_strip_layer <- function(x) sub("^.*?__", "", x)
