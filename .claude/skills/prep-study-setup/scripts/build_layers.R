#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## build_layers.R — read study_config.json, align every omic layer to the
## metadata, audit the 4-visit design, and save raw SummarizedExperiments.
##
## Usage (from the project root):
##   Rscript .claude/skills/prep-study-setup/scripts/build_layers.R [config]
## Outputs:
##   results/data/<layer>_raw_se.rds
##   results/data/raw_mae.rds              (MultiAssayExperiment of all layers)
##   results/setup/design_audit.txt
##   results/setup/sample_overlap.csv      (which layers each sample has)
## -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
})
source("R/lomics_utils.R")

args <- commandArgs(trailingOnly = TRUE)
cfg  <- lomics_read_config(if (length(args)) args[1] else "study_config.json")
meta <- lomics_read_metadata(cfg)

audit_lines <- c(sprintf("Study: %s", cfg$study_name),
                 sprintf("Visits (ordered): %s", paste(cfg$visits, collapse = " -> ")),
                 sprintf("Metadata samples: %d, subjects: %d, groups: %s",
                         nrow(meta), length(unique(meta$subject_id)),
                         paste(unique(meta$group), collapse = ", ")),
                 "")
se_list <- list()
for (i in seq_len(nrow(cfg$layers))) {
  L <- cfg$layers[i, ]
  x <- lomics_read_matrix(L$file)
  al <- lomics_align(x, meta, L$name)
  se <- lomics_build_se(al$x, al$meta, assay_name = "raw",
                        feature_type = L$feature_type)
  S4Vectors::metadata(se)$scale <- L$scale
  saveRDS(se, lomics_out(cfg, "data", paste0(L$name, "_raw_se.rds")))
  se_list[[L$name]] <- se

  au <- lomics_design_audit(al$meta, al$x, cfg$visits)
  audit_lines <- c(audit_lines,
                   sprintf("== %s (%s scale) ==", L$name, L$scale),
                   utils::capture.output(print(au$samples_by_group_visit)),
                   utils::capture.output(print(au$follow_up, row.names = FALSE)),
                   utils::capture.output(print(au$missingness, row.names = FALSE)),
                   "")
  if (L$scale == "counts" && any(al$x < 0, na.rm = TRUE)) {
    warning(L$name, " is declared as counts but has negative values.")
  }
  if (L$scale == "log2" && max(al$x, na.rm = TRUE) > 50) {
    warning(L$name, " is declared log2 but max value > 50; is it really logged?")
  }
}

## Sample overlap across layers (needed for integration)
all_ids <- meta$sample_id
overlap <- data.frame(sample_id = all_ids,
                      subject_id = meta$subject_id,
                      visit = meta$visit,
                      group = meta$group)
for (nm in names(se_list)) overlap[[nm]] <- all_ids %in% colnames(se_list[[nm]])
overlap$n_layers <- rowSums(overlap[, names(se_list), drop = FALSE])
utils::write.csv(overlap, lomics_out(cfg, "setup", "sample_overlap.csv"),
                 row.names = FALSE)
audit_lines <- c(audit_lines, "== Samples with all layers, by visit ==",
                 utils::capture.output(print(table(
                   overlap$visit[overlap$n_layers == length(se_list)]))))

## MultiAssayExperiment of raw layers (colData = subject-visit level metadata)
mae <- MultiAssayExperiment(
  experiments = ExperimentList(lapply(se_list, function(se) {
    SummarizedExperiment(assays = assays(se), rowData = rowData(se))
  })),
  colData = S4Vectors::DataFrame(meta, row.names = meta$sample_id)
)
saveRDS(mae, lomics_out(cfg, "data", "raw_mae.rds"))

writeLines(audit_lines, lomics_out(cfg, "setup", "design_audit.txt"))
cat(audit_lines, sep = "\n")
lomics_session_log(cfg, "setup")
