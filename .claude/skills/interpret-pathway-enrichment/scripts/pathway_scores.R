#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## pathway_scores.R — sample-level pathway scores (ssGSEA or GSVA) for one
## layer, saved as a NEW processed layer so the longitudinal-differential and
## longomicstraj-trajectories scripts can be run on pathways exactly as on
## features. Optionally also runs preranked fgsea on every DE contrast.
##
## Usage (from the project root):
##   Rscript .claude/skills/interpret-pathway-enrichment/scripts/pathway_scores.R \
##       <config> <layer> <gene_sets> [method] [min_size] [max_size]
##   gene_sets  path to a .gmt file  OR  "msigdbr:<collection>" e.g. msigdbr:H
##              (for metabolomics supply a .gmt whose members match your
##               metabolite IDs, e.g. built from KEGG / SMPDB / Reactome)
##   method     ssgsea (default) | gsva
##
## Feature IDs in the layer must match the gene-set member IDs (e.g. HGNC
## symbols). Map Ensembl IDs / UniProt accessions / Olink IDs first.
##
## Outputs:
##   results/data/pathways_<layer>_processed_se.rds   (assay "values")
##   results/pathways/<layer>/coverage.csv            set size before/after matching
##   results/pathways/<layer>/fgsea_<contrast>.csv    if fgsea + DE results exist
## Then run, e.g.:
##   Rscript .claude/skills/stats-longitudinal-de/scripts/longitudinal_de.R study_config.json pathways_<layer>
## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lomics_utils.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: pathway_scores.R <config> <layer> <gene_sets> [method] [min] [max]")
cfg      <- lomics_read_config(args[1])
layer    <- args[2]
gs_spec  <- args[3]
method   <- if (length(args) >= 4) args[4] else "ssgsea"
min_size <- if (length(args) >= 5) as.integer(args[5]) else 10L
max_size <- if (length(args) >= 6) as.integer(args[6]) else 500L

read_gmt <- function(f) {
  lines <- strsplit(readLines(f, warn = FALSE), "\t")
  stats::setNames(lapply(lines, function(l) unique(l[-(1:2)])), vapply(lines, `[`, "", 1))
}
gene_sets <- if (startsWith(gs_spec, "msigdbr:")) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) stop("Install msigdbr or supply a .gmt file.")
  coll <- sub("^msigdbr:", "", gs_spec)
  ## msigdbr >= 10 uses `collection`; older versions use `category`
  fa <- names(formals(msigdbr::msigdbr))
  m <- if ("collection" %in% fa) msigdbr::msigdbr(species = "Homo sapiens", collection = coll)
       else msigdbr::msigdbr(species = "Homo sapiens", category = coll)
  split(m$gene_symbol, m$gs_name)
} else read_gmt(gs_spec)

se <- lomics_load_layer(cfg, layer, "processed")
x  <- assay(se, "values")
matched <- lapply(gene_sets, intersect, rownames(x))
cov <- data.frame(set = names(gene_sets), size = lengths(gene_sets), matched = lengths(matched))
outdir <- dirname(lomics_out(cfg, "pathways", layer, "x"))
utils::write.csv(cov, file.path(outdir, "coverage.csv"), row.names = FALSE)
keep <- cov$matched >= min_size & cov$matched <= max_size
message(sprintf("%s: %d / %d sets with %d-%d matched members (median matched fraction %.2f)",
                layer, sum(keep), length(keep), min_size, max_size,
                stats::median(cov$matched / pmax(cov$size, 1))))
if (sum(keep) == 0) stop("No gene sets pass the size filter — check that feature IDs match the gene-set IDs.")
sets <- matched[keep]

if (exists("ssgseaParam", asNamespace("GSVA"))) {            # GSVA >= 1.50
  par <- if (method == "ssgsea") GSVA::ssgseaParam(x, sets, minSize = min_size, maxSize = max_size)
         else GSVA::gsvaParam(x, sets, minSize = min_size, maxSize = max_size)
  scores <- GSVA::gsva(par, verbose = FALSE)
} else {                                                       # older GSVA
  scores <- GSVA::gsva(x, sets, method = method, min.sz = min_size, max.sz = max_size, verbose = FALSE)
}
## ssGSEA scores are on an arbitrary scale; z-score each pathway across
## samples so LongOmicsTraj thresholds (delta_min) are comparable to log2 data.
if (method == "ssgsea") scores <- t(scale(t(scores)))

out <- SummarizedExperiment(
  assays = list(values = scores),
  colData = colData(se),
  rowData = S4Vectors::DataFrame(feature_id = rownames(scores), feature_type = "pathway",
                                 n_matched = lengths(sets)[rownames(scores)],
                                 row.names = rownames(scores)))
S4Vectors::metadata(out)$source <- list(layer = layer, gene_sets = gs_spec, method = method)
saveRDS(out, lomics_out(cfg, "data", paste0("pathways_", layer, "_processed_se.rds")))
message("Saved pathway layer: pathways_", layer, " (", nrow(scores), " pathways)")

## Preranked GSEA on every DE contrast of the source layer (optional)
de_dir <- file.path(cfg$output_dir, "de", layer)
if (dir.exists(de_dir) && requireNamespace("fgsea", quietly = TRUE)) {
  for (f in list.files(de_dir, pattern = "^(change|step|interaction)_.*\\.csv$", full.names = TRUE)) {
    de <- utils::read.csv(f)
    stat <- stats::setNames(de$t, de$feature)
    stat <- stat[is.finite(stat)]
    set.seed(42)
    res <- fgsea::fgsea(pathways = sets, stats = stat, minSize = min_size, maxSize = max_size)
    res$leadingEdge <- vapply(res$leadingEdge, paste, "", collapse = ";")
    utils::write.csv(res[order(res$padj), ], file.path(outdir, paste0("fgsea_", basename(f))), row.names = FALSE)
  }
} else if (dir.exists(de_dir)) {
  message("fgsea not installed — skipping preranked GSEA on DE contrasts.")
}
lomics_session_log(cfg, paste0("pathways_", layer))
