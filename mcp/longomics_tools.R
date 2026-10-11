## =============================================================================
## longomics_tools.R — custom MCP tools for the 4-visit multi-omics project.
##
## Served by mcptools:   Rscript -e "mcptools::mcp_server(tools = 'mcp/longomics_tools.R')"
## (see .mcp.json). mcptools sources this file and uses the list returned by
## its LAST expression as the tool set.
##
## The plain R functions (lomics_tool_*) can also be called directly in RStudio,
## which is how they were tested.
## =============================================================================

.lomics_root <- function() normalizePath(Sys.getenv("LOMICS_PROJECT", getwd()), mustWork = FALSE)

.lomics_trunc <- function(txt, max_chars = 20000) {
  txt <- paste(txt, collapse = "\n")
  if (nchar(txt) > max_chars) paste0(substr(txt, 1, max_chars), "\n... [truncated]") else txt
}

.lomics_cfg <- function() {
  old <- setwd(.lomics_root()); on.exit(setwd(old))
  source("R/lomics_utils.R", local = TRUE)
  cfg <- lomics_read_config("study_config.json")
  cfg$output_dir <- file.path(.lomics_root(), cfg$output_dir)
  cfg
}

## 1. Evaluate R code and return printed output, messages, warnings, errors.
lomics_tool_run_r <- function(code) {
  out <- character(); notes <- character()
  res <- withCallingHandlers(
    tryCatch({
      exprs <- parse(text = code)
      out <- utils::capture.output({
        for (e in exprs) {                       # print each visible result, like the console
          val <- withVisible(eval(e, envir = globalenv()))
          if (val$visible) print(val$value)
        }
      })
      "ok"
    }, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(w) { notes <<- c(notes, paste("Warning:", conditionMessage(w))); invokeRestart("muffleWarning") },
    message = function(m) { notes <<- c(notes, sub("\n$", "", conditionMessage(m))); invokeRestart("muffleMessage") })
  .lomics_trunc(c(out, notes, if (res != "ok") res))
}

## 2. Which pipeline steps have outputs on disk.
lomics_tool_status <- function() {
  cfg <- .lomics_cfg(); od <- cfg$output_dir
  layers <- cfg$layers$name
  derived <- sub("_processed_se\\.rds$", "", list.files(file.path(od, "data"), pattern = "^pathways_.*_processed_se\\.rds$"))
  chk <- function(...) file.exists(file.path(od, ...))
  rows <- lapply(c(layers, derived), function(l) data.frame(
    layer = l,
    raw = chk("data", paste0(l, "_raw_se.rds")),
    processed = chk("data", paste0(l, "_processed_se.rds")),
    de = dir.exists(file.path(od, "de", l)),
    trajectories = chk("trajectories", l, "consensus.csv")))
  st <- do.call(rbind, rows)
  runs <- unlist(lapply(c("dimred", "tests", "regression", "survival", "clustering",
                          "feature_selection", "classification", "integration", "pathways"), function(cat) {
    dd <- list.dirs(file.path(od, cat), recursive = FALSE, full.names = FALSE)
    if (length(dd)) sprintf("%s: %s", cat, paste(dd, collapse = ", "))
  }))
  extra <- c("", "Completed analysis runs (results/<category>/<run>):", if (length(runs)) runs else "none yet")
  .lomics_trunc(c(sprintf("Project: %s | visits: %s", .lomics_root(), paste(cfg$visits, collapse = " -> ")),
                  utils::capture.output(print(st, row.names = FALSE)), extra))
}

## 3. Skills and their scripts, discovered from .claude/skills/*/scripts/*.R
.lomics_skill_scripts <- function() {
  f <- Sys.glob(file.path(.lomics_root(), ".claude", "skills", "*", "scripts", "*.R"))
  data.frame(skill = basename(dirname(dirname(f))), script = basename(f),
             path = sub(paste0("^", .lomics_root(), "/"), "", f), stringsAsFactors = FALSE)
}

## Split an argument string like a shell would: respects "double" or 'single'
## quotes, so --formula="fev1 ~ age + sex" stays one argument.
.lomics_split_args <- function(x) {
  if (!nzchar(trimws(x))) return(character())
  tok <- regmatches(x, gregexpr("(?:[^\\s\"']+|\"[^\"]*\"|'[^']*')+", x, perl = TRUE))[[1]]
  gsub("[\"']", "", tok)
}

lomics_tool_list_skills <- function() {
  sc <- .lomics_skill_scripts()
  desc <- vapply(unique(sc$skill), function(sk) {
    md <- readLines(file.path(.lomics_root(), ".claude", "skills", sk, "SKILL.md"), warn = FALSE)
    d <- sub("^description:\\s*", "", grep("^description:", md, value = TRUE)[1])
    if (nchar(d) > 160) paste0(substr(d, 1, 157), "...") else d
  }, "")
  .lomics_trunc(c("Skills with runnable scripts (read the skill's SKILL.md for arguments):",
                  sprintf("- %s [%s]: %s", sc$skill, sc$script, desc[sc$skill])))
}

## Run a skill script as a separate Rscript process from the project root
## (keeps RStudio free and reproduces exactly what a batch run would do).
lomics_tool_run_script <- function(skill, args = "", script = "") {
  sc <- .lomics_skill_scripts(); sc <- sc[sc$skill == skill, ]
  if (!nrow(sc)) return(paste("Unknown skill. Available:", paste(unique(.lomics_skill_scripts()$skill), collapse = ", ")))
  if (nzchar(script)) sc <- sc[sc$script == script, ]
  if (nrow(sc) != 1) return(paste("Skill", skill, "has several scripts; pass script = one of:",
                                  paste(.lomics_skill_scripts()$script[.lomics_skill_scripts()$skill == skill], collapse = ", ")))
  old <- setwd(.lomics_root()); on.exit(setwd(old))
  a <- .lomics_split_args(args)
  rscript <- file.path(R.home("bin"), "Rscript")
  out <- suppressWarnings(system2(rscript, shQuote(c(sc$path, a)), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status"); if (is.null(status)) status <- 0
  .lomics_trunc(c(sprintf("[%s/%s] exit status %s", skill, sc$script, status),
                  sprintf("args: %s", paste(a, collapse = " | ")), utils::tail(out, 150)))
}

## 4. Query trajectory signatures. `topology` accepts wildcards: "up_*_*",
##    "*_down", "flat_flat_flat". Uses significance-gated labels if gated = TRUE.
lomics_tool_query_ots <- function(layer = "", group = "", topology = "", gated = TRUE, max_rows = 50) {
  cfg <- .lomics_cfg()
  f <- file.path(cfg$output_dir, "trajectories", "ots_all_layers.csv")
  if (!file.exists(f)) return("No trajectories yet — run the 'trajectories' step.")
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  col <- if (isTRUE(gated) && "topology_gated" %in% names(d) && any(!is.na(d$topology_gated))) "topology_gated" else "topology"
  if (nzchar(layer)) d <- d[d$layer == layer, ]
  if (nzchar(group)) d <- d[d$group == group, ]
  if (nzchar(topology)) {
    rx <- paste0("^", gsub("\\*", "[a-z]+", topology), "$")
    d <- d[!is.na(d[[col]]) & grepl(rx, d[[col]]), ]
  }
  .lomics_trunc(c(sprintf("%d matches (label column: %s)", nrow(d), col),
                  utils::capture.output(print(utils::head(d, max_rows), row.names = FALSE)),
                  "", "Counts by layer x group:",
                  utils::capture.output(print(table(d$layer, d$group)))))
}

## 5. Everything known about one feature: DE statistics across contrasts + OTS.
lomics_tool_feature <- function(feature, layer) {
  cfg <- .lomics_cfg()
  de_dir <- file.path(cfg$output_dir, "de", layer)
  rows <- list()
  for (f in list.files(de_dir, pattern = "\\.csv$", full.names = TRUE)) {
    d <- utils::read.csv(f); r <- d[d$feature == feature, , drop = FALSE]
    if (nrow(r)) rows[[basename(f)]] <- data.frame(
      contrast = sub("\\.csv$", "", basename(f)),
      logFC = if ("logFC" %in% names(r)) signif(r$logFC, 3) else NA,
      P = signif(r$P.Value, 3), FDR = signif(r$adj.P.Val, 3))
  }
  ots_f <- file.path(cfg$output_dir, "trajectories", "ots_all_layers.csv")
  ots <- if (file.exists(ots_f)) { o <- utils::read.csv(ots_f); o[o$layer == layer & o$feature == feature, ] } else NULL
  .lomics_trunc(c(sprintf("Feature %s in %s", feature, layer),
                  if (length(rows)) utils::capture.output(print(do.call(rbind, rows), row.names = FALSE)) else "No DE results.",
                  "", "Trajectory signatures:",
                  if (!is.null(ots) && nrow(ots)) utils::capture.output(print(ots, row.names = FALSE)) else "None."))
}

## 6. Read a results text/CSV file (paths relative to the project root).
lomics_tool_read_result <- function(path, max_rows = 40) {
  p <- file.path(.lomics_root(), path)
  if (!file.exists(p)) return(paste("Not found:", path))
  if (grepl("\\.csv$", p)) {
    d <- utils::read.csv(p, check.names = FALSE)
    return(.lomics_trunc(c(sprintf("%d rows x %d cols", nrow(d), ncol(d)),
                           utils::capture.output(print(utils::head(d, max_rows), row.names = FALSE)))))
  }
  .lomics_trunc(readLines(p, warn = FALSE))
}

## ---------------------------------------------------------------------------
## Tool registry (last expression = value returned to mcptools)
## ---------------------------------------------------------------------------
if (requireNamespace("ellmer", quietly = TRUE)) {
  list(
    ellmer::tool(lomics_tool_run_r, name = "run_r",
      description = "Evaluate R code in the global environment of the R process and return printed output, messages, warnings and errors. Use for inspecting objects, quick checks and small computations; use run_skill_script for full analyses.",
      arguments = list(code = ellmer::type_string("R code to evaluate."))),
    ellmer::tool(lomics_tool_status, name = "project_status",
      description = "Show the study design from study_config.json and which preparation steps exist for each omic layer, and which analysis runs (dimred, tests, regression, survival, clustering, feature selection, classification, pathways, integration) have results.",
      arguments = list()),
    ellmer::tool(lomics_tool_list_skills, name = "list_skills",
      description = "List every analysis skill that has a runnable R script (data preparation, dimensionality reduction, statistical testing, regression, longitudinal, survival, clustering, feature selection, classification, pathways, integration) with a one-line description.",
      arguments = list()),
    ellmer::tool(lomics_tool_run_script, name = "run_skill_script",
      description = "Run a skill's R script as a separate Rscript process from the project root and return its exit status and last 150 output lines. args is the argument string exactly as on the command line, e.g. '--outcome=binary --layers=proteomics --visit=V1' or '--mode=formula --formula=\"fev1 ~ age + PROT21\"'. Read the skill's SKILL.md first for the arguments.",
      arguments = list(
        skill = ellmer::type_string("Skill folder name, e.g. ml-classification."),
        args = ellmer::type_string("Command-line arguments.", required = FALSE),
        script = ellmer::type_string("Script file name, only if the skill has several.", required = FALSE))),
    ellmer::tool(lomics_tool_query_ots, name = "query_trajectories",
      description = "Query Omics Trajectory Signatures (LongOmicsTraj consensus) across layers. Filter by layer, group (treatment arm) and topology pattern with * wildcards, e.g. 'up_*_down'. gated = TRUE uses labels where features without a significant within-group time effect are 'no_time_effect'.",
      arguments = list(
        layer = ellmer::type_string("Omic layer name, or empty for all.", required = FALSE),
        group = ellmer::type_string("Group / arm, or empty for all.", required = FALSE),
        topology = ellmer::type_string("Topology pattern, e.g. 'flat_up_up' or 'up_*_*'.", required = FALSE),
        gated = ellmer::type_boolean("Use significance-gated labels.", required = FALSE),
        max_rows = ellmer::type_integer("Maximum rows to print.", required = FALSE))),
    ellmer::tool(lomics_tool_feature, name = "feature_report",
      description = "Report every longitudinal differential statistic (all contrasts) and trajectory signature for one feature in one layer.",
      arguments = list(
        feature = ellmer::type_string("Feature ID (gene, protein, metabolite or pathway)."),
        layer = ellmer::type_string("Omic layer name."))),
    ellmer::tool(lomics_tool_read_result, name = "read_result",
      description = "Read a result file (CSV preview or text) by path relative to the project root, e.g. 'results/de/de_summary.csv' or 'results/qc/proteomics_qc_summary.txt'.",
      arguments = list(
        path = ellmer::type_string("Path relative to the project root."),
        max_rows = ellmer::type_integer("Rows to show for CSV files.", required = FALSE)))
  )
} else {
  warning("ellmer is not installed; returning no MCP tools. install.packages('ellmer')")
  list()
}
