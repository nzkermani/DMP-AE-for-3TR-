#!/usr/bin/env Rscript
## -----------------------------------------------------------------------------
## multilevel_splsda.R — paired (multilevel) sparse PLS-DA between two visits,
## keepX tuning, backward feature elimination to a small panel, and validation
## at three levels of honesty.
##
## Reproduces the analysis in Vassileva et al., Respir Med 2026 (Methods 2.5)
## and adds the checks that design needs:
##   A  training AUROC of the final panel (paper: "training performance")
##   B  leave-one-subject-out CV of the FIXED panel (paper's LOSO-CV). The
##      panel was selected on all subjects, so B is still optimistic.
##   C  nested CV: tuning + elimination repeated inside every outer fold
##      (subject-grouped K-fold x repeats, or LOSO; skip with --nested=FALSE)
##   D  within-subject label-flip permutation null of A and B (--n_perm)
##
## Usage (from the project root):
##   Rscript .claude/skills/ml-multilevel-splsda/scripts/multilevel_splsda.R \
##     --config=study_config.json --layer=breath_voc --visits=BL,M4 \
##     [--ncomp=2] [--grid=5,10,15,20,30,50,75,100] [--zscore=TRUE] \
##     [--drop_max=0.05] [--drop_type=relative|absolute] [--target_auc=0.9|NA] \
##     [--elim_auc=train|loso] [--nested=TRUE] [--nested_grid=5,10,20,50] \
##     [--nested_folds=5|loso] [--nested_repeats=5] \
##     [--n_perm=20] [--panel=f1,f2,...] [--plot_visits=BL,M1,M4] \
##     [--colour_by=group] [--engine=auto|mixOmics|internal] [--seed=]
##   --panel skips selection and evaluates a given feature set.
##
## Outputs: <output_dir>/splsda/<v1>vs<v2>_<layer>/
##   tuning.csv, elimination_log.csv, panel.csv, performance.csv,
##   predictions_train_loso.csv, predictions_nested.csv, nested_panels.csv,
##   nested_repeats.csv,
##   permutation.csv, sample_plot.(pdf|png), roc.(pdf|png),
##   panel_boxplots.(pdf|png), summary.txt
## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
source("R/lomics_utils.R"); source("R/lomics_spls.R")
a <- lomics_args(list(config = "study_config.json", layer = NULL, visits = NULL, ncomp = "2",
                      grid = "5,10,15,20,30,50,75,100", zscore = "TRUE", drop_max = "0.05",
                      drop_type = "relative", target_auc = "0.9", elim_auc = "train",
                      nested = "TRUE", nested_grid = "5,10,20,50", nested_folds = "5", nested_repeats = "5", n_perm = "20", panel = NULL,
                      plot_visits = NULL, colour_by = "group", engine = "auto", seed = NULL, tag = NULL))
cfg <- lomics_read_config(a$config)
layer <- if (is.null(a$layer)) cfg$layers$name[1] else a$layer
vis <- if (is.null(a$visits)) cfg$visits[c(1, length(cfg$visits))] else strsplit(a$visits, ",")[[1]]
stopifnot(length(vis) == 2)
num <- function(s) as.numeric(strsplit(s, ",")[[1]])
ncomp <- as.integer(a$ncomp); grid <- num(a$grid); ngrid <- num(a$nested_grid)
target <- suppressWarnings(as.numeric(a$target_auc))
seed <- if (is.null(a$seed)) (if (is.null(cfg$seed)) 1 else cfg$seed) else as.integer(a$seed)
engine <- a$engine
out_dir <- dirname(lomics_out(cfg, "splsda", sprintf("%svs%s_%s%s", vis[1], vis[2], layer,
                                                     if (is.null(a$tag)) "" else paste0("_", a$tag)), "x"))
unlink(list.files(out_dir, full.names = TRUE))   # no stale files from earlier runs

## ---------------------------------------------------------------- data
se <- lomics_load_layer(cfg, layer, "processed")
x_all <- assay(se, "values"); meta_all <- as.data.frame(colData(se))
meta_all$visit <- as.character(meta_all$visit)
if (as.logical(a$zscore)) {                       # paper: z-score every fragment (caret center/scale)
  x_all <- t(scale(t(x_all)))
}
rd <- as.data.frame(rowData(se))
annot_cols <- intersect(c("compound", "mz", "rt_s", "class", "formula"), colnames(rd))

sel <- meta_all$visit %in% vis
paired_subj <- names(which(table(meta_all$subject_id[sel]) == 2))
sel <- sel & meta_all$subject_id %in% paired_subj
meta <- meta_all[sel, ]; meta <- meta[order(meta$subject_id, meta$visit), ]
X <- t(x_all[, meta$sample_id])
y <- factor(meta$visit, levels = vis); subject <- meta$subject_id
Xw <- lomics_within_variation(X, subject)
eng_used <- if (engine == "auto") (if (requireNamespace("mixOmics", quietly = TRUE)) "mixOmics" else "internal") else engine
lines <- c(sprintf("Multilevel sPLS-DA: %s vs %s | layer %s | %d subjects with both visits (%d samples) | %d features",
                   vis[1], vis[2], layer, length(paired_subj), nrow(X), ncol(X)),
           sprintf("engine: %s | ncomp: %d | z-score: %s | seed: %d", eng_used, ncomp, a$zscore, seed))

pick_args <- list(ncomp = ncomp, drop_max = as.numeric(a$drop_max), drop_type = a$drop_type,
                  target_auc = target, auc_on = a$elim_auc, engine = engine)
select_panel <- function(Xw, y, subject, grid) {
  do.call(lomics_spls_select, c(list(Xw = Xw, y = y, subject = subject, grid = grid), pick_args))
}

## ------------------------------------------------- selection on all data
set.seed(seed)
if (is.null(a$panel)) {
  s <- select_panel(Xw, y, subject, grid)
  panel <- s$panel
  utils::write.csv(s$tuning$table, file.path(out_dir, "tuning.csv"), row.names = FALSE)
  utils::write.csv(s$elimination, file.path(out_dir, "elimination_log.csv"), row.names = FALSE)
  lines <- c(lines, sprintf("Tuned keepX: %s (LOSO AUC); start panel %d features -> final panel %d after backward elimination (AUC on %s, drop_max %s %s, target %s)",
                            paste(s$tuning$keepX, collapse = ","), length(s$start_panel), length(panel),
                            a$elim_auc, a$drop_max, a$drop_type, a$target_auc))
} else {
  panel <- strsplit(a$panel, ",")[[1]]
  missing <- setdiff(panel, colnames(Xw)); if (length(missing)) stop("Panel features not found: ", paste(missing, collapse = ", "))
  lines <- c(lines, sprintf("Fixed panel supplied (%d features); no selection.", length(panel)))
}
nc <- min(ncomp, length(panel))
Xp <- Xw[, panel, drop = FALSE]
fit <- lomics_spls_fit(Xp, y, nc, rep(length(panel), nc), engine)

## A. training ROC ---------------------------------------------------------
sc_train <- lomics_spls_score(fit, Xp)
ci_train <- lomics_auc_ci(y, sc_train, vis[2])
## B. LOSO of the fixed panel ----------------------------------------------
sc_loso <- lomics_spls_loso(Xp, y, subject, nc, NULL, engine)
ci_loso <- lomics_auc_ci(y, sc_loso, vis[2])
utils::write.csv(data.frame(sample_id = meta$sample_id, subject_id = subject, visit = y,
                            score_train = sc_train, score_loso_fixed_panel = sc_loso),
                 file.path(out_dir, "predictions_train_loso.csv"), row.names = FALSE)

## C. nested CV: selection repeated inside every outer fold ---------------
##    --nested_folds=loso (one subject out) or K (subject-grouped K-fold,
##    repeated --nested_repeats times). With strong selection, LOSO is
##    biased DOWNWARDS: a held-out non-responder makes its feature look
##    stronger in the rest, so it gets picked and the held-out subject is
##    misclassified. Repeated K-fold reduces this; report both if they differ.
ci_nest <- c(NA, NA, NA); nest_rep <- NULL
if (as.logical(a$nested) && is.null(a$panel)) {
  us <- unique(subject)
  if (identical(a$nested_folds, "loso")) {
    fold_sets <- list(stats::setNames(seq_along(us), us))
  } else {
    K <- as.integer(a$nested_folds); set.seed(seed + 2)
    fold_sets <- lapply(seq_len(as.integer(a$nested_repeats)), function(r)
      stats::setNames(sample(rep_len(seq_len(K), length(us))), us))
  }
  panels <- list(); sc_all <- list()
  for (r in seq_along(fold_sets)) {
    fs <- fold_sets[[r]]; sc_r <- rep(NA_real_, nrow(Xw))
    for (k in sort(unique(fs))) {
      te <- subject %in% names(fs)[fs == k]; tr <- !te
      sn <- select_panel(Xw[tr, , drop = FALSE], y[tr], subject[tr], ngrid)
      ncn <- min(ncomp, length(sn$panel))
      fn <- lomics_spls_fit(Xw[tr, sn$panel, drop = FALSE], y[tr], ncn, rep(length(sn$panel), ncn), engine)
      sc_r[te] <- lomics_spls_score(fn, Xw[te, sn$panel, drop = FALSE])
      panels[[length(panels) + 1]] <- data.frame(repeat_id = r, fold = k,
        held_out = paste(names(fs)[fs == k], collapse = ";"), n_features = length(sn$panel),
        features = paste(sn$panel, collapse = ";"))
    }
    sc_all[[r]] <- sc_r
  }
  nest_rep <- do.call(rbind, lapply(seq_along(sc_all), function(r) {
    pa_r <- lomics_paired_accuracy(y, sc_all[[r]], subject, vis[2])
    data.frame(repeat_id = r, auc = lomics_auc(y, sc_all[[r]], vis[2]), paired_acc = pa_r[["acc"]])
  }))
  sc_nest <- sc_all[[1]]
  ci_nest <- if (length(sc_all) == 1) lomics_auc_ci(y, sc_nest, vis[2]) else
    c(min(nest_rep$auc), mean(nest_rep$auc), max(nest_rep$auc))
  np <- do.call(rbind, panels)
  utils::write.csv(np, file.path(out_dir, "nested_panels.csv"), row.names = FALSE)
  utils::write.csv(nest_rep, file.path(out_dir, "nested_repeats.csv"), row.names = FALSE)
  utils::write.csv(data.frame(sample_id = meta$sample_id, subject_id = subject, visit = y,
                              do.call(cbind, stats::setNames(sc_all, paste0("score_nested_rep", seq_along(sc_all))))),
                   file.path(out_dir, "predictions_nested.csv"), row.names = FALSE)
  freq <- table(unlist(strsplit(np$features, ";")))
  lines <- c(lines, sprintf("Nested CV (%s%s): panel size median %d (range %d-%d); final-panel features re-selected in %s of %d fold fits",
                            if (identical(a$nested_folds, "loso")) "leave-one-subject-out" else paste0(a$nested_folds, "-fold grouped"),
                            if (identical(a$nested_folds, "loso")) "" else paste0(" x ", a$nested_repeats, " repeats"),
                            as.integer(stats::median(np$n_features)), min(np$n_features), max(np$n_features),
                            paste(sprintf("%s %d", panel, as.integer(ifelse(is.na(freq[panel]), 0, freq[panel]))), collapse = ", "),
                            nrow(np)))
}

## D. permutation null (flip visit labels within random subjects) -----------
n_perm <- as.integer(a$n_perm); perm <- NULL
if (n_perm > 0 && is.null(a$panel)) {
  set.seed(seed + 1)
  perm <- do.call(rbind, lapply(seq_len(n_perm), function(i) {
    flip <- sample(c(TRUE, FALSE), length(unique(subject)), replace = TRUE); names(flip) <- unique(subject)
    yp <- y; f <- flip[subject]
    yp[f] <- ifelse(y[f] == vis[1], vis[2], vis[1])
    sp <- select_panel(Xw, yp, subject, grid)
    ncp <- min(ncomp, length(sp$panel)); Xpp <- Xw[, sp$panel, drop = FALSE]
    scp <- lomics_spls_loso(Xpp, yp, subject, ncp, NULL, engine)
    data.frame(perm = i, n_features = length(sp$panel), auc_train = sp$train_auc,
               auc_loso_fixed_panel = lomics_auc(yp, scp, vis[2]),
               paired_acc_loso_fixed_panel = lomics_paired_accuracy(yp, scp, subject, vis[2])[["acc"]])
  }))
  utils::write.csv(perm, file.path(out_dir, "permutation.csv"), row.names = FALSE)
}
emp_p <- function(obs, null) if (is.null(null)) NA else (1 + sum(null >= obs)) / (1 + length(null))

pa <- function(sc) if (all(is.na(sc))) c(acc = NA, ci_low = NA, ci_high = NA, n = NA) else lomics_paired_accuracy(y, sc, subject, vis[2])
pa_tr <- pa(sc_train); pa_lo <- pa(sc_loso)
pa_ne <- if (is.null(nest_rep)) pa(NA) else if (nrow(nest_rep) == 1) pa(sc_nest) else
  c(acc = mean(nest_rep$paired_acc), ci_low = min(nest_rep$paired_acc), ci_high = max(nest_rep$paired_acc), n = length(unique(subject)))
perf <- data.frame(
  estimate = c("A_training (selection on all data)", "B_LOSO_fixed_panel (panel selected on all data)",
               if (identical(a$nested_folds, "loso")) "C_nested_LOSO (selection inside each fold)" else
                 sprintf("C_nested_%sfold_x%s (selection inside each fold; mean [min-max] over repeats)", a$nested_folds, a$nested_repeats)),
  auc = c(ci_train[2], ci_loso[2], ci_nest[2]), auc_ci_low = c(ci_train[1], ci_loso[1], ci_nest[1]),
  auc_ci_high = c(ci_train[3], ci_loso[3], ci_nest[3]),
  paired_acc = c(pa_tr[["acc"]], pa_lo[["acc"]], pa_ne[["acc"]]),
  paired_acc_ci_low = c(pa_tr[["ci_low"]], pa_lo[["ci_low"]], pa_ne[["ci_low"]]),
  paired_acc_ci_high = c(pa_tr[["ci_high"]], pa_lo[["ci_high"]], pa_ne[["ci_high"]]),
  perm_null_mean_auc = c(if (!is.null(perm)) mean(perm$auc_train) else NA,
                         if (!is.null(perm)) mean(perm$auc_loso_fixed_panel) else NA, NA),
  perm_p_auc = c(emp_p(ci_train[2], perm$auc_train), emp_p(ci_loso[2], perm$auc_loso_fixed_panel), NA),
  n_perm = if (is.null(perm)) 0 else nrow(perm))
utils::write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

## ------------------------------------------------------- panel table
wd <- lapply(split(seq_len(nrow(Xp)), subject), function(ix) {
  i1 <- ix[y[ix] == vis[1]]; i2 <- ix[y[ix] == vis[2]]; X[i2, panel] - X[i1, panel]
})
dmat <- matrix(unlist(wd), ncol = length(panel), byrow = TRUE, dimnames = list(names(wd), panel))
ptab <- data.frame(feature = panel, rd[panel, annot_cols, drop = FALSE],
                   loading_comp1 = round(fit$loadings[panel, 1], 3),
                   mean_change_z = round(colMeans(dmat), 3),
                   direction = ifelse(colMeans(dmat) > 0, "increase", "decrease"),
                   n_up = colSums(dmat > 0), n_down = colSums(dmat < 0),
                   wilcoxon_p_descriptive = signif(apply(dmat, 2, function(d) stats::wilcox.test(d)$p.value), 3),
                   row.names = NULL, check.names = FALSE)
utils::write.csv(ptab, file.path(out_dir, "panel.csv"), row.names = FALSE)

## ------------------------------------------------------- plots
pv <- if (is.null(a$plot_visits)) vis else strsplit(a$plot_visits, ",")[[1]]
grp <- if (a$colour_by %in% colnames(meta_all)) a$colour_by else "group"
plot_samples <- function() {
  V <- fit$variates
  if (ncol(V) < 2) V <- cbind(V, 0)
  cols <- c("#2C7BB6", "#E66101")[as.integer(y)]
  graphics::plot(V[, 1], V[, 2], col = cols, pch = 16, xlab = "Component 1", ylab = "Component 2",
                 main = sprintf("Multilevel sPLS-DA (%d features)", length(panel)))
  graphics::legend("topright", legend = vis, col = c("#2C7BB6", "#E66101"), pch = 16, bty = "n")
}
plot_roc <- function() {
  r1 <- pROC::roc(as.integer(y == vis[2]), sc_train, direction = "<", quiet = TRUE)
  r2 <- pROC::roc(as.integer(y == vis[2]), sc_loso, direction = "<", quiet = TRUE)
  graphics::plot(1 - r1$specificities, r1$sensitivities, type = "l", lwd = 2, col = "grey40",
                 xlab = "1 - specificity", ylab = "Sensitivity", main = sprintf("%s vs %s", vis[1], vis[2]))
  graphics::lines(1 - r2$specificities, r2$sensitivities, lwd = 2, col = "#E66101")
  leg <- c(sprintf("training AUC %.2f", ci_train[2]), sprintf("LOSO fixed panel %.2f", ci_loso[2]))
  cl <- c("grey40", "#E66101")
  if (exists("sc_nest")) {
    r3 <- pROC::roc(as.integer(y == vis[2]), sc_nest, direction = "<", quiet = TRUE)
    graphics::lines(1 - r3$specificities, r3$sensitivities, lwd = 2, col = "#5E3C99", lty = 2)
    leg <- c(leg, sprintf("nested CV%s %.2f", if (identical(a$nested_folds, "loso")) "" else " (repeat 1)", lomics_auc(y, sc_nest, vis[2]))); cl <- c(cl, "#5E3C99")
  }
  graphics::abline(0, 1, lty = 3); graphics::legend("bottomright", leg, col = cl, lwd = 2, bty = "n")
}
plot_box <- function() {
  mm <- meta_all[meta_all$visit %in% pv, ]; mm$visit <- factor(mm$visit, levels = pv)
  k <- length(panel); nr <- ceiling(k / 3)
  op <- graphics::par(mfrow = c(nr, min(3, k)), mar = c(3, 4, 3, 1))
  gl <- factor(mm[[grp]]); gcols <- grDevices::hcl.colors(nlevels(gl), "Dark 3")
  for (f in panel) {
    v <- x_all[f, mm$sample_id]
    ttl <- if ("compound" %in% annot_cols && rd[f, "compound"] != "unidentified") rd[f, "compound"] else f
    graphics::boxplot(v ~ mm$visit, outline = FALSE, main = ttl, ylab = "z-score", xlab = "", col = "grey92",
                      ylim = range(v, na.rm = TRUE), cex.main = 0.9)
    for (s in unique(mm$subject_id)) {
      ix <- which(mm$subject_id == s); ix <- ix[order(mm$visit[ix])]
      graphics::lines(as.integer(mm$visit[ix]), v[ix], col = grDevices::adjustcolor(gcols[gl[ix[1]]], 0.6))
    }
    graphics::points(jitter(as.integer(mm$visit), 0.3), v, pch = 16, cex = 0.5, col = gcols[gl])
    if (f == panel[1]) graphics::legend("topright", legend = levels(gl), col = gcols, lwd = 2, bty = "n", cex = 0.7)
  }
  graphics::par(op)
}
save_plot <- function(fun, name, w, h) {
  grDevices::pdf(file.path(out_dir, paste0(name, ".pdf")), width = w, height = h); fun(); invisible(grDevices::dev.off())
  grDevices::png(file.path(out_dir, paste0(name, ".png")), width = w * 150, height = h * 150, res = 150); fun(); invisible(grDevices::dev.off())
}
save_plot(plot_samples, "sample_plot", 6, 5)
if (requireNamespace("pROC", quietly = TRUE)) save_plot(plot_roc, "roc", 5.5, 5.5)
save_plot(plot_box, "panel_boxplots", 10, 3.3 * ceiling(length(panel) / 3))

lines <- c(lines, "", "Performance:", utils::capture.output(print(perf, row.names = FALSE, digits = 3)),
           "", "Panel:", utils::capture.output(print(ptab, row.names = FALSE)),
           "", "Read A and B as optimistic: the panel was chosen on the same subjects.",
           "C (nested) is the honest estimate; D shows what A and B look like with no real effect.")
writeLines(lines, file.path(out_dir, "summary.txt"))
lomics_session_log(cfg, sprintf("splsda_%svs%s_%s", vis[1], vis[2], layer))
cat(lines, sep = "\n")
