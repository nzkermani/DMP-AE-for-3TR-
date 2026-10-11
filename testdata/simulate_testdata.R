## Synthetic 4-visit, 3-omic study used to test every skill script.
## Usage: Rscript testdata/simulate_testdata.R testdata
##
## Planted signals (so each method has something known to recover):
##   * Treatment arm: first 10% of features in each layer rise at V2-V3 and
##     return at V4 (trajectory, longitudinal DE, interaction).
##   * Subject-level outcomes depend on baseline (V1) values of
##     PROT20-PROT24 and GENE50-GENE54:
##       responder (binary), fev1 (continuous), time_to_exacerbation /
##       exacerbation (right-censored survival).
##   * Two latent baseline "endotypes" shift GENE100-GENE139 and MET30-MET39
##     (sample clustering, dimension reduction).
set.seed(1)
out <- commandArgs(TRUE)[1]; if (is.na(out)) out <- "testdata"
n_subj <- 60
subj <- sprintf("S%02d", seq_len(n_subj)); visits <- c("V1", "V2", "V3", "V4")
meta <- expand.grid(subject_id = subj, visit = visits, stringsAsFactors = FALSE)
meta$group <- ifelse(as.integer(sub("S", "", meta$subject_id)) <= n_subj / 2, "Placebo", "Treatment")
meta$sample_id <- paste(meta$subject_id, meta$visit, sep = "_")
meta <- meta[-c(5, 30, 77, 150, 201), ]                 # a few missed visits
subj_info <- data.frame(subject_id = subj,
                        age = round(runif(n_subj, 20, 70)),
                        sex = sample(c("F", "M"), n_subj, TRUE),
                        endotype = sample(c("E1", "E2"), n_subj, TRUE))
meta <- merge(meta, subj_info, by = "subject_id")
meta$batch <- sample(c("B1", "B2"), nrow(meta), TRUE)
tvec <- c(V1 = 0, V2 = 1, V3 = 2, V4 = 3)

mk <- function(nf, prefix, endo_rows = integer(), counts = FALSE, missing = 0) {
  subj_eff <- matrix(rnorm(nf * n_subj, sd = .5), nf, n_subj, dimnames = list(NULL, subj))
  m <- sapply(seq_len(nrow(meta)), function(i) {
    s <- meta$subject_id[i]; t <- tvec[meta$visit[i]]
    n_trt <- nf / 10
    eff <- c(rep(0.8, n_trt) * (meta$group[i] == "Treatment") * c(0, 1, 1, 0.2)[t + 1],
             rep(0, nf - n_trt))
    eff[endo_rows] <- eff[endo_rows] + 1.2 * (meta$endotype[i] == "E2")
    5 + subj_eff[, s] + eff + rnorm(nf, sd = .4)
  })
  rownames(m) <- paste0(prefix, seq_len(nf)); colnames(m) <- meta$sample_id
  if (counts) m <- matrix(rpois(length(m), 2^m), nrow(m), dimnames = dimnames(m))
  if (missing > 0) m[sample(length(m), missing * length(m))] <- NA
  m
}
tx  <- mk(400, "GENE", endo_rows = 100:139, counts = TRUE)
pr  <- mk(150, "PROT", missing = .05)
met <- mk(100, "MET", endo_rows = 30:39, missing = .1)

## Subject-level outcomes from baseline values
v1 <- meta[meta$visit == "V1", ]
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
score <- z(colMeans(pr[20:24, v1$sample_id], na.rm = TRUE)) +
         z(colMeans(log2(tx[50:54, v1$sample_id] + 1)))
names(score) <- v1$subject_id
score <- score[subj]; score[is.na(score)] <- 0              # subjects without V1
oc <- data.frame(subject_id = subj,
                 responder = ifelse(plogis(1.5 * score) > runif(n_subj), "Yes", "No"),
                 fev1 = round(2.8 + 0.35 * score + rnorm(n_subj, sd = .3), 2))
haz <- 0.05 * exp(0.9 * score)
t_event <- rexp(n_subj, haz); t_cens <- runif(n_subj, 5, 40)
oc$time_to_exacerbation <- round(pmin(t_event, t_cens), 1)
oc$exacerbation <- as.integer(t_event <= t_cens)
meta <- merge(meta, oc, by = "subject_id")
meta <- meta[order(meta$subject_id, meta$visit), c("sample_id", setdiff(names(meta), "sample_id"))]

dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(meta, file.path(out, "metadata.csv"), row.names = FALSE)
write.csv(tx,  file.path(out, "transcriptomics_counts.csv"))
write.csv(pr,  file.path(out, "proteomics_log2.csv"))
write.csv(2^met, file.path(out, "metabolomics_raw.csv"))
cat("Wrote", nrow(meta), "samples from", n_subj, "subjects to", out, "\n")
