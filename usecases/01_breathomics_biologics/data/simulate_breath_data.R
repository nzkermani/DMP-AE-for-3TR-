#!/usr/bin/env Rscript
## =============================================================================
## simulate_breath_data.R — synthetic 3TR-ABC-style exhaled-breath GC-MS study
##
## Mimics the design of Vassileva et al., Respir Med 262 (2026) 109067
## ("The impact of benralizumab and dupilumab on exhaled VOC profiles of severe
## asthma patients"), for which patient-level data are not public.
##
##   * 25 severe-asthma patients, 2 sites (Amsterdam, Copenhagen)
##   * biologic: dupilumab (n = 20), benralizumab (n = 5)
##   * visits: BL (baseline), M1 (1 month, 15 patients), M4 (4 months)
##   * clinical data matching the summary statistics in the paper's Tables 1 & 3
##   * GC-MS feature table AFTER XCMS peak picking/alignment: 4270 ion fragments
##     (features x samples, raw intensities), of which 1062 are flagged as
##     known contaminants (siloxanes, phthalates, Tenax artefacts)
##   * 4 GC-MS analysis batches (aligned with site, crossed with visit)
##
## PLANTED SIGNALS (known truth, used to check the pipeline):
##   time effect BL -> M4 (half at M1), identifying fragments as in Table 2:
##     MZ052_RT0496 methyl isobutyl ketone  down
##     MZ039_RT0933 2-ethyl-1-hexanol       up
##     MZ098_RT1193 5-methyloctadecane      down
##     MZ104_RT1080 2-methylundecane        down
##     MZ141_RT1143 1,1'-oxybis(decane)     down
##   plus 15 weaker time-responsive "background" fragments (unnamed)
##   baseline trait associations (constant within subject, so invisible to the
##   multilevel model but visible in baseline correlations):
##     2-ethyl-1-hexanol   ~ + baseline ACQ-5
##     1,1'-oxybis(decane) ~ - baseline blood eosinophils
##   no planted difference between biologics in VOC change.
##
## Everything here is SYNTHETIC. Numbers produced by the notebook describe this
## simulation, not the 3TR-ABC cohort.
##
## Usage (from the project root):
##   Rscript usecases/01_breathomics_biologics/data/simulate_breath_data.R [out_dir] [seed]
## =============================================================================

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "usecases/01_breathomics_biologics/data"
seed    <- if (length(args) >= 2) as.integer(args[2]) else 2026L
set.seed(seed)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)
n <- 25
subj <- sprintf("ABC%03d", seq_len(n))

## ------------------------------------------------------------- patient level
site     <- sample(rep(c("Amsterdam", "Copenhagen"), c(14, 11)))
biologic <- sample(rep(c("dupilumab", "benralizumab"), c(20, 5)))
sex      <- sample(rep(c("female", "male"), c(8, 17)))
age      <- round(clip(rnorm(n, 61, 16), 22, 85))
onset    <- round(pmin(clip(rnorm(n, 43, 22), 3, 80), age - 1))
bmi      <- round(clip(rnorm(n, 26, 4), 18, 38), 1)
smoking  <- sample(rep(c("never", "former", "passive"), c(13, 11, 1)))
pack_yrs <- ifelse(smoking == "former", round(rlnorm(n, log(11), 0.9), 1), NA)
exac     <- sample(c(1, 2, 2, 2, 3, 3, 3, 3, 4, 5), n, replace = TRUE)
ics_dose <- round(clip(rnorm(n, 1610, 370), 800, 2400) / 100) * 100
lama     <- sample(rep(c("yes", "no"), c(12, 13)))
ltra     <- sample(rep(c("yes", "no"), c(17, 8)))
ocs_m    <- sample(rep(c("yes", "no"), c(4, 21)))
ocs_dose <- ifelse(ocs_m == "yes", round(clip(rnorm(n, 16, 14), 5, 40) / 2.5) * 2.5, NA)
fvc_pp   <- round(clip(rnorm(n, 95, 12), 65, 125), 1)
fev1_fvc <- round(clip(rnorm(n, 63, 12.5), 35, 85), 1)

## clinical trajectories (BL, M1, M4); dupilumab vs benralizumab biology
acq_bl  <- clip(rnorm(n, 2.89, 1.02), 0.6, 5.4)
acq_m4  <- clip(acq_bl - 1.5 + rnorm(n, 0, 0.9), 0, 5)
aqlq_bl <- clip(rnorm(n, 4.61, 0.94), 2, 6.8)
aqlq_m4 <- clip(aqlq_bl + 1.1 + rnorm(n, 0, 0.8), 2, 7)
fev_bl  <- clip(rnorm(n, 77, 18), 35, 120)
fev_m1  <- clip(fev_bl + 9 + rnorm(n, 0, 7), 35, 130)
fev_m4  <- clip(fev_bl + 7 + rnorm(n, 0, 6), 35, 130)
dupi    <- biologic == "dupilumab"
eos_bl  <- exp(rnorm(n, log(0.30), 0.85))
eos_m1  <- ifelse(dupi, eos_bl * exp(rnorm(n, 0.30, 0.5)), runif(n, 0, 0.03))
eos_m4  <- ifelse(dupi, eos_bl * exp(rnorm(n, 0.20, 0.5)), runif(n, 0, 0.04))
ige_bl  <- exp(rnorm(n, log(80), 1.4))
ige_m1  <- ige_bl * ifelse(dupi, 0.8, 1) * exp(rnorm(n, 0, 0.2))
ige_m4  <- ige_bl * ifelse(dupi, 0.5, 1) * exp(rnorm(n, 0, 0.25))
feno_bl <- exp(rnorm(n, log(39), 0.7))
feno_m1 <- feno_bl * ifelse(dupi, 0.65, 0.95) * exp(rnorm(n, 0, 0.3))
feno_m4 <- feno_bl * ifelse(dupi, 0.6, 0.95) * exp(rnorm(n, 0, 0.3))
miss_feno <- sample(n, 2); miss_ige <- sample(n, 2)

## which patients have the 1-month breath sample (15 of 25, unrelated to clinic)
has_m1 <- seq_len(n) %in% sample(n, 15)

## --------------------------------------------------------------- sample level
rows <- list()
for (i in seq_len(n)) {
  for (v in c("BL", "M1", "M4")) {
    if (v == "M1" && !has_m1[i]) next
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = paste0(subj[i], "_", v), subject_id = subj[i], visit = v,
      biologic = biologic[i], site = site[i], sex = sex[i], age = age[i],
      age_onset = onset[i], bmi = bmi[i], smoking = smoking[i], pack_years = pack_yrs[i],
      exacerbations_prev_year = exac[i],
      acq5 = round(switch(v, BL = acq_bl[i], M1 = NA, M4 = acq_m4[i]), 2),
      aqlq = round(switch(v, BL = aqlq_bl[i], M1 = NA, M4 = aqlq_m4[i]), 2),
      blood_eos = round(switch(v, BL = eos_bl[i], M1 = eos_m1[i], M4 = eos_m4[i]), 2),
      ige = if (i %in% miss_ige) NA else round(switch(v, BL = ige_bl[i], M1 = ige_m1[i], M4 = ige_m4[i])),
      feno = if (i %in% miss_feno) NA else round(switch(v, BL = feno_bl[i], M1 = feno_m1[i], M4 = feno_m4[i]), 1),
      fev1_pp = round(switch(v, BL = fev_bl[i], M1 = fev_m1[i], M4 = fev_m4[i]), 1),
      fvc_pp = fvc_pp[i], fev1_fvc = fev1_fvc[i],
      ics_dose = ics_dose[i], lama = lama[i], ltra = ltra[i],
      ocs_maintenance = ocs_m[i], ocs_dose = ocs_dose[i],
      stringsAsFactors = FALSE)
  }
}
meta <- do.call(rbind, rows)
## clinically meaningful ACQ response (decrease >= 0.5) — handy as an outcome
resp <- stats::setNames(ifelse(acq_bl - acq_m4 >= 0.5, "responder", "non_responder"), subj)
meta$acq_responder <- resp[meta$subject_id]

## GC-MS analysis batch: by analysis date, aligned with site, crossed with visit
meta$gcms_batch <- NA_character_
for (s in unique(meta$site)) {
  idx <- which(meta$site == s)
  b <- if (s == "Amsterdam") c("B1", "B2") else c("B3", "B4")
  meta$gcms_batch[idx] <- sample(rep(b, length.out = length(idx)))
}
S <- nrow(meta)

## -------------------------------------------------------- feature definition
P <- 4270; n_contam <- 1062
mz <- sample(37:300, P, replace = TRUE)
rt <- round(runif(P, 300, 1500))
planted <- data.frame(
  compound = c("methyl isobutyl ketone", "2-ethyl-1-hexanol", "5-methyloctadecane",
               "2-methylundecane", "1,1'-oxybis(decane)"),
  class = c("Ketones", "Fatty alcohols", "Alkanes", "Alkanes", "Ethers"),
  formula = c("C6H12O", "C8H18O", "C19H40", "C12H26", "C20H42O"),
  mz = c(52, 39, 98, 104, 141), rt = c(496, 933, 1193, 1080, 1143),
  effect_m4 = c(-0.80, 0.80, -0.75, -0.75, -0.75),   # log-units, BL -> M4
  stringsAsFactors = FALSE)
mz[1:5] <- planted$mz; rt[1:5] <- planted$rt
feat <- sprintf("MZ%03d_RT%04d", mz, rt)
while (anyDuplicated(feat)) {                         # make IDs unique
  d <- which(duplicated(feat)); rt[d] <- rt[d] + 1L
  feat <- sprintf("MZ%03d_RT%04d", mz, rt)
}
contam_idx <- sample(6:P, n_contam)
bg_idx     <- sample(setdiff(6:P, contam_idx), 15)   # weaker background responders
contam_names <- c("hexamethylcyclotrisiloxane", "octamethylcyclotetrasiloxane",
                  "decamethylcyclopentasiloxane", "dodecamethylcyclohexasiloxane",
                  "diethyl phthalate", "dibutyl phthalate", "benzaldehyde (Tenax artefact)",
                  "acetophenone (Tenax artefact)", "phenol (Tenax artefact)")
annot <- data.frame(feature = feat, mz = mz, rt_s = rt, compound = "unidentified",
                    class = NA_character_, formula = NA_character_,
                    contaminant = FALSE, stringsAsFactors = FALSE)
annot[1:5, c("compound", "class", "formula")] <- planted[, c("compound", "class", "formula")]
annot$contaminant[contam_idx] <- TRUE
annot$compound[contam_idx] <- sample(contam_names, n_contam, replace = TRUE)
annot$class[contam_idx] <- ifelse(grepl("siloxane", annot$compound[contam_idx]), "Siloxanes",
                           ifelse(grepl("phthalate", annot$compound[contam_idx]), "Phthalates", "Sorbent artefacts"))

## --------------------------------------------------------- intensity model
## log intensity = feature mean + subject trait (low-rank + idiosyncratic)
##               + batch shift + visit effects + noise
mu <- rnorm(P, 9, 1.6); mu[contam_idx] <- mu[contam_idx] + 1.5
subj_ix <- match(meta$subject_id, subj)
K <- 6
Fs <- matrix(rnorm(n * K), n, K)                       # subject latent traits
Ls <- matrix(rnorm(P * K, 0, 0.25), P, K)
subj_eff <- Fs %*% t(Ls) + matrix(rnorm(n * P, 0, 0.45), n, P)   # n x P
batches <- sort(unique(meta$gcms_batch))
batch_eff <- matrix(rnorm(length(batches) * P, 0, 0.5), length(batches), P,
                    dimnames = list(batches, NULL))
batch_eff[, contam_idx] <- batch_eff[, contam_idx] * 3   # contaminants are batchy
noise <- matrix(rnorm(S * P, 0, 0.35), S, P)
L <- matrix(mu, S, P, byrow = TRUE) + subj_eff[subj_ix, ] +
     batch_eff[meta$gcms_batch, ] + noise

tfrac <- c(BL = 0, M1 = 0.5, M4 = 1)[meta$visit]
## heterogeneous response: each patient responds on each planted VOC with
## probability 0.75 (independently across VOCs), so single VOCs separate the
## visits only moderately and the panel does better than any one VOC
resp_vp <- matrix(stats::rbinom(n * 5, 1, 0.8), n, 5)
for (j in 1:5) L[, j] <- L[, j] + tfrac * planted$effect_m4[j] * resp_vp[subj_ix, j]
bg_eff <- sample(c(-1, 1), 15, replace = TRUE) * runif(15, 0.1, 0.2)
for (k in seq_along(bg_idx)) L[, bg_idx[k]] <- L[, bg_idx[k]] + tfrac * bg_eff[k]

## baseline traits (constant within subject)
z <- function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE)
L[, 2] <- L[, 2] + 0.65 * z(acq_bl)[subj_ix]           # 2-ethyl-1-hexanol ~ ACQ
L[, 5] <- L[, 5] - 0.90 * z(log(eos_bl))[subj_ix]      # oxybis(decane) ~ -eos

X <- exp(L)
## values below the detection limit -> 0 (as in an XCMS table without fill-in)
lod <- stats::quantile(X, 0.03)
X[X < lod] <- 0
X <- round(t(X))                                       # features x samples
dimnames(X) <- list(feat, meta$sample_id)

## ------------------------------------------------------------------- write
utils::write.csv(meta, file.path(out_dir, "metadata.csv"), row.names = FALSE)
utils::write.csv(data.frame(feature = feat, X, check.names = FALSE),
                 file.path(out_dir, "breath_gcms_raw.csv"), row.names = FALSE)
utils::write.csv(annot, file.path(out_dir, "feature_annotation.csv"), row.names = FALSE)
truth <- rbind(
  data.frame(feature = feat[1:5], compound = planted$compound, type = "time_effect",
             effect_log_M4_vs_BL = planted$effect_m4),
  data.frame(feature = feat[bg_idx], compound = "unidentified", type = "time_effect_background",
             effect_log_M4_vs_BL = round(bg_eff, 3)),
  data.frame(feature = feat[c(2, 5)], compound = planted$compound[c(2, 5)],
             type = c("baseline_trait: + ACQ-5", "baseline_trait: - blood eosinophils"),
             effect_log_M4_vs_BL = NA))
utils::write.csv(truth, file.path(out_dir, "planted_truth.csv"), row.names = FALSE)
message(sprintf("Wrote %d samples (%d subjects) x %d fragments (%d contaminants) to %s [seed %d]",
                S, n, P, n_contam, out_dir, seed))
