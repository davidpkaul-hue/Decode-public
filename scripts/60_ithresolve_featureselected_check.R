#!/usr/bin/env Rscript
################################################################################
# Script 60: Sensitivity check of the ITHresolveGBM differentiation score
# ----------------------------------------------------------------------------
# The manuscript's differentiation analysis (Figure 4) projects the ITHresolveGBM
# component profiles with the file `profiles_full`. The authors' own reference
# script (apply_deconvolute.Rmd, line 55) recommends a different file for this
# purpose:
#
#   "profiles_after_feature_selection_with_sc ... These profiles are better for
#    gaining cell state resolution ... If you are interested in cell types and
#    microenviorement composition more you can try profiles_full"
#   "IMPORTANT! Make sure the number of CpGs overlapping is 80K or more"
#
# A differentiation score IS a cell-state readout, so this reruns the projection
# with the feature-selected profiles and compares every number the manuscript
# quotes. No classifier run is needed; the purified calls are reused.
################################################################################
suppressPackageStartupMessages({ library(EpiDISH) })

FS   <- "data/external/ithresolve/profiles_after_feature_selection_with_sc"
BETA <- "/mnt/volume-hel1-1/siewert/betamatrix.rds"
OUT  <- "data/processed/ithresolve_featureselected_diffscore.csv"

simplify <- function(x) {
  out <- rep("other", length(x))
  out[grepl("RTK1", x)] <- "RTK I"
  out[grepl("RTK2", x)] <- "RTK II"
  out[grepl("MES|mesenchymal", x, ignore.case = TRUE)] <- "MES"
  out[grepl("PNC|primitive neuronal", x, ignore.case = TRUE)] <- "PNC"
  out
}

message("Loading feature-selected profiles ...")
prof <- read.table(FS, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message(sprintf("  %d probes x 6 components", nrow(prof)))

message("Loading beta matrix ...")
beta <- readRDS(BETA)
common <- intersect(prof$probe, rownames(beta))
message(sprintf("  overlap with our arrays: %d CpGs (authors require >= 80,000)", length(common)))
stopifnot(length(common) >= 80000)

beta <- beta[common, ]
ref  <- as.matrix(prof[match(common, prof$probe), c("immune","neuron","glia","V4","V5","V6")])
rownames(ref) <- common
keep <- !(rowSums(is.na(beta)) > 0 | rowSums(is.na(ref)) > 0)
beta <- beta[keep, ]; ref <- ref[keep, ]
message(sprintf("  after NA removal: %d CpGs x %d samples", nrow(beta), ncol(beta)))

message("EpiDISH RPC (maxit = 2000), same call as the authors' reference script ...")
fr <- as.data.frame(epidish(beta.m = beta, ref.m = ref, method = "RPC", maxit = 2000)$estF)
colnames(fr) <- c("Immune","Neuron","Glia","Stem","Diff2","Diff1")   # V4, V5, V6
fr$Idat <- rownames(fr)
fr$Purity <- fr$Stem + fr$Diff1 + fr$Diff2
fr$DiffScore <- fr$Stem - (fr$Diff1 + fr$Diff2)      # V4 - (V5 + V6), their code
write.csv(fr, OUT, row.names = FALSE)
message("  wrote ", OUT)

# ---------------------------------------------------------------- comparison
old <- read.csv("data/processed/siewert_ithresolve_full.csv", stringsAsFactors = FALSE)
meta <- read.csv("data/processed/siewert_cohort_results.csv", stringsAsFactors = FALSE, check.names = FALSE)
meta <- meta[grepl("Glioblastoma, IDH", meta$Max_family), ]
meta$orig <- simplify(meta$Max_subclass)
source("scripts/_purification_config.R")
calls <- read.csv(PURIFIED_CALLS, stringsAsFactors = FALSE)
calls$pur <- simplify(calls$pred_sub)

report <- function(score, idat, label) {
  d <- data.frame(Idat = idat, S = score, stringsAsFactors = FALSE)
  d$orig <- meta$orig[match(d$Idat, meta$Idat)]
  d$pur  <- calls$pur[match(d$Idat, calls$Sample)]
  d <- d[!is.na(d$S) & !is.na(d$orig) & !is.na(d$pur), ]
  d$T <- paste(d$orig, "->", d$pur)
  k <- c("RTK I -> RTK I", "MES -> RTK I", "MES -> RTK II", "RTK II -> RTK II")
  d <- d[d$T %in% k, ]; d$T <- factor(d$T, levels = k)
  cat("\n---", label, "---\n")
  for (g in k) cat(sprintf("  %-18s n=%3d  median %+.3f\n", g, sum(d$T == g), median(d$S[d$T == g])))
  cat(sprintf("  Kruskal-Wallis p = %.2e\n", kruskal.test(S ~ T, data = d)$p.value))
  d
}
a <- report(old$DiffScore, old$Idat, "profiles_full (manuscript)")
b <- report(fr$DiffScore,  fr$Idat,  "profiles_after_feature_selection_with_sc")

m <- merge(data.frame(Idat = old$Idat, full = old$DiffScore),
           data.frame(Idat = fr$Idat,  fs   = fr$DiffScore), by = "Idat")
cat(sprintf("\nPer-sample agreement: Pearson r = %.3f, Spearman rho = %.3f (n = %d)\n",
            cor(m$full, m$fs), cor(m$full, m$fs, method = "spearman"), nrow(m)))
cat(sprintf("Sign agreement (stem-like vs differentiated): %.1f%%\n",
            100 * mean((m$full > 0) == (m$fs > 0))))
