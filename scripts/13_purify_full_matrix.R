#!/usr/bin/env Rscript
################################################################################
# DECODE: In-Silico Purification of FULL Beta Matrix
################################################################################
# The Heidelberg classifier requires ~10,000 specific CpG features.
# Our previous purification only covered 128k CpGs (ITH profile overlap).
# This script purifies ALL 721k CpGs using chunked processing to fit in 8GB RAM.
#
# Method: Extended InfiniumPurify (Silverbush et al.)
#   For each CpG: Y = asin(2*beta - 1)
#   Fit: Y ~ immune + neuron + glia (per-sample fractions from ITH deconvolution)
#   Purified: sin(Y - b1*immune - b2*neuron - b3*glia) + 1) / 2
#
# Input:  /mnt/volume-hel1-1/siewert/betamatrix.rds (721k x 436)
#         data/processed/siewert_ithresolve_full.csv (ITH fractions)
# Output: /mnt/volume-hel1-1/siewert/betamatrix_purified_full.rds
#         /mnt/volume-hel1-1/siewert/methylset_purified_full_v1.rds
#         /mnt/volume-hel1-1/siewert/methylset_purified_full_v2.rds
################################################################################

if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

suppressPackageStartupMessages({
  library(minfi)
})

message("================================================================")
message("In-Silico Purification: Full Beta Matrix (721k CpGs)")
message("================================================================\n")

# ============================================================================
# Load ITH deconvolution fractions (immune, neuron, glia per sample)
# ============================================================================

message("Loading ITH deconvolution results...")
ith <- read.csv("data/processed/siewert_ithresolve_full.csv", stringsAsFactors = FALSE)
message(sprintf("  %d samples with ITH fractions", nrow(ith)))

# ============================================================================
# Load full beta matrix
# ============================================================================

message("\nLoading full beta matrix...")
beta_full <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
message(sprintf("  %d CpGs x %d samples", nrow(beta_full), ncol(beta_full)))

# Align samples
samples <- intersect(colnames(beta_full), ith$Idat)
beta_full <- beta_full[, samples]
ith <- ith[match(samples, ith$Idat), ]
message(sprintf("  Aligned: %d samples", length(samples)))

# Non-malignant fractions for regression
A <- ith$Immune
B <- ith$Neuron
C <- ith$Glia

# ============================================================================
# Chunked purification to stay within 8GB RAM
# ============================================================================

.myasin <- function(x) asin(pmin(pmax(2 * x - 1, -1), 1))

# Design matrix (nsamples x 4)
X_design <- cbind(1, A, B, C)
XtX_inv <- solve(t(X_design) %*% X_design)
XtX_inv_Xt <- XtX_inv %*% t(X_design)  # 4 x nsamples — precompute once

n_probes <- nrow(beta_full)
chunk_size <- 30000  # Process 30k CpGs at a time
n_chunks <- ceiling(n_probes / chunk_size)

message(sprintf("\nPurifying %d CpGs in %d chunks of %d (in-place)...", n_probes, n_chunks, chunk_size))

# Purify IN PLACE to avoid allocating a second 2.5GB matrix
for (ch in seq_len(n_chunks)) {
  start_idx <- (ch - 1) * chunk_size + 1
  end_idx <- min(ch * chunk_size, n_probes)
  chunk_rows <- start_idx:end_idx

  # Work on chunk
  beta_chunk <- beta_full[chunk_rows, , drop = FALSE]
  na_mask <- is.na(beta_chunk)
  beta_chunk[na_mask] <- 0.5

  Y_chunk <- .myasin(beta_chunk)
  beta_hat <- XtX_inv_Xt %*% t(Y_chunk)
  effects <- outer(beta_hat[2, ], A) + outer(beta_hat[3, ], B) + outer(beta_hat[4, ], C)
  Y_purified <- Y_chunk - effects

  bp <- (sin(Y_purified) + 1) / 2
  bp[bp < 0] <- 0
  bp[bp > 1] <- 1
  bp[na_mask] <- NA

  # Write back in place
  beta_full[chunk_rows, ] <- bp

  message(sprintf("  Chunk %d/%d: CpGs %d-%d done", ch, n_chunks, start_idx, end_idx))
  rm(beta_chunk, Y_chunk, beta_hat, effects, Y_purified, bp, na_mask)
  gc(verbose = FALSE)
}

# beta_full is now the purified matrix
beta_purified <- beta_full
rm(beta_full)
gc(verbose = FALSE)

message(sprintf("\nPurified matrix: %d CpGs x %d samples", nrow(beta_purified), ncol(beta_purified)))

# Check how many NAs remain
n_na <- sum(is.na(beta_purified))
message(sprintf("  NA values: %d (%.1f%%)", n_na, n_na / length(beta_purified) * 100))

# ============================================================================
# Save purified beta matrix
# ============================================================================

out_beta <- "/mnt/volume-hel1-1/siewert/betamatrix_purified_full.rds"
message(sprintf("\nSaving purified beta matrix..."))
saveRDS(beta_purified, out_beta)
message(sprintf("  Saved: %s (%.0f MB)", out_beta, file.size(out_beta) / 1e6))

# ============================================================================
# Convert to MethylSets (v1 and v2 separately)
# ============================================================================

message("\nCreating MethylSet objects...")

sentrix <- as.numeric(sub("_.*", "", colnames(beta_purified)))
is_v2 <- sentrix >= 207000000000

TOTAL <- 10000

for (arr in c("v1", "v2")) {
  idx <- if (arr == "v1") !is_v2 else is_v2
  b <- beta_purified[, idx]
  n <- ncol(b)
  message(sprintf("  %s: %d samples", arr, n))

  pdata <- data.frame(Sample_Name = colnames(b), row.names = colnames(b))

  ann <- if (arr == "v1") {
    c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
  } else {
    c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  }

  # Handle NAs for MethylSet: replace with 0
  b_clean <- b
  b_clean[is.na(b_clean)] <- 0.5

  mset <- MethylSet(
    Meth = as.matrix(round(b_clean * TOTAL)),
    Unmeth = as.matrix(round((1 - b_clean) * TOTAL)),
    annotation = ann,
    colData = DataFrame(pdata)
  )

  outpath <- sprintf("/mnt/volume-hel1-1/siewert/methylset_purified_full_%s.rds", arr)
  saveRDS(mset, outpath)
  message(sprintf("  Saved: %s (%.0f MB)", outpath, file.size(outpath) / 1e6))

  rm(b, b_clean, mset)
  gc(verbose = FALSE)
}

message("\n================================================================")
message("Done! Full purified matrix + MethylSets ready for classifier.")
message("================================================================")
