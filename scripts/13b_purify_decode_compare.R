# 13b: DECODE-fraction in-silico purification, 3-component vs 4-component.
# Compares purified matrices (3-comp immune/neuron/glia vs 4-comp +stromal),
# focus on MES samples. Saves the primary 3-component matrix for classification.
suppressWarnings(suppressMessages(library(stats)))

myasin <- function(x) asin(pmin(pmax(2 * x - 1, -1), 1))

message("Loading DECODE fractions + MES labels...")
dec <- read.csv("/mnt/volume-hel1-1/data_processed/siewert_cohort_results.csv", stringsAsFactors = FALSE)
ith <- read.csv("/mnt/volume-hel1-1/data_processed/siewert_ithresolve_full.csv", stringsAsFactors = FALSE)

message("Loading raw beta matrix (2.3 GB)...")
beta <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
message(sprintf("  %d CpGs x %d samples", nrow(beta), ncol(beta)))

# Align samples to DECODE fractions
samp <- intersect(colnames(beta), dec$Idat)
beta <- beta[, samp]
di <- match(samp, dec$Idat)
immune <- dec$Immune_total[di]
neuron <- dec$Neurons[di]
glia   <- dec$Astrocytes[di] + dec$Oligodendrocytes[di]
stromal <- dec$Endothelial[di] + dec$Fibroblasts[di]
mes <- ith$is_MES[match(samp, ith$Idat)]
mes[is.na(mes)] <- FALSE
message(sprintf("  Aligned: %d samples (%d MES, %d non-MES)", length(samp), sum(mes), sum(!mes)))

# Design matrices + projectors
X3 <- cbind(1, immune, neuron, glia)
X4 <- cbind(1, immune, neuron, glia, stromal)
P3 <- solve(t(X3) %*% X3) %*% t(X3)   # 4 x n
P4 <- solve(t(X4) %*% X4) %*% t(X4)   # 5 x n

n_probes <- nrow(beta); chunk <- 30000; nch <- ceiling(n_probes / chunk)
message(sprintf("Purifying %d CpGs in %d chunks (3-comp + 4-comp)...", n_probes, nch))

# Accumulate per-sample abs diff (3-comp vs 4-comp), counting non-NA
sum_abs <- numeric(length(samp)); max_abs <- numeric(length(samp)); cnt <- numeric(length(samp))

for (ch in seq_len(nch)) {
  i0 <- (ch - 1) * chunk + 1; i1 <- min(ch * chunk, n_probes); rows <- i0:i1
  bc <- beta[rows, , drop = FALSE]
  na <- is.na(bc); bc[na] <- 0.5
  Y <- myasin(bc)

  b3 <- P3 %*% t(Y)                                   # 4 x m
  e3 <- outer(b3[2, ], immune) + outer(b3[3, ], neuron) + outer(b3[4, ], glia)
  p3 <- (sin(Y - e3) + 1) / 2; p3[p3 < 0] <- 0; p3[p3 > 1] <- 1

  b4 <- P4 %*% t(Y)                                   # 5 x m
  e4 <- outer(b4[2, ], immune) + outer(b4[3, ], neuron) + outer(b4[4, ], glia) + outer(b4[5, ], stromal)
  p4 <- (sin(Y - e4) + 1) / 2; p4[p4 < 0] <- 0; p4[p4 > 1] <- 1

  d <- abs(p3 - p4); d[na] <- NA
  sum_abs <- sum_abs + colSums(d, na.rm = TRUE)
  max_abs <- pmax(max_abs, apply(d, 2, function(x) max(x, na.rm = TRUE)))
  cnt <- cnt + colSums(!na)

  p3[na] <- NA
  beta[rows, ] <- p3                                  # overwrite raw -> becomes 3-comp purified
  if (ch %% 5 == 0) message(sprintf("  chunk %d/%d", ch, nch))
  rm(bc, Y, b3, e3, p3, b4, e4, p4, d, na); gc(verbose = FALSE)
}

mean_abs <- sum_abs / cnt
cat("\n================ 3-component vs 4-component purified beta ================\n")
report <- function(lbl, idx) cat(sprintf("  %-10s n=%3d | mean|Δ|=%.4f  median|Δ|=%.4f  max|Δ|=%.4f\n",
  lbl, sum(idx), mean(mean_abs[idx]), median(mean_abs[idx]), max(max_abs[idx])))
report("ALL", rep(TRUE, length(samp)))
report("MES", mes)
report("non-MES", !mes)
cat("\nInterpretation: mean|Δ| = mean absolute difference in purified beta per sample\n")
cat("(0.00 = adding stromal changes nothing; larger = stromal matters, esp. if MES > non-MES)\n")

# Persist the sensitivity comparison. It answers the obvious reviewer question
# ("why was stroma not regressed out too?") with numbers, and previously existed
# only as console output that was gone once the run finished.
sens <- do.call(rbind, lapply(
  list(ALL = rep(TRUE, length(samp)), MES = mes, `non-MES` = !mes),
  function(idx) data.frame(n = sum(idx),
                           mean_abs_delta   = mean(mean_abs[idx]),
                           median_abs_delta = median(mean_abs[idx]),
                           max_abs_delta    = max(max_abs[idx]))))
sens <- cbind(group = rownames(sens), sens)
write.csv(sens, "data/processed/purification_3comp_vs_4comp.csv", row.names = FALSE)
cat("\nSaved: data/processed/purification_3comp_vs_4comp.csv\n")

message("\nSaving primary 3-component purified matrix...")
saveRDS(beta, "/mnt/volume-hel1-1/siewert/betamatrix_purified_decode3.rds")
message("Done: betamatrix_purified_decode3.rds")

# ============================================================================
# MethylSet export for the Heidelberg classifier
# ============================================================================
# The classifier (mnp.v12epicv2) takes a minfi MethylSet, not a beta matrix, so
# the purified betas are written back as synthetic Meth/Unmeth intensities that
# sum to TOTAL. Split by array version because v1 and v2 need different
# annotations. This block reproduces the June 2026 export that produced
# methylset_purified_decode3_v1.rds / _v2.rds, the files Christin Siewert
# classified; it previously lived only as a loose helper on the data volume and
# was therefore missing from the public repository.
# Existing files are not overwritten unless OVERWRITE is set to TRUE.

suppressPackageStartupMessages(library(minfi))

OVERWRITE <- FALSE
TOTAL <- 10000

sentrix <- as.numeric(sub("_.*", "", colnames(beta)))
is_v2 <- sentrix >= 207000000000   # EPIC v2 sentrix IDs start at 207...

message("\nCreating MethylSet objects...")
for (arr in c("v1", "v2")) {
  outpath <- sprintf("/mnt/volume-hel1-1/siewert/methylset_purified_decode3_%s.rds", arr)
  if (file.exists(outpath) && !OVERWRITE) {
    message(sprintf("  %s: exists, skipping (set OVERWRITE <- TRUE to rebuild)", basename(outpath)))
    next
  }
  idx <- if (arr == "v1") !is_v2 else is_v2
  b <- beta[, idx]
  message(sprintf("  %s: %d samples", arr, ncol(b)))

  ann <- if (arr == "v1") {
    c(array = "IlluminaHumanMethylationEPIC",   annotation = "ilm10b4.hg19")
  } else {
    c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  }

  # MethylSet cannot carry NAs, but the purified matrix has plenty: the raw NA
  # pattern is restored after the regression (p3[na] <- NA above), so 4.7% of all
  # values are missing, every sample is affected, and one carries 416k NAs. They
  # are therefore handed to the classifier as 0.5, i.e. as fully intermediate
  # methylation. This is a real modelling choice, not a formality -- see Methods.
  b_clean <- b
  b_clean[is.na(b_clean)] <- 0.5

  mset <- MethylSet(
    Meth = as.matrix(round(b_clean * TOTAL)),
    Unmeth = as.matrix(round((1 - b_clean) * TOTAL)),
    annotation = ann,
    colData = DataFrame(data.frame(Sample_Name = colnames(b), row.names = colnames(b)))
  )
  saveRDS(mset, outpath)
  message(sprintf("  Saved: %s (%.0f MB)", outpath, file.size(outpath) / 1e6))
  rm(b, b_clean, mset); gc(verbose = FALSE)
}
