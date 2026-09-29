# ============================================================================
# Script 53: In-silico purification with microglia as a separate component
# ----------------------------------------------------------------------------
# Motivation: the canonical purification (script 13b -> betamatrix_purified_
# decode3.rds) regresses on immune (six blood-derived types), neurons and glia
# (astrocytes + oligodendrocytes). Microglia entered NONE of these terms, even
# though they are the largest single immune population in MES tumors
# (14.5% vs 5.0% in RTK) and are only 4% explained by the three existing
# regressors (R^2 = 0.038).
#
# This script adds microglia as its own regressor (4 components). A separate
# term is preferable to folding microglia into the immune term: microglia and
# blood-derived macrophages have distinct methylomes, and a shared coefficient
# would distort the immune term. Empirically the separate term also changes the
# purified matrix less (mean|delta| 0.0137 vs 0.0229).
#
# Output: full-CpG purified matrix + MethylSets (v1/v2) for the real Heidelberg
# classifier run at the Institute of Neuropathology, Charite.
# ============================================================================

suppressPackageStartupMessages({ library(minfi) })
myasin <- function(x) asin(pmin(pmax(2 * x - 1, -1), 1))

OUT_RDS <- "/mnt/volume-hel1-1/siewert/betamatrix_purified_microglia4.rds"
EXPORT  <- "/mnt/volume-hel1-1/siewert/purified_methylsets_microglia4"
dir.create(EXPORT, showWarnings = FALSE)

message("Loading DECODE fractions ...")
dec <- read.csv("data/processed/siewert_cohort_results.csv", stringsAsFactors = FALSE)

message("Loading raw beta matrix (2.2 GB) ...")
beta <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
samp <- intersect(colnames(beta), dec$Idat)
beta <- beta[, samp]
di <- match(samp, dec$Idat)

immune    <- dec$Immune_total[di]                                  # six blood-derived types
microglia <- dec$Microglia[di]                                     # NEW: own component
neuron    <- dec$Neurons[di]
glia      <- dec$Astrocytes[di] + dec$Oligodendrocytes[di]
stopifnot(!any(is.na(c(immune, microglia, neuron, glia))))
message(sprintf("  %d CpGs x %d samples aligned", nrow(beta), ncol(beta)))

X <- cbind(1, immune, microglia, neuron, glia)
P <- solve(t(X) %*% X) %*% t(X)

n_probes <- nrow(beta); chunk <- 30000; nch <- ceiling(n_probes / chunk)
message(sprintf("Purifying %d CpGs in %d chunks (4 components) ...", n_probes, nch))
for (ch in seq_len(nch)) {
  i0 <- (ch - 1) * chunk + 1; i1 <- min(ch * chunk, n_probes); rows <- i0:i1
  bc <- beta[rows, , drop = FALSE]
  na <- is.na(bc); bc[na] <- 0.5
  Y  <- myasin(bc)
  b  <- P %*% t(Y)
  e  <- outer(b[2, ], immune) + outer(b[3, ], microglia) +
        outer(b[4, ], neuron) + outer(b[5, ], glia)
  p <- (sin(Y - e) + 1) / 2; p[p < 0] <- 0; p[p > 1] <- 1
  p[na] <- NA
  beta[rows, ] <- p
  if (ch %% 5 == 0) message(sprintf("  chunk %d/%d", ch, nch))
}
message("Saving purified matrix ...")
saveRDS(beta, OUT_RDS)
message(sprintf("  %s (%.0f MB)", OUT_RDS, file.size(OUT_RDS) / 1e6))

# --- MethylSet export, identical format to the previous Christin runs --------
sentrix <- as.numeric(sub("_.*", "", colnames(beta)))
is_v2 <- sentrix >= 207000000000
TOTAL <- 10000
for (arr in c("v1", "v2")) {
  idx <- if (arr == "v1") !is_v2 else is_v2
  if (!any(idx)) next
  bb <- beta[, idx, drop = FALSE]
  bb[is.na(bb)] <- 0.5
  ann <- if (arr == "v1") c(array = "IlluminaHumanMethylationEPIC",   annotation = "ilm10b4.hg19")
         else            c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  ms <- MethylSet(Meth = round(bb * TOTAL), Unmeth = round((1 - bb) * TOTAL),
                  annotation = ann,
                  colData = DataFrame(data.frame(Sample_Name = colnames(bb),
                                                 row.names = colnames(bb))))
  fp <- file.path(EXPORT, sprintf("methylset_purified_microglia4_%s.rds", arr))
  saveRDS(ms, fp)
  message(sprintf("  Saved %s: %d samples (%.0f MB)", fp, ncol(bb), file.size(fp) / 1e6))
  rm(bb, ms); gc(verbose = FALSE)
}

man <- data.frame(Sample = samp, immune_blood = immune, microglia = microglia,
                  neuron = neuron, glia = glia,
                  platform = ifelse(is_v2, "EPICv2", "EPICv1"))
write.csv(man, file.path(EXPORT, "manifest_microglia4.csv"), row.names = FALSE)
write.csv(man, "data/processed/purification_microglia4_manifest.csv", row.names = FALSE)
message("\nDONE. Hand the two MethylSets in ", EXPORT, " to Christin for mnp.v12epicv2.")
