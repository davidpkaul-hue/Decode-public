#!/usr/bin/env Rscript
################################################################################
# Sensitivity check: DECODE-fraction vs ITHresolveGBM-fraction purification
################################################################################
# The Discussion argues that the reclassification does not depend on which
# deconvolution supplies the immune/neuron/glia fractions. This quantifies it by
# correlating the two purified matrices:
#   betamatrix_purified_decode3.rds  (script 13b, DECODE fractions)
#   betamatrix_purified_full.rds     (script 13,  ITHresolveGBM fractions)
#
# Both matrices are 2.3 GB, so they are compared on a random subset of CpGs;
# with 50,000 probes x 436 samples the estimate rests on ~21 million values and
# is far more precise than the two reported digits.
#
# Result of the run behind the manuscript (seed 42):
#   Pearson r = 0.987, mean |delta| = 0.0368, per-sample median r = 0.989
#   (range 0.935 to 0.997)
# Output: data/processed/purified_decode_vs_ith.csv
################################################################################

set.seed(42)
N_CPG <- 50000

message("Loading DECODE-purified matrix ...")
a <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix_purified_decode3.rds")
idx <- sort(sample(nrow(a), N_CPG))
cn <- colnames(a); rn <- rownames(a)[idx]
a <- a[idx, ]; gc(verbose = FALSE)

message("Loading ITH-purified matrix ...")
b <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix_purified_full.rds")
b <- b[rn, cn]; gc(verbose = FALSE)

ok <- is.finite(a) & is.finite(b)
r_all  <- cor(a[ok], b[ok])
mad_all <- mean(abs(a[ok] - b[ok]))
per <- sapply(seq_len(ncol(a)), function(j) {
  o <- ok[, j]
  if (sum(o) < 100) NA_real_ else cor(a[o, j], b[o, j])
})

message(sprintf("\nvalues compared : %d", sum(ok)))
message(sprintf("Pearson r       : %.4f", r_all))
message(sprintf("mean |delta|    : %.4f", mad_all))
message(sprintf("per sample r    : median %.4f, range %.4f to %.4f",
                median(per, na.rm = TRUE), min(per, na.rm = TRUE), max(per, na.rm = TRUE)))

write.csv(data.frame(n_cpg = N_CPG, n_values = sum(ok), pearson_r = r_all,
                     mean_abs_delta = mad_all,
                     per_sample_r_median = median(per, na.rm = TRUE),
                     per_sample_r_min = min(per, na.rm = TRUE),
                     per_sample_r_max = max(per, na.rm = TRUE)),
          "data/processed/purified_decode_vs_ith.csv", row.names = FALSE)
message("Saved: data/processed/purified_decode_vs_ith.csv")
