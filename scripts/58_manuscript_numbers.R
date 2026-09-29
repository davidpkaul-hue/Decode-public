#!/usr/bin/env Rscript
################################################################################
# Every number the Results text quotes about purification and reclassification
################################################################################
# These were computed ad hoc until now, which meant that switching purification
# variants left no way to regenerate them. This script derives all of them from
# the configured purification in one place, so the manuscript text can be
# checked against a single run.
#
# Sections match the Results paragraphs they feed.
################################################################################

suppressPackageStartupMessages({ library(matrixStats) })
source("scripts/_purification_config.R")

meta <- read.csv("data/processed/siewert_cohort_results.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)
meta <- meta[meta$is_GBM %in% c(TRUE, "TRUE"), ]
lab <- function(x) ifelse(grepl("RTK1", x), "RTK I",
                  ifelse(grepl("RTK2", x), "RTK II",
                  ifelse(grepl("mesenchymal|MES_", x, ignore.case = TRUE), "MES",
                  ifelse(grepl("primitive neuronal|PNC", x, ignore.case = TRUE), "PNC", "other"))))
meta$orig <- lab(meta$Max_subclass)
cal <- read.csv(PURIFIED_CALLS, stringsAsFactors = FALSE)
d <- merge(meta[, c("Idat","orig","Max_subclass_score")], cal, by.x = "Idat", by.y = "Sample")
d$pur <- d$subtype_pur

cat(sprintf("\n=== Cohort: %d GBM, %d carrying one of the four subtype labels ===\n",
            nrow(d), sum(d$orig != "other")))
cat("\n--- Cross-tabulation original vs purified ---\n")
print(table(original = d$orig, purified = d$pur))

cat("\n--- Reclassification per original subtype ---\n")
for (s in c("MES","RTK I","RTK II","PNC")) {
  a <- d[d$orig == s, ]; if (!nrow(a)) next
  ch <- sum(a$pur != s)
  cat(sprintf("  %-7s n=%3d | stable %3d (%.1f%%) | reclassified %3d (%.1f%%)\n",
              s, nrow(a), nrow(a)-ch, 100*(nrow(a)-ch)/nrow(a), ch, 100*ch/nrow(a)))
}
mes <- d[d$orig == "MES", ]
cat(sprintf("\n  MES destinations: RTK II %d | RTK I %d | stays MES %d | other %d | PNC %d\n",
    sum(mes$pur=="RTK II"), sum(mes$pur=="RTK I"), sum(mes$pur=="MES"),
    sum(mes$pur=="other"), sum(mes$pur=="PNC")))
rec <- mes[mes$pur != "MES", ]
cat(sprintf("  of the %d reclassified MES, %d (%.0f%%) go to RTK II\n",
    nrow(rec), sum(rec$pur=="RTK II"), 100*mean(rec$pur=="RTK II")))

# --- calibration threshold subset (0.84 is the classifier's own cutoff) ------
hi <- mes[!is.na(mes$Max_subclass_score) & mes$Max_subclass_score >= 0.84, ]
cat(sprintf("\n  MES above the 0.84 calibration threshold: %d, of which %d (%.1f%%) reclassified\n",
    nrow(hi), sum(hi$pur != "MES"), 100*mean(hi$pur != "MES")))

# --- low-confidence non-GBM calls -------------------------------------------
lo <- d[d$pur == "other", ]
cat(sprintf("\n--- Low-confidence non-GBM calls after purification: %d ---\n", nrow(lo)))
cat(sprintf("  among the %d subtyped samples: %d\n", sum(d$orig != "other"), sum(lo$orig != "other")))
print(table(lo$pred_sub))
cat(sprintf("  score range %.2f - %.2f (median %.2f)\n",
    min(lo$score_sub), max(lo$score_sub), median(lo$score_sub)))
cat("  origin: "); print(table(lo$orig))

# --- score comparison --------------------------------------------------------
cat("\n--- Subclass score before vs after purification (paired Wilcoxon) ---\n")
for (s in c("MES","RTK I","RTK II","PNC")) {
  a <- d[d$orig == s & d$pur != "other" & !is.na(d$Max_subclass_score), ]
  if (nrow(a) < 3) next
  p <- suppressWarnings(wilcox.test(a$Max_subclass_score, a$score_sub, paired = TRUE)$p.value)
  cat(sprintf("  %-7s n=%3d | median %.3f -> %.3f | p = %.3g\n",
      s, nrow(a), median(a$Max_subclass_score), median(a$score_sub), p))
}
cat(sprintf("\n  calls below 0.5 after purification: %d\n", sum(d$score_sub < 0.5)))

# --- failed probes explain the low-confidence calls --------------------------
cat("\n--- Failed probes (NA fraction in the ORIGINAL matrix) ---\n")
b <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
ids <- intersect(d$Idat, colnames(b))
naf <- colMeans(is.na(b[, ids, drop = FALSE]))
rm(b); gc(verbose = FALSE)
d$na_frac <- naf[match(d$Idat, names(naf))]
g1 <- d$na_frac[d$pur == "other"]; g2 <- d$na_frac[d$pur != "other"]
w  <- suppressWarnings(wilcox.test(g1, g2))
auc <- w$statistic / (length(g1)*length(g2)); auc <- max(auc, 1-auc)
cat(sprintf("  low-confidence: median %.1f%% failed | GBM calls: median %.1f%%\n",
            100*median(g1), 100*median(g2)))
cat(sprintf("  Wilcoxon p = %.2g | AUC %.2f\n", w$p.value, auc))
cat(sprintf("  Spearman rho(NA fraction, purified score) = %.2f (p = %.2g)\n",
    cor(d$na_frac, d$score_sub, method = "spearman", use = "complete.obs"),
    suppressWarnings(cor.test(d$na_frac, d$score_sub, method = "spearman"))$p.value))
cat(sprintf("  highest NA fraction in the cohort: %.1f%% (%s, purified call %s, score %.2f)\n",
    100*max(d$na_frac, na.rm = TRUE), d$Idat[which.max(d$na_frac)],
    d$pur[which.max(d$na_frac)], d$score_sub[which.max(d$na_frac)]))

write.csv(d[, c("Idat","orig","pur","Max_subclass_score","score_sub","pred_sub","na_frac")],
          "data/processed/manuscript_numbers_table.csv", row.names = FALSE)
message("\nWrote data/processed/manuscript_numbers_table.csv")
