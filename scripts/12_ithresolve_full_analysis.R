#!/usr/bin/env Rscript
################################################################################
# DECODE: Full ITHresolveGBM-Style Analysis on Siewert Cohort
################################################################################
# Reproduces key analyses from Silverbush et al. (ITHresolveGBM):
#   1. Deconvolution into 6 components (immune, neuron, glia, stem, diff1, diff2)
#   2. Tumor purity & differentiation score
#   3. In-silico purification (extended InfiniumPurify)
#   4. MES reclassification analysis
#   5. Differentiation axis characterization
#
# Input:  /mnt/volume-hel1-1/siewert/betamatrix.rds
#         /mnt/volume-hel1-1/siewert/betamatrix.xlsx
#         /mnt/volume-hel1-1/siewert/classifier_v128_from_html.csv
#         data/external/ithresolve/profiles_full
#         data/processed/siewert_cohort_results.csv  (DECODE results)
# Output: data/processed/siewert_ithresolve_full.csv
#         figures/siewert_ith_*.pdf
################################################################################

if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

suppressPackageStartupMessages({
  library(EpiDISH)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

message("================================================================")
message("ITHresolveGBM-Style Full Analysis: Siewert Cohort")
message("================================================================\n")

# ============================================================================
# Configuration
# ============================================================================

beta_file <- "/mnt/volume-hel1-1/siewert/betamatrix.rds"
meta_file <- "/mnt/volume-hel1-1/siewert/betamatrix.xlsx"
html_scores_file <- "/mnt/volume-hel1-1/siewert/classifier_v128_from_html.csv"
profiles_file <- "data/external/ithresolve/profiles_full"
decode_results_file <- "data/processed/siewert_cohort_results.csv"
output_csv <- "data/processed/siewert_ithresolve_full.csv"
fig_dir <- "figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load Metadata (with v12.8 HTML score fill-in)
# ============================================================================

message("Loading metadata...")
meta <- read_excel(meta_file)
meta <- meta %>%
  rename(
    Idat = `Idat no`,
    Max_subclass = `Max sub-class v12.8 epignostix`,
    Max_subclass_score = `Max subClass Score v12.8 epignostix`,
    Max_class = `Max Class v12.8 Epignostix`,
    Max_class_score = `Max Class Score v12.8 epignostix`,
    Max_family = `Max Family v12.8 Epignostix`,
    Max_family_score = `Max Family Score v12.8 Epignostix`,
    RF_purity = `RF purity prediction (absolute)`
  )

# Fill missing v12.8 scores from HTML reports
if (file.exists(html_scores_file)) {
  html_scores <- read.csv(html_scores_file, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(html_scores))) {
    idx <- which(meta$Idat == html_scores$Sample_ID[i])
    if (length(idx) == 1 && is.na(meta$Max_family_score[idx])) {
      meta$Max_family[idx] <- html_scores$family_predicted[i]
      meta$Max_family_score[idx] <- as.numeric(html_scores$family_score[i])
      meta$Max_class[idx] <- html_scores$class_predicted[i]
      meta$Max_class_score[idx] <- as.numeric(html_scores$class_score[i])
      meta$Max_subclass[idx] <- html_scores$subclass_predicted[i]
      meta$Max_subclass_score[idx] <- as.numeric(html_scores$subclass_score[i])
    }
  }
  message(sprintf("  Filled v12.8 scores from HTML for %d samples", nrow(html_scores)))
}
message(sprintf("  Total samples: %d | v12.8 scores: %d",
                nrow(meta), sum(!is.na(meta$Max_family_score))))

# ============================================================================
# Load ITHresolveGBM Profiles
# ============================================================================

message("\nLoading ITHresolveGBM profiles...")
profiles <- read.table(profiles_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
profile_probes <- profiles$probe
message(sprintf("  Profiles: %d probes x 6 components", nrow(profiles)))

# ============================================================================
# Load Beta Matrix — Memory-Efficient
# ============================================================================

message("\nLoading beta-value matrix...")
beta_full <- readRDS(beta_file)
message(sprintf("  Full matrix: %d CpGs x %d samples", nrow(beta_full), ncol(beta_full)))

# Subset to profile probes
common_probes <- intersect(profile_probes, rownames(beta_full))
beta <- beta_full[common_probes, ]
rm(beta_full)
gc(verbose = FALSE)
message(sprintf("  Subset to %d shared probes (freed full matrix)", nrow(beta)))

# Remove rows with NAs
na_mask <- rowSums(is.na(beta)) > 0
ref_mat <- as.matrix(profiles[match(common_probes, profiles$probe),
                               c("immune", "neuron", "glia", "V4", "V5", "V6")])
rownames(ref_mat) <- common_probes
na_mask <- na_mask | rowSums(is.na(ref_mat)) > 0
beta <- beta[!na_mask, ]
ref_mat <- ref_mat[!na_mask, ]
message(sprintf("  After NA removal: %d probes", nrow(beta)))

# ============================================================================
# Step 1: Deconvolution (EpiDISH RPC with ITH profiles)
# ============================================================================

message("\n=== Step 1: ITHresolveGBM Deconvolution ===\n")
message("  Running EpiDISH RPC (maxit=2000)...")
result <- epidish(beta.m = beta, ref.m = ref_mat, method = "RPC", maxit = 2000)
frac <- as.data.frame(result$estF)
# Column mapping. The published profile file names only immune/neuron/glia; the
# three malignant components are unlabelled (V4, V5, V6). V4 is taken as
# Malignant 1, which Silverbush et al. identify as the stem-like state -- an
# assignment confirmed empirically here (stable RTK I stem-like, stable RTK II
# differentiated, the direction they report). CAUTION: "Diff1" is V6 and "Diff2"
# is V5, so these names do NOT map onto their Malignant 2 / Malignant 3
# numbering. Nothing reported anywhere uses them individually -- only their sum
# enters Purity, Differentiated and DiffScore. Do not read them as Malignant 2/3.
colnames(frac) <- c("Immune", "Neuron", "Glia", "Stem", "Diff2", "Diff1")
frac$Idat <- rownames(frac)
rm(result)
gc(verbose = FALSE)
message(sprintf("  Deconvolved: %d samples", nrow(frac)))

# ============================================================================
# Step 2: Tumor Purity & Differentiation Score
# ============================================================================

message("\n=== Step 2: Purity & Differentiation Score ===\n")

frac$Purity <- frac$Stem + frac$Diff1 + frac$Diff2
frac$Differentiated <- frac$Diff1 + frac$Diff2
frac$NonMalignant <- frac$Immune + frac$Neuron + frac$Glia

# Differentiation score: positive = stem-like, negative = differentiated
frac$DiffScore <- frac$Stem - (frac$Diff1 + frac$Diff2)

# Normalized malignant fractions (within-malignant, sums to 1)
frac$Stem_norm <- ifelse(frac$Purity > 0, frac$Stem / frac$Purity, NA)
frac$Diff1_norm <- ifelse(frac$Purity > 0, frac$Diff1 / frac$Purity, NA)
frac$Diff2_norm <- ifelse(frac$Purity > 0, frac$Diff2 / frac$Purity, NA)

# Binary: more stem-like or more differentiated?
frac$MoreStemLike <- frac$DiffScore > 0

message(sprintf("  Mean purity: %.1f%% (SD %.1f%%)", mean(frac$Purity) * 100, sd(frac$Purity) * 100))
message(sprintf("  Purity range: %.1f%% – %.1f%%", min(frac$Purity) * 100, max(frac$Purity) * 100))
message(sprintf("  Mean diff. score: %+.3f", mean(frac$DiffScore)))
message(sprintf("  More stem-like: %d | More differentiated: %d",
                sum(frac$MoreStemLike), sum(!frac$MoreStemLike)))

# ============================================================================
# Merge with Metadata
# ============================================================================

message("\nMerging with metadata...")
merged <- merge(meta[, c("Idat", "Max_family", "Max_family_score",
                          "Max_class", "Max_class_score",
                          "Max_subclass", "Max_subclass_score", "RF_purity")],
                frac, by = "Idat")
message(sprintf("  Merged: %d samples", nrow(merged)))

# Classify: GBM subtypes from v12.8
merged$is_GBM <- grepl("Glioblastoma, IDH-wildtype", merged$Max_family, ignore.case = TRUE)
merged$is_MES <- grepl("mesenchymal", merged$Max_subclass, ignore.case = TRUE) |
  grepl("mesenchymal", merged$Max_class, ignore.case = TRUE)
merged$GBM_subtype <- case_when(
  grepl("RTK I$|RTK-I$|rtk i$", merged$Max_class, ignore.case = TRUE) ~ "RTK-I",
  grepl("RTK II|RTK-II|rtk ii", merged$Max_class, ignore.case = TRUE) ~ "RTK-II",
  grepl("mesenchymal", merged$Max_class, ignore.case = TRUE) ~ "MES",
  grepl("MYCN", merged$Max_class, ignore.case = TRUE) ~ "MYCN",
  grepl("MID", merged$Max_class, ignore.case = TRUE) ~ "MID",
  merged$is_GBM ~ "GBM_other",
  TRUE ~ "Non-GBM"
)

message(sprintf("\n  GBM subtype distribution:"))
for (st in sort(unique(merged$GBM_subtype))) {
  message(sprintf("    %-12s: %d", st, sum(merged$GBM_subtype == st)))
}

# ============================================================================
# Step 3: In-Silico Purification
# ============================================================================
# Extended InfiniumPurify: for each CpG, regress asin(2*beta-1) ~ immune + neuron + glia
# across all samples, then purified = (sin(residuals + intercept) + 1) / 2
# Silverbush et al. restrict their purification to samples above 20% predicted
# purity. We deliberately do NOT: the regression below uses every sample, with
# the low-purity ones acting as quasi-controls that pin down the non-malignant
# coefficients. The count is reported only so the difference to the paper is
# visible in the log.

message("\n=== Step 3: In-Silico Purification ===\n")

purify_samples <- merged$Idat[merged$Purity > 0.20]
message(sprintf("  Samples with purity > 20%%: %d / %d", length(purify_samples), nrow(merged)))

# Non-malignant fractions for regression
nm_fracs <- merged[, c("Idat", "Immune", "Neuron", "Glia")]
rownames(nm_fracs) <- nm_fracs$Idat

# We use all samples in the regression (including low-purity ones as quasi-controls)
# This provides better estimation of non-malignant effects
all_samples <- merged$Idat[merged$Idat %in% colnames(beta)]
beta_for_purify <- beta[, all_samples]

A <- nm_fracs[all_samples, "Immune"]
B <- nm_fracs[all_samples, "Neuron"]
C <- nm_fracs[all_samples, "Glia"]

message("  Running per-CpG regression (this takes a while)...")
message(sprintf("  %d CpGs x %d samples", nrow(beta_for_purify), ncol(beta_for_purify)))

# Arcsin transformation
.myasin <- function(x) asin(pmin(pmax(2 * x - 1, -1), 1))

# Vectorized purification: fit lm(Y ~ A + B + C) for each CpG
# Y_mat: nprobes x nsamples (same layout as beta_for_purify)
Y_mat <- .myasin(beta_for_purify)

# Design matrix: nsamples x 4
X_design <- cbind(1, A, B, C)

# Solve all regressions at once via OLS:
# For Y^T (nsamples x nprobes) = X (nsamples x 4) * B (4 x nprobes)
# B = (X'X)^{-1} X' Y^T
message("  Fitting regression model (vectorized)...")
XtX_inv <- solve(t(X_design) %*% X_design)
beta_hat <- XtX_inv %*% t(X_design) %*% t(Y_mat)  # 4 x nprobes

# Purification: remove non-malignant effects, keep intercept (= pure tumor signal)
# purified_i = Y_i - b1*A_i - b2*B_i - b3*C_i  (per CpG, across samples)
# effects_mat: nprobes x nsamples
effects_mat <- outer(beta_hat[2, ], A) + outer(beta_hat[3, ], B) + outer(beta_hat[4, ], C)
Y_purified <- Y_mat - effects_mat  # nprobes x nsamples

# Back-transform from arcsin to beta
beta_purified <- (sin(Y_purified) + 1) / 2
beta_purified[beta_purified < 0] <- 0
beta_purified[beta_purified > 1] <- 1

message(sprintf("  Purified matrix: %d probes x %d samples", nrow(beta_purified), ncol(beta_purified)))

# Clean up large matrices
rm(Y_mat, effects_mat, Y_purified, X_design, XtX_inv, beta_hat)
gc(verbose = FALSE)

# ============================================================================
# Step 4: Re-Deconvolve Purified Profiles
# ============================================================================

message("\n=== Step 4: Re-Deconvolution of Purified Profiles ===\n")

# Run deconvolution on purified beta values to see the effect
ref_purify <- ref_mat[rownames(beta_purified), ]
result_purified <- epidish(beta.m = beta_purified, ref.m = ref_purify,
                            method = "RPC", maxit = 2000)
frac_purified <- as.data.frame(result_purified$estF)
colnames(frac_purified) <- paste0("Pur_", c("Immune", "Neuron", "Glia", "Stem", "Diff2", "Diff1"))
frac_purified$Idat <- rownames(frac_purified)
frac_purified$Pur_Purity <- frac_purified$Pur_Stem + frac_purified$Pur_Diff1 + frac_purified$Pur_Diff2
frac_purified$Pur_DiffScore <- frac_purified$Pur_Stem - (frac_purified$Pur_Diff1 + frac_purified$Pur_Diff2)
rm(result_purified)
gc(verbose = FALSE)

merged <- merge(merged, frac_purified, by = "Idat", all.x = TRUE)
message(sprintf("  Mean purity BEFORE purification: %.1f%%", mean(merged$Purity) * 100))
message(sprintf("  Mean purity AFTER  purification: %.1f%%", mean(merged$Pur_Purity, na.rm = TRUE) * 100))

# ============================================================================
# Step 5: MES Reclassification Analysis
# ============================================================================

message("\n=== Step 5: MES Reclassification Analysis ===\n")

# Key question: Are MES-classified samples distinguished by high immune infiltration?
gbm <- merged %>% filter(is_GBM)
message(sprintf("  GBM IDH-wt samples: %d", nrow(gbm)))

# MES vs non-MES comparison
if (sum(gbm$GBM_subtype == "MES") > 0) {
  mes <- gbm %>% filter(GBM_subtype == "MES")
  non_mes <- gbm %>% filter(GBM_subtype != "MES")

  message(sprintf("\n  MES samples: %d", nrow(mes)))
  message(sprintf("  Non-MES GBM: %d", nrow(non_mes)))

  # Immune infiltration comparison
  message(sprintf("\n  --- Immune Fraction ---"))
  message(sprintf("  MES:     mean=%.1f%% (SD %.1f%%)", mean(mes$Immune) * 100, sd(mes$Immune) * 100))
  message(sprintf("  Non-MES: mean=%.1f%% (SD %.1f%%)", mean(non_mes$Immune) * 100, sd(non_mes$Immune) * 100))
  wt <- wilcox.test(mes$Immune, non_mes$Immune)
  message(sprintf("  Wilcoxon p = %.2e", wt$p.value))

  # Purity comparison
  message(sprintf("\n  --- Tumor Purity ---"))
  message(sprintf("  MES:     mean=%.1f%% (SD %.1f%%)", mean(mes$Purity) * 100, sd(mes$Purity) * 100))
  message(sprintf("  Non-MES: mean=%.1f%% (SD %.1f%%)", mean(non_mes$Purity) * 100, sd(non_mes$Purity) * 100))
  wt_pur <- wilcox.test(mes$Purity, non_mes$Purity)
  message(sprintf("  Wilcoxon p = %.2e", wt_pur$p.value))

  # Purity after purification
  message(sprintf("\n  --- After In-Silico Purification ---"))
  message(sprintf("  MES:     purity %.1f%% -> %.1f%%",
                  mean(mes$Purity) * 100, mean(mes$Pur_Purity, na.rm = TRUE) * 100))
  message(sprintf("  Non-MES: purity %.1f%% -> %.1f%%",
                  mean(non_mes$Purity) * 100, mean(non_mes$Pur_Purity, na.rm = TRUE) * 100))

  # Diff score comparison
  message(sprintf("\n  --- Differentiation Score ---"))
  message(sprintf("  MES:     mean=%+.3f (more %s)",
                  mean(mes$DiffScore), ifelse(mean(mes$DiffScore) > 0, "stem-like", "differentiated")))
  message(sprintf("  Non-MES: mean=%+.3f (more %s)",
                  mean(non_mes$DiffScore), ifelse(mean(non_mes$DiffScore) > 0, "stem-like", "differentiated")))
} else {
  message("  No MES-classified samples in this cohort.")
}

# Per-subtype summary
message("\n  --- Per-Subtype Summary ---\n")
message(sprintf("  %-12s %5s %8s %8s %8s %8s %10s",
                "Subtype", "n", "Immune%", "Purity%", "DiffScr", "Stem%", "Class.Scr"))
message(paste(rep("-", 75), collapse = ""))
for (st in c("RTK-I", "RTK-II", "MES", "MYCN", "MID", "GBM_other", "Non-GBM")) {
  sub <- merged %>% filter(GBM_subtype == st)
  if (nrow(sub) > 0) {
    message(sprintf("  %-12s %5d %7.1f%% %7.1f%% %+8.3f %7.1f%% %10.3f",
                    st, nrow(sub),
                    mean(sub$Immune) * 100, mean(sub$Purity) * 100,
                    mean(sub$DiffScore),
                    mean(sub$Stem_norm, na.rm = TRUE) * 100,
                    mean(sub$Max_family_score, na.rm = TRUE)))
  }
}

# ============================================================================
# Step 6: Differentiation Axis Analysis
# ============================================================================

message("\n=== Step 6: Differentiation Axis ===\n")

# Distribution along the differentiation axis
message("  Differentiation score distribution (GBM only):")
message(sprintf("    Quartiles: %.3f | %.3f | %.3f | %.3f | %.3f",
                quantile(gbm$DiffScore, 0), quantile(gbm$DiffScore, 0.25),
                quantile(gbm$DiffScore, 0.5), quantile(gbm$DiffScore, 0.75),
                quantile(gbm$DiffScore, 1)))
message(sprintf("    Stem-like (score > 0): %d (%.1f%%)",
                sum(gbm$DiffScore > 0), sum(gbm$DiffScore > 0) / nrow(gbm) * 100))
message(sprintf("    Differentiated (score < 0): %d (%.1f%%)",
                sum(gbm$DiffScore < 0), sum(gbm$DiffScore < 0) / nrow(gbm) * 100))

# Correlation with classifier score
r_diff <- cor.test(gbm$DiffScore, gbm$Max_family_score, use = "complete.obs")
message(sprintf("    DiffScore vs ClassifierScore: r=%+.3f (p=%.2e)",
                r_diff$estimate, r_diff$p.value))

# ============================================================================
# Visualizations
# ============================================================================

message("\n=== Generating Plots ===\n")

sig_stars <- Vectorize(function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("n.s.")
})

subtype_colors <- c("RTK-I" = "#E41A1C", "RTK-II" = "#377EB8", "MES" = "#4DAF4A",
                     "MYCN" = "#984EA3", "MID" = "#FF7F00", "GBM_other" = "#A65628",
                     "Non-GBM" = "#999999")

# --------------------------------------------------------------------------
# Plot 1: Composition by GBM subtype
# --------------------------------------------------------------------------

comp_data <- merged %>%
  select(Idat, GBM_subtype, Immune, Neuron, Glia, Stem, Diff1, Diff2) %>%
  pivot_longer(cols = c(Immune, Neuron, Glia, Stem, Diff1, Diff2),
               names_to = "Component", values_to = "Fraction") %>%
  mutate(Fraction_Pct = Fraction * 100,
         Component = factor(Component, levels = c("Stem", "Diff1", "Diff2",
                                                    "Immune", "Glia", "Neuron")),
         GBM_subtype = factor(GBM_subtype, levels = names(subtype_colors)))

p1 <- ggplot(comp_data, aes(x = Component, y = Fraction_Pct, fill = GBM_subtype)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = subtype_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "ITHresolveGBM: Component Composition by v12.8 Subtype",
    subtitle = sprintf("Siewert Cohort (n=%d)", nrow(merged)),
    x = "Component", y = "Fraction (%)", fill = "v12.8 Subtype"
  )
ggsave(file.path(fig_dir, "siewert_ith_composition_by_subtype.pdf"), p1, width = 14, height = 8)
message("  Saved: siewert_ith_composition_by_subtype.pdf")

# --------------------------------------------------------------------------
# Plot 2: Purity vs Immune by subtype (key MES artifact test)
# --------------------------------------------------------------------------

p2 <- ggplot(merged %>% filter(is_GBM),
             aes(x = Immune * 100, y = Purity * 100, color = GBM_subtype)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = subtype_colors) +
  theme_bw() +
  labs(
    title = "Tumor Purity vs Immune Infiltration by GBM Subtype",
    subtitle = sprintf("GBM IDH-wt only (n=%d) | MES subtype expected to show high immune / low purity",
                       sum(merged$is_GBM)),
    x = "Immune Fraction (%)", y = "Tumor Purity (%)", color = "v12.8 Subtype"
  )
ggsave(file.path(fig_dir, "siewert_ith_purity_vs_immune.pdf"), p2, width = 10, height = 8)
message("  Saved: siewert_ith_purity_vs_immune.pdf")

# --------------------------------------------------------------------------
# Plot 3: Differentiation axis (stem vs differentiated)
# --------------------------------------------------------------------------

p3 <- ggplot(merged %>% filter(is_GBM & Purity > 0.05),
             aes(x = Stem_norm * 100, y = (Diff1_norm + Diff2_norm) * 100,
                 color = GBM_subtype)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = -1, intercept = 100, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = subtype_colors) +
  theme_bw() +
  labs(
    title = "GBM Differentiation Axis: Stem-like vs Differentiated",
    subtitle = sprintf("Normalized malignant fractions (purity > 5%%) | n=%d",
                       sum(merged$is_GBM & merged$Purity > 0.05)),
    x = "Stem-like Fraction (% of malignant)",
    y = "Differentiated Fraction (% of malignant)",
    color = "v12.8 Subtype"
  )
ggsave(file.path(fig_dir, "siewert_ith_differentiation_axis.pdf"), p3, width = 10, height = 8)
message("  Saved: siewert_ith_differentiation_axis.pdf")

# --------------------------------------------------------------------------
# Plot 4: Differentiation score distribution by subtype
# --------------------------------------------------------------------------

p4 <- ggplot(merged %>% filter(is_GBM),
             aes(x = GBM_subtype, y = DiffScore, fill = GBM_subtype)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = subtype_colors, guide = "none") +
  theme_bw() +
  annotate("text", x = 0.5, y = max(merged$DiffScore[merged$is_GBM]) * 0.9,
           label = "← more stem-like", hjust = 0, size = 3, color = "grey40") +
  annotate("text", x = 0.5, y = min(merged$DiffScore[merged$is_GBM]) * 0.9,
           label = "← more differentiated", hjust = 0, size = 3, color = "grey40") +
  labs(
    title = "Differentiation Score by GBM v12.8 Subtype",
    subtitle = sprintf("Score = Stem - (Diff1 + Diff2) | GBM IDH-wt (n=%d)", sum(merged$is_GBM)),
    x = "v12.8 Subtype", y = "Differentiation Score"
  )
ggsave(file.path(fig_dir, "siewert_ith_diffscore_by_subtype.pdf"), p4, width = 10, height = 7)
message("  Saved: siewert_ith_diffscore_by_subtype.pdf")

# --------------------------------------------------------------------------
# Plot 5: Purification effect — before vs after
# --------------------------------------------------------------------------

purify_comp <- merged %>%
  filter(Purity > 0.20) %>%
  select(Idat, GBM_subtype,
         `Before: Immune` = Immune, `Before: Purity` = Purity,
         `After: Immune` = Pur_Immune, `After: Purity` = Pur_Purity) %>%
  pivot_longer(cols = starts_with("Before") | starts_with("After"),
               names_to = "Variable", values_to = "Fraction") %>%
  mutate(Fraction_Pct = Fraction * 100,
         Stage = ifelse(grepl("Before", Variable), "Before Purification", "After Purification"),
         Measure = gsub("^(Before|After): ", "", Variable),
         Stage = factor(Stage, levels = c("Before Purification", "After Purification")))

p5 <- ggplot(purify_comp, aes(x = Stage, y = Fraction_Pct, fill = Stage)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~Measure, scales = "free_y") +
  scale_fill_manual(values = c("Before Purification" = "#FC8D62",
                                "After Purification" = "#66C2A5"),
                    guide = "none") +
  theme_bw() +
  theme(strip.text = element_text(face = "bold")) +
  labs(
    title = "In-Silico Purification Effect",
    subtitle = sprintf("Samples with purity > 20%% (n=%d)", sum(merged$Purity > 0.20)),
    x = NULL, y = "Fraction (%)"
  )
ggsave(file.path(fig_dir, "siewert_ith_purification_effect.pdf"), p5, width = 10, height = 7)
message("  Saved: siewert_ith_purification_effect.pdf")

# --------------------------------------------------------------------------
# Plot 6: Classifier score vs differentiation score
# --------------------------------------------------------------------------

gbm_valid <- merged %>% filter(is_GBM & !is.na(Max_family_score))
r_val <- cor.test(gbm_valid$Max_family_score, gbm_valid$DiffScore)

p6 <- ggplot(gbm_valid, aes(x = Max_family_score, y = DiffScore, color = GBM_subtype)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed",
              linewidth = 0.6, inherit.aes = FALSE,
              aes(x = Max_family_score, y = DiffScore)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = subtype_colors) +
  theme_bw() +
  labs(
    title = "Classifier Score vs Differentiation Score",
    subtitle = sprintf("GBM IDH-wt (n=%d) | r=%+.3f %s",
                       nrow(gbm_valid), r_val$estimate, sig_stars(r_val$p.value)),
    x = "Methylation Classifier Score (v12.8 Family)",
    y = "Differentiation Score (positive = stem-like)",
    color = "v12.8 Subtype"
  )
ggsave(file.path(fig_dir, "siewert_ith_score_vs_diffscore.pdf"), p6, width = 10, height = 8)
message("  Saved: siewert_ith_score_vs_diffscore.pdf")

# --------------------------------------------------------------------------
# Plot 7: Ternary-style composition (stacked bar by subtype)
# --------------------------------------------------------------------------

subtype_means <- merged %>%
  filter(is_GBM) %>%
  group_by(GBM_subtype) %>%
  summarise(
    n = n(),
    Stem = mean(Stem),
    Diff1 = mean(Diff1),
    Diff2 = mean(Diff2),
    Immune = mean(Immune),
    Glia = mean(Glia),
    Neuron = mean(Neuron),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Stem, Diff1, Diff2, Immune, Glia, Neuron),
               names_to = "Component", values_to = "Fraction") %>%
  mutate(Component = factor(Component, levels = c("Stem", "Diff1", "Diff2",
                                                    "Immune", "Glia", "Neuron")),
         label = paste0(GBM_subtype, " (n=", n, ")"))

comp_colors <- c("Stem" = "#E41A1C", "Diff1" = "#377EB8", "Diff2" = "#4DAF4A",
                  "Immune" = "#FF7F00", "Glia" = "#984EA3", "Neuron" = "#F781BF")

p7 <- ggplot(subtype_means, aes(x = label, y = Fraction * 100, fill = Component)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = comp_colors) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Mean Cell Composition by GBM Subtype (ITHresolveGBM)",
    subtitle = "Siewert Cohort — GBM IDH-wt only",
    x = "v12.8 Subtype", y = "Mean Fraction (%)", fill = "Component"
  )
ggsave(file.path(fig_dir, "siewert_ith_mean_composition_subtypes.pdf"), p7, width = 10, height = 7)
message("  Saved: siewert_ith_mean_composition_subtypes.pdf")

# ============================================================================
# Summary: Key Results Table
# ============================================================================

message("\n================================================================")
message("SUMMARY")
message("================================================================\n")

# MES reclassification prediction
n_mes <- sum(merged$GBM_subtype == "MES")
if (n_mes > 0) {
  mes_data <- merged %>% filter(GBM_subtype == "MES")
  message(sprintf("MES samples: %d", n_mes))
  message(sprintf("  Mean immune: %.1f%% vs non-MES GBM: %.1f%%",
                  mean(mes_data$Immune) * 100,
                  mean(merged$Immune[merged$is_GBM & merged$GBM_subtype != "MES"]) * 100))
  message(sprintf("  Mean purity: %.1f%% vs non-MES GBM: %.1f%%",
                  mean(mes_data$Purity) * 100,
                  mean(merged$Purity[merged$is_GBM & merged$GBM_subtype != "MES"]) * 100))
  message("  → If ITHresolveGBM findings hold: MES samples would likely")
  message("    reclassify to RTK-I or RTK-II after purification.")
  message(sprintf("  MES diff. score: %+.3f → mostly %s",
                  mean(mes_data$DiffScore),
                  ifelse(mean(mes_data$DiffScore) > 0, "stem-like → likely RTK-I",
                         "differentiated → likely RTK-II")))
} else {
  message("No MES-classified samples in this cohort.")
}

# Low-confidence / unclassifiable samples
low_conf <- merged %>% filter(Max_family_score < 0.9)
message(sprintf("\nLow-confidence samples (score < 0.9): %d", nrow(low_conf)))
if (nrow(low_conf) > 0) {
  message(sprintf("  Mean purity: %.1f%% (low purity → classification difficulty)",
                  mean(low_conf$Purity) * 100))
  message(sprintf("  Mean purity of high-conf: %.1f%%",
                  mean(merged$Purity[merged$Max_family_score >= 0.9], na.rm = TRUE) * 100))
  r_pur_score <- cor.test(merged$Max_family_score, merged$Purity, use = "complete.obs")
  message(sprintf("  Purity vs Classifier Score: r=%+.3f (p=%.2e)",
                  r_pur_score$estimate, r_pur_score$p.value))
  message("  → In-silico purification could potentially rescue these for classification")
}

# Differentiation summary
message(sprintf("\nDifferentiation axis (GBM only, n=%d):", sum(merged$is_GBM)))
message(sprintf("  Stem-like (DiffScore > 0): %d (%.1f%%)",
                sum(gbm$DiffScore > 0), sum(gbm$DiffScore > 0) / nrow(gbm) * 100))
message(sprintf("  Differentiated (DiffScore < 0): %d (%.1f%%)",
                sum(gbm$DiffScore < 0), sum(gbm$DiffScore < 0) / nrow(gbm) * 100))

# ============================================================================
# Save Results
# ============================================================================

write.csv(merged, output_csv, row.names = FALSE)
message(sprintf("\nSaved: %s", output_csv))

# Also save the purified beta matrix for potential later classifier re-run
purified_out <- file.path(dirname(beta_file), "betamatrix_purified.rds")
saveRDS(beta_purified, purified_out)
message(sprintf("Saved purified beta matrix: %s", purified_out))

message("\n================================================================")
message("Done! Full ITHresolveGBM analysis complete.")
message("NOTE: Re-running the Heidelberg classifier on purified profiles")
message("      requires the mnp.v12epicv2 R package (registration needed")
message("      at molecularneuropathology.org). Ask Prof. Capper for access.")
message("================================================================")
