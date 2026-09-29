#!/usr/bin/env Rscript
################################################################################
# DECODE vs ITHresolveGBM Comparison
################################################################################
# Purpose: Apply ITHresolveGBM (Silverbush et al.) deconvolution to the Siewert
#          cohort and compare with DECODE results
#
# ITHresolveGBM uses 6 components: immune, neuron, glia and three unlabelled
#          malignant columns V4, V5, V6. V4 = Malignant 1 (stem-like per
#          Silverbush et al.). The Diff1/Diff2 names below are V6/V5 and do NOT
#          correspond to their Malignant 2 / Malignant 3 numbering; only the sum
#          of the two is ever used in a reported quantity.
# DECODE uses 13 cell types with finer non-malignant resolution
#
# Both methods use EpiDISH RPC — the difference is the reference profiles.
#
# Input:  /mnt/volume-hel1-1/siewert/betamatrix.rds
#         data/external/ithresolve/profiles_full
#         data/external/ithresolve/profiles_after_feature_selection_with_sc
#         data/processed/siewert_cohort_results.csv  (DECODE results)
# Output: data/processed/siewert_ithresolve_results.csv
#         figures/siewert_ithresolve_*.pdf
################################################################################

if (file.exists("scripts/R45_patch.R")) {
  source("scripts/R45_patch.R")
}

suppressPackageStartupMessages({
  library(EpiDISH)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

message("==============================================================")
message("DECODE vs ITHresolveGBM: Siewert Cohort Comparison")
message("==============================================================\n")

# ============================================================================
# Configuration
# ============================================================================

beta_file <- "/mnt/volume-hel1-1/siewert/betamatrix.rds"
profiles_full_file <- "data/external/ithresolve/profiles_full"
profiles_fs_file <- "data/external/ithresolve/profiles_after_feature_selection_with_sc"
decode_results_file <- "data/processed/siewert_cohort_results.csv"
output_csv <- "data/processed/siewert_ithresolve_results.csv"
fig_dir <- "figures"

# ============================================================================
# Load ITHresolveGBM Profiles
# ============================================================================

message("Loading ITHresolveGBM reference profiles...")
profiles_full <- read.table(profiles_full_file, header = TRUE, sep = "\t",
                            stringsAsFactors = FALSE)
profiles_fs <- read.table(profiles_fs_file, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
message(sprintf("  Full profiles: %d probes x %d components", nrow(profiles_full), ncol(profiles_full) - 1))
message(sprintf("  Feature-selected profiles: %d probes x %d components", nrow(profiles_fs), ncol(profiles_fs) - 1))

# ============================================================================
# Load Beta Matrix — Memory-Efficient: Only Load Needed Probes
# ============================================================================
# The full beta matrix is 2.2GB (721K x 436). On 8GB RAM, we must subset
# immediately. We load the full matrix, extract needed probes, and free it.

message("\nIdentifying needed probes...")
all_needed_probes <- unique(c(profiles_full$probe, profiles_fs$probe))
message(sprintf("  Total unique probes needed: %d", length(all_needed_probes)))

message("Loading beta-value matrix...")
beta_full <- readRDS(beta_file)
message(sprintf("  Full matrix: %d CpGs x %d samples", nrow(beta_full), ncol(beta_full)))

# Subset to needed probes immediately
common_all <- intersect(all_needed_probes, rownames(beta_full))
message(sprintf("  Probes found in beta matrix: %d", length(common_all)))
beta_subset <- beta_full[common_all, ]
rm(beta_full)
gc(verbose = FALSE)
message(sprintf("  Subset matrix: %d probes x %d samples (freed full matrix)",
                nrow(beta_subset), ncol(beta_subset)))

# ============================================================================
# Run ITHresolveGBM Deconvolution (Full Profiles)
# ============================================================================

message("\n--- ITHresolveGBM with FULL profiles ---")

common_full <- intersect(profiles_full$probe, rownames(beta_subset))
message(sprintf("  Overlapping probes: %d (minimum required: 80,000)", length(common_full)))

# Build reference matrix
ref_full <- as.matrix(profiles_full[match(common_full, profiles_full$probe),
                                     c("immune", "neuron", "glia", "V4", "V5", "V6")])
rownames(ref_full) <- common_full
beta_common_full <- beta_subset[common_full, ]

# Remove rows with NAs
na_rows <- rowSums(is.na(beta_common_full)) > 0 | rowSums(is.na(ref_full)) > 0
message(sprintf("  Removing %d rows with NAs", sum(na_rows)))
ref_full <- ref_full[!na_rows, ]
beta_common_full <- beta_common_full[!na_rows, ]
message(sprintf("  Final probes for deconvolution: %d", nrow(ref_full)))

# Run EpiDISH RPC
message("  Running EpiDISH RPC (maxit=2000)...")
result_full <- epidish(beta.m = beta_common_full, ref.m = ref_full,
                       method = "RPC", maxit = 2000)
fractions_full <- as.data.frame(result_full$estF)
colnames(fractions_full) <- c("ITH_Immune", "ITH_Neuron", "ITH_Glia",
                               "ITH_Stem", "ITH_Diff2", "ITH_Diff1")
fractions_full$Idat <- rownames(fractions_full)
message(sprintf("  Deconvolved: %d samples", nrow(fractions_full)))
rm(ref_full, beta_common_full, result_full)
gc(verbose = FALSE)

# Derived columns
fractions_full$ITH_Purity <- fractions_full$ITH_Stem + fractions_full$ITH_Diff1 +
                              fractions_full$ITH_Diff2
fractions_full$ITH_Differentiated <- fractions_full$ITH_Diff1 + fractions_full$ITH_Diff2
fractions_full$ITH_NonMalignant <- fractions_full$ITH_Immune + fractions_full$ITH_Neuron +
                                    fractions_full$ITH_Glia

# ============================================================================
# Run ITHresolveGBM Deconvolution (Feature-Selected Profiles)
# ============================================================================

message("\n--- ITHresolveGBM with FEATURE-SELECTED profiles ---")

common_fs <- intersect(profiles_fs$probe, rownames(beta_subset))
message(sprintf("  Overlapping probes: %d", length(common_fs)))

ref_fs <- as.matrix(profiles_fs[match(common_fs, profiles_fs$probe),
                                 c("immune", "neuron", "glia", "V4", "V5", "V6")])
rownames(ref_fs) <- common_fs
beta_common_fs <- beta_subset[common_fs, ]

na_rows_fs <- rowSums(is.na(beta_common_fs)) > 0 | rowSums(is.na(ref_fs)) > 0
message(sprintf("  Removing %d rows with NAs", sum(na_rows_fs)))
ref_fs <- ref_fs[!na_rows_fs, ]
beta_common_fs <- beta_common_fs[!na_rows_fs, ]
message(sprintf("  Final probes for deconvolution: %d", nrow(ref_fs)))

message("  Running EpiDISH RPC (maxit=2000)...")
result_fs <- epidish(beta.m = beta_common_fs, ref.m = ref_fs,
                     method = "RPC", maxit = 2000)
fractions_fs <- as.data.frame(result_fs$estF)
colnames(fractions_fs) <- c("ITHfs_Immune", "ITHfs_Neuron", "ITHfs_Glia",
                              "ITHfs_Stem", "ITHfs_Diff2", "ITHfs_Diff1")
fractions_fs$Idat <- rownames(fractions_fs)
fractions_fs$ITHfs_Purity <- fractions_fs$ITHfs_Stem + fractions_fs$ITHfs_Diff1 +
                              fractions_fs$ITHfs_Diff2

rm(beta_subset, beta_common_fs, ref_fs, result_fs, profiles_full, profiles_fs)
gc(verbose = FALSE)

# ============================================================================
# Load DECODE Results & Merge
# ============================================================================

message("\nLoading DECODE results...")
decode <- read.csv(decode_results_file, stringsAsFactors = FALSE)
message(sprintf("  DECODE samples: %d", nrow(decode)))

# Merge all three
merged <- merge(decode, fractions_full, by = "Idat", all.x = TRUE)
merged <- merge(merged, fractions_fs, by = "Idat", all.x = TRUE)
message(sprintf("  Merged: %d samples", nrow(merged)))

# ============================================================================
# Comparison: DECODE vs ITHresolveGBM
# ============================================================================

message("\n==============================================================")
message("COMPARISON: DECODE vs ITHresolveGBM (Full Profiles)")
message("==============================================================\n")

# Tumor purity comparison
r_purity <- cor.test(merged$GBM_malignant, merged$ITH_Purity, use = "complete.obs")
r_purity_s <- cor.test(merged$GBM_malignant, merged$ITH_Purity, method = "spearman",
                       use = "complete.obs")
message(sprintf("Tumor purity: DECODE GBM_malignant vs ITH Purity"))
message(sprintf("  Pearson r = %.3f (p = %.2e)", r_purity$estimate, r_purity$p.value))
message(sprintf("  Spearman rho = %.3f (p = %.2e)", r_purity_s$estimate, r_purity_s$p.value))

# Component-level comparison
message("\n--- Component Correlations (Pearson r) ---\n")
comparisons <- list(
  c("GBM_malignant", "ITH_Purity", "Tumor (DECODE) vs Purity (ITH)"),
  c("Brain_total", "ITH_Glia", "Brain (DECODE) vs Glia (ITH)"),
  c("Neurons", "ITH_Neuron", "Neurons (DECODE) vs Neuron (ITH)"),
  c("Immune_total", "ITH_Immune", "Immune (DECODE) vs Immune (ITH)"),
  c("Microglia", "ITH_Immune", "Microglia (DECODE) vs Immune (ITH)"),
  c("Astrocytes", "ITH_Glia", "Astrocytes (DECODE) vs Glia (ITH)"),
  c("Oligodendrocytes", "ITH_Glia", "Oligo (DECODE) vs Glia (ITH)")
)

message(sprintf("%-45s %8s %10s %8s %10s", "Comparison", "r", "p", "rho", "p_s"))
message(paste(rep("-", 85), collapse = ""))
for (comp in comparisons) {
  if (comp[1] %in% names(merged) && comp[2] %in% names(merged)) {
    rt <- cor.test(merged[[comp[1]]], merged[[comp[2]]], use = "complete.obs")
    st <- cor.test(merged[[comp[1]]], merged[[comp[2]]], method = "spearman",
                   use = "complete.obs")
    sig <- ifelse(rt$p.value < 0.001, "***",
                  ifelse(rt$p.value < 0.01, "**",
                         ifelse(rt$p.value < 0.05, "*", "")))
    message(sprintf("%-45s %+8.3f %10.2e %+8.3f %10.2e %s",
                    comp[3], rt$estimate, rt$p.value, st$estimate, st$p.value, sig))
  }
}

# Score correlations comparison
message("\n--- Classifier Score Correlations ---\n")
valid <- merged %>% filter(!is.na(ClassifierScore))
message(sprintf("Samples with classifier score: %d\n", nrow(valid)))

score_vars <- c("GBM_malignant", "ITH_Purity", "ITH_Stem", "ITH_Differentiated",
                "Brain_total", "ITH_Glia", "ITH_Neuron",
                "Immune_total", "ITH_Immune",
                "Stromal_total")
message(sprintf("%-25s %8s %10s %8s %10s", "Variable", "r", "p", "rho", "p_s"))
message(paste(rep("-", 65), collapse = ""))
for (v in score_vars) {
  if (v %in% names(valid)) {
    rt <- cor.test(valid$ClassifierScore, valid[[v]], use = "complete.obs")
    st <- cor.test(valid$ClassifierScore, valid[[v]], method = "spearman",
                   use = "complete.obs")
    sig <- ifelse(rt$p.value < 0.001, "***",
                  ifelse(rt$p.value < 0.01, "**",
                         ifelse(rt$p.value < 0.05, "*", "")))
    message(sprintf("%-25s %+8.3f %10.2e %+8.3f %10.2e %s",
                    v, rt$estimate, rt$p.value, st$estimate, st$p.value, sig))
  }
}

# ============================================================================
# Summary Statistics
# ============================================================================

message("\n--- ITHresolveGBM Mean Fractions (All Samples) ---\n")
ith_vars <- c("ITH_Immune", "ITH_Neuron", "ITH_Glia", "ITH_Stem",
              "ITH_Diff1", "ITH_Diff2", "ITH_Purity")
message(sprintf("%-20s %8s %8s %8s %8s", "Component", "Mean%", "SD%", "Min%", "Max%"))
message(paste(rep("-", 56), collapse = ""))
for (v in ith_vars) {
  vals <- merged[[v]] * 100
  message(sprintf("%-20s %7.1f%% %7.1f%% %7.1f%% %7.1f%%",
                  v, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
                  min(vals, na.rm = TRUE), max(vals, na.rm = TRUE)))
}

# ============================================================================
# Visualizations
# ============================================================================

message("\nGenerating comparison plots...")

# Helper
sig_stars <- Vectorize(function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("n.s.")
})

# --------------------------------------------------------------------------
# Plot 1: DECODE vs ITH — Tumor purity scatter
# --------------------------------------------------------------------------

r_val <- cor(merged$GBM_malignant, merged$ITH_Purity, use = "complete.obs")
p1 <- ggplot(merged, aes(x = GBM_malignant * 100, y = ITH_Purity * 100)) +
  geom_point(size = 1.2, alpha = 0.5, color = "#E41A1C") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "grey50", linetype = "dotted") +
  theme_bw() +
  labs(
    title = "Tumor Purity: DECODE vs ITHresolveGBM",
    subtitle = sprintf("n=%d | Pearson r = %.3f | Dashed = regression, dotted = identity",
                       nrow(merged), r_val),
    x = "DECODE: GBM Malignant Fraction (%)",
    y = "ITHresolveGBM: Tumor Purity (Stem + Diff1 + Diff2) (%)"
  )
ggsave(file.path(fig_dir, "siewert_ithresolve_purity_comparison.pdf"), p1, width = 9, height = 8)
message("  Saved: siewert_ithresolve_purity_comparison.pdf")

# --------------------------------------------------------------------------
# Plot 2: Side-by-side composition comparison
# --------------------------------------------------------------------------

# Map both methods to comparable categories
decode_cat <- merged %>%
  transmute(
    Idat = Idat,
    Method = "DECODE",
    Tumor = GBM_malignant,
    Brain = Brain_total,
    Immune = Immune_total,
    Stromal = Stromal_total
  )

ith_cat <- merged %>%
  transmute(
    Idat = Idat,
    Method = "ITHresolveGBM",
    Tumor = ITH_Purity,
    Brain = ITH_Glia + ITH_Neuron,
    Immune = ITH_Immune,
    Stromal = 0  # ITH has no stromal component
  )

combined <- bind_rows(decode_cat, ith_cat) %>%
  pivot_longer(cols = c(Tumor, Brain, Immune, Stromal),
               names_to = "Category", values_to = "Fraction") %>%
  mutate(Category = factor(Category, levels = c("Tumor", "Brain", "Immune", "Stromal")),
         Fraction_Pct = Fraction * 100)

p2 <- ggplot(combined, aes(x = Method, y = Fraction_Pct, fill = Method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~Category, scales = "free_y") +
  scale_fill_manual(values = c("DECODE" = "#377EB8", "ITHresolveGBM" = "#E41A1C")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold")) +
  labs(
    title = "Cell Composition: DECODE vs ITHresolveGBM",
    subtitle = sprintf("Siewert Cohort (n=%d) | ITH has no Stromal component", nrow(merged)),
    x = NULL, y = "Fraction (%)", fill = NULL
  )
ggsave(file.path(fig_dir, "siewert_ithresolve_composition_comparison.pdf"), p2, width = 10, height = 8)
message("  Saved: siewert_ithresolve_composition_comparison.pdf")

# --------------------------------------------------------------------------
# Plot 3: Score correlations — both methods side by side
# --------------------------------------------------------------------------

valid_plot <- merged %>% filter(!is.na(ClassifierScore))

score_comp <- valid_plot %>%
  transmute(
    Idat = Idat,
    ClassifierScore = ClassifierScore,
    `DECODE: Tumor` = GBM_malignant,
    `ITH: Tumor Purity` = ITH_Purity,
    `DECODE: Brain` = Brain_total,
    `ITH: Glia+Neuron` = ITH_Glia + ITH_Neuron,
    `DECODE: Immune` = Immune_total,
    `ITH: Immune` = ITH_Immune
  ) %>%
  pivot_longer(cols = -c(Idat, ClassifierScore),
               names_to = "Variable", values_to = "Fraction") %>%
  mutate(
    Fraction_Pct = Fraction * 100,
    Method = ifelse(grepl("DECODE", Variable), "DECODE", "ITHresolveGBM"),
    Category = gsub("^(DECODE|ITH): ", "", Variable)
  )

# Compute per-facet stats
score_stats <- score_comp %>%
  group_by(Variable) %>%
  summarise(
    r = cor(ClassifierScore, Fraction, use = "complete.obs"),
    p = cor.test(ClassifierScore, Fraction)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r=%+.3f %s", r, sig_stars(p)))

y_pos <- score_comp %>%
  group_by(Variable) %>%
  summarise(y_pos = max(Fraction_Pct, na.rm = TRUE) * 0.95, .groups = "drop")
score_stats <- merge(score_stats, y_pos)
score_stats$x_pos <- min(valid_plot$ClassifierScore, na.rm = TRUE) +
  diff(range(valid_plot$ClassifierScore, na.rm = TRUE)) * 0.02

p3 <- ggplot(score_comp, aes(x = ClassifierScore, y = Fraction_Pct)) +
  geom_point(aes(color = Method), size = 0.8, alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed",
              linewidth = 0.5) +
  geom_text(data = score_stats,
            aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 1, size = 3, family = "mono") +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("DECODE" = "#377EB8", "ITHresolveGBM" = "#E41A1C")) +
  theme_bw() +
  theme(strip.text = element_text(size = 9, face = "bold"),
        legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 8)) +
  labs(
    title = "Classifier Score vs Cell Fractions: DECODE vs ITHresolveGBM",
    subtitle = sprintf("Siewert Cohort (n=%d)", nrow(valid_plot)),
    x = "Methylation Classifier Score (v12.8 Family)",
    y = "Fraction (%)",
    caption = "*** p<0.001, ** p<0.01, * p<0.05, n.s. not significant",
    color = NULL
  )
ggsave(file.path(fig_dir, "siewert_ithresolve_score_comparison.pdf"), p3, width = 12, height = 10)
message("  Saved: siewert_ithresolve_score_comparison.pdf")

# --------------------------------------------------------------------------
# Plot 4: ITH malignant subtypes vs classifier score
# --------------------------------------------------------------------------

ith_malig <- valid_plot %>%
  transmute(
    Idat = Idat,
    ClassifierScore = ClassifierScore,
    `Stem-like (V4)` = ITH_Stem,
    `Differentiated 1 (V6)` = ITH_Diff1,
    `Differentiated 2 (V5)` = ITH_Diff2,
    `Total Purity` = ITH_Purity
  ) %>%
  pivot_longer(cols = -c(Idat, ClassifierScore),
               names_to = "Component", values_to = "Fraction") %>%
  mutate(Fraction_Pct = Fraction * 100)

malig_stats <- ith_malig %>%
  group_by(Component) %>%
  summarise(
    r = cor(ClassifierScore, Fraction, use = "complete.obs"),
    p = cor.test(ClassifierScore, Fraction)$p.value,
    rho = cor(ClassifierScore, Fraction, method = "spearman", use = "complete.obs"),
    p_s = cor.test(ClassifierScore, Fraction, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r=%+.3f %s | rho=%+.3f %s", r, sig_stars(p), rho, sig_stars(p_s)))

y_pos_m <- ith_malig %>%
  group_by(Component) %>%
  summarise(y_pos = max(Fraction_Pct, na.rm = TRUE) * 0.95, .groups = "drop")
malig_stats <- merge(malig_stats, y_pos_m)
malig_stats$x_pos <- min(valid_plot$ClassifierScore, na.rm = TRUE) +
  diff(range(valid_plot$ClassifierScore, na.rm = TRUE)) * 0.02

p4 <- ggplot(ith_malig, aes(x = ClassifierScore, y = Fraction_Pct)) +
  geom_point(size = 1, alpha = 0.4, color = "#E41A1C") +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "black",
              linewidth = 0.6) +
  geom_text(data = malig_stats,
            aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 1, size = 3, family = "mono") +
  facet_wrap(~Component, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(strip.text = element_text(size = 10, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8)) +
  labs(
    title = "ITHresolveGBM: Malignant Subtypes vs Classifier Score",
    subtitle = sprintf("Siewert Cohort (n=%d) | V4 = stem-like, V5+V6 = differentiated",
                       nrow(valid_plot)),
    x = "Methylation Classifier Score (v12.8 Family)",
    y = "Fraction (%)",
    caption = "Stats: Pearson r and Spearman rho. *** p<0.001, ** p<0.01, * p<0.05, n.s. not significant"
  )
ggsave(file.path(fig_dir, "siewert_ithresolve_malignant_subtypes.pdf"), p4, width = 11, height = 9)
message("  Saved: siewert_ithresolve_malignant_subtypes.pdf")

# --------------------------------------------------------------------------
# Plot 5: Concordance heatmap — correlation matrix between all components
# --------------------------------------------------------------------------

cor_vars <- c("GBM_malignant", "Astrocytes", "Oligodendrocytes", "Microglia",
              "Neurons", "Immune_total", "Endothelial", "Fibroblasts",
              "ITH_Stem", "ITH_Diff1", "ITH_Diff2", "ITH_Immune",
              "ITH_Neuron", "ITH_Glia")
cor_labels <- c("DECODE: Tumor", "DECODE: Astrocytes", "DECODE: Oligo",
                "DECODE: Microglia", "DECODE: Neurons", "DECODE: Immune",
                "DECODE: Endothelial", "DECODE: Fibroblasts",
                "ITH: Stem-like", "ITH: Diff1", "ITH: Diff2",
                "ITH: Immune", "ITH: Neuron", "ITH: Glia")

cor_data <- merged[, cor_vars]
cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
colnames(cor_mat) <- cor_labels
rownames(cor_mat) <- cor_labels

# Melt for ggplot heatmap
cor_melt <- as.data.frame(as.table(cor_mat))
names(cor_melt) <- c("Var1", "Var2", "Correlation")

p5 <- ggplot(cor_melt, aes(x = Var1, y = Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 2.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 12)) +
  labs(
    title = "Cross-Method Correlation: DECODE vs ITHresolveGBM Components",
    subtitle = sprintf("Siewert Cohort (n=%d) | Pearson correlations", nrow(merged)),
    x = NULL, y = NULL, fill = "r"
  )
ggsave(file.path(fig_dir, "siewert_ithresolve_correlation_heatmap.pdf"), p5, width = 13, height = 11)
message("  Saved: siewert_ithresolve_correlation_heatmap.pdf")

# ============================================================================
# Save Results
# ============================================================================

# Save ITH results merged with DECODE
out_cols <- c("Idat", "ClassifierScore", "is_GBM", "Max_family",
              "GBM_malignant", "Brain_total", "Immune_total", "Stromal_total",
              "ITH_Immune", "ITH_Neuron", "ITH_Glia",
              "ITH_Stem", "ITH_Diff1", "ITH_Diff2", "ITH_Purity",
              "ITHfs_Stem", "ITHfs_Diff1", "ITHfs_Diff2", "ITHfs_Purity")
out_cols <- out_cols[out_cols %in% names(merged)]
write.csv(merged[, out_cols], output_csv, row.names = FALSE)
message(sprintf("\nSaved: %s", output_csv))

message("\n==============================================================")
message("Done! DECODE vs ITHresolveGBM comparison complete.")
message("==============================================================")
