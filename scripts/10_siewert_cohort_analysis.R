#!/usr/bin/env Rscript
################################################################################
# DECODE: Siewert / Charité-Internal GBM Cohort Analysis
################################################################################
# Purpose: Deconvolution of the Charité-internal GBM cohort (Christin Siewert)
#          Beta-value matrix provided directly (no IDAT processing needed)
#          Correlation of cell fractions with methylation classifier scores
#
# Input:  /mnt/volume-hel1-1/siewert/betamatrix.rds              (beta-value matrix)
#         /mnt/volume-hel1-1/siewert/betamatrix.xlsx              (sample metadata)
#         /mnt/volume-hel1-1/siewert/classifier_v128_from_html.csv (v12.8 scores from HTML reports)
#         data/processed/reference_oneVsAll.rds
# Output: data/processed/siewert_cohort_results.csv
#         figures/siewert_*.pdf
################################################################################

if (file.exists("scripts/R45_patch.R")) {
  source("scripts/R45_patch.R")
}

suppressPackageStartupMessages({
  library(EpiDISH)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
})

message("==================================================")
message("DECODE: Siewert / Charité-Internal Cohort Analysis")
message("==================================================\n")

# ============================================================================
# Configuration
# ============================================================================

beta_file <- "/mnt/volume-hel1-1/siewert/betamatrix.rds"
meta_file <- "/mnt/volume-hel1-1/siewert/betamatrix.xlsx"
html_scores_file <- "/mnt/volume-hel1-1/siewert/classifier_v128_from_html.csv"
reference_file <- "data/processed/reference_oneVsAll.rds"
output_csv <- "data/processed/siewert_cohort_results.csv"
fig_dir <- "figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load Reference Matrix
# ============================================================================

message("Loading reference matrix...")
ref <- readRDS(reference_file)
sig <- ref$signature
marker_probes <- rownames(sig)
message(sprintf("  Signature: %d probes x %d cell types", nrow(sig), ncol(sig)))

# ============================================================================
# Load Metadata
# ============================================================================

message("Loading metadata...")
meta <- read_excel(meta_file)
message(sprintf("  Total entries: %d", nrow(meta)))

# Clean column names for easier access
meta <- meta %>%
  rename(
    Idat = `Idat no`,
    Diagnosis_text = `Aktuell beste Diagnose / Abschließende Diagnose`,
    Age = `Patient age at OP`,
    Sex = `Incoming gender`,
    Tumor_pct = `Tumor cell %`,
    EPIC_date = `EPIC scan date`,
    Genotype_check = `Genotype check EPIC`,
    p_detection = `Mean detection p-value`,
    Quality = `Technical Quality`,
    Tumor_location = `Standard tumor location`,
    Max_subclass = `Max sub-class v12.8 epignostix`,
    Max_subclass_score = `Max subClass Score v12.8 epignostix`,
    Max_class = `Max Class v12.8 Epignostix`,
    Max_class_score = `Max Class Score v12.8 epignostix`,
    Max_family = `Max Family v12.8 Epignostix`,
    Max_family_score = `Max Family Score v12.8 Epignostix`,
    MGMT_status = `mgmt status v12.8 Epignostix`,
    MGMT_methylation = `mgmt est. methylation v12.8 Epignostix`,
    IDH1_result = `IDH1 mutation result`,
    IDH1_exchange = `IDH1 aa exchange`,
    CNV_Load = `CNV_Load`,
    RF_purity = `RF purity prediction (absolute)`
  )

# ============================================================================
# Fill Missing v12.8 Scores from HTML Classifier Reports
# ============================================================================
# 63 samples originally had no v12.8 scores in the Excel metadata.
# Christin Siewert provided the HTML classifier reports for these IDATs,
# which were parsed to extract the v12.8 family/class/subclass scores.

if (file.exists(html_scores_file)) {
  html_scores <- read.csv(html_scores_file, stringsAsFactors = FALSE)
  message(sprintf("  Loaded v12.8 HTML scores for %d samples", nrow(html_scores)))

  n_filled <- 0
  for (i in seq_len(nrow(html_scores))) {
    idx <- which(meta$Idat == html_scores$Sample_ID[i])
    if (length(idx) == 1 && is.na(meta$Max_family_score[idx])) {
      meta$Max_family[idx] <- html_scores$family_predicted[i]
      meta$Max_family_score[idx] <- as.numeric(html_scores$family_score[i])
      meta$Max_class[idx] <- html_scores$class_predicted[i]
      meta$Max_class_score[idx] <- as.numeric(html_scores$class_score[i])
      meta$Max_subclass[idx] <- html_scores$subclass_predicted[i]
      meta$Max_subclass_score[idx] <- as.numeric(html_scores$subclass_score[i])
      n_filled <- n_filled + 1
    }
  }
  message(sprintf("  Filled %d/%d missing v12.8 scores from HTML reports", n_filled,
                  sum(is.na(meta$Max_family_score) | html_scores$Sample_ID %in% meta$Idat)))
} else {
  message("  WARN: HTML scores file not found, will use v11b4 fallback for missing v12.8")
}

message(sprintf("  Samples with v12.8 family score: %d / %d",
                sum(!is.na(meta$Max_family_score)), nrow(meta)))
message(sprintf("  Classifier families: %s",
                paste(names(table(meta$Max_family)), collapse = ", ")))

# ============================================================================
# Load Beta Matrix & Subset to Marker Probes
# ============================================================================

message("Loading beta-value matrix (this may take a moment)...")
beta_full <- readRDS(beta_file)
message(sprintf("  Full matrix: %d CpGs x %d samples", nrow(beta_full), ncol(beta_full)))

# Subset to marker probes immediately to save memory
common <- intersect(marker_probes, rownames(beta_full))
message(sprintf("  Marker probes found: %d / %d", length(common), length(marker_probes)))

beta <- beta_full[common, ]
rm(beta_full)
gc(verbose = FALSE)

# Align reference to available probes
sig_common <- sig[common, ]

# ============================================================================
# Match Samples: Beta Matrix <-> Metadata
# ============================================================================

message("\nMatching samples to metadata...")
beta_samples <- colnames(beta)
meta_idats <- meta$Idat

matched <- intersect(beta_samples, meta_idats)
message(sprintf("  Beta matrix samples: %d", length(beta_samples)))
message(sprintf("  Metadata entries:    %d", length(meta_idats)))
message(sprintf("  Matched:             %d", length(matched)))

unmatched_beta <- setdiff(beta_samples, meta_idats)
if (length(unmatched_beta) > 0) {
  message(sprintf("  WARN: %d samples in beta matrix without metadata", length(unmatched_beta)))
}

# ============================================================================
# Run Deconvolution (EpiDISH RPC)
# ============================================================================

message("\nRunning EpiDISH RPC deconvolution...")

# Use all samples in beta matrix (even unmatched - we can still deconvolve)
result <- epidish(beta.m = beta, ref.m = sig_common, method = "RPC")
fractions <- as.data.frame(result$estF)
fractions$Idat <- rownames(fractions)

message(sprintf("  Deconvolved: %d samples", nrow(fractions)))

# Free memory
rm(beta, result)
gc(verbose = FALSE)

# ============================================================================
# Merge with Metadata
# ============================================================================

message("\nMerging results with metadata...")
merged <- merge(meta, fractions, by = "Idat", all.y = TRUE)
message(sprintf("  Merged: %d samples", nrow(merged)))

# Add derived columns
cell_types <- c("GBM_malignant", "Astrocytes", "Oligodendrocytes", "Microglia",
                "Neurons", "Bcell", "CD4T", "CD8T", "Mono", "NK", "Neu",
                "Endothelial", "Fibroblasts")

merged$Immune_total <- merged$Bcell + merged$CD4T + merged$CD8T +
                       merged$Mono + merged$NK + merged$Neu
merged$Brain_total <- merged$Astrocytes + merged$Oligodendrocytes +
                      merged$Microglia + merged$Neurons
merged$Stromal_total <- merged$Endothelial + merged$Fibroblasts

# Use v12.8 classifier score (all 436 samples now have v12.8 from Excel + HTML reports)
merged$ClassifierScore <- merged$Max_family_score
n_v128 <- sum(!is.na(merged$ClassifierScore))
n_missing <- sum(is.na(merged$ClassifierScore))
message(sprintf("  Classifier scores: %d v12.8, %d missing", n_v128, n_missing))
if (n_missing > 0) {
  message("  WARN: Some samples still missing v12.8 score — check HTML parsing")
}

# ============================================================================
# Identify GBM IDH-wt Samples
# ============================================================================

message("\n=== Sample Classification ===\n")

# Use v12.8 classifier family for GBM identification (all samples have v12.8 now)
merged$is_GBM <- grepl("Glioblastoma, IDH-wildtype", merged$Max_family, ignore.case = TRUE)

message(sprintf("  GBM IDH-wt (by classifier):  %d", sum(merged$is_GBM, na.rm = TRUE)))
message(sprintf("  Other / non-GBM:              %d", sum(!merged$is_GBM, na.rm = TRUE)))

# Quality filter
if (!is.null(merged$p_detection)) {
  poor_quality <- !is.na(merged$p_detection) & merged$p_detection > 0.05
  message(sprintf("  Poor quality (p > 0.05):      %d", sum(poor_quality)))
}

# Genotype duplicates
if (!is.null(merged$Genotype_check)) {
  n_unique_genotypes <- length(unique(na.omit(merged$Genotype_check)))
  message(sprintf("  Unique genotype groups:       %d", n_unique_genotypes))
}

# ============================================================================
# Summary Statistics
# ============================================================================

message("\n=== Mean Cell Fractions (All Samples) ===\n")
message(sprintf("Total samples: %d (GBM IDH-wt: %d, Non-GBM: %d)\n",
                nrow(merged), sum(merged$is_GBM, na.rm = TRUE),
                sum(!merged$is_GBM, na.rm = TRUE)))

message(sprintf("%-20s %8s %8s %8s %8s", "Cell Type", "Mean%", "SD%", "Min%", "Max%"))
message(paste(rep("-", 56), collapse = ""))
for (ct in cell_types) {
  vals <- merged[[ct]] * 100
  message(sprintf("%-20s %7.1f%% %7.1f%% %7.1f%% %7.1f%%",
                  ct, mean(vals), sd(vals), min(vals), max(vals)))
}

# Category summary
message("\n--- Biological Categories ---")
for (cat in c("Brain_total", "Immune_total", "Stromal_total", "GBM_malignant")) {
  vals <- merged[[cat]] * 100
  message(sprintf("%-20s %7.1f%% (SD %.1f%%)", cat, mean(vals), sd(vals)))
}

# ============================================================================
# Correlation Analysis: ClassifierScore vs Cell Fractions
# ============================================================================

message("\n=== Correlation: Classifier Score vs Cell Fractions (All Samples) ===\n")

all_cell_types <- c("GBM_malignant", "Brain_total", "Immune_total", "Stromal_total",
                    cell_types)

message(sprintf("%-20s %8s %10s %8s %10s", "Cell Type", "r", "p", "rho", "p_spear"))
message(paste(rep("-", 60), collapse = ""))

valid <- merged %>% filter(!is.na(ClassifierScore))
message(sprintf("Samples with classifier score: %d\n", nrow(valid)))
for (ct in all_cell_types) {
  if (ct %in% names(valid) && nrow(valid) > 3) {
    r_test <- cor.test(valid$ClassifierScore, valid[[ct]], use = "complete.obs")
    s_test <- cor.test(valid$ClassifierScore, valid[[ct]],
                       method = "spearman", use = "complete.obs")
    sig <- ifelse(r_test$p.value < 0.001, "***",
                  ifelse(r_test$p.value < 0.01, "**",
                         ifelse(r_test$p.value < 0.05, "*", "")))
    message(sprintf("%-20s %+8.3f %10.2e %+8.3f %10.2e %s",
                    ct, r_test$estimate, r_test$p.value,
                    s_test$estimate, s_test$p.value, sig))
  }
}

# ============================================================================
# Visualizations
# ============================================================================

message("\nGenerating plots...")

plot_color <- "#E41A1C"

# Helper functions
fmt_p <- Vectorize(function(p) {
  if (p < 0.001) return("p < 0.001")
  sprintf("p = %.3f", p)
})
sig_stars <- Vectorize(function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  return("n.s.")
})

# Use all samples; filter only for valid classifier score in scatter plots
all_plot <- merged %>% filter(!is.na(ClassifierScore))
message(sprintf("  Samples for scatter plots (with score): %d", nrow(all_plot)))

# --------------------------------------------------------------------------
# Plot 1: Cell type composition boxplot (all samples)
# --------------------------------------------------------------------------

plot_data <- merged %>%
  select(Idat, all_of(cell_types)) %>%
  pivot_longer(cols = all_of(cell_types),
               names_to = "CellType", values_to = "Fraction") %>%
  mutate(Fraction_Pct = Fraction * 100)

p1 <- ggplot(plot_data, aes(x = reorder(CellType, -Fraction_Pct, FUN = median),
                             y = Fraction_Pct)) +
  geom_boxplot(fill = plot_color, alpha = 0.6, outlier.size = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(
    title = "Cell Type Composition: Siewert Charité Cohort",
    subtitle = sprintf("All samples (n=%d)", nrow(merged)),
    x = "Cell Type", y = "Fraction (%)"
  )
ggsave(file.path(fig_dir, "siewert_celltype_composition.pdf"), p1, width = 12, height = 7)
message("  Saved: siewert_celltype_composition.pdf")

# --------------------------------------------------------------------------
# Plot 2: Score vs GBM malignant (transformed axes + stats)
# --------------------------------------------------------------------------

r_p <- cor.test(all_plot$ClassifierScore, all_plot$GBM_malignant)
r_s <- cor.test(all_plot$ClassifierScore, all_plot$GBM_malignant, method = "spearman")
stats_text <- sprintf(
  "n=%d | Pearson r = %+.3f, %s %s | Spearman rho = %+.3f, %s %s",
  nrow(all_plot),
  r_p$estimate, fmt_p(r_p$p.value), sig_stars(r_p$p.value),
  r_s$estimate, fmt_p(r_s$p.value), sig_stars(r_s$p.value)
)

p2 <- ggplot(all_plot, aes(x = ClassifierScore, y = GBM_malignant * 100)) +
  geom_point(size = 1.5, alpha = 0.6, color = plot_color) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = plot_color) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0, size = 8, family = "mono")) +
  labs(
    title = "Classifier Score vs GBM Malignant Cell Fraction",
    subtitle = sprintf("Siewert Charité Cohort (all samples, n=%d)", nrow(all_plot)),
    x = "Methylation Classifier Score",
    y = "GBM Malignant Fraction (%)",
    caption = stats_text
  )
ggsave(file.path(fig_dir, "siewert_score_vs_GBM_malignant.pdf"), p2, width = 10, height = 8)
message("  Saved: siewert_score_vs_GBM_malignant.pdf")

# --------------------------------------------------------------------------
# Plot 3: Faceted cell types vs score (transformed + per-facet stats)
# --------------------------------------------------------------------------

cell_types_plot <- c("Astrocytes", "Microglia", "Oligodendrocytes", "Neurons",
                     "Bcell", "CD4T", "CD8T", "Mono", "NK", "Neu",
                     "Endothelial", "Fibroblasts")

plot_facet <- all_plot %>%
  select(Idat, ClassifierScore, all_of(cell_types_plot)) %>%
  pivot_longer(cols = all_of(cell_types_plot),
               names_to = "CellType", values_to = "Fraction") %>%
  mutate(Fraction_Pct = Fraction * 100)

facet_stats <- plot_facet %>%
  group_by(CellType) %>%
  summarise(
    r = cor(ClassifierScore, Fraction, use = "complete.obs"),
    p = cor.test(ClassifierScore, Fraction)$p.value,
    rho = cor(ClassifierScore, Fraction, method = "spearman", use = "complete.obs"),
    p_s = cor.test(ClassifierScore, Fraction, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r=%+.2f %s | rho=%+.2f %s",
                          r, sig_stars(p), rho, sig_stars(p_s)))

y_positions <- plot_facet %>%
  group_by(CellType) %>%
  summarise(y_pos = max(Fraction_Pct, na.rm = TRUE) * 0.95, .groups = "drop")
facet_labels <- merge(facet_stats, y_positions)
facet_labels$x_pos <- min(all_plot$ClassifierScore, na.rm = TRUE) +
  diff(range(all_plot$ClassifierScore, na.rm = TRUE)) * 0.02

p3 <- ggplot(plot_facet, aes(x = ClassifierScore, y = Fraction_Pct)) +
  geom_point(size = 0.8, alpha = 0.4, color = plot_color) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.5,
              color = plot_color) +
  geom_text(data = facet_labels,
            aes(x = x_pos, y = y_pos, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 2.8, family = "mono") +
  facet_wrap(~CellType, scales = "free_y", ncol = 4) +
  theme_bw() +
  theme(strip.text = element_text(size = 9),
        plot.caption = element_text(hjust = 0, size = 8)) +
  labs(
    title = "Classifier Score vs Individual Cell Type Fractions",
    subtitle = sprintf("Siewert Charité Cohort (all samples, n=%d)", nrow(all_plot)),
    x = "Methylation Classifier Score",
    y = "Cell Fraction (%)",
    caption = "Stats: Pearson r and Spearman rho. *** p<0.001, ** p<0.01, * p<0.05, n.s. not significant"
  )
ggsave(file.path(fig_dir, "siewert_score_vs_celltypes_faceted.pdf"), p3, width = 14, height = 10)
message("  Saved: siewert_score_vs_celltypes_faceted.pdf")

# --------------------------------------------------------------------------
# Plot 4: Biological categories vs score
# --------------------------------------------------------------------------

categories <- c("Brain", "Immune", "Stromal", "Tumor", "Non-Tumor")
cat_colors <- c("Brain" = "#377EB8", "Immune" = "#4DAF4A",
                "Stromal" = "#984EA3", "Tumor" = "#E41A1C",
                "Non-Tumor" = "#FF7F00")

all_plot$Brain <- all_plot$Astrocytes + all_plot$Neurons + all_plot$Oligodendrocytes
all_plot$Immune <- all_plot$Bcell + all_plot$CD4T + all_plot$CD8T +
                   all_plot$Mono + all_plot$NK + all_plot$Neu + all_plot$Microglia
all_plot$Stromal <- all_plot$Endothelial + all_plot$Fibroblasts
all_plot$Tumor <- all_plot$GBM_malignant
all_plot$`Non-Tumor` <- 1 - all_plot$GBM_malignant

plot_cat <- all_plot %>%
  select(Idat, ClassifierScore, all_of(categories)) %>%
  pivot_longer(cols = all_of(categories),
               names_to = "Category", values_to = "Fraction") %>%
  mutate(Fraction_Pct = Fraction * 100,
         Category = factor(Category, levels = categories))

cat_stats <- plot_cat %>%
  group_by(Category) %>%
  summarise(
    r = cor(ClassifierScore, Fraction, use = "complete.obs"),
    p = cor.test(ClassifierScore, Fraction)$p.value,
    rho = cor(ClassifierScore, Fraction, method = "spearman", use = "complete.obs"),
    p_s = cor.test(ClassifierScore, Fraction, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("r=%+.2f %s | rho=%+.2f %s",
                          r, sig_stars(p), rho, sig_stars(p_s)))

y_pos <- plot_cat %>%
  group_by(Category) %>%
  summarise(y_pos = max(Fraction_Pct, na.rm = TRUE) * 0.95, .groups = "drop")
cat_stats <- merge(cat_stats, y_pos)
cat_stats$x_pos <- min(all_plot$ClassifierScore, na.rm = TRUE) +
  diff(range(all_plot$ClassifierScore, na.rm = TRUE)) * 0.02

p4 <- ggplot(plot_cat, aes(x = ClassifierScore, y = Fraction_Pct)) +
  geom_point(aes(color = Category), size = 1.2, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed",
              linewidth = 0.7) +
  geom_text(data = cat_stats,
            aes(x = x_pos, y = y_pos, label = label),
            hjust = 0, vjust = 1, size = 3.2, family = "mono") +
  facet_wrap(~Category, scales = "free_y", ncol = 3) +
  scale_color_manual(values = cat_colors, guide = "none") +
  theme_bw() +
  theme(strip.text = element_text(size = 11, face = "bold"),
        plot.caption = element_text(hjust = 0, size = 8, family = "mono")) +
  labs(
    title = "Biological Categories vs Classifier Score",
    subtitle = sprintf("Siewert Charité Cohort (n=%d) | Non-Tumor = Brain + Immune + Stromal",
                       nrow(all_plot)),
    x = "Methylation Classifier Score",
    y = "Category Fraction (%)",
    caption = "Stats: Pearson r and Spearman rho. *** p<0.001, ** p<0.01, * p<0.05, n.s. not significant"
  )
ggsave(file.path(fig_dir, "siewert_gbm_categories_vs_score.pdf"), p4, width = 12, height = 9)
message("  Saved: siewert_gbm_categories_vs_score.pdf")

# --------------------------------------------------------------------------
# Plot 5: GBM subclass composition
# --------------------------------------------------------------------------

gbm_only <- merged %>% filter(is_GBM)
if ("Max_subclass" %in% names(gbm_only) && sum(!is.na(gbm_only$Max_subclass)) > 0) {
  subclass_plot <- gbm_only %>%
    filter(!is.na(Max_subclass)) %>%
    select(Idat, Max_subclass, all_of(cell_types)) %>%
    pivot_longer(cols = all_of(cell_types),
                 names_to = "CellType", values_to = "Fraction") %>%
    mutate(Fraction_Pct = Fraction * 100,
           Subclass = gsub("glioblastoma, IDH-wildtype, ", "", Max_subclass,
                           ignore.case = TRUE))

  p5 <- ggplot(subclass_plot, aes(x = CellType, y = Fraction_Pct, fill = Subclass)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Cell Type Composition by GBM Methylation Subclass",
      subtitle = sprintf("Siewert Charité GBM IDH-wt (n=%d)", nrow(gbm_only)),
      x = "Cell Type", y = "Fraction (%)", fill = "Subclass"
    )
  ggsave(file.path(fig_dir, "siewert_gbm_subclass_composition.pdf"), p5, width = 14, height = 8)
  message("  Saved: siewert_gbm_subclass_composition.pdf")
}

# ============================================================================
# Save Results
# ============================================================================

write.csv(merged, output_csv, row.names = FALSE)
message(sprintf("\nSaved: %s", output_csv))

message("\n==================================================")
message(sprintf("Done! Processed %d samples (%d GBM IDH-wt).",
                nrow(merged), sum(merged$is_GBM, na.rm = TRUE)))
message("==================================================")
