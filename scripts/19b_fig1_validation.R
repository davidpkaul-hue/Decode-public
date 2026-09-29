#!/usr/bin/env Rscript
################################################################################
# DECODE: Unified Fig 1 (Validation + Cross-Method Benchmark)
################################################################################
# Purpose: One 4-panel figure for the manuscript that combines internal
#          validation (confusion matrix + mixture scatter) with cross-method
#          DECODE only. The GIMiCC benchmark was removed from the manuscript
#          on 2026-09-02 together with the cross-method section; the panels are
#          otherwise unchanged.
# Inputs:  data/processed/validation_results.rds                    (Panel a)
#          data/processed/validation_insilico_true_fractions.csv    (Panel b)
#          data/processed/validation_insilico_predicted_our.csv     (Panel b)
#          (Panel c is computed from the two mixture files above)
#          data/processed/siewert_cohort_results.csv                (Panel d)

# Output:  figures/Fig1_validation.pdf + .png
################################################################################

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(gridExtra)
})

theme_manuscript <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

col_decode <- "#2166AC"   # blue

# ============================================================================
# Panel a: Confusion matrix (DECODE held-out controls)
# ============================================================================
val <- readRDS("data/processed/validation_results.rds")
ctrl <- val$control_samples$proportions / 100

celltypes_order <- c("Bcell", "CD4T", "CD8T", "Mono", "NK", "Neu",
                     "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia",
                     "Endothelial", "Fibroblasts", "GBM_malignant")
# Written out rather than abbreviated: Capper asked for the cell type names to
# be readable without a key in the legend (review 28.09.2026).
ct_labels <- c("B cells", "CD4+ T cells", "CD8+ T cells", "Monocytes", "NK cells",
               "Neutrophils", "Neurons", "Oligodendrocytes", "Astrocytes",
               "Microglia", "Endothelial", "Fibroblasts", "GBM malignant")

ctrl_df <- data.frame(
  Sample = rep(ct_labels, each = 13),
  CellType = rep(ct_labels, 13),
  Fraction = as.vector(t(ctrl)))
ctrl_df$Sample <- factor(ctrl_df$Sample, levels = rev(ct_labels))
ctrl_df$CellType <- factor(ctrl_df$CellType, levels = ct_labels)

p_a <- ggplot(ctrl_df, aes(x = CellType, y = Sample, fill = Fraction)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(round(Fraction * 100) >= 1, sprintf("%.0f", Fraction * 100), ""),
                color = Fraction > 0.5),
            size = 2.4, show.legend = FALSE) +
  scale_color_manual(values = c(`FALSE` = "grey20", `TRUE` = "white")) +
  scale_fill_gradient(low = "white", high = col_decode, limits = c(0, 1),
                      labels = scales::percent_format(accuracy = 1)) +
  theme_manuscript +
  labs(tag = "a", x = "Predicted cell type", y = "True cell type", fill = "Fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        legend.key.height = unit(0.6, "cm"),
        legend.position = "right")

# ============================================================================
# Panel b: In-silico mixture scatter (DECODE)
# ============================================================================
mix_true   <- read.csv("data/processed/validation_insilico_true_fractions.csv", stringsAsFactors = FALSE)
mix_decode <- read.csv("data/processed/validation_insilico_predicted_our.csv", stringsAsFactors = FALSE)

mix_true_long   <- pivot_longer(mix_true,   -Sample, names_to = "CellType", values_to = "True")
mix_decode_long <- pivot_longer(mix_decode, -Sample, names_to = "CellType", values_to = "Predicted")

mix_true_long$True        <- mix_true_long$True / 100
mix_decode_long$Predicted <- mix_decode_long$Predicted / 100

df_decode <- merge(mix_true_long, mix_decode_long, by = c("Sample", "CellType"))
mix_df <- df_decode

mae_decode <- mean(abs(df_decode$True - df_decode$Predicted), na.rm = TRUE)
r_decode   <- cor(df_decode$True, df_decode$Predicted, use = "complete")

p_b <- ggplot(mix_df, aes(x = True, y = Predicted)) +
  geom_point(size = 0.9, alpha = 0.6, stroke = 0, colour = col_decode) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  theme_manuscript +
  labs(tag = "b", x = "True fraction", y = "Predicted fraction") +
  annotate("text", x = 0.52, y = 0.03,
           label = sprintf("r = %.3f, MAE = %.1f%%", r_decode, mae_decode * 100),
           size = 2.8, hjust = 0.5, color = col_decode) +
  # Square-root axes: 272 of the 299 points sit below a true fraction of 0.1,
  # where a linear axis piles them onto each other. Square root spreads that
  # range, keeps the identity line straight and, unlike a log or logit axis,
  # still shows the 43 true and 21 predicted zeros.
  scale_x_sqrt(limits = c(0, 1), breaks = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 1),
               expand = c(0.02, 0)) +
  scale_y_sqrt(limits = c(0, 1), breaks = c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 1),
               expand = c(0.02, 0)) +
  coord_equal() +
  theme(legend.position = "none")

# ============================================================================
# Panel c: Per-cell-type error, split into bias and scatter
# ----------------------------------------------------------------------------
# The mean ABSOLUTE error hides two things. First, the error is systematic in
# sign (91-100% of mixtures agree for the large biases), so a signed axis
# carries more information. Second, the biases pair up: oligodendrocytes +1.98 against
# astrocytes -1.91 (sum +0.07, i.e. the pooled glial compartment is unbiased and
# only its split is wrong), and CD4+ T +2.09 against CD8+ T -1.17 and NK -1.32,
# the same confusion panel a shows as 21% CD8 into CD4. Sorted by |bias| so the
# two pairs sit adjacent and their opposite directions are visible at a glance.
# Computed from the mixture files directly, not from the GIMiCC benchmark table.
# ============================================================================
err_df <- merge(mix_true_long, mix_decode_long, by = c("Sample", "CellType"))
err_df$Err <- 100 * (err_df$Predicted - err_df$True)          # percentage points

err_sum <- err_df %>%
  group_by(CellType) %>%
  # 95% CI of the mean (t, n = 23 mixtures), not the SD. The SD is dominated by the
  # mixture design: the error grows with the true fraction (|error| vs true r up to
  # 0.98), so a large SD here means "the mixtures differ", not "the estimate is
  # noisy", and it visually contradicted the text. The CI shows what the panel
  # claims: whether each bias is reliably different from zero.
  summarise(Bias = mean(Err),
            CI = qt(0.975, n() - 1) * sd(Err) / sqrt(n()), .groups = "drop") %>%
  as.data.frame()
ct_label_map <- setNames(ct_labels, celltypes_order)
err_sum$CellType <- ct_label_map[err_sum$CellType]
err_sum <- err_sum[order(abs(err_sum$Bias)), ]
err_sum$CellType <- factor(err_sum$CellType, levels = err_sum$CellType)

p_c <- ggplot(err_sum, aes(x = Bias, y = CellType)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = Bias - CI, xmax = Bias + CI),
                orientation = "y", width = 0.25, colour = "grey55", linewidth = 0.4) +
  geom_point(size = 2, colour = col_decode) +
  theme_manuscript +
  labs(tag = "c", x = "Mean error (percentage points, 95% CI)", y = NULL) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 8))

# ============================================================================
# Panel d: Tumor purity vs classifier score (Siewert cohort)
# ============================================================================
sie <- read.csv("data/processed/siewert_cohort_results.csv", stringsAsFactors = FALSE)

# Restrict to GBM IDH-wildtype, the cohort the caption reports. This previously
# ran on all 436 arrays including the 18 non-GBM.
sie_a <- sie[grepl("Glioblastoma, IDH", sie$Max_family), ]
message(sprintf("Panel d: %d GBM IDH-wildtype samples", nrow(sie_a)))

dfd <- data.frame(Purity = sie_a$GBM_malignant, Score = sie_a$ClassifierScore)
dfd <- dfd[complete.cases(dfd), ]
r_d_decode <- cor(sie_a$GBM_malignant, sie_a$ClassifierScore, use = "complete.obs")

p_d <- ggplot(dfd, aes(x = Purity, y = Score)) +
  geom_point(size = 0.8, alpha = 0.4, colour = col_decode) +
  theme_manuscript +
  coord_cartesian(ylim = c(0, 1)) +   # classifier scores are bounded [0,1]
  labs(tag = "d", x = "Estimated tumor purity", y = "Classifier score") +
  annotate("text", x = 0.02, y = 0.06,
           label = sprintf("r = %.3f", r_d_decode),
           size = 3, hjust = 0, color = col_decode) +
  theme(legend.position = "none")

# ============================================================================
# Compose 2x2
# ============================================================================
pdf("figures/Fig1_validation.pdf", width = 12, height = 10)
grid.arrange(p_a, p_b, p_c, p_d, ncol = 2, widths = c(1.1, 1), heights = c(1, 1))
dev.off()
png("figures/Fig1_validation.png", width = 12, height = 10, units = "in", res = 300)
grid.arrange(p_a, p_b, p_c, p_d, ncol = 2, widths = c(1.1, 1), heights = c(1, 1))
dev.off()

cat("Saved: figures/Fig1_validation.pdf\n")
cat("Saved: figures/Fig1_validation.png\n")

# Mirror PNG into manuscript figures folder (used by masterarbeit*.md, supplementary.md, progress_report.md)
manuscript_fig_dir <- "manuscript/figures"
if (dir.exists(manuscript_fig_dir)) {
  file.copy("figures/Fig1_validation.png",
            file.path(manuscript_fig_dir, "Fig1_validation.png"),
            overwrite = TRUE)
  cat("Mirrored: ", file.path(manuscript_fig_dir, "Fig1_validation.png"), "\n", sep = "")
}
