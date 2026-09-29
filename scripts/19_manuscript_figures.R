#!/usr/bin/env Rscript
################################################################################
# DECODE: Manuscript Figures (Nature Communications style)
################################################################################
# Creates 6 composite main figures for the manuscript.
# Fig 1: Validation (mixture recovery + score vs tumor + DECODE vs ITH)
# Fig 2: Cell type composition by subtype (stacked + boxplots)
# Fig 3: Reclassification after purification (cross-tab + flow + scores)
# Fig 4: Microenvironment (pre-computed from script 15e)
# Fig 5: t-SNE before/after (pre-computed from script 15b)
# Fig 6: MGMT purification (probability + changes + clinical)
################################################################################

if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

source("scripts/_purification_config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

set.seed(42)

message("================================================================")
message("DECODE: Manuscript Figures")
message("================================================================\n")

# Common theme
theme_manuscript <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        plot.tag = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9))

subtype_colors <- c("RTK I" = "#E41A1C", "RTK II" = "#377EB8",
                     "MES" = "#4DAF4A", "PNC" = "#FF7F00")

# ============================================================================
# Load all data
# ============================================================================

meta <- read.csv("data/processed/siewert_cohort_results.csv", stringsAsFactors = FALSE)
gbm <- meta[grepl("Glioblastoma, IDH", meta$Max_family), ]

ith <- read.csv("data/processed/siewert_ithresolve_full.csv", stringsAsFactors = FALSE)

gbm$Subtype <- NA_character_
gbm$Subtype[grepl("RTK1", gbm$Max_subclass)] <- "RTK I"
gbm$Subtype[grepl("RTK2", gbm$Max_subclass)] <- "RTK II"
gbm$Subtype[grepl("mesenchymal", gbm$Max_subclass, ignore.case = TRUE)] <- "MES"
gbm$Subtype[grepl("primitive neuronal", gbm$Max_subclass, ignore.case = TRUE)] <- "PNC"
gbm <- gbm[!is.na(gbm$Subtype), ]

gbm$Immune <- ith$Immune[match(gbm$Idat, ith$Idat)]
gbm$Neuron <- ith$Neuron[match(gbm$Idat, ith$Idat)]
gbm$Glia <- ith$Glia[match(gbm$Idat, ith$Idat)]
gbm$ITH_Purity <- ith$Purity[match(gbm$Idat, ith$Idat)]

# Purified classifier (DECODE-fraction re-run, Heidelberg mnp.v12epicv2)
# Canonical calls from data/processed/siewert_purified_decode3_calls.csv
pur_calls <- read.csv(PURIFIED_CALLS, stringsAsFactors = FALSE)
pur_labels <- data.frame(
  Idat = pur_calls$Sample,
  pur_class = pur_calls$pred_sub,
  stringsAsFactors = FALSE
)
pur_labels$Pur_Subtype <- NA_character_
pur_labels$Pur_Subtype[pur_labels$pur_class == "GBM_RTK1"] <- "RTK I"
pur_labels$Pur_Subtype[pur_labels$pur_class == "GBM_RTK2"] <- "RTK II"
pur_labels$Pur_Subtype[pur_labels$pur_class %in% c("GBM_MES_TYP", "GBM_MES_ATYP")] <- "MES"
pur_labels$Pur_Subtype[pur_labels$pur_class == "GBM_PNC"] <- "PNC"

gbm$Pur_Subtype <- pur_labels$Pur_Subtype[match(gbm$Idat, pur_labels$Idat)]

message(sprintf("GBM samples: %d (with purified class: %d)", nrow(gbm), sum(!is.na(gbm$Pur_Subtype))))

########################################################################
# FIGURE 1: DECODE Validation -> moved to scripts/19b_fig1_validation.R
# (4-panel unified figure: confusion matrix + mixtures + per-CT MAE + clinical)
########################################################################

# Run separately:  Rscript scripts/19b_fig1_validation.R

# Common cell-type labels still used by Figure 2 below:
celltypes_order <- c("Bcell", "CD4T", "CD8T", "Mono", "NK", "Neu",
                      "Neurons", "Oligodendrocytes", "Astrocytes", "Microglia",
                      "Endothelial", "Fibroblasts", "GBM_malignant")
ct_labels <- c("B cells", "CD4+ T", "CD8+ T", "Mono", "NK", "Neu",
               "Neurons", "Oligo", "Astro", "Microglia",
               "Endoth.", "Fibro", "GBM")

########################################################################
# FIGURE 2: Cell type composition by subtype
########################################################################

message("--- Figure 2: Cell type composition ---")

# Astrocytes and oligodendrocytes are pooled into a single glial compartment:
# the FANS-derived reference panels define astrocytes by exclusion and cannot
# separate the two reliably in bulk tissue (see manuscript Limitations).
gbm$Glia <- gbm$Astrocytes + gbm$Oligodendrocytes
# Rare immune types are pooled exactly as in Table 4: individually they sit near
# the detection limit of the deconvolution and their ordering is not interpreted.
gbm$OtherImmune <- gbm$CD4T + gbm$Bcell + gbm$NK + gbm$Neu

celltypes <- c("GBM_malignant", "Mono", "Microglia", "CD8T", "OtherImmune",
               "Neurons", "Glia", "Endothelial", "Fibroblasts")
celltype_labels <- c("Tumor", "Mono/Macro", "Microglia", "CD8+ T", "Other immune",
                      "Neurons", "Glia", "Endothelial", "Fibroblasts")
celltype_colors <- c(
  "Tumor" = "#E31A1C", "Mono/Macro" = "#FF7F00", "Microglia" = "#FDBF6F",
  "CD8+ T" = "#6A3D9A", "Other immune" = "#CAB2D6", "Neurons" = "#33A02C",
  # Endothelial was grey and sat directly on the brown glial band, which are hard
  # to separate at similar lightness; sky blue keeps the three top bands distinct.
  "Glia" = "#B15928", "Endothelial" = "#56B4E9",
  "Fibroblasts" = "#999999")

# Panel A: Stacked barplot
means <- gbm %>%
  group_by(Subtype) %>%
  summarise(across(all_of(celltypes), ~mean(.x, na.rm = TRUE)), .groups = "drop")
means_long <- means %>%
  pivot_longer(cols = all_of(celltypes), names_to = "CellType", values_to = "Fraction")
means_long$CellLabel <- celltype_labels[match(means_long$CellType, celltypes)]
means_long$CellLabel <- factor(means_long$CellLabel, levels = rev(celltype_labels))
means_long$Subtype <- factor(means_long$Subtype, levels = c("RTK I", "RTK II", "MES", "PNC"))

p2a <- ggplot(means_long, aes(x = Subtype, y = Fraction, fill = CellLabel)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  theme_manuscript +
  labs(tag = "a", x = "", y = "Mean fraction", fill = "Cell type") +
  theme(legend.position = "right", legend.key.size = unit(0.35, "cm"))

# Panel B: Key cell types boxplot
key_types <- c("GBM_malignant", "Mono", "Microglia", "CD8T")
key_labels <- c("Tumor", "Mono/Macro", "Microglia", "CD8+ T")
box_long <- gbm %>%
  select(Idat, Subtype, all_of(key_types)) %>%
  pivot_longer(cols = all_of(key_types), names_to = "CellType", values_to = "Fraction")
box_long$CellLabel <- key_labels[match(box_long$CellType, key_types)]
box_long$CellLabel <- factor(box_long$CellLabel, levels = key_labels)
box_long$Subtype <- factor(box_long$Subtype, levels = c("RTK I", "RTK II", "MES", "PNC"))

p2b <- ggplot(box_long, aes(x = Subtype, y = Fraction, fill = Subtype)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~CellLabel, scales = "free_y", nrow = 1) +
  # Square-root axis: a handful of extreme CD8+ T and microglia samples otherwise
  # stretch those panels to ~45% and flatten the boxes that carry the message.
  # sqrt keeps every point visible and maps zero to zero.
  scale_y_sqrt(labels = scales::percent_format(), breaks = c(0, 0.02, 0.05, 0.1, 0.2, 0.4)) +
  scale_fill_manual(values = subtype_colors) +
  theme_manuscript +
  labs(tag = "b", x = "", y = "Fraction (square-root scale)") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold", size = 10))

# Panel C: Score vs Tumor fraction scatter (logit y-axis to spread high scores)
# Subclass score, not class or family: this panel supports the claim about the
# reliability of the SUBTYPE call, and the subclass score is the only level that
# measures confidence in the RTK/MES assignment. (Table 8 deliberately uses the
# family score instead, as the purity proxy least affected by subtype ambiguity.)
# The classifier never reports more than 0.9999. Scores are clipped at 0.999 for
# display so the upper tail does not stretch the logit axis; the axis break is
# labelled accordingly so the pile-up is visible.
gbm$score_plot <- pmin(gbm$Max_subclass_score, 0.999)
p2c <- ggplot(gbm, aes(x = GBM_malignant, y = score_plot, color = Subtype)) +
  geom_point(size = 1.2, alpha = 0.5) +
  scale_color_manual(values = subtype_colors) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(trans = "logit",
                     breaks = c(0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999),
                     labels = c("0.1", "0.3", "0.5", "0.7", "0.9", "0.95", "0.99",
                                "\u2265 0.999")) +
  theme_manuscript +
  labs(tag = "c", x = "DECODE tumor fraction", y = "Classifier score (logit scale)", color = "Subtype") +
  theme(legend.position = c(0.85, 0.25), legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"))

fig2_layout <- function() {
  grid.arrange(
    arrangeGrob(p2a, p2b, ncol = 1, heights = c(1, 0.8)),
    arrangeGrob(p2c, nullGrob(), ncol = 1, heights = c(1, 0.05)),
    ncol = 2, widths = c(1.3, 1)
  )
}
pdf("figures/Fig2_composition.pdf", width = 14, height = 8); fig2_layout(); dev.off()
# PNG straight into the manuscript figure folder (this is the file the docx uses)
png("manuscript/figures/Fig2_composition.png", width = 14, height = 8,
    units = "in", res = 150); fig2_layout(); dev.off()
message("  Saved: figures/Fig2_composition.pdf + manuscript/figures/Fig2_composition.png")

########################################################################
# FIGURE 3: Reclassification after purification
########################################################################

message("--- Figure 3: Reclassification ---")

gbm_pur <- gbm[!is.na(gbm$Pur_Subtype), ]

# Panel A: Cross-tabulation heatmap.
# Low-confidence non-GBM calls get their own column so that rows sum to the
# subtype totals reported in Table 5 (otherwise 7 MES, 8 RTK II and 3 RTK I
# samples silently disappear).
ct_src <- gbm[!is.na(gbm$Subtype), ]
ct_src$PurCol <- ifelse(is.na(ct_src$Pur_Subtype), "non-GBM", ct_src$Pur_Subtype)
ct <- as.data.frame(table(Original = ct_src$Subtype, Purified = ct_src$PurCol))
ct$Original <- factor(ct$Original, levels = c("RTK I", "RTK II", "MES", "PNC"))
ct$Purified <- factor(ct$Purified, levels = c("RTK I", "RTK II", "MES", "PNC", "non-GBM"))

# The cells where the purified call equals the original one carry the message of
# the panel, so they are outlined rather than left to be found by eye. The
# non-GBM column is separated by a rule because it is not a subtype: those calls
# are the low-confidence artifacts discussed in the text.
ct$stable <- as.character(ct$Original) == as.character(ct$Purified)
# Fill is the share of the ORIGINAL class, not the raw count: on a count scale the
# two large cells (169, 149) dominate and the 82% RTK I stability (62 of 76) reads
# as pale noise. The absolute numbers stay printed in the cells, so nothing is lost.
ct$RowFrac <- ct$Freq / ave(ct$Freq, ct$Original, FUN = sum)
sep_x <- length(levels(ct$Purified)) - 0.5   # zwischen PNC und non-GBM

p3a <- ggplot(ct, aes(x = Purified, y = Original, fill = RowFrac)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_tile(data = subset(ct, stable), fill = NA, color = "grey15", linewidth = 1.1) +
  geom_text(aes(label = Freq), size = 5, fontface = "bold") +
  annotate("segment", x = sep_x, xend = sep_x, y = 0.4, yend = length(levels(ct$Original)) + 0.6,
           linewidth = 0.8, colour = "grey35", linetype = "22") +
  scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1),
                      labels = function(x) paste0(round(100 * x), "%")) +
  theme_manuscript +
  labs(tag = "a", x = "After purification", y = "Before purification",
       fill = "Share of\noriginal class") +
  theme(panel.grid = element_blank()) +
  coord_equal(clip = "off")

# Panel C: Score improvement boxplots (class-level score, original vs purified)
# Subclass level throughout: this compares confidence in the SUBTYPE call before
# and after purification, so the subclass score is the matching quantity. The
# class score would mix in the (near-ceiling) confidence that the sample is a GBM
# at all. Only Table 8 uses the family score, and for a different purpose.
gbm_pur$orig_score <- gbm_pur$Max_subclass_score
gbm_pur$pur_score <- pur_calls$score_sub[match(gbm_pur$Idat, pur_calls$Sample)]

# Paired test behind the sentence in the Results: the class-level score of the
# originally MES-classified samples before vs after purification, same samples,
# so signed-rank rather than rank-sum. Previously computed ad hoc and reported
# without a script.
paired <- do.call(rbind, lapply(c("MES", "RTK I", "RTK II", "PNC"), function(st) {
  d <- gbm_pur[gbm_pur$Subtype == st & !is.na(gbm_pur$orig_score) & !is.na(gbm_pur$pur_score), ]
  if (nrow(d) < 3) return(NULL)
  data.frame(Subtype = st, n = nrow(d),
             median_original = median(d$orig_score),
             median_purified = median(d$pur_score),
             p_paired_wilcoxon = suppressWarnings(
               wilcox.test(d$orig_score, d$pur_score, paired = TRUE)$p.value))
}))
print(paired, row.names = FALSE)
write.csv(paired, "data/processed/classifier_score_paired_tests.csv", row.names = FALSE)

score_df <- data.frame(
  Stage = rep(c("Original", "Purified"), each = nrow(gbm_pur)),
  Score = c(gbm_pur$orig_score, gbm_pur$pur_score),
  Subtype = rep(gbm_pur$Subtype, 2)
)
score_df$Stage <- factor(score_df$Stage, levels = c("Original", "Purified"))
# Same subtype order as panel a and every table; without this the facets sort
# alphabetically and the figure contradicts itself.
score_df$Subtype <- factor(score_df$Subtype, levels = c("RTK I", "RTK II", "MES", "PNC"))

p3c <- ggplot(score_df, aes(x = Stage, y = Score, fill = Stage)) +
  # PNC has n = 3: a box with whiskers and outliers implies structure that three
  # values cannot carry, so those are drawn as the individual points they are.
  geom_boxplot(data = subset(score_df, Subtype != "PNC"),
               alpha = 0.7, outlier.size = 0.5) +
  geom_point(data = subset(score_df, Subtype == "PNC"), aes(colour = Stage),
             size = 1.8, show.legend = FALSE,
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  facet_wrap(~Subtype, nrow = 1) +
  scale_fill_manual(values = c("Original" = "#D6604D", "Purified" = "#4393C3")) +
  scale_colour_manual(values = c("Original" = "#D6604D", "Purified" = "#4393C3")) +
  theme_manuscript +
  labs(tag = "b", x = "", y = "Subclass-level classifier score") +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

# Two panels since 2026-09-02: the MES-flow bars were the MES row of panel a
# replotted, and were cited nowhere in the text.
fig3_layout <- function() {
  grid.arrange(p3a, p3c, ncol = 1, heights = c(1, 0.8))
}
pdf("figures/Fig3_reclassification.pdf", width = 10, height = 9); fig3_layout(); dev.off()
# PNG straight into the manuscript figure folder (this is the file the docx uses),
# so the figure no longer has to be converted from the PDF by hand.
png("manuscript/figures/Fig3_reclassification.png", width = 10, height = 9,
    units = "in", res = 150); fig3_layout(); dev.off()
message("  Saved: figures/Fig3_reclassification.pdf + manuscript/figures/Fig3_reclassification.png")

gc(verbose = FALSE)

########################################################################
# FIGURE 4: Microenvironment (pre-computed)
########################################################################

message("--- Figure 4: Microenvironment ---")
message("  -> Pre-computed from script 15e:")
message("  -> tsne_microenv_combat_subtype.pdf (panel a)")
message("  -> tsne_microenv_combat_immune_types.pdf (panel b)")
message("  -> tsne_microenv_combat_purity.pdf (panel c)")
message("  Assemble manually or use pre-computed PDFs directly.")

########################################################################
# FIGURE 5: t-SNE before/after (pre-computed)
########################################################################

message("--- Figure 5: t-SNE before/after purification ---")
message("  -> Pre-computed from script 15b:")
message("  -> tsne_combat_comparison_subtype.pdf (panel a)")
message("  -> tsne_combat_mes_migration.pdf (panel b)")
message("  Assemble manually or use pre-computed PDFs directly.")

########################################################################
# MGMT figure removed 2026-08-30
########################################################################
# This script used to build a fourth panel set (Fig6_mgmt.pdf) from
# siewert_mgmt_purification.csv (script 16) and mgmt_pyro_metrics.csv
# (script 30). MGMT was dropped from the manuscript in July 2026 after the
# pyrosequencing validation showed no clinical benefit, so the figure has no
# place in the paper and the dependency on scripts 16 and 30 was the only
# thing keeping them in the thesis pipeline. The MGMT analysis itself lives
# on in scripts 16 and 30.

message("\n================================================================")
message("Manuscript figures generated: Fig 2 and Fig 3 (Fig 1 comes from script 19b)")
message("Fig4 (microenvironment) and Fig5 (t-SNE) use pre-computed PDFs")
message("================================================================")
