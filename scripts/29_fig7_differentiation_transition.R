#!/usr/bin/env Rscript
################################################################################
# Manuscript Figure 5 (file Fig6_differentiation_score.png): ITHresolveGBM
# differentiation score across the
# original -> purified subtype transition (DECODE-fraction purification).
################################################################################
# Reframed 2026-06-12 (David's decision): the earlier two-trajectory contrast
# (MES->RTK I vs MES->RTK II) no longer holds under the corrected DECODE-fraction
# purification, where the MES->RTK I subgroup collapses to n=2. Instead we show
# the differentiation-score GRADIENT across the four transition groups:
#   RTK I->RTK I (stem-like)  >  MES->RTK I (n=5, sits at zero)
#   >  MES->RTK II  >  RTK II->RTK II (differentiated).
# Message: after purification, mesenchymal tumours collapse almost entirely onto
# the differentiated RTK-II end of the axis (MES->RTK II n=165 vs MES->RTK I n=2).
#
# DiffScore = Malignant 1 - (Malignant 2 + 3) (ITHresolveGBM), from
# siewert_ithresolve_full.csv. Purified class from the Heidelberg classifier
# (DECODE-fraction re-run), via siewert_purified_decode3_calls.csv.
################################################################################
suppressPackageStartupMessages(library(ggplot2))

simplify <- function(x) {
  out <- rep("other", length(x))
  out[grepl("RTK1", x)] <- "RTK I"
  out[grepl("RTK2", x)] <- "RTK II"
  out[grepl("MES|mesenchymal", x, ignore.case = TRUE)] <- "MES"
  out[grepl("PNC|primitive neuronal", x, ignore.case = TRUE)] <- "PNC"
  out
}

# --- Differentiation scores + original class ---
ith  <- read.csv("data/processed/siewert_ithresolve_full.csv", stringsAsFactors = FALSE)
meta <- read.csv("data/processed/siewert_cohort_results.csv", stringsAsFactors = FALSE, check.names = FALSE)
# Restrict to GBM IDH-wildtype first. Without this filter, simplify() also matches
# non-GBM entities whose subclass string happens to contain RTK/mesenchymal --
# 4 control brain tissues and 2 paediatric high-grade gliomas were being pulled
# into the transition groups, inflating MES->RTK II from 165 to 169 and
# RTK I->RTK I from 60 to 62 relative to Table 5.
meta <- meta[grepl("Glioblastoma, IDH", meta$Max_family), ]
meta$orig <- simplify(meta$Max_subclass)

# --- Purified class from the DECODE-fraction Heidelberg re-run ---
source("scripts/_purification_config.R")
calls <- read.csv(PURIFIED_CALLS, stringsAsFactors = FALSE)
calls$pur <- simplify(calls$pred_sub)

# --- Merge ---
# Renormalised over the three malignant components, as the contamination
# experiment already does (Supplementary Methods). The raw difference
# Stem - (Diff1 + Diff2) scales with the malignant fraction, and that fraction
# differs systematically between the transition groups (median ITHresolveGBM
# purity 33% in MES against 62% in RTK), so the raw score confounds
# differentiation state with tumour content. Dividing by the malignant total
# puts every sample on the same scale.
d <- data.frame(Idat = ith$Idat,
                DiffScore = (ith$Stem - (ith$Diff1 + ith$Diff2)) / ith$Purity,
                stringsAsFactors = FALSE)
d$orig <- meta$orig[match(d$Idat, meta$Idat)]
d$pur  <- calls$pur[match(d$Idat, calls$Sample)]
d <- d[!is.na(d$DiffScore) & !is.na(d$orig) & !is.na(d$pur), ]

d$Transition <- paste(d$orig, "->", d$pur)
# Ordered stem-like to differentiated by the renormalised medians. On the raw
# score MES -> RTK II sat between the two RTK groups; renormalised it is the
# most differentiated group of the four, so the order changes here.
keep <- c("RTK I -> RTK I", "MES -> RTK I", "RTK II -> RTK II", "MES -> RTK II")
d <- d[d$Transition %in% keep, ]
d$Transition <- factor(d$Transition, levels = keep)

ns  <- table(d$Transition)
med <- tapply(d$DiffScore, d$Transition, median)
message("n per group: ", paste(sprintf("%s=%d", levels(d$Transition), ns), collapse = " | "))
message("median DiffScore: ", paste(sprintf("%s=%.2f", names(med), med), collapse = " | "))

# --- Gradient statistic: Kruskal-Wallis across the four ordered groups ---
kw <- kruskal.test(DiffScore ~ Transition, data = d)
message(sprintf("Kruskal-Wallis across groups: p = %.2e", kw$p.value))
# Purified MES (MES->RTK II) vs the differentiated anchor (RTK II->RTK II)
w_diff <- wilcox.test(d$DiffScore[d$Transition == "MES -> RTK II"],
                      d$DiffScore[d$Transition == "RTK II -> RTK II"])
message(sprintf("MES->RTK II vs RTK II->RTK II: p = %.2e", w_diff$p.value))

# --- Plot: gradient boxplot (stem-like -> differentiated) ---
cols <- c("RTK I -> RTK I" = "#E41A1C", "MES -> RTK I" = "#F4A6A6",
          "MES -> RTK II" = "#A6CEE3", "RTK II -> RTK II" = "#377EB8")
xlabs <- sprintf("%s\n(n=%d)", levels(d$Transition), as.integer(ns))

p7 <- ggplot(d, aes(Transition, DiffScore, fill = Transition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.size = 0.6, alpha = 0.9, width = 0.6) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = xlabs) +
  # No "negligible" label any more. Renormalised, the MES -> RTK I group sits at
  # exactly zero and its size is already on the x axis as (n=5), so the word
  # added nothing the plot did not show, and the caption had to explain it.
  # Das Label sass auf dem Whisker der ersten Gruppe. Jetzt oberhalb aller Daten,
  # mit etwas Luft nach oben, damit nichts ueberlappt.
  annotate("text", x = 0.5, y = max(d$DiffScore, na.rm = TRUE) +
             0.10 * diff(range(d$DiffScore, na.rm = TRUE)),
           label = sprintf("Kruskal-Wallis p = %.1e", kw$p.value),
           size = 4, hjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major.x = element_blank()) +
  # Kein Untertitel: die Orientierung der Achse (stem-like links, differenziert
  # rechts) steht in der Abbildungslegende im Manuskript.
  labs(x = NULL,
       y = "Differentiation score, renormalized\n(Malignant 1 - (Malignant 2+3)) / malignant total")

ggsave("figures/Fig6_differentiation_score.png", p7, width = 9, height = 7, dpi = 300)
ggsave("manuscript/figures/Fig6_differentiation_score.png", p7, width = 9, height = 7, dpi = 300)
message("Saved: Fig6_differentiation_score.png (figures/ + manuscript/figures/)")
