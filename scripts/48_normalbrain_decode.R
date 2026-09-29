#!/usr/bin/env Rscript
# ============================================================================
# Script 48: DECODE auf oeffentlichem Normalhirn (Laurens-Test, 2026-07-20)
# ----------------------------------------------------------------------------
# Frage: Ueberschaetzt DECODE den Normalhirn-Anteil (insbesondere Astrozyten)?
# In den MOSAIC-Tumorproben gab DECODE ~32% Normalhirn an, Single-cell nur ~5%.
# Laurens' tumorfreie Hirnproben zeigen, dass Single-cell Normalhirn sehr wohl
# findet (65%). Also muss DECODE selbst geprueft werden: auf echtem Normalhirn,
# wo die Wahrheit ungefaehr bekannt ist.
#
# Datensatz: GSE74193 (Jaffe et al.), bulk DLPFC, 450k, 278 erwachsene
# Kontrollen. Auswahlkriterien aus der GEO-Series-Matrix rekonstruiert (21.09.2026):
#   group == "Control"  (schliesst die 225 Schizophrenie-Proben aus)
#   dropsample == FALSE (QC-Flag der Autoren)
#   age >= 18           (schliesst die 172 fetalen/paediatrischen Proben aus)
# Das ergibt exakt 278 Proben, Alter 18,05 bis 85,20. KEINE Altersobergrenze --
# die frueher notierte Spanne "18-85" beschreibt den Bereich, nicht den Filter.
# Die Eingaben wurden am 21.09.2026 aus GSE74193 neu gebaut, weil sie in einem
# Session-Scratchpad lagen und verlorengegangen waren; der Lauf reproduziert die
# alten Fraktionen exakt (groesste Abweichung 0).
# Beta-Werte auf die 2225 Referenz-Marker
# reduziert (gestreamt, s. scratchpad/extract.sh).
#
# Zwei Pruefungen:
#   (1) Spezifitaet: GBM_malignant sollte auf Normalhirn ~0% sein.
#   (2) Plausibilitaet: Neurone/Astrozyten/Oligo vs Literatur und vs Laurens'
#       Single-cell-Normalhirn.
# Output: data/processed/normalbrain_decode_fractions.csv
#         figures/normalbrain_decode_vs_singlecell.png
# ============================================================================

suppressPackageStartupMessages({ library(EpiDISH); library(ggplot2) })

SP <- "/mnt/volume-hel1-1/public_normalbrain"   # dauerhaft, war bis 21.09.2026 ein Session-Scratchpad

ref <- readRDS("data/processed/reference_oneVsAll.rds")
if (is.list(ref) && !is.matrix(ref)) ref <- ref[[1]]
message(sprintf("Referenz: %d Probes x %d Zelltypen", nrow(ref), ncol(ref)))

meta <- read.csv(file.path(SP, "gse74193_adult_controls.csv"), stringsAsFactors = FALSE)
raw  <- read.csv(file.path(SP, "gse74193_markers.csv"), row.names = 1, check.names = FALSE)

# nur die Beta-Spalten der ausgewaehlten erwachsenen Kontrollen
keep <- paste0(meta$sample, "_Beta")
keep <- keep[keep %in% colnames(raw)]
B <- as.matrix(raw[, keep, drop = FALSE])
colnames(B) <- sub("_Beta$", "", colnames(B))
message(sprintf("Beta-Matrix: %d Probes x %d Proben", nrow(B), ncol(B)))

# Marker-Schnittmenge
common <- intersect(rownames(ref), rownames(B))
message(sprintf("gemeinsame Marker: %d von %d (%.1f%%)",
                length(common), nrow(ref), 100 * length(common) / nrow(ref)))
B <- B[common, , drop = FALSE]; R <- ref[common, , drop = FALSE]

# Proben mit zu vielen NAs raus, Rest zeilenweise imputieren
ok <- colMeans(is.na(B)) < 0.10
message(sprintf("Proben mit <10%% NA: %d von %d", sum(ok), length(ok)))
B <- B[, ok, drop = FALSE]
if (anyNA(B)) {
  rm_ <- rowMeans(B, na.rm = TRUE)
  idx <- which(is.na(B), arr.ind = TRUE)
  B[idx] <- rm_[idx[, 1]]
}

res <- epidish(beta.m = B, ref.m = R, method = "RPC")$estF
fr  <- as.data.frame(100 * res)
fr$sample <- rownames(fr)
fr <- merge(fr, meta[, c("sample", "age", "sex")], by = "sample")
write.csv(fr, "data/processed/normalbrain_decode_fractions.csv", row.names = FALSE)

ct <- colnames(res)
S <- data.frame(
  celltype = ct,
  median   = round(apply(fr[, ct], 2, median), 1),
  q25      = round(apply(fr[, ct], 2, quantile, 0.25), 1),
  q75      = round(apply(fr[, ct], 2, quantile, 0.75), 1))
S <- S[order(-S$median), ]
cat("\n=== DECODE auf Normalhirn (GSE74193 DLPFC, n =", ncol(B), ") ===\n")
print(S, row.names = FALSE)

brain <- c("Neurons", "Astrocytes", "Oligodendrocytes")
cat(sprintf("\nNormalhirn gesamt (Neurone+Astro+Oligo): median %.1f%%\n",
            median(rowSums(fr[, brain]))))
cat(sprintf("GBM_malignant (sollte ~0 sein):          median %.1f%%  (max %.1f%%)\n",
            median(fr$GBM_malignant), max(fr$GBM_malignant)))
immune <- c("Bcell", "CD4T", "CD8T", "Mono", "NK", "Neu", "Microglia")
cat(sprintf("Immun gesamt:                            median %.1f%%\n",
            median(rowSums(fr[, immune]))))

# --- Vergleich mit Laurens' Single-cell-Normalhirn ---------------------------
sc <- read.csv("/mnt/volume-hel1-1/laurens_zaschke/aging_normal/AgingGB_normal_annotation_level_1.csv",
               row.names = 1, check.names = FALSE)
scp <- 100 * sweep(sc, 1, rowSums(sc), "/")
cmp <- data.frame(
  celltype = c("Neurone", "Astrozyten", "Oligodendrozyten", "Myeloid/Mikroglia", "Tumor"),
  DECODE = c(median(fr$Neurons), median(fr$Astrocytes), median(fr$Oligodendrocytes),
             median(fr$Microglia + fr$Mono), median(fr$GBM_malignant)),
  Singlecell = c(median(scp$Neuron), median(scp$Astrocyte), median(scp$Oligodendrocyte),
                 median(scp$Myeloid), median(scp$Malignant)))
cmp[, 2:3] <- round(cmp[, 2:3], 1)
cat("\n=== DECODE (bulk DLPFC) vs Single-cell (tumorfreies Hirn, Laurens) ===\n")
print(cmp, row.names = FALSE)

L <- reshape(cmp, varying = c("DECODE", "Singlecell"), v.names = "pct",
             timevar = "Methode", times = c("DECODE (Methylierung)", "Single-cell (Laurens)"),
             direction = "long")
L$celltype <- factor(L$celltype, levels = cmp$celltype)
p <- ggplot(L, aes(celltype, pct, fill = Methode)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f", pct)), position = position_dodge(0.8),
            vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c("DECODE (Methylierung)" = "#08519c",
                               "Single-cell (Laurens)" = "#d95f02")) +
  labs(x = NULL, y = "Anteil (%)", fill = NULL,
       title = "Normalhirn: DECODE-Dekonvolution vs Single-cell",
       subtitle = sprintf("DECODE auf %d erwachsenen DLPFC-Kontrollen (GSE74193), Single-cell auf %d tumorfreien Hirnproben",
                          ncol(B), nrow(sc))) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave("figures/normalbrain_decode_vs_singlecell.png", p, width = 9, height = 5, dpi = 150)
message("wrote figures/normalbrain_decode_vs_singlecell.png")
