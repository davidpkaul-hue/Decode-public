#!/usr/bin/env Rscript
# ============================================================================
# Script 49: Gegenpruefung von Script 48 an einem zweiten Normalhirn-Datensatz
# ----------------------------------------------------------------------------
# Script 48 (GSE74193, 450k, DLPFC, n=278) ergab: DECODE ueberschaetzt
# Astrozyten massiv (38,8% vs ~13% Single-cell / ~15-25% Literatur), bei sonst
# plausiblem Gesamtbild und guter Tumor-Spezifitaet.
#
# Diese Replikation nutzt einen unabhaengigen Datensatz auf ANDERER PLATTFORM:
#   GSE284764, Infinium EPIC, praefrontaler Cortex, Newcastle Brain Bank,
#   111 Kontrollen + 131 Alzheimer.
# Damit laesst sich pruefen, ob die Astrozyten-Inflation (a) reproduzierbar und
# (b) kein 450k-Plattformartefakt ist. AD-Faelle werden getrennt ausgewertet,
# da Astrogliose dort einen ECHTEN Astrozyten-Anstieg erzeugen kann.
# Output: data/processed/normalbrain_replication_fractions.csv
#         figures/normalbrain_replication.png
# ============================================================================

suppressPackageStartupMessages({ library(EpiDISH); library(ggplot2) })

RAW <- "/mnt/volume-hel1-1/public_normalbrain/GSE284764_betas.RData"
SUB <- "/mnt/volume-hel1-1/public_normalbrain/GSE284764_markers.rds"

ref <- readRDS("data/processed/reference_oneVsAll.rds")
if (is.list(ref) && !is.matrix(ref)) ref <- ref[[1]]

# --- Beta-Matrix laden und sofort auf die Marker reduzieren ------------------
if (!file.exists(SUB)) {
  message("lade grosse Beta-Matrix (einmalig) ...")
  e <- new.env()
  load(RAW, envir = e)   # load() erkennt die gzip-Kompression der RData selbst
  obj <- ls(e)
  message("Objekte in RData: ", paste(obj, collapse = ", "))
  B <- get(obj[which.max(sapply(obj, function(o) length(get(o, envir = e))))], envir = e)
  if (is.data.frame(B)) B <- as.matrix(B)
  message(sprintf("volle Matrix: %d x %d", nrow(B), ncol(B)))
  # EPICv2-Suffixe abstreifen, falls vorhanden
  rownames(B) <- sub("_.*$", "", rownames(B))
  keep <- intersect(rownames(ref), rownames(B))
  B <- B[keep, , drop = FALSE]
  saveRDS(B, SUB); rm(e); gc()
  message("Marker-Subset gespeichert: ", SUB)
} else {
  B <- readRDS(SUB)
}
message(sprintf("Marker-Matrix: %d Probes x %d Proben", nrow(B), ncol(B)))

# --- Metadaten aus der Series-Matrix ----------------------------------------
hdr <- "/tmp/h_764.txt"
h <- readLines(hdr)
unq <- function(v) gsub('^"|"$', "", v)
rowv <- function(ln) unq(strsplit(ln, "\t")[[1]][-1])
gsm <- rowv(grep("^!Sample_geo_accession", h, value = TRUE)[1])
ch  <- grep("^!Sample_characteristics_ch1", h, value = TRUE)
getch <- function(key) sub("^[^:]*: ", "", rowv(ch[grep(key, ch, fixed = TRUE)][1]))
meta <- data.frame(gsm = gsm, dx = getch("ad diagnosis"),
                   sex = getch("Sex"), stringsAsFactors = FALSE)
message("Diagnosen: ", paste(names(table(meta$dx)), table(meta$dx), collapse = ", "))

# Spalten der Beta-Matrix auf GSM abbilden
cn <- colnames(B)
if (!any(cn %in% meta$gsm)) {
  message("Hinweis: Spaltennamen sind keine GSM-IDs (", paste(head(cn, 2), collapse = ", "),
          ") -> Reihenfolge der Series-Matrix wird angenommen")
  if (ncol(B) == nrow(meta)) colnames(B) <- meta$gsm else
    stop("Spaltenzahl passt nicht zur Metadatenzahl, Zuordnung unklar")
}
meta <- meta[match(colnames(B), meta$gsm), ]

# --- Deconvolution -----------------------------------------------------------
common <- intersect(rownames(ref), rownames(B))
message(sprintf("gemeinsame Marker: %d von %d (%.1f%%)",
                length(common), nrow(ref), 100 * length(common) / nrow(ref)))
Bm <- B[common, , drop = FALSE]; R <- ref[common, , drop = FALSE]
ok <- colMeans(is.na(Bm)) < 0.10
Bm <- Bm[, ok, drop = FALSE]; meta <- meta[ok, ]
if (anyNA(Bm)) { rm_ <- rowMeans(Bm, na.rm = TRUE)
  idx <- which(is.na(Bm), arr.ind = TRUE); Bm[idx] <- rm_[idx[, 1]] }

fr <- as.data.frame(100 * epidish(beta.m = Bm, ref.m = R, method = "RPC")$estF)
fr$gsm <- rownames(fr); fr$dx <- meta$dx
write.csv(fr, "data/processed/normalbrain_replication_fractions.csv", row.names = FALSE)

ct <- colnames(R)
show <- function(lab, d) {
  cat("\n---", lab, "(n =", nrow(d), ") ---\n")
  s <- data.frame(celltype = ct, median = round(apply(d[, ct], 2, median), 1))
  print(s[order(-s$median), ], row.names = FALSE)
  cat(sprintf("Normalhirn (Neurone+Astro+Oligo): %.1f%% | GBM_malignant: %.1f%%\n",
              median(rowSums(d[, c("Neurons", "Astrocytes", "Oligodendrocytes")])),
              median(d$GBM_malignant)))
}
ctl <- fr[fr$dx == "Ctl", ]; ad <- fr[fr$dx != "Ctl", ]
show("GSE284764 Kontrollen (EPIC, PFC)", ctl)
show("GSE284764 Alzheimer (EPIC, PFC)", ad)

# --- Vergleich der beiden Datensaetze + Single-cell ---------------------------
d48 <- read.csv("data/processed/normalbrain_decode_fractions.csv")
sc <- read.csv("/mnt/volume-hel1-1/laurens_zaschke/aging_normal/AgingGB_normal_annotation_level_1.csv",
               row.names = 1, check.names = FALSE)
scp <- 100 * sweep(sc, 1, rowSums(sc), "/")
cmp <- data.frame(
  celltype = rep(c("Neurone", "Astrozyten", "Oligodendrozyten", "Tumor"), 3),
  Quelle = rep(c("DECODE 450k (GSE74193, n=278)",
                 sprintf("DECODE EPIC (GSE284764 Ctl, n=%d)", nrow(ctl)),
                 "Single-cell (Laurens, n=10)"), each = 4),
  pct = c(median(d48$Neurons), median(d48$Astrocytes), median(d48$Oligodendrocytes), median(d48$GBM_malignant),
          median(ctl$Neurons), median(ctl$Astrocytes), median(ctl$Oligodendrocytes), median(ctl$GBM_malignant),
          median(scp$Neuron), median(scp$Astrocyte), median(scp$Oligodendrocyte), median(scp$Malignant)))
cmp$pct <- round(cmp$pct, 1)
cmp$celltype <- factor(cmp$celltype, levels = c("Neurone", "Astrozyten", "Oligodendrozyten", "Tumor"))
cat("\n=== Zusammenfassung ===\n"); print(cmp, row.names = FALSE)

p <- ggplot(cmp, aes(celltype, pct, fill = Quelle)) +
  geom_col(position = position_dodge(0.8), width = 0.75, colour = "grey30", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f", pct)), position = position_dodge(0.8),
            vjust = -0.35, size = 3.1) +
  annotate("rect", xmin = 1.55, xmax = 2.45, ymin = -Inf, ymax = Inf, alpha = 0.07, fill = "red") +
  scale_fill_manual(values = c("#08519c", "#3182bd", "#d95f02")) +
  labs(x = NULL, y = "Anteil (%)", fill = NULL,
       title = "Normalhirn: DECODE auf zwei unabhaengigen Datensaetzen vs Single-cell",
       subtitle = "Astrozyten (rot hinterlegt) werden auf beiden Plattformen ueberschaetzt") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(nrow = 2))
ggsave("figures/normalbrain_replication.png", p, width = 9.5, height = 5.4, dpi = 150)
message("wrote figures/normalbrain_replication.png")
