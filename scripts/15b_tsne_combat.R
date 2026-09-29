#!/usr/bin/env Rscript
################################################################################
# DECODE: t-SNE with ComBat batch correction (EPIC v1/v2) for comparison
# Produces manuscript FIGURE 8d,e (files Fig5a/b_tsne_*.png -- the file names predate
# the current figure numbering).
################################################################################
# Applies ComBat (sva) to purified beta values to remove v1/v2 batch effect,
# then runs t-SNE on all samples combined.
# Compare with split-by-array approach from script 15.
################################################################################

if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(Rtsne)
  library(sva)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})

set.seed(42)

message("================================================================")
message("DECODE: t-SNE with ComBat Batch Correction")
message("================================================================\n")

# ============================================================================
# Load metadata
# ============================================================================

meta <- read.csv("data/processed/siewert_cohort_results.csv", stringsAsFactors = FALSE)
gbm_ids <- meta$Idat[grepl("Glioblastoma, IDH", meta$Max_family)]

assign_gbm_class <- function(subclass) {
  cls <- rep(NA_character_, length(subclass))
  cls[grepl("RTK1", subclass)] <- "RTK I"
  cls[grepl("RTK2", subclass)] <- "RTK II"
  cls[grepl("mesenchymal", subclass, ignore.case = TRUE)] <- "MES"
  cls[grepl("primitive neuronal", subclass, ignore.case = TRUE)] <- "PNC"
  cls
}

# Array version
sentrix <- as.numeric(gsub("_.*", "", gbm_ids))
array_version <- ifelse(sentrix >= 207000000000, "EPIC v2", "EPIC v1")
names(array_version) <- gbm_ids

# Purified classifier labels (DECODE-fraction re-run, Heidelberg mnp.v12epicv2)
# Canonical calls from data/processed/siewert_purified_decode3_calls.csv
source("scripts/_purification_config.R")
pur_calls <- read.csv(PURIFIED_CALLS, stringsAsFactors = FALSE)
pur_labels <- data.frame(
  Idat = pur_calls$Sample,
  pur_gbm_class = NA_character_,
  stringsAsFactors = FALSE
)
pur_labels$pur_gbm_class[pur_calls$pred_sub == "GBM_RTK1"] <- "RTK I"
pur_labels$pur_gbm_class[pur_calls$pred_sub == "GBM_RTK2"] <- "RTK II"
pur_labels$pur_gbm_class[pur_calls$pred_sub %in% c("GBM_MES_TYP", "GBM_MES_ATYP")] <- "MES"
pur_labels$pur_gbm_class[pur_calls$pred_sub == "GBM_PNC"] <- "PNC"
pur_labels$pur_gbm_class[is.na(pur_labels$pur_gbm_class)] <- "non-GBM"

message(sprintf("GBM IDH-wt samples: %d (v1: %d, v2: %d)",
                length(gbm_ids),
                sum(array_version == "EPIC v1"),
                sum(array_version == "EPIC v2")))

# ============================================================================
# Load purified beta matrix, filter probes
# ============================================================================

message("\nLoading purified beta matrix...")
# DECODE-fraction purification (3 components: immune/neuron/glia) -- 2026-06-12
# Same cache as script 15e: everything down to the t-SNE loads two 2.2 GB
# matrices and runs Rtsne three times, but the plots need only the three small
# coordinate frames. Reused unless an input is newer or DECODE_REBUILD_TSNE is set.
TSNE_CACHE <- "data/processed/tsne_combat_cache.rds"
REBUILD_TSNE <- as.logical(Sys.getenv("DECODE_REBUILD_TSNE", "FALSE"))
cache_ok <- !REBUILD_TSNE && file.exists(TSNE_CACHE) &&
  file.mtime(TSNE_CACHE) > max(file.mtime(PURIFIED_BETA),
                               file.mtime("/mnt/volume-hel1-1/siewert/betamatrix.rds"),
                               file.mtime("data/processed/siewert_cohort_results.csv"))
if (cache_ok) {
  message("Reusing cached t-SNE coordinates (", TSNE_CACHE, ")")
  .cache <- readRDS(TSNE_CACHE)
  df_orig <- .cache$orig; df_nobc <- .cache$nobc; df_comb <- .cache$comb
} else {

beta_pur <- readRDS(PURIFIED_BETA)
gbm_in_beta <- intersect(gbm_ids, colnames(beta_pur))
beta_pur <- beta_pur[, gbm_in_beta]

# Remove sex chromosomes
anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
xy_probes <- rownames(anno[anno$chr %in% c("chrX", "chrY"), ])
beta_pur <- beta_pur[setdiff(rownames(beta_pur), xy_probes), ]
rm(anno); gc(verbose = FALSE)
message(sprintf("  Autosomal probes: %d", nrow(beta_pur)))

# Select top 10k most variable probes
message("Selecting top 10,000 most variable CpGs...")
n_probes <- nrow(beta_pur)
chunk_size <- 50000
vars <- numeric(n_probes)
for (i in seq(1, n_probes, chunk_size)) {
  end <- min(i + chunk_size - 1, n_probes)
  chunk <- beta_pur[i:end, ]
  for (r in seq_len(nrow(chunk))) {
    na_idx <- is.na(chunk[r, ])
    if (any(na_idx)) chunk[r, na_idx] <- mean(chunk[r, !na_idx], na.rm = TRUE)
  }
  vars[i:end] <- apply(chunk, 1, var, na.rm = TRUE)
}
top_idx <- order(vars, decreasing = TRUE)[1:10000]
beta_pur <- beta_pur[top_idx, ]
beta_pur[is.na(beta_pur)] <- 0.5
rm(vars); gc(verbose = FALSE)
message(sprintf("  Selected 10,000 probes"))

# ============================================================================
# Apply ComBat batch correction
# ============================================================================

message("\nApplying ComBat batch correction (v1 vs v2)...")
batch <- array_version[gbm_in_beta]

# ComBat needs a model matrix — we preserve biological variation (subtype)
orig_subclass <- meta$Max_subclass[match(gbm_in_beta, meta$Idat)]
orig_class <- assign_gbm_class(orig_subclass)
orig_class[is.na(orig_class)] <- "unknown"

mod <- model.matrix(~ factor(orig_class))

# ComBat on M-values (logit-transformed beta) is recommended
# But beta values work too for visualization purposes
beta_combat <- ComBat(
  dat = beta_pur,
  batch = batch,
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
)

# Clip to [0, 1] range (ComBat can produce values outside)
beta_combat[beta_combat < 0] <- 0
beta_combat[beta_combat > 1] <- 1

message("  ComBat done.")

# ============================================================================
# Run t-SNEs: uncorrected vs ComBat-corrected
# ============================================================================

message("\nRunning t-SNE on purified (no batch correction)...")
tsne_nocombat <- Rtsne(t(beta_pur), dims = 2, perplexity = 30, max_iter = 1000,
                       check_duplicates = FALSE, pca = TRUE, pca_center = TRUE,
                       num_threads = 4)

message("Running t-SNE on purified + ComBat corrected...")
tsne_combat <- Rtsne(t(beta_combat), dims = 2, perplexity = 30, max_iter = 1000,
                     check_duplicates = FALSE, pca = TRUE, pca_center = TRUE,
                     num_threads = 4)

rm(beta_pur, beta_combat); gc(verbose = FALSE)

# ============================================================================
# Also load and run original (unpurified) for 3-panel comparison
# ============================================================================

message("Loading original beta matrix for comparison...")
beta_orig <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
beta_orig <- beta_orig[setdiff(rownames(beta_orig), xy_probes), gbm_in_beta]

# Top 10k most variable
vars_orig <- numeric(nrow(beta_orig))
for (i in seq(1, nrow(beta_orig), chunk_size)) {
  end <- min(i + chunk_size - 1, nrow(beta_orig))
  chunk <- beta_orig[i:end, ]
  for (r in seq_len(nrow(chunk))) {
    na_idx <- is.na(chunk[r, ])
    if (any(na_idx)) chunk[r, na_idx] <- mean(chunk[r, !na_idx], na.rm = TRUE)
  }
  vars_orig[i:end] <- apply(chunk, 1, var, na.rm = TRUE)
}
top_orig <- order(vars_orig, decreasing = TRUE)[1:10000]
beta_orig <- beta_orig[top_orig, ]
beta_orig[is.na(beta_orig)] <- 0.5
rm(vars_orig); gc(verbose = FALSE)

message("Running t-SNE on original beta values...")
tsne_orig <- Rtsne(t(beta_orig), dims = 2, perplexity = 30, max_iter = 1000,
                   check_duplicates = FALSE, pca = TRUE, pca_center = TRUE,
                   num_threads = 4)
rm(beta_orig); gc(verbose = FALSE)

# ============================================================================
# Build plot data
# ============================================================================

orig_class_vec <- assign_gbm_class(meta$Max_subclass[match(gbm_in_beta, meta$Idat)])
# Die eine Probe ohne konfidenten Subtyp heisst "Other", nicht NA, damit die
# Legende im Manuskript keine NA-Kategorie zeigt (wie in Script 15e).
orig_class_vec[is.na(orig_class_vec)] <- "Other"
pur_class_vec <- pur_labels$pur_gbm_class[match(gbm_in_beta, pur_labels$Idat)]
array_vec <- array_version[gbm_in_beta]
was_mes <- grepl("mesenchymal", meta$Max_subclass[match(gbm_in_beta, meta$Idat)], ignore.case = TRUE)

make_df <- function(tsne, stage) {
  data.frame(tSNE1 = tsne$Y[, 1], tSNE2 = tsne$Y[, 2],
             Original_Subtype = orig_class_vec,
             Purified_Subtype = pur_class_vec,
             Array = array_vec, was_MES = was_mes,
             Stage = stage)
}

df_orig <- make_df(tsne_orig, "1. Original")
df_nobc <- make_df(tsne_nocombat, "2. Purified (no batch corr.)")
df_comb <- make_df(tsne_combat, "3. Purified + ComBat")

saveRDS(list(orig = df_orig, nobc = df_nobc, comb = df_comb), TSNE_CACHE)
message("Cached t-SNE coordinates: ", TSNE_CACHE)
}

df_all <- rbind(df_orig, df_nobc, df_comb)
df_all$Stage <- factor(df_all$Stage)

# ============================================================================
# Plots
# ============================================================================

message("\n--- Generating Figures ---")

subtype_colors <- c("RTK I" = "#E41A1C", "RTK II" = "#377EB8",
                     "MES" = "#4DAF4A", "PNC" = "#FF7F00", "non-GBM" = "#999999",
                     "Other" = "#CCCCCC")
array_colors <- c("EPIC v1" = "#E66101", "EPIC v2" = "#5E3C99")

facet_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank())

# --- 3-panel: ALL panels coloured by the ORIGINAL classification ---
# Until 2026-08-30 the purified panels were coloured by the PURIFIED call. Green
# then vanished after purification, but that is relabelling (165 MES become
# RTK II), not geometry, so the panel could not show what it was cited for --
# that the MES cluster dissolves. Keeping the original label fixed makes the
# dispersal of the green points the actual evidence. The purified-label view is
# already given separately by the migration panel below.
df_all$Subtype <- df_all$Original_Subtype

p1 <- ggplot(df_all, aes(x = tSNE1, y = tSNE2, color = Subtype)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~Stage, scales = "free", nrow = 1) +
  scale_color_manual(values = subtype_colors) +
  facet_theme +
  labs(title = "t-SNE: Original vs Purified vs Purified+ComBat",
       subtitle = "Original colored by original class (with MES), purified panels by purified class (MES dissolved)",
       color = "GBM Subtype")
ggsave("figures/tsne_combat_comparison_subtype.pdf", p1, width = 18, height = 6)

# Manuscript version: the ORIGINAL embedding only. The three-panel version
# showed the purified state twice (with and without ComBat) and the ComBat one
# again in the migration panel below, so the figure carried the same state three
# times. What the manuscript needs from this plot is the "before": that MES
# forms a cluster of its own at all. The batch-correction control is stated in
# the text instead, and the "after" is the migration panel, which shows the same
# embedding plus where each formerly-MES sample went.
p1_orig <- ggplot(df_orig, aes(x = tSNE1, y = tSNE2, color = Original_Subtype)) +
  geom_point(size = 1.1, alpha = 0.75) +
  scale_color_manual(values = subtype_colors) +
  facet_theme +
  labs(x = "tSNE1", y = "tSNE2", color = "Original subtype")
# Written only as a diagnostic single panel; the manuscript version is the
# before/after pair built at the end of this script.
ggsave("figures/tsne_original_subtype_single.png", p1_orig, width = 9, height = 8, dpi = 150)
message("  Saved: figures/tsne_combat_comparison_subtype.pdf")

# --- 3-panel: colored by array version (batch effect diagnostic) ---
p2 <- ggplot(df_all, aes(x = tSNE1, y = tSNE2, color = Array)) +
  geom_point(size = 1, alpha = 0.6) +
  facet_wrap(~Stage, scales = "free", nrow = 1) +
  scale_color_manual(values = array_colors) +
  facet_theme +
  labs(title = "t-SNE: Batch Effect Diagnostic (EPIC v1 vs v2)",
       subtitle = "Does ComBat remove the v1/v2 split after purification?",
       color = "Array Version")
ggsave("figures/tsne_combat_comparison_array.pdf", p2, width = 18, height = 6)
message("  Saved: figures/tsne_combat_comparison_array.pdf")

# --- 3-panel: all colored by ORIGINAL subtype (track MES migration) ---
p3 <- ggplot(df_all, aes(x = tSNE1, y = tSNE2, color = Original_Subtype)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~Stage, scales = "free", nrow = 1) +
  scale_color_manual(values = subtype_colors) +
  facet_theme +
  labs(title = "t-SNE: MES Migration Tracking",
       subtitle = "All panels colored by original classification - where do MES samples move?",
       color = "Original Subtype")
ggsave("figures/tsne_combat_comparison_migration.pdf", p3, width = 18, height = 6)
message("  Saved: figures/tsne_combat_comparison_migration.pdf")

# --- ComBat only: MES migration (triangles) ---
# Shape is its own aesthetic with its own legend. Until 2026-09-02 the two point
# layers shared a single colour legend while only one of them fixed the shape,
# so each subtype inherited the shape of whichever layer it happened to occur
# in: PNC exists only among the non-MES samples and therefore appeared as a
# circle while the rest were triangles, encoding nothing. The explanation of the
# shapes also lived in the subtitle, which the manuscript version strips.
df_comb$Origin <- factor(ifelse(df_comb$was_MES, "originally MES", "other samples"),
                         levels = c("originally MES", "other samples"))
df_p4 <- df_comb[order(df_comb$Origin, decreasing = TRUE), ]   # triangles drawn last

p4 <- ggplot(df_p4, aes(x = tSNE1, y = tSNE2, color = Purified_Subtype,
                        shape = Origin, size = Origin, alpha = Origin)) +
  geom_point() +
  scale_color_manual(values = subtype_colors) +
  scale_shape_manual(values = c("originally MES" = 17, "other samples" = 16)) +
  scale_size_manual(values = c("originally MES" = 2.5, "other samples" = 1.2), guide = "none") +
  scale_alpha_manual(values = c("originally MES" = 0.85, "other samples" = 0.3), guide = "none") +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.6, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(color = "grey30", size = 2.6, alpha = 1))) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(0.1, "cm")) +
  labs(title = "Purified + ComBat: MES Migration",
       subtitle = "All points coloured by post-purification class.",
       color = "Purified subtype", shape = NULL)
ggsave("figures/tsne_combat_mes_migration.pdf", p4, width = 9, height = 8)

# --- Manuscript panel b -----------------------------------------------------
# Colour carries the ORIGINAL subtype here, exactly as in panel a, so one colour
# means one thing across the pair. Colouring by the purified call (as p4 above
# does, deliberately, as a diagnostic) makes the manuscript panel partly
# circular: the ex-MES samples sit in the RTK II region and are drawn RTK II blue
# *because* they were reclassified there, so the colour cannot show anything
# else. The informative half is the geometry, and the purified call belongs on
# the shape, where it adds something the colour does not already say.
df_p4b <- df_comb
df_p4b$Moved <- factor(ifelse(as.character(df_p4b$Original_Subtype) !=
                              as.character(df_p4b$Purified_Subtype),
                              "reclassified", "unchanged"),
                       levels = c("reclassified", "unchanged"))
df_p4b <- df_p4b[order(df_p4b$Moved, decreasing = TRUE), ]   # triangles drawn last

p4b <- ggplot(df_p4b, aes(x = tSNE1, y = tSNE2, color = Original_Subtype,
                          shape = Moved, size = Moved, alpha = Moved)) +
  geom_point() +
  scale_color_manual(values = subtype_colors) +
  scale_shape_manual(values = c("reclassified" = 17, "unchanged" = 16)) +
  scale_size_manual(values = c("reclassified" = 2.2, "unchanged" = 1.3), guide = "none") +
  scale_alpha_manual(values = c("reclassified" = 0.9, "unchanged" = 0.45), guide = "none") +
  guides(color = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.6, alpha = 1)),
         shape = guide_legend(order = 2, override.aes = list(color = "grey30", size = 2.6, alpha = 1))) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(0.1, "cm")) +
  labs(x = "tSNE1", y = "tSNE2",
       color = "Original subtype", shape = "Purified call")
# --- Manuscript figure: before and after side by side -----------------------
# The two embeddings are a before/after pair and were previously stacked as two
# separate images, which is the hardest arrangement to compare across. They are
# now one figure with the panel letters drawn in, so the caption still refers to
# c and d. Both keep their own legend, because they are coloured by different
# things: the left by the ORIGINAL subtype, the right by the PURIFIED one.
tag_theme <- theme(plot.tag = element_text(size = 15, face = "bold"),
                   plot.title = element_text(size = 12, hjust = 0.5),
                   legend.position = "bottom", legend.box = "vertical",
                   legend.key.size = unit(0.35, "cm"),
                   legend.text = element_text(size = 9),
                   legend.title = element_text(size = 9))

p_before <- p1_orig + tag_theme + labs(tag = "a", title = "Before purification") +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 2.4, alpha = 1)))
p_after  <- p4b + labs(title = NULL, subtitle = NULL) + tag_theme +
  labs(tag = "b", title = "After purification and ComBat") +
  guides(color = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 16, size = 2.4, alpha = 1)),
         shape = guide_legend(order = 2, nrow = 1, override.aes = list(color = "grey30", size = 2.4, alpha = 1)))

ggsave("manuscript/figures/Fig5cd_tsne_before_after.png",
       gridExtra::arrangeGrob(p_before, p_after, ncol = 2),
       width = 14, height = 7, dpi = 300)   # >>> Manuscript Figure 7c,d
message("  Saved: manuscript/figures/Fig5cd_tsne_before_after.png (before/after pair)")
message("  Saved: figures/tsne_combat_mes_migration.pdf")

message("\n================================================================")
message("Done!")
message("================================================================")
