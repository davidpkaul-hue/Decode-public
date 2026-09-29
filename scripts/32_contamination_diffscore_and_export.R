# ============================================================================
# Script 32: Contamination experiment - differentiation-score check + export
# ----------------------------------------------------------------------------
# (A) DIFFERENTIATION CHECK: apply the ITHresolveGBM differentiation score to
#     the contaminated RTK mixtures. If tumor biology is untouched (only the
#     immune admixture grows), the malignant-renormalised differentiation score
#     stays FLAT while the classifier label flips -> proof that contamination
#     changes only the label, not the tumor cell state.
# (B) EXPORT (Stage 2): build full-CpG (721k) contaminated MethylSets (v1/v2)
#     at a titration of immune fractions for the real Heidelberg v12.8 call
#     (to be run by Christin), plus a manifest.
# ============================================================================

suppressPackageStartupMessages({
  library(EpiDISH); library(minfi); library(ggplot2); library(matrixStats)
})
if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

fig_dir <- "figures"; dir.create(fig_dir, showWarnings = FALSE)
export_dir <- "/mnt/volume-hel1-1/siewert/contamination_methylsets"
dir.create(export_dir, showWarnings = FALSE)

# --- MES-like immune microenvironment profile (same weights as script 31) ----
mes_w <- c(Mono = 10.2, Microglia = 14.5, CD8T = 4.4, CD4T = 1.8,
           Bcell = 1.5, NK = 0.3, Neu = 3.4); mes_w <- mes_w / sum(mes_w)
imm <- readRDS("data/processed/immune_profiles_fullcpg.rds")[, names(mes_w)]
immune_profile <- as.numeric(imm %*% mes_w); names(immune_profile) <- rownames(imm)

# --- samples to contaminate (the 12 from script 31) + metadata --------------
cres <- read.csv("data/processed/contamination_results.csv")
rtk_ids <- unique(cres$sample)
orig_sub <- setNames(cres$orig_subtype[match(rtk_ids, cres$sample)], rtk_ids)

message("Loading cohort beta matrix ...")
beta <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
beta <- beta[, rtk_ids, drop = FALSE]              # only the 12 samples we need
gc(verbose = FALSE)

# ============================================================================
# (A) Differentiation-score titration
# ============================================================================
message("\n(A) ITHresolveGBM differentiation score across contamination ...")
prof <- read.table("data/external/ithresolve/profiles_full",
                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
ith_probes <- intersect(prof$probe, intersect(rownames(beta), names(immune_profile)))
ref_ith <- as.matrix(prof[match(ith_probes, prof$probe),
                          c("immune","neuron","glia","V4","V5","V6")])
rownames(ref_ith) <- ith_probes
ip_ith <- immune_profile[ith_probes]
b_ith  <- beta[ith_probes, , drop = FALSE]
# drop probes with any NA (ref, immune profile, or any sample) to avoid singular fits
ok <- complete.cases(ref_ith) & !is.na(ip_ith) & rowSums(is.na(b_ith)) == 0
message(sprintf("  ITH probes usable: %d of %d", sum(ok), length(ok)))
ref_ith <- ref_ith[ok, ]; ip_ith <- ip_ith[ok]; b_ith <- b_ith[ok, , drop = FALSE]
ith_probes <- ith_probes[ok]

grid_A <- seq(0, 0.6, by = 0.05)
# build all mixtures (probes x [samples*fracs]) then one deconvolution
mix_cols <- list(); labs <- list()
for (id in rtk_ids) for (f in grid_A) {
  mix_cols[[paste0(id, "|", f)]] <- (1 - f) * b_ith[, id] + f * ip_ith
  labs[[paste0(id, "|", f)]] <- data.frame(sample = id, orig = orig_sub[id], f = f)
}
mixM <- do.call(cbind, mix_cols)
est <- epidish(mixM, ref_ith, method = "RPC")$estF     # cols: immune neuron glia V4 V5 V6
lab <- do.call(rbind, labs)
# Differentiation score, malignant-renormalised: (Stem - Diff1 - Diff2) / malig_sum.
# The division is load-bearing, not cosmetic. Adding immune signal shrinks all three
# malignant components proportionally, so the RAW score would drift towards zero for
# purely arithmetic reasons and the "differentiation state is untouched" argument
# would be circular. Renormalising to the malignant compartment makes the score
# report its composition rather than its size. Note this differs from the raw score
# in script 12; the Methods section states both.
stem <- est[, "V4"]; diff1 <- est[, "V6"]; diff2 <- est[, "V5"]
malig <- stem + diff1 + diff2
lab$diff_score <- (stem - diff1 - diff2) / malig
lab$decode_immune <- cres$decode_immune[match(paste0(lab$sample, lab$f),
                                              paste0(cres$sample, cres$f))]
write.csv(lab, "data/processed/contamination_diffscore.csv", row.names = FALSE)

# flatness stats: change in diff score from f=0 to max f, per sample
d0 <- lab$diff_score[lab$f == 0]; dmax <- lab$diff_score[lab$f == max(grid_A)]
message(sprintf("  Differentiation score change f=0 -> f=%.2f: median |Δ| = %.3f (range %.3f-%.3f)",
        max(grid_A), median(abs(dmax - d0)), min(abs(dmax - d0)), max(abs(dmax - d0))))

pA <- ggplot(lab, aes(100 * f, diff_score, group = sample, colour = orig)) +
  geom_line(alpha = 0.7) + geom_point(size = 1) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey50") +
  labs(x = "Added immune fraction f (%)",
       y = "Malignant-renormalised differentiation score\n(stem - differentiated)",
       colour = "Original",
       title = "Tumor differentiation state is unchanged by contamination",
       subtitle = "Flat lines = malignant cell state untouched while the classifier label flips") +
  theme_bw()
ggsave(file.path(fig_dir, "contamination_diffscore.png"), pA, width = 8, height = 5, dpi = 150)
message("  Wrote figures/contamination_diffscore.png")

# ============================================================================
# (B) Full-CpG contaminated MethylSets for the real classifier (Christin)
# ============================================================================
message("\n(B) Building full-CpG contaminated MethylSets for Christin ...")
common <- intersect(rownames(beta), names(immune_profile))
beta <- beta[common, , drop = FALSE]
ip   <- immune_profile[common]

grid_B <- seq(0, 0.5, by = 0.1)                    # 6 levels incl. f=0 control
cols <- list(); manifest <- list()
for (id in rtk_ids) for (f in grid_B) {
  nm <- sprintf("%s__f%02d", id, round(f * 100))
  cols[[nm]] <- (1 - f) * beta[, id] + f * ip
  manifest[[nm]] <- data.frame(
    pseudo_sample = nm, orig_sample = id, orig_subtype = orig_sub[id], frac = f,
    expected_decode_immune = cres$decode_immune[match(paste0(id, f),
                                                paste0(cres$sample, cres$f))])
}
mixfull <- do.call(cbind, cols)
manifest <- do.call(rbind, manifest)
write.csv(manifest, "data/processed/contamination_manifest.csv", row.names = FALSE)
message("  ", ncol(mixfull), " pseudo-samples (", length(rtk_ids), " RTK x ",
        length(grid_B), " fractions)")

sentrix <- as.numeric(sub("_.*", "", colnames(mixfull)))
is_v2 <- sentrix >= 207000000000
TOTAL <- 10000
for (arr in c("v1", "v2")) {
  idx <- if (arr == "v1") !is_v2 else is_v2
  if (!any(idx)) { message("  ", arr, ": 0 samples, skipped"); next }
  b <- mixfull[, idx, drop = FALSE]
  b[is.na(b)] <- 0.5
  ann <- if (arr == "v1") c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b4.hg19")
         else            c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
  mset <- MethylSet(Meth = round(b * TOTAL), Unmeth = round((1 - b) * TOTAL),
                    annotation = ann,
                    colData = DataFrame(data.frame(Sample_Name = colnames(b),
                                                   row.names = colnames(b))))
  outpath <- file.path(export_dir, sprintf("methylset_contaminated_%s.rds", arr))
  saveRDS(mset, outpath)
  message(sprintf("  Saved %s: %d samples (%.0f MB)", outpath, ncol(b), file.size(outpath)/1e6))
  rm(b, mset); gc(verbose = FALSE)
}
message("\nDONE. Manifest: data/processed/contamination_manifest.csv")
message("Hand the two MethylSets in ", export_dir, " to Christin for mnp.v12epicv2.")
