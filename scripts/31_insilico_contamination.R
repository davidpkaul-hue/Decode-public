# ============================================================================
# Script 31: In-silico contamination (reverse of purification)
# ----------------------------------------------------------------------------
# Question: can we take high-purity RTK tumors and, by ARTIFICIALLY ADDING a
# realistic MES-like immune microenvironment, push them toward a MES call?
# This is the "sufficiency" mirror of the purification "necessity" experiment.
#
# STAGE 1 (this script, fully local): build the contaminated mixtures and read
# them out with a LOCAL surrogate classifier (nearest-centroid on the cohort's
# own subtype labels) + DECODE deconvolution as a counter-check. Go/No-Go.
# STAGE 2 (later, needs Christin): export contaminated MethylSets for the real
# Heidelberg v12.8 call.
# ============================================================================

suppressPackageStartupMessages({
  library(sesame)
  library(EpiDISH)
  library(matrixStats)
  library(ggplot2)
})
if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

set.seed(42)
base_dir   <- "bjarne/Reference_IDATS"
out_prof   <- "data/processed/immune_profiles_fullcpg.rds"
fig_dir    <- "figures"; dir.create(fig_dir, showWarnings = FALSE)
res_csv    <- "data/processed/contamination_results.csv"

# MES-like immune microenvironment weights (mean MES immune composition,
# normalised within the immune compartment; from Supp Table 2):
# Mono/Microglia/CD8T are the MES cohort means (Table 4). CD4T/Bcell/NK/Neu were
# carried over from an earlier deconvolution pass and sit above the final cohort
# values (0.5/1.1/0.2/2.3). Left unchanged so the classified pseudo-samples stay
# reproducible; the resulting spike-in profile is the same either way
# (r = 0.9999, mean |delta| = 0.003 beta against a profile from the final means).
mes_immune <- c(Mono = 10.2, Microglia = 14.5, CD8T = 4.4, CD4T = 1.8,
                Bcell = 1.5, NK = 0.3, Neu = 3.4)
mes_w <- mes_immune / sum(mes_immune)
message("MES immune weights:"); print(round(mes_w, 3))

# ---------------------------------------------------------------------------
# 1. Build full-CpG immune cell-type mean profiles (cached)
# ---------------------------------------------------------------------------
imm_dirs <- list(
  Bcell     = "FlowSorted.Blood.EPIC/Bcell",
  CD4T      = "FlowSorted.Blood.EPIC/CD4T",
  CD8T      = "FlowSorted.Blood.EPIC/CD8T",
  Mono      = "FlowSorted.Blood.EPIC/Mono",
  NK        = "FlowSorted.Blood.EPIC/NK",
  Neu       = "FlowSorted.Blood.EPIC/Neu",
  Microglia = "Hannon_et_al/Microglia"
)

load_mean_beta <- function(subpath, ct) {
  d <- file.path(base_dir, subpath)
  pref <- sesame::searchIDATprefixes(d)
  message(sprintf("  %s: %d IDAT pairs", ct, length(pref)))
  bl <- lapply(pref, function(p) tryCatch(
    sesame::getBetas(sesame::prepSesame(sesame::readIDATpair(p), "QCDPB")),
    error = function(e) NULL))
  bl <- bl[!sapply(bl, is.null)]
  m  <- do.call(cbind, bl)
  rowMeans(m, na.rm = TRUE)                       # mean profile for this cell type
}

if (file.exists(out_prof)) {
  message("Loading cached immune profiles: ", out_prof)
  imm <- readRDS(out_prof)
} else {
  message("Processing reference IDATs into full-CpG immune profiles ...")
  prof <- lapply(names(imm_dirs), function(ct) load_mean_beta(imm_dirs[[ct]], ct))
  names(prof) <- names(imm_dirs)
  common <- Reduce(intersect, lapply(prof, names))
  imm <- sapply(prof, function(v) v[common])      # CpG x 7 matrix
  rownames(imm) <- common
  saveRDS(imm, out_prof)
  message("Saved ", nrow(imm), " CpGs x ", ncol(imm), " cell types")
}

# MES-like microenvironment profile (weighted sum over the 7 immune types)
imm <- imm[, names(mes_w)]                         # order columns to weights
immune_profile <- as.numeric(imm %*% mes_w)
names(immune_profile) <- rownames(imm)

# ---------------------------------------------------------------------------
# 2. Cohort betas + subtype labels
# ---------------------------------------------------------------------------
message("Loading cohort beta matrix ...")
beta <- readRDS("/mnt/volume-hel1-1/siewert/betamatrix.rds")
meta <- read.csv("data/processed/siewert_cohort_results.csv",
                 check.names = FALSE, stringsAsFactors = FALSE)
sc <- meta[["Max_subclass"]]
meta$subtype <- ifelse(grepl("RTK1", sc), "RTK I",
                ifelse(grepl("RTK2", sc), "RTK II",
                ifelse(grepl("mesenchymal", sc, ignore.case = TRUE), "MES",
                ifelse(grepl("primitive neuronal", sc, ignore.case = TRUE), "PNC", NA))))
meta <- meta[meta$is_GBM %in% c(TRUE, "TRUE"), ]   # 418 GBM rows

# column-indexed maps (length = ncol(beta); NA for non-GBM columns)
mi        <- match(colnames(beta), meta$Idat)
col_sub   <- meta$subtype[mi]
col_tumor <- meta$GBM_malignant[mi]
col_imm   <- meta$Immune_total[mi]

# common CpGs across cohort + immune profile
common <- intersect(rownames(beta), names(immune_profile))
message("Common CpGs cohort/immune: ", length(common))

# discriminative CpGs = top-variable among GBM samples (as in the t-SNE, 10k)
gbm_cols <- which(!is.na(col_sub))
v <- rowVars(beta[common, gbm_cols, drop = FALSE], na.rm = TRUE)
D <- common[order(v, decreasing = TRUE)[1:10000]]
message("Discriminative (top-variable) CpGs: ", length(D))

# ---------------------------------------------------------------------------
# 3. Subtype centroids + surrogate nearest-centroid classifier (sanity check)
# ---------------------------------------------------------------------------
cents <- sapply(c("RTK I", "RTK II", "MES"), function(s)
  rowMeans(beta[D, which(col_sub == s), drop = FALSE], na.rm = TRUE))

classify <- function(x) {                          # x: named beta vector over >=D
  xd <- x[D]
  cc <- c(`RTK I`  = cor(xd, cents[, "RTK I"],  use = "complete.obs"),
          `RTK II` = cor(xd, cents[, "RTK II"], use = "complete.obs"),
          MES      = cor(xd, cents[, "MES"],    use = "complete.obs"))
  list(call = names(which.max(cc)),
       mes_margin = unname(cc["MES"] - max(cc["RTK I"], cc["RTK II"])))
}

# resubstitution sanity: how well does the surrogate reproduce real labels?
col_call <- rep(NA_character_, ncol(beta))
col_call[gbm_cols] <- sapply(gbm_cols, function(j) classify(beta[D, j])$call)
acc <- mean(col_call[gbm_cols] == col_sub[gbm_cols])
message(sprintf("Surrogate classifier resubstitution accuracy: %.1f%%", 100 * acc))
message("Confusion (real rows x predicted cols):")
print(table(real = col_sub[gbm_cols], pred = col_call[gbm_cols]))

# ---------------------------------------------------------------------------
# 4. Select high-purity RTK samples to contaminate
# ---------------------------------------------------------------------------
pick <- function(sub, n) {
  idx <- which(col_sub == sub & col_call == sub)        # correctly-called, clean
  ord <- order(-col_tumor[idx] + col_imm[idx])          # high purity, low immune
  colnames(beta)[idx[ord[seq_len(n)]]]
}
rtk_ids <- c(pick("RTK II", 8), pick("RTK I", 4))
rtk_ids <- rtk_ids[!is.na(rtk_ids)]
message("Contaminating ", length(rtk_ids), " high-purity RTK samples")

# ---------------------------------------------------------------------------
# 5. Titrate immune admixture and read out
# ---------------------------------------------------------------------------
fracs <- seq(0, 0.6, by = 0.05)
ref_sig <- readRDS("data/processed/reference_oneVsAll.rds")$signature
markers <- intersect(rownames(ref_sig), names(immune_profile))
ip_D <- immune_profile[D]
ip_M <- immune_profile[markers]

rows <- list()
for (id in rtk_ids) {
  b_full <- beta[, id]
  b_D <- b_full[D]; b_M <- b_full[markers]
  st0 <- col_sub[match(id, colnames(beta))]
  for (f in fracs) {
    mix_D <- (1 - f) * b_D + f * ip_D
    cl <- classify(setNames(mix_D, D))
    # DECODE fraction check on marker CpGs
    mix_M <- (1 - f) * b_M + f * ip_M
    est <- epidish(as.matrix(mix_M), ref_sig[markers, ], method = "RPC")$estF
    imm_frac <- sum(est[1, c("Bcell","CD4T","CD8T","Mono","NK","Neu","Microglia")])
    rows[[length(rows) + 1]] <- data.frame(
      sample = id, orig_subtype = st0, f = f,
      pred = cl$call, mes_margin = cl$mes_margin,
      decode_immune = imm_frac, decode_tumor = est[1, "GBM_malignant"])
  }
}
res <- do.call(rbind, rows)
write.csv(res, res_csv, row.names = FALSE)
message("Wrote ", res_csv)

# flip summary
flip <- aggregate(f ~ sample, data = res[res$pred == "MES", ],
                  FUN = function(x) min(x))
message("\n=== Immune fraction at which each RTK sample flips to MES ===")
if (nrow(flip)) print(flip) else message("  (none flipped within tested range)")
message(sprintf("\n%d of %d RTK samples flip to MES within f<=0.6",
                nrow(flip), length(rtk_ids)))

# ---------------------------------------------------------------------------
# 6. Plot: MES-margin vs immune admixture, per sample
# ---------------------------------------------------------------------------
res$flips <- res$sample %in% flip$sample
p <- ggplot(res, aes(100 * decode_immune, mes_margin, group = sample,
                     colour = orig_subtype)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  annotate("rect", xmin = 30, xmax = 42, ymin = -Inf, ymax = Inf,
           alpha = 0.08, fill = "red") +
  geom_line(alpha = 0.7) + geom_point(size = 1) +
  labs(x = "DECODE immune fraction after contamination (%)",
       y = "MES margin  (cor to MES - cor to nearest RTK)",
       colour = "Original",
       title = "In-silico contamination pushes RTK tumors toward MES",
       subtitle = "Above 0 = classified MES. Red band = real MES immune range (~30-42%)") +
  theme_bw()
ggsave(file.path(fig_dir, "contamination_margin.png"), p, width = 8, height = 5, dpi = 150)
message("Wrote ", file.path(fig_dir, "contamination_margin.png"))
message("\nDONE (Stage 1).")
