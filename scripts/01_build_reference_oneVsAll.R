#!/usr/bin/env Rscript
################################################################################
# DECODE: Build Reference Matrix with oneVsAllTtest Marker Selection
################################################################################
# Purpose: Build reference matrix using oneVsAllTtest for marker selection
#          - Keeps Oligodendrocytes and Astrocytes separate
#          - Holds out 1 sample per cell type as control
#          - Selects 500 markers per cell type
# Output:  data/processed/reference_oneVsAll.rds (reference)
#          data/processed/control_samples.rds (held-out controls)
# Author:  davidpkaul-hue
# Date:    2026-02-02
################################################################################

suppressPackageStartupMessages({
  library(sesame)
  library(EpiDISH)
})

message("======================================")
message("DECODE: Build Reference (oneVsAllTtest)")
message("======================================\n")

# ============================================================================
# Configuration
# ============================================================================

base_dir <- "bjarne/Reference_IDATS"
output_ref <- "data/processed/reference_oneVsAll.rds"
output_ctrl <- "data/processed/control_samples.rds"
n_markers_per_ct <- 500

dir.create(dirname(output_ref), showWarnings = FALSE, recursive = TRUE)

# Reference directories - Oligo and Astro SEPARATE
ref_dirs <- list(
  # Blood cells
  Bcell = file.path(base_dir, "FlowSorted.Blood.EPIC/Bcell"),
  CD4T = file.path(base_dir, "FlowSorted.Blood.EPIC/CD4T"),
  CD8T = file.path(base_dir, "FlowSorted.Blood.EPIC/CD8T"),
  Mono = file.path(base_dir, "FlowSorted.Blood.EPIC/Mono"),
  NK = file.path(base_dir, "FlowSorted.Blood.EPIC/NK"),
  Neu = file.path(base_dir, "FlowSorted.Blood.EPIC/Neu"),

  # Brain cells - SEPARATE Oligo and Astro
  Neurons = file.path(base_dir, "Hannon_et_al/Neurons"),
  Oligodendrocytes = file.path(base_dir, "Hannon_et_al/Oligodendrocytes"),
  Astrocytes = file.path(base_dir, "Hannon_et_al/Astrocytes"),
  Microglia = file.path(base_dir, "Hannon_et_al/Microglia"),

  # Stromal
  Endothelial = file.path(base_dir, "Moss_et_al"),
  Fibroblasts = file.path(base_dir, "Fibroblasts_Dura"),

  # Tumor
  GBM_malignant = file.path(base_dir, "Shen_et_al/diagnostic_score_gbm_12.8")
)

# ============================================================================
# Load all IDATs
# ============================================================================

message("Loading reference IDAT files...\n")

load_idats <- function(idat_dir, cell_type) {
  message(sprintf("  Loading %s...", cell_type))

  if (!dir.exists(idat_dir)) {
    warning(sprintf("    Directory not found: %s", idat_dir))
    return(NULL)
  }

  prefixes <- sesame::searchIDATprefixes(idat_dir)
  if (length(prefixes) == 0) {
    warning(sprintf("    No IDATs found"))
    return(NULL)
  }

  beta_list <- lapply(prefixes, function(pfx) {
    tryCatch({
      sesame::getBetas(sesame::prepSesame(sesame::readIDATpair(pfx), "QCDPB"))
    }, error = function(e) NULL)
  })

  beta_list <- beta_list[!sapply(beta_list, is.null)]
  if (length(beta_list) == 0) return(NULL)

  beta_mat <- do.call(cbind, beta_list)
  message(sprintf("    Loaded %d samples", ncol(beta_mat)))
  return(beta_mat)
}

all_betas <- list()
all_phenotypes <- c()

for (ct in names(ref_dirs)) {
  beta_mat <- load_idats(ref_dirs[[ct]], ct)
  if (!is.null(beta_mat)) {
    all_betas[[ct]] <- beta_mat
    all_phenotypes <- c(all_phenotypes, rep(ct, ncol(beta_mat)))
  }
}

# ============================================================================
# Combine and find common probes
# ============================================================================

message("\nCombining samples...")
common_probes <- Reduce(intersect, lapply(all_betas, rownames))
message(sprintf("  Common probes: %d", length(common_probes)))

combined_beta <- do.call(cbind, lapply(all_betas, function(x) x[common_probes, ]))
message(sprintf("  Combined: %d probes x %d samples", nrow(combined_beta), ncol(combined_beta)))

# ============================================================================
# Hold out 1 sample per cell type as control
# ============================================================================

message("\nHolding out 1 sample per cell type as control...")

set.seed(42)  # Reproducible
cell_types <- unique(all_phenotypes)

ctrl_idx <- c()
for (ct in cell_types) {
  ct_idx <- which(all_phenotypes == ct)
  hold_out <- sample(ct_idx, 1)
  ctrl_idx <- c(ctrl_idx, hold_out)
  message(sprintf("  %s: holding out sample %d", ct, hold_out))
}

# Split into reference and control
ref_idx <- setdiff(1:ncol(combined_beta), ctrl_idx)

ref_beta <- combined_beta[, ref_idx]
ref_pheno <- all_phenotypes[ref_idx]

ctrl_beta <- combined_beta[, ctrl_idx]
ctrl_pheno <- all_phenotypes[ctrl_idx]

message(sprintf("\n  Reference samples: %d", ncol(ref_beta)))
message(sprintf("  Control samples: %d", ncol(ctrl_beta)))

# ============================================================================
# oneVsAllTtest: Select markers
# ============================================================================

message("\nSelecting markers with oneVsAllTtest (500 per cell type)...\n")

oneVsAllTtest <- function(beta, phenotype, n_markers = 500, p_threshold = 0.01) {
  cell_types <- unique(phenotype)
  all_markers <- list()
  marker_info <- list()

  for (ct in cell_types) {
    message(sprintf("  Processing %s...", ct))

    idx_ct <- which(phenotype == ct)
    idx_other <- which(phenotype != ct)

    if (length(idx_ct) < 2 || length(idx_other) < 2) {
      warning(sprintf("    Skipping %s: not enough samples", ct))
      next
    }

    # T-test for each probe
    pvals <- apply(beta, 1, function(x) {
      tryCatch({
        t.test(x[idx_ct], x[idx_other])$p.value
      }, error = function(e) NA)
    })

    # Effect size (mean difference)
    effect <- rowMeans(beta[, idx_ct, drop = FALSE], na.rm = TRUE) -
              rowMeans(beta[, idx_other, drop = FALSE], na.rm = TRUE)

    # Select significant probes
    valid <- !is.na(pvals) & pvals < p_threshold
    if (sum(valid) == 0) {
      warning(sprintf("    No significant markers for %s", ct))
      next
    }

    # Sort by absolute effect size
    sig_probes <- names(pvals)[valid]
    sig_effects <- abs(effect[sig_probes])
    top_markers <- names(sort(sig_effects, decreasing = TRUE))[1:min(n_markers, length(sig_effects))]

    all_markers[[ct]] <- top_markers
    marker_info[[ct]] <- data.frame(
      probe = top_markers,
      effect = effect[top_markers],
      pval = pvals[top_markers]
    )

    message(sprintf("    Selected %d markers", length(top_markers)))
  }

  unique_markers <- unique(unlist(all_markers))
  message(sprintf("\n  Total unique markers: %d", length(unique_markers)))

  return(list(
    markers = unique_markers,
    markers_per_ct = all_markers,
    marker_info = marker_info
  ))
}

marker_result <- oneVsAllTtest(ref_beta, ref_pheno, n_markers = n_markers_per_ct)
markers <- marker_result$markers

# ============================================================================
# Build signature matrix (mean per cell type, marker probes only)
# ============================================================================

message("\nBuilding signature matrix...")

# Subset to markers
ref_beta_markers <- ref_beta[markers, ]
ref_beta_markers <- ref_beta_markers[complete.cases(ref_beta_markers), ]
message(sprintf("  Using %d markers (after removing NAs)", nrow(ref_beta_markers)))

# Calculate mean per cell type
signature <- sapply(cell_types, function(ct) {
  idx <- which(ref_pheno == ct)
  if (length(idx) == 1) {
    ref_beta_markers[, idx]
  } else {
    rowMeans(ref_beta_markers[, idx], na.rm = TRUE)
  }
})

message(sprintf("  Signature matrix: %d probes x %d cell types",
                nrow(signature), ncol(signature)))

# ============================================================================
# Save reference
# ============================================================================

message("\nSaving reference matrix...")

reference <- list(
  signature = signature,
  beta = ref_beta_markers,
  phenotype = ref_pheno,
  cell_types = cell_types,
  markers = markers,
  markers_per_ct = marker_result$markers_per_ct,
  marker_info = marker_result$marker_info,
  n_markers_per_ct = n_markers_per_ct,
  method = "oneVsAllTtest"
)

saveRDS(reference, output_ref)
message(sprintf("  Saved reference: %s", output_ref))

# ============================================================================
# Save control samples
# ============================================================================

message("\nSaving control samples...")

# Also subset controls to marker probes
ctrl_beta_markers <- ctrl_beta[markers, ]
ctrl_beta_markers <- ctrl_beta_markers[rownames(ref_beta_markers), ]

controls <- list(
  beta = ctrl_beta_markers,
  phenotype = ctrl_pheno,
  cell_types = cell_types
)

saveRDS(controls, output_ctrl)
message(sprintf("  Saved controls: %s", output_ctrl))

# ============================================================================
# Summary
# ============================================================================

message("\n======================================")
message("Summary")
message("======================================")
message(sprintf("Reference samples: %d", ncol(ref_beta_markers)))
message(sprintf("Control samples: %d (1 per cell type)", ncol(ctrl_beta_markers)))
message(sprintf("Marker probes: %d", nrow(signature)))
message(sprintf("Cell types: %d", length(cell_types)))

message("\nCell type distribution (reference):")
for (ct in cell_types) {
  n <- sum(ref_pheno == ct)
  message(sprintf("  - %s: %d samples", ct, n))
}

message("\nMarkers per cell type:")
for (ct in names(marker_result$markers_per_ct)) {
  n <- length(marker_result$markers_per_ct[[ct]])
  message(sprintf("  - %s: %d markers", ct, n))
}

message("\nDone!")
