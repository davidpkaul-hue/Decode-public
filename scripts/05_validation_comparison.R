#!/usr/bin/env Rscript
################################################################################
# DECODE: Validation Comparison - Bjarne vs Our Pipeline
################################################################################
# Purpose: Compare two deconvolution approaches:
#   1. Bjarne's pipeline (methylDeConv, 11 cell types)
#   2. Our pipeline (EpiDISH, 13 cell types)
#
# Validation methods:
#   A. Real samples: Bjarne's 10 GBM samples (ground truth = Bjarne's results)
#   B. In-silico mixtures: Known proportions from reference samples
#
# Output:
#   - data/processed/validation_comparison.rds
#   - figures/validation_*.pdf
################################################################################

suppressPackageStartupMessages({
  library(sesame)
  library(EpiDISH)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

message("==========================================")
message("DECODE: Validation Comparison")
message("==========================================\n")

# ============================================================================
# Configuration
# ============================================================================

our_ref_file <- "data/processed/reference_oneVsAll.rds"
bjarne_csv_dir <- "bjarne/cell_fractions_bjarne"
bjarne_idat_dir <- "bjarne/query_IDATS"
output_file <- "data/processed/validation_comparison.rds"
fig_dir <- "figures"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load Our Reference
# ============================================================================

message("Loading our reference matrix...")
if (!file.exists(our_ref_file)) {
  stop("Reference file not found. Run 01_build_reference_oneVsAll.R first.")
}

our_ref <- readRDS(our_ref_file)
our_sig <- our_ref$signature
message(sprintf("  Our signature: %d probes x %d cell types",
                nrow(our_sig), ncol(our_sig)))
message(sprintf("  Cell types: %s", paste(colnames(our_sig), collapse = ", ")))

# ============================================================================
# PART A: Real Sample Comparison (Bjarne's 10 GBM samples)
# ============================================================================

message("\n==========================================")
message("PART A: Real Sample Comparison")
message("==========================================\n")

# Load Bjarne's results
message("Loading Bjarne's results...")
csv_files <- list.files(bjarne_csv_dir, pattern = "*.csv", full.names = TRUE)

if (length(csv_files) == 0) {
  message("  No Bjarne CSV files found. Skipping real sample comparison.")
  bjarne_results <- NULL
  our_results_real <- NULL
} else {
  bjarne_results <- do.call(rbind, lapply(csv_files, function(f) {
    sample_id <- gsub("_cellfractions.csv", "", basename(f))
    df <- read.csv(f)
    df$Sample <- sample_id
    df
  }))
  message(sprintf("  Loaded %d samples from Bjarne", nrow(bjarne_results)))

  # Process same IDATs with our pipeline
  message("\nProcessing IDATs with our pipeline...")
  prefixes <- sesame::searchIDATprefixes(bjarne_idat_dir)

  beta_query <- do.call(cbind, lapply(prefixes, function(pfx) {
    getBetas(prepSesame(readIDATpair(pfx), "QCDPB"))
  }))
  message(sprintf("  Processed %d samples", ncol(beta_query)))

  # Align probes
  common_probes <- intersect(rownames(beta_query), rownames(our_sig))
  message(sprintf("  Common probes: %d", length(common_probes)))

  beta_query <- beta_query[common_probes, ]
  sig_common <- our_sig[common_probes, ]

  # Run our deconvolution
  message("Running EpiDISH-RPC deconvolution...")
  res <- epidish(beta.m = beta_query, ref.m = sig_common, method = "RPC")
  our_results_real <- as.data.frame(res$estF)
  our_results_real$Sample <- colnames(beta_query)

  # Calculate comparison metrics
  message("\n--- Comparison: Bjarne vs Our Pipeline ---\n")

  # Map cell type names
  cell_map <- c(
    "malignant" = "GBM_malignant",
    "astrocytes" = "Astrocytes",
    "neurons" = "Neurons",
    "oligodendrocytes" = "Oligodendrocytes",
    "microglia" = "Microglia",
    "Bcell" = "Bcell",
    "CD4T" = "CD4T",
    "CD8T" = "CD8T",
    "Mono" = "Mono",
    "NK" = "NK",
    "Neu" = "Neu",
    "endothel" = "Endothelial"
  )

  comparison_real <- data.frame()

  for (bj_name in names(cell_map)) {
    our_name <- cell_map[bj_name]

    if (bj_name %in% names(bjarne_results) && our_name %in% names(our_results_real)) {
      bj_vals <- c()
      our_vals <- c()

      for (sample in bjarne_results$Sample) {
        bj_vals <- c(bj_vals, bjarne_results[bjarne_results$Sample == sample, bj_name])
        our_vals <- c(our_vals, our_results_real[our_results_real$Sample == sample, our_name])
      }

      cor_val <- cor(bj_vals, our_vals, use = "complete.obs")
      mean_diff <- mean(abs(bj_vals - our_vals), na.rm = TRUE)

      comparison_real <- rbind(comparison_real, data.frame(
        CellType = bj_name,
        Correlation = cor_val,
        MeanAbsDiff = mean_diff,
        stringsAsFactors = FALSE
      ))
    }
  }

  comparison_real <- comparison_real[order(-comparison_real$Correlation), ]

  message("Correlation per cell type:")
  for (i in 1:nrow(comparison_real)) {
    message(sprintf("  %-20s r=%.3f  MAD=%.1f%%",
                    comparison_real$CellType[i],
                    comparison_real$Correlation[i],
                    comparison_real$MeanAbsDiff[i] * 100))
  }

  overall_cor <- mean(comparison_real$Correlation, na.rm = TRUE)
  overall_mad <- mean(comparison_real$MeanAbsDiff, na.rm = TRUE)
  message(sprintf("\n  OVERALL: r=%.3f, MAD=%.1f%%", overall_cor, overall_mad * 100))
}

# ============================================================================
# PART B: In-Silico Validation
# ============================================================================

message("\n==========================================")
message("PART B: In-Silico Validation")
message("==========================================\n")

# Load HELD-OUT control samples (NOT used in reference building!)
# These samples were excluded during reference matrix construction
# to provide an independent validation set.
ctrl_file <- "data/processed/control_samples.rds"

if (!file.exists(ctrl_file)) {
  message("Control samples not found. Skipping in-silico validation.")
  insilico_our <- NULL
  insilico_bjarne <- NULL
} else {
  ctrl <- readRDS(ctrl_file)
  ctrl_beta <- ctrl$beta
  ctrl_celltypes <- ctrl$cell_types  # Note: cell_types not celltypes

  message(sprintf("Loaded %d HELD-OUT control samples (independent from reference!)", ncol(ctrl_beta)))
  message(sprintf("Cell types: %s", paste(unique(ctrl_celltypes), collapse = ", ")))
  message("  NOTE: These samples were excluded during reference building (1 per cell type)")
  message("        In-silico mixtures use ONLY these held-out samples, not reference data.")

  # ----------------------------------------------------------------------------
  # Create in-silico mixtures (structured, same as script 04)
  # ----------------------------------------------------------------------------

  message("\nCreating structured in-silico mixtures...")

  # Get common probes between control samples and reference
  common_probes <- intersect(rownames(ctrl_beta), rownames(our_sig))
  ctrl_beta_common <- ctrl_beta[common_probes, ]

  cell_types <- unique(ctrl_celltypes)
  n_ct <- length(cell_types)

  # Create cell type profiles (one beta vector per cell type)
  ctrl_profiles <- list()
  for (ct in cell_types) {
    ct_idx <- which(ctrl_celltypes == ct)
    ctrl_profiles[[ct]] <- ctrl_beta_common[, ct_idx[1]]  # Use held-out sample
  }

  mixtures <- list()
  true_fractions <- list()
  mixture_names <- c()

  # --- Test 1: Equal mixture (all cell types at 1/n_ct each) ---
  message(sprintf("  Creating equal mixture (%.1f%% each)...", 100/n_ct))

  mix_equal <- rep(0, length(common_probes))
  props_equal <- rep(1/n_ct, n_ct)
  names(props_equal) <- cell_types

  for (ct in cell_types) {
    mix_equal <- mix_equal + (1/n_ct) * ctrl_profiles[[ct]]
  }

  mixtures[["Equal_Mix"]] <- mix_equal
  true_fractions[["Equal_Mix"]] <- props_equal
  mixture_names <- c(mixture_names, "Equal_Mix")

  # --- Test 2: 50% Dominant mixtures (one per cell type) ---
  message(sprintf("  Creating %d dominant mixtures (50%% + rest equal)...", n_ct))

  for (dominant in cell_types) {
    mix_name <- paste0(dominant, "_50pct")

    # 50% dominant, rest equally distributed
    mix <- 0.5 * ctrl_profiles[[dominant]]
    other_prop <- 0.5 / (n_ct - 1)

    props <- rep(other_prop, n_ct)
    names(props) <- cell_types
    props[dominant] <- 0.5

    for (other in cell_types[cell_types != dominant]) {
      mix <- mix + other_prop * ctrl_profiles[[other]]
    }

    mixtures[[mix_name]] <- mix
    true_fractions[[mix_name]] <- props
    mixture_names <- c(mixture_names, mix_name)
  }

  # --- Test 3: Tumor gradient (0/20/40/60/80/100% GBM, rest equally distributed) ---
  message("  Creating 6 tumor-gradient mixtures (0/20/40/60/80/100% GBM)...")

  tumor_levels <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
  non_tumor    <- cell_types[cell_types != "GBM_malignant"]

  for (tf in tumor_levels) {
    mix_name <- sprintf("Tumor_%dpct", round(tf * 100))
    other_prop <- (1 - tf) / length(non_tumor)

    mix <- tf * ctrl_profiles[["GBM_malignant"]]
    for (other in non_tumor) {
      mix <- mix + other_prop * ctrl_profiles[[other]]
    }

    props <- setNames(rep(other_prop, n_ct), cell_types)
    props["GBM_malignant"] <- tf

    mixtures[[mix_name]] <- mix
    true_fractions[[mix_name]] <- props
    mixture_names <- c(mixture_names, mix_name)
  }

  # --- Test 4: Three-way Oligo/Astro/GBM mixtures (biologically realistic GBM tissue) ---
  message("  Creating 3 three-way Oligo/Astro/GBM mixtures...")

  three_way_configs <- list(
    "GliaTumor_33_33_33"  = c(Oligodendrocytes = 1/3,  Astrocytes = 1/3,  GBM_malignant = 1/3),
    "GliaTumor_25_25_50"  = c(Oligodendrocytes = 0.25, Astrocytes = 0.25, GBM_malignant = 0.50),
    "GliaTumor_45_45_10"  = c(Oligodendrocytes = 0.45, Astrocytes = 0.45, GBM_malignant = 0.10)
  )

  for (cfg_name in names(three_way_configs)) {
    cfg <- three_way_configs[[cfg_name]]
    mix <- rep(0, length(common_probes))
    for (ct in names(cfg)) mix <- mix + cfg[ct] * ctrl_profiles[[ct]]

    props <- setNames(rep(0, n_ct), cell_types)
    props[names(cfg)] <- cfg

    mixtures[[cfg_name]] <- mix
    true_fractions[[cfg_name]] <- props
    mixture_names <- c(mixture_names, cfg_name)
  }

  n_mixtures <- length(mixtures)

  # Convert to matrix
  mix_beta <- do.call(cbind, mixtures)
  colnames(mix_beta) <- mixture_names
  rownames(mix_beta) <- common_probes

  # Create true fractions dataframe
  true_frac_df <- do.call(rbind, lapply(mixture_names, function(m) {
    data.frame(
      Sample = m,
      as.list(true_fractions[[m]]),
      stringsAsFactors = FALSE
    )
  }))

  message(sprintf("  Total: %d mixtures (%d probes)", n_mixtures, nrow(mix_beta)))
  message(sprintf("    - 1 equal mixture"))
  message(sprintf("    - %d x 50%% dominant", n_ct))
  message(sprintf("    - %d x 70%% dominant", n_ct))

  # ----------------------------------------------------------------------------
  # Run OUR pipeline on in-silico mixtures
  # ----------------------------------------------------------------------------

  message("\nRunning OUR pipeline (EpiDISH-RPC) on in-silico mixtures...")

  sig_insilico <- our_sig[common_probes, ]
  res_our <- epidish(beta.m = mix_beta, ref.m = sig_insilico, method = "RPC")
  pred_our <- as.data.frame(res_our$estF)
  pred_our$Sample <- rownames(pred_our)

  # Calculate errors for our pipeline
  errors_our <- data.frame()

  for (sample in mixture_names) {
    for (ct in cell_types) {
      true_val <- true_fractions[[sample]][ct]
      pred_val <- pred_our[pred_our$Sample == sample, ct]

      if (length(pred_val) > 0 && !is.na(pred_val) && !is.na(true_val)) {
        errors_our <- rbind(errors_our, data.frame(
          Sample = sample,
          CellType = ct,
          True = true_val,
          Predicted = pred_val,
          Error = pred_val - true_val,
          AbsError = abs(pred_val - true_val),
          Pipeline = "Our (EpiDISH)",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  message(sprintf("  Our pipeline - Mean Absolute Error: %.2f%%",
                  mean(errors_our$AbsError) * 100))

  # ----------------------------------------------------------------------------
  # Run BJARNE-style pipeline on in-silico mixtures
  # ----------------------------------------------------------------------------

  message("\nRunning BJARNE-style pipeline (methylDeConv-RPC) on in-silico mixtures...")

  # Bjarne uses 11 cell types (no Endothelial, Fibroblasts)
  bjarne_celltypes <- setdiff(colnames(our_sig), c("Endothelial", "Fibroblasts"))
  sig_bjarne <- our_sig[common_probes, bjarne_celltypes]

  # Bjarne uses methylDeConv::RPC, but EpiDISH-RPC is algorithmically equivalent
  res_bjarne <- epidish(beta.m = mix_beta, ref.m = sig_bjarne, method = "RPC")
  pred_bjarne <- as.data.frame(res_bjarne$estF)
  pred_bjarne$Sample <- rownames(pred_bjarne)

  # Calculate errors for Bjarne-style pipeline
  errors_bjarne <- data.frame()

  for (sample in mixture_names) {
    for (ct in bjarne_celltypes) {
      true_val <- true_fractions[[sample]][ct]
      pred_val <- pred_bjarne[pred_bjarne$Sample == sample, ct]

      if (length(pred_val) > 0 && !is.na(pred_val) && !is.na(true_val)) {
        errors_bjarne <- rbind(errors_bjarne, data.frame(
          Sample = sample,
          CellType = ct,
          True = true_val,
          Predicted = pred_val,
          Error = pred_val - true_val,
          AbsError = abs(pred_val - true_val),
          Pipeline = "Bjarne-style (11 types)",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  message(sprintf("  Bjarne-style pipeline - Mean Absolute Error: %.2f%%",
                  mean(errors_bjarne$AbsError) * 100))

  # Combine errors
  all_errors <- rbind(errors_our, errors_bjarne)

  insilico_our <- pred_our
  insilico_bjarne <- pred_bjarne
}

# ============================================================================
# Summary and Visualization
# ============================================================================

message("\n==========================================")
message("Summary")
message("==========================================\n")

if (!is.null(bjarne_results)) {
  message("PART A - Real Sample Comparison:")
  message(sprintf("  Samples: %d", nrow(bjarne_results)))
  message(sprintf("  Overall correlation: r = %.3f", overall_cor))
  message(sprintf("  Mean absolute difference: %.1f%%", overall_mad * 100))
}

if (exists("all_errors")) {
  message("\nPART B - In-Silico Validation:")

  summary_by_pipeline <- all_errors %>%
    group_by(Pipeline) %>%
    summarise(
      MAE = mean(AbsError) * 100,
      Correlation = cor(True, Predicted),
      .groups = "drop"
    )

  for (i in 1:nrow(summary_by_pipeline)) {
    message(sprintf("  %s: MAE=%.2f%%, r=%.3f",
                    summary_by_pipeline$Pipeline[i],
                    summary_by_pipeline$MAE[i],
                    summary_by_pipeline$Correlation[i]))
  }

  # Per cell type comparison
  message("\nPer cell type (Our pipeline):")
  summary_by_ct_our <- all_errors %>%
    filter(Pipeline == "Our (EpiDISH)") %>%
    group_by(CellType) %>%
    summarise(
      MAE = mean(AbsError) * 100,
      Correlation = cor(True, Predicted),
      .groups = "drop"
    ) %>%
    arrange(MAE)

  for (i in 1:nrow(summary_by_ct_our)) {
    message(sprintf("  %-20s MAE=%.2f%%, r=%.3f",
                    summary_by_ct_our$CellType[i],
                    summary_by_ct_our$MAE[i],
                    summary_by_ct_our$Correlation[i]))
  }
}

# ============================================================================
# Visualizations
# ============================================================================

message("\nGenerating visualizations...")

# Plot 1: Real sample comparison (if available)
if (!is.null(bjarne_results) && !is.null(our_results_real)) {

  # Create comparison data for plotting
  plot_data_real <- data.frame()

  for (bj_name in names(cell_map)) {
    our_name <- cell_map[bj_name]

    if (bj_name %in% names(bjarne_results) && our_name %in% names(our_results_real)) {
      for (sample in bjarne_results$Sample) {
        bj_val <- bjarne_results[bjarne_results$Sample == sample, bj_name]
        our_val <- our_results_real[our_results_real$Sample == sample, our_name]

        plot_data_real <- rbind(plot_data_real, data.frame(
          Sample = sample,
          CellType = bj_name,
          Bjarne = bj_val,
          Our = our_val,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  p1 <- ggplot(plot_data_real, aes(x = Bjarne, y = Our)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    facet_wrap(~CellType, scales = "free") +
    theme_bw() +
    labs(
      title = "Real Sample Comparison: Bjarne vs Our Pipeline",
      subtitle = sprintf("n=%d samples, overall r=%.3f", nrow(bjarne_results), overall_cor),
      x = "Bjarne's Pipeline (fraction)",
      y = "Our Pipeline (fraction)"
    )

  ggsave(file.path(fig_dir, "validation_real_samples.pdf"), p1, width = 12, height = 10)
  message("  Saved: validation_real_samples.pdf")
}

# Plot 2: In-silico validation (if available)
if (exists("all_errors")) {

  p2 <- ggplot(all_errors, aes(x = True, y = Predicted, color = Pipeline)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    facet_wrap(~CellType, scales = "free") +
    theme_bw() +
    labs(
      title = "In-Silico Validation: True vs Predicted Fractions",
      subtitle = sprintf("n=%d mixtures", n_mixtures),
      x = "True Fraction",
      y = "Predicted Fraction"
    ) +
    scale_color_manual(values = c("Our (EpiDISH)" = "#3498DB",
                                   "Bjarne-style (11 types)" = "#E74C3C"))

  ggsave(file.path(fig_dir, "validation_insilico.pdf"), p2, width = 14, height = 10)
  message("  Saved: validation_insilico.pdf")

  # Plot 3: Error distribution by pipeline
  p3 <- ggplot(all_errors, aes(x = Pipeline, y = AbsError * 100, fill = Pipeline)) +
    geom_boxplot(alpha = 0.7) +
    theme_bw() +
    labs(
      title = "In-Silico Validation: Error Distribution",
      x = "",
      y = "Absolute Error (%)"
    ) +
    scale_fill_manual(values = c("Our (EpiDISH)" = "#3498DB",
                                  "Bjarne-style (11 types)" = "#E74C3C")) +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "validation_error_boxplot.pdf"), p3, width = 6, height = 5)
  message("  Saved: validation_error_boxplot.pdf")
}

# ============================================================================
# Save Results
# ============================================================================

message("\nSaving results...")

results <- list(
  # Part A: Real samples
  bjarne_results = bjarne_results,
  our_results_real = our_results_real,
  comparison_real = if(exists("comparison_real")) comparison_real else NULL,

  # Part B: In-silico
  insilico_our = if(exists("insilico_our")) insilico_our else NULL,
  insilico_bjarne = if(exists("insilico_bjarne")) insilico_bjarne else NULL,
  insilico_errors = if(exists("all_errors")) all_errors else NULL,
  insilico_true = if(exists("true_frac_df")) true_frac_df else NULL,

  # Metadata
  n_real_samples = if(!is.null(bjarne_results)) nrow(bjarne_results) else 0,
  n_insilico_mixtures = if(exists("n_mixtures")) n_mixtures else 0
)

saveRDS(results, output_file)
message(sprintf("  Saved: %s", output_file))

# ============================================================================
# Export detailed CSV files
# ============================================================================

message("\nExporting detailed CSV files...")

# 1. Real sample comparison (Bjarne vs Our)
if (!is.null(bjarne_results) && !is.null(our_results_real)) {

  detailed_real <- data.frame()

  for (bj_name in names(cell_map)) {
    our_name <- cell_map[bj_name]

    if (bj_name %in% names(bjarne_results) && our_name %in% names(our_results_real)) {
      for (sample in bjarne_results$Sample) {
        bj_val <- bjarne_results[bjarne_results$Sample == sample, bj_name]
        our_val <- our_results_real[our_results_real$Sample == sample, our_name]

        detailed_real <- rbind(detailed_real, data.frame(
          Sample = sample,
          CellType = bj_name,
          Bjarne_Fraction = round(bj_val, 4),
          Our_Fraction = round(our_val, 4),
          Difference = round(our_val - bj_val, 4),
          Abs_Difference = round(abs(our_val - bj_val), 4),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Sort by sample and cell type
  detailed_real <- detailed_real[order(detailed_real$Sample, detailed_real$CellType), ]

  write.csv(detailed_real, "data/processed/validation_real_samples_detailed.csv", row.names = FALSE)
  message("  Saved: data/processed/validation_real_samples_detailed.csv")

  # Also create wide format (one row per sample)
  wide_bjarne <- bjarne_results
  names(wide_bjarne) <- paste0(names(wide_bjarne), "_Bjarne")
  names(wide_bjarne)[names(wide_bjarne) == "Sample_Bjarne"] <- "Sample"

  wide_our <- our_results_real[, c("Sample", setdiff(names(our_results_real), "Sample"))]
  names(wide_our)[-1] <- paste0(names(wide_our)[-1], "_Our")

  wide_combined <- merge(wide_bjarne, wide_our, by = "Sample")
  write.csv(wide_combined, "data/processed/validation_real_samples_wide.csv", row.names = FALSE)
  message("  Saved: data/processed/validation_real_samples_wide.csv")
}

# 2. In-silico validation details
if (exists("all_errors") && nrow(all_errors) > 0) {

  # Detailed errors
  insilico_detailed <- all_errors %>%
    mutate(
      True_Percent = round(True * 100, 2),
      Predicted_Percent = round(Predicted * 100, 2),
      Error_Percent = round(Error * 100, 2),
      AbsError_Percent = round(AbsError * 100, 2)
    ) %>%
    select(Sample, CellType, Pipeline, True_Percent, Predicted_Percent, Error_Percent, AbsError_Percent) %>%
    arrange(Pipeline, Sample, CellType)

  write.csv(insilico_detailed, "data/processed/validation_insilico_detailed.csv", row.names = FALSE)
  message("  Saved: data/processed/validation_insilico_detailed.csv")

  # True fractions used for mixtures
  true_frac_percent <- true_frac_df
  for (col in setdiff(names(true_frac_percent), "Sample")) {
    true_frac_percent[[col]] <- round(true_frac_percent[[col]] * 100, 2)
  }
  write.csv(true_frac_percent, "data/processed/validation_insilico_true_fractions.csv", row.names = FALSE)
  message("  Saved: data/processed/validation_insilico_true_fractions.csv")

  # Predicted fractions (Our pipeline)
  pred_our_percent <- insilico_our
  for (col in setdiff(names(pred_our_percent), "Sample")) {
    pred_our_percent[[col]] <- round(pred_our_percent[[col]] * 100, 2)
  }
  write.csv(pred_our_percent, "data/processed/validation_insilico_predicted_our.csv", row.names = FALSE)
  message("  Saved: data/processed/validation_insilico_predicted_our.csv")

  # Predicted fractions (Bjarne-style pipeline)
  pred_bjarne_percent <- insilico_bjarne
  for (col in setdiff(names(pred_bjarne_percent), "Sample")) {
    pred_bjarne_percent[[col]] <- round(pred_bjarne_percent[[col]] * 100, 2)
  }
  write.csv(pred_bjarne_percent, "data/processed/validation_insilico_predicted_bjarne.csv", row.names = FALSE)
  message("  Saved: data/processed/validation_insilico_predicted_bjarne.csv")
}

message("\n==========================================")
message("Done!")
message("==========================================")
