#!/usr/bin/env Rscript
################################################################################
# DECODE: In-Silico Validation
################################################################################
# Purpose: Validate deconvolution accuracy using held-out control samples
#          Creates synthetic mixtures with known proportions and tests recovery
# Input:   data/processed/reference_oneVsAll.rds
#          data/processed/control_samples.rds
# Output:  data/processed/validation_results.rds
#          figures/validation_*.pdf
# Author:  davidpkaul-hue
# Date:    2026-02-02
################################################################################

suppressPackageStartupMessages({
  library(EpiDISH)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

message("======================================")
message("DECODE: In-Silico Validation")
message("======================================\n")

# ============================================================================
# Configuration
# ============================================================================

reference_file <- "data/processed/reference_oneVsAll.rds"
control_file <- "data/processed/control_samples.rds"
output_dir <- "data/processed"
figures_dir <- "figures"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# Load Data
# ============================================================================

message("Loading reference and control samples...")

ref <- readRDS(reference_file)
ctrl <- readRDS(control_file)

signature <- ref$signature
ctrl_beta <- ctrl$beta
ctrl_pheno <- ctrl$phenotype
cell_types <- colnames(signature)
n_ct <- length(cell_types)

message(sprintf("  Reference: %d probes x %d cell types", nrow(signature), ncol(signature)))
message(sprintf("  Controls: %d probes x %d samples (held-out)", nrow(ctrl_beta), ncol(ctrl_beta)))

# Get control profiles per cell type
ctrl_profiles <- list()
for (i in 1:ncol(ctrl_beta)) {
  ct <- ctrl_pheno[i]
  ctrl_profiles[[ct]] <- ctrl_beta[, i]
}

# ============================================================================
# Test 1: Control Sample Deconvolution
# ============================================================================

message("\n--- Test 1: Control Sample Deconvolution ---")
message("Testing if held-out controls are correctly identified...\n")

result_ctrl <- epidish(beta.m = ctrl_beta, ref.m = signature, method = "RPC")
props_ctrl <- round(result_ctrl$estF * 100, 1)

ctrl_results <- data.frame(
  CellType = ctrl_pheno,
  SelfPrediction = sapply(1:length(ctrl_pheno), function(i) props_ctrl[i, ctrl_pheno[i]]),
  TopPredicted = colnames(props_ctrl)[apply(props_ctrl, 1, which.max)]
)
ctrl_results$Correct <- ctrl_results$CellType == ctrl_results$TopPredicted

message("Results:")
for (i in 1:nrow(ctrl_results)) {
  status <- ifelse(ctrl_results$Correct[i], "OK", "WRONG")
  message(sprintf("  %s: %.1f%% [%s]",
                  ctrl_results$CellType[i],
                  ctrl_results$SelfPrediction[i],
                  status))
}

message(sprintf("\nAccuracy: %d/%d (%.0f%%)",
                sum(ctrl_results$Correct), nrow(ctrl_results),
                100 * mean(ctrl_results$Correct)))
message(sprintf("Mean self-prediction: %.1f%%", mean(ctrl_results$SelfPrediction)))

# ============================================================================
# Test 2: Equal Mixture (all cell types at 7.7%)
# ============================================================================

message("\n--- Test 2: Equal Mixture (7.7% each) ---")

# Create equal mixture and add to a matrix with a dummy for RPC stability
mix_equal <- rep(0, nrow(ctrl_beta))
for (ct in cell_types) {
  mix_equal <- mix_equal + (1/n_ct) * ctrl_profiles[[ct]]
}

# RPC needs multiple samples, so we add a small perturbation as second sample
mix_equal_matrix <- cbind(mix_equal, mix_equal + rnorm(length(mix_equal), 0, 0.001))
colnames(mix_equal_matrix) <- c("Equal_Mix", "Equal_Mix_dup")

result_equal <- epidish(beta.m = mix_equal_matrix, ref.m = signature, method = "RPC")
props_equal <- round(result_equal$estF * 100, 1)[1, ]  # Take first sample

message("\nExpected vs Observed:")
equal_results <- data.frame(
  CellType = cell_types,
  Expected = round(100/n_ct, 1),
  Observed = as.numeric(props_equal[cell_types]),
  Diff = as.numeric(props_equal[cell_types]) - round(100/n_ct, 1)
)

for (i in 1:nrow(equal_results)) {
  message(sprintf("  %s: Exp %.1f%% Obs %.1f%% (Diff %+.1f)",
                  equal_results$CellType[i],
                  equal_results$Expected[i],
                  equal_results$Observed[i],
                  equal_results$Diff[i]))
}

message(sprintf("\nMean absolute error: %.1f%%", mean(abs(equal_results$Diff))))

# ============================================================================
# Test 3: 50% Dominant Mixtures
# ============================================================================

message("\n--- Test 3: 50% Dominant Mixtures ---")

dominant_results <- data.frame(
  CellType = character(),
  Expected = numeric(),
  Observed = numeric(),
  Error = numeric(),
  TopPredicted = character(),
  Correct = logical()
)

mixtures_50 <- list()

for (dominant in cell_types) {
  # Create mixture: 50% dominant, rest equal
  mix <- 0.5 * ctrl_profiles[[dominant]]
  other_prop <- 0.5 / (n_ct - 1)
  for (other in cell_types[cell_types != dominant]) {
    mix <- mix + other_prop * ctrl_profiles[[other]]
  }
  mixtures_50[[dominant]] <- mix
}

mix_matrix_50 <- do.call(cbind, mixtures_50)
colnames(mix_matrix_50) <- paste0(cell_types, "_50pct")

result_50 <- epidish(beta.m = mix_matrix_50, ref.m = signature, method = "RPC")
props_50 <- round(result_50$estF * 100, 1)

message("\nResults:")
for (dominant in cell_types) {
  mix_name <- paste0(dominant, "_50pct")
  obs <- props_50[mix_name, dominant]
  err <- obs - 50
  top_pred <- names(which.max(props_50[mix_name, ]))
  correct <- top_pred == dominant
  status <- ifelse(correct, "OK", paste0("WRONG (", top_pred, ")"))

  message(sprintf("  %s: Exp 50.0%% Obs %.1f%% (Err %+.1f) [%s]",
                  dominant, obs, err, status))

  dominant_results <- rbind(dominant_results, data.frame(
    CellType = dominant,
    Expected = 50,
    Observed = obs,
    Error = err,
    TopPredicted = top_pred,
    Correct = correct
  ))
}

message(sprintf("\nAccuracy: %d/%d (%.0f%%)",
                sum(dominant_results$Correct), nrow(dominant_results),
                100 * mean(dominant_results$Correct)))
message(sprintf("Mean absolute error: %.1f%%", mean(abs(dominant_results$Error))))

# ============================================================================
# Test 4: Custom 3-Way Mixtures (Oligo/Astro/GBM)
# ============================================================================

message("\n--- Test 4: 3-Way Mixtures (Oligo/Astro/GBM) ---")

three_way_configs <- list(
  "Equal_33_33_33" = c(Oligodendrocytes = 1/3, Astrocytes = 1/3, GBM_malignant = 1/3),
  "GBM50_Oligo25_Astro25" = c(Oligodendrocytes = 0.25, Astrocytes = 0.25, GBM_malignant = 0.50),
  "GBM10_Oligo45_Astro45" = c(Oligodendrocytes = 0.45, Astrocytes = 0.45, GBM_malignant = 0.10)
)

three_way_results <- list()

for (config_name in names(three_way_configs)) {
  props_config <- three_way_configs[[config_name]]

  mix <- rep(0, nrow(ctrl_beta))
  for (ct in names(props_config)) {
    mix <- mix + props_config[ct] * ctrl_profiles[[ct]]
  }

  # RPC needs multiple samples
  mix_matrix <- cbind(mix, mix + rnorm(length(mix), 0, 0.001))
  result <- epidish(beta.m = mix_matrix, ref.m = signature, method = "RPC")
  obs <- round(result$estF * 100, 1)[1, ]

  three_way_results[[config_name]] <- list(
    expected = props_config * 100,
    observed = obs[names(props_config)]
  )

  message(sprintf("\n%s:", config_name))
  for (ct in names(props_config)) {
    exp <- props_config[ct] * 100
    obs_val <- obs[ct]
    message(sprintf("  %s: Exp %.1f%% Obs %.1f%% (Diff %+.1f)",
                    ct, exp, obs_val, obs_val - exp))
  }
}

# ============================================================================
# Generate Plots
# ============================================================================

message("\n--- Generating Plots ---")

# Plot 1: Control sample self-prediction
p1 <- ggplot(ctrl_results, aes(x = reorder(CellType, -SelfPrediction), y = SelfPrediction,
                                fill = Correct)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Control Sample Validation",
       subtitle = "Self-prediction percentage for held-out samples",
       x = "Cell Type", y = "Self-Prediction (%)", fill = "Correct") +
  scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
  ylim(0, 105)

ggsave(file.path(figures_dir, "validation_control_samples.pdf"), p1, width = 10, height = 6)
message("  Saved: validation_control_samples.pdf")

# Plot 2: 50% dominant mixture accuracy
p2 <- ggplot(dominant_results, aes(x = reorder(CellType, -Observed), y = Observed,
                                    fill = Correct)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "50% Dominant Mixture Validation",
       subtitle = "Observed percentage when expected is 50%",
       x = "Dominant Cell Type", y = "Observed (%)", fill = "Correct") +
  scale_fill_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
  ylim(0, 70)

ggsave(file.path(figures_dir, "validation_50pct_mixtures.pdf"), p2, width = 10, height = 6)
message("  Saved: validation_50pct_mixtures.pdf")

# Plot 3: Equal mixture deviation
p3 <- ggplot(equal_results, aes(x = reorder(CellType, abs(Diff)), y = Diff,
                                 fill = Diff > 0)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Equal Mixture Validation (7.7% each)",
       subtitle = "Deviation from expected proportion",
       x = "Cell Type", y = "Deviation (%)", fill = "Over-estimated") +
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C"))

ggsave(file.path(figures_dir, "validation_equal_mixture.pdf"), p3, width = 10, height = 6)
message("  Saved: validation_equal_mixture.pdf")

# ============================================================================
# Save Results
# ============================================================================

message("\n--- Saving Results ---")

validation_results <- list(
  control_samples = list(
    proportions = props_ctrl,
    summary = ctrl_results,
    accuracy = mean(ctrl_results$Correct),
    mean_self_prediction = mean(ctrl_results$SelfPrediction)
  ),
  equal_mixture = list(
    proportions = props_equal,
    summary = equal_results,
    mean_absolute_error = mean(abs(equal_results$Diff))
  ),
  dominant_50pct = list(
    proportions = props_50,
    summary = dominant_results,
    accuracy = mean(dominant_results$Correct),
    mean_absolute_error = mean(abs(dominant_results$Error))
  ),
  three_way = three_way_results
)

saveRDS(validation_results, file.path(output_dir, "validation_results.rds"))
message(sprintf("  Saved: %s", file.path(output_dir, "validation_results.rds")))

# ============================================================================
# Summary
# ============================================================================

message("\n======================================")
message("VALIDATION SUMMARY")
message("======================================")
message(sprintf("Control samples:     %d/%d correct (%.0f%%), mean self-pred %.1f%%",
                sum(ctrl_results$Correct), nrow(ctrl_results),
                100 * mean(ctrl_results$Correct),
                mean(ctrl_results$SelfPrediction)))
message(sprintf("Equal mixture:       Mean absolute error %.1f%%",
                mean(abs(equal_results$Diff))))
message(sprintf("50%% dominant:        %d/%d correct (%.0f%%), mean abs error %.1f%%",
                sum(dominant_results$Correct), nrow(dominant_results),
                100 * mean(dominant_results$Correct),
                mean(abs(dominant_results$Error))))

message("\nDone!")
