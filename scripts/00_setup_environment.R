#!/usr/bin/env Rscript
################################################################################
# DECODE: Environment Setup Script
################################################################################
# Purpose: Initialize R environment, install required packages, and verify setup
# Author: Prof. Dr. David Kaul
# Date: 2026-01-23 (updated 2026-05-24: added EpiDISH/sesame/InfiniumPurify/GIMiCC etc.)
################################################################################

message("======================================")
message("DECODE Environment Setup")
message("======================================\n")

# ============================================================================
# 1. Check R Version
# ============================================================================

message("Checking R version...")
r_version <- R.version.string
message(paste("  R version:", r_version))

if (getRversion() < "4.3.0") {
  warning("R version 4.3.0 or higher is recommended!")
}

# ============================================================================
# 2. Initialize renv (if not already initialized)
# ============================================================================

message("\nChecking renv status...")

if (!requireNamespace("renv", quietly = TRUE)) {
  message("  Installing renv...")
  install.packages("renv")
}

library(renv)

# Initialize renv if not already done
if (!file.exists("renv.lock")) {
  message("  Initializing renv...")
  renv::init()
} else {
  message("  renv already initialized")
  message("  Run renv::restore() to install project packages")
}

# ============================================================================
# 3. Define Required Packages
# ============================================================================

required_packages <- c(
  # Bioconductor core
  "BiocManager",

  # Methylation processing / annotation
  "minfi",
  "sesame",                                              # IDAT -> beta values
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",       # EPIC v1 annotation
  "FlowSorted.Blood.EPIC",                               # Blood reference

  # Deconvolution / purification
  "EpiDISH",                                             # RPC deconvolution (core algorithm)
  "InfiniumPurify",                                      # Purification baseline

  # Genomics utilities
  "BiocParallel",
  "Biostrings",
  "GenomeInfoDb",
  "GenomicRanges",
  "SummarizedExperiment",
  "missMethyl",
  "sva",                                                 # ComBat batch correction

  # Data processing
  "data.table",
  "dplyr",
  "tidyr",
  "readr",
  "readxl",                                              # Excel metadata

  # Visualization
  "ggplot2",
  "RColorBrewer",
  "gridExtra",
  "png",
  "uwot",                                                # UMAP
  "Rtsne",                                               # t-SNE

  # Statistical analysis / modeling
  "MASS",
  "glmnet",

  # Reporting
  "knitr",
  "rmarkdown"
)

# Packages that need GitHub installation (not on CRAN/Bioconductor)
github_packages <- list(
)

# ============================================================================
# 4. Package Installation Function
# ============================================================================

# Bioconductor packages (installed via BiocManager rather than install.packages)
bioc_packages <- c(
  "minfi", "sesame", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "FlowSorted.Blood.EPIC", "EpiDISH", "BiocParallel", "Biostrings",
  "GenomeInfoDb", "GenomicRanges", "SummarizedExperiment", "missMethyl", "sva"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("  Installing", pkg, "..."))

    # GitHub-only packages
    if (pkg %in% names(github_packages)) {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        install.packages("remotes")
      }
      remotes::install_github(github_packages[[pkg]], upgrade = "never")
    } else if (pkg %in% bioc_packages) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  } else {
    message(paste("  ✓", pkg, "already installed"))
  }
}

# ============================================================================
# 5. Install Required Packages (if --install flag is provided)
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

if ("--install" %in% args) {
  message("\nInstalling required packages...")

  # Install BiocManager first
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  # Install all required packages
  invisible(lapply(required_packages, install_if_missing))

  message("\n✓ Package installation complete!")

  # Snapshot the environment
  message("\nCreating renv snapshot...")
  renv::snapshot()

} else {
  message("\nTo install packages, run:")
  message("  Rscript scripts/00_setup_environment.R --install")
}

# ============================================================================
# 6. Verify Installation
# ============================================================================

message("\nVerifying package installation...")

verification_results <- data.frame(
  Package = character(),
  Installed = logical(),
  Version = character(),
  stringsAsFactors = FALSE
)

for (pkg in required_packages) {
  is_installed <- requireNamespace(pkg, quietly = TRUE)

  if (is_installed) {
    version <- as.character(packageVersion(pkg))
  } else {
    version <- "Not installed"
  }

  verification_results <- rbind(
    verification_results,
    data.frame(
      Package = pkg,
      Installed = is_installed,
      Version = version,
      stringsAsFactors = FALSE
    )
  )
}

print(verification_results)

# Summary
installed_count <- sum(verification_results$Installed)
total_count <- nrow(verification_results)

message(sprintf("\n%d/%d packages installed", installed_count, total_count))

if (installed_count == total_count) {
  message("✓ All required packages are installed!")
} else {
  warning("Some packages are missing. Run with --install flag.")
}

# ============================================================================
# 7. R 4.5 Compatibility Patch (if needed)
# ============================================================================

message("\nChecking R 4.5 compatibility...")

if (getRversion() >= "4.5.0") {
  message("  R 4.5+ detected - applying GenomicRanges patch...")

  patch_code <- '
# R 4.5 Patch for GenomicRanges::seqinfo
if (getRversion() >= "4.5.0") {
  library(GenomicRanges)
  setMethod("seqinfo", "RangedSummarizedExperiment",
            function(x) {
              GenomeInfoDb::seqinfo(SummarizedExperiment::rowRanges(x))
            })
  message("✓ GenomicRanges patch applied")
}
'

  # Write patch to a separate file
  writeLines(patch_code, "scripts/R45_patch.R")
  message("  ✓ Patch code written to scripts/R45_patch.R")
  message("  Include this at the start of analysis scripts")
}

# ============================================================================
# 8. Create Directory Structure
# ============================================================================

message("\nVerifying directory structure...")

required_dirs <- c(
  "data/raw",
  "data/processed",
  "figures",
  "scripts",
  "reports",
  "renv"
)

for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(paste("  Created:", dir))
  } else {
    message(paste("  ✓", dir, "exists"))
  }
}

# ============================================================================
# 9. Final Message
# ============================================================================

message("\n======================================")
message("Setup Complete!")
message("======================================")
message("\nNext steps:")
message("  1. Run: renv::restore()")
message("  2. Place reference IDAT files in bjarne/Reference_IDATS/")
message("  3. Run: Rscript scripts/01_build_reference_oneVsAll.R")
message("\nFor questions, contact: dpk194@student.bham.ac.uk")
