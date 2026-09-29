#!/usr/bin/env Rscript
################################################################################
# Build the public thesis repository from this working repository
################################################################################
# The working repo carries far more than the thesis: MOSAIC, HLA, contamination
# round 2, MGMT, the Berlin and Knoll cohorts. The public repo referenced under
# Code availability must contain only what reproduces the Masterarbeit.
#
# This is a WHITELIST, not a delete list: nothing is published unless it is named
# here. That is the safer direction for a cohort under GDPR, and it means the
# export can be regenerated and re-inspected at any time.
#
# Usage:  Rscript tools/export_public.R [target-directory]
# Default target: ../Decode-public-export
################################################################################

suppressPackageStartupMessages(library(tools))

TARGET <- commandArgs(trailingOnly = TRUE)
TARGET <- if (length(TARGET)) TARGET[1] else "../Decode-public-export"

# ---------------------------------------------------------------------------
# 1. Scripts that the thesis depends on
# ---------------------------------------------------------------------------
# Derived by tracing every manuscript figure backwards through the files each
# script reads and writes (2026-08-30). Grouped by the role they play.
SCRIPTS <- c(
  # environment and reference panel
  "00_setup_environment.R",
  "01_build_reference_oneVsAll.R",
  # internal validation and in-silico mixtures
  "04_validation_insilico.R",
  "05_validation_comparison.R",
  # cohort deconvolution
  "10_siewert_cohort_analysis.R",
  # ITHresolveGBM projection and differentiation score
  "11_ithresolve_comparison.R",
  "12_ithresolve_full_analysis.R",
  # in-silico purification
  "13_purify_full_matrix.R",           # superseded ITH-fraction version, kept because 13c compares against it
  "13b_purify_decode_compare.R",       # 3-component version, superseded but kept: it is the sensitivity comparison
  "13c_compare_purified_ith_decode.R", # sensitivity: DECODE vs ITH fractions
  "53_purify_microglia4.R",            # THE purification behind every reported result (4 components)
  "14b_build_purified_calls.R",        # canonical call table, 3-component run
  "14c_build_purified_calls_microglia4.R", # canonical call table, 4-component run
  "56_microglia4_reclassification.R",  # 4- against 3-component reclassification
  "58_manuscript_numbers.R",           # derives every purification number quoted in the Results
  "_purification_config.R",            # which purification the analyses run on - one switch
  # figures and downstream analyses
  "15b_tsne_combat.R",
  "19_manuscript_figures.R",
  "19b_fig1_validation.R",
  "29_fig7_differentiation_transition.R",
  # in-silico contamination
  "31_insilico_contamination.R",
  "32_contamination_diffscore_and_export.R",
  "33_contamination_classifier_analysis.R",
  # normal-brain validation of the glial fractions
  "48_normalbrain_decode.R",
  "49_normalbrain_replication.R",
  # sensitivity check of the ITHresolveGBM profile file (Methods)
  "60_ithresolve_featureselected_check.R"
)

# ---------------------------------------------------------------------------
# 2. Result tables, with an explicit column whitelist
# ---------------------------------------------------------------------------
# Sample identifiers are replaced by running numbers; the mapping is written
# next to the working repo and is NOT part of the export. Clinical variables
# (age, sex, free-text diagnosis, case and NGS identifiers) are never exported:
# Table 1 of the manuscript therefore cannot be rebuilt from the public data.
CELLTYPES <- c("Bcell","CD4T","CD8T","Mono","NK","Neu","Neurons",
               "Oligodendrocytes","Astrocytes","Microglia","Endothelial",
               "Fibroblasts","GBM_malignant")

TABLES <- list(
  list(file = "siewert_cohort_results.csv", id = "Idat",
       keep = c(CELLTYPES, "Immune_total","Brain_total","Stromal_total",
                "Max_family","Max_family_score","Max_class","Max_class_score",
                "Max_subclass","Max_subclass_score","RF_purity","CNV_Load",
                "focal.CN.Amplifications","focal.CN.Deletions")),
  # Input to 14c/56, not produced by any script here: the Heidelberg classifier
  # runs externally at the Charite. Without it the purification chain can be read
  # but not re-run, so it is exported (pseudonymised) rather than only its result.
  list(file = "classifier_scores_microglia4.csv", id = "Sample", keep = NULL),
  list(file = "siewert_purified_decode3_calls.csv", id = "Sample", keep = NULL),
  list(file = "siewert_purified_microglia4_calls.csv", id = "Sample", keep = NULL),
  list(file = "microglia4_reclassification.csv", id = "Idat", keep = NULL),
  list(file = "manuscript_numbers_table.csv", id = "Idat", keep = NULL),
  list(file = "siewert_ithresolve_full.csv", id = "Idat", keep = NULL),
  list(file = "purified_na_per_sample.csv", id = "Idat", keep = NULL),
  list(file = "contamination_manifest.csv", id = c("pseudo_sample","orig_sample"), keep = NULL),
  list(file = "contamination_classifier_calls.csv", id = "Sample", keep = NULL),
  list(file = "contamination_diffscore.csv", id = "sample", keep = NULL),
  # no sample identifiers in these
  list(file = "classifier_score_paired_tests.csv", id = NULL, keep = NULL),
  list(file = "purification_3comp_vs_4comp.csv", id = NULL, keep = NULL),
  list(file = "purified_decode_vs_ith.csv", id = NULL, keep = NULL),
  # public data, GSM accessions may stay as they are
  list(file = "normalbrain_decode_fractions.csv", id = NULL, keep = NULL),
  list(file = "normalbrain_replication_fractions.csv", id = NULL, keep = NULL),
  list(file = "ithresolve_featureselected_diffscore.csv", id = "Idat", keep = NULL)
)

# The working repo's .gitignore is deliberately NOT copied: it carries
# "data/processed/**", which would ignore exactly the derived result tables this
# export exists to publish. Older tables stayed tracked while new ones silently
# did not, so the public repo is given its own rules instead (written below).
OTHER <- c("LICENSE.md", "renv.lock", "DECODE.Rproj")  # README: siehe unten
# masterarbeit.md is deliberately NOT exported (decision 29.09.2026). The Data
# availability statement promises the code and the processed results behind
# every figure and number, not the text of the dissertation. Birmingham does
# not archive the dissertation either, and publishing the full text here would
# put it online ahead of the planned paper. The figures and the build scripts
# stay, so the manuscript remains reproducible from the repository.
MANUSCRIPT <- c("build_docx.sh", "color_span.lua",
                "caption_fontsize.py", "reference_black.docx")

SENTRIX <- "20[0-9]{10}_R[0-9]{2}C[0-9]{2}"

# ---------------------------------------------------------------------------
# 3. Build the identifier map from every sentrix ID the tables contain
# ---------------------------------------------------------------------------
message("Collecting sample identifiers ...")
ids <- character()
for (tb in TABLES) {
  p <- file.path("data/processed", tb$file)
  if (!file.exists(p)) next
  d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in tb$id) if (!is.null(col) && col %in% names(d))
    ids <- c(ids, grep(SENTRIX, d[[col]], value = TRUE))
}
ids <- sort(unique(sub("__.*$", "", ids)))
idmap <- setNames(sprintf("S%03d", seq_along(ids)), ids)
message(sprintf("  %d distinct samples -> S001 to S%03d", length(ids), length(ids)))

deident <- function(x) {
  # handles both bare IDs and the pseudo-sample form "<id>__f30"
  vapply(as.character(x), function(v) {
    m <- regmatches(v, regexpr(SENTRIX, v))
    if (!length(m)) return(v)
    sub(m, idmap[[m]], v, fixed = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

# ---------------------------------------------------------------------------
# 4. Assemble the export
# ---------------------------------------------------------------------------
unlink(TARGET, recursive = TRUE)
for (d in c("scripts", "data/processed", "manuscript/figures", "tools"))
  dir.create(file.path(TARGET, d), recursive = TRUE, showWarnings = FALSE)

missing <- SCRIPTS[!file.exists(file.path("scripts", SCRIPTS))]
if (length(missing)) stop("Missing scripts: ", paste(missing, collapse = ", "))
invisible(file.copy(file.path("scripts", SCRIPTS), file.path(TARGET, "scripts"), overwrite = TRUE))
message(sprintf("Scripts: %d of %d in the working repo", length(SCRIPTS),
                length(list.files("scripts", pattern = "\\.R$"))))

for (tb in TABLES) {
  p <- file.path("data/processed", tb$file)
  if (!file.exists(p)) { message("  missing, skipped: ", tb$file); next }
  d <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in tb$id) if (!is.null(col) && col %in% names(d)) d[[col]] <- deident(d[[col]])
  if (!is.null(tb$keep)) {
    keep <- c(tb$id, intersect(tb$keep, names(d)))
    d <- d[, keep, drop = FALSE]
  }
  write.csv(d, file.path(TARGET, "data/processed", tb$file), row.names = FALSE)
}
message(sprintf("Tables: %d", length(TABLES)))

md <- readLines("manuscript/masterarbeit.md", warn = FALSE)
figs <- unique(sub(".*figures/", "", regmatches(md, regexpr("figures/[^)]+", md))))
figs <- figs[nzchar(figs)]
invisible(file.copy(file.path("manuscript/figures", figs),
                    file.path(TARGET, "manuscript/figures"), overwrite = TRUE))
invisible(file.copy(file.path("manuscript", MANUSCRIPT), file.path(TARGET, "manuscript"), overwrite = TRUE))
invisible(file.copy(OTHER[file.exists(OTHER)], TARGET, overwrite = TRUE))

# .gitignore written for the public repository, not inherited from the working one
writeLines(c(
  "# Rules for the public repository. The derived result tables under",
  "# data/processed are part of the publication and are deliberately tracked.",
  "",
  "# Never any primary data",
  "*.idat",
  "*.rds",
  "*.RData",
  "*.Rda",
  "data/raw/",
  "data/external/",
  "",
  "# Local identifier map, must never be published",
  ".public_id_map.csv",
  "",
  "# R working files",
  ".Rproj.user/",
  ".Rhistory",
  ".RData",
  ".Ruserdata",
  "renv/library/",
  "",
  "# Build products",
  "manuscript/*.docx",
  "!manuscript/reference_black.docx",   # Build-Vorlage, kein Build-Produkt
  "Rplots.pdf",
  "",
  "# System",
  ".DS_Store",
  "Thumbs.db"
), file.path(TARGET, ".gitignore"))
invisible(file.copy("tools/export_public.R", file.path(TARGET, "tools"), overwrite = TRUE))
message(sprintf("Manuscript: %d figures + %d files", length(figs), length(MANUSCRIPT)))

# ---------------------------------------------------------------------------
# 5. Refuse to finish if anything identifiable slipped through
# ---------------------------------------------------------------------------
message("\nScanning the export for sample identifiers ...")
files <- list.files(TARGET, recursive = TRUE, full.names = TRUE)
text <- files[grepl("\\.(R|md|csv|sh|py|lua|Rproj|lock|json)$", files)]
hits <- character()
for (f in text) {
  ln <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
  if (any(grepl(SENTRIX, ln))) hits <- c(hits, f)
}
if (length(hits)) {
  stop("Sample identifiers found in:\n  ", paste(hits, collapse = "\n  "),
       "\nExport aborted; fix the source or extend the whitelist.")
}
message("  clean: no sentrix identifier in any exported text file")

# ---------------------------------------------------------------------------
# 6. Describe the scope of the export in the published README
# ---------------------------------------------------------------------------
# The public README is written for readers of the repository, not for us, so it
# is maintained as its own file (manuscript/README_public.md) and copied over
# the working README here. It states the scope, the two load-bearing
# limitations, and what the data folder deliberately leaves out.
readme_src <- "manuscript/README_public.md"
stopifnot(file.exists(readme_src))
writeLines(c(readLines(readme_src, warn = FALSE), "",
             sprintf("Generated by `tools/export_public.R` on %s.", format(Sys.Date()))),
           file.path(TARGET, "README.md"))
message("README: taken from ", readme_src)

writeLines(paste0(names(idmap), ",", idmap),
           "data/processed/.public_id_map.csv")
message("Identifier map kept locally at data/processed/.public_id_map.csv (never exported)")
# ---------------------------------------------------------------------------
# Regression guard: did anything that was exported before stop being exported?
# ---------------------------------------------------------------------------
# The export is a whitelist, so a file drops out silently the moment its entry
# here is lost. That happened on 2026-09-06: the entry for script 60 was wiped
# by an accidental `git reset --hard`, and the next export would have removed a
# script the Methods cite, without saying so. The manifest is tracked in git, so
# an unintended removal also surfaces as a diff before anyone commits it.
MANIFEST <- "tools/export_manifest.txt"
now <- sort(list.files(TARGET, recursive = TRUE, all.files = TRUE, no.. = TRUE))

if (file.exists(MANIFEST)) {
  before <- readLines(MANIFEST, warn = FALSE)
  gone  <- setdiff(before, now)
  added <- setdiff(now, before)
  if (length(gone)) {
    message("\n", strrep("!", 72))
    message(sprintf("WARNING: %d file(s) exported previously are NOT in this export:", length(gone)))
    for (g in gone) message("  - ", g)
    message("A whitelist entry above is missing, or the source file is gone.")
    message("If the removal is intended, commit the updated manifest.")
    message(strrep("!", 72))
  }
  if (length(added)) {
    message(sprintf("\nNew in this export: %d file(s)", length(added)))
    for (a in utils::head(added, 10)) message("  + ", a)
    if (length(added) > 10) message(sprintf("  ... and %d more", length(added) - 10))
  }
  if (!length(gone) && !length(added)) message("\nSame set of files as the last run (contents may still have changed).")
} else {
  message("\nNo previous manifest; writing the first one.")
}
writeLines(now, MANIFEST)
message(sprintf("Manifest: %s (%d files)", MANIFEST, length(now)))

message(sprintf("\nExport written to %s", normalizePath(TARGET)))
