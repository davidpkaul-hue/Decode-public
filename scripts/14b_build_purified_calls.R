#!/usr/bin/env Rscript
################################################################################
# Build the canonical purified-classifier call table
################################################################################
# The purified MethylSets (script 13b) were classified externally with the
# Heidelberg Brain Tumor Classifier (mnp.v12epicv2, v12.8) at the Institute of
# Neuropathology, Charite - Universitaetsmedizin Berlin. This script turns the
# two returned score objects into the single table every downstream analysis
# reads (scripts 15b, 18, 19, 29, 35a).
#
# It previously had no script at all: the table was assembled ad hoc, so the
# step between "classifier was run" and "figures were made" was missing from the
# repository.
#
# Input:  <RESULTS>/classifier_scores_purified_v1.rds  (240 EPIC v1 samples)
#         <RESULTS>/classifier_scores_purified_v2.rds  (196 EPIC v2 samples)
#         each a list of three data frames (sub, family, class), every one with
#         Sample / predicted / maxscore plus the full score matrix
# Output: data/processed/siewert_purified_decode3_calls.csv
################################################################################

RESULTS <- "/mnt/volume-hel1-1/siewert/christin_decode3_results/christin-new_run"
OUT     <- "data/processed/siewert_purified_decode3_calls.csv"

read_batch <- function(file, batch) {
  x <- readRDS(file)
  stopifnot(identical(x$sub$Sample, x$family$Sample), identical(x$sub$Sample, x$class$Sample))
  data.frame(
    Sample    = x$sub$Sample,
    pred_sub  = x$sub$predicted,     score_sub = x$sub$maxscore,
    pred_fam  = x$family$predicted,  score_fam = x$family$maxscore,
    pred_cls  = x$class$predicted,   score_cls = x$class$maxscore,
    batch     = batch,
    stringsAsFactors = FALSE
  )
}

calls <- rbind(
  read_batch(file.path(RESULTS, "classifier_scores_purified_v1.rds"), "v1_EPIC"),
  read_batch(file.path(RESULTS, "classifier_scores_purified_v2.rds"), "v2_EPICv2")
)

# Readable subtype label. Everything the classifier does not call as one of the
# four GBM subclasses is collapsed to "other"; those are the low-confidence
# non-GBM calls discussed in the Results.
calls$subtype_pur <- c(
  GBM_RTK1     = "RTK I",
  GBM_RTK2     = "RTK II",
  GBM_MES_TYP  = "MES",
  GBM_MES_ATYP = "MES",
  GBM_PNC      = "PNC"
)[calls$pred_sub]
calls$subtype_pur[is.na(calls$subtype_pur)] <- "other"

message(sprintf("Samples: %d (v1 %d, v2 %d)", nrow(calls),
                sum(calls$batch == "v1_EPIC"), sum(calls$batch == "v2_EPICv2")))
print(table(calls$subtype_pur))

write.csv(calls, OUT, row.names = FALSE)
message("Saved: ", OUT)
