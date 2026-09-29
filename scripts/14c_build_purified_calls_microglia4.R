#!/usr/bin/env Rscript
################################################################################
# Build the canonical purified-classifier call table, 4-component purification
################################################################################
# Same job as script 14b, for the purification that removes microglia as their
# own regressor (script 53). Christin returned one CSV this time rather than the
# two RDS score objects of the June run, so the reader is simpler.
#
# Why this became the main analysis: the 3-component purification regressed out
# immune + neuron + glia, where the immune term was Immune_total = the six
# blood-borne types. Microglia sit in the brain compartment of that table and
# were therefore in none of the three terms, so they were never removed -
# although they are the largest single immune population in MES tumours.
#
# Input:  data/processed/classifier_scores_microglia4.csv (436 samples, v12.8)
# Output: data/processed/siewert_purified_microglia4_calls.csv
#         same columns as siewert_purified_decode3_calls.csv, so every
#         downstream script reads it unchanged.
################################################################################

IN  <- "data/processed/classifier_scores_microglia4.csv"
OUT <- "data/processed/siewert_purified_microglia4_calls.csv"

raw <- read.csv(IN, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Sample","Array","Subclass","Score_sub","Family","Score_fam",
                "Class","Score_class") %in% names(raw)))

calls <- data.frame(
  Sample    = raw$Sample,
  pred_sub  = raw$Subclass, score_sub = raw$Score_sub,
  pred_fam  = raw$Family,   score_fam = raw$Score_fam,
  pred_cls  = raw$Class,    score_cls = raw$Score_class,
  batch     = ifelse(raw$Array == "EPICv1", "v1_EPIC", "v2_EPICv2"),
  stringsAsFactors = FALSE
)

# Identical label mapping to script 14b: anything the classifier does not call
# as one of the four GBM subclasses is collapsed to "other".
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

# Guard: the 3-component table is the reference for sample identity, not for
# calls. If the sample sets ever diverge, every downstream join would silently
# drop rows, so fail loudly here instead.
ref <- read.csv("data/processed/siewert_purified_decode3_calls.csv", stringsAsFactors = FALSE)
miss <- setdiff(ref$Sample, calls$Sample)
if (length(miss)) stop(length(miss), " samples missing versus the 3-component table")
message("Sample set identical to the 3-component table.")

write.csv(calls, OUT, row.names = FALSE)
message("Saved: ", OUT)
