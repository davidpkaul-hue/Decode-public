# ============================================================================
# Script 56: Reclassification after the 4-component purification (microglia)
# ----------------------------------------------------------------------------
# The result-carrying purification regressed out immune + neuron + glia, where
# the immune term came from Immune_total = the six BLOOD-BORNE cell types.
# Microglia sit in the brain compartment of that table, so they were in none of
# the three terms and were never actually removed - even though microglia are
# the largest single immune population in MES tumours. Script 53 added them as
# their own regressor (4 components) and Christin ran the real v12.8 on the
# result.
#
# This script reads her calls and asks the question the thesis rests on: does
# the MES reclassification hold when microglia are genuinely removed?
#
# The 3-component numbers it is compared against (from the DECODE re-run):
#   MES 181 originally -> 174 reclassified (96.1%), 165 to RTK II, 2 to RTK I,
#   7 stay MES; RTK II stable 149/157 (94.9%); RTK I 60/76 (78.9%).
#
# Column names in Christin's CSV are auto-detected, since the format has varied
# between runs (RDS list in June, CSV here).
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2) })
if (file.exists("scripts/R45_patch.R")) source("scripts/R45_patch.R")

CALLS <- "data/processed/classifier_scores_microglia4.csv"
if (!file.exists(CALLS)) stop("Not found: ", CALLS,
  "\nPut Christin's classifier_scores_microglia4.csv there first.")
raw <- read.csv(CALLS, stringsAsFactors = FALSE, check.names = FALSE)
message(sprintf("Read %d rows, %d columns: %s", nrow(raw), ncol(raw),
                paste(names(raw), collapse = ", ")))

pick <- function(pats, what) {
  for (p in pats) { h <- grep(p, names(raw), ignore.case = TRUE, value = TRUE)
                    if (length(h)) return(h[1]) }
  stop("Could not find the ", what, " column among: ", paste(names(raw), collapse = ", "))
}
c_id  <- pick(c("^sample$", "sample", "^id$", "idat", "sentrix"), "sample id")
c_sub <- pick(c("pred_sub", "subclass$", "^subclass", "predicted.*sub", "mcf.*sub"), "subclass call")
c_scr <- pick(c("score_sub", "subclass.*score", "score.*subclass", "^score"), "subclass score")
message(sprintf("Using: id = %s | subclass = %s | score = %s", c_id, c_sub, c_scr))

lab <- function(x) ifelse(grepl("RTK1|RTK I\\b", x), "RTK I",
                  ifelse(grepl("RTK2|RTK II",   x), "RTK II",
                  ifelse(grepl("mesenchymal|MES", x, ignore.case = TRUE), "MES",
                  ifelse(grepl("primitive neuronal|PNC", x, ignore.case = TRUE), "PNC", "other"))))
new <- data.frame(Sample = as.character(raw[[c_id]]),
                  pred_sub = raw[[c_sub]],
                  score_sub = suppressWarnings(as.numeric(raw[[c_scr]])),
                  subtype_mg4 = lab(raw[[c_sub]]), stringsAsFactors = FALSE)

# --- original (unpurified) labels and the 3-component purified ones ----------
meta <- read.csv("data/processed/siewert_cohort_results.csv", check.names = FALSE,
                 stringsAsFactors = FALSE)
meta <- meta[meta$is_GBM %in% c(TRUE, "TRUE"), ]
sc <- meta$Max_subclass
meta$subtype_orig <- lab(sc)
old <- read.csv("data/processed/siewert_purified_decode3_calls.csv", stringsAsFactors = FALSE)

d <- merge(meta[, c("Idat", "subtype_orig")], new, by.x = "Idat", by.y = "Sample")
d <- merge(d, old[, c("Sample", "subtype_pur")], by.x = "Idat", by.y = "Sample", all.x = TRUE)
names(d)[names(d) == "subtype_pur"] <- "subtype_3comp"
message(sprintf("\n%d GBM samples matched to Christin's 4-component calls", nrow(d)))

cat("\n=== Original -> 4-component purified ===\n")
print(table(original = d$subtype_orig, purified_4comp = d$subtype_mg4))

cat("\n=== Side by side with the 3-component purification ===\n")
for (s in c("MES", "RTK I", "RTK II")) {
  a <- d[d$subtype_orig == s, ]
  if (!nrow(a)) next
  ch3 <- sum(a$subtype_3comp != s, na.rm = TRUE); ch4 <- sum(a$subtype_mg4 != s)
  cat(sprintf("  %-7s n=%3d | 3-comp reclassified %3d (%.1f%%) | 4-comp %3d (%.1f%%)\n",
      s, nrow(a), ch3, 100*ch3/nrow(a), ch4, 100*ch4/nrow(a)))
}

mes <- d[d$subtype_orig == "MES", ]
cat(sprintf("\n=== Where the %d originally-MES tumours land ===\n", nrow(mes)))
print(rbind(`3-comp` = table(factor(mes$subtype_3comp, levels = c("RTK I","RTK II","MES","PNC","other"))),
            `4-comp` = table(factor(mes$subtype_mg4,   levels = c("RTK I","RTK II","MES","PNC","other")))))

cat("\n=== Do the two purifications agree per sample? ===\n")
ok <- !is.na(d$subtype_3comp)
cat(sprintf("  identical call: %d / %d (%.1f%%)\n", sum(d$subtype_mg4[ok] == d$subtype_3comp[ok]),
            sum(ok), 100*mean(d$subtype_mg4[ok] == d$subtype_3comp[ok])))
disc <- d[ok & d$subtype_mg4 != d$subtype_3comp, c("Idat","subtype_orig","subtype_3comp","subtype_mg4","score_sub")]
if (nrow(disc)) { cat(sprintf("  %d discordant:\n", nrow(disc)))
  print(head(disc[order(-disc$score_sub), ], 20), row.names = FALSE) }

cat("\n=== Confidence ===\n")
cat(sprintf("  median subclass score: 3-comp %.3f | 4-comp %.3f\n",
            median(old$score_sub, na.rm = TRUE), median(d$score_sub, na.rm = TRUE)))
cat(sprintf("  calls below 0.5: 3-comp %d | 4-comp %d\n",
            sum(old$score_sub < 0.5, na.rm = TRUE), sum(d$score_sub < 0.5, na.rm = TRUE)))

write.csv(d, "data/processed/microglia4_reclassification.csv", row.names = FALSE)
message("\nWrote data/processed/microglia4_reclassification.csv")
