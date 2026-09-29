################################################################################
# Which purification the manuscript analyses run on - ONE switch
################################################################################
# Every downstream analysis (figures 2, 3, 4, 5, 8 and tables 5, 6) reads the
# purified beta matrix and the purified classifier calls. Both come in two
# variants, and this file decides which one they get, so the choice sits in one
# auditable place instead of being spread over six scripts.
#
#   "microglia4" (default, since 2026-08-31) - 4 components:
#       immune(6 blood types) + microglia + neuron + glia
#   "decode3" - 3 components: immune(6 blood types) + neuron + glia
#
# Why the default changed: the 3-component immune term came from Immune_total,
# which holds the six BLOOD-BORNE cell types only. Microglia sit in the brain
# compartment of that table, so they were in none of the three terms and were
# never actually removed - although they are the largest single immune
# population in MES tumours (14.5 percent versus 5.0 in RTK). The methods text
# claimed they were removed; with the 4-component variant that is now true.
#
# The reclassification does not weaken, it sharpens: MES 174/181 reclassified
# (96.1 percent) under 3 components, 179/181 (98.9 percent) under 4, with the
# median subclass score rising from 0.973 to 0.980 and low-confidence calls
# falling from 58 to 43. The two agree on 402 of 418 samples.
#
# To go back, set PURIFICATION to "decode3" and rerun scripts 19, 18, 29, 15b,
# 15e, 51. Nothing else needs touching.
################################################################################

PURIFICATION <- Sys.getenv("DECODE_PURIFICATION", "microglia4")

PURIFIED <- switch(PURIFICATION,
  microglia4 = list(
    beta  = "/mnt/volume-hel1-1/siewert/betamatrix_purified_microglia4.rds",
    calls = "data/processed/siewert_purified_microglia4_calls.csv"),
  decode3 = list(
    beta  = "/mnt/volume-hel1-1/siewert/betamatrix_purified_decode3.rds",
    calls = "data/processed/siewert_purified_decode3_calls.csv"),
  stop("Unknown PURIFICATION: ", PURIFICATION)
)
PURIFIED_BETA  <- PURIFIED$beta
PURIFIED_CALLS <- PURIFIED$calls
if (!file.exists(PURIFIED_CALLS)) stop("Calls file missing: ", PURIFIED_CALLS)
message(sprintf("[purification] %s\n  beta  : %s\n  calls : %s",
                PURIFICATION, PURIFIED_BETA, PURIFIED_CALLS))
