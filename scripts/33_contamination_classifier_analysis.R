# ============================================================================
# Script 33: Contamination experiment - REAL Heidelberg classifier read-out
# Produces manuscript FIGURE 6 (file manuscript/figures/Fig7_contamination.png -- the
# file name predates the current figure numbering).
# ----------------------------------------------------------------------------
# Christin ran mnp.v12epicv2 (v12.8) on the 72 in-silico contaminated pseudo-
# samples (12 high-purity RTK x 6 immune fractions f=0..0.5). This script joins
# her calls to the manifest, finds the immune fraction at which each RTK flips
# to MES, compares it to the real MES immune range, and renders the figure.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

calls <- read.csv("data/processed/contamination_classifier_calls.csv",
                  stringsAsFactors = FALSE)
mani  <- read.csv("data/processed/contamination_manifest.csv",
                  stringsAsFactors = FALSE)

# --- join calls to manifest (pseudo_sample == Sample) -----------------------
d <- merge(calls, mani, by.x = "Sample", by.y = "pseudo_sample", all.x = TRUE)
# NOTE: only GBM_MES_TYP counts as MES here, whereas scripts 14 and 19 count
# GBM_MES_TYP and GBM_MES_ATYP. It makes no difference to any number reported:
# none of the 72 pseudo-samples was ever called GBM_MES_ATYP. Kept as-is so the
# published figure stays reproducible, but harmonise before reusing this code.
d$is_MES <- d$Subclass == "GBM_MES_TYP"
# short subclass family for plotting
d$call <- ifelse(grepl("MES", d$Subclass), "MES",
          ifelse(grepl("RTK1", d$Subclass), "RTK I",
          ifelse(grepl("RTK2", d$Subclass), "RTK II", "other")))
d <- d[order(d$orig_sample, d$frac), ]

# --- per-sample flip fraction (first f called MES) --------------------------
samps <- unique(d$orig_sample)
flip <- do.call(rbind, lapply(samps, function(s) {
  ds <- d[d$orig_sample == s, ]
  mes_f <- ds$frac[ds$is_MES]
  data.frame(orig_sample = s,
             orig_subtype = ds$orig_subtype[1],
             call_f0 = ds$call[ds$frac == 0],
             flip_frac = if (length(mes_f)) min(mes_f) else NA_real_,
             flip_decode_immune = if (length(mes_f))
                 ds$expected_decode_immune[ds$frac == min(mes_f)] else NA_real_)
}))

cat("\n=== Per-sample flip to MES (real Heidelberg v12.8) ===\n")
print(flip, row.names = FALSE)

n_flip <- sum(!is.na(flip$flip_frac))
cat(sprintf("\n%d of %d RTK samples flip to MES within f<=0.5\n", n_flip, nrow(flip)))
cat(sprintf("Flip immune fraction (added f):   median %.0f%%  (range %.0f-%.0f%%)\n",
    100*median(flip$flip_frac, na.rm=TRUE),
    100*min(flip$flip_frac, na.rm=TRUE), 100*max(flip$flip_frac, na.rm=TRUE)))
cat(sprintf("Flip DECODE immune fraction:      median %.1f%%  (range %.1f-%.1f%%)\n",
    100*median(flip$flip_decode_immune, na.rm=TRUE),
    100*min(flip$flip_decode_immune, na.rm=TRUE),
    100*max(flip$flip_decode_immune, na.rm=TRUE)))
cat("Real MES immune reference (cohort): ~33.3% mean, ~24-42% IQR (DECODE)\n")

write.csv(flip, "data/processed/contamination_flip_summary.csv", row.names = FALSE)
write.csv(d[, c("Sample","orig_sample","orig_subtype","frac","expected_decode_immune",
                "Subclass","Score_sub","call","is_MES")],
          "data/processed/contamination_classifier_joined.csv", row.names = FALSE)

# --- Figure: MES subclass score vs added immune fraction, per sample --------
d$mes_score <- ifelse(d$call == "MES", d$Score_sub, NA)
# use Score_sub as the MES-ness proxy only when MES; else plot 0 baseline via call colour
p <- ggplot(d, aes(100*frac, Score_sub, group = orig_sample, colour = orig_subtype)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0) +
  geom_line(alpha = 0.5) +
  geom_point(aes(shape = call), size = 2) +
  scale_shape_manual(values = c("MES"=17, "RTK I"=16, "RTK II"=16, "other"=4)) +
  labs(x = "Added immune fraction f (%)",
       y = "Heidelberg subclass score (confidence of called subclass)",
       colour = "Original subtype", shape = "v12.8 call",
       title = "Real Heidelberg classifier flips RTK to MES under in-silico contamination",
       subtitle = sprintf("%d/%d high-purity RTK tumors become MES by adding immune signal; triangles = MES call", n_flip, nrow(flip))) +
  theme_bw()
ggsave("figures/contamination_classifier_flip.png", p, width = 9, height = 5.5, dpi = 150)
cat("\nWrote figures/contamination_classifier_flip.png\n")

# --- Figure 2: fraction of samples called MES vs f (dose-response) ----------
agg <- aggregate(is_MES ~ frac, data = d, FUN = function(x) 100*mean(x))
names(agg)[2] <- "pct_MES"
cat("\n=== % of RTK samples called MES at each fraction ===\n"); print(agg, row.names=FALSE)
p2 <- ggplot(agg, aes(100*frac, pct_MES)) +
  geom_line(colour = "grey30") + geom_point(size = 3, colour = "firebrick") +
  scale_y_continuous(limits = c(0,100)) +
  labs(x = "Added immune fraction f (%)", y = "RTK tumors called MES (%)",
       title = "Dose-response: immune admixture drives the MES call",
       subtitle = "12 high-purity RTK tumors, real Heidelberg v12.8") +
  theme_bw()
ggsave("figures/contamination_doseresponse.png", p2, width = 7, height = 4.5, dpi = 150)
cat("Wrote figures/contamination_doseresponse.png\n")

# ---------------------------------------------------------------------------
# Manuscript Figure 8 (two panels via gridExtra) -> manuscript/figures/
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(gridExtra))
# The manuscript reports the ten tumors whose f = 0 control reproduced the original
# subtype (Methods). The two that failed it (least confident original call; 12.1%
# failed probes) are excluded here; both also reached a MES call, so the headline
# numbers above, computed on all twelve, do not depend on the exclusion.
ctrl0 <- d[d$frac == 0, c("orig_sample", "orig_subtype", "call")]
unclean <- ctrl0$orig_sample[as.character(ctrl0$call) != as.character(ctrl0$orig_subtype)]
stopifnot(length(unclean) == 2)
d <- d[!d$orig_sample %in% unclean, ]
stopifnot(length(unique(d$orig_sample)) == 10)
agg <- aggregate(is_MES ~ frac, data = d, FUN = function(x) 100*mean(x))
names(agg)[2] <- "pct_MES"
cat("\n=== Manuscript Fig 6, ten clean-control tumors: % MES per fraction ===\n")
print(agg, row.names = FALSE)
# Panel a: dose-response (% MES vs added immune fraction)
pa <- ggplot(agg, aes(100*frac, pct_MES)) +
  geom_line(colour = "grey40") +
  geom_point(size = 3, colour = "firebrick") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  labs(x = "Added immune fraction (%)", y = "RTK tumors called MES (%)",
       title = "a") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12))

# Panel b: per-sample called-subclass score trajectory, MES calls as triangles
d$orig2 <- factor(d$orig_subtype, levels = c("RTK I", "RTK II"))
pb <- ggplot(d, aes(100*frac, Score_sub, group = orig_sample, colour = orig2)) +
  annotate("rect", xmin = 24, xmax = 42, ymin = -Inf, ymax = Inf,
           alpha = 0.10, fill = "grey50") +
  geom_line(alpha = 0.45) +
  # Triangles (MES calls) drawn larger than circles so the flip stands out.
  geom_point(aes(shape = call, size = call == "MES")) +
  scale_size_manual(values = c("FALSE" = 1.8, "TRUE" = 3.2), guide = "none") +
  # RTK I and RTK II share the circle -- they were two identical legend entries,
  # the subtype is carried by colour anyway. The cross (non-GBM call) was missing
  # from the legend entirely and only explained in the caption.
  scale_shape_manual(values = c("MES"=17, "RTK I"=16, "RTK II"=16, "other"=4),
                     breaks = c("RTK II", "MES", "other"),
                     labels = c("RTK", "MES", "non-GBM")) +
  # Manuskript-Subtyp-Palette (wie Fig 2, 3, 7a); gruen ist dort MES und darf
  # hier nicht fuer RTK I stehen.
  scale_colour_manual(values = c("RTK I" = "#E41A1C", "RTK II" = "#377EB8")) +
  guides(colour = guide_legend(order = 1),
         shape = guide_legend(order = 2, override.aes = list(size = c(1.8, 3.2)))) +
  scale_x_continuous(breaks = seq(0, 50, 10)) +
  labs(x = "Added immune fraction (%)",
       y = "Heidelberg subclass score", colour = "Unmixed\ntumor", shape = "Mixture\ncalled",
       title = "b") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 12))

g <- arrangeGrob(pa, pb, ncol = 2, widths = c(1, 1.42))
ggsave("manuscript/figures/Fig7_contamination.png", g, width = 11, height = 4.6, dpi = 150)
cat("Wrote manuscript/figures/Fig7_contamination.png\n")
