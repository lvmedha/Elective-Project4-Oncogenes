# 11_myc_cnv_analysis.R
# MYC-focused follow-up to R/09: why does AlphaMissense call some MYC
# missense variants "likely benign" even though the cell line carrying them
# is highly MYC-dependent (very negative Chronos)?
#
# Hypothesis: those lines are MYC-amplified (or otherwise overexpress MYC).
# Dependency is set by gene dosage (CN, transcription, translocation), not
# by the residue change. AlphaMissense scores the residue, so it correctly
# calls many variants "likely benign" structurally; the missense event is
# a passenger on top of an already-amplified locus.
#
# Inputs:
#   ~/Documents/depmap_cnv.csv                 — gene-level relative CN
#   ~/Documents/depmap_CGE.csv                 — gene-level Chronos
#   ~/Documents/depmap_meta.csv                — lineage / disease
#   results/ped_gof_snv/09_am_chronos_per_event.tsv
#
# Outputs (results/<RUN_NAME>/):
#   11_myc_cnv_per_line.tsv         one row per cohort MYC-mutant cell line
#   11_myc_cnv_summary.txt          counts, full-DepMap correlation, puzzling cases
#   11_myc_cnv_vs_chronos.pdf       scatter: full DepMap + cohort overlay
#   11_myc_cnv_per_line_bars.pdf    per-line CN bar, colored by AM class

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
RUN_NAME <- "ped_gof_snv"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CN_CSV     <- file.path(DOC_DIR, "depmap_cnv.csv")
CGE_CSV    <- file.path(DOC_DIR, "depmap_CGE.csv")
META_CSV   <- file.path(DOC_DIR, "depmap_meta.csv")
EVENTS_TSV <- file.path(OUT_DIR, "09_am_chronos_per_event.tsv")

for (f in c(CN_CSV, CGE_CSV, META_CSV, EVENTS_TSV)) {
  if (!file.exists(f)) stop("Missing: ", f)
}

# ---- 1. Load only the MYC column from the 400 MB CN matrix ---------------
cn_header <- fread(CN_CSV, nrows = 0)
myc_cn_col <- grep("^MYC \\(", colnames(cn_header), value = TRUE)
stopifnot(length(myc_cn_col) == 1L)
cn_keys <- intersect(c("ModelID"), colnames(cn_header))
stopifnot(length(cn_keys) == 1L)
cn <- fread(CN_CSV, select = c(cn_keys, myc_cn_col))
setnames(cn, c(cn_keys, myc_cn_col), c("ModelID", "MYC_CN"))
cat(sprintf("MYC CN loaded for %d DepMap lines (column='%s')\n",
            nrow(cn), myc_cn_col))

# ---- 2. Load only the MYC column from the CGE matrix ---------------------
cge_header <- fread(CGE_CSV, nrows = 0)
myc_cge_col <- grep("^MYC \\(", colnames(cge_header), value = TRUE)
stopifnot(length(myc_cge_col) == 1L)
cge_myc <- fread(CGE_CSV,
                 select = c(colnames(cge_header)[1], myc_cge_col))
setnames(cge_myc, c("ModelID", "MYC_Chronos"))

# ---- 3. Lineage metadata --------------------------------------------------
meta <- fread(META_CSV,
              select = c("ModelID", "OncotreeLineage",
                         "OncotreePrimaryDisease", "StrippedCellLineName"))

# ---- 4. Cohort MYC events from R/09 ---------------------------------------
ev <- fread(EVENTS_TSV)
myc_ev <- ev[gene_symbol == "MYC"]
cat(sprintf("Cohort MYC events: %d (across %d cell lines)\n",
            nrow(myc_ev), uniqueN(myc_ev$ACH)))

myc_ev <- merge(myc_ev, cn, by.x = "ACH", by.y = "ModelID", all.x = TRUE)
myc_ev <- merge(myc_ev, meta, by.x = "ACH", by.y = "ModelID", all.x = TRUE)

# Per-line summary: prefer most-pathogenic AM class if a line has multiple
# variants. RAJI is the only line with a duplicated locus (T73I as SNV +
# MNV), and only one of those calls has AM, so coalescing is the right move.
am_priority <- function(classes) {
  if (any(grepl("likely pathogenic", classes, fixed = TRUE))) {
    "likely pathogenic"
  } else if (any(grepl("ambiguous", classes, fixed = TRUE))) {
    "ambiguous"
  } else if (any(grepl("likely benign", classes, fixed = TRUE))) {
    "likely benign"
  } else {
    "no AM score"
  }
}

per_line <- myc_ev[, .(
  variants = paste(sort(unique(HGVSp_Short)), collapse = ", "),
  am_classes_raw = paste(sort(unique(am_class_lab)), collapse = ", "),
  am_class_top = am_priority(unique(am_class_lab)),
  am_max_pathogenicity = suppressWarnings({
    m <- max(am_pathogenicity, na.rm = TRUE)
    if (is.finite(m)) m else NA_real_
  }),
  Chronos = unique(GeneEffect)[1],
  MYC_CN = unique(MYC_CN)[1],
  lineage = unique(OncotreeLineage)[1],
  disease = unique(OncotreePrimaryDisease)[1],
  CellLine = unique(StrippedCellLineName)[1]
), by = ACH]

per_line[, am_class_top := factor(
  am_class_top,
  levels = c("likely pathogenic", "ambiguous", "likely benign", "no AM score")
)]

fwrite(per_line, file.path(OUT_DIR, "11_myc_cnv_per_line.tsv"), sep = "\t")

# ---- 5. Full DepMap background ------------------------------------------
bg <- merge(cge_myc, cn, by = "ModelID", all = FALSE)
bg <- bg[!is.na(MYC_Chronos) & !is.na(MYC_CN)]
bg[, is_cohort := ModelID %in% per_line$ACH]
cat(sprintf("Full DepMap background: %d lines with both CN and Chronos\n",
            nrow(bg)))

# ---- 6. Summary text ------------------------------------------------------
sink(file.path(OUT_DIR, "11_myc_cnv_summary.txt"))
cat("11_myc_cnv_analysis.R summary\n")
cat("=============================\n\n")
cat(sprintf("Full DepMap (lines with MYC CN and MYC Chronos): n=%d\n",
            nrow(bg)))
cat(sprintf("  median MYC_CN     : %.2f  (IQR %.2f to %.2f, range %.2f to %.2f)\n",
            median(bg$MYC_CN), quantile(bg$MYC_CN, 0.25),
            quantile(bg$MYC_CN, 0.75),
            min(bg$MYC_CN), max(bg$MYC_CN)))
cat(sprintf("  median MYC_Chronos: %.2f  (IQR %.2f to %.2f)\n",
            median(bg$MYC_Chronos), quantile(bg$MYC_Chronos, 0.25),
            quantile(bg$MYC_Chronos, 0.75)))
ct <- suppressWarnings(cor.test(bg$MYC_CN, bg$MYC_Chronos,
                                method = "spearman", exact = FALSE))
cat(sprintf("  Spearman MYC_CN vs MYC_Chronos: rho=%.3f, p=%.2g\n",
            unname(ct$estimate), ct$p.value))
cat("  (negative rho means: more amplification -> more negative Chronos = ",
    "more dependent)\n\n", sep = "")

cat("Cohort MYC-mutant cell lines (sorted by Chronos, most dependent first):\n\n")
print(per_line[order(Chronos),
               .(CellLine, lineage, disease,
                 MYC_CN = round(MYC_CN, 2),
                 Chronos = round(Chronos, 2),
                 am_class_top, am_max_pathogenicity = round(am_max_pathogenicity, 2),
                 variants)])

cat("\nCohort regression: Chronos ~ log2(MYC_CN)\n")
fit_dt_s <- per_line[!is.na(MYC_CN) & !is.na(Chronos) & MYC_CN > 0]
if (nrow(fit_dt_s) >= 3L) {
  lm_s <- lm(Chronos ~ log2(MYC_CN), data = fit_dt_s)
  rho_s <- suppressWarnings(cor(fit_dt_s$MYC_CN, fit_dt_s$Chronos,
                                method = "spearman"))
  pear_s <- suppressWarnings(cor(log2(fit_dt_s$MYC_CN), fit_dt_s$Chronos,
                                 method = "pearson"))
  cat(sprintf("  n             : %d\n", nrow(fit_dt_s)))
  cat(sprintf("  intercept     : %+0.3f\n", unname(coef(lm_s)[1])))
  cat(sprintf("  slope (log2CN): %+0.3f  (Chronos change per CN doubling)\n",
              unname(coef(lm_s)[2])))
  cat(sprintf("  R^2           : %.3f\n", summary(lm_s)$r.squared))
  cat(sprintf("  p (slope!=0)  : %.3g\n",
              coef(summary(lm_s))[2, "Pr(>|t|)"]))
  cat(sprintf("  Pearson r     : %+0.3f  (on log2(MYC_CN))\n", pear_s))
  cat(sprintf("  Spearman rho  : %+0.3f  (rank-based, scale-free)\n", rho_s))
  cat("  Negative slope/rho confirms the same direction as the full-DepMap\n")
  cat("  trend: more amplification -> more negative Chronos = more dependent.\n")
} else {
  cat(sprintf("  Skipped (n=%d, need >= 3 lines with CN and Chronos).\n",
              nrow(fit_dt_s)))
}

cat("\n--- The puzzle: AM=likely_benign yet Chronos<-2 (highly dependent) ---\n")
puzzling <- per_line[Chronos < -2 & am_class_top == "likely benign"]
if (nrow(puzzling)) {
  print(puzzling[order(Chronos),
                 .(CellLine, lineage, disease,
                   MYC_CN = round(MYC_CN, 2),
                   Chronos = round(Chronos, 2),
                   am_max_pathogenicity = round(am_max_pathogenicity, 2),
                   variants)])
  cat(sprintf("\nMedian MYC_CN among these puzzling lines: %.2f\n",
              median(puzzling$MYC_CN, na.rm = TRUE)))
  cat(sprintf("Median MYC_CN in full DepMap: %.2f\n",
              median(bg$MYC_CN)))
  cat("Interpretation:\n")
  cat("  If MYC_CN is elevated (>~1.3) in the puzzling lines, the dependency\n")
  cat("  is being driven by amplification, not the residue change. The\n")
  cat("  AM-likely-benign call is structurally correct (those residues are\n")
  cat("  not destabilising) but biologically uninformative for an oncogene\n")
  cat("  whose activation mechanism is gene dosage / overexpression.\n")
  cat("  Burkitt-style B-cell lymphoma lines may have near-diploid CN but\n")
  cat("  high MYC mRNA via IG-MYC translocation -- CN alone cannot detect\n")
  cat("  that case (mRNA analysis would).\n")
} else {
  cat("None — every Chronos<-2 line in the cohort has at least one\n")
  cat("AM=likely_pathogenic call.\n")
}
sink()

# ---- 7. Plots --------------------------------------------------------------
base_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

xmax <- max(c(bg$MYC_CN, per_line$MYC_CN), na.rm = TRUE)

# Cohort-only regression: log2(MYC_CN) -> Chronos. Fit in log2 space because
# CN is naturally multiplicative (a doubling is the meaningful unit) and the
# x-axis is on a log2 scale. Spearman rho is reported alongside R^2 to mirror
# the non-parametric stat used for the full DepMap above.
fit_dt <- per_line[!is.na(MYC_CN) & !is.na(Chronos) & MYC_CN > 0]
n_fit <- nrow(fit_dt)
lm_subtitle <- ""
if (n_fit >= 3L) {
  lm_coh   <- lm(Chronos ~ log2(MYC_CN), data = fit_dt)
  intc_coh <- unname(coef(lm_coh)[1])
  slope_coh <- unname(coef(lm_coh)[2])
  r2_coh   <- summary(lm_coh)$r.squared
  pval_coh <- coef(summary(lm_coh))[2, "Pr(>|t|)"]
  rho_coh  <- suppressWarnings(cor(fit_dt$MYC_CN, fit_dt$Chronos,
                                   method = "spearman"))
  # Plain ASCII for the PDF device on Windows (Greek/superscript glyphs
  # don't survive mbcsToSbcs conversion in the default PDF font).
  lm_subtitle <- sprintf(
    "Cohort fit (black line): Chronos = %.2f + %.2f * log2(MYC_CN); R^2=%.2f, p=%.2g, Spearman rho=%.2f (n=%d)",
    intc_coh, slope_coh, r2_coh, pval_coh, rho_coh, n_fit)
}

p1 <- ggplot(bg, aes(MYC_CN, MYC_Chronos)) +
  geom_point(color = "grey80", alpha = 0.45, size = 0.85) +
  geom_hline(yintercept = -1, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_point(data = per_line,
             aes(MYC_CN, Chronos, color = am_class_top),
             size = 3.2, alpha = 0.9, inherit.aes = FALSE)

if (n_fit >= 3L) {
  p1 <- p1 + geom_smooth(
    data = fit_dt,
    aes(MYC_CN, Chronos),
    method = "lm", formula = y ~ log2(x),
    se = TRUE, color = "grey15", fill = "grey55",
    linewidth = 0.55, alpha = 0.18, inherit.aes = FALSE)
}

p1 <- p1 +
  geom_text(data = per_line,
            aes(MYC_CN, Chronos, label = CellLine),
            size = 2.6, hjust = -0.15, vjust = 0.4,
            check_overlap = TRUE, inherit.aes = FALSE) +
  scale_x_continuous(trans = "log2",
                     breaks = c(0.5, 1, 2, 4, 8, 16, 32),
                     limits = c(0.4, xmax * 1.6)) +
  scale_color_brewer(palette = "Set1", drop = FALSE,
                     na.value = "grey50") +
  labs(
    title = "MYC copy number vs Chronos dependency",
    subtitle = paste0(
      "Grey = full DepMap (n=", nrow(bg), ");  ",
      "colored = cohort MYC-mutant lines (n=", nrow(per_line), "). ",
      "Dashed: CN=1 (diploid) and Chronos=-1 (essentiality threshold).",
      if (nzchar(lm_subtitle)) paste0("\n", lm_subtitle) else ""
    ),
    x = "MYC relative copy number (log2 axis; 1 = diploid)",
    y = "MYC Chronos gene effect",
    color = "AM class\n(per line, max-priority)"
  ) +
  base_theme

ggsave(file.path(OUT_DIR, "11_myc_cnv_vs_chronos.pdf"),
       p1, width = 10, height = 6.5)

# Per-line bars: sort by Chronos so the most-dependent lines are at the top
per_line_b <- copy(per_line)
setorder(per_line_b, -Chronos)
per_line_b[, line_label := paste0(
  CellLine, "  (Chr ", sprintf("%+0.1f", Chronos), ")"
)]
per_line_b[, line_label := factor(line_label, levels = line_label)]

p2 <- ggplot(per_line_b,
             aes(line_label, MYC_CN, fill = am_class_top)) +
  geom_col(width = 0.78) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", MYC_CN)),
            hjust = -0.15, size = 2.9) +
  coord_flip() +
  scale_y_continuous(limits = c(0, max(per_line_b$MYC_CN, na.rm = TRUE) * 1.18)) +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    na.value = "grey60") +
  labs(
    title = "MYC copy number per cohort mutant cell line",
    subtitle = paste0(
      "Sorted by Chronos (most dependent at bottom). ",
      "Dashed line = CN=1 (diploid baseline). Bars colored by AM class."
    ),
    x = NULL,
    y = "MYC relative copy number (1 = diploid; >1 = amplified)",
    fill = "AM class"
  ) +
  base_theme

ggsave(file.path(OUT_DIR, "11_myc_cnv_per_line_bars.pdf"),
       p2, width = 9.5, height = 0.45 * nrow(per_line_b) + 2)

cat("\nWrote:\n",
    " - ", file.path(OUT_DIR, "11_myc_cnv_per_line.tsv"), "\n",
    " - ", file.path(OUT_DIR, "11_myc_cnv_summary.txt"), "\n",
    " - ", file.path(OUT_DIR, "11_myc_cnv_vs_chronos.pdf"), "\n",
    " - ", file.path(OUT_DIR, "11_myc_cnv_per_line_bars.pdf"), "\n",
    sep = "")
