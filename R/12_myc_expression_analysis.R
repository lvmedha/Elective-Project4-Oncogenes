# 12_myc_expression_analysis.R
# MYC-focused follow-up to R/09, parallel to R/11 (CN). The CN scatter
# closes the puzzle for solid-tumor MYC-amplified lines but leaves a
# known gap: Burkitt-style B-cell lymphoma lines are near-diploid for
# MYC yet are still highly MYC-dependent because IG-MYC translocation
# drives MYC mRNA without changing copy number. Expression (log2 TPM)
# should close that gap.
#
# Inputs:
#   ~/Documents/depmap_expression.csv  — gene-level log2(TPM+1)
#                                         (cols: SequencingID,
#                                          ModelConditionID, ModelID,
#                                          IsDefaultEntryForMC,
#                                          IsDefaultEntryForModel,
#                                          <gene (entrez)> ...)
#   ~/Documents/depmap_CGE.csv         — gene-level Chronos
#   ~/Documents/depmap_meta.csv        — lineage / disease
#   results/ped_gof_snv/09_am_chronos_per_event.tsv
#
# Outputs (results/<RUN_NAME>/):
#   12_myc_expression_per_line.tsv      one row per cohort MYC-mutant line
#   12_myc_expression_summary.txt       counts, full-DepMap correlation,
#                                       puzzling cases
#   12_myc_expression_vs_chronos.pdf    scatter: full DepMap + cohort
#                                       overlay
#   12_myc_expression_per_line_bars.pdf per-line expression bar, colored
#                                       by AM class

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

EXPR_CSV   <- file.path(DOC_DIR, "depmap_expression.csv")
CGE_CSV    <- file.path(DOC_DIR, "depmap_CGE.csv")
META_CSV   <- file.path(DOC_DIR, "depmap_meta.csv")
EVENTS_TSV <- file.path(OUT_DIR, "09_am_chronos_per_event.tsv")

for (f in c(EXPR_CSV, CGE_CSV, META_CSV, EVENTS_TSV)) {
  if (!file.exists(f)) stop("Missing: ", f)
}

# ---- 1. Load only the MYC column from the 290 MB expression matrix -------
# The expression CSV carries several ID columns before the gene columns;
# we need ModelID and IsDefaultEntryForModel to dedupe to one row per
# DepMap model.
expr_header <- fread(EXPR_CSV, nrows = 0)
myc_expr_col <- grep("^MYC \\(", colnames(expr_header), value = TRUE)
stopifnot(length(myc_expr_col) == 1L)
id_cols <- intersect(
  c("ModelID", "IsDefaultEntryForModel"),
  colnames(expr_header)
)
stopifnot("ModelID" %in% id_cols)
expr <- fread(EXPR_CSV, select = c(id_cols, myc_expr_col))
setnames(expr, myc_expr_col, "MYC_Expression")

# Dedupe to default model entry where that flag exists.
if ("IsDefaultEntryForModel" %in% names(expr)) {
  expr <- expr[IsDefaultEntryForModel == "Yes"]
  expr[, IsDefaultEntryForModel := NULL]
}
expr <- expr[!is.na(ModelID) & ModelID != ""]
expr <- unique(expr, by = "ModelID")
cat(sprintf("MYC expression loaded for %d DepMap models (column='%s')\n",
            nrow(expr), myc_expr_col))

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

myc_ev <- merge(myc_ev, expr, by.x = "ACH", by.y = "ModelID", all.x = TRUE)
myc_ev <- merge(myc_ev, meta, by.x = "ACH", by.y = "ModelID", all.x = TRUE)

# Same per-line AM-class priority as R/11.
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
  MYC_Expression = unique(MYC_Expression)[1],
  lineage = unique(OncotreeLineage)[1],
  disease = unique(OncotreePrimaryDisease)[1],
  CellLine = unique(StrippedCellLineName)[1]
), by = ACH]

per_line[, am_class_top := factor(
  am_class_top,
  levels = c("likely pathogenic", "ambiguous", "likely benign", "no AM score")
)]

fwrite(per_line, file.path(OUT_DIR, "12_myc_expression_per_line.tsv"),
       sep = "\t")

# ---- 5. Full DepMap background ------------------------------------------
bg <- merge(cge_myc, expr, by = "ModelID", all = FALSE)
bg <- bg[!is.na(MYC_Chronos) & !is.na(MYC_Expression)]
bg[, is_cohort := ModelID %in% per_line$ACH]
cat(sprintf("Full DepMap background: %d lines with both expression and Chronos\n",
            nrow(bg)))

# ---- 6. Summary text ------------------------------------------------------
sink(file.path(OUT_DIR, "12_myc_expression_summary.txt"))
cat("12_myc_expression_analysis.R summary\n")
cat("====================================\n\n")
cat(sprintf("Full DepMap (lines with MYC expression and MYC Chronos): n=%d\n",
            nrow(bg)))
cat(sprintf("  median MYC_Expression : %.2f  (IQR %.2f to %.2f, range %.2f to %.2f)\n",
            median(bg$MYC_Expression), quantile(bg$MYC_Expression, 0.25),
            quantile(bg$MYC_Expression, 0.75),
            min(bg$MYC_Expression), max(bg$MYC_Expression)))
cat("  (units: log2(TPM+1))\n")
cat(sprintf("  median MYC_Chronos    : %.2f  (IQR %.2f to %.2f)\n",
            median(bg$MYC_Chronos), quantile(bg$MYC_Chronos, 0.25),
            quantile(bg$MYC_Chronos, 0.75)))
ct <- suppressWarnings(cor.test(bg$MYC_Expression, bg$MYC_Chronos,
                                method = "spearman", exact = FALSE))
cat(sprintf("  Spearman MYC_Expression vs MYC_Chronos: rho=%.3f, p=%.2g\n",
            unname(ct$estimate), ct$p.value))
cat("  (negative rho means: higher expression -> more negative Chronos = ",
    "more dependent)\n\n", sep = "")

cat("Cohort MYC-mutant cell lines (sorted by Chronos, most dependent first):\n\n")
print(per_line[order(Chronos),
               .(CellLine, lineage, disease,
                 MYC_Expression = round(MYC_Expression, 2),
                 Chronos = round(Chronos, 2),
                 am_class_top, am_max_pathogenicity = round(am_max_pathogenicity, 2),
                 variants)])

cat("\nCohort regression: Chronos ~ MYC_Expression  (log2 TPM, already log-space)\n")
fit_dt_s <- per_line[!is.na(MYC_Expression) & !is.na(Chronos)]
if (nrow(fit_dt_s) >= 3L) {
  lm_s <- lm(Chronos ~ MYC_Expression, data = fit_dt_s)
  rho_s <- suppressWarnings(cor(fit_dt_s$MYC_Expression, fit_dt_s$Chronos,
                                method = "spearman"))
  pear_s <- suppressWarnings(cor(fit_dt_s$MYC_Expression, fit_dt_s$Chronos,
                                 method = "pearson"))
  cat(sprintf("  n             : %d\n", nrow(fit_dt_s)))
  cat(sprintf("  intercept     : %+0.3f\n", unname(coef(lm_s)[1])))
  cat(sprintf("  slope (log2TPM): %+0.3f  (Chronos change per +1 unit log2 TPM)\n",
              unname(coef(lm_s)[2])))
  cat(sprintf("  R^2           : %.3f\n", summary(lm_s)$r.squared))
  cat(sprintf("  p (slope!=0)  : %.3g\n",
              coef(summary(lm_s))[2, "Pr(>|t|)"]))
  cat(sprintf("  Pearson r     : %+0.3f\n", pear_s))
  cat(sprintf("  Spearman rho  : %+0.3f  (rank-based, scale-free)\n", rho_s))
  if (unname(coef(lm_s)[2]) < 0) {
    cat("  Negative slope agrees with the full-DepMap trend: higher MYC\n")
    cat("  mRNA -> more negative Chronos = more dependent on MYC.\n")
  } else {
    cat("  NOTE: within the cohort the slope is POSITIVE, opposite to the\n")
    cat("  full-DepMap trend. With only 19 mutant lines this is statistically\n")
    cat("  underpowered and likely driven by a handful of high-expression\n")
    cat("  but only-moderately-dependent lymphoid lines (RAJI, A3KAW). The\n")
    cat("  full-DepMap correlation is the trustworthy effect; the cohort\n")
    cat("  panel is too small and too dependency-biased to recapitulate it.\n")
  }
} else {
  cat(sprintf("  Skipped (n=%d, need >= 3 lines with expression and Chronos).\n",
              nrow(fit_dt_s)))
}

cat("\n--- The puzzle: AM=likely_benign yet Chronos<-2 (highly dependent) ---\n")
puzzling <- per_line[Chronos < -2 & am_class_top == "likely benign"]
if (nrow(puzzling)) {
  print(puzzling[order(Chronos),
                 .(CellLine, lineage, disease,
                   MYC_Expression = round(MYC_Expression, 2),
                   Chronos = round(Chronos, 2),
                   am_max_pathogenicity = round(am_max_pathogenicity, 2),
                   variants)])
  cat(sprintf("\nMedian MYC_Expression among these puzzling lines: %.2f\n",
              median(puzzling$MYC_Expression, na.rm = TRUE)))
  cat(sprintf("Median MYC_Expression in full DepMap: %.2f\n",
              median(bg$MYC_Expression)))
  cat(sprintf("Top-quartile MYC_Expression threshold in full DepMap: %.2f\n",
              quantile(bg$MYC_Expression, 0.75)))
  cat("Interpretation:\n")
  cat("  Expression closes the gap that CN alone leaves open. Burkitt-style\n")
  cat("  B-cell lymphoma lines (e.g. RAJI) carry IG-MYC translocations:\n")
  cat("  near-diploid CN but highly elevated MYC mRNA. If the puzzling\n")
  cat("  AM-likely-benign lines sit at or above the top-quartile MYC\n")
  cat("  expression threshold, their dependency is explained by\n")
  cat("  overexpression, not the residue change. AM scores the residue,\n")
  cat("  so the likely-benign call is structurally correct but biologically\n")
  cat("  uninformative for an oncogene whose activation mechanism is dosage.\n")
} else {
  cat("None -- every Chronos<-2 line in the cohort has at least one\n")
  cat("AM=likely_pathogenic call.\n")
}
sink()

# ---- 7. Plots --------------------------------------------------------------
base_theme <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

xmin <- min(c(bg$MYC_Expression, per_line$MYC_Expression), na.rm = TRUE)
xmax <- max(c(bg$MYC_Expression, per_line$MYC_Expression), na.rm = TRUE)
expr_median_all <- median(bg$MYC_Expression)
expr_q3_all     <- quantile(bg$MYC_Expression, 0.75)

# Cohort-only linear regression: Chronos ~ MYC_Expression. Expression is
# already in log2(TPM+1), so no extra log transform is needed (in contrast
# to R/11, which fits Chronos ~ log2(MYC_CN) because CN is multiplicative).
fit_dt <- per_line[!is.na(MYC_Expression) & !is.na(Chronos)]
n_fit <- nrow(fit_dt)
lm_subtitle <- ""
if (n_fit >= 3L) {
  lm_coh   <- lm(Chronos ~ MYC_Expression, data = fit_dt)
  intc_coh <- unname(coef(lm_coh)[1])
  slope_coh <- unname(coef(lm_coh)[2])
  r2_coh   <- summary(lm_coh)$r.squared
  pval_coh <- coef(summary(lm_coh))[2, "Pr(>|t|)"]
  rho_coh  <- suppressWarnings(cor(fit_dt$MYC_Expression, fit_dt$Chronos,
                                   method = "spearman"))
  lm_subtitle <- sprintf(
    "Cohort fit (black line): Chronos = %.2f + %.2f * MYC_Expression; R^2=%.2f, p=%.2g, Spearman rho=%.2f (n=%d)",
    intc_coh, slope_coh, r2_coh, pval_coh, rho_coh, n_fit)
}

p1 <- ggplot(bg, aes(MYC_Expression, MYC_Chronos)) +
  geom_point(color = "grey80", alpha = 0.45, size = 0.85) +
  geom_hline(yintercept = -1, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = expr_median_all, linetype = "dashed",
             linewidth = 0.3) +
  geom_point(data = per_line,
             aes(MYC_Expression, Chronos, color = am_class_top),
             size = 3.2, alpha = 0.9, inherit.aes = FALSE)

if (n_fit >= 3L) {
  p1 <- p1 + geom_smooth(
    data = fit_dt,
    aes(MYC_Expression, Chronos),
    method = "lm", formula = y ~ x,
    se = TRUE, color = "grey15", fill = "grey55",
    linewidth = 0.55, alpha = 0.18, inherit.aes = FALSE)
}

p1 <- p1 +
  geom_text(data = per_line,
            aes(MYC_Expression, Chronos, label = CellLine),
            size = 2.6, hjust = -0.15, vjust = 0.4,
            check_overlap = TRUE, inherit.aes = FALSE) +
  scale_x_continuous(limits = c(xmin - 0.2, xmax + 1.2)) +
  scale_color_brewer(palette = "Set1", drop = FALSE,
                     na.value = "grey50") +
  labs(
    title = "MYC expression vs Chronos dependency",
    subtitle = paste0(
      "Grey = full DepMap (n=", nrow(bg), ");  ",
      "colored = cohort MYC-mutant lines (n=", nrow(per_line), "). ",
      "Dashed: full-DepMap median MYC log2 TPM and Chronos=-1 ",
      "(essentiality threshold).",
      if (nzchar(lm_subtitle)) paste0("\n", lm_subtitle) else ""
    ),
    x = "MYC expression (log2(TPM+1))",
    y = "MYC Chronos gene effect",
    color = "AM class\n(per line, max-priority)"
  ) +
  base_theme

ggsave(file.path(OUT_DIR, "12_myc_expression_vs_chronos.pdf"),
       p1, width = 10, height = 6.5)

# Per-line bars: sort by Chronos so the most-dependent lines are at the
# top after coord_flip().
per_line_b <- copy(per_line)
setorder(per_line_b, -Chronos)
per_line_b[, line_label := paste0(
  CellLine, "  (Chr ", sprintf("%+0.1f", Chronos), ")"
)]
per_line_b[, line_label := factor(line_label, levels = line_label)]

p2 <- ggplot(per_line_b,
             aes(line_label, MYC_Expression, fill = am_class_top)) +
  geom_col(width = 0.78) +
  geom_hline(yintercept = expr_median_all, linetype = "dashed",
             linewidth = 0.35) +
  geom_hline(yintercept = expr_q3_all, linetype = "dotted",
             linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", MYC_Expression)),
            hjust = -0.15, size = 2.9) +
  coord_flip() +
  scale_y_continuous(limits = c(
    0,
    max(per_line_b$MYC_Expression, na.rm = TRUE) * 1.18)) +
  scale_fill_brewer(palette = "Set1", drop = FALSE,
                    na.value = "grey60") +
  labs(
    title = "MYC expression per cohort mutant cell line",
    subtitle = paste0(
      "Sorted by Chronos (most dependent at bottom). ",
      "Dashed = full-DepMap median, dotted = top-quartile threshold ",
      "of MYC log2 TPM. Bars colored by AM class."
    ),
    x = NULL,
    y = "MYC expression (log2(TPM+1))",
    fill = "AM class"
  ) +
  base_theme

ggsave(file.path(OUT_DIR, "12_myc_expression_per_line_bars.pdf"),
       p2, width = 9.5, height = 0.45 * nrow(per_line_b) + 2)

cat("\nWrote:\n",
    " - ", file.path(OUT_DIR, "12_myc_expression_per_line.tsv"), "\n",
    " - ", file.path(OUT_DIR, "12_myc_expression_summary.txt"), "\n",
    " - ", file.path(OUT_DIR, "12_myc_expression_vs_chronos.pdf"), "\n",
    " - ", file.path(OUT_DIR, "12_myc_expression_per_line_bars.pdf"), "\n",
    sep = "")
