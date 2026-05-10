# 10_prism_mutation_sensitivity.R
# PRISM secondary readout: compare log2 AUC (drug response) between cell lines
# that carry a panel oncogene mutation vs lines without that mutation in the
# same gene ("WT for that gene" in our callset mapping).
#
# Addresses practicum item (4): targeted search for selective drug responses
# linked to activating events — including on-target inhibitors and collateral
# sensitivities (drugs whose annotated target is not the mutated gene).
#
# Inputs (defaults under ~/Documents):
#   PRISMOncologyReferenceSeqLog2AUCMatrix.csv   rows = ACH ModelID
#   PRISMOncologyReferenceSeqCompoundList.csv     SampleID = PRC-* column key
#   results/ped_gof_snv/06_mutations_with_ACH.tsv
#   data/target_genes.tsv
#
# Convention: log2 AUC is from the reference matrix; lower values typically
# indicate stronger growth inhibition / sensitivity (same framing as PedDep).
#
# Outputs (results/ped_gof_snv/):
#   10_prism_mut_vs_wt_drug_tests.tsv   all genes x drugs (min n thresholds)
#   10_prism_collateral_top_hits.tsv    FDR + ranked collateral-only drugs
#   10_prism_volcano_<GENE>.pdf        volcano; drug labels if -log10(p) > 5

# Volcano plots label drug names only when -log10(p) > LABEL_MIN_NEGLOG10P
# (ggrepel if available; else geom_text).
LABEL_MIN_NEGLOG10P <- 5

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
DATA_DIR <- file.path(PROJ_DIR, "data")
RUN_NAME <- "ped_gof_snv"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)

PRISM_AUC  <- file.path(DOC_DIR, "PRISMOncologyReferenceSeqLog2AUCMatrix.csv")
PRISM_CMP  <- file.path(DOC_DIR, "PRISMOncologyReferenceSeqCompoundList.csv")

# Lines labelled mutant if Mut_Status is one of these (excludes Silent/WT)
MUT_STATUSES <- c("Hotspot", "Missense_Other", "InframeIndel", "Truncating")

# Minimum PRISM lines per arm for Wilcoxon
MIN_MUT <- 3L
MIN_WT  <- 10L

# Genes with at least this many mutant lines in PRISM get a volcano PDF
MIN_MUT_VOLCANO <- 8L

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(PRISM_AUC) || !file.exists(PRISM_CMP)) {
  stop("Missing PRISM files under ", DOC_DIR,
       "\nExpected:\n  ", basename(PRISM_AUC),
       "\n  ", basename(PRISM_CMP))
}

target_genes <- fread(file.path(DATA_DIR, "target_genes.tsv"), select = "gene")$gene

mut <- fread(file.path(OUT_DIR, "06_mutations_with_ACH.tsv"))
mut <- mut[!is.na(ACH) & nzchar(ACH)]
mut <- mut[Mut_Status %in% MUT_STATUSES & gene_symbol %in% target_genes]
mut_lines_by_gene <- mut[, .(ModelID = unique(ACH)), by = gene_symbol]

# ---- Load PRISM long format -------------------------------------------------
pr <- fread(PRISM_AUC)
if (names(pr)[1L] == "" || is.na(names(pr)[1L])) {
  setnames(pr, 1L, "ModelID")
} else if (!"ModelID" %in% names(pr)) {
  setnames(pr, 1L, "ModelID")
}
pr[, ModelID := as.character(ModelID)]

cmp <- fread(PRISM_CMP)
setnames(cmp, "SampleID", "PRC_id", skip_absent = FALSE)

drug_cols <- setdiff(names(pr), "ModelID")
long <- melt(pr, id.vars = "ModelID", variable.name = "PRC_id",
             value.name = "log2AUC")
long[, PRC_id := as.character(PRC_id)]
long[, log2AUC := as.numeric(log2AUC)]

long <- merge(long,
              unique(cmp[, .(PRC_id, CompoundName, TargetOrMechanism,
                             GeneSymbolOfTargets)]),
              by = "PRC_id", all.x = TRUE)

all_lines <- unique(pr$ModelID)
cat("PRISM matrix: ", length(all_lines), " ACH lines, ",
    length(drug_cols), " drug columns.\n", sep = "")

# Whole-genome target list for a compound (semicolon-separated)
target_vec <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character())
  unlist(strsplit(s, ";", fixed = TRUE), use.names = FALSE)
}

targets_include_gene <- function(target_str, g) {
  tv <- target_vec(target_str)
  if (!length(tv)) return(FALSE)
  any(tv == g)
}

# Collateral = annotated targets exist and do NOT include mutated gene g
is_collateral_drug <- function(target_str, g) {
  tv <- target_vec(target_str)
  if (!length(tv)) return(NA)
  !any(tv == g)
}

# ---- Per-gene, per-drug Wilcoxon -------------------------------------------
res_list <- list()
volcano_genes <- character()

for (g in target_genes) {
  mlines <- mut_lines_by_gene[gene_symbol == g, ModelID]
  mlines <- mlines[mlines %in% all_lines]
  wlines <- setdiff(all_lines, mlines)
  if (length(mlines) < MIN_MUT || length(wlines) < MIN_WT) next

  sub <- long[ModelID %in% c(mlines, wlines)]
  sub[, arm := ifelse(ModelID %in% mlines, "mut", "wt")]

  by_drug <- sub[, {
    x <- log2AUC[arm == "mut"]
    y <- log2AUC[arm == "wt"]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    n_m <- length(x)
    n_w <- length(y)
    if (n_m < MIN_MUT || n_w < MIN_WT) {
      list(median_mut = NA_real_, median_wt = NA_real_,
           delta_mut_minus_wt = NA_real_, p_wilcox = NA_real_,
           n_mut = n_m, n_wt = n_w)
    } else {
      wt <- suppressWarnings(wilcox.test(x, y))
      med_m <- median(x)
      med_w <- median(y)
      list(median_mut = med_m, median_wt = med_w,
           delta_mut_minus_wt = med_m - med_w,
           p_wilcox = unname(wt$p.value),
           n_mut = n_m, n_wt = n_w)
    }
  }, by = .(PRC_id, CompoundName, GeneSymbolOfTargets, TargetOrMechanism)]

  by_drug[, mutated_gene := g]
  by_drug[, drug_targets_mutated_gene := vapply(
    GeneSymbolOfTargets, targets_include_gene, logical(1L), g = g)]
  by_drug[, collateral := vapply(
    GeneSymbolOfTargets, is_collateral_drug, logical(1L), g = g)]

  res_list[[g]] <- by_drug

  if (length(mlines) >= MIN_MUT_VOLCANO) {
    volcano_genes <- c(volcano_genes, g)
  }
}

if (!length(res_list)) {
  stop("No genes passed PRISM sample-size filters; check ACH overlap with ",
       "PRISM matrix.")
}

res <- rbindlist(res_list, fill = TRUE)
res[, fdr_within_gene := p.adjust(p_wilcox, method = "BH"),
    by = mutated_gene]

fwrite(res, file.path(OUT_DIR, "10_prism_mut_vs_wt_drug_tests.tsv"), sep = "\t")

# Collateral-only ranked table (exclude NA collateral; require finite p)
coll <- res[collateral == TRUE & is.finite(p_wilcox) & is.finite(delta_mut_minus_wt)]
coll[, abs_delta := abs(delta_mut_minus_wt)]
setorder(coll, mutated_gene, p_wilcox, -abs_delta)
coll[, abs_delta := NULL]
coll[, rank_within_gene := seq_len(.N), by = mutated_gene]
fwrite(coll[rank_within_gene <= 25L],
       file.path(OUT_DIR, "10_prism_collateral_top_hits.tsv"),
       sep = "\t")

cat("Wrote ", nrow(res), " gene x drug test rows.\n", sep = "")
cat("Collateral top-hits preview (first gene):\n")
print(head(coll[mutated_gene == coll[1, mutated_gene]],
            8)[, .(mutated_gene, CompoundName, delta_mut_minus_wt,
                   p_wilcox, fdr_within_gene, GeneSymbolOfTargets)])

# ---- Volcano plots (genes with enough mutant lines) -------------------------
for (g in volcano_genes) {
  dg <- res[mutated_gene == g & is.finite(p_wilcox) &
              is.finite(delta_mut_minus_wt)]
  if (nrow(dg) < 5L) next
  dg[, neglog10p := -log10(pmax(p_wilcox, 1e-15))]
  dg[, hit := fdr_within_gene < 0.25]  # lenient for exploration
  dg[, compound_lab := fifelse(is.na(CompoundName) | !nzchar(CompoundName),
                               as.character(PRC_id), CompoundName)]
  dg_lbl <- dg[neglog10p > LABEL_MIN_NEGLOG10P]

  p <- ggplot(dg, aes(delta_mut_minus_wt, neglog10p)) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    geom_point(aes(color = drug_targets_mutated_gene, shape = collateral),
               alpha = 0.85, size = 2.2) +
    scale_color_manual(
      values = c(`TRUE` = "#D7263D", `FALSE` = "#555555"),
      labels = c(`TRUE` = "Drug target includes gene",
                 `FALSE` = "Drug target ≠ gene"),
      name = NULL) +
    scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16, `NA` = 1),
                       labels = c(`TRUE` = "collateral",
                                  `FALSE` = "not collateral",
                                  `NA` = "unknown target"),
                       name = NULL) +
    labs(
      title = paste0("PRISM: mutant vs non-mutant ", g),
      subtitle = paste0(
        "delta = median(log2 AUC) mut − WT for ", g,
        "; negative ⇒ more sensitive when mutated; ",
        "labels if -log10(p) > ", LABEL_MIN_NEGLOG10P
      ),
      x = expression(Delta ~ "median log2 AUC (mut - WT)"),
      y = expression(-log[10] ~ italic(p) ~ "(Wilcoxon)")
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  if (nrow(dg_lbl) > 0L) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = dg_lbl,
        aes(label = compound_lab),
        size = 2.6,
        max.overlaps = 40,
        segment.color = "grey45",
        segment.size = 0.25,
        box.padding = 0.45,
        point.padding = 0.2,
        min.segment.length = 0,
        inherit.aes = TRUE
      )
    } else {
      message("Install package 'ggrepel' for nicer drug labels; using geom_text.")
      p <- p + geom_text(
        data = dg_lbl,
        aes(label = compound_lab),
        size = 2.4,
        check_overlap = TRUE,
        inherit.aes = TRUE
      )
    }
  }

  outfile <- file.path(OUT_DIR, paste0("10_prism_volcano_", g, ".pdf"))
  ggsave(outfile, p, width = 8.5, height = 5.5)
  cat("Wrote ", outfile, "\n", sep = "")
}

cat("\nDone. Interpret collateral hits cautiously (many drugs, BH-FDR is ",
    "within gene only).\n", sep = "")
