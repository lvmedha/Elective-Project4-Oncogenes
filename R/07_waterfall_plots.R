# 07_waterfall_plots.R
# Build the per-gene waterfall plots Declan asked for:
#   - one plot per target gene
#   - x = cell lines sorted by their CRISPR Chronos score for that gene
#   - y = Chronos gene effect
#   - color = Mut_Status (Hotspot / Missense_Other / InframeIndel /
#     Truncating / WT)
#   - dashed line at -1 (essentiality threshold)
# Two PDFs are written:
#   - 07_waterfalls_pancancer.pdf    (one panel per gene, all lineages)
#   - 07_waterfalls_facet_lineage.pdf (per gene, facet by OncotreePrimaryDisease)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
OUT_DIR  <- file.path(PROJ_DIR, "results")

target_genes <- c("KRAS","NRAS","HRAS","BRAF","PIK3CA",
                  "MAP2K1","EGFR","ERBB2","ALK","CTNNB1")

# ---- 1. Load mutation table with ACH attached -----------------------------
mut <- fread(file.path(OUT_DIR, "06_mutations_with_ACH.tsv"))
mut <- mut[!is.na(ACH) & gene_symbol %in% target_genes]

# Collapse to one Mut_Status per (cell line x gene), prefer the most-activating
# event if a line has multiple (Hotspot > Missense_Other > InframeIndel >
# Truncating).
status_priority <- c(Hotspot=5, Missense_Other=4, InframeIndel=3,
                     Truncating=2, Silent=1)
mut[, status_rank := status_priority[Mut_Status]]
setorder(mut, ACH, gene_symbol, -status_rank, -VAF)
mut_per_line_gene <- unique(mut, by = c("ACH", "gene_symbol"))

# ---- 2. Pull just our 10 genes from the wide CGE file ---------------------
# CGE has ~18 500 gene columns (~440 MB on disk). We only need 10 of them,
# so:
#   1. read just the header to learn what the column names look like
#   2. pick the matching columns -- DepMap names them "HUGO (entrez_id)",
#      so we strip the "(...)" part to compare against our gene list
#   3. fread() with `select = ...` to load only those 11 columns (1 ID +
#      10 genes) -- this drops memory use from ~1 GB to a few MB.
cge_header <- fread(file.path(DOC_DIR, "depmap_CGE.csv"), nrows = 0)
gene_cols  <- colnames(cge_header)[-1]
hugo_only  <- sub(" \\(.*", "", gene_cols)
sel <- gene_cols[hugo_only %in% target_genes]
cat("Found", length(sel), "matching gene columns in CGE.\n")
print(sel)

cge <- fread(file.path(DOC_DIR, "depmap_CGE.csv"),
             select = c(colnames(cge_header)[1], sel))
setnames(cge, colnames(cge_header)[1], "ModelID")
# Rename the gene columns to bare HUGO symbols so they match `target_genes`
# and our mutation table later.
old <- copy(colnames(cge))
new <- c("ModelID", sub(" \\(.*", "", old[-1]))
setnames(cge, old, new)

# Pivot from wide (one column per gene) to long (one row per cell-line x
# gene combination). Long format is what ggplot expects and makes the
# downstream join with the mutation table a simple two-key merge.
cge_long <- melt(cge, id.vars = "ModelID",
                 variable.name = "gene_symbol",
                 value.name    = "GeneEffect")
cge_long[, gene_symbol := as.character(gene_symbol)]
cat("CGE long rows:", nrow(cge_long), "\n")

# ---- 3. Add metadata + mutation info --------------------------------------
meta <- fread(file.path(DOC_DIR, "depmap_meta.csv"),
              select = c("ModelID","OncotreeLineage","OncotreePrimaryDisease",
                         "AgeCategory","StrippedCellLineName"))
plot_dt <- merge(cge_long, meta, by = "ModelID", all.x = TRUE)
plot_dt <- merge(plot_dt,
                 mut_per_line_gene[, .(ACH, gene_symbol, Mut_Status,
                                       HGVSp_Short, VAF)],
                 by.x = c("ModelID", "gene_symbol"),
                 by.y = c("ACH", "gene_symbol"),
                 all.x = TRUE)

plot_dt[is.na(Mut_Status), Mut_Status := "WT"]
plot_dt[, Mut_Status := factor(Mut_Status,
                               levels = c("Hotspot","Missense_Other",
                                          "InframeIndel","Truncating",
                                          "Silent","WT"))]

# Drop rows where the gene is missing in this line (NA gene effect)
plot_dt <- plot_dt[!is.na(GeneEffect)]

# ---- 4. Make the per-gene waterfall plot ----------------------------------
pal <- c(Hotspot        = "#D7263D",
         Missense_Other = "#F46036",
         InframeIndel   = "#9D4EDD",
         Truncating     = "#3A506B",
         Silent         = "#999999",
         WT             = "#E0E0E0")

waterfall_one_gene <- function(g, df) {
  d <- df[gene_symbol == g]
  # The waterfall ordering: sort cell lines by GeneEffect (most negative
  # first, i.e. most-dependent on the left), then assign an integer x-position
  # via .I (data.table's row-number variable). aes(x = x, ...) then draws
  # them left-to-right in that order.
  setorder(d, GeneEffect)
  d[, x := .I]
  ggplot(d, aes(x = x, y = GeneEffect, fill = Mut_Status)) +
    geom_col(width = 1) +
    # Reference lines: dashed at -1 (DepMap's standard essentiality
    # threshold), solid at 0 (no effect).
    geom_hline(yintercept = -1, linetype = "dashed",
               color = "black", linewidth = 0.4) +
    geom_hline(yintercept =  0, linewidth = 0.3, color = "black") +
    # drop = FALSE keeps every Mut_Status level in the legend even if the
    # gene has zero cell lines in that category, so legends are consistent
    # across the 10 plots.
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(title = paste0(g, "  ( n=", nrow(d), " cell lines )"),
         x = NULL, y = "Chronos gene effect", fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "right")
}

pdf(file.path(OUT_DIR, "07_waterfalls_pancancer.pdf"),
    width = 9, height = 4.5, onefile = TRUE)
for (g in target_genes) {
  print(waterfall_one_gene(g, plot_dt))
}
invisible(dev.off())
cat("\nWrote:", file.path(OUT_DIR, "07_waterfalls_pancancer.pdf"), "\n")

# ---- 5. Faceted version: same plot but split by OncotreePrimaryDisease ---
# Only show diseases with at least 5 lines for readability.
plot_dt[, disease := ifelse(is.na(OncotreePrimaryDisease) |
                            OncotreePrimaryDisease == "",
                            "Unknown", OncotreePrimaryDisease)]

waterfall_facet <- function(g, df, min_n = 5) {
  d <- df[gene_symbol == g]
  # DepMap has ~70 OncotreePrimaryDisease values; if we faceted by all of
  # them most panels would have 1-2 cell lines and be unreadable. Lump any
  # disease with < min_n lines into a single "Other" panel.
  big_groups <- d[, .N, by = disease][N >= min_n, disease]
  d[, disease_grp := ifelse(disease %in% big_groups, disease, "Other")]
  d[, disease_grp := factor(disease_grp,
                            levels = c(sort(setdiff(unique(d$disease_grp),
                                                    "Other")), "Other"))]
  # x positions restart per facet (1..N within each disease_grp) so each
  # panel renders its own waterfall with `scales = "free_x"`.
  setorder(d, disease_grp, GeneEffect)
  d[, x := seq_len(.N), by = disease_grp]
  ggplot(d, aes(x = x, y = GeneEffect, fill = Mut_Status)) +
    geom_col(width = 1) +
    geom_hline(yintercept = -1, linetype = "dashed",
               color = "black", linewidth = 0.4) +
    geom_hline(yintercept =  0, linewidth = 0.3, color = "black") +
    scale_fill_manual(values = pal, drop = FALSE) +
    facet_wrap(~ disease_grp, scales = "free_x") +
    labs(title = paste0(g, "  -  facet by OncotreePrimaryDisease"),
         x = NULL, y = "Chronos gene effect", fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          plot.title = element_text(face = "bold"),
          strip.text = element_text(size = 7),
          legend.position = "bottom")
}

pdf(file.path(OUT_DIR, "07_waterfalls_facet_lineage.pdf"),
    width = 11, height = 8, onefile = TRUE)
for (g in target_genes) {
  print(waterfall_facet(g, plot_dt))
}
invisible(dev.off())
cat("Wrote:", file.path(OUT_DIR, "07_waterfalls_facet_lineage.pdf"), "\n")

# ---- 6. Quick effect-size summary per gene -------------------------------
# A first-pass quantification of oncogene addiction: per gene, compare the
# median Chronos score in HOTSPOT-mutated cell lines vs WT lines. A more
# negative delta means stronger self-dependency in the mutated set
# (= the addiction signal Declan asked us to look for).
# This is a cheap median-difference statistic only -- formal Wilcoxon /
# Cliff's-delta tests will go in a later script.
eff <- plot_dt[, .(median_mut = median(GeneEffect[Mut_Status == "Hotspot"],
                                       na.rm = TRUE),
                   median_wt  = median(GeneEffect[Mut_Status == "WT"],
                                       na.rm = TRUE),
                   n_mut = sum(Mut_Status == "Hotspot"),
                   n_wt  = sum(Mut_Status == "WT")),
               by = gene_symbol]
eff[, delta_median := median_mut - median_wt]
setorder(eff, delta_median)
cat("\nQuick effect sizes (median Hotspot - median WT, more negative = more dependent in mutants):\n")
print(eff)
fwrite(eff, file.path(OUT_DIR, "07_effect_sizes.tsv"), sep = "\t")
