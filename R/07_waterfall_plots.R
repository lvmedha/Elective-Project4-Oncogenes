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
DATA_DIR <- file.path(PROJ_DIR, "data")
RUN_NAME <- "ped_gof_snv"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load target gene table from data/target_genes.tsv (built by
# R/00_build_target_genes.R from the SJPedPanel paper's GoF + pediatric
# solid-tumor filter). Also pull `tier` so we can group plots by tier.
target_table <- fread(file.path(DATA_DIR, "target_genes.tsv"),
                      select = c("gene","tier","score"))
setorder(target_table, tier, -score, gene)
target_genes <- target_table$gene
cat(sprintf("Plotting waterfalls for %d target genes (Tier counts: %s)\n",
            length(target_genes),
            paste(target_table[, .N, by = tier][order(tier),
                                                 paste0(tier, "=", N)],
                  collapse = " ")))

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

# ---- 2. Pull just our target genes from the wide CGE file -----------------
# CGE has ~18 500 gene columns (~440 MB on disk). We only need a few dozen.
#   1. read just the header to learn what the column names look like
#   2. pick the matching columns -- DepMap names them "HUGO (entrez_id)",
#      so we strip the "(...)" part to compare against our gene list
#   3. fread() with `select = ...` to load only those columns -- drops
#      memory use from ~1 GB to a few MB.
# Also handle the H3F3A -> H3-3A rename: DepMap may use either symbol.
cge_header <- fread(file.path(DOC_DIR, "depmap_CGE.csv"), nrows = 0)
gene_cols  <- colnames(cge_header)[-1]
hugo_only  <- sub(" \\(.*", "", gene_cols)
# Map HGNC-renamed legacy symbols (H3F3A -> H3-3A) so DepMap CGE's
# current symbols still match. Only inject aliases for targets we
# actually have in our list.
hgnc_aliases <- c("H3F3A" = "H3-3A", "H3F3B" = "H3-3B",
                  "HIST1H3B" = "H3C2", "HIST1H3C" = "H3C3")
target_aliases <- c(target_genes,
                    unname(hgnc_aliases[target_genes]))
target_aliases <- target_aliases[!is.na(target_aliases)]
sel <- gene_cols[hugo_only %in% target_aliases]
cat("Found", length(sel), "matching gene columns in CGE",
    "(out of", length(target_genes), "targets).\n")
missing_in_cge <- setdiff(target_genes, hugo_only)
missing_in_cge <- setdiff(missing_in_cge,
                          c(if ("H3F3A" %in% target_genes &
                                "H3-3A" %in% hugo_only) "H3F3A" else NULL))
missing_in_cge <- setdiff(missing_in_cge,
                          c(if ("HIST1H3B" %in% target_genes &
                                "H3C2" %in% hugo_only) "HIST1H3B" else NULL))
missing_in_cge <- setdiff(missing_in_cge,
                          c(if ("HIST1H3C" %in% target_genes &
                                "H3C3" %in% hugo_only) "HIST1H3C" else NULL))
if (length(missing_in_cge)) {
  cat("Targets missing in DepMap CGE (will be skipped):",
      paste(missing_in_cge, collapse = ","), "\n")
}
print(sel)

cge <- fread(file.path(DOC_DIR, "depmap_CGE.csv"),
             select = c(colnames(cge_header)[1], sel))
setnames(cge, colnames(cge_header)[1], "ModelID")
# Rename the gene columns to bare HUGO symbols so they match `target_genes`
# and our mutation table later. Also normalise the histone-rename so
# downstream code sees the legacy symbol the SJPedPanel paper uses.
old <- copy(colnames(cge))
new <- c("ModelID", sub(" \\(.*", "", old[-1]))
new[new == "H3-3A"] <- "H3F3A"
new[new == "H3-3B"] <- "H3F3B"
new[new == "H3C2"] <- "HIST1H3B"
new[new == "H3C3"] <- "HIST1H3C"
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

# Lookup "tier" by gene so titles can show it.
gene_tier <- setNames(target_table$tier, target_table$gene)

waterfall_one_gene <- function(g, df) {
  d <- df[gene_symbol == g]
  if (!nrow(d)) {
    tier_lab <- if (g %in% names(gene_tier)) paste0(" [Tier ", gene_tier[g], "]") else ""
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No DepMap Chronos (CGE) rows for this symbol.\nCheck hgnc_aliases / DepMap column names.",
                 size = 3.2, hjust = 0.5, vjust = 0.5) +
        labs(title = paste0(g, tier_lab, "  ( n=0 cell lines )"),
             x = NULL, y = "Chronos gene effect") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold"),
              panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())
    )
  }
  # The waterfall ordering: sort cell lines by GeneEffect (most negative
  # first, i.e. most-dependent on the left), then assign an integer x-position
  # via .I (data.table's row-number variable). aes(x = x, ...) then draws
  # them left-to-right in that order.
  setorder(d, GeneEffect)
  d[, x := .I]
  tier_lab <- if (g %in% names(gene_tier)) paste0(" [Tier ", gene_tier[g], "]") else ""
  ggplot(d, aes(x = x, y = GeneEffect, fill = Mut_Status)) +
    geom_col(width = 1) +
    # Reference lines: dashed at -1 (DepMap's standard essentiality
    # threshold), solid at 0 (no effect).
    geom_hline(yintercept = -1, linetype = "dashed",
               color = "black", linewidth = 0.4) +
    geom_hline(yintercept =  0, linewidth = 0.3, color = "black") +
    # drop = FALSE keeps every Mut_Status level in the legend even if the
    # gene has zero cell lines in that category, so legends are consistent
    # across all plots.
    scale_fill_manual(values = pal, drop = FALSE) +
    labs(title = paste0(g, tier_lab,
                        "  ( n=", nrow(d), " cell lines )"),
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
n_plotted <- 0L
for (g in target_genes) {
  p <- waterfall_one_gene(g, plot_dt)
  print(p)
  n_plotted <- n_plotted + 1L
}
invisible(dev.off())
cat(sprintf("\nWrote: %s (%d / %d genes (one page each))\n",
            file.path(OUT_DIR, "07_waterfalls_pancancer.pdf"),
            n_plotted, length(target_genes)))

# ---- 5. Faceted version: same plot but split by OncotreePrimaryDisease ---
# Only show diseases with at least 5 lines for readability.
plot_dt[, disease := ifelse(is.na(OncotreePrimaryDisease) |
                            OncotreePrimaryDisease == "",
                            "Unknown", OncotreePrimaryDisease)]

waterfall_facet <- function(g, df, min_n = 5) {
  d <- df[gene_symbol == g]
  if (!nrow(d)) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = "No DepMap Chronos (CGE) rows for this symbol.",
                 size = 3.2, hjust = 0.5, vjust = 0.5) +
        labs(title = paste0(g,
                            if (g %in% names(gene_tier))
                              paste0(" [Tier ", gene_tier[g], "]") else "",
                            "  -  facet by OncotreePrimaryDisease"),
             x = NULL, y = "Chronos gene effect") +
        theme_minimal(base_size = 9) +
        theme(plot.title = element_text(face = "bold"),
              panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())
    )
  }
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
    labs(title = paste0(g,
                        if (g %in% names(gene_tier))
                          paste0(" [Tier ", gene_tier[g], "]") else "",
                        "  -  facet by OncotreePrimaryDisease"),
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
n_plotted <- 0L
for (g in target_genes) {
  p <- waterfall_facet(g, plot_dt)
  print(p)
  n_plotted <- n_plotted + 1L
}
invisible(dev.off())
cat(sprintf("Wrote: %s (%d / %d genes (one page each))\n",
            file.path(OUT_DIR, "07_waterfalls_facet_lineage.pdf"),
            n_plotted, length(target_genes)))

# ---- 6. Quick effect-size summary per gene -------------------------------
# A first-pass quantification of oncogene addiction: per gene, compare the
# median Chronos score in HOTSPOT-mutated cell lines vs WT lines. Also
# include a "MutAny" comparison (Hotspot OR Missense_Other OR InframeIndel)
# vs WT, because for many of the new fusion-driven Tier-B/C oncogenes
# nobody has curated Tier-0 codons -- those genes will always have n_mut=0
# under the strict Hotspot bucket.
eff <- plot_dt[, {
  is_hot <- Mut_Status == "Hotspot"
  is_any <- Mut_Status %in% c("Hotspot","Missense_Other","InframeIndel")
  is_wt  <- Mut_Status == "WT"
  list(
    median_hot      = median(GeneEffect[is_hot], na.rm = TRUE),
    median_any_mut  = median(GeneEffect[is_any], na.rm = TRUE),
    median_wt       = median(GeneEffect[is_wt],  na.rm = TRUE),
    n_hot           = sum(is_hot),
    n_any_mut       = sum(is_any),
    n_wt            = sum(is_wt)
  )
}, by = gene_symbol]
eff <- merge(eff, target_table, by.x = "gene_symbol", by.y = "gene", all.x = TRUE)
eff[, delta_hot_vs_wt     := median_hot     - median_wt]
eff[, delta_anymut_vs_wt  := median_any_mut - median_wt]
setorder(eff, tier, delta_anymut_vs_wt)
cat("\nQuick effect sizes (more negative = more dependent in mutants):\n")
print(eff)
fwrite(eff, file.path(OUT_DIR, "07_effect_sizes.tsv"), sep = "\t")
