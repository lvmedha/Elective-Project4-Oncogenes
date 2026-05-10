# 09_am_dependency_figures.R
# Join local VEP + AlphaMissense (cohort-wide VCF) to:
#   - panel mutations with DepMap ACH (06_mutations_with_ACH.tsv)
#   - Chronos gene-effect (CGE) for the altered oncogene in that cell line
# and build manuscript-style comparisons: AM pathogenicity vs dependency,
# and AM class vs dependency, for missense events in target genes.
#
# Prerequisites:
#   - results/cohort_full/09_vep_full_alphamissense.vcf.gz (VEP 115 + AM plugin)
#   - results/ped_gof_snv/06_mutations_with_ACH.tsv
#   - ~/Documents/depmap_CGE.csv, depmap_meta.csv (same as 07_waterfall_plots.R)
#
# Outputs (results/ped_gof_snv/):
#   09_am_chronos_per_event.tsv   one row per ACH x variant x gene (deduped)
#   09_am_chronos_summary.txt     counts + Spearman correlation (missense + AM)
#   09_am_vs_chronos_scatter.pdf              pooled (all genes)
#   09_am_vs_chronos_scatter_faceted_by_gene.pdf  one panel per gene (AM points)
#   09_am_class_vs_chronos_boxplot.pdf      pooled by AM class
#   09_am_class_vs_chronos_faceted_by_gene.pdf    AM class x Chronos, facet = gene
# Part 2 triage (unique variants, no Chronos): R/11_am_missense_triage_figures.R

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
DATA_DIR <- file.path(PROJ_DIR, "data")
RUN_NAME <- "ped_gof_snv"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
COHORT_DIR <- file.path(PROJ_DIR, "results", "cohort_full")
VCF_AM <- file.path(COHORT_DIR, "09_vep_full_alphamissense.vcf.gz")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(VCF_AM)) {
  stop("Missing ", VCF_AM,
       "\nRun VEP + AlphaMissense on 08_vep_input_full.vcf first ",
       "(see scripts/run_vep_alphamissense.sh).")
}

target_table <- fread(file.path(DATA_DIR, "target_genes.tsv"),
                      select = c("gene", "tier", "score"))
target_genes <- target_table$gene
tg_set <- unique(target_genes)

# ---- Parse VCF CSQ into variant x gene AlphaMissense ------------------------
# CSQ Format (VEP 115): ...|SYMBOL|...|BIOTYPE|...|am_class|am_pathogenicity
# Indices 1-based in R after strsplit(..., "|") : SYMBOL=4, BIOTYPE=8,
# Consequence=2, am_class=24, am_pathogenicity=25

extract_info_csq <- function(info) {
  parts <- strsplit(info, ";", fixed = TRUE)[[1]]
  hit <- parts[startsWith(parts, "CSQ=")]
  if (length(hit) == 0L) return(NA_character_)
  sub("^CSQ=", "", hit[1L])
}

parse_vcf_am_lut <- function(path, genes_allow) {
  con <- gzfile(path, "r")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("VCF: unexpected EOF before #CHROM")
    if (startsWith(line, "#CHROM")) break
  }
  rows <- list()
  n <- 0L
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (!nzchar(line)) next
    f <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(f) < 8L) next
    chrom <- f[1L]
    pos   <- as.integer(f[2L])
    ref   <- f[4L]
    alt   <- f[5L]
    info  <- f[8L]
    csq_raw <- extract_info_csq(info)
    if (is.na(csq_raw)) next
    for (rec in strsplit(csq_raw, ",", fixed = FALSE)[[1]]) {
      g <- strsplit(rec, "|", fixed = TRUE)[[1]]
      if (length(g) < 25L) next
      cons <- g[2L]
      sym  <- g[4L]
      bio  <- g[8L]
      amc  <- g[24L]
      amp  <- g[25L]
      if (!grepl("missense", cons, fixed = TRUE)) next
      if (bio != "protein_coding") next
      if (!sym %in% genes_allow) next
      amp_num <- if (nzchar(amp)) suppressWarnings(as.numeric(amp)) else NA_real_
      n <- n + 1L
      rows[[n]] <- list(
        chrom = chrom, pos = pos, ref = ref, alt = alt,
        gene_symbol = sym,
        am_class = if (nzchar(amc)) amc else NA_character_,
        am_pathogenicity = amp_num
      )
    }
  }
  if (n == 0L) {
    return(data.table(
      chrom = character(), pos = integer(), ref = character(),
      alt = character(), gene_symbol = character(),
      am_class = character(), am_pathogenicity = numeric()
    ))
  }
  dt_all <- rbindlist(rows)
  setorder(dt_all, chrom, pos, ref, alt, gene_symbol, -am_pathogenicity)
  unique(dt_all, by = c("chrom", "pos", "ref", "alt", "gene_symbol"))
}

cat("Parsing AlphaMissense from VCF (target genes only)...\n")
am_lut <- parse_vcf_am_lut(VCF_AM, tg_set)
cat("  AM lookup rows (unique locus x gene):", nrow(am_lut), "\n")

# ---- Mutations + locus keys -------------------------------------------------
mut <- fread(file.path(OUT_DIR, "06_mutations_with_ACH.tsv"))
mut[, c("chrom", "pos", "ref", "alt") :=
      tstrsplit(Marker_hg38, ".", fixed = TRUE, keep = 1:4)]
mut[, pos := as.integer(pos)]
mut[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]

mut <- mut[gene_symbol %in% tg_set]
mut_miss <- mut[Variant_Class == "Missense"]

m2 <- merge(mut_miss, am_lut,
            by = c("chrom", "pos", "ref", "alt", "gene_symbol"),
            all.x = TRUE)
cat("Missense rows in panel:", nrow(mut_miss),
    "  with AM join:", sum(!is.na(m2$am_pathogenicity)), "\n")

# ---- CGE long (same column selection as 07) ---------------------------------
hgnc_aliases <- c(
  "H3F3A" = "H3-3A", "H3F3B" = "H3-3B",
  "HIST1H3B" = "H3C2", "HIST1H3C" = "H3C3"
)
cge_header <- fread(file.path(DOC_DIR, "depmap_CGE.csv"), nrows = 0)
gene_cols  <- colnames(cge_header)[-1]
hugo_only  <- sub(" \\(.*", "", gene_cols)
target_aliases <- c(target_genes, unname(hgnc_aliases[target_genes]))
target_aliases <- target_aliases[!is.na(target_aliases)]
sel <- gene_cols[hugo_only %in% target_aliases]
cge <- fread(file.path(DOC_DIR, "depmap_CGE.csv"),
             select = c(colnames(cge_header)[1], sel))
setnames(cge, colnames(cge_header)[1], "ModelID")
old <- copy(colnames(cge))
new <- c("ModelID", sub(" \\(.*", "", old[-1]))
new[new == "H3-3A"] <- "H3F3A"
new[new == "H3-3B"] <- "H3F3B"
new[new == "H3C2"] <- "HIST1H3B"
new[new == "H3C3"] <- "HIST1H3C"
setnames(cge, old, new)
cge_long <- melt(cge, id.vars = "ModelID",
                 variable.name = "gene_symbol", value.name = "GeneEffect")
cge_long[, gene_symbol := as.character(gene_symbol)]

# One observation per cell line x variant (Chronos is per gene per line)
m2 <- m2[!is.na(ACH) & nzchar(ACH)]
plot_dt <- merge(m2, cge_long,
                 by.x = c("ACH", "gene_symbol"),
                 by.y = c("ModelID", "gene_symbol"),
                 all.x = TRUE)
plot_dt <- plot_dt[!is.na(GeneEffect)]

plot_dt[, am_class_lab := fifelse(
  is.na(am_class) | !nzchar(am_class),
  "no AM score",
  gsub("_", " ", am_class, fixed = TRUE))]
plot_dt[, am_class_lab := factor(
  am_class_lab,
  levels = c("likely pathogenic", "ambiguous", "likely benign", "no AM score"))]

# Dedupe: same ACH + variant + gene
plot_u <- unique(plot_dt, by = c("ACH", "gene_symbol", "Marker_hg38"))
plot_u[, gene_f := factor(gene_symbol, levels = target_genes)]

fwrite(plot_u,
       file.path(OUT_DIR, "09_am_chronos_per_event.tsv"),
       sep = "\t")

# ---- Summary stats ---------------------------------------------------------
sink(file.path(OUT_DIR, "09_am_chronos_summary.txt"))
cat("09_am_dependency_figures.R summary\n")
cat("Rows (deduped ACH x variant x gene) with CGE:", nrow(plot_u), "\n")
cat("With numeric AlphaMissense:", sum(!is.na(plot_u$am_pathogenicity)), "\n")
cat("\nPer Mut_Status:\n")
print(plot_u[, .N, by = Mut_Status][order(-N)])
cat("\nPer gene (n points):\n")
print(plot_u[, .N, by = gene_symbol][order(-N)])

pw <- plot_u[!is.na(am_pathogenicity)]
if (nrow(pw) >= 5L) {
  ct <- suppressWarnings(cor.test(pw$am_pathogenicity, pw$GeneEffect,
                                  method = "spearman",
                                  exact = FALSE))
  cat("\nSpearman rho (AM pathogenicity vs Chronos), missense + AM, n=",
      nrow(pw), ":\n", sep = "")
  print(ct)
  cat(
    "\nInterpretation: more negative Chronos = stronger dependency.\n",
    "Negative rho means higher AM pathogenicity tracks with stronger ",
    "dependency.\n", sep = "")
} else {
  cat("\nToo few points with AM for stable Spearman (need >= 5).\n")
}
sink()

# ---- Figures ---------------------------------------------------------------
base <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# 1) Scatter (x requires a missense AM score)
plot_sc <- plot_u[!is.na(am_pathogenicity)]
if (nrow(plot_sc) >= 2L) {
  p1 <- ggplot(plot_sc, aes(am_pathogenicity, GeneEffect)) +
    geom_hline(yintercept = -1, linetype = "dashed", linewidth = 0.35) +
    geom_point(aes(color = Mut_Status), alpha = 0.75, size = 2.2) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "AlphaMissense vs Chronos (altered oncogene)",
      subtitle = paste0(
        "Missense + AM score; ", RUN_NAME,
        "; ACH lines with CGE (n=", nrow(plot_sc), ")"
      ),
      x = "AlphaMissense pathogenicity (0–1)",
      y = "Chronos gene effect (dependency)",
      color = "Mut_Status"
    ) +
    base
  if (nrow(plot_sc) >= 8L) {
    p1 <- p1 + geom_smooth(method = "lm", se = TRUE, linewidth = 0.35,
                            alpha = 0.12, color = "grey35")
  }
  ggsave(file.path(OUT_DIR, "09_am_vs_chronos_scatter.pdf"),
         p1, width = 8, height = 5.5)
} else {
  cat("Too few points with AM for scatter PDF; skipped.\n")
}

# 2) Boxplot by AM class (only rows with known class; add no-score separately)
plot_box <- copy(plot_u)
plot_box[, am_box := as.character(am_class_lab)]
p2 <- ggplot(plot_box, aes(am_box, GeneEffect)) +
  geom_hline(yintercept = -1, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  geom_jitter(aes(color = Mut_Status), width = 0.12, height = 0,
              alpha = 0.7, size = 1.8) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Chronos by AlphaMissense class",
    x = "AlphaMissense class",
    y = "Chronos gene effect",
    color = "Mut_Status"
  ) +
  base +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(OUT_DIR, "09_am_class_vs_chronos_boxplot.pdf"),
       p2, width = 8, height = 5.5)

ncol_f <- 4L

# 3) Scatter faceted by gene (every gene with ≥1 AM-scored point; free scales)
plot_sc <- plot_u[!is.na(am_pathogenicity)]
if (nrow(plot_sc) >= 1L) {
  ng <- uniqueN(plot_sc$gene_symbol)
  p3 <- ggplot(plot_sc, aes(am_pathogenicity, GeneEffect)) +
    geom_hline(yintercept = -1, linetype = "dashed", linewidth = 0.25) +
    geom_point(aes(color = Mut_Status), alpha = 0.8, size = 1.8) +
    facet_wrap(~ gene_f, scales = "free", ncol = ncol_f, drop = TRUE) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "AlphaMissense vs Chronos (faceted by gene)",
      subtitle = paste0(
        RUN_NAME,
        "; one panel per gene with ≥1 missense + AM; free x/y scales"
      ),
      x = "AlphaMissense pathogenicity",
      y = "Chronos gene effect",
      color = "Mut_Status"
    ) +
    base +
    theme(strip.text = element_text(size = 7))
  ggsave(
    file.path(OUT_DIR, "09_am_vs_chronos_scatter_faceted_by_gene.pdf"),
    p3,
    width = 11,
    height = 2 + 2.5 * ceiling(ng / ncol_f)
  )
} else {
  cat("No AM-scored points; skipping faceted-by-gene scatter.\n")
}

# 4) AM-class boxplot faceted by gene (shows within-gene AM vs dependency spread)
plot_box2 <- copy(plot_u)
plot_box2[, am_box := as.character(am_class_lab)]
p4 <- ggplot(plot_box2, aes(am_box, GeneEffect)) +
  geom_hline(yintercept = -1, linetype = "dashed", linewidth = 0.25) +
  geom_boxplot(outlier.shape = NA, coef = 0, linewidth = 0.35) +
  geom_jitter(aes(color = Mut_Status), width = 0.1, height = 0,
              alpha = 0.65, size = 1.2) +
  facet_wrap(~ gene_f, scales = "free_y", ncol = ncol_f, drop = TRUE) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Chronos by AlphaMissense class (faceted by gene)",
    x = "AlphaMissense class",
    y = "Chronos gene effect",
    color = "Mut_Status"
  ) +
  base +
  theme(
    strip.text = element_text(size = 7),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 6)
  )
ng2 <- uniqueN(plot_box2$gene_symbol)
ggsave(
  file.path(OUT_DIR, "09_am_class_vs_chronos_faceted_by_gene.pdf"),
  p4,
  width = 11,
  height = 2 + 2.6 * ceiling(ng2 / ncol_f)
)

cat("\nWrote:\n",
    " - ", file.path(OUT_DIR, "09_am_chronos_per_event.tsv"), "\n",
    " - ", file.path(OUT_DIR, "09_am_chronos_summary.txt"), "\n",
    " - ", file.path(OUT_DIR, "09_am_vs_chronos_scatter.pdf"),
    " (if enough AM points)\n",
    " - ", file.path(OUT_DIR, "09_am_vs_chronos_scatter_faceted_by_gene.pdf"),
    "\n",
    " - ", file.path(OUT_DIR, "09_am_class_vs_chronos_boxplot.pdf"), "\n",
    " - ", file.path(OUT_DIR, "09_am_class_vs_chronos_faceted_by_gene.pdf"),
    "\n",
    sep = "")
