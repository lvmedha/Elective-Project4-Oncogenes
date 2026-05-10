# 11_am_missense_triage_figures.R
# Part 2 (practicum): transparent AlphaMissense + hotspot triage for panel
# missense variants; tables + figures.
#
# Prerequisites: same as 09 — VCF with AM, 06_mutations_with_ACH.tsv
#
# Triage rules (evaluated in assign_triage):
#   strong_candidate — Tier-0 Hotspot OR recurrent (≥ MIN_PID_RECURRENT patients)
#                      OR AM class likely_pathogenic
#   likely_weak      — AM likely_benign OR (no class but pathogenicity ≤
#                      AM_NUM_WEAK_CUTOFF)
#   uncertain        — everything else (incl. missing AM)
#
# Outputs (results/<RUN_NAME>/):
#   11_am_missense_triage.tsv
#   11_am_triage_README.txt
#   11_am_pathogenicity_by_mut_status.pdf
#   11_am_triage_counts_stacked.pdf
#   11_am_triage_by_gene_heatmap.pdf
#   11_am_class_vs_mut_status_heatmap.pdf

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

MIN_PID_RECURRENT <- 2L
AM_NUM_WEAK_CUTOFF <- 0.2

PROJ_DIR   <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
DATA_DIR   <- file.path(PROJ_DIR, "data")
RUN_NAME   <- "ped_gof_snv"
OUT_DIR    <- file.path(PROJ_DIR, "results", RUN_NAME)
COHORT_DIR <- file.path(PROJ_DIR, "results", "cohort_full")
VCF_AM     <- file.path(COHORT_DIR, "09_vep_full_alphamissense.vcf.gz")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(VCF_AM)) {
  stop("Missing ", VCF_AM, "\nRun VEP + AlphaMissense first.")
}

target_table <- fread(file.path(DATA_DIR, "target_genes.tsv"),
                      select = c("gene", "tier", "score"))
target_genes <- target_table$gene
tg_set <- unique(target_genes)

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

status_priority <- c(
  Hotspot = 5L, Missense_Other = 4L, InframeIndel = 3L,
  Truncating = 2L, Silent = 1L
)

assign_triage <- function(is_hot, n_pat, am_class, amp) {
  n <- length(is_hot)
  out <- rep("uncertain", n)
  weak_cls <- !is.na(am_class) & am_class == "likely_benign"
  weak_num <- is.na(am_class) & is.finite(amp) & amp <= AM_NUM_WEAK_CUTOFF
  strong_rec <- n_pat >= MIN_PID_RECURRENT
  strong_am <- !is.na(am_class) & am_class == "likely_pathogenic"
  is_hot_l <- as.logical(is_hot)
  out[is_hot_l] <- "strong_candidate"
  take_weak <- !is_hot_l & (weak_cls | weak_num)
  out[take_weak] <- "likely_weak"
  take_strong <- !is_hot_l & !take_weak & (strong_rec | strong_am)
  out[take_strong] <- "strong_candidate"
  out
}

cat("Parsing AlphaMissense from VCF...\n")
am_lut <- parse_vcf_am_lut(VCF_AM, tg_set)

mut <- fread(file.path(OUT_DIR, "06_mutations_with_ACH.tsv"))
mut[, c("chrom", "pos", "ref", "alt") :=
      tstrsplit(Marker_hg38, ".", fixed = TRUE, keep = 1:4)]
mut[, pos := as.integer(pos)]
mut[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]
mut <- mut[gene_symbol %in% tg_set & Variant_Class == "Missense"]
mut[, status_rank := as.integer(status_priority[Mut_Status])]

m2 <- merge(mut, am_lut,
            by = c("chrom", "pos", "ref", "alt", "gene_symbol"),
            all.x = TRUE)

vars <- m2[, .(
  n_patients = uniqueN(PID),
  Mut_Status = Mut_Status[which.max(status_rank)],
  am_pathogenicity = {
    x <- am_pathogenicity[is.finite(am_pathogenicity)]
    if (length(x)) x[1L] else NA_real_
  },
  am_class = {
    z <- unique(am_class[nzchar(am_class) & !is.na(am_class)])
    if (length(z)) z[1L] else NA_character_
  }
), by = .(gene_symbol, Marker_hg38, HGVSp_Short)]

vars[, is_hotspot := Mut_Status == "Hotspot"]
vars[, triage_bucket := assign_triage(is_hotspot, n_patients, am_class,
                                      am_pathogenicity)]
vars[, gene_f := factor(gene_symbol, levels = target_genes)]
vars[, triage_f := factor(
  triage_bucket,
  levels = c("strong_candidate", "uncertain", "likely_weak"),
  labels = c("Strong candidate", "Uncertain", "Likely weak / tolerated")
)]
vars[, mut_lab := factor(
  Mut_Status,
  levels = c("Hotspot", "Missense_Other")
)]

vars[, has_am_score := is.finite(am_pathogenicity)]

fwrite(vars, file.path(OUT_DIR, "11_am_missense_triage.tsv"), sep = "\t")

readme <- c(
  "11_am_missense_triage — methods (auto-generated)",
  "==============================================",
  "",
  paste0("MIN_PID_RECURRENT = ", MIN_PID_RECURRENT,
         " (distinct patients at same Marker_hg38 → strong_candidate if not weak)."),
  paste0("AM_NUM_WEAK_CUTOFF = ", AM_NUM_WEAK_CUTOFF,
         " (if AM class missing, pathogenicity ≤ this → likely_weak)."),
  "",
  "Priority:",
  "  1) Hotspot (Tier-0 curated codon) → strong_candidate",
  "  2) likely_benign, or missing class with low pathogenicity → likely_weak",
  "  3) likely_pathogenic OR recurrent (≥ MIN_PID_RECURRENT patients) → strong_candidate",
  "  4) else → uncertain",
  "",
  "AlphaMissense predicts damaging vs benign missense broadly; it is not a",
  "dedicated oncogenic-activation classifier — combine with hotspots and recurrence.",
  "",
  "Row unit: unique (gene_symbol, Marker_hg38) protein-defining variant in cohort."
)
writeLines(readme, file.path(OUT_DIR, "11_am_triage_README.txt"))

base <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

v_am <- vars[is.finite(am_pathogenicity) &
               Mut_Status %in% c("Hotspot", "Missense_Other")]
if (nrow(v_am) >= 3L) {
  p1 <- ggplot(v_am, aes(mut_lab, am_pathogenicity)) +
    geom_violin(alpha = 0.35, color = NA, aes(fill = mut_lab)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, coef = 0) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.45, size = 1.6) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    labs(
      title = "AlphaMissense pathogenicity: Tier-0 hotspot vs other missense",
      subtitle = paste0(
        "Unique missense variants (n=", nrow(v_am), "); AM is not activation-specific"
      ),
      x = NULL,
      y = "AlphaMissense pathogenicity (0–1)"
    ) +
    base
  ggsave(file.path(OUT_DIR, "11_am_pathogenicity_by_mut_status.pdf"),
         p1, width = 7, height = 4.8)
} else {
  cat("Too few AM-scored variants for violin; skip 11_am_pathogenicity_by_mut_status.pdf\n")
}

tc <- vars[, .N, by = triage_f][order(triage_f)]
if (nrow(tc)) {
  p2 <- ggplot(tc, aes(triage_f, N, fill = triage_f)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = N), vjust = -0.35, size = 3.8) +
    scale_fill_manual(values = c(
      "Strong candidate" = "#2E86AB",
      "Uncertain" = "#A8ADAC",
      "Likely weak / tolerated" = "#F4A261"
    ), guide = "none") +
    labs(
      title = "Missense variant triage (cohort-unique loci)",
      subtitle = paste0(
        "Rules in 11_am_triage_README.txt; total variants = ", nrow(vars)
      ),
      x = NULL,
      y = "Count of unique variants"
    ) +
    base +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
  ggsave(file.path(OUT_DIR, "11_am_triage_counts_stacked.pdf"),
         p2, width = 7.5, height = 4.8)
}

tb <- vars[, .N, by = .(gene_f, triage_f)]
if (nrow(tb)) {
  p3 <- ggplot(tb, aes(triage_f, gene_f, fill = N)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = N), color = "grey15", size = 3) +
    scale_fill_gradient(low = "#f0f4f8", high = "#1b4965", name = "Variants") +
    labs(
      title = "Triage bucket by gene (unique missense variants)",
      x = NULL,
      y = NULL
    ) +
    base +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
  ggsave(file.path(OUT_DIR, "11_am_triage_by_gene_heatmap.pdf"),
         p3, width = 8, height = 9)
}

vars[, am_class_plot := fifelse(
  is.na(am_class) | !nzchar(am_class), "no AM class",
  gsub("_", " ", am_class, fixed = TRUE))]
vars[, am_class_plot := factor(
  am_class_plot,
  levels = c("likely pathogenic", "ambiguous", "likely benign", "no AM class")
)]

cm <- vars[, .N, by = .(mut_lab, am_class_plot)]
cm <- cm[!is.na(mut_lab)]
if (nrow(cm)) {
  p4 <- ggplot(cm, aes(am_class_plot, mut_lab, fill = N)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = N), size = 3.2) +
    scale_fill_gradient(low = "#fef9f3", high = "#c1121f", name = "Variants") +
    labs(
      title = "AlphaMissense class vs mutation tier (unique variants)",
      subtitle = "Hotspot = Tier-0 curated activating codons",
      x = "AM class (VEP plugin)",
      y = "Mutation tier"
    ) +
    base +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
  ggsave(file.path(OUT_DIR, "11_am_class_vs_mut_status_heatmap.pdf"),
         p4, width = 8.5, height = 4.5)
}

cat("\nWrote:\n",
    " - ", file.path(OUT_DIR, "11_am_missense_triage.tsv"), "\n",
    " - ", file.path(OUT_DIR, "11_am_triage_README.txt"), "\n",
    " - PDFs 11_am_*.pdf in ", OUT_DIR, "\n",
    sep = "")
