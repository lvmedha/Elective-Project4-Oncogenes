# 08_make_vep_input.R
# Convert the filtered/tiered 42-gene mutation table into VEP-ready inputs
# so we can run a local VEP with --plugin AlphaMissense and pull AM scores
# back into the tiering step (R/05).
#
# Inputs  (default):
#   results/<RUN_NAME>/05_mutations_tiered.tsv   (post-filter, post-tier;
#                                                 1 row per PID x gene)
#
# Outputs (written to results/<RUN_NAME>/):
#   08_vep_input.vcf                  proper VCFv4.2 (sortable, bcftools-able)
#   08_vep_input.tsv                  VEP "default" format (CHROM POS . REF ALT . . .)
#   08_vep_input_README.txt           run instructions (local VEP + AM plugin)
#
# Genome: GRCh38 / hg38. Chromosome names are written WITHOUT the "chr" prefix
# to match the Ensembl VEP cache (most common). If your local cache uses RefSeq
# style ("chr1", "chrX") add `--chr_synonyms` or rename with bcftools annotate.
#
# Normalisation: the source St. Jude markers are already in VCF anchor-base
# representation (e.g. "AA -> AGCTCCAA" for an insertion). For STRICT left-
# alignment + minimal-representation normalisation, run bcftools after this
# script:
#   bcftools norm -f /path/to/GRCh38.fa -c w 08_vep_input.vcf \
#     -Oz -o 08_vep_input.norm.vcf.gz
# Then point VEP at the normalised .vcf.gz instead. This script does not call
# bcftools because it's a Linux/WSL tool; the unnormalised VCF still works
# with VEP for the variants in this dataset (all source coords are GRCh38).

suppressPackageStartupMessages({
  library(data.table)
})

PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
RUN_NAME <- "all_42genes"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)

# Source of "filtered" mutations. Switch to 03_vaf_annotated.tsv if you want
# the broader pre-tier coding set, or to 05_mutations_tiered_long.tsv if you
# want one row per (PID, variant) including silent calls.
IN_FILE  <- file.path(OUT_DIR, "05_mutations_tiered.tsv")
stopifnot(file.exists(IN_FILE))

mut <- fread(IN_FILE)
cat(sprintf("Loaded %d mutation rows from %s\n", nrow(mut), basename(IN_FILE)))

# -----------------------------------------------------------------------------
# 1. Parse Marker_hg38 ("chr12.25245350.C.T") into VCF columns and dedup to
#    unique genomic variants. We aggregate per-PID context as INFO so it can
#    be carried through VEP without losing the link back to the tiered table.
# -----------------------------------------------------------------------------
mut[, c("chrom_raw", "pos", "ref", "alt") :=
      tstrsplit(Marker_hg38, ".", fixed = TRUE, keep = 1:4)]
mut[, pos := as.integer(pos)]

# Strip "chr" prefix -> Ensembl GRCh38 cache convention
mut[, chrom := sub("^chr", "", chrom_raw)]

# Validate allele alphabets (A/C/G/T/N only, no IUPAC ambiguity)
bad_alleles <- mut[!grepl("^[ACGTNacgtn]+$", ref) |
                   !grepl("^[ACGTNacgtn]+$", alt)]
if (nrow(bad_alleles)) {
  cat("WARNING:", nrow(bad_alleles),
      "rows have non-ACGTN alleles -- dropping:\n")
  print(bad_alleles[, .(Marker_hg38, ref, alt)])
  mut <- mut[grepl("^[ACGTN]+$", ref, ignore.case = TRUE) &
             grepl("^[ACGTN]+$", alt, ignore.case = TRUE)]
}
mut[, ref := toupper(ref)]
mut[, alt := toupper(alt)]

# Aggregate per-variant: list of carrier PIDs, dominant gene/HGVSp, status mix
per_var <- mut[, .(
  gene       = paste(sort(unique(gene_symbol)), collapse = ","),
  hgvsp      = paste(sort(unique(HGVSp_Short)), collapse = ","),
  mut_status = paste(sort(unique(Mut_Status)),  collapse = ","),
  n_pids     = uniqueN(PID),
  vclass     = paste(sort(unique(Variant_Class)), collapse = ",")
), by = .(chrom, pos, ref, alt)]

cat(sprintf("Deduplicated to %d unique genomic variants (from %d patient rows)\n",
            nrow(per_var), nrow(mut)))

# -----------------------------------------------------------------------------
# 2. Sort by chromosome (1..22, X, Y, MT) and position.
# -----------------------------------------------------------------------------
chrom_order <- c(as.character(1:22), "X", "Y", "MT", "M")
per_var[, chrom_idx := match(chrom, chrom_order)]
if (anyNA(per_var$chrom_idx)) {
  cat("WARNING: unknown chromosomes (will sort to end):",
      paste(unique(per_var[is.na(chrom_idx), chrom]), collapse = ", "), "\n")
}
setorder(per_var, chrom_idx, pos, ref, alt, na.last = TRUE)
per_var[, chrom_idx := NULL]

# -----------------------------------------------------------------------------
# 3. Write VCFv4.2.
# -----------------------------------------------------------------------------
vcf_path <- file.path(OUT_DIR, "08_vep_input.vcf")

# VCF INFO encoder: replace VCF-illegal characters in free-text fields.
vcf_escape <- function(x) {
  x <- gsub(";", "|", x, fixed = TRUE)
  x <- gsub(" ", "_", x, fixed = TRUE)
  x <- gsub(",", "/", x, fixed = TRUE)
  x
}

contigs_used <- unique(per_var$chrom)
contigs_used <- contigs_used[order(match(contigs_used, chrom_order))]

vcf_header <- c(
  "##fileformat=VCFv4.2",
  sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
  "##source=Elective-Project4-Oncogenes/R/08_make_vep_input.R",
  "##reference=GRCh38",
  sprintf("##contig=<ID=%s>", contigs_used),
  '##INFO=<ID=GENE,Number=1,Type=String,Description="Target oncogene (gene_symbol from R/03)">',
  '##INFO=<ID=HGVSP,Number=1,Type=String,Description="HGVSp_Short from R/03 (comma-separated if multi-tx)">',
  '##INFO=<ID=MUT_STATUS,Number=1,Type=String,Description="Tiered Mut_Status from R/05">',
  '##INFO=<ID=VCLASS,Number=1,Type=String,Description="Variant_Class from R/03">',
  '##INFO=<ID=N_PIDS,Number=1,Type=Integer,Description="Number of patients carrying this variant in 05_mutations_tiered.tsv">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

per_var[, info := sprintf(
  "GENE=%s;HGVSP=%s;MUT_STATUS=%s;VCLASS=%s;N_PIDS=%d",
  vcf_escape(gene),
  vcf_escape(hgvsp),
  vcf_escape(mut_status),
  vcf_escape(vclass),
  n_pids
)]

vcf_body <- per_var[, sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t%s",
                              chrom, pos, ref, alt, info)]

writeLines(c(vcf_header, vcf_body), vcf_path)
cat("Wrote", vcf_path, "\n")

# -----------------------------------------------------------------------------
# 4. Write VEP "default" format (space-separated, no header). VEP autodetects
#    this when fed via `--input_file ... --format ensembl`. Matches the
#    `vep_input` column already used in R/03_annotate_variants.R.
# -----------------------------------------------------------------------------
default_path <- file.path(OUT_DIR, "08_vep_input.tsv")
default_body <- per_var[, sprintf("%s %d . %s %s . . .", chrom, pos, ref, alt)]
writeLines(default_body, default_path)
cat("Wrote", default_path, "\n")

# -----------------------------------------------------------------------------
# 5. Tiny README with the local-VEP + AlphaMissense command.
# -----------------------------------------------------------------------------
readme_path <- file.path(OUT_DIR, "08_vep_input_README.txt")
readme <- c(
  "08_vep_input.vcf / 08_vep_input.tsv -- VEP-ready inputs",
  "========================================================",
  "",
  sprintf("Source : %s", IN_FILE),
  sprintf("Variants: %d unique (chrom, pos, ref, alt) on GRCh38",
          nrow(per_var)),
  "Genome  : GRCh38 / hg38, Ensembl chromosome naming (no 'chr' prefix)",
  "",
  "Files:",
  "  08_vep_input.vcf  -- VCFv4.2 with INFO fields GENE / HGVSP / MUT_STATUS / VCLASS / N_PIDS",
  "  08_vep_input.tsv  -- VEP default format (CHROM POS . REF ALT . . .)",
  "",
  "Optional strict normalisation (recommended; run inside WSL):",
  "  bcftools norm -f /path/to/GRCh38.fa -c w 08_vep_input.vcf \\",
  "    -Oz -o 08_vep_input.norm.vcf.gz",
  "  tabix -p vcf 08_vep_input.norm.vcf.gz",
  "",
  "Local VEP + AlphaMissense plugin (run inside the alpha_tools conda env):",
  "  vep \\",
  "    --input_file  08_vep_input.vcf \\",
  "    --output_file 08_vep_alphamissense.tsv \\",
  "    --species homo_sapiens --assembly GRCh38 \\",
  "    --cache --offline --dir_cache /path/to/vep_cache \\",
  "    --fasta /path/to/GRCh38.fa \\",
  "    --tab --everything --pick \\",
  "    --plugin AlphaMissense,file=/path/to/AlphaMissense_hg38.tsv.gz",
  "",
  "Notes:",
  "  - The AlphaMissense TSV must be bgzipped + tabix-indexed:",
  "      tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz",
  "  - --pick keeps one consequence per variant; drop it for all transcripts.",
  "  - --everything turns on canonical, mane, hgvs, protein, biotype, etc.",
  "  - If your VEP cache uses 'chr1' style names, add `--chr_synonyms` or",
  "    rename contigs first (bcftools annotate --rename-chrs)."
)
writeLines(readme, readme_path)
cat("Wrote", readme_path, "\n")

cat("\nDone. Inputs ready in", OUT_DIR, "\n")
