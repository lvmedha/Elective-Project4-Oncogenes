# 08_make_vep_input.R
# Convert the WHOLE VAF dataset into VEP-ready inputs so a local VEP run
# with --plugin AlphaMissense (+ LOFTEE + gnomAD on the St. Jude cluster)
# can score every coding variant in the cohort and we can pull AM scores
# back into the downstream analysis.
#
# This is panel-agnostic: every coding variant in the cohort gets an
# AlphaMissense score, regardless of whether the gene is in the 42-gene
# target list. AlphaMissense only scores missense substitutions; VEP will
# still annotate non-missense variants normally and the AM column will
# simply be empty for those rows. The resulting cohort-wide table can be
# joined back onto any panel-specific analysis (e.g.
# results/ped_gof_snv/03_vaf_annotated.tsv) by (chrom, pos, ref, alt).
#
# Inputs:
#   ~/Documents/AllMarkers_VAF_long.tsv       (raw long-format VAF table)
#
# Outputs (results/cohort_full/):
#   08_vep_input_full.vcf                     VCFv4.2 (Part 3.2 of the
#                                             cluster protocol). chrom
#                                             names use "chr" prefix so
#                                             they match the cluster's
#                                             VEP cache convention. INFO
#                                             carries N_PIDS;VCLASS for
#                                             traceability.
#   08_vep_input_full.ensembl                 VEP "default" format (Part
#                                             3.1). Tab-separated:
#                                             chrom start end allele
#                                             strand identifier
#                                             with proper indel encoding
#                                             (-/INS, DEL/-, etc.).
#   08_vep_input_full_README.txt              cluster-workflow notes.
#
# Genome: GRCh38 / hg38. Chromosome names INCLUDE the "chr" prefix per
# the cluster protocol (the cluster's VEP cache uses chr1..chr22, chrX,
# chrY, chrM). If you ever target an Ensembl-style cache instead, strip
# the prefix or pass `--synonyms /path/to/chr_synonyms.txt` to vep.

suppressPackageStartupMessages({
  library(data.table)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
OUT_DIR  <- file.path(PROJ_DIR, "results", "cohort_full")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

IN_FILE  <- file.path(DOC_DIR, "AllMarkers_VAF_long.tsv")
stopifnot(file.exists(IN_FILE))

vaf <- fread(IN_FILE)
cat(sprintf("Loaded %d rows from %s\n", nrow(vaf), basename(IN_FILE)))

# Drop SVs (no single chrom/pos/ref/alt). Keep SNV + Indel + MNV.
keep_types <- c("SNV", "Indel", "MNV")
vaf_var <- vaf[Marker.Type %in% keep_types]
cat(sprintf("After keeping only %s: %d rows\n",
            paste(keep_types, collapse = "/"), nrow(vaf_var)))

# Parse Marker_hg38 ("chr12.25245347.C.T") into chrom/pos/ref/alt.
vaf_var[, c("chrom_raw", "pos", "ref", "alt") :=
          tstrsplit(Marker_hg38, ".", fixed = TRUE, keep = 1:4)]
vaf_var[, pos := as.integer(pos)]

# Cluster protocol uses "chr"-prefixed names. Normalise so every row
# starts with "chr" exactly once.
vaf_var[, chrom := ifelse(startsWith(chrom_raw, "chr"),
                          chrom_raw,
                          paste0("chr", chrom_raw))]

# Validate allele alphabets (A/C/G/T/N only, no IUPAC ambiguity codes).
bad_alleles <- vaf_var[!grepl("^[ACGTNacgtn]+$", ref) |
                       !grepl("^[ACGTNacgtn]+$", alt)]
if (nrow(bad_alleles)) {
  cat(sprintf("WARNING: dropping %d rows with non-ACGTN alleles\n",
              nrow(bad_alleles)))
  vaf_var <- vaf_var[grepl("^[ACGTN]+$", ref, ignore.case = TRUE) &
                     grepl("^[ACGTN]+$", alt, ignore.case = TRUE)]
}
vaf_var[, ref := toupper(ref)]
vaf_var[, alt := toupper(alt)]

# Aggregate per unique (chrom, pos, ref, alt).
per_var <- vaf_var[, .(
  n_pids = uniqueN(PID),
  vclass = paste(sort(unique(Marker.Type)), collapse = ",")
), by = .(chrom, pos, ref, alt)]
cat(sprintf("Unique genomic variants: %d (from %d patient rows)\n",
            nrow(per_var), nrow(vaf_var)))

# Sort by chromosome (chr1..chr22, chrX, chrY, chrM) then position.
chrom_order <- paste0("chr", c(as.character(1:22), "X", "Y", "M", "MT"))
per_var[, chrom_idx := match(chrom, chrom_order)]
if (anyNA(per_var$chrom_idx)) {
  cat("WARNING: unknown chromosomes (will sort to end):",
      paste(unique(per_var[is.na(chrom_idx), chrom]), collapse = ", "), "\n")
}
setorder(per_var, chrom_idx, pos, ref, alt, na.last = TRUE)
per_var[, chrom_idx := NULL]

# -----------------------------------------------------------------------------
# Output 1: VCFv4.2 (cluster protocol Part 3.2). The .lsf job feeds this in
# directly: `vep -i 08_vep_input_full.vcf` -- no --format flag needed.
# -----------------------------------------------------------------------------
vcf_path <- file.path(OUT_DIR, "08_vep_input_full.vcf")

contigs_used <- unique(per_var$chrom)
contigs_used <- contigs_used[order(match(contigs_used, chrom_order))]

vcf_header <- c(
  "##fileformat=VCFv4.2",
  sprintf("##fileDate=%s", format(Sys.Date(), "%Y%m%d")),
  "##source=Elective-Project4-Oncogenes/R/08_make_vep_input.R",
  "##reference=GRCh38",
  sprintf("##contig=<ID=%s>", contigs_used),
  '##INFO=<ID=N_PIDS,Number=1,Type=Integer,Description="Number of patients carrying this variant in AllMarkers_VAF_long.tsv">',
  '##INFO=<ID=VCLASS,Number=1,Type=String,Description="Marker.Type from source VAF (SNV/Indel/MNV)">',
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
)

per_var[, info := sprintf("N_PIDS=%d;VCLASS=%s", n_pids, vclass)]
vcf_body <- per_var[, sprintf("%s\t%d\t.\t%s\t%s\t.\tPASS\t%s",
                              chrom, pos, ref, alt, info)]
writeLines(c(vcf_header, vcf_body), vcf_path)
cat("Wrote", vcf_path, "\n")

# -----------------------------------------------------------------------------
# Output 2: VEP "default" / Ensembl format (cluster protocol Part 3.1).
# Tab-separated, no header, columns: chrom start end allele strand identifier.
# Indel encoding follows the protocol's Python snippet exactly:
#
#   SNV       : start = end = pos             allele = REF/ALT
#   simple del: start = pos+1, end = pos+|REF|-1   allele = REF[1:]/-
#   simple ins: start = pos,   end = pos+1         allele = -/ALT[1:]
#   complex   : start = pos,   end = pos+|REF|-1   allele = REF/ALT
#
# Use this if you ever pass --format ensembl to vep instead of feeding the VCF.
# -----------------------------------------------------------------------------
encode_vep_default <- function(pos, ref, alt) {
  pos <- as.integer(pos)
  rl  <- nchar(ref)
  al  <- nchar(alt)

  is_snv <- rl == 1L & al == 1L
  is_del <- rl >  al & substr(alt, 1L, 1L) == substr(ref, 1L, 1L) & al == 1L
  is_ins <- al >  rl & substr(ref, 1L, 1L) == substr(alt, 1L, 1L) & rl == 1L

  start  <- integer(length(pos))
  end    <- integer(length(pos))
  allele <- character(length(pos))

  start[is_snv] <- pos[is_snv]
  end  [is_snv] <- pos[is_snv]
  allele[is_snv] <- paste0(ref[is_snv], "/", alt[is_snv])

  start[is_del] <- pos[is_del] + 1L
  end  [is_del] <- pos[is_del] + rl[is_del] - 1L
  allele[is_del] <- paste0(substring(ref[is_del], 2L), "/-")

  start[is_ins] <- pos[is_ins]
  end  [is_ins] <- pos[is_ins] + 1L
  allele[is_ins] <- paste0("-/", substring(alt[is_ins], 2L))

  cx <- !is_snv & !is_del & !is_ins
  start[cx] <- pos[cx]
  end  [cx] <- pos[cx] + rl[cx] - 1L
  allele[cx] <- paste0(ref[cx], "/", alt[cx])

  list(start = start, end = end, allele = allele)
}

enc <- encode_vep_default(per_var$pos, per_var$ref, per_var$alt)
per_var[, c("vep_start", "vep_end", "vep_allele") :=
          list(enc$start, enc$end, enc$allele)]
per_var[, vep_id := sprintf("var%d", .I)]

ensembl_path <- file.path(OUT_DIR, "08_vep_input_full.ensembl")
ensembl_body <- per_var[, sprintf("%s\t%d\t%d\t%s\t+\t%s",
                                  chrom, vep_start, vep_end,
                                  vep_allele, vep_id)]
writeLines(ensembl_body, ensembl_path)
cat("Wrote", ensembl_path, "\n")

# -----------------------------------------------------------------------------
# Output 3: in-folder README documenting the cluster workflow.
# -----------------------------------------------------------------------------
readme_path <- file.path(OUT_DIR, "08_vep_input_full_README.txt")
readme <- c(
  "08_vep_input_full.* -- VEP-ready inputs (whole VAF)",
  "===================================================",
  "",
  sprintf("Source : %s", IN_FILE),
  sprintf("Variants: %d unique (chrom, pos, ref, alt) on GRCh38",
          nrow(per_var)),
  "Genome  : GRCh38 / hg38, UCSC chromosome naming (chr1..chr22, chrX, chrY)",
  "Scope   : whole cohort, panel-agnostic. AlphaMissense will only score",
  "          missense substitutions; non-missense rows will get a normal VEP",
  "          annotation but an empty AM score, which is expected.",
  "",
  "Files (this folder):",
  "  08_vep_input_full.vcf       VCFv4.2 (use with `vep -i ... .vcf`)",
  "  08_vep_input_full.ensembl   VEP default tab format (Part 3.1):",
  "                              chrom start end allele strand identifier",
  "                              (use with `vep -i ... .ensembl --format ensembl`)",
  "",
  "Joining AM scores back into the project",
  "---------------------------------------",
  "Once an annotated VCF is produced (vep + --plugin AlphaMissense), parse",
  "the CSQ INFO field, emit a TSV keyed by (chrom, pos, ref, alt) carrying",
  "am_pathogenicity / am_class / LoF / gene / consequence, and join onto",
  "results/ped_gof_snv/03_vaf_annotated.tsv (and any future panel TSV) by",
  "(chrom, pos, ref, alt) -- same key as this VCF.",
  "",
  "Notes",
  "-----",
  "  - chrom naming: this file uses chr1..chr22, chrX, chrY (UCSC-style)",
  "    to match a UCSC-style VEP cache. For an Ensembl-style cache,",
  "    strip the chr prefix or pass --synonyms to vep (chr_synonyms.txt path).",
  "  - For variants where the AM TSV is missing a row (synonymous, intronic,",
  "    indels longer than 1 nt), VEP writes the row with am_pathogenicity",
  "    blank -- expected, not an error."
)
writeLines(readme, readme_path)
cat("Wrote", readme_path, "\n")

cat("\nDone. Inputs ready in", OUT_DIR, "\n")
