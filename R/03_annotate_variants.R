# 03_annotate_variants.R
# "Fix the format" of AllMarkers_VAF_long.tsv:
#
#   - keep only SNV/Indel/MNV markers that fall in any of our 10 target
#     oncogene loci
#   - submit them to Ensembl VEP REST (GRCh38) to get gene_symbol,
#     consequence, HGVSp, protein position, amino-acid change
#   - keep canonical-transcript annotations only
#   - write a tidy annotated TSV to results/03_vaf_annotated.tsv
#
# NB: "AllMarkers_VAF_long.tsv" is keyed by St. Jude tumor PIDs, not by
# DepMap ModelIDs. This script only fixes the *variant-level* format.
# Patient-to-cell-line linkage is a separate problem (see 99_TODO.md).

suppressPackageStartupMessages({
  library(data.table)
  library(httr2)
  library(jsonlite)
})

DOC_DIR <- "C:/Users/mvijayan/Documents"
OUT_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes/results"
CACHE_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes/data/cache"
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. The 10 target oncogenes and their canonical hg38 loci (with 2-kb pad).
#    Coordinates verified against UCSC/Ensembl GRCh38 gene records.
# -----------------------------------------------------------------------------
target_genes <- data.table(
  gene  = c("KRAS",  "NRAS",  "HRAS",  "BRAF",   "PIK3CA",
            "MAP2K1","EGFR",  "ERBB2", "ALK",    "CTNNB1"),
  chrom = c("chr12", "chr1",  "chr11", "chr7",   "chr3",
            "chr15", "chr7",  "chr17", "chr2",   "chr3"),
  start = c(25205246,114704464,533488,140719327,179148115,
            66386654,55019017, 39688094,29192774,41194741),
  end   = c(25250929,114716894,535567,140924929,179240096,
            66495020,55211628, 39728660,29921586,41260096)
)
# Pad each gene's locus by 2 kb on each side so we don't miss splice-region
# variants and 5'/3'-UTR mutations that sit just outside the CDS.
PAD <- 2000L
target_genes[, c("start", "end") := .(start - PAD, end + PAD)]
setkey(target_genes, chrom, start, end)

# -----------------------------------------------------------------------------
# 2. Load VAF, parse Marker_hg38, restrict to SNV/Indel/MNV in target loci.
# -----------------------------------------------------------------------------
vaf <- fread(file.path(DOC_DIR, "AllMarkers_VAF_long.tsv"))
vaf <- vaf[Marker.Type %in% c("SNV", "Indel", "MNV")]

vaf[, c("chrom", "pos", "ref", "alt") :=
      tstrsplit(Marker_hg38, ".", fixed = TRUE, keep = 1:4)]
vaf[, pos := as.integer(pos)]

cat("VAF rows after restricting to SNV/Indel/MNV:", nrow(vaf), "\n")

# foverlaps() does an interval-overlap join. It expects BOTH tables to have
# (chrom, start, end) -- target_genes already does, but a VAF row is a
# single point, so we make a degenerate 1-bp window where start == end == pos.
# type = "within" means we keep VAF rows whose [start,end] is fully INSIDE a
# target_genes [start,end]. nomatch = NULL drops non-overlapping rows.
vaf[, start := pos]
vaf[, end   := pos]
setkey(vaf, chrom, start, end)

vaf_in_targets <- foverlaps(vaf, target_genes,
                            by.x = c("chrom", "start", "end"),
                            type = "within", nomatch = NULL)
vaf_in_targets[, gene_target := gene]   # gene = column from target_genes
setnames(vaf_in_targets, c("i.start", "i.end"), c("vaf_start", "vaf_end"),
         skip_absent = TRUE)

cat("Markers in target-gene loci:", nrow(vaf_in_targets), "\n")
cat("Unique variants:", uniqueN(vaf_in_targets$Marker_hg38), "\n")

# Per-gene tally before VEP (informative)
cat("\nPer-target-gene marker counts (pre-VEP, locus-based):\n")
print(vaf_in_targets[, .(n_markers = .N,
                         n_unique  = uniqueN(Marker_hg38),
                         n_pids    = uniqueN(PID)), by = gene_target])

# -----------------------------------------------------------------------------
# 3. Build the unique VEP-input list. VEP expects 'CHROM POS . REF ALT . . .'
#    with chromosome names WITHOUT the 'chr' prefix.
# -----------------------------------------------------------------------------
# We dedup before sending: 100s of patients may carry the same KRAS G12D, but
# VEP's answer for that variant is the same every time. Dedup -> ~10x fewer
# API calls -> faster + nicer to Ensembl's servers.
unique_variants <- unique(
  vaf_in_targets[, .(chrom, pos, ref, alt, Marker_hg38)]
)
# Ensembl uses GRCh38 chromosome names WITHOUT the 'chr' prefix
# (e.g. "12" not "chr12"), so strip it before building the VEP input string.
# The format VEP expects is the standard VCF body: CHROM POS ID REF ALT QUAL FILTER INFO
# (the 4 trailing dots are placeholder ID/QUAL/FILTER/INFO fields).
unique_variants[, vep_chrom := sub("^chr", "", chrom)]
unique_variants[, vep_input := sprintf("%s %d . %s %s . . .",
                                       vep_chrom, pos, ref, alt)]

cat("\nUnique variants to send to VEP:", nrow(unique_variants), "\n")

# -----------------------------------------------------------------------------
# 4. Call Ensembl VEP REST in batches of 200.
#    Endpoint: POST https://rest.ensembl.org/vep/human/region
#    Cached to data/cache/vep_<hash>.json so re-runs are free.
# -----------------------------------------------------------------------------
vep_endpoint <- "https://rest.ensembl.org/vep/human/region"

vep_query <- function(variants_chunk) {
  # The flags request optional VEP outputs: canonical-transcript flag,
  # full HGVS strings, protein-level annotation, exon/intron numbers,
  # and MANE Select transcript identifiers.
  body <- list(
    variants    = variants_chunk,
    canonical   = 1,
    hgvs        = 1,
    protein     = 1,
    numbers     = 1,
    mane        = 1
  )
  body_json <- toJSON(body, auto_unbox = TRUE)
  # Cache by hash of the request body. If the same chunk is queried again
  # (e.g. when re-running this script), we read the response from disk
  # instead of hitting the API. ~10s vs ~3s on a warm cache.
  cache_key <- digest::digest(body_json, algo = "md5")
  cache_file <- file.path(CACHE_DIR, paste0("vep_", cache_key, ".json"))
  if (file.exists(cache_file)) return(fromJSON(cache_file,
                                               simplifyVector = FALSE))
  # POST the JSON body to the REST endpoint; retry up to 4 times with a 5s
  # backoff to absorb transient network or rate-limit failures.
  resp <- request(vep_endpoint) |>
    req_method("POST") |>
    req_headers(`Content-Type` = "application/json",
                Accept = "application/json") |>
    req_body_raw(body_json) |>
    req_retry(max_tries = 4, backoff = ~ 5) |>
    req_timeout(120) |>
    req_perform()
  result <- resp_body_json(resp)
  writeLines(toJSON(result, auto_unbox = TRUE, pretty = FALSE),
             cache_file)
  result
}

# digest is a tiny package -- install if missing
if (!requireNamespace("digest", quietly = TRUE)) {
  install.packages("digest", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(digest))

# Ensembl REST allows up to 200 variants per /vep/region POST. We split into
# chunks of 200 just to stay polite and resilient (smaller batches mean
# smaller retries on failure). With 172 variants we end up with one batch.
batch_size <- 200
batches <- split(unique_variants$vep_input,
                 ceiling(seq_along(unique_variants$vep_input) / batch_size))
cat("\nQuerying VEP in", length(batches), "batches of up to",
    batch_size, "variants each ...\n")

all_responses <- vector("list", length(batches))
for (i in seq_along(batches)) {
  cat("  batch", i, "/", length(batches), "  (n =",
      length(batches[[i]]), ") ... ")
  all_responses[[i]] <- vep_query(batches[[i]])
  cat("ok\n")
}

# -----------------------------------------------------------------------------
# 5. Flatten the JSON. VEP returns a list-per-variant; each variant has
#    transcript_consequences -> list. We pull the canonical row (or MANE Select)
#    per variant.
# -----------------------------------------------------------------------------
# Helper: VEP JSON fields are missing whenever they don't apply (e.g.
# protein_start is absent for intronic variants). na_or() returns a safe
# typed default in those cases so we can build a tidy data.table without
# blowing up on NULLs.
na_or <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0) return(default)
  x[[1]]
}

# VEP's response is a list-of-variants; each variant has 0..N transcript
# consequences (one per overlapping transcript). We emit ONE data.table row
# per (variant x transcript) so that downstream we can pick the right
# transcript per gene (see the dedup step further down).
flatten_vep <- function(per_variant) {
  out <- list()
  for (v in per_variant) {
    tcs <- v$transcript_consequences
    if (is.null(tcs) || length(tcs) == 0) next
    for (t in tcs) {
      out[[length(out) + 1]] <- data.table(
        vep_input     = v$input,
        most_severe   = na_or(v$most_severe_consequence, NA_character_),
        gene_symbol   = na_or(t$gene_symbol,   NA_character_),
        gene_id       = na_or(t$gene_id,       NA_character_),
        transcript_id = na_or(t$transcript_id, NA_character_),
        mane_select   = na_or(t$mane_select,   NA_character_),
        canonical     = as.integer(na_or(t$canonical, 0L)),
        biotype       = na_or(t$biotype,       NA_character_),
        consequence   = paste(unlist(t$consequence_terms), collapse = ","),
        protein_start = as.integer(na_or(t$protein_start, NA_integer_)),
        protein_end   = as.integer(na_or(t$protein_end,   NA_integer_)),
        amino_acids   = na_or(t$amino_acids,   NA_character_),
        codons        = na_or(t$codons,        NA_character_),
        hgvsc         = na_or(t$hgvsc,         NA_character_),
        hgvsp         = na_or(t$hgvsp,         NA_character_)
      )
    }
  }
  rbindlist(out, fill = TRUE)
}

vep_dt <- rbindlist(lapply(all_responses, flatten_vep), fill = TRUE)
cat("\nVEP returned", nrow(vep_dt),
    "transcript consequences across",
    uniqueN(vep_dt$vep_input), "variants.\n")

# Restrict to our 10 target genes only. NRAS-locus variants will produce
# rows for both NRAS and the overlapping CSDE1 gene; this filter drops CSDE1.
vep_dt <- vep_dt[gene_symbol %in% target_genes$gene]
# We still might have multiple rows per (variant, gene): one per Ensembl
# transcript (canonical, RefSeq, alt isoforms). We want ONE row representing
# the standard transcript for each gene, with this preference order:
#   1. MANE Select (= NCBI/Ensembl jointly-curated reference transcript)
#   2. Ensembl canonical
#   3. anything else
# setorder with negative columns sorts descending, so MANE-tagged rows
# (has_mane = 1) come before non-MANE (0); then unique() keeps the FIRST
# row per group => effectively "best transcript" per gene.
vep_dt[, has_mane := as.integer(nzchar(mane_select))]
setorder(vep_dt, vep_input, gene_symbol, -has_mane, -canonical)
vep_dt <- unique(vep_dt, by = c("vep_input", "gene_symbol"))
vep_dt[, has_mane := NULL]
cat("After dedup (1 row per variant x target gene):",
    nrow(vep_dt), "rows across",
    uniqueN(vep_dt$vep_input), "variants.\n")

# -----------------------------------------------------------------------------
# 6. Build a compact HGVSp_Short ('p.G12D' style) from amino_acids + protein_start
#    and parse the consequence into a coarse Variant_Class.
# -----------------------------------------------------------------------------
aa3to1 <- c(Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Glu="E", Gln="Q",
            Gly="G", His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F",
            Pro="P", Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V", Sec="U",
            Pyl="O", Ter="*", Stop="*")

# VEP encodes the protein change as "REF/ALT" in the amino_acids field
# (e.g. "G/D" for KRAS G12D). Split on "/" to get separate REF and ALT AAs,
# then build a compact HGVSp_Short like "p.G12D" using protein_start as the
# residue number. Variants without a protein change (intronic, UTR) get NA.
vep_dt[, c("aa_ref", "aa_alt") := tstrsplit(amino_acids, "/", fixed = TRUE)]
vep_dt[, HGVSp_Short := fifelse(
  !is.na(aa_ref) & !is.na(aa_alt) & !is.na(protein_start),
  sprintf("p.%s%d%s",
          aa_ref,
          protein_start,
          aa_alt),
  NA_character_
)]

# VEP returns a comma-separated list of consequence terms for each variant
# (e.g. "missense_variant,splice_region_variant"). We collapse this multi-
# label output into a single coarse class using the priority below: the
# FIRST matching branch wins, so e.g. a missense+splice variant becomes
# "Missense". The order encodes biological severity / interpretability.
classify_consequence <- function(cs) {
  if (is.na(cs) || cs == "") return(NA_character_)
  csv <- strsplit(cs, ",", fixed = TRUE)[[1]]
  if ("missense_variant"             %in% csv) return("Missense")
  if ("stop_gained"                  %in% csv) return("Nonsense")
  if ("frameshift_variant"           %in% csv) return("Frameshift")
  if ("stop_lost"                    %in% csv) return("StopLost")
  if ("start_lost"                   %in% csv) return("StartLost")
  if ("inframe_insertion"            %in% csv) return("InframeIns")
  if ("inframe_deletion"             %in% csv) return("InframeDel")
  if ("splice_acceptor_variant"      %in% csv ||
      "splice_donor_variant"         %in% csv) return("Splice")
  if ("synonymous_variant"           %in% csv) return("Silent")
  if ("protein_altering_variant"     %in% csv) return("InframeOther")
  if ("5_prime_UTR_variant"          %in% csv ||
      "3_prime_UTR_variant"          %in% csv) return("UTR")
  if ("intron_variant"               %in% csv) return("Intron")
  if ("upstream_gene_variant"        %in% csv ||
      "downstream_gene_variant"      %in% csv) return("Flanking")
  csv[1]
}
vep_dt[, Variant_Class := vapply(consequence, classify_consequence, "")]

# -----------------------------------------------------------------------------
# 7. Join VEP back to the per-PID VAF rows; keep only target genes.
# -----------------------------------------------------------------------------
# Two-step join:
#   ann   = unique variants enriched with their VEP annotation
#   final = expand back to one row per (PID, variant) -- the per-patient grain.
# allow.cartesian = TRUE because variants in overlapping genes may produce
# multiple matching rows; the gene_symbol == gene_target filter below
# resolves them to one annotation per locus.
ann <- merge(unique_variants, vep_dt, by = "vep_input",
             all.x = FALSE, allow.cartesian = TRUE)

final <- merge(vaf_in_targets[, .(PID, Sample, Marker_hg38, Marker.Type,
                                  nM, nT, VAF, gene_target)],
               ann[, .(Marker_hg38, gene_symbol, transcript_id, mane_select,
                       canonical, biotype, consequence, Variant_Class,
                       protein_start, amino_acids, HGVSp_Short, hgvsc, hgvsp)],
               by = "Marker_hg38", allow.cartesian = TRUE)

# Keep only rows where VEP's gene_symbol matches the gene we got by locus
# overlap. This removes e.g. CSDE1-annotated rows that come along when an
# NRAS variant happens to fall inside the CSDE1 gene body too.
final <- final[gene_symbol == gene_target]
setcolorder(final, c("PID","Sample","gene_target","gene_symbol",
                     "Marker_hg38","Marker.Type","Variant_Class",
                     "HGVSp_Short","amino_acids","protein_start",
                     "consequence","mane_select","transcript_id","canonical",
                     "VAF","nM","nT","hgvsc","hgvsp","biotype"))

# -----------------------------------------------------------------------------
# 8. Per-gene QC summary
# -----------------------------------------------------------------------------
cat("\nPer-target-gene Variant_Class summary (post-VEP):\n")
print(final[!is.na(gene_symbol),
            .N, keyby = .(gene_symbol, Variant_Class)])

cat("\nMissense-only counts per gene (n_pids):\n")
print(final[Variant_Class == "Missense",
            .(n_variants = .N, n_unique = uniqueN(HGVSp_Short),
              n_pids = uniqueN(PID)), by = gene_symbol][order(-n_pids)])

# -----------------------------------------------------------------------------
# 9. Write outputs
# -----------------------------------------------------------------------------
fwrite(final,    file.path(OUT_DIR, "03_vaf_annotated.tsv"), sep = "\t")
fwrite(vep_dt,   file.path(OUT_DIR, "03_vep_per_variant.tsv"), sep = "\t")
cat("\nWrote:\n",
    " - ", file.path(OUT_DIR, "03_vaf_annotated.tsv"), "\n",
    " - ", file.path(OUT_DIR, "03_vep_per_variant.tsv"), "\n", sep = "")
