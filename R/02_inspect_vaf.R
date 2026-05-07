# 02_inspect_vaf.R
# Deeper look at the VAF file:
#   - what tumor types are in PID prefixes?
#   - does anything in the SJ ID space line up with DepMap models?
#   - how many markers fall in canonical oncogene loci?

suppressPackageStartupMessages({
  library(data.table)
})

DOC_DIR <- "C:/Users/mvijayan/Documents"
OUT_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes/results"

vaf  <- fread(file.path(DOC_DIR, "AllMarkers_VAF_long.tsv"))
meta <- fread(file.path(DOC_DIR, "depmap_meta.csv"))

# 1. PID prefix distribution -- the SJ prefix usually encodes tumor type.
cat("\n--- PID prefix distribution (first 5 letters) ---\n")
vaf[, prefix5 := substr(PID, 1, 5)]
print(vaf[, .(n_rows = .N, n_pids = uniqueN(PID)), by = prefix5][order(-n_pids)])

cat("\n--- PID prefix distribution (first 6 letters) ---\n")
vaf[, prefix6 := substr(PID, 1, 6)]
print(vaf[, .(n_rows = .N, n_pids = uniqueN(PID)), by = prefix6][order(-n_pids)])

# 2. Sample column inside VAF
cat("\n--- Sample column distribution ---\n")
print(vaf[, .N, by = Sample][order(-N)][1:20])

# 3. Look for cross-references inside meta. Is anything ever an SJxxx ID?
cat("\n--- looking for SJ-style IDs anywhere in meta ---\n")
char_cols <- names(meta)[sapply(meta, is.character)]
for (col in char_cols) {
  hits <- sum(grepl("^SJ", meta[[col]]))
  if (hits > 0) cat(sprintf("  meta$%s has %d rows starting with SJ\n", col, hits))
}

cat("\n--- ModelSubtypeFeatures or PatientSubtypeFeatures sample values ---\n")
print(head(unique(meta$ModelSubtypeFeatures), 5))
print(head(unique(meta$PatientSubtypeFeatures), 5))

# 4. A few canonical oncogene loci on hg38. We just want to know whether
#    *any* of the 10K markers land in them. This tells us whether the
#    file is going to be informative for the RAS-pathway analysis.
oncogene_loci <- data.table(
  gene  = c("KRAS",  "HRAS",  "NRAS",  "BRAF",   "PIK3CA", "EGFR",
            "ALK",   "MYCN",  "TP53",  "IDH1",   "FGFR1",  "MAP2K1",
            "ERBB2", "MET",   "MYC",   "CDKN2A", "RB1",    "PTEN",
            "CTNNB1","SMARCB1"),
  chrom = c("chr12","chr11","chr1", "chr7",   "chr3", "chr7",
            "chr2", "chr2", "chr17","chr2",   "chr8", "chr15",
            "chr17","chr7", "chr8", "chr9",   "chr13","chr10",
            "chr3", "chr22"),
  start = c(25205246, 533488,  114704464, 140719327, 179148115, 55019017,
            29192774, 15940550,7668402,  208236227, 38411138, 66386654,
            39688094, 116672196,127736233,21967751, 48303257, 87863438,
            41194741, 23786966),
  end   = c(25250929, 535567,  114716894, 140924929, 179240096, 55211628,
            29921586, 15947007,7687550,  208266074, 38468834, 66495020,
            39728660, 116798386,127741434,21995301, 48535063, 87971930,
            41260096, 23834540)
)

# Marker_hg38 looks like "chr12.25245347.C.T" -- one string per variant
# encoding chrom.pos.ref.alt with a literal period as the separator.
# tstrsplit() is the data.table vectorized version of strsplit; keep = 1:2
# extracts only the first two fields, returned as a list of column vectors.
# We drop SVs (structural variants) because they don't have a single
# position/ref/alt and can't be VEP'd the way SNVs/Indels can.
snv <- vaf[Marker.Type %in% c("SNV", "Indel", "MNV")]
snv[, c("chrom", "pos") := tstrsplit(Marker_hg38, ".", fixed = TRUE,
                                     keep = 1:2)]
snv[, pos := as.integer(pos)]

cat("\n--- markers per oncogene locus ---\n")
for (i in seq_len(nrow(oncogene_loci))) {
  g  <- oncogene_loci[i]
  hit <- snv[chrom == g$chrom & pos >= g$start & pos <= g$end]
  cat(sprintf("  %-7s %s:%d-%d  ->  %d markers, %d unique PIDs\n",
              g$gene, g$chrom, g$start, g$end,
              nrow(hit), uniqueN(hit$PID)))
}

# 5. Pull the markers that fall in any of those loci so we can eyeball them.
# The 'on = .(... pos >= start, pos <= end)' syntax is a data.table NON-EQUI
# JOIN: it joins each snv row to the oncogene_loci row whose [start, end]
# interval contains the snv's pos. nomatch = NULL drops snv rows that don't
# fall in any locus (i.e. inner join behaviour).
all_hits <- snv[oncogene_loci, on = .(chrom = chrom, pos >= start, pos <= end),
                nomatch = NULL,
                .(PID, Marker_hg38, Marker.Type, Sample, VAF, gene)]

cat(sprintf("\nTotal markers in those %d oncogenes: %d (across %d PIDs)\n",
            nrow(oncogene_loci), nrow(all_hits), uniqueN(all_hits$PID)))
cat("\n--- per-gene PID counts in those 20 loci ---\n")
print(all_hits[, .(n_markers = .N, n_pids = uniqueN(PID)), by = gene][order(-n_pids)])

cat("\n--- first 30 oncogene-locus hits ---\n")
print(head(all_hits, 30))

fwrite(all_hits,
       file.path(OUT_DIR, "02_oncogene_locus_hits.csv"))
cat("\nWrote:", file.path(OUT_DIR, "02_oncogene_locus_hits.csv"), "\n")
