# 02_inspect_vaf.R
# Deeper look at the VAF file:
#   - what tumor types are in PID prefixes?
#   - does anything in the SJ ID space line up with DepMap models?
#   - how many markers fall in canonical oncogene loci?

suppressPackageStartupMessages({
  library(data.table)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
DATA_DIR <- file.path(PROJ_DIR, "data")
# Per-run output folder. See R/01_profile_inputs.R for the layout rationale.
RUN_NAME <- "all_42genes"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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

# 4. Target oncogene loci (hg38) come from data/target_genes.tsv -- the
#    single source of truth for the rest of the pipeline. Built from the
#    SJPedPanel paper's GoF + pediatric-solid-tumor filter
#    (R/00_build_target_genes.R + data/oncogene_shortlist_sjpedpanel.tsv).
TARGETS_FILE <- file.path(DATA_DIR, "target_genes.tsv")
if (!file.exists(TARGETS_FILE)) {
  stop("Missing ", TARGETS_FILE,
       ". Run `Rscript R/00_build_target_genes.R` first.")
}
oncogene_loci <- fread(TARGETS_FILE,
                        select = c("gene","chrom","start","end","tier","score"))
cat(sprintf("Loaded %d target genes from %s\n",
            nrow(oncogene_loci), TARGETS_FILE))

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
cat("\n--- per-gene PID counts in those loci ---\n")
print(all_hits[, .(n_markers = .N, n_pids = uniqueN(PID)),
               by = gene][order(-n_pids)])

# Per-tier rollup so we can sanity-check that every tier has SOME signal
cat("\n--- per-TIER PID counts ---\n")
all_hits_tier <- merge(all_hits,
                       oncogene_loci[, .(gene, tier)],
                       by = "gene", all.x = TRUE)
print(all_hits_tier[, .(n_markers = .N, n_pids = uniqueN(PID),
                         n_genes_with_hits = uniqueN(gene)),
                     by = tier][order(tier)])

cat("\n--- first 30 oncogene-locus hits ---\n")
print(head(all_hits, 30))

fwrite(all_hits,
       file.path(OUT_DIR, "02_oncogene_locus_hits.csv"))
cat("\nWrote:", file.path(OUT_DIR, "02_oncogene_locus_hits.csv"), "\n")
