# 04_hotspot_summary.R
# Inspect the annotated VAF: which protein changes are most recurrent in
# each of the 10 target genes? Then (manually) compare to the canonical
# activating-hotspot list to build the Tier-0 set.

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes/results"
ann <- fread(file.path(OUT_DIR, "03_vaf_annotated.tsv"))

cat("Annotated rows:", nrow(ann), "\n")
cat("Unique PIDs in annotated table:", uniqueN(ann$PID), "\n\n")

# Top recurrent HGVSp per gene (missense only)
mis <- ann[Variant_Class == "Missense" & !is.na(HGVSp_Short)]
top <- mis[, .(n_pids = uniqueN(PID),
               n_rows = .N,
               example_PIDs = paste(head(unique(PID), 3), collapse=";")),
           by = .(gene_symbol, HGVSp_Short)]
setorder(top, gene_symbol, -n_pids)

cat("Top recurrent missense changes per gene:\n")
for (g in unique(top$gene_symbol)) {
  cat("\n=== ", g, " ===\n", sep = "")
  print(top[gene_symbol == g][1:min(.N, 10)])
}

fwrite(top, file.path(OUT_DIR, "04_recurrent_missense_per_gene.tsv"),
       sep = "\t")
cat("\nWrote:", file.path(OUT_DIR, "04_recurrent_missense_per_gene.tsv"), "\n")

# Check how many missense changes hit a known activating residue. We use
# a coarse, well-known list. (Full Tier-0 curation goes in next step.)
known_codons <- list(
  KRAS   = c(12, 13, 61, 117, 146),
  NRAS   = c(12, 13, 61),
  HRAS   = c(12, 13, 61),
  BRAF   = c(469, 581, 594, 600, 597, 601),
  PIK3CA = c(88, 110, 542, 545, 546, 1043, 1047),
  MAP2K1 = c(53, 56, 57, 67, 124, 211),
  EGFR   = c(719, 746, 747, 858, 861),
  ERBB2  = c(310, 755, 777, 842, 869),
  ALK    = c(1174, 1196, 1245, 1275),
  CTNNB1 = c(32, 33, 34, 35, 37, 41, 45)
)

ann[, codon := protein_start]
# For each row, check whether its (gene, codon) is a known activating
# residue. mapply walks the two columns in lockstep; the closure looks up
# the gene's codon vector in `known_codons` and tests membership.
# (Used here for the QC summary; the production tier assignment in 05
# uses a faster vectorised paste()-trick on the curated CSV.)
ann[, hotspot_codon := mapply(function(g, c) {
  if (is.na(g) || is.na(c)) return(FALSE)
  c %in% known_codons[[g]]
}, gene_symbol, codon)]

cat("\nKnown-codon hits per gene:\n")
print(ann[Variant_Class == "Missense",
          .(n_unique_HGVSp = uniqueN(HGVSp_Short),
            n_unique_HGVSp_at_hotspot = uniqueN(HGVSp_Short[hotspot_codon]),
            n_pids = uniqueN(PID),
            n_pids_at_hotspot = uniqueN(PID[hotspot_codon])),
          by = gene_symbol][order(-n_pids_at_hotspot)])

# Also: tumor-type mix among hotspot-hit PIDs
ann[, prefix6 := substr(PID, 1, 6)]
cat("\nTumor-type mix in PIDs carrying a known hotspot in any of the 10 genes:\n")
print(ann[Variant_Class == "Missense" & hotspot_codon == TRUE,
          .(n_pids = uniqueN(PID)), by = prefix6][order(-n_pids)])
