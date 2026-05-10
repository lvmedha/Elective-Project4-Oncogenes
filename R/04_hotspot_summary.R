# 04_hotspot_summary.R
# Inspect the annotated VAF: which protein changes are most recurrent in
# each target gene? Compare to the curated activating-hotspot list
# (data/hotspots_tier0.csv) to summarize Tier-0 vs other-missense calls
# before 05_apply_hotspot_tiers.R does the formal assignment.

suppressPackageStartupMessages({
  library(data.table)
})

PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
DATA_DIR <- file.path(PROJ_DIR, "data")
RUN_NAME <- "ped_gof_snv"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
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

# Check how many missense changes hit a known activating residue.
# Source the (gene, codon) hotspot list from data/hotspots_tier0.csv so it
# stays in sync with what 05_apply_hotspot_tiers.R uses for formal
# Mut_Status assignment. Genes in the target set without curated entries
# (e.g. fusion partners, amplification-driven oncogenes) just get FALSE for
# every variant -- they'll all flow into Missense_Other in step 05.
hotspots <- fread(file.path(DATA_DIR, "hotspots_tier0.csv"))
hotspots[, codon := as.integer(codon)]
cat("Tier-0 hotspot file: ", nrow(hotspots),
    " (gene, codon) entries across ",
    uniqueN(hotspots$gene), " genes\n", sep = "")

ann[, codon := protein_start]
# Two-key lookup ("is this gene+codon in the hotspot file?"): flatten both
# keys into a single string ("KRAS 12") and use %in% set membership --
# faster and shorter than a join.
ann[, hotspot_codon := paste(gene_symbol, codon) %in%
                       paste(hotspots$gene, hotspots$codon)]

cat("\nKnown-codon hits per gene:\n")
print(ann[Variant_Class == "Missense",
          .(n_unique_HGVSp = uniqueN(HGVSp_Short),
            n_unique_HGVSp_at_hotspot = uniqueN(HGVSp_Short[hotspot_codon]),
            n_pids = uniqueN(PID),
            n_pids_at_hotspot = uniqueN(PID[hotspot_codon])),
          by = gene_symbol][order(-n_pids_at_hotspot)])

# Also: tumor-type mix among hotspot-hit PIDs
ann[, prefix6 := substr(PID, 1, 6)]
cat("\nTumor-type mix in PIDs carrying a known hotspot in any target gene:\n")
print(ann[Variant_Class == "Missense" & hotspot_codon == TRUE,
          .(n_pids = uniqueN(PID)), by = prefix6][order(-n_pids)])
