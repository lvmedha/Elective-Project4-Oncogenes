# 05_apply_hotspot_tiers.R
# Take the annotated VAF (03) and the curated Tier-0 hotspot list, and emit
# a per-PID-x-gene mutation table with a Mut_Status column.
#
# Tiers used here (we'll add AlphaMissense later for Tier 1/2/3):
#   - Hotspot       : missense at a canonical activating codon (Tier-0 file)
#   - Missense_Other: missense not at a Tier-0 codon
#   - InframeIndel  : inframe insertion / deletion / protein_altering
#   - Truncating    : nonsense / frameshift / splice / start- or stop-lost
#   - Silent        : synonymous
#   - Other         : UTR, intron, flanking, non-coding (we drop these)

suppressPackageStartupMessages({
  library(data.table)
})

PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
RUN_NAME <- "ped_gof_snv"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ann      <- fread(file.path(OUT_DIR, "03_vaf_annotated.tsv"))
hotspots <- fread(file.path(PROJ_DIR, "data/hotspots_tier0.csv"))

# Drop non-coding consequences entirely
keep_classes <- c("Missense", "Nonsense", "Frameshift", "InframeDel",
                  "InframeIns", "InframeOther", "Splice", "StartLost",
                  "StopLost", "Silent")
muts <- ann[Variant_Class %in% keep_classes]
cat("Coding/exon variants kept:", nrow(muts), "across",
    uniqueN(muts$PID), "PIDs\n")

# Tag hotspot codon membership.
# We need a TWO-key lookup ("is this gene+codon in the hotspot file?"). The
# trick below flattens both keys into a single string ("KRAS 12") and uses
# %in% set membership -- much faster and shorter than a join for this size.
hotspots[, codon := as.integer(codon)]
muts[, codon := as.integer(protein_start)]
muts[, hotspot_codon := paste(gene_symbol, codon) %in%
                       paste(hotspots$gene, hotspots$codon)]

# fcase() is data.table's vectorised switch/case (like SQL CASE WHEN).
# Each (condition, value) pair is checked in order and the FIRST match wins,
# so the order of branches matters. Unmatched rows fall through to `default`.
muts[, Mut_Status := fcase(
  Variant_Class == "Missense"     & hotspot_codon  , "Hotspot",
  Variant_Class == "Missense"     & !hotspot_codon , "Missense_Other",
  Variant_Class %in% c("InframeDel","InframeIns","InframeOther"),
                                                    "InframeIndel",
  Variant_Class %in% c("Nonsense","Frameshift","Splice","StartLost",
                       "StopLost"),                  "Truncating",
  Variant_Class == "Silent",                         "Silent",
  default = "Other"
)]

cat("\nMut_Status x gene distribution:\n")
print(dcast(muts, gene_symbol ~ Mut_Status,
            value.var = "PID", fun.aggregate = uniqueN, fill = 0L))

# Per-gene per-PID collapse. A patient can have several variants in the
# same gene (e.g. a hotspot AND a passenger). For the addiction analysis we
# want one row per (PID, gene) representing the strongest event:
# Hotspot > Missense_Other > InframeIndel > Truncating > Silent.
#
# Trick: build a NAMED VECTOR mapping status -> rank, then index into it
# with the column to get a numeric rank column without writing a loop or
# a long if/else. setorder() then puts the strongest event first per
# (PID, gene), and unique(by = ...) keeps the FIRST row of each group.
status_priority <- c(Hotspot=5, Missense_Other=4, InframeIndel=3,
                     Truncating=2, Silent=1, Other=0)
muts[, status_rank := status_priority[Mut_Status]]
setorder(muts, PID, gene_symbol, -status_rank, -VAF)
mut_per_pid_gene <- unique(muts, by = c("PID", "gene_symbol"))

cat("\nPer-PID x per-gene unique events:", nrow(mut_per_pid_gene), "\n")

cat("\nTop few rows of the cleaned table:\n")
print(head(mut_per_pid_gene[, .(PID, gene_symbol, HGVSp_Short, Variant_Class,
                                Mut_Status, VAF)], 10))

cat("\nN PIDs per (gene x Mut_Status):\n")
print(mut_per_pid_gene[, .N, by = .(gene_symbol, Mut_Status)
                      ][order(gene_symbol, -N)])

fwrite(mut_per_pid_gene, file.path(OUT_DIR, "05_mutations_tiered.tsv"),
       sep = "\t")
fwrite(muts,             file.path(OUT_DIR, "05_mutations_tiered_long.tsv"),
       sep = "\t")
cat("\nWrote:\n",
    " - ", file.path(OUT_DIR, "05_mutations_tiered.tsv"),
    "  (1 row per PID x gene)\n",
    " - ", file.path(OUT_DIR, "05_mutations_tiered_long.tsv"),
    "  (1 row per variant)\n", sep = "")
