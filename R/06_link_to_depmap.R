# 06_link_to_depmap.R
# Bring it all together:
#   - load the SJ-PID -> DepMap ACH mapping (samples.txt)
#   - attach ACH ModelIDs to the tiered mutation table from 05
#   - check how many cell lines we end up with per gene x Mut_Status
#   - check coverage against the CGE row IDs

suppressPackageStartupMessages({
  library(data.table)
})

DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
OUT_DIR  <- file.path(PROJ_DIR, "results")

# ---- 1. Load samples.txt --------------------------------------------------
samp <- fread(file.path(DOC_DIR, "samples.txt"),
              header = FALSE,
              col.names = c("ACH", "SJ_full", "SJ_clone", "CellLineName"))
cat("samples.txt rows:", nrow(samp), "\n")
cat("unique ACH IDs:", uniqueN(samp$ACH), "\n")
cat("unique SJ_clone keys:", uniqueN(samp$SJ_clone), "\n\n")

# samples.txt's third column looks like "SJRHB098835_C1". The piece before
# "_C<N>" is the SJ patient ID; "C<N>" is the clone (also encoded in VAF's
# `Sample` field as e.g. "C1_WGS"). We split them with two regexes so the
# join key in step 4 is (PID, Sample) and matches both files.
samp[, PID    := sub("_C[0-9]+$", "", SJ_clone)]
samp[, Sample := sub("^.*_(C[0-9]+)$", "\\1", SJ_clone)]
cat("First 5 rows of decomposed mapping:\n")
print(head(samp[, .(PID, Sample, ACH, CellLineName)], 5))

# ---- 2. Load tiered mutations + CGE row ids -------------------------------
mut <- fread(file.path(OUT_DIR, "05_mutations_tiered.tsv"))
cat("\nMutation table rows:", nrow(mut),
    " (", uniqueN(mut$PID), "PIDs )\n")

# VAF stores Sample as e.g. "C1_WGS" but samples.txt stores it as "C1".
# Strip the trailing assay suffix so the two column values can join.
mut[, sample_clone := sub("_(WGS|RNA|WES)$", "", Sample)]

# Two-key join (PID, sample_clone). all.x = TRUE keeps every mutation row,
# even those that don't have a cell-line equivalent (those rows get NA in
# the ACH column and are filtered/flagged downstream).
mut2 <- merge(mut,
              samp[, .(PID, sample_clone = Sample, ACH, CellLineName)],
              by = c("PID", "sample_clone"),
              all.x = TRUE)

cat("\nPID+sample_clone keys that mapped to an ACH:",
    sum(!is.na(mut2$ACH)),
    "/", nrow(mut2), "rows\n")
cat("Distinct mutated cell lines (ACH):",
    uniqueN(mut2$ACH[!is.na(mut2$ACH)]), "\n")

# ---- 3. Cross-check against the CGE row IDs -------------------------------
cge_first_col <- fread(file.path(DOC_DIR, "depmap_CGE.csv"),
                       select = 1, col.names = "ModelID")
cat("\nCGE rows:", nrow(cge_first_col), "\n")
cat("Mutated ACH that also have CGE scores:",
    sum(unique(mut2$ACH) %in% cge_first_col$ModelID,
        na.rm = TRUE), "\n")

# ---- 4. Per-gene cell-line counts now that we have ACH IDs ----------------
have_cge <- cge_first_col$ModelID
mut2[, in_CGE := !is.na(ACH) & ACH %in% have_cge]

cat("\nPer-gene mutated CELL LINE counts (i.e. with both mutation AND CGE):\n")
print(mut2[in_CGE == TRUE, .N, by = .(gene_symbol, Mut_Status)
           ][order(gene_symbol, -N)])

cat("\nPer-gene Hotspot-only cell-line counts:\n")
print(mut2[in_CGE & Mut_Status == "Hotspot",
           .(n_lines = uniqueN(ACH)), by = gene_symbol][order(-n_lines)])

# ---- 5. How many of those mutated cell lines are pediatric solid? --------
meta <- fread(file.path(DOC_DIR, "depmap_meta.csv"))
ped_solid_lineages <- c("Bone", "Peripheral Nervous System", "Soft Tissue",
                        "CNS/Brain", "Kidney", "Liver", "Eye")
ped_solid <- meta[AgeCategory == "Pediatric" &
                  OncotreeLineage %in% ped_solid_lineages, ModelID]
cat("\nPediatric-solid-tumor cell lines (in meta):", length(ped_solid), "\n")
cat("...of which have CGE:", sum(ped_solid %in% have_cge), "\n")
cat("...of which carry a mutation in any of our 10 genes:",
    uniqueN(mut2[ACH %in% ped_solid]$ACH), "\n")

cat("\nPer-gene mutated cell-line counts WITHIN pediatric solid tumors:\n")
print(mut2[in_CGE & ACH %in% ped_solid,
           .(n_lines = uniqueN(ACH)), by = .(gene_symbol, Mut_Status)
           ][order(gene_symbol, -n_lines)])

# ---- 6. Save the joined per-cell-line mutation table ----------------------
fwrite(mut2, file.path(OUT_DIR, "06_mutations_with_ACH.tsv"), sep = "\t")
cat("\nWrote:", file.path(OUT_DIR, "06_mutations_with_ACH.tsv"), "\n")
