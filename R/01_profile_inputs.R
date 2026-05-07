# 01_profile_inputs.R
# Quick profile of the three input datasets so we know exactly what we have
# before we start joining anything.
#
# Inputs (from Documents/):
#   - depmap_meta.csv          (cell-line metadata, keyed on ModelID)
#   - depmap_CGE.csv           (CRISPR Chronos gene effect, wide format)
#   - AllMarkers_VAF_long.tsv  (mutation calls, long format)

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE))
    install.packages("data.table", repos = "https://cloud.r-project.org")
  library(data.table)
})

# OUT_DIR is per-run so we can keep yesterday's pilot outputs frozen under
# results/pilot_10genes/ while regenerating the 42-gene run under
# results/all_42genes/. Change RUN_NAME to point a fresh re-run at a new
# subfolder.
DOC_DIR  <- "C:/Users/mvijayan/Documents"
PROJ_DIR <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
RUN_NAME <- "all_42genes"
OUT_DIR  <- file.path(PROJ_DIR, "results", RUN_NAME)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n================ depmap_meta.csv ================\n")
meta <- fread(file.path(DOC_DIR, "depmap_meta.csv"))
cat("dim:", dim(meta), "\n")
cat("ncol:", ncol(meta), "\n\n")
cat("colnames:\n"); print(colnames(meta))
cat("\nfirst 3 ModelIDs:\n"); print(head(meta$ModelID, 3))

cat("\n--- AgeCategory counts ---\n")
print(meta[, .N, by = AgeCategory][order(-N)])

cat("\n--- OncotreeLineage counts (top 25) ---\n")
print(meta[, .N, by = OncotreeLineage][order(-N)][1:25])

cat("\n--- AgeCategory x OncotreeLineage (Pediatric only, top 25) ---\n")
print(meta[AgeCategory == "Pediatric", .N, by = OncotreeLineage][order(-N)][1:25])

cat("\n================ AllMarkers_VAF_long.tsv ================\n")
vaf <- fread(file.path(DOC_DIR, "AllMarkers_VAF_long.tsv"))
cat("dim:", dim(vaf), "\n")
cat("colnames:\n"); print(colnames(vaf))
cat("\nfirst 5 rows:\n"); print(head(vaf, 5))

cat("\n--- unique PIDs ---\n")
cat("n unique PID:", length(unique(vaf$PID)), "\n")
cat("first 10 PIDs:\n"); print(head(unique(vaf$PID), 10))
cat("first 10 Samples:\n"); print(head(unique(vaf$Sample), 10))

cat("\n--- Marker.Type counts ---\n")
print(vaf[, .N, by = Marker.Type])

# Does PID match any DepMap ModelID directly?
cat("\n--- PID vs DepMap ModelID overlap ---\n")
cat("vaf PIDs that match meta$ModelID:",
    sum(unique(vaf$PID) %in% meta$ModelID), "\n")
cat("vaf PIDs that match meta$StrippedCellLineName:",
    sum(unique(vaf$PID) %in% meta$StrippedCellLineName), "\n")
cat("vaf PIDs that match meta$CCLEName:",
    sum(unique(vaf$PID) %in% meta$CCLEName), "\n")
cat("vaf PIDs that match meta$SangerModelID:",
    sum(unique(vaf$PID) %in% meta$SangerModelID), "\n")

# Try a fuzzy partial match too
cat("\nany PIDs that look like ACH-XXXXXX?:",
    sum(grepl("^ACH-", unique(vaf$PID))), "\n")

cat("\n================ depmap_CGE.csv ================\n")
# CGE is 440 MB and ~18 500 columns wide. We never need all of it loaded
# at once during profiling, so we read in two cheap passes:
#   pass 1: nrows = 0 -> just the header line, gives us the column names
cge_header <- fread(file.path(DOC_DIR, "depmap_CGE.csv"), nrows = 0)
cat("ncol:", ncol(cge_header), "\n")
cat("first 5 columns:\n"); print(head(colnames(cge_header), 5))
cat("last 5 columns:\n"); print(tail(colnames(cge_header), 5))

#   pass 2: select = first column only, nrows = 5 -> first 5 row IDs
first_col <- fread(file.path(DOC_DIR, "depmap_CGE.csv"),
                   select = colnames(cge_header)[1], nrows = 5)
cat("\nfirst 5 rows of column 1 (", colnames(cge_header)[1], "):\n", sep = "")
print(first_col)

# Same trick to count rows without ever loading the gene matrix.
cge_rows <- fread(file.path(DOC_DIR, "depmap_CGE.csv"),
                  select = colnames(cge_header)[1])
cat("\nCGE rows (cell lines):", nrow(cge_rows), "\n")
cat("how many overlap with meta$ModelID:",
    sum(cge_rows[[1]] %in% meta$ModelID), "\n")

# Save profile to a small text file for reference
sink(file.path(OUT_DIR, "01_profile_summary.txt"))
cat("Profile run at: ", as.character(Sys.time()), "\n\n")
cat("META:    ", nrow(meta), "rows x", ncol(meta), "cols\n")
cat("VAF:     ", nrow(vaf), "rows x", ncol(vaf), "cols\n")
cat("CGE:     ", nrow(cge_rows), "rows x", ncol(cge_header), "cols\n\n")
cat("VAF unique PIDs:", length(unique(vaf$PID)), "\n")
cat("VAF PID -> meta$ModelID matches:",
    sum(unique(vaf$PID) %in% meta$ModelID), "\n")
cat("VAF PID -> meta$CCLEName matches:",
    sum(unique(vaf$PID) %in% meta$CCLEName), "\n")
cat("VAF PID -> meta$StrippedCellLineName:",
    sum(unique(vaf$PID) %in% meta$StrippedCellLineName), "\n")
cat("CGE rows -> meta$ModelID matches:",
    sum(cge_rows[[1]] %in% meta$ModelID), "\n")
sink()

cat("\nWrote:", file.path(OUT_DIR, "01_profile_summary.txt"), "\n")
