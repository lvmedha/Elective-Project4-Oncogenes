# 00_build_target_genes.R
# Build the canonical target_genes table for the rest of the pipeline.
#
# Source of the gene list:
#   - data/oncogene_shortlist_sjpedpanel.tsv  (Tiers A/B/C, 42 genes; output of
#     filter_oncogenes_sjpedpanel.R parsing the SJPedPanel supp. tables, Karol
#     et al. 2024, CCR-24-1063).
#
# This script:
#   1. Reads the shortlist.
#   2. Looks up each gene's GRCh38 coordinates via the Ensembl REST API
#      (POST /lookup/symbol/homo_sapiens, batched; per-gene file cached
#      under data/cache/ so re-runs are free).
#   3. Writes data/target_genes.tsv, the SINGLE SOURCE OF TRUTH used by
#      02..07 in place of the previous hard-coded 10-gene blocks.
#
# Output columns:
#   gene chrom start end strand ensembl_id biotype tier score evidence_tags

suppressPackageStartupMessages({
  for (p in c("data.table","httr2","jsonlite","digest")) {
    if (!requireNamespace(p, quietly = TRUE))
      install.packages(p, repos = "https://cloud.r-project.org")
  }
  library(data.table)
  library(httr2)
  library(jsonlite)
  library(digest)
})

PROJ_DIR  <- "C:/Users/mvijayan/Documents/Elective-Project4-Oncogenes"
DATA_DIR  <- file.path(PROJ_DIR, "data")
CACHE_DIR <- file.path(DATA_DIR, "cache")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

SHORTLIST <- file.path(DATA_DIR, "oncogene_shortlist_sjpedpanel.tsv")
OUT       <- file.path(DATA_DIR, "target_genes.tsv")

stopifnot(file.exists(SHORTLIST))
sl <- fread(SHORTLIST)
cat("Shortlist rows:", nrow(sl), "\n")
cat("Tiers:\n"); print(sl[, .N, by = tier])

# --- Ensembl REST: POST /lookup/symbol/homo_sapiens (batch up to 1000) ------
# Response shape: { "GENE": { id, seq_region_name, start, end, strand, biotype, ... }, ... }
# We cache by hash of the request body so re-runs are free.
ensembl_lookup_symbols <- function(symbols, server = "https://rest.ensembl.org") {
  body <- list(symbols = symbols)
  body_json <- toJSON(body, auto_unbox = TRUE)
  cache_key <- digest(body_json, algo = "md5")
  cache_file <- file.path(CACHE_DIR, paste0("lookup_", cache_key, ".json"))
  if (file.exists(cache_file)) {
    return(fromJSON(cache_file, simplifyVector = FALSE))
  }
  resp <- request(paste0(server, "/lookup/symbol/homo_sapiens")) |>
    req_method("POST") |>
    req_headers(`Content-Type` = "application/json",
                Accept         = "application/json") |>
    req_body_raw(body_json) |>
    req_retry(max_tries = 4, backoff = ~ 5) |>
    req_timeout(120) |>
    req_perform()
  result <- resp_body_json(resp)
  writeLines(toJSON(result, auto_unbox = TRUE, pretty = FALSE), cache_file)
  result
}

# HGNC re-named several histone genes in 2019; Ensembl follows the new
# symbols. Map legacy aliases (used by the SJPedPanel paper and most
# tumor-mutation databases) to the current symbol so the lookup hits.
aliases <- c(
  "H3F3A"  = "H3-3A",
  "H3F3B"  = "H3-3B",
  "HIST1H3B" = "H3C2",
  "HIST1H3C" = "H3C3"
)
sl[, lookup_symbol := ifelse(gene %in% names(aliases),
                              aliases[gene], gene)]

cat("\nLooking up coordinates via Ensembl REST...\n")
res <- ensembl_lookup_symbols(unique(sl$lookup_symbol))

# Reverse-map: keys in `res` are the looked-up symbols; the gene column we
# write back uses the original (legacy-friendly) name.
found_lookup <- names(res)
sl[, found := lookup_symbol %in% found_lookup]
missed <- sl[found == FALSE, gene]
if (length(missed)) {
  cat("WARNING: Ensembl returned no record for:",
      paste(missed, collapse = ","), "\n")
  cat("        These will be dropped from the target list.\n")
}

coords <- rbindlist(lapply(sl[found == TRUE]$gene, function(g) {
  lk <- aliases[g]; if (is.na(lk)) lk <- g
  r <- res[[lk]]
  data.table(
    gene       = g,
    chrom      = paste0("chr", r$seq_region_name),
    start      = as.integer(r$start),
    end        = as.integer(r$end),
    strand     = as.integer(r$strand),
    ensembl_id = r$id,
    biotype    = r$biotype
  )
}), fill = TRUE)

# Merge with shortlist metadata (tier, score, evidence_tags)
keep_cols <- intersect(c("gene","tier","score","evidence_tags",
                         "ma2018","grobner2018","in_sjped_panel",
                         "hotspots_only_n_panels","fusion_panel_cov",
                         "n_s4c_solid_gof","n_s4d_solid","n_s4e_focal_gain",
                         "subtypes_s4c","subtypes_s4d","subtypes_s4e"),
                       names(sl))
target <- merge(coords, sl[, ..keep_cols], by = "gene", all.x = TRUE)

# Sort: Tier A first, then by score desc within tier
target[, tier := factor(tier, levels = c("A","B","C","D"))]
setorder(target, tier, -score, gene)

# Sanity: drop any rows that didn't get a chrom (shouldn't happen)
target <- target[!is.na(chrom) & !is.na(start) & !is.na(end)]

cat(sprintf("\nFinal target_genes: %d genes (Tier A=%d, B=%d, C=%d)\n",
            nrow(target),
            sum(target$tier == "A", na.rm = TRUE),
            sum(target$tier == "B", na.rm = TRUE),
            sum(target$tier == "C", na.rm = TRUE)))

cat("\nTier A:\n"); print(target[tier == "A", .(gene, chrom, start, end, score)])
cat("\nTier B:\n"); print(target[tier == "B", .(gene, chrom, start, end, score)])
cat("\nTier C:\n"); print(target[tier == "C", .(gene, chrom, start, end, score)])

fwrite(target, OUT, sep = "\t")
cat("\nWrote:", OUT, "\n")
