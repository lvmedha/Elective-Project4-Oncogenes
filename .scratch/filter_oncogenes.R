## Filter the SJPedPanel supplementary tables (Karol et al. 2024, CCR-24-1063)
## for oncogenes that are
##   (a) gain-of-function relevant, and
##   (b) recurrent in pediatric SOLID / BRAIN tumors.
##
## Filter rules
##   PED_SOLID_RELEVANT  == TRUE  iff the gene appears in S4C/S4D/S4E with
##                                  Cancer category in {Solid tumor, Brain tumor}.
##                                  (S4 = empirical 113-case cohort.)
##
##   GOF_RELEVANT        == TRUE  iff at least one of:
##                                  - In S4C with class in {missense, in-frame
##                                    ins/del, exon skipping, splice_region}
##                                    AND cancer category in solid/brain.
##                                  - In S4D with cancer category in solid/brain
##                                    (every SV/fusion qualifying for the
##                                    diagnostic yield of a solid/brain case
##                                    is GoF by curation).
##                                  - In S4E with Type of change == GAIN AND
##                                    Size of CNA == focal AND cancer category
##                                    in solid/brain.
##                                  - Marked "Hotspots or SVs only" in ANY of
##                                    the 6 non-SJPedPanel panels in S7B
##                                    (panel designers chose hotspot-only
##                                    coverage => activating oncogene).
##                                  - Listed as a representative oncogenic
##                                    fusion gene in S7D AND covered in
##                                    SJPedPanel.
##
##   NOT_TSG             == TRUE  unless gene is in a curated canonical-TSG
##                                  list AND only has S4E GAIN evidence
##                                  (TSGs are LoF, not GoF).
##
## Outputs (under .scratch/)
##   oncogene_evidence.tsv  one row per gene, all evidence columns
##   oncogene_shortlist.tsv filtered + ranked recommendation
##   oncogene_filter.log    human-readable log

suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
})

XLSX <- "C:/Users/mvijayan/Downloads/ccr-24-1063_supplementary_tables_1_suppts1.xlsx"
OUT_DIR <- ".scratch"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_lines <- character()
log <- function(...) {
  msg <- paste0(...)
  log_lines <<- c(log_lines, msg)
  cat(msg, "\n", sep = "")
}

read_raw <- function(sheet) {
  suppressMessages(read_excel(
    XLSX, sheet = sheet, col_names = FALSE,
    col_types = "text", .name_repair = "minimal"
  ))
}

# -----------------------------------------------------------------------------
# 1. S8 - 183 pediatric cancer driver genes (Ma 2018 + Grobner 2018)
# -----------------------------------------------------------------------------
log("=== S8: 183 pediatric cancer driver genes ===")
s8 <- read_raw("S8 -183 drivers panel comp.")
s8_body <- s8[3:nrow(s8), , drop = FALSE]
s8dt <- data.table(
  gene          = s8_body[[2]],
  ma2018        = !is.na(s8_body[[3]]) & nzchar(s8_body[[3]]),
  ma2018_n      = suppressWarnings(as.numeric(s8_body[[4]])),
  ma2018_pct    = suppressWarnings(as.numeric(s8_body[[5]])),
  grobner2018   = !is.na(s8_body[[6]]) & nzchar(s8_body[[6]]),
  grobner_n     = suppressWarnings(as.numeric(s8_body[[7]])),
  grobner_pct   = suppressWarnings(as.numeric(s8_body[[8]])),
  in_sjped_panel= s8_body[[9]] == "1",
  oncokb        = s8_body[[16]] == "1"
)
s8dt <- s8dt[!is.na(gene) & nzchar(gene)]
log(sprintf("  %d genes; Ma 2018 = %d, Grobner 2018 = %d, on SJPedPanel = %d",
            nrow(s8dt), sum(s8dt$ma2018), sum(s8dt$grobner2018),
            sum(s8dt$in_sjped_panel)))

# -----------------------------------------------------------------------------
# 2. S7B - hotspots-only flag across the 6 NON-SJPedPanel panels
#    SJPedPanel covers full CDS for all 357 genes, so its own
#    hotspots-only flag is uninformative. Other panels' designers
#    chose to selectively sequence only mutational hotspots / SV
#    breakpoints for genes whose oncogenic role is activating /
#    fusion-driven => strong GoF prior.
#
#    Column indices (1-based) of "Hotspots or SVs only" by panel:
#       SJPedPanel        12 (skip)
#       Foundation 1 Heme 17
#       Foundation 1 CDX  22
#       OncoKids          27
#       Thermo Fisher OCAv3 32
#       MSK-IMPACT        38
#       CHOP CHMP/CSTP    43
# -----------------------------------------------------------------------------
log("=== S7B: hotspots-only flag across 6 commercial / institutional panels ===")
s7b <- read_raw("S7B - Panel gene comparsions")
s7b_body <- s7b[3:nrow(s7b), , drop = FALSE]
hs_cols <- c(17, 22, 27, 32, 38, 43)
hs_panels <- c("FoundationOneHeme", "FoundationOneCDX", "OncoKids",
               "ThermoFisherOCAv3", "MSK-IMPACT", "CHOP CHMP/CSTP")
hs_mat <- sapply(hs_cols, function(j) s7b_body[[j]] == "1")
hs_mat[is.na(hs_mat)] <- FALSE
n_panels_hs_only <- rowSums(hs_mat)
hs_panel_names <- apply(hs_mat, 1, function(v) {
  paste(hs_panels[v], collapse = ",")
})
s7b_dt <- data.table(
  gene                       = s7b_body[[2]],
  in_oncokb_cancer_gene_list = s7b_body[[4]] == "Y",
  oncokb_annotated           = s7b_body[[5]] == "Y",
  sjped_covered_S7B          = s7b_body[[6]] == "1",
  hotspots_only_n_panels     = n_panels_hs_only,
  hotspots_only_panels       = hs_panel_names,
  hotspots_only_any          = n_panels_hs_only > 0
)
s7b_dt <- s7b_dt[!is.na(gene) & nzchar(gene)]
log(sprintf("  %d genes; with hotspots-only in >=1 panel = %d (>=2 panels = %d)",
            nrow(s7b_dt), sum(s7b_dt$hotspots_only_any),
            sum(s7b_dt$hotspots_only_n_panels >= 2)))

# -----------------------------------------------------------------------------
# 3. S7D - Reported pediatric oncogenic fusions
# -----------------------------------------------------------------------------
log("=== S7D: 73 representative pediatric fusion oncogenes ===")
s7d <- read_raw("S7D - Reported Fusions Coverage")
s7d_body <- s7d[2:nrow(s7d), , drop = FALSE]
s7d_dt <- data.table(
  gene             = s7d_body[[2]],
  fusion_partners  = s7d_body[[3]],
  n_partners       = suppressWarnings(as.integer(s7d_body[[4]])),
  fusion_panel_cov = s7d_body[[5]] == "Y"
)
s7d_dt <- s7d_dt[!is.na(gene) & nzchar(gene)]
log(sprintf("  %d genes; covered in SJPedPanel = %d",
            nrow(s7d_dt), sum(s7d_dt$fusion_panel_cov, na.rm = TRUE)))

# -----------------------------------------------------------------------------
# 4. S4C - SNV/InDel diagnostic yield: solid+brain GoF-class events
# -----------------------------------------------------------------------------
log("=== S4C: SNV/InDel diagnostic yield (113-case cohort) ===")
s4c <- read_raw("S4C-Diagnostic yield SNVInDel")
s4c_body <- s4c[2:nrow(s4c), , drop = FALSE]
colnames(s4c_body) <- as.character(unlist(s4c[1, ]))
s4c_dt <- as.data.table(s4c_body)
setnames(s4c_dt, c("Cancer category", "Cancer subtype", "Gene", "Class",
                   "AA Change", "Medal"),
         c("cancer_cat", "cancer_subtype", "gene", "class", "aa", "medal"))
s4c_dt[, solid_or_brain := cancer_cat %in% c("Solid tumor", "Brain tumor")]
s4c_dt[, gof_class := grepl("missense|in-?frame|proteinIns|proteinDel|exon[ _]?skip",
                             class, ignore.case = TRUE)]
s4c_solid_gof <- s4c_dt[solid_or_brain == TRUE & gof_class == TRUE,
                        .(n_s4c_solid_gof = .N,
                          aa_examples_s4c = paste(unique(aa[!is.na(aa)]), collapse = ","),
                          subtypes_s4c    = paste(sort(unique(cancer_subtype)),
                                                  collapse = ",")),
                        by = gene]
log(sprintf("  %d solid/brain GoF-class events across %d genes",
            sum(s4c_solid_gof$n_s4c_solid_gof), nrow(s4c_solid_gof)))

# -----------------------------------------------------------------------------
# 5. S4D - SV/fusion diagnostic yield: solid+brain
# -----------------------------------------------------------------------------
log("=== S4D: SV/fusion diagnostic yield ===")
s4d <- read_raw("S4D - Diagnostic yield SV")
s4d_body <- s4d[2:nrow(s4d), , drop = FALSE]
colnames(s4d_body) <- as.character(unlist(s4d[1, ]))
s4d_dt <- as.data.table(s4d_body)
setnames(s4d_dt, c("Cancer category", "Cancer subtype", "Gene/Fusion Pair",
                   "Fusion/SV"),
         c("cancer_cat", "cancer_subtype", "pair", "type"))
s4d_dt[, solid_or_brain := cancer_cat %in% c("Solid tumor", "Brain tumor")]
explode_partners <- function(x) {
  parts <- unique(unlist(strsplit(x, "::|--|/")))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}
s4d_long <- s4d_dt[solid_or_brain == TRUE, {
  parts <- explode_partners(pair)
  list(gene = parts, subtype = cancer_subtype, type = type)
}, by = .I]
s4d_solid <- s4d_long[, .(n_s4d_solid = .N,
                           subtypes_s4d = paste(sort(unique(subtype)),
                                                collapse = ",")),
                       by = gene]
log(sprintf("  %d solid/brain SV/fusion gene-mentions across %d unique genes",
            nrow(s4d_long), nrow(s4d_solid)))

# -----------------------------------------------------------------------------
# 6. S4E - CNV: focal GAIN in solid/brain
# -----------------------------------------------------------------------------
log("=== S4E: CNV diagnostic yield ===")
s4e <- read_raw("S4E - Diagnostic yield CNV")
s4e_body <- s4e[3:nrow(s4e), , drop = FALSE]
colnames(s4e_body) <- as.character(unlist(s4e[2, ]))
s4e_dt <- as.data.table(s4e_body)
setnames(s4e_dt, c("Cancer category", "Cancer subtype", "Gene/Locus",
                   "Type of change", "Size of CNA"),
         c("cancer_cat", "cancer_subtype", "locus", "change", "size"))
s4e_dt[, solid_or_brain := cancer_cat %in% c("Solid tumor", "Brain tumor")]
s4e_dt[, gain := toupper(change) == "GAIN"]
s4e_dt[, focal := grepl("Gene-centric|Focal", size, ignore.case = TRUE)]
s4e_solid_gain <- s4e_dt[solid_or_brain == TRUE & gain == TRUE & focal == TRUE,
                          .(n_s4e_focal_gain = .N,
                            subtypes_s4e = paste(sort(unique(cancer_subtype)),
                                                 collapse = ",")),
                          by = .(gene = locus)]
log(sprintf("  %d solid/brain focal-GAIN events across %d loci",
            sum(s4e_solid_gain$n_s4e_focal_gain), nrow(s4e_solid_gain)))

# -----------------------------------------------------------------------------
# 7. Universe and merge
# -----------------------------------------------------------------------------
# Drop non-gene rows in S4 (e.g. "Intergenic") and split combined CNV
# loci like "PDGFRA/KIT" into their component genes.
not_a_gene <- c("Intergenic", "intergenic", "")
explode_combined <- function(dt, gene_col = "gene") {
  if (!nrow(dt)) return(dt)
  dt[, (gene_col) := strsplit(get(gene_col), "/")]
  dt <- dt[, lapply(.SD, function(x) if (is.list(x)) unlist(x) else
                     rep(x, lengths(dt[[gene_col]])))]
  dt
}
s4c_solid_gof <- s4c_solid_gof[!gene %in% not_a_gene]
s4d_solid     <- s4d_solid[!gene %in% not_a_gene]
s4e_solid_gain <- s4e_solid_gain[!gene %in% not_a_gene]
# Split slash-combined CNV loci ("PDGFRA/KIT") -> separate rows
if (any(grepl("/", s4e_solid_gain$gene))) {
  rows_with_slash <- grep("/", s4e_solid_gain$gene)
  expanded <- s4e_solid_gain[rows_with_slash][,
    .(gene = unlist(strsplit(gene, "/"))),
    by = setdiff(names(s4e_solid_gain), "gene")]
  s4e_solid_gain <- rbind(
    s4e_solid_gain[-rows_with_slash],
    expanded[, names(s4e_solid_gain), with = FALSE]
  )[, .(n_s4e_focal_gain = sum(n_s4e_focal_gain),
        subtypes_s4e     = paste(unique(subtypes_s4e), collapse = ",")),
     by = gene]
}

universe <- unique(c(s4c_solid_gof$gene, s4d_solid$gene, s4e_solid_gain$gene,
                      s7d_dt[fusion_panel_cov == TRUE]$gene))
universe <- universe[!is.na(universe) & nzchar(universe) & !universe %in% not_a_gene]
log(sprintf("=== Universe (any solid/brain or pediatric-fusion evidence): %d genes ===",
            length(universe)))

ev <- data.table(gene = universe)
ev <- merge(ev, s8dt, by = "gene", all.x = TRUE)
ev <- merge(ev, s7b_dt, by = "gene", all.x = TRUE)
ev <- merge(ev, s7d_dt, by = "gene", all.x = TRUE)
ev <- merge(ev, s4c_solid_gof, by = "gene", all.x = TRUE)
ev <- merge(ev, s4d_solid, by = "gene", all.x = TRUE)
ev <- merge(ev, s4e_solid_gain, by = "gene", all.x = TRUE)

num_zero <- function(x) { x[is.na(x)] <- 0; x }
for (cc in c("n_s4c_solid_gof", "n_s4d_solid", "n_s4e_focal_gain",
             "hotspots_only_n_panels", "n_partners",
             "ma2018_n", "ma2018_pct", "grobner_n", "grobner_pct")) {
  if (cc %in% names(ev)) ev[[cc]] <- num_zero(ev[[cc]])
}
yn <- function(x) ifelse(is.na(x) | !x, "", "Y")
for (cc in c("ma2018", "grobner2018", "in_sjped_panel", "oncokb",
             "in_oncokb_cancer_gene_list", "oncokb_annotated",
             "sjped_covered_S7B", "hotspots_only_any", "fusion_panel_cov")) {
  if (cc %in% names(ev)) ev[[cc]] <- yn(ev[[cc]])
}

# -----------------------------------------------------------------------------
# 8. Apply filter
# -----------------------------------------------------------------------------
ev[, ped_solid_relevant := (n_s4c_solid_gof + n_s4d_solid + n_s4e_focal_gain) > 0]
ev[, gof_relevant := (n_s4c_solid_gof > 0 |
                       n_s4d_solid    > 0 |
                       n_s4e_focal_gain > 0 |
                       hotspots_only_any == "Y" |
                       fusion_panel_cov == "Y")]
ev[, both := ped_solid_relevant & gof_relevant]

canonical_tsgs <- c(
  "TP53","RB1","CDKN2A","CDKN2B","NF1","NF2","ATRX","DAXX",
  "ARID1A","ARID1B","ARID2","SMARCA4","SMARCB1","BAP1",
  "PTEN","APC","BRCA1","BRCA2","ATM","CHEK2","PALB2",
  "MSH2","MLH1","MSH6","PMS2","KEAP1","STK11","FBXW7",
  "FAT1","FAT3","TET2","DNMT3A","ASXL1","BCOR","CREBBP",
  "EP300","KMT2C","KMT2D","SETD2","WT1","DICER1","SUFU",
  "PTCH1","PTCH2","NOTCH1","NOTCH2","CIC","FUBP1","TSC1","TSC2",
  "VHL","SMAD4","SMAD2","CDH1","RUNX1","ETV6","PAX5",
  "IKZF1","IKZF3","NCOR1","NCOR2","KDM6A","KDM5A","RNF43",
  "PHF6","ZFHX3","ZFHX4","RAD21","STAG2","SPEN","NSD2",
  "MGA","KAT6A","KAT6B","CDKN1B","CDKN1A","RAD51C","BRIP1",
  "MRE11","NBN","MSH3","ELF3","CDC73","MEN1","NF1","NF2",
  "POLE","POLD1","BLM","WRN"
)
ev[, looks_tsg := gene %in% canonical_tsgs]
# Keep if both criteria met AND not canonical TSG. Drop singletons with
# zero independent corroboration (only 1 S4 event AND no panel signal AND
# not a pediatric driver gene) -- usually fusion partners or noise.
ev[, kept := both & !looks_tsg]
ev[, weak_singleton := (n_s4c_solid_gof + n_s4d_solid + n_s4e_focal_gain) <= 1 &
                       hotspots_only_n_panels == 0 &
                       fusion_panel_cov != "Y" &
                       ma2018 != "Y" &
                       grobner2018 != "Y"]
ev[weak_singleton == TRUE, kept := FALSE]

# Score for ranking
ev[, score := 0]
ev[, score := score + 2 * (ma2018 == "Y")]
ev[, score := score + 1 * (grobner2018 == "Y")]
ev[, score := score + pmin(n_s4c_solid_gof, 10) * 1.0]
ev[, score := score + pmin(n_s4d_solid,    10) * 1.5]
ev[, score := score + pmin(n_s4e_focal_gain, 10) * 1.0]
ev[, score := score + 0.75 * hotspots_only_n_panels]
ev[, score := score + 1.5 * (fusion_panel_cov == "Y")]

# What evidence triggered each gene? (compact tag)
ev[, evidence_tags := {
  tags <- character(.N)
  tags <- ifelse(n_s4c_solid_gof > 0,
                 paste0(tags, ifelse(nchar(tags), ";", ""),
                        "S4C(", n_s4c_solid_gof, ")"), tags)
  tags <- ifelse(n_s4d_solid > 0,
                 paste0(tags, ifelse(nchar(tags), ";", ""),
                        "S4D(", n_s4d_solid, ")"), tags)
  tags <- ifelse(n_s4e_focal_gain > 0,
                 paste0(tags, ifelse(nchar(tags), ";", ""),
                        "S4E_GAIN(", n_s4e_focal_gain, ")"), tags)
  tags <- ifelse(hotspots_only_any == "Y",
                 paste0(tags, ifelse(nchar(tags), ";", ""),
                        "S7B_hotspots_only(", hotspots_only_n_panels, "p)"), tags)
  tags <- ifelse(fusion_panel_cov == "Y",
                 paste0(tags, ifelse(nchar(tags), ";", ""), "S7D_fusion"), tags)
  tags <- ifelse(ma2018 == "Y",
                 paste0(tags, ifelse(nchar(tags), ";", ""), "S8_Ma2018"), tags)
  tags <- ifelse(grobner2018 == "Y",
                 paste0(tags, ifelse(nchar(tags), ";", ""), "S8_Grobner2018"), tags)
  tags
}]

# -----------------------------------------------------------------------------
# Tier classification
#   A: pediatric driver list (Ma OR Grobner) AND empirical solid/brain
#      event in this cohort (S4C or S4D or S4E) AND GoF-flavored panel
#      design (hotspots_only in >= 2 panels OR S7D-fusion).
#   B: empirical solid/brain event AND (pediatric driver OR GoF-flavored
#      panel design); covers strong fusion oncogenes (EWSR1, ZFTA, ...)
#      not in S8 and gene-level GAIN drivers.
#   C: empirical solid/brain event but no other support (often partner
#      genes of fusions; biology to be confirmed individually).
# -----------------------------------------------------------------------------
ev[, has_s4_evidence := (n_s4c_solid_gof + n_s4d_solid + n_s4e_focal_gain) > 0]
ev[, in_peddep_drivers := ma2018 == "Y" | grobner2018 == "Y"]
ev[, gof_panel_signal  := hotspots_only_n_panels >= 2 | fusion_panel_cov == "Y"]
ev[, tier := fcase(
  has_s4_evidence & in_peddep_drivers & gof_panel_signal, "A",
  has_s4_evidence & (in_peddep_drivers | gof_panel_signal), "B",
  has_s4_evidence, "C",
  default = "D"
)]
setorder(ev, -kept, tier, -score, gene)

# Pretty subtypes (truncate)
trunc_str <- function(x, n = 35) ifelse(is.na(x), "",
                                         ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x))

# -----------------------------------------------------------------------------
# 9. Pilot-set sanity check  (pull from the broader S7B/S7D/S8 universe so
#    we see pilot genes even if they have no S4 solid/brain hit)
# -----------------------------------------------------------------------------
PILOT <- c("KRAS","NRAS","HRAS","BRAF","PIK3CA","MAP2K1","EGFR","ERBB2","ALK","CTNNB1")
log("=== Pilot-set sanity check (your current 10 genes) ===")
all_supp <- merge(s8dt[, .(gene, ma2018, grobner2018, in_sjped_panel, oncokb)],
                  s7b_dt[, .(gene, hotspots_only_n_panels, fusion_panel_cov_s7b = NA)],
                  by = "gene", all = TRUE)
all_supp <- merge(all_supp, s7d_dt[, .(gene, fusion_panel_cov)],
                  by = "gene", all.x = TRUE)
all_supp <- merge(all_supp, s4c_solid_gof[, .(gene, n_s4c_solid_gof, subtypes_s4c)],
                  by = "gene", all.x = TRUE)
all_supp <- merge(all_supp, s4d_solid[, .(gene, n_s4d_solid)], by = "gene", all.x = TRUE)
all_supp <- merge(all_supp, s4e_solid_gain[, .(gene, n_s4e_focal_gain)],
                  by = "gene", all.x = TRUE)
for (cc in c("hotspots_only_n_panels","n_s4c_solid_gof","n_s4d_solid","n_s4e_focal_gain")) {
  all_supp[[cc]] <- num_zero(all_supp[[cc]])
}
all_supp[, ma2018          := yn(ma2018)]
all_supp[, grobner2018     := yn(grobner2018)]
all_supp[, fusion_panel_cov := yn(fusion_panel_cov)]
ps_show <- all_supp[gene %in% PILOT][order(match(gene, PILOT)),
                  .(gene, ma2018, grobner2018,
                    hsP = hotspots_only_n_panels,
                    fus = fusion_panel_cov,
                    S4C = n_s4c_solid_gof, S4D = n_s4d_solid, S4E = n_s4e_focal_gain,
                    subtypes_s4c = trunc_str(subtypes_s4c, 35))]
print(ps_show)
log("  (Genes missing from this table have no evidence in S4C/S4D/S4E
       solid+brain rows -- usually because their solid-tumor mutations
       are predominantly in adult cancers or the cohort just didn't
       sample them in solid/brain cases.)")

# -----------------------------------------------------------------------------
# 10. Save and print final shortlist
# -----------------------------------------------------------------------------
fwrite(ev, file.path(OUT_DIR, "oncogene_evidence.tsv"), sep = "\t")

shortlist <- ev[kept == TRUE]
fwrite(shortlist, file.path(OUT_DIR, "oncogene_shortlist.tsv"), sep = "\t")

log(sprintf("=== Final shortlist after GoF + pediatric-solid filter: %d genes ===",
            nrow(shortlist)))
sl <- shortlist[, .(gene, score,
                     ma2018, grobner2018,
                     in_sjped = in_sjped_panel,
                     hsP = hotspots_only_n_panels,
                     fus = fusion_panel_cov,
                     S4C = n_s4c_solid_gof,
                     S4D = n_s4d_solid,
                     S4E = n_s4e_focal_gain,
                     subtypes_s4c  = trunc_str(subtypes_s4c, 30),
                     subtypes_s4d  = trunc_str(subtypes_s4d, 25),
                     subtypes_s4e  = trunc_str(subtypes_s4e, 20),
                     evidence_tags = trunc_str(evidence_tags, 90))]
print(sl, nrows = 200)

writeLines(log_lines, file.path(OUT_DIR, "oncogene_filter.log"))
cat("\nWrote:\n",
    " - .scratch/oncogene_evidence.tsv  (universe with all evidence)\n",
    " - .scratch/oncogene_shortlist.tsv (filtered + ranked)\n",
    " - .scratch/oncogene_filter.log\n", sep = "")
