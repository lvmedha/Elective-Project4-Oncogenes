# Project 4 — Step 1 status

## What's done

The `AllMarkers_VAF_long.tsv` has now been re-shaped from
`PID + chr.pos.ref.alt` rows into a tidy gene/protein-level mutation
table for **10 target oncogenes**:

`KRAS, NRAS, HRAS, BRAF, PIK3CA, MAP2K1, EGFR, ERBB2, ALK, CTNNB1`

### Pipeline (see `R/`)

1. `01_profile_inputs.R` — sanity-checks all three input files.
2. `02_inspect_vaf.R` — explores VAF; counts markers in canonical
   oncogene loci to choose the 10 target genes.
3. `03_annotate_variants.R` — submits all coding variants in those
   loci to **Ensembl VEP REST API (GRCh38)**; restricts to MANE Select
   / canonical transcripts; emits `HugoSymbol`, `HGVSp_Short`,
   `Variant_Class`.
4. `04_hotspot_summary.R` — recurrence summary per gene.
5. `05_apply_hotspot_tiers.R` — joins the curated Tier-0 hotspot file
   `data/hotspots_tier0.csv` and assigns `Mut_Status`
   (`Hotspot / Missense_Other / InframeIndel / Truncating / Silent`).

### Outputs (`results/`)

- `03_vaf_annotated.tsv` — every coding variant in the 10 target
  genes with VEP annotation.
- `04_recurrent_missense_per_gene.tsv` — per-gene counts of
  recurrent protein changes (useful for sanity-checking & for
  curating Tier-0).
- `05_mutations_tiered.tsv` — **the file to use for the waterfall
  plot.** One row per `(PID, gene_symbol)` carrying a
  hotspot/missense/inframe/truncating event.

### Per-gene PID counts (Tier-0 hotspot vs other missense)

| Gene    | Hotspot | Missense_Other | InframeIndel | Truncating |
|---------|--------:|---------------:|-------------:|-----------:|
| KRAS    | 126     | 9              | 2            | 0          |
| NRAS    | 89      | 1              | 0            | 0          |
| PIK3CA  | 71      | 17             | 3            | 0          |
| BRAF    | 57      | 7              | 2            | 1          |
| CTNNB1  | 18      | 4              | 2            | 2          |
| ALK     | 11      | 1              | 0            | 0          |
| HRAS    | 10      | 5              | 0            | 0          |
| MAP2K1  | 9       | 2              | 0            | 0          |
| ERBB2   | 6       | 13             | 1            | 0          |
| EGFR    | 4       | 4              | 7            | 0          |

Biology checks pass: BRAF V600E in 50 patients, KRAS G12D 38, NRAS
Q61K 25, PIK3CA H1047R 22, ALK F1174L/R1275Q in `SJNBL` neuroblastoma
patients, EGFR L858R/T790M in `SJLUAD`. The Tier-0 codon list captures
> 75% of missense hits in 9 / 10 genes (ERBB2 is the exception — only
32% — worth reviewing).

## Blockers / open questions for Declan

### 1. Patient ↔ cell-line linkage (the big one)

`AllMarkers_VAF_long.tsv` is keyed by St. Jude tumor patient IDs
(`SJST*`, `SJOS*`, `SJNBL*`, …), **not** DepMap `ACH-` ModelIDs.
Of the 1,095 unique PIDs in the VAF, **0** match anything in
`depmap_meta.csv` directly.

The PedDep CRISPR file `depmap_CGE.csv` is keyed on `ACH-` IDs and
its 1,208 cell lines all match `depmap_meta.csv` cleanly.

**This means we cannot currently compare mutated-vs-WT *cell-line*
dependency**, because mutations are on tumors and CRISPR scores are
on cell lines. Possibilities to discuss:

- (a) Is there a separate PID → ACH mapping file for the PedDep cell
  collection that we should be using?
- (b) Should we treat this VAF as the *tumor cohort projection* layer
  (Figure 4) and use a different DepMap/PedDep mutation file for the
  cell-line analysis (Figure 3)?
- (c) If (b), the obvious file is DepMap's `OmicsSomaticMutations.csv`
  (already keyed on `ModelID`, already has `HugoSymbol` and
  `HGVSp_Short` and `LikelyLoF`/`Hotspot` flags).

### 2. ERBB2 hotspot list looks short

Only 6/19 ERBB2 missense PIDs hit a Tier-0 codon. Worth checking
whether we should add ERBB2 R678Q (recurrent in our data, 2 PIDs)
and any kinase-domain residues from the recent ERBB2 hotspot
literature.

### 3. AlphaMissense layer

Not yet integrated. Once 1 is resolved we can layer
`AlphaMissense_aa_substitutions.tsv.gz` onto the
`(transcript_id, protein_start, aa_alt)` key in
`results/03_vaf_annotated.tsv`. Then we'll have the AM tier alongside
the Tier-0 tier and can do the AM-evaluation analysis we planned.
