# CLAUDE.md

Persistent project context for Claude / Cursor agent sessions in this repo.
For human-facing docs read [`README.md`](./README.md) and
[`results/STATUS.md`](./results/STATUS.md) first — they are the source of
truth.

---

## What this project is

Course project analyzing oncogene missense variants across pediatric tumor
cohorts and DepMap cell lines, plus an **AlphaMissense + VEP** annotation
layer. Two analytical layers:

1. **R pipeline** (`R/00..R/10`) — pure R + Ensembl VEP REST. Produces
   tiered mutation tables, waterfall plots, AM-vs-Chronos figures.
2. **AlphaMissense / local VEP layer** — optional, runs in WSL/conda
   (`alpha_tools` env). The cluster job produces
   `results/cohort_full/09_vep_full_alphamissense.vcf.gz`, which `R/09`
   joins back into the R outputs.

Default `RUN_NAME` in every numbered script is `"ped_gof_snv"` →
outputs land under `results/ped_gof_snv/`. The 10-gene pilot is
preserved at `results/pilot_10genes/` for comparison; don't overwrite it.

---

## R pipeline ordering

| Step | Script | Reads | Writes (under `results/<RUN_NAME>/`) |
| --- | --- | --- | --- |
| 00 | `R/00_build_target_genes.R` | `data/oncogene_shortlist_ped_gof_snvs.tsv` | `data/target_genes.tsv` |
| 01 | `R/01_profile_inputs.R` | VAF table, DepMap CSVs | sanity-check stdout |
| 02 | `R/02_inspect_vaf.R` | VAF table | VAF QC |
| 03 | `R/03_annotate_variants.R` | VAF + Ensembl VEP REST | `03_vaf_annotated.tsv` |
| 04 | `R/04_hotspot_summary.R` | `03_vaf_annotated.tsv`, `data/hotspots_tier0.csv` | `04_recurrent_missense_per_gene.tsv` |
| 05 | `R/05_apply_hotspot_tiers.R` | 04 + `data/hotspots_tier0.csv` | `05_mutations_tiered.tsv` |
| 06 | `R/06_link_to_depmap.R` | 05 + `samples.txt` | `06_mutations_with_ACH.tsv` |
| 07 | `R/07_waterfall_plots.R` | 06 + DepMap CGE | `07_waterfalls_*.pdf`, `07_effect_sizes.tsv` |
| 08 | `R/08_make_vep_input.R` (optional) | 03 | `results/cohort_full/08_vep_input_full.{vcf,ensembl}` |
| 09 | `R/09_am_dependency_figures.R` | `09_vep_full_alphamissense.vcf.gz` + 06 + DepMap CGE | `09_am_chronos_*.{tsv,pdf}` |
| 10 | `R/10_prism_mutation_sensitivity.R` (optional) | PRISM CSVs + 06 | `10_prism_*.tsv`, `10_prism_volcano_*.pdf` |
| 11 | `R/11_myc_cnv_analysis.R` | `~/Documents/depmap_cnv.csv` + 09 events | `11_myc_cnv_*.{tsv,txt,pdf}` |

Steps 09–10 need files step 08 (or the cluster) produces. Step 10 needs
PRISM CSVs in `~/Documents/`. Step 11 is a MYC-only follow-up to 09 that
joins DepMap MYC copy number to the cohort's MYC missense events to
explain why some AlphaMissense-likely-benign variants sit on highly
MYC-dependent lines (answer: amplification or IG-MYC translocation).
The earlier triage script (`R/11_am_missense_triage_figures.R`) was
removed in May 2026 — do not recreate it under that name. The new `11`
slot is the MYC CNV analysis.

---

## How to run things

```powershell
# Run any single step
Rscript R/07_waterfall_plots.R

# Smoke-test the cluster VEP+AM job
bash scripts/run_vep_alphamissense_smoke.sh
```

All scripts are idempotent; rerunning overwrites the previous run's
PDFs/TSVs under the same `RUN_NAME`. To freeze a run, change `RUN_NAME`
at the top of the relevant script (or all of 01..10) **before** re-running.

R packages needed: `data.table`, `httr2`, `jsonlite`, `digest`, `readxl`,
`ggplot2`, `RColorBrewer` (for 07's palette), optional `ggrepel` (for
10's volcano labels).

---

## Conventions

### R style

- `suppressPackageStartupMessages({ library(...) })` at the top of every
  numbered script. Always `library(data.table); library(ggplot2)` at
  minimum; add `RColorBrewer` when using `brewer.pal()` directly.
- Heavy data wrangling: **`data.table`** (`setorder`, `fread`, `:=`,
  by-reference). Avoid `dplyr` / tidyverse — none of the existing scripts
  use it and we're keeping deps minimal.
- File paths: build with `file.path(OUT_DIR, "<name>.<ext>")`, never
  string-concat with `/`. `OUT_DIR <- file.path(PROJ_DIR, "results", RUN_NAME)`.
- Reproducibility: every script defines `PROJ_DIR`, `DATA_DIR`, `OUT_DIR`,
  `RUN_NAME` at the top. Keep that pattern.
- Plots: `theme_minimal(base_size = 11)`, bold `plot.title`. Mutation-tier
  fills use the Set1 palette via `RColorBrewer::brewer.pal(9, "Set1")` —
  see the `pal` block in `R/07_waterfall_plots.R`. Keep `Hotspot = red`,
  `Missense_Other = blue` (don't return to red/orange — too similar).
- AM-faceted scatter (`R/09`): x-axis is locked to `[0, 1]` via
  `scale_x_continuous(limits = c(0, 1))` so panels are comparable. Don't
  switch back to `scales = "free"`.

### Mutation tier vocabulary

`Mut_Status` factor levels, always in this order:
`Hotspot > Missense_Other > InframeIndel > Truncating > Silent > WT`.
Status priority for "most-activating event wins" dedup is encoded as
`status_priority <- c(Hotspot=5, Missense_Other=4, InframeIndel=3,
Truncating=2, Silent=1)`.

### HGNC alias handling

DepMap uses current HGNC names; this panel keeps legacy names (e.g.
`H3F3A`, `HIST1H3B`, `HIST1H3C`). Aliases are normalized in `R/07` and
`R/09` via `hgnc_aliases <- c("H3F3A"="H3-3A", "H3F3B"="H3-3B",
"HIST1H3B"="H3C2", "HIST1H3C"="H3C3")`. Add new aliases there if you
extend the panel.

---

## What's tracked vs ignored

Only the **panel-definition** files in `data/` are tracked:
`target_genes.tsv`, `oncogene_shortlist_ped_gof_snvs.tsv`,
`hotspots_tier0.csv`. Everything else under `data/` is `.gitignore`d
(including `data/cache/` VEP REST cache). The AlphaMissense hg38 table
(`AlphaMissense_hg38.tsv.bgz`) and `/vep_cache/` are also ignored — they
are multi-GB.

The cloned AlphaMissense source repo (`/alphamissense/`) is **hard-ignored**
because Windows-native git can't read WSL symlinks on `/mnt/c`. Don't
re-symlink it into the worktree.

`.scratch/` is a sandbox: `*.out`/`*.log` ignored, `*.sh` tracked. Use it
for ad-hoc shell scripts you don't want in `scripts/`.

---

## External inputs (not in repo)

Paths default to `C:/Users/mvijayan/Documents/`. Adjust `DOC_DIR` at the
top of the relevant R script if your layout differs.

| File | Used by |
| --- | --- |
| `AllMarkers_VAF_long.tsv` (or project-specific path) | 01–06 |
| `depmap_CGE.csv`, `depmap_meta.csv` | 06, 07, 09 |
| `samples.txt` (PID ↔ ACH map) | 06 |
| `results/cohort_full/09_vep_full_alphamissense.vcf.gz` | 09 |
| `PRISMOncology*.csv` | 10 |

---

## Platform notes (Windows / WSL)

- This is a Windows host. R runs natively on Windows. VEP +
  AlphaMissense **must** run in WSL2 / Ubuntu — bioconda doesn't ship
  Windows wheels for the genomics stack.
- Do not commit WSL-side symlinks back into `/mnt/c/...`; Windows git
  treats them as text blobs.
- Don't add Python to the R-only pipeline. The optional AM tooling lives
  in `setup_alpha_tools.sh` / `setup_env.sh` and stays isolated in the
  `alpha_tools` conda env.

---

## Working-with-this-repo etiquette

- Edit existing scripts in place. The numbered ordering is the project's
  contract — adding new steps means picking the next free number, not
  renumbering existing ones.
- When you change a plot script, **regenerate the PDF** so the
  `results/<RUN_NAME>/` outputs stay in sync with the code. Don't commit
  code-only edits that leave stale PDFs behind unless the user asks.
- Don't create documentation files (`*.md`) proactively. `README.md` and
  `results/STATUS.md` are the documentation; update them when behavior
  changes, don't fork them.
- Don't run `git commit` / `git push` unless explicitly asked.
- After substantive R edits, lint-check with `ReadLints` for the touched
  files.

---

## Recent decisions worth remembering

- **May 2026** — removed `R/11_am_missense_triage_figures.R` and all its
  outputs (`results/ped_gof_snv/11_am_*`). Don't reintroduce.
- **May 2026** — switched waterfall palette in `R/07` from a custom
  red/orange scheme to `RColorBrewer::brewer.pal(9, "Set1")` so Hotspot
  (red) and Missense_Other (blue) are visually distinct.
- **May 2026** — locked the AM-vs-Chronos faceted scatter's x-axis to
  `[0, 1]` so panels are comparable across genes.
- **May 2026** — added
  `09_chronos_density_with_am_rug_faceted_by_gene.pdf`: per-gene density
  of the full DepMap Chronos distribution with the cohort's AM-scored
  mutant lines as a rug.
- **May 2026** — added `R/11_myc_cnv_analysis.R` (MYC-only): joins
  `~/Documents/depmap_cnv.csv` to the cohort's MYC events. Closes the
  AlphaMissense-vs-Chronos story for MYC: AM-likely-benign + highly
  dependent lines are explained by MYC amplification (CN > ~1.5 in
  solid-tumor carriers) or by IG-MYC translocation (B-cell lymphoma
  lines, near-diploid CN but high mRNA — CN cannot detect this). Full
  DepMap correlation: Spearman ρ(MYC_CN, MYC_Chronos) = −0.44, n=872.
