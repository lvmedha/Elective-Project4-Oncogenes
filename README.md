# Elective-Project4-Oncogenes

Course project analyzing oncogene missense variants across pediatric
tumor cohorts and the PedDep cell-line panel. Combines:

- An **R pipeline** (under [`R/`](./R)) that takes a long-format VAF
  table, restricts to **24 target oncogenes** (curated pediatric
  gain-of-function SNV/indel panel — see
  [`data/oncogene_shortlist_ped_gof_snvs.tsv`](./data/oncogene_shortlist_ped_gof_snvs.tsv),
  built into [`data/target_genes.tsv`](./data/target_genes.tsv) by
  `R/00_build_target_genes.R`). A **legacy 42-gene** SJPedPanel-derived
  list remains available as
  [`data/oncogene_shortlist_sjpedpanel.tsv`](./data/oncogene_shortlist_sjpedpanel.tsv)
  (switch `SHORTLIST` in `00_build_target_genes.R` to use it). The
  pipeline annotates every coding variant via the **Ensembl VEP** REST API, applies a
  curated Tier-0 hotspot list, and produces tiered mutation tables and
  waterfall plots.
- An **AlphaMissense layer** (Google DeepMind's
  [AlphaMissense](https://github.com/google-deepmind/alphamissense))
  for proteome-wide missense pathogenicity scoring of the same variants.
- A local **Ensembl VEP** install (via bioconda) as a fallback to the
  REST API when batch sizes / network are an issue.

For a step-by-step description of what the R pipeline produces and the
current scientific status, see [`results/STATUS.md`](./results/STATUS.md).

---

## Do I need Python / WSL / AlphaMissense?

**No — not to run the core 00–07 pipeline.** Those steps are self-contained
and need **R + an internet connection** (Ensembl VEP REST). **Steps 09–11**
read **AlphaMissense** from a **cohort VCF** produced on the cluster (VEP +
`--plugin AlphaMissense`); you can regenerate that VCF using
[`scripts/run_vep_alphamissense.sh`](./scripts/run_vep_alphamissense.sh)
after `R/08_make_vep_input.R`. Python, WSL, and conda are **optional** and
only needed for local VEP / AlphaMissense tooling (see below).

The optional **step 08** (`R/08_make_vep_input.R`) is also pure R — it
re-emits the cohort as a `.vcf` / `.ensembl` pair that a cluster VEP
+ AlphaMissense run can consume. The cluster job itself runs outside
this repo (St. Jude HPC), so the only thing you need locally is R.

---

## TL;DR — R-only quick start (no Python required)

This path reproduces tables and figures under `results/ped_gof_snv/`
(default `RUN_NAME` in the numbered `R/` scripts). **External files** you
must supply locally are listed under **External inputs** below (VAF table,
DepMap CSVs, optional PRISM matrices).

```r
# From R, in the project root:
install.packages(c("data.table","httr2","jsonlite","digest","readxl","ggplot2"))
# Optional: nicer PRISM volcano labels
install.packages("ggrepel")

source("R/00_build_target_genes.R")    # build data/target_genes.tsv (24 genes)
source("R/01_profile_inputs.R")
source("R/02_inspect_vaf.R")
source("R/03_annotate_variants.R")     # calls Ensembl VEP REST
source("R/04_hotspot_summary.R")
source("R/05_apply_hotspot_tiers.R")
source("R/06_link_to_depmap.R")
source("R/07_waterfall_plots.R")

# Optional: cohort VCF for cluster VEP + AlphaMissense
# source("R/08_make_vep_input.R")      # -> results/cohort_full/08_vep_input_full.*

# After 09_vep_full_alphamissense.vcf.gz exists under results/cohort_full/:
source("R/09_am_dependency_figures.R")   # AM vs Chronos PDFs + 09_am_chronos_*.tsv

# Optional: PRISM (expects CSVs under ~/Documents; see R/10 header)
source("R/10_prism_mutation_sensitivity.R")
```

Panel files **`data/target_genes.tsv`**, **`data/oncogene_shortlist_ped_gof_snvs.tsv`**,
and **`data/hotspots_tier0.csv`** are **tracked in git**; other `data/*`
(e.g. `cache/`, large TSVs) stays local. For the **legacy** 42-gene panel,
set `SHORTLIST` in `00_build_target_genes.R` to `oncogene_shortlist_sjpedpanel.tsv`
and add that file locally (not tracked by default).

Core **00–07** need only **R + network**. **09** needs the annotated VCF;
**10** needs PRISM files on disk.

---

## (Optional) TL;DR — AlphaMissense / local-VEP setup

Only follow this if you want to add AlphaMissense scores or run VEP
offline. On Windows this **must** be done inside **WSL2 / Ubuntu** (see
the WSL note below). All commands below assume an Ubuntu shell.

```bash
# From a fresh WSL/Ubuntu shell, in this repo's root:
bash setup_alpha_tools.sh    # creates the alpha_tools conda env (VEP + Python 3.11)
bash setup_env.sh --all      # also installs AlphaMissense (chains into setup_alpha_tools.sh)

conda activate alpha_tools
python -c "import alphamissense; from alphamissense.model import config; \
           print('OK', config.model_config().model.num_recycle)"
```

If `conda activate` reports "command not found", source it once with
`source ~/mambaforge/etc/profile.d/conda.sh` (or your conda's path) and
try again.

---

## Repo layout

```
.
├── R/                                # R pipeline 00–10 (see Running section)
│   ├── 00_build_target_genes.R
│   ├── 01_profile_inputs.R … 07_waterfall_plots.R
│   ├── 08_make_vep_input.R           # cohort VCF / Ensembl input
│   ├── 09_am_dependency_figures.R    # AM vs Chronos (needs annotated VCF)
│   └── 10_prism_mutation_sensitivity.R
├── data/                             # tracked: panel TSV/CSV only (see .gitignore)
│   ├── target_genes.tsv
│   ├── oncogene_shortlist_ped_gof_snvs.tsv
│   ├── hotspots_tier0.csv
│   └── cache/                        # local VEP REST cache (not in git)
├── scripts/
│   ├── run_vep_alphamissense.sh      # cluster / local VEP + AM
│   └── run_vep_alphamissense_smoke.sh
├── results/
│   ├── STATUS.md
│   ├── cohort_full/                  # whole-cohort VEP input + AM VCF
│   ├── pilot_10genes/
│   └── ped_gof_snv/                  # default RUN_NAME outputs
├── environment.yml
├── setup_alpha_tools.sh
├── setup_env.sh
└── README.md
```

`results/` is split per-run. **`RUN_NAME`** in each script is currently
`"ped_gof_snv"`. Change it to freeze another run, then re-run **01..07**
(and **08–10** as needed).

The R pipeline is the project's core analysis path; AlphaMissense and
the local VEP install are auxiliary layers that the R outputs feed into.

---

## Environment setup (optional — AlphaMissense / local VEP only)

> Skip this entire section if you just want to run the R pipeline. The
> R scripts call Ensembl VEP over REST and do not need any of the
> environments below.

There are two complementary environments:

| Env | Built by | Contents | Used for |
| --- | --- | --- | --- |
| `alpha_tools` (conda) | `setup_alpha_tools.sh` | Python 3.11, ensembl-vep (+BLAST/HMMER/BioPerl/etc.), open-mpi | Local VEP, hosting AlphaMissense |
| `alphamissense/venv` (pip venv) | `setup_env.sh` | Python 3.11 venv, AlphaMissense + JAX | Optional, separate AlphaMissense install |

**Recommended** path is the conda env (`alpha_tools`) because:

- It already provides Python 3.11 and `jackhmmer`, so AlphaMissense can
  be installed *into it* without `apt install python3.11-venv` and
  without `sudo`.
- It is the env the R pipeline's local-VEP fallback also uses, so
  there's just one thing to `conda activate`.

The pip-venv path (`setup_env.sh`) is preserved as a fallback that
mirrors AlphaMissense's official upstream install instructions.

### Prereq: WSL2 (Windows only)

AlphaMissense and `ensembl-vep` are Linux-only. On Windows you must
run **inside WSL2 / Ubuntu**. From an **Administrator** PowerShell:

```powershell
wsl --install
```

Reboot when prompted, launch **Ubuntu** from the Start menu, create a
Linux username and password. Open the project from inside WSL with:

```bash
cd /mnt/c/Users/<you>/Documents/Elective-Project4-Oncogenes
```

> **WSL filesystem caveat.** The `/mnt/c/...` mount does not support
> POSIX permission bits, so cloning git repos or creating symlinks
> *inside* `/mnt/c/...` from WSL fails with `chmod ... Operation not
> permitted` and `Function not implemented`. The setup scripts therefore
> clone AlphaMissense into the WSL **native** filesystem
> (`~/alphamissense`) and Python's editable install lives there, not
> inside the project folder. The project's R/data/results files stay on
> `/mnt/c/...` so Windows tooling and Git can see them normally.

### Option A — `setup_alpha_tools.sh` (recommended)

```bash
bash setup_alpha_tools.sh                  # create alpha_tools from environment.yml
bash setup_alpha_tools.sh --lock           # use environment.lock.yml (pinned, exact)
bash setup_alpha_tools.sh --with-vep-cache # also fetch the GRCh38 VEP cache (~25 GB)
bash setup_alpha_tools.sh --force          # recreate env even if it exists
```

The script:

1. Detects/sources an existing `mambaforge` / `miniforge3` / `miniconda3` /
   `anaconda3` install, or downloads and installs Miniforge if none exists.
2. Creates the `alpha_tools` conda env from `environment.yml` (high-level)
   or `environment.lock.yml` (exact, ~290 pinned packages).
3. Verifies `vep` and `vep_install` are on PATH inside the env.
4. (Optional) Runs `vep_install` to download the GRCh38 cache into
   `data/vep/` for offline VEP runs.

To install AlphaMissense **into the same env** afterwards (no separate
venv, no sudo, no apt):

```bash
conda activate alpha_tools
conda install -n alpha_tools pip   # alpha_tools doesn't ship pip by default
git clone https://github.com/google-deepmind/alphamissense.git ~/alphamissense
cd ~/alphamissense
# AlphaMissense's requirements.txt pins jaxlib==0.4.14 which has been
# yanked from PyPI; bump to 0.4.18 (oldest available, dm-haiku 0.0.10
# is compatible).
sed -i 's/^jax==0\.4\.14$/jax==0.4.18/'      requirements.txt
sed -i 's/^jaxlib==0\.4\.14$/jaxlib==0.4.18/' requirements.txt
pip install -r requirements.txt
pip install -e .
```

Quick smoke test (from any directory):

```bash
python -c "
import alphamissense
from alphamissense.model import config
from alphamissense.data import pipeline_missense
print('alphamissense OK')
print('default num_recycle:', config.model_config().model.num_recycle)
"
```

> The shipped `test/test_installation.py` hardcodes
> `/usr/bin/jackhmmer` and will fail until that path is patched (or
> until you symlink the env's `jackhmmer` to `/usr/bin/jackhmmer` with
> `sudo`). The above import smoke test does not need `jackhmmer`.

### Option B — `setup_env.sh` (separate pip venv)

The original installer creates a standalone `alphamissense/venv`
inside the cloned repo. This **requires `sudo`** (for `apt-get install
python3.11-venv aria2 hmmer git`) and creates duplicates of Python and
HMMER that `alpha_tools` would already provide.

```bash
bash setup_env.sh                # AlphaMissense venv only
bash setup_env.sh --with-data    # also fetch ~5–10 GB precomputed predictions
bash setup_env.sh --all          # also chain into setup_alpha_tools.sh
bash setup_env.sh --help
```

Activate with:

```bash
cd alphamissense
source venv/bin/activate
```

Both scripts are idempotent — re-running is safe and skips work that's
already complete.

### Optional: auto-activation in WSL

If you want every shell that's `cd`'d into this project to auto-activate
`alpha_tools` (so you don't have to `conda activate` manually), append
this block to `~/.bashrc` inside WSL:

```bash
# >>> alpha_tools project hook (Elective-Project4-Oncogenes) >>>
__alpha_tools_project_root="/mnt/c/Users/<you>/Documents/Elective-Project4-Oncogenes"
__alpha_tools_maybe_activate() {
    case "$PWD" in
        "$__alpha_tools_project_root"|"$__alpha_tools_project_root"/*)
            if [ "$CONDA_DEFAULT_ENV" != "alpha_tools" ]; then
                conda activate alpha_tools >/dev/null 2>&1
            fi
            ;;
    esac
}
__alpha_tools_maybe_activate
PROMPT_COMMAND="__alpha_tools_maybe_activate;${PROMPT_COMMAND:-}"
# <<< alpha_tools project hook <<<
```

After that, opening any WSL shell (Cursor terminal, Windows Terminal,
plain `wsl`, ...) and `cd`-ing into the project will activate the env
silently.

---

## (Optional) Precomputed AlphaMissense predictions

DeepMind publishes precomputed AlphaMissense scores for every possible
human missense substitution in a public Google Cloud Storage bucket.
For most downstream oncogene analyses you only need these files — you
do **not** need to run the model yourself.

The easiest way is `bash setup_env.sh --with-data`, which downloads
the two most commonly used files via `aria2` into `./data/`:

| File | Contents |
| --- | --- |
| `AlphaMissense_aa_substitutions.tsv.gz` | Pathogenicity scores for all possible single-amino-acid substitutions in the human proteome. |
| `AlphaMissense_hg38.tsv.gz` | Same predictions mapped to GRCh38 (hg38) coordinates. |

To download manually instead:

```text
https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
```

Full bucket: <https://console.cloud.google.com/storage/browser/dm_alphamissense>

---

## (Optional) Genetic databases

If you intend to *run* the AlphaMissense data pipeline (not just use
the precomputed predictions), you will additionally need the genetic
sequence databases (BFD, MGnify, UniRef90). Follow the instructions in
the [AlphaFold repository](https://github.com/google-deepmind/alphafold)
to download them. These are large (hundreds of GB) and only required
for inference from raw sequences.

---

## Running the R pipeline

This is the project's core analysis path and **does not require Python,
WSL, or conda** — see the R-only TL;DR above for the minimum setup.

The numbered scripts in [`R/`](./R) form a linear pipeline. Each script
reads from `data/` (and paths documented in **External inputs**) and writes
into `results/<RUN_NAME>/` (and `results/cohort_full/` for step 08).

```r
source("R/00_build_target_genes.R")
source("R/01_profile_inputs.R")
source("R/02_inspect_vaf.R")
source("R/03_annotate_variants.R")
source("R/04_hotspot_summary.R")
source("R/05_apply_hotspot_tiers.R")
source("R/06_link_to_depmap.R")
source("R/07_waterfall_plots.R")
# source("R/08_make_vep_input.R")              # optional: cohort VEP input
# After results/cohort_full/09_vep_full_alphamissense.vcf.gz exists:
source("R/09_am_dependency_figures.R")
# source("R/10_prism_mutation_sensitivity.R")   # optional: PRISM CSVs
```

### External inputs (not in this repo)

Paths are configured at the top of each script (default **`C:/Users/mvijayan/Documents`**).

| Input | Used by |
| --- | --- |
| Long-format VAF / marker table (`AllMarkers_VAF_long.tsv` or project-specific path) | `01`–`06` |
| `depmap_CGE.csv`, `depmap_meta.csv` | `06`, `07`, `09` |
| `samples.txt` (PID ↔ ACH mapping) | `06` |
| `results/cohort_full/09_vep_full_alphamissense.vcf.gz` | `09` |
| `PRISMOncologyReferenceSeqLog2AUCMatrix.csv`, `PRISMOncologyReferenceSeqCompoundList.csv` | `10` |

### Target oncogenes — single source of truth

`data/target_genes.tsv` is the **single source of truth** for the
gene set the rest of the pipeline scans. By default,
`R/00_build_target_genes.R` reads
`data/oncogene_shortlist_ped_gof_snvs.tsv` (**24 genes**, all tier **A**
in that file, with per-gene `score` / `evidence_tags`), then adds GRCh38
coordinates via Ensembl REST. The panel is tuned for **SNV/InDel/MNV
callsets** (no fusion partners): RAS/MAPK, PI3K–AKT–mTOR, RTKs, WNT,
selected oncohistones, **FLT3**, etc.

**24 genes (alphabetical):**  
`ACVR1, AKT1, ALK, BRAF, CTNNB1, EGFR, FGFR1, FGFR4, FLT3, HRAS, H3F3A, HIST1H3B, HIST1H3C, IDH1, KIT, KRAS, MTOR, MYC, MYCN, NRAS, PDGFRA, PIK3CA, PTPN11, SMO`

**Legacy 42-gene list:** `data/oncogene_shortlist_sjpedpanel.tsv` (Tier
A/B/C from SJPedPanel supplementary filtering — Karol et al. 2024) plus
`data/filter_oncogenes_sjpedpanel.R`. To rebuild that list after table
updates, run `Rscript data/filter_oncogenes_sjpedpanel.R`, point
`SHORTLIST` in `00_build_target_genes.R` at that TSV, then
`Rscript R/00_build_target_genes.R`.

Required R packages: `data.table`, `httr2`, `jsonlite`, `digest`,
`ggplot2`. Optional: `ggrepel` (volcano labels in `10`). Install with
`install.packages(c("data.table","httr2","jsonlite","digest","ggplot2"))`.

`R/03_annotate_variants.R` calls the public Ensembl VEP REST endpoint
by default and caches responses to `data/cache/`. If you have the local
VEP install (via `setup_alpha_tools.sh`), you can swap it in for offline
or large-batch runs.

All scripts write into `results/<RUN_NAME>/` (default `ped_gof_snv`).
The 10-gene pilot is preserved at `results/pilot_10genes/`. To re-run
with a different gene set, change `RUN_NAME` and/or the shortlist TSV
in `00_build_target_genes.R`, rebuild `target_genes.tsv`, then re-run
01..07.

See [`results/STATUS.md`](./results/STATUS.md) for the side-by-side
comparison of frozen runs and outstanding scientific questions.

---

## References

- Cheng J. *et al.* **Accurate proteome-wide missense variant effect
  prediction with AlphaMissense.** *Science* (2023).
  DOI: [10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492)
- AlphaMissense source code:
  <https://github.com/google-deepmind/alphamissense>
- AlphaMissense predictions bucket:
  <https://console.cloud.google.com/storage/browser/dm_alphamissense>
- Ensembl Variant Effect Predictor (VEP):
  <https://www.ensembl.org/info/docs/tools/vep/index.html>
- bioconda `ensembl-vep` recipe:
  <https://bioconda.github.io/recipes/ensembl-vep/README.html>
