#!/usr/bin/env bash
# Local WSL2: run Ensembl VEP (conda env alpha_tools) on the cohort VCF from
# R/08_make_vep_input.R with the AlphaMissense plugin and a GRCh38 cache.
#
# Defaults expect a local mirror under the repo (same layout as MastersStudy/Data):
#   ./vep_cache/homo_sapiens/112_GRCh38/...
#   ./AlphaMissense_hg38.tsv.bgz  (+ .tbi)
# Override with VEP_CACHE_DIR / AM_TSV if yours live elsewhere (e.g. /mnt/z/...).
#
#   cd /mnt/c/Users/mvijayan/Documents/Elective-Project4-Oncogenes
#   bash scripts/run_vep_alphamissense.sh
#
#   VEP_FORKS=8 OUT_VCF=/tmp/test.vcf.gz bash scripts/run_vep_alphamissense.sh
#
# Quick check on first N variants: scripts/run_vep_alphamissense_smoke.sh
#
set -euo pipefail

echo "run_vep_alphamissense.sh: starting" >&2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# --- user-tunable defaults ---------------------------------------------------
: "${VEP_CACHE_DIR:=${PROJECT_ROOT}/vep_cache}"
: "${CACHE_VERSION:=112}"
: "${ASSEMBLY:=GRCh38}"
: "${SPECIES:=homo_sapiens}"

: "${AM_TSV:=${PROJECT_ROOT}/AlphaMissense_hg38.tsv.bgz}"
# If the table lives under ./data/ instead: AM_TSV="${PROJECT_ROOT}/data/AlphaMissense_hg38.tsv.bgz" ...
: "${CHR_SYNONYMS:=${VEP_CACHE_DIR}/homo_sapiens/${CACHE_VERSION}_${ASSEMBLY}/chr_synonyms.txt}"

: "${IN_VCF:=${PROJECT_ROOT}/results/cohort_full/08_vep_input_full.vcf}"
: "${OUT_VCF:=${PROJECT_ROOT}/results/cohort_full/09_vep_full_alphamissense.vcf.gz}"

: "${VEP_FORKS:=4}"
: "${STATS_FILE:=${PROJECT_ROOT}/results/cohort_full/09_vep_full_alphamissense_stats.txt}"
: "${WARNING_FILE:=${PROJECT_ROOT}/results/cohort_full/09_vep_full_alphamissense_warnings.txt}"

# Optional: reference FASTA (not required for most cache-only runs). Uncomment
# in the vep invocation below if VEP asks for --fasta.
# : "${REF_FASTA:=${VEP_CACHE_DIR}/homo_sapiens/115_${ASSEMBLY}/Homo_sapiens.GRCh38.dna.toplevel.fa.gz}"

# --- conda (same search order as setup_alpha_tools / .scratch helpers) --------
CONDA_SH=""
for cand in "${HOME}/mambaforge/etc/profile.d/conda.sh" \
            "${HOME}/miniforge3/etc/profile.d/conda.sh" \
            "${HOME}/miniconda3/etc/profile.d/conda.sh" \
            "${HOME}/anaconda3/etc/profile.d/conda.sh"; do
  if [[ -f "$cand" ]]; then
    CONDA_SH="$cand"
    break
  fi
done
if [[ -z "$CONDA_SH" ]]; then
  echo "ERROR: conda.sh not found under ~/mambaforge, ~/miniforge3, ~/miniconda3, or ~/anaconda3" >&2
  exit 1
fi
# shellcheck disable=SC1090
source "$CONDA_SH"
if ! conda activate alpha_tools 2>/dev/null; then
  echo "ERROR: conda env 'alpha_tools' not found. Create it with: bash setup_alpha_tools.sh" >&2
  exit 1
fi

command -v vep >/dev/null || { echo "ERROR: vep not on PATH inside alpha_tools" >&2; exit 1; }

# --- sanity checks ------------------------------------------------------------
for p in "$IN_VCF" "$AM_TSV" "$CHR_SYNONYMS"; do
  if [[ ! -f "$p" ]]; then
    echo "ERROR: missing file: $p" >&2
    exit 1
  fi
done
if [[ ! -d "${VEP_CACHE_DIR}/homo_sapiens/${CACHE_VERSION}_${ASSEMBLY}" ]]; then
  echo "ERROR: cache dir not found: ${VEP_CACHE_DIR}/homo_sapiens/${CACHE_VERSION}_${ASSEMBLY}" >&2
  exit 1
fi

command -v tabix >/dev/null || {
  echo "ERROR: tabix not on PATH (htslib). The AlphaMissense VEP plugin calls tabix on the TSV." >&2
  exit 1
}
if ! tabix -H "$AM_TSV" >/dev/null 2>&1; then
  echo "ERROR: tabix cannot read AlphaMissense file (missing or incompatible .tbi next to the .bgz/.gz?):"
  echo "       $AM_TSV"
  echo "       Expected index: ${AM_TSV}.tbi"
  echo "       If you need to build it, see AlphaMissense.pm (tabix -s 1 -b 2 -e 2 -f -S 1 ...)." >&2
  exit 1
fi

mkdir -p "$(dirname "$OUT_VCF")" "$(dirname "$STATS_FILE")" "$(dirname "$WARNING_FILE")"

# Absolute path for plugin (recommended in AlphaMissense.pm)
AM_TSV_ABS="$(readlink -f "$AM_TSV" 2>/dev/null || echo "$AM_TSV")"

echo "VEP cache dir     : $VEP_CACHE_DIR"
echo "Cache version     : $CACHE_VERSION  assembly: $ASSEMBLY"
echo "Input VCF         : $IN_VCF"
echo "Output VCF        : $OUT_VCF"
echo "AlphaMissense TSV : $AM_TSV_ABS"
echo "synonyms file     : $CHR_SYNONYMS"
echo "Forks             : $VEP_FORKS"
echo
echo "Running VEP (no annotated VCF until this finishes; errors print below) ..."
echo

vep \
  --input_file "$IN_VCF" \
  --output_file "$OUT_VCF" \
  --fork "$VEP_FORKS" \
  --vcf \
  --compress gzip \
  --offline \
  --cache \
  --species "$SPECIES" \
  --assembly "$ASSEMBLY" \
  --dir "$VEP_CACHE_DIR" \
  --cache_version "$CACHE_VERSION" \
  --synonyms "$CHR_SYNONYMS" \
  --force_overwrite \
  --warning_file "$WARNING_FILE" \
  --stats_file "$STATS_FILE" \
  --plugin "AlphaMissense,file=${AM_TSV_ABS}"

echo
echo "Done. Main output: $OUT_VCF"
echo "Stats: $STATS_FILE | Warnings: $WARNING_FILE"
