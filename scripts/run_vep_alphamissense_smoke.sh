#!/usr/bin/env bash
# Smoke test: same VEP + AlphaMissense setup as run_vep_alphamissense.sh, but
# only the first N variant rows from the cohort VCF (full ## / #CHROM headers
# preserved). Use before the full ~4k-variant run.
# VEP cache + AlphaMissense file: defaults are repo-local (see run_vep_alphamissense.sh).
#
#   cd /mnt/c/Users/mvijayan/Documents/Elective-Project4-Oncogenes
#   bash scripts/run_vep_alphamissense_smoke.sh
#
#   SMOKE_N=100 VEP_FORKS=2 bash scripts/run_vep_alphamissense_smoke.sh
#
# Defaults to VEP_FORKS=1 for the delegated run (fewer issues with tabix-backed
# plugins); set VEP_FORKS explicitly to parallelize.
set -euo pipefail

echo "run_vep_alphamissense_smoke.sh: starting" >&2

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${SMOKE_N:=25}"
: "${SOURCE_VCF:=${PROJECT_ROOT}/results/cohort_full/08_vep_input_full.vcf}"

if ! [[ "${SMOKE_N}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: SMOKE_N must be a positive integer, got: ${SMOKE_N}" >&2
  exit 1
fi

if [[ ! -f "$SOURCE_VCF" ]]; then
  echo "ERROR: source VCF not found: $SOURCE_VCF" >&2
  exit 1
fi

SMOKE_DIR="${PROJECT_ROOT}/results/cohort_full/smoke"
mkdir -p "$SMOKE_DIR"

SMOKE_VCF="${SMOKE_DIR}/08_vep_input_first${SMOKE_N}_variants.vcf"
# Do not use `grep … | head` here: with `set -o pipefail` grep exits 141 (SIGPIPE)
# when head stops reading, and the script can exit with almost no message.
awk -v n="$SMOKE_N" '
  /^##/      { print; next }
  /^#CHROM/  { print; next }
  /^#/       { next }
  { if (count < n) { print; count++ } }
' "$SOURCE_VCF" > "$SMOKE_VCF"

n_var="$(grep -cv '^#' "$SMOKE_VCF" || true)"
if [[ "$n_var" -lt 1 ]]; then
  echo "ERROR: smoke VCF has no variant rows (check SOURCE_VCF / SMOKE_N)." >&2
  exit 1
fi

echo "Smoke test: first ${SMOKE_N} variant line(s) requested; ${n_var} row(s) in ${SMOKE_VCF}"
echo

# Single-process VEP avoids rare fork + tabix-plugin issues; override with VEP_FORKS=4 ...
: "${VEP_FORKS:=1}"

IN_VCF="$SMOKE_VCF" \
OUT_VCF="${SMOKE_DIR}/09_vep_smoke_n${SMOKE_N}_alphamissense.vcf.gz" \
STATS_FILE="${SMOKE_DIR}/09_vep_smoke_n${SMOKE_N}_stats.txt" \
WARNING_FILE="${SMOKE_DIR}/09_vep_smoke_n${SMOKE_N}_warnings.txt" \
VEP_FORKS="$VEP_FORKS" \
bash "${SCRIPT_DIR}/run_vep_alphamissense.sh"
