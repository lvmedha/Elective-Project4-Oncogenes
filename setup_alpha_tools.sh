#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_alpha_tools.sh
#
# Reproducible setup for the `alpha_tools` conda env used by this project.
# Provides a LOCAL Ensembl VEP install (via bioconda) as an optional fallback
# for R/03_annotate_variants.R, which by default uses the Ensembl REST API.
#
# What this script does:
#   1. Make sure a mamba/conda runtime exists (installs Miniforge if not).
#   2. Create the `alpha_tools` env from environment.yml
#      (or environment.lock.yml with --lock for exact reproducibility).
#   3. Optionally download the VEP GRCh38 cache (~25 GB) with --with-vep-cache.
#
# Usage:
#   bash setup_alpha_tools.sh                  # create env from environment.yml
#   bash setup_alpha_tools.sh --lock           # use environment.lock.yml (pinned)
#   bash setup_alpha_tools.sh --with-vep-cache # also fetch GRCh38 VEP cache
#   bash setup_alpha_tools.sh --force          # recreate env even if it exists
#   bash setup_alpha_tools.sh --help
#
# After it finishes, activate the env with:
#   conda activate alpha_tools
# (you may need: source ~/mambaforge/etc/profile.d/conda.sh   first)
# ---------------------------------------------------------------------------

set -euo pipefail

# ----- arg parsing ---------------------------------------------------------
USE_LOCK=0
WITH_VEP_CACHE=0
FORCE=0
for arg in "$@"; do
    case "$arg" in
        --lock)            USE_LOCK=1 ;;
        --with-vep-cache)  WITH_VEP_CACHE=1 ;;
        --force)           FORCE=1 ;;
        -h|--help)
            sed -n '2,22p' "$0"
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Run with --help for usage." >&2
            exit 1
            ;;
    esac
done

# ----- helpers -------------------------------------------------------------
log()  { printf '\n\033[1;34m[alpha_tools]\033[0m %s\n' "$*"; }
warn() { printf '\n\033[1;33m[warn]\033[0m        %s\n' "$*"; }
die()  { printf '\n\033[1;31m[error]\033[0m       %s\n' "$*" >&2; exit 1; }

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_NAME="alpha_tools"
SPEC_FILE="environment.yml"
[[ $USE_LOCK -eq 1 ]] && SPEC_FILE="environment.lock.yml"

# ----- sanity checks -------------------------------------------------------
if [[ "$(uname -s)" != "Linux" ]]; then
    die "This script targets Linux (Ensembl VEP / bioconda). On Windows, run inside WSL2."
fi

if [[ ! -f "$PROJECT_DIR/$SPEC_FILE" ]]; then
    die "Could not find $SPEC_FILE in $PROJECT_DIR. Run this from the project root."
fi

# ----- 1. ensure conda/mamba is available ---------------------------------
# Order of preference: existing $CONDA_EXE > mamba > conda > install Miniforge.
CONDA_BIN=""
MAMBA_BIN=""

if command -v mamba >/dev/null 2>&1; then
    MAMBA_BIN="$(command -v mamba)"
fi
if command -v conda >/dev/null 2>&1; then
    CONDA_BIN="$(command -v conda)"
fi

# Pick up the user's existing mambaforge / miniforge / miniconda if they
# haven't been initialized in this shell (common with `bash -lc`).
if [[ -z "$CONDA_BIN" ]]; then
    for cand in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/anaconda3" "$HOME/miniconda3" /opt/conda; do
        if [[ -x "$cand/bin/conda" ]]; then
            log "Found conda at $cand (not yet initialized in this shell)."
            # shellcheck disable=SC1091
            source "$cand/etc/profile.d/conda.sh"
            CONDA_BIN="$cand/bin/conda"
            [[ -x "$cand/bin/mamba" ]] && MAMBA_BIN="$cand/bin/mamba"
            break
        fi
    done
fi

if [[ -z "$CONDA_BIN" ]]; then
    log "No conda/mamba found. Installing Miniforge to ~/miniforge3 ..."
    INSTALLER="/tmp/miniforge.sh"
    URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
    if command -v curl >/dev/null 2>&1; then
        curl -L "$URL" -o "$INSTALLER"
    elif command -v wget >/dev/null 2>&1; then
        wget -O "$INSTALLER" "$URL"
    else
        die "Need curl or wget to download Miniforge."
    fi
    bash "$INSTALLER" -b -p "$HOME/miniforge3"
    rm -f "$INSTALLER"
    # shellcheck disable=SC1091
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    CONDA_BIN="$HOME/miniforge3/bin/conda"
    MAMBA_BIN="$HOME/miniforge3/bin/mamba"
    log "Initializing conda for bash (one-time)."
    "$CONDA_BIN" init bash >/dev/null
fi

log "Using conda: $CONDA_BIN"
[[ -n "$MAMBA_BIN" ]] && log "Using mamba: $MAMBA_BIN"

# Prefer mamba for env creation (much faster solver); fall back to conda.
SOLVER="$CONDA_BIN"
[[ -n "$MAMBA_BIN" ]] && SOLVER="$MAMBA_BIN"

# ----- 2. create / update the env -----------------------------------------
ENV_EXISTS=0
if "$CONDA_BIN" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    ENV_EXISTS=1
fi

if [[ $ENV_EXISTS -eq 1 && $FORCE -eq 0 ]]; then
    log "Env '$ENV_NAME' already exists. Skipping create. (Pass --force to recreate.)"
elif [[ $ENV_EXISTS -eq 1 && $FORCE -eq 1 ]]; then
    log "Removing existing env '$ENV_NAME' (--force) ..."
    "$CONDA_BIN" env remove -n "$ENV_NAME" -y
    log "Creating '$ENV_NAME' from $SPEC_FILE ..."
    "$SOLVER" env create -f "$PROJECT_DIR/$SPEC_FILE"
else
    log "Creating '$ENV_NAME' from $SPEC_FILE ..."
    "$SOLVER" env create -f "$PROJECT_DIR/$SPEC_FILE"
fi

# ----- 3. quick sanity check ----------------------------------------------
log "Verifying VEP is on PATH inside the env ..."
"$CONDA_BIN" run -n "$ENV_NAME" bash -c '
    set -e
    which vep
    vep --help | head -3
    which vep_install
'

# ----- 4. (optional) VEP cache --------------------------------------------
if [[ $WITH_VEP_CACHE -eq 1 ]]; then
    CACHE_DIR="$PROJECT_DIR/data/vep"
    mkdir -p "$CACHE_DIR"
    log "Downloading + indexing VEP GRCh38 cache into $CACHE_DIR (this is ~25 GB and slow) ..."
    "$CONDA_BIN" run -n "$ENV_NAME" \
        vep_install -a cf -s homo_sapiens -y GRCh38 -c "$CACHE_DIR" --CONVERT --NO_HTSLIB --NO_UPDATE
    log "VEP cache ready at $CACHE_DIR"
fi

# ----- done ----------------------------------------------------------------
cat <<EOF

\033[1;32m[alpha_tools] All done.\033[0m

To activate:

    source ~/mambaforge/etc/profile.d/conda.sh   # or your conda dir
    conda activate alpha_tools

To rebuild from the exact lockfile:

    bash setup_alpha_tools.sh --lock --force

To also download the VEP GRCh38 cache (~25 GB):

    bash setup_alpha_tools.sh --with-vep-cache
EOF
