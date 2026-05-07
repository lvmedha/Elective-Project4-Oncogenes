#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# setup_env.sh
#
# One-shot installer for the AlphaMissense environment used in this project.
# Tested on Ubuntu 22.04 / WSL2-Ubuntu. Should work on any Debian-based distro.
#
# Note: this only sets up Google DeepMind's AlphaMissense (a python venv at
# alphamissense/venv). For the SEPARATE Ensembl-VEP / `alpha_tools` conda env
# used as a local fallback by R/03_annotate_variants.R, see setup_alpha_tools.sh
# (or pass --all here to chain both).
#
# Usage:
#   bash setup_env.sh                # AlphaMissense env only
#   bash setup_env.sh --with-data    # also download precomputed AlphaMissense
#                                    # predictions (~5-10 GB, into ./data)
#   bash setup_env.sh --all          # also run setup_alpha_tools.sh after
#   bash setup_env.sh --help
#
# After it finishes, activate the AlphaMissense env with:
#   source alphamissense/venv/bin/activate
# ---------------------------------------------------------------------------

set -euo pipefail

# ----- arg parsing ---------------------------------------------------------
WITH_DATA=0
ALSO_ALPHA_TOOLS=0
for arg in "$@"; do
    case "$arg" in
        --with-data) WITH_DATA=1 ;;
        --all)       ALSO_ALPHA_TOOLS=1 ;;
        -h|--help)
            sed -n '2,23p' "$0"
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
log()  { printf '\n\033[1;34m[setup]\033[0m %s\n' "$*"; }
warn() { printf '\n\033[1;33m[warn]\033[0m  %s\n' "$*"; }
die()  { printf '\n\033[1;31m[error]\033[0m %s\n' "$*" >&2; exit 1; }

# ----- sanity checks -------------------------------------------------------
if [[ "$(uname -s)" != "Linux" ]]; then
    die "AlphaMissense requires Linux. On Windows, run this script inside WSL2 (Ubuntu)."
fi

if ! command -v sudo >/dev/null 2>&1; then
    die "'sudo' is required to install system packages. Please install sudo first."
fi

# ----- 1. system packages --------------------------------------------------
log "Installing system packages (python3.11-venv, aria2, hmmer, git)..."
sudo apt-get update -y
sudo apt-get install -y python3.11-venv aria2 hmmer git

# ----- 2. clone AlphaMissense ---------------------------------------------
if [[ ! -d alphamissense ]]; then
    log "Cloning google-deepmind/alphamissense..."
    git clone https://github.com/google-deepmind/alphamissense.git
else
    log "alphamissense/ already exists - skipping clone."
fi

cd alphamissense

# ----- 3. python venv ------------------------------------------------------
if [[ ! -d venv ]]; then
    log "Creating Python virtual environment in alphamissense/venv ..."
    python3 -m venv ./venv
else
    log "venv/ already exists - reusing it."
fi

log "Upgrading pip and installing AlphaMissense + dependencies..."
./venv/bin/pip install --upgrade pip
./venv/bin/pip install -r requirements.txt
./venv/bin/pip install -e .

# ----- 4. installation test -----------------------------------------------
log "Running the AlphaMissense installation test..."
if ./venv/bin/python test/test_installation.py; then
    log "Installation test passed."
else
    warn "Installation test failed. Check the output above before continuing."
fi

cd ..

# ----- 5. (optional) precomputed predictions ------------------------------
if [[ $WITH_DATA -eq 1 ]]; then
    log "Downloading precomputed AlphaMissense predictions into ./data ..."
    mkdir -p data

    BASE_URL="https://storage.googleapis.com/dm_alphamissense"
    FILES=(
        "AlphaMissense_aa_substitutions.tsv.gz"
        "AlphaMissense_hg38.tsv.gz"
    )

    for f in "${FILES[@]}"; do
        if [[ -f "data/$f" ]]; then
            log "  data/$f already present - skipping."
            continue
        fi
        log "  downloading $f ..."
        aria2c -x 8 -s 8 -d data -o "$f" "$BASE_URL/$f"
    done
fi

# ----- 6. (optional) chain into setup_alpha_tools.sh ----------------------
if [[ $ALSO_ALPHA_TOOLS -eq 1 ]]; then
    if [[ -x ./setup_alpha_tools.sh || -f ./setup_alpha_tools.sh ]]; then
        log "Running setup_alpha_tools.sh ..."
        bash ./setup_alpha_tools.sh
    else
        warn "--all was passed but setup_alpha_tools.sh was not found in $(pwd). Skipping."
    fi
fi

# ----- done ----------------------------------------------------------------
cat <<EOF

\033[1;32mAll done.\033[0m

To start using the AlphaMissense environment:

    cd alphamissense
    source venv/bin/activate

To deactivate when finished:

    deactivate

To also build the VEP / alpha_tools conda env (one-time):

    bash setup_alpha_tools.sh
EOF
