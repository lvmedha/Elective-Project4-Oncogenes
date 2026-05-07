#!/usr/bin/env bash
# Helper to run a python script inside the alpha_tools conda env from WSL.
# Usage: bash run_in_alpha_tools.sh <script.py> [extra args...]
set -euo pipefail

for cand in "$HOME/mambaforge/etc/profile.d/conda.sh" \
            "$HOME/miniforge3/etc/profile.d/conda.sh" \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh"; do
    if [ -f "$cand" ]; then
        # shellcheck disable=SC1090
        source "$cand"
        break
    fi
done

conda activate alpha_tools
exec python "$@"
