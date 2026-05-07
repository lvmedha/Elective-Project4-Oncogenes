#!/usr/bin/env bash
set -euo pipefail
python3 -m pip install --quiet --user --break-system-packages openpyxl pandas
python3 .scratch/inspect_xlsx.py
