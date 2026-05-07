"""Inspect the supplementary Excel: list sheet names, sizes, and the first
few rows + headers of each sheet so we can understand the structure.
"""
import sys
import pandas as pd

XLSX = "/mnt/c/Users/mvijayan/Downloads/ccr-24-1063_supplementary_tables_1_suppts1.xlsx"

xl = pd.ExcelFile(XLSX, engine="openpyxl")
print("=== Sheets ===")
for s in xl.sheet_names:
    print(repr(s))

for s in xl.sheet_names:
    print("\n\n=========================================")
    print(f"Sheet: {s!r}")
    print("=========================================")
    # Read with no header so we see any title rows
    df = pd.read_excel(XLSX, sheet_name=s, header=None, engine="openpyxl")
    print(f"shape: {df.shape}")
    nshow = min(15, len(df))
    with pd.option_context("display.max_columns", 30, "display.width", 220, "display.max_colwidth", 60):
        print(df.head(nshow).to_string())
