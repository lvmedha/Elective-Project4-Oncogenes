#!/usr/bin/env bash
echo "--- python3 ---"
which python3
python3 --version
echo "--- pip3 ---"
which pip3 || echo "no pip3"
pip3 --version 2>&1 | head -3 || true
echo "--- python3 -m pip ---"
python3 -m pip --version 2>&1 | head -3 || true
echo "--- modules already installed ---"
python3 -c "
mods = ['openpyxl','pandas','xlrd','tablib']
for m in mods:
    try:
        __import__(m)
        print(m, 'OK')
    except Exception as e:
        print(m, 'MISSING', e)
"
echo "--- apt python3-openpyxl/pandas ---"
dpkg -l 2>/dev/null | awk '$2 ~ /^python3-(openpyxl|pandas|pip)$/ {print}'
