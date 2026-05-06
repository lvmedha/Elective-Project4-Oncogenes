# Elective-Project4-Oncogenes

Course project analyzing oncogene missense variants using
[AlphaMissense](https://github.com/google-deepmind/alphamissense).

This README documents the steps needed to reproduce the AlphaMissense
environment used in this project.

---

## Environment Setup (AlphaMissense)

AlphaMissense is built for **Linux**. On Windows, use **WSL2 (Ubuntu)**.
On macOS or native Linux, skip the WSL step.

### 1. (Windows only) Install WSL2 + Ubuntu

In an **Administrator** PowerShell:

```powershell
wsl --install
```

Reboot when prompted, then launch **Ubuntu** from the Start menu and
create your Linux username and password. All remaining steps are run
**inside the WSL/Ubuntu shell**.

### 2. Run the automated setup script (recommended)

The repo includes [`setup_env.sh`](./setup_env.sh) which performs every
install step for you (system packages, cloning AlphaMissense, creating
the Python venv, installing dependencies, and running the install test).

From the repo root:

```bash
bash setup_env.sh
```

To also download the precomputed AlphaMissense predictions
(~5–10 GB, into `./data`):

```bash
bash setup_env.sh --with-data
```

When the script finishes, activate the environment with:

```bash
cd alphamissense
source venv/bin/activate
```

The script is idempotent — re-running it is safe and will skip steps
that are already complete.

### 3. Manual setup (fallback)

If you prefer to run the steps yourself, or the script fails on your
system, the equivalent manual commands are:

```bash
sudo apt update
sudo apt install -y python3.11-venv aria2 hmmer git

git clone https://github.com/google-deepmind/alphamissense.git
cd alphamissense

python3 -m venv ./venv
venv/bin/pip install --upgrade pip
venv/bin/pip install -r requirements.txt
venv/bin/pip install -e .

venv/bin/python test/test_installation.py
```

`hmmer` provides `jackhmmer`, which AlphaMissense uses to build the
multiple sequence alignments. `aria2` is used for fast database
downloads.

---

## (Optional) Download precomputed predictions

DeepMind publishes precomputed AlphaMissense scores for every possible
human missense substitution in a public Google Cloud Storage bucket.
For most downstream oncogene analyses you only need these files — you
do **not** need to run the model yourself.

The easiest way is `bash setup_env.sh --with-data` (see above), which
downloads the two most commonly used files via `aria2`.

To download manually instead, the files are at:

```text
https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
```

You can browse the full bucket here:
<https://console.cloud.google.com/storage/browser/dm_alphamissense>

| File | Contents |
| --- | --- |
| `AlphaMissense_aa_substitutions.tsv.gz` | Pathogenicity scores for all possible single-amino-acid substitutions in the human proteome. |
| `AlphaMissense_hg38.tsv.gz` | Same predictions mapped to the GRCh38 (hg38) genome coordinates. |

---

## (Optional) Genetic databases

If you intend to *run* the AlphaMissense data pipeline (not just use
the precomputed predictions), you will additionally need the genetic
sequence databases (BFD, MGnify, UniRef90). Follow the instructions in
the [AlphaFold repository](https://github.com/google-deepmind/alphafold)
to download them. These are large (hundreds of GB).

---

## References

- Cheng J. *et al.* **Accurate proteome-wide missense variant effect
  prediction with AlphaMissense.** *Science* (2023).
  DOI: [10.1126/science.adg7492](https://doi.org/10.1126/science.adg7492)
- AlphaMissense source code:
  <https://github.com/google-deepmind/alphamissense>
- AlphaMissense predictions bucket:
  <https://console.cloud.google.com/storage/browser/dm_alphamissense>
