# Elective-Project4-Oncogenes

Course project analyzing oncogene missense variants using
[AlphaMissense](https://github.com/google-deepmind/alphamissense).

This README documents the steps needed to reproduce the AlphaMissense
environment used in this project.

---

## Environment Setup (AlphaMissense)

AlphaMissense is built for **Linux**. On Windows, use **WSL2 (Ubuntu)**.
On macOS or native Linux, skip step 1.

### 1. (Windows only) Install WSL2 + Ubuntu

In an **Administrator** PowerShell:

```powershell
wsl --install
```

Reboot when prompted, then launch **Ubuntu** from the Start menu and
create your Linux username and password. All remaining steps are run
**inside the WSL/Ubuntu shell**.

### 2. Install system dependencies

```bash
sudo apt update
sudo apt install -y python3.11-venv aria2 hmmer git
```

`hmmer` provides `jackhmmer`, which AlphaMissense uses to build the
multiple sequence alignments. `aria2` is used for fast database
downloads.

### 3. Clone AlphaMissense

```bash
git clone https://github.com/google-deepmind/alphamissense.git
cd alphamissense
```

### 4. Create the Python virtual environment

```bash
python3 -m venv ./venv
venv/bin/pip install --upgrade pip
venv/bin/pip install -r requirements.txt
venv/bin/pip install -e .
```

### 5. Verify the installation

```bash
venv/bin/python test/test_installation.py
```

If the test passes, the environment is ready.

### 6. Reactivating the environment in a new shell

```bash
cd alphamissense
source venv/bin/activate
```

---

## (Optional) Download precomputed predictions

DeepMind publishes precomputed AlphaMissense scores for every possible
human missense substitution in a public Google Cloud Storage bucket.
For most downstream oncogene analyses you only need these files — you
do **not** need to run the model yourself.

Browse the bucket here:
<https://console.cloud.google.com/storage/browser/dm_alphamissense>

To download the two most commonly used files (~several GB each) using
`gsutil` (install instructions: <https://cloud.google.com/storage/docs/gsutil_install>):

```bash
mkdir -p data
gsutil -m cp gs://dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz ./data/
gsutil -m cp gs://dm_alphamissense/AlphaMissense_hg38.tsv.gz             ./data/
```

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
