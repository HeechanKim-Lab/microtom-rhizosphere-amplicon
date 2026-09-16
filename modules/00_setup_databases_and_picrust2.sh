#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 00_download_classifiers.sh
# Target: Reference Classifiers & PICRUSt2 (v2.6.3) Environment Setup
# Platform: macOS (osx-64 via Rosetta 2 translation)
# Description: Downloads pre-trained Naive Bayes classifiers for SILVA 138 &
#              UNITE v10, builds the picrust2 conda prefix, and patches default
#              reference database assets.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

echo "=== [$(date +'%T')] Initializing Database and Classifier Setup Module ==="

# ==============================================================================
# 1. Reference Taxonomic Classifiers (SILVA 138 & UNITE v10)
# ==============================================================================
mkdir -p classifiers

# 1.1 SILVA 138 99% Naive Bayes Classifier (16S Bacteria)
if [[ ! -f "classifiers/silva-138-99-nb-classifier.qza" ]]; then
  echo "[+] Downloading SILVA 138 99% Naive Bayes Classifier..."
  wget -O classifiers/silva-138-99-nb-classifier.qza \
    "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
else
  echo "[*] SILVA 138 classifier already exists. Skipping download."
fi

# 1.2 UNITE v10 99% Naive Bayes Classifier (ITS Fungi)
if [[ ! -f "classifiers/unite_ver10_99_04.04.2024-Q2-2024.5.qza" ]]; then
  echo "[+] Downloading UNITE v10 99% Naive Bayes Classifier..."
  wget -O classifiers/unite_ver10_99_04.04.2024-Q2-2024.5.qza \
    "https://github.com/colinbrislawn/unite-train/releases/download/v10.0-v04.04.2024-qiime2-2024.2/unite_ver10_99_04.04.2024-Q2-2024.5.qza"
else
  echo "[*] UNITE v10 classifier already exists. Skipping download."
fi

# ==============================================================================
# 2. PICRUSt2 (v2.6.3) Environment Setup & Asset Patching
# ==============================================================================
# Set up dedicated x86_64 Conda prefix for PICRUSt2 if not already present
if ! conda info --envs | grep -q "picrust2"; then
  echo "[+] Creating dedicated picrust2 conda environment (CONDA_SUBDIR=osx-64)..."
  CONDA_SUBDIR=osx-64 conda create -n picrust2 \
    -c bioconda \
    -c conda-forge \
    python=3.9 \
    r-base r-castor \
    hmmer epa-ng pplacer gappa glpk sepp \
    biom-format ete3 dendropy h5py jinja2 joblib numpy pandas scipy statsmodels typeguard wget \
    -y

  echo "[+] Configuring osx-64 pinning for picrust2 prefix..."
  conda config --env --set subdir osx-64

  # Temporarily activate to execute pip installation and default file patching
  eval "$(conda shell.bash hook)"
  conda activate picrust2

  echo "[+] Installing PICRUSt2 v2.6.3 core via pip..."
  pip install https://github.com/picrust/picrust2/archive/v2.6.3.tar.gz

  echo "[+] Fetching and patching default reference database assets..."
  git clone -b v2.6.3 --depth 1 https://github.com/picrust/picrust2.git /tmp/picrust2_src
  cp -r /tmp/picrust2_src/picrust2/default_files "$CONDA_PREFIX/lib/python3.9/site-packages/picrust2/"
  rm -rf /tmp/picrust2_src

  echo "[+] Verifying critical binary dependencies in system PATH..."
  which hmmalign epa-ng pplacer gappa Rscript

  conda deactivate
else
  echo "[*] picrust2 environment already detected. Skipping creation."
fi

echo "=== [$(date +'%T')] Classifier Downloads and PICRUSt2 Setup Complete ==="