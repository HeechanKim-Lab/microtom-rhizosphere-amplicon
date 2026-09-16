#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 01_remediation_and_import.sh
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# Description: Conda environment initialization, library duplication handling,
#              directory scaffolding, and automated Casava 1.8 batch ingestion.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

# ==============================================================================
# 0. Environment Setup (Run once if environment is not already present)
# ==============================================================================
# CONDA_NO_PLUGINS=true CONDA_SUBDIR=osx-64 conda env create --solver=classic \
#   -n qiime2-amp --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
# conda activate qiime2-amp
# conda config --env --set subdir osx-64

echo "=== [$(date +'%T')] Initializing Remediation and Import Module ==="
echo "Active QIIME 2 Binary: $(which qiime)"

# ==============================================================================
# 1. Directory Scaffolding
# ==============================================================================
mkdir -p 01_imported 02_trimmed 03_denoised 04_merged 05_taxonomy 06_phylo 07_diversity 08_exported

# ==============================================================================
# 2. Re-Sequencing Concatenation & Failed Run Quarantine
# ==============================================================================

#### Handling Duplication of Bacteria_GJ2 (B_J3-4: S7 + S79)
if [[ -f "Bacteria_GJ2/B_J3-4_S7_L001_R1_001.fastq.gz" && -f "Bacteria_GJ2/B_J3-4_S79_L001_R1_001.fastq.gz" ]]; then
  echo "[+] Concatenating re-sequenced pairs for Bacteria_GJ2: B_J3-4 (S7 + S79)..."
  mkdir -p Bacteria_GJ2/raw_runs
  
  cat Bacteria_GJ2/B_J3-4_S7_L001_R1_001.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R1_001.fastq.gz > Bacteria_GJ2/B_J3-4_merged_R1.fastq.gz
  cat Bacteria_GJ2/B_J3-4_S7_L001_R2_001.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R2_001.fastq.gz > Bacteria_GJ2/B_J3-4_merged_R2.fastq.gz
  
  mv Bacteria_GJ2/B_J3-4_S7_* Bacteria_GJ2/B_J3-4_S79_* Bacteria_GJ2/raw_runs/
  mv Bacteria_GJ2/B_J3-4_merged_R1.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R1_001.fastq.gz
  mv Bacteria_GJ2/B_J3-4_merged_R2.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R2_001.fastq.gz
else
  echo "[*] Bacteria_GJ2 B_J3-4 already remediated or raw files archived. Skipping."
fi

#### Handling Duplication of Bacteria_GB3 (B_M4-1: S45 Quarantined, S78 + S96 Merged)
if [[ -f "Bacteria_GB3/B_M4-1_S45_L001_R1_001.fastq.gz" || -f "Bacteria_GB3/B_M4-1_S78_L001_R1_001.fastq.gz" ]]; then
  echo "[+] Quarantining failed run S45 and concatenating S78 + S96 for Bacteria_GB3: B_M4-1..."
  mkdir -p Bacteria_GB3/raw_runs
  
  # Archive failed under-sequenced run
  if [[ -f "Bacteria_GB3/B_M4-1_S45_L001_R1_001.fastq.gz" ]]; then
    mv Bacteria_GB3/B_M4-1_S45_* Bacteria_GB3/raw_runs/
  fi
  
  # Concatenate supplementary runs
  if [[ -f "Bacteria_GB3/B_M4-1_S78_L001_R1_001.fastq.gz" && -f "Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz" ]]; then
    cat Bacteria_GB3/B_M4-1_S78_L001_R1_001.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz > Bacteria_GB3/B_M4-1_merged_R1.fastq.gz
    cat Bacteria_GB3/B_M4-1_S78_L001_R2_001.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R2_001.fastq.gz > Bacteria_GB3/B_M4-1_merged_R2.fastq.gz
    
    mv Bacteria_GB3/B_M4-1_S78_* Bacteria_GB3/B_M4-1_S96_* Bacteria_GB3/raw_runs/
    mv Bacteria_GB3/B_M4-1_merged_R1.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz
    mv Bacteria_GB3/B_M4-1_merged_R2.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R2_001.fastq.gz
  fi
else
  echo "[*] Bacteria_GB3 B_M4-1 already remediated or raw files archived. Skipping."
fi

# ==============================================================================
# 3. Automated Batch Ingestion into QIIME 2 Artifacts
# ==============================================================================
# Cohort locked to 5 verified batches (GJ3 omitted per July 29 cohort audit)
BATCHES=("GB1" "GB2" "GB3" "GJ1" "GJ2")

#### Import 16S Bacterial Batches
for batch in "${BATCHES[@]}"; do
  echo "--- Ingesting Bacteria_${batch} ---"
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "Bacteria_${batch}" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "01_imported/demux_16S_${batch}.qza"
done

#### Import ITS Fungal Batches
for batch in "${BATCHES[@]}"; do
  echo "--- Ingesting Fungi_${batch} ---"
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "Fungi_${batch}" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "01_imported/demux_ITS_${batch}.qza"
done

echo "=== [$(date +'%T')] Batch Ingestion Complete. Outputs saved to 01_imported/ ==="