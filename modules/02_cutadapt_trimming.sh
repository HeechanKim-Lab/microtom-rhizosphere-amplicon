#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 02_cutadapt_trimming.sh
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# Description: Demultiplexed paired-end primer and heterogeneity spacer excision
#              via Cutadapt with degenerate wildcard matching.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

THREADS=10
BATCHES=("GB1" "GB2" "GB3" "GJ1" "GJ2")

echo "=== [$(date +'%T')] Initializing Cutadapt Primer Trimming Module ==="

mkdir -p 02_trimmed_16S 02_trimmed_ITS

# ==============================================================================
# 1. Bacterial 16S rRNA (V3-V4) Primer & Spacer Trimming
# ==============================================================================
# 341F: CCTACGGGNGGCWGCAG (+ 0-3 bp N-spacers)
# 805R: GACTACHVGGGTATCTAATCC (+ 0-3 bp N-spacers)
for batch in "${BATCHES[@]}"; do
  echo "=================================================="
  echo "Trimming primers and N-spacers for Bacteria_${batch}..."
  echo "=================================================="

  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "01_imported/demux_16S_${batch}.qza" \
    --p-front-f CCTACGGGNGGCWGCAG \
    --p-front-r GACTACHVGGGTATCTAATCC \
    --p-error-rate 0.15 \
    --p-match-adapter-wildcards \
    --p-match-read-wildcards \
    --p-discard-untrimmed \
    --p-cores "${THREADS}" \
    --o-trimmed-sequences "02_trimmed_16S/trimmed_16S_${batch}.qza"

  # Generate diagnostic summary visualizer
  qiime demux summarize \
    --i-data "02_trimmed_16S/trimmed_16S_${batch}.qza" \
    --o-visualization "02_trimmed_16S/trimmed_16S_${batch}.qzv"
done

# ==============================================================================
# 2. Fungal ITS2 Primer Trimming
# ==============================================================================
# ITS3: GCATCGATGAAGAACGCAGC
# ITS4: TCCTCCGCTTATTGATATGC
for batch in "${BATCHES[@]}"; do
  echo "=================================================="
  echo "Trimming ITS primers for Fungi_${batch}..."
  echo "=================================================="

  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "01_imported/demux_ITS_${batch}.qza" \
    --p-front-f GCATCGATGAAGAACGCAGC \
    --p-front-r TCCTCCGCTTATTGATATGC \
    --p-error-rate 0.15 \
    --p-match-adapter-wildcards \
    --p-match-read-wildcards \
    --p-discard-untrimmed \
    --p-cores "${THREADS}" \
    --o-trimmed-sequences "02_trimmed_ITS/trimmed_ITS_${batch}.qza"

  # Generate diagnostic summary visualizer
  qiime demux summarize \
    --i-data "02_trimmed_ITS/trimmed_ITS_${batch}.qza" \
    --o-visualization "02_trimmed_ITS/trimmed_ITS_${batch}.qzv"
done

echo "=== [$(date +'%T')] Cutadapt Primer Trimming Complete. Outputs saved to 02_trimmed_16S/ and 02_trimmed_ITS/ ==="