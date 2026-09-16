#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 03_dada2_denoise.sh
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# Description: Within-run parametric error rate learning, quality filtering,
#              dereplication, paired-end merging, and bimera removal via DADA2.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

THREADS=10
BATCHES=("GB1" "GB2" "GB3" "GJ1" "GJ2")

echo "=== [$(date +'%T')] Initializing DADA2 Denoising Module ==="

mkdir -p 03_denoised_16S 03_denoised_ITS

# ==============================================================================
# 1. Bacterial 16S DADA2 Denoising (Positional Truncation)
# ==============================================================================
# Parameters: Trunc F=260, R=200; maxEE F=2.0, R=3.0 (retains >20 bp overlap for 460 bp amplicon)
for batch in "${BATCHES[@]}"; do
  echo "=================================================="
  echo "Running DADA2 for Bacteria_${batch}..."
  echo "=================================================="

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "02_trimmed_16S/trimmed_16S_${batch}.qza" \
    --p-trunc-len-f 260 \
    --p-trunc-len-r 200 \
    --p-max-ee-f 2.0 \
    --p-max-ee-r 3.0 \
    --p-n-threads "${THREADS}" \
    --o-table "03_denoised_16S/table_16S_${batch}.qza" \
    --o-representative-sequences "03_denoised_16S/repseq_16S_${batch}.qza" \
    --o-denoising-stats "03_denoised_16S/stats_16S_${batch}.qza"

  # Tabulate denoising stats visualizer
  qiime metadata tabulate \
    --m-input-file "03_denoised_16S/stats_16S_${batch}.qza" \
    --o-visualization "03_denoised_16S/stats_16S_${batch}.qzv"
done

# ==============================================================================
# 2. Fungal ITS2 DADA2 Denoising (Zero Truncation)
# ==============================================================================
# Parameters: Trunc F=0, R=0 (preserves biological length polymorphism); maxEE F=2.0, R=3.0
for batch in "${BATCHES[@]}"; do
  echo "=================================================="
  echo "Running DADA2 for Fungi_${batch}..."
  echo "=================================================="

  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "02_trimmed_ITS/trimmed_ITS_${batch}.qza" \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-max-ee-f 2.0 \
    --p-max-ee-r 3.0 \
    --p-n-threads "${THREADS}" \
    --o-table "03_denoised_ITS/table_ITS_${batch}.qza" \
    --o-representative-sequences "03_denoised_ITS/repseq_ITS_${batch}.qza" \
    --o-denoising-stats "03_denoised_ITS/stats_ITS_${batch}.qza"

  # Tabulate denoising stats visualizer
  qiime metadata tabulate \
    --m-input-file "03_denoised_ITS/stats_ITS_${batch}.qza" \
    --o-visualization "03_denoised_ITS/stats_ITS_${batch}.qzv"
done

echo "=== [$(date +'%T')] DADA2 Denoising Complete. Outputs saved to 03_denoised_16S/ and 03_denoised_ITS/ ==="