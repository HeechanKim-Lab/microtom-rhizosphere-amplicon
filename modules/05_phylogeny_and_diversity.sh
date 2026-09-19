#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 05_phylogeny_and_diversity.sh
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Sensitivity Testing (Full) & Primary Experimental Subsetting (Standard)
# Description: Rarefaction normalization, phylogenetic/non-phylogenetic alpha
#              and beta diversity metric generation, and PERMANOVA audits
#              evaluating Treatment Group contrasts and Sequencing Batch effects.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

DEPTH_16S=27000
DEPTH_ITS=8600
METADATA_FILE="metadata_P_generation.tsv"

echo "=== [$(date +'%T')] Initializing Diversity Metrics and PERMANOVA Audit Module ==="

# Verify metadata presence
if [[ ! -f "${METADATA_FILE}" ]]; then
  echo "[-] ERROR: Metadata file '${METADATA_FILE}' not found in working directory." >&2
  exit 1
fi

mkdir -p 07_diversity_16S/full 07_diversity_16S/standard 07_diversity_16S/permanova_audit \
         07_diversity_ITS/full 07_diversity_ITS/standard 07_diversity_ITS/batch_audit

# ==============================================================================
# 1. Core Diversity Metric Calculation: Bacterial 16S rRNA
# ==============================================================================
# Rarefaction Depth: 27,000 reads/sample (retains 100% of biological libraries)

#### A. FULL Dataset (Standard + Extra x/y Samples; Sensitivity Baseline)
echo "--- Calculating 16S Core Diversity Metrics (Full Dataset) ---"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --p-sampling-depth "${DEPTH_16S}" \
  --m-metadata-file "${METADATA_FILE}" \
  --output-dir 07_diversity_16S/full

#### B. STANDARD Dataset (Strictly Curated Baseline)
echo "--- Calculating 16S Core Diversity Metrics (Standard Dataset) ---"
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --i-table 05_taxonomy_16S/table_16S_standard.qza \
  --p-sampling-depth "${DEPTH_16S}" \
  --m-metadata-file "${METADATA_FILE}" \
  --output-dir 07_diversity_16S/standard

# ==============================================================================
# 2. Core Diversity Metric Calculation: Fungal ITS2
# ==============================================================================
# Rarefaction Depth: 8,600 reads/sample (Phylogeny omitted for non-coding ITS2)

#### A. FULL Dataset (With Extra x/y Samples)
echo "--- Calculating ITS Core Diversity Metrics (Full Dataset) ---"
qiime diversity core-metrics \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --p-sampling-depth "${DEPTH_ITS}" \
  --m-metadata-file "${METADATA_FILE}" \
  --output-dir 07_diversity_ITS/full

#### B. STANDARD Dataset (Without Extra x/y Samples)
echo "--- Calculating ITS Core Diversity Metrics (Standard Dataset) ---"
qiime diversity core-metrics \
  --i-table 05_taxonomy_ITS/table_ITS_standard.qza \
  --p-sampling-depth "${DEPTH_ITS}" \
  --m-metadata-file "${METADATA_FILE}" \
  --output-dir 07_diversity_ITS/standard

# ==============================================================================
# 3. PERMANOVA Audit: Primary Treatment Group Contrasts
# ==============================================================================

#### A. Bacteria (16S) Treatment Significance
echo "--- Running 16S Treatment PERMANOVA (Full & Standard) ---"
qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_16S/full/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Treatment_Group \
  --o-visualization 07_diversity_16S/permanova_audit/PERMANOVA_treatment_full.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_16S/standard/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Treatment_Group \
  --o-visualization 07_diversity_16S/permanova_audit/PERMANOVA_treatment_standard.qzv

#### B. Fungi (ITS) Treatment Significance
echo "--- Running ITS Treatment PERMANOVA (Full & Standard) ---"
qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_ITS/full/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Treatment_Group \
  --o-visualization 07_diversity_ITS/PERMANOVA_ITS_full.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_ITS/standard/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Treatment_Group \
  --o-visualization 07_diversity_ITS/PERMANOVA_ITS_standard.qzv

# ==============================================================================
# 4. PERMANOVA Audit: Batch Effect Technical Variance Evaluation
# ==============================================================================

#### A. Bacteria (16S) Batch Effect Significance
echo "--- Running 16S Batch Effect PERMANOVA (Full & Standard) ---"
qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_16S/full/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Sequencing_Batch \
  --o-visualization 07_diversity_16S/permanova_audit/PERMANOVA_batch_full.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_16S/standard/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Sequencing_Batch \
  --o-visualization 07_diversity_16S/permanova_audit/PERMANOVA_batch_standard.qzv

#### B. Fungi (ITS) Batch Effect Significance
echo "--- Running ITS Batch Effect PERMANOVA (Full & Standard) ---"
qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_ITS/full/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Sequencing_Batch \
  --o-visualization 07_diversity_ITS/batch_audit/ITS_batch_effect_full.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_ITS/standard/bray_curtis_distance_matrix.qza \
  --m-metadata-file "${METADATA_FILE}" \
  --m-metadata-column Sequencing_Batch \
  --o-visualization 07_diversity_ITS/batch_audit/ITS_batch_effect_standard.qzv

echo "=== [$(date +'%T')] Core Diversity Metrics and PERMANOVA Audits Complete ==="