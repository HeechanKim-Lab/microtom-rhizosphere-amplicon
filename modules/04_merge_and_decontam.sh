#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 04_merge_and_decontam.sh
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# Description: Feature table and sequence consolidation, 16S phylogenetic tree
#              reconstruction (MAFFT + FastTree), taxonomic classification,
#              host organelle/contaminant filtering, and metadata-driven
#              standard cohort stratification.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

THREADS=8

echo "=== [$(date +'%T')] Initializing Merging, Phylogeny, and Decontamination Module ==="

mkdir -p 04_merged_16S 04_merged_ITS \
         05_taxonomy_16S 05_taxonomy_ITS \
         06_phylo_16S

# ==============================================================================
# 1. Feature Table & Sequence Consolidation (Verified 5 Batches)
# ==============================================================================

#### A. 16S Bacteria Merging
echo "--- Merging 16S Feature Tables and Rep-Seqs ---"
qiime feature-table merge \
  --i-tables 03_denoised_16S/table_16S_GB1.qza \
             03_denoised_16S/table_16S_GB2.qza \
             03_denoised_16S/table_16S_GB3.qza \
             03_denoised_16S/table_16S_GJ1.qza \
             03_denoised_16S/table_16S_GJ2.qza \
  --o-merged-table 04_merged_16S/table_16S_merged.qza

qiime feature-table merge-seqs \
  --i-data 03_denoised_16S/repseq_16S_GB1.qza \
           03_denoised_16S/repseq_16S_GB2.qza \
           03_denoised_16S/repseq_16S_GB3.qza \
           03_denoised_16S/repseq_16S_GJ1.qza \
           03_denoised_16S/repseq_16S_GJ2.qza \
  --o-merged-data 04_merged_16S/repseq_16S_merged.qza

qiime feature-table summarize \
  --i-table 04_merged_16S/table_16S_merged.qza \
  --o-visualization 04_merged_16S/table_16S_merged.qzv

#### B. ITS Fungi Merging
echo "--- Merging ITS Feature Tables and Rep-Seqs ---"
qiime feature-table merge \
  --i-tables 03_denoised_ITS/table_ITS_GB1.qza \
             03_denoised_ITS/table_ITS_GB2.qza \
             03_denoised_ITS/table_ITS_GB3.qza \
             03_denoised_ITS/table_ITS_GJ1.qza \
             03_denoised_ITS/table_ITS_GJ2.qza \
  --o-merged-table 04_merged_ITS/table_ITS_merged.qza

qiime feature-table merge-seqs \
  --i-data 03_denoised_ITS/repseq_ITS_GB1.qza \
           03_denoised_ITS/repseq_ITS_GB2.qza \
           03_denoised_ITS/repseq_ITS_GB3.qza \
           03_denoised_ITS/repseq_ITS_GJ1.qza \
           03_denoised_ITS/repseq_ITS_GJ2.qza \
  --o-merged-data 04_merged_ITS/repseq_ITS_merged.qza

qiime feature-table summarize \
  --i-table 04_merged_ITS/table_ITS_merged.qza \
  --o-visualization 04_merged_ITS/table_ITS_merged.qzv

# ==============================================================================
# 2. 16S Multiple Sequence Alignment & FastTree Phylogenetic Inference
# ==============================================================================
echo "--- Constructing 16S Rooted Phylogenetic Tree ---"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 04_merged_16S/repseq_16S_merged.qza \
  --p-n-threads "${THREADS}" \
  --o-alignment 06_phylo_16S/aligned_repseq_16S.qza \
  --o-masked-alignment 06_phylo_16S/masked_aligned_repseq_16S.qza \
  --o-tree 06_phylo_16S/unrooted_tree_16S.qza \
  --o-rooted-tree 06_phylo_16S/rooted_tree_16S.qza

# ==============================================================================
# 3. Taxonomic Classification & Decontamination Filtering
# ==============================================================================

#### A. 16S Bacteria (SILVA 138 & Host Organelle Negative Exclusion)
echo "--- Classifying and Filtering 16S Bacteria ---"
qiime feature-classifier classify-sklearn \
  --i-classifier classifiers/silva-138-99-nb-classifier.qza \
  --i-reads 04_merged_16S/repseq_16S_merged.qza \
  --p-n-jobs "${THREADS}" \
  --o-classification 05_taxonomy_16S/taxonomy_16S.qza

qiime taxa filter-table \
  --i-table 04_merged_16S/table_16S_merged.qza \
  --i-taxonomy 05_taxonomy_16S/taxonomy_16S.qza \
  --p-exclude mitochondria,chloroplast,Archaea,Eukaryota,Cyanobacteria,Rickettsiales,Rickettsiellales,Unassigned \
  --o-filtered-table 05_taxonomy_16S/table_16S_clean.qza

qiime feature-table summarize \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --o-visualization 05_taxonomy_16S/table_16S_clean.qzv

#### B. ITS Fungi (UNITE v10 & Strict k__Fungi Positive Inclusion)
echo "--- Classifying and Filtering ITS Fungi ---"
qiime feature-classifier classify-sklearn \
  --i-classifier classifiers/unite_ver10_99_04.04.2024-Q2-2024.5.qza \
  --i-reads 04_merged_ITS/repseq_ITS_merged.qza \
  --p-n-jobs "${THREADS}" \
  --o-classification 05_taxonomy_ITS/taxonomy_ITS.qza

qiime taxa filter-table \
  --i-table 04_merged_ITS/table_ITS_merged.qza \
  --i-taxonomy 05_taxonomy_ITS/taxonomy_ITS.qza \
  --p-include k__Fungi \
  --o-filtered-table 05_taxonomy_ITS/table_ITS_clean.qza

qiime feature-table summarize \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --o-visualization 05_taxonomy_ITS/table_ITS_clean.qzv

# ==============================================================================
# 4. Generate Metadata Mapping File
# ==============================================================================
echo "--- Generating Sample Metadata ---"
if [[ -f "modules/make_metadata.py" ]]; then
  python modules/make_metadata.py
elif [[ -f "make_metadata.py" ]]; then
  python make_metadata.py
fi

# ==============================================================================
# 5. Metadata-Driven Cohort Stratification (Standard Only)
# ==============================================================================
if [[ -f "metadata_P_generation.tsv" ]]; then
  echo "--- Subsetting Standard Sample Cohorts ---"

  #### A. Filter 16S Bacteria (Standard Only)
  qiime feature-table filter-samples \
    --i-table 05_taxonomy_16S/table_16S_clean.qza \
    --m-metadata-file metadata_P_generation.tsv \
    --p-where "[Sample_Type]='Standard'" \
    --o-filtered-table 05_taxonomy_16S/table_16S_standard.qza

  qiime feature-table summarize \
    --i-table 05_taxonomy_16S/table_16S_standard.qza \
    --o-visualization 05_taxonomy_16S/table_16S_standard.qzv

  #### B. Filter ITS Fungi (Standard Only)
  qiime feature-table filter-samples \
    --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
    --m-metadata-file metadata_P_generation.tsv \
    --p-where "[Sample_Type]='Standard'" \
    --o-filtered-table 05_taxonomy_ITS/table_ITS_standard.qza

  qiime feature-table summarize \
    --i-table 05_taxonomy_ITS/table_ITS_standard.qza \
    --o-visualization 05_taxonomy_ITS/table_ITS_standard.qzv
else
  echo "[!] Warning: metadata_P_generation.tsv not found. Skipping standard subsetting."
fi

echo "=== [$(date +'%T')] Merging, Phylogeny, and Filtering Complete ==="