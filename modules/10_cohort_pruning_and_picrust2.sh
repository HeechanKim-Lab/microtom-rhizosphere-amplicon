#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 10_cohort_pruning_and_picrust2.sh
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# Platform: macOS (osx-64 via Rosetta 2)
# Description: Purges legacy/ambiguous GJ3 artifacts, merges the 5 verified
#              batches, executes strict taxonomic decontamination, reconstructs
#              16S phylogeny, computes core diversity metrics, and runs PICRUSt2.
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

THREADS=8
DEPTH_16S=27000
DEPTH_ITS=8600
METADATA_FILE="metadata_P_generation.tsv"

echo "=== [$(date +'%T')] Initializing Cohort Pruning, Re-Analysis, and PICRUSt2 Pipeline ==="

# ==============================================================================
# 1. Purge Legacy/GJ3 Files & Re-Structure Workspace
# ==============================================================================
echo "--- Step 1: Purging GJ3 Legacy Artifacts & Rebuilding Hierarchy ---"

# 1.1 Remove raw, imported, and denoised GJ3 files
rm -rf Bacteria_GJ3 Fungi_GJ3
rm -rf 04_merged_16S 04_merged_ITS 07_diversity_16S 07_diversity_ITS 08_exported
rm -f 01_imported/*GJ3* 02_trimmed_16S/*GJ3* 02_trimmed_ITS/*GJ3* 03_denoised_16S/*GJ3* 03_denoised_ITS/*GJ3*
rm -f 05_taxonomy_16S/table_* 05_taxonomy_ITS/table_*

# 1.2 Create fresh directory structure
mkdir -p 04_merged_16S 04_merged_ITS \
         05_taxonomy_16S 05_taxonomy_ITS \
         06_phylo_16S \
         07_diversity_16S/core_metrics 07_diversity_ITS/core_metrics \
         09_qc_reports \
         picrust2_inputs

# ==============================================================================
# 2. Merge Feature Tables & Representative Sequences (Verified 5 Batches)
# ==============================================================================
echo "--- Step 2: Merging Verified 5-Batch Datasets (N=110) ---"

#### A. 16S Bacteria Merging
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

#### B. ITS Fungi Merging
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

# ==============================================================================
# 3. Strict Taxonomic Contaminant Filtering
# ==============================================================================
echo "--- Step 3: Taxonomic Decontamination Filtering ---"

# Filter 16S Bacteria (Negative exclusion against SILVA 138)
qiime taxa filter-table \
  --i-table 04_merged_16S/table_16S_merged.qza \
  --i-taxonomy 05_taxonomy_16S/taxonomy_16S.qza \
  --p-exclude mitochondria,chloroplast,Archaea,Eukaryota,Cyanobacteria,Rickettsiales,Rickettsiellales,Unassigned \
  --o-filtered-table 05_taxonomy_16S/table_16S_clean.qza

# Filter ITS Fungi (Positive inclusion strictly within Kingdom Fungi)
qiime taxa filter-table \
  --i-table 04_merged_ITS/table_ITS_merged.qza \
  --i-taxonomy 05_taxonomy_ITS/taxonomy_ITS.qza \
  --p-include k__Fungi \
  --o-filtered-table 05_taxonomy_ITS/table_ITS_clean.qza

# ==============================================================================
# 4. Reconstruct 16S Phylogenetic Tree (MAFFT + FastTree)
# ==============================================================================
echo "--- Step 4: De Novo Phylogenetic Alignment & Reconstruction ---"

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 04_merged_16S/repseq_16S_merged.qza \
  --p-n-threads "${THREADS}" \
  --o-alignment 06_phylo_16S/aligned_repseq_16S.qza \
  --o-masked-alignment 06_phylo_16S/masked_aligned_repseq_16S.qza \
  --o-tree 06_phylo_16S/unrooted_tree_16S.qza \
  --o-rooted-tree 06_phylo_16S/rooted_tree_16S.qza

# ==============================================================================
# 5. Generate Senior QC Deliverables & Tracking Visualizers
# ==============================================================================
echo "--- Step 5: Generating Senior Quality Control Visualizers ---"

# 5.1 DADA2 Read Tracking Visualizations
if [[ -f "modules/merge_16S_dada2_stats.py" ]]; then
  python modules/merge_16S_dada2_stats.py
elif [[ -f "merge_16S_dada2_stats.py" ]]; then
  python merge_16S_dada2_stats.py
fi

if [[ -f "09_qc_reports/merged_16S_dada2_stats.tsv" ]]; then
  qiime metadata tabulate \
    --m-input-file 09_qc_reports/merged_16S_dada2_stats.tsv \
    --o-visualization 09_qc_reports/16S_dada2_filtering_stats.qzv
fi

if [[ -f "modules/merge_ITS_dada2_stats.py" ]]; then
  python modules/merge_ITS_dada2_stats.py
elif [[ -f "merge_ITS_dada2_stats.py" ]]; then
  python merge_ITS_dada2_stats.py
fi

if [[ -f "09_qc_reports/merged_ITS_dada2_stats.tsv" ]]; then
  qiime metadata tabulate \
    --m-input-file 09_qc_reports/merged_ITS_dada2_stats.tsv \
    --o-visualization 09_qc_reports/ITS_dada2_filtering_stats.qzv
fi

# 5.2 Post-Filtering Sequencing Depth Summaries
qiime feature-table summarize \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --o-visualization 09_qc_reports/16S_sequencing_depth_summary.qzv

qiime feature-table summarize \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --o-visualization 09_qc_reports/ITS_sequencing_depth_summary.qzv

# 5.3 Alpha Rarefaction Saturation Curves
qiime diversity alpha-rarefaction \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --p-max-depth "${DEPTH_16S}" \
  --m-metadata-file "${METADATA_FILE}" \
  --o-visualization 09_qc_reports/16S_alpha_rarefaction_curves.qzv

qiime diversity alpha-rarefaction \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --p-max-depth 20000 \
  --m-metadata-file "${METADATA_FILE}" \
  --o-visualization 09_qc_reports/ITS_alpha_rarefaction_curves.qzv

# ==============================================================================
# 6. Core Diversity Metric Calculations (Locked N=110 Cohort)
# ==============================================================================
echo "--- Step 6: Calculating 16S and ITS Core Diversity Metrics ---"

# 6.1 16S Phylogenetic Core Diversity (Depth = 27,000)
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --p-sampling-depth "${DEPTH_16S}" \
  --m-metadata-file "${METADATA_FILE}" \
  --output-dir 07_diversity_16S/core_metrics

# 6.2 ITS Core Diversity (Depth = 8,600)
qiime diversity core-metrics \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --p-sampling-depth "${DEPTH_ITS}" \
  --m-metadata-file "${METADATA_FILE}" \
  --output-dir 07_diversity_ITS/core_metrics

# ==============================================================================
# 7. PICRUSt2 Functional Metagenome Inference
# ==============================================================================
echo "--- Step 7: Exporting Artifacts and Executing PICRUSt2 ---"

# 7.1 Export FASTA Sequences and BIOM Feature Table
qiime tools export \
  --input-path 04_merged_16S/repseq_16S_merged.qza \
  --output-path picrust2_inputs

qiime tools export \
  --input-path 05_taxonomy_16S/table_16S_clean.qza \
  --output-path picrust2_inputs

# 7.2 Run PICRUSt2 Full Pipeline
# Note: Ensure the 'picrust2' conda prefix is active if installed separately from qiime2-amp.
if command -v picrust2_pipeline.py &> /dev/null; then
  picrust2_pipeline.py \
    -s picrust2_inputs/dna-sequences.fasta \
    -i picrust2_inputs/feature-table.biom \
    -o picrust2_out_16S \
    -p "${THREADS}" \
    --verbose
else
  echo "[!] WARNING: 'picrust2_pipeline.py' not found in current PATH."
  echo "    Activate the dedicated environment: 'conda activate picrust2'"
  echo "    Then re-run Step 7."
fi

echo "=== [$(date +'%T')] Cohort Pruning, Diversity Calculations, and PICRUSt2 Execution Complete ==="