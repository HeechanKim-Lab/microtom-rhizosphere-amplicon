#!/usr/bin/env bash
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Target Assays: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# ==============================================================================

set -euo pipefail
IFS=$'\n\t'

# Initialize Conda for subshell scripting
eval "$(conda shell.bash hook)"

################# Environment 1: QIIME 2 (qiime2-amp) #################
echo "Activating qiime2-amp environment..."
conda activate qiime2-amp
conda config --env --set subdir osx-64
which qiime



################# Directory Setup & Automated QIIME 2 Import #################

#### Handling Duplication of Bacteria_GJ2
# 1. Create a backup folder for unmerged files
mkdir -p Bacteria_GJ2/raw_runs
# 2. Concatenate forward (R1) and reverse (R2) reads
cat Bacteria_GJ2/B_J3-4_S7_L001_R1_001.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R1_001.fastq.gz > Bacteria_GJ2/B_J3-4_merged_R1.fastq.gz
cat Bacteria_GJ2/B_J3-4_S7_L001_R2_001.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R2_001.fastq.gz > Bacteria_GJ2/B_J3-4_merged_R2.fastq.gz
# 3. Move originals to backup and rename merged pair to standard Casava format
mv Bacteria_GJ2/B_J3-4_S7_* Bacteria_GJ2/B_J3-4_S79_* Bacteria_GJ2/raw_runs/
mv Bacteria_GJ2/B_J3-4_merged_R1.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R1_001.fastq.gz
mv Bacteria_GJ2/B_J3-4_merged_R2.fastq.gz Bacteria_GJ2/B_J3-4_S79_L001_R2_001.fastq.gz

#### Handling Duplication of Bacteria_GB3
# 1. Create backup directory for unmerged files
mkdir -p Bacteria_GB3/raw_runs
# 2. Move failed run S45 to backup
mv Bacteria_GB3/B_M4-1_S45_* Bacteria_GB3/raw_runs/
# 3. Concatenate forward (R1) and reverse (R2) reads for S78 and S96
cat Bacteria_GB3/B_M4-1_S78_L001_R1_001.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz > Bacteria_GB3/B_M4-1_merged_R1.fastq.gz
cat Bacteria_GB3/B_M4-1_S78_L001_R2_001.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R2_001.fastq.gz > Bacteria_GB3/B_M4-1_merged_R2.fastq.gz
# 4. Move unmerged S78 & S96 to backup
mv Bacteria_GB3/B_M4-1_S78_* Bacteria_GB3/B_M4-1_S96_* Bacteria_GB3/raw_runs/
# 5. Rename merged files to standard Casava format
mv Bacteria_GB3/B_M4-1_merged_R1.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz
mv Bacteria_GB3/B_M4-1_merged_R2.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R2_001.fastq.gz

#### Import 16S Bacterial Batches (5 Batches)
mkdir -p 01_imported 02_trimmed_16S 02_trimmed_ITS 03_denoised_16S 03_denoised_ITS \
         04_merged_16S 04_merged_ITS 05_taxonomy_16S 05_taxonomy_ITS 06_phylo_16S \
         07_diversity_16S 07_diversity_ITS 09_qc_reports picrust2_inputs

for batch in GB1 GB2 GB3 GJ1 GJ2; do
  echo "Importing Bacteria_${batch}..."
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "Bacteria_${batch}" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "01_imported/demux_16S_${batch}.qza"
done

#### Import ITS Fungal Batches (5 Batches)
for batch in GB1 GB2 GB3 GJ1 GJ2; do
  echo "Importing Fungi_${batch}..."
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "Fungi_${batch}" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "01_imported/demux_ITS_${batch}.qza"
done



################# 16S Bacterial Primer Trimming #################
mkdir -p 02_trimmed_16S

for batch in GB1 GB2 GB3 GJ1 GJ2; do
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
    --p-cores 10 \
    --o-trimmed-sequences "02_trimmed_16S/trimmed_16S_${batch}.qza"
    
  qiime demux summarize \
    --i-data "02_trimmed_16S/trimmed_16S_${batch}.qza" \
    --o-visualization "02_trimmed_16S/trimmed_16S_${batch}.qzv"
done



################# ITS Fungal Primer Trimming #################
mkdir -p 02_trimmed_ITS

for batch in GB1 GB2 GB3 GJ1 GJ2; do
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
    --p-cores 10 \
    --o-trimmed-sequences "02_trimmed_ITS/trimmed_ITS_${batch}.qza"
  
  qiime demux summarize \
    --i-data "02_trimmed_ITS/trimmed_ITS_${batch}.qza" \
    --o-visualization "02_trimmed_ITS/trimmed_ITS_${batch}.qzv"
done



################# DADA2 Denoising #################
#### Bacteria (16S)
mkdir -p 03_denoised_16S

for batch in GB1 GB2 GB3 GJ1 GJ2; do
  echo "=================================================="
  echo "Running DADA2 for Bacteria_${batch}..."
  echo "=================================================="
  
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "02_trimmed_16S/trimmed_16S_${batch}.qza" \
    --p-trunc-len-f 260 \
    --p-trunc-len-r 200 \
    --p-max-ee-f 2.0 \
    --p-max-ee-r 3.0 \
    --p-n-threads 10 \
    --o-table "03_denoised_16S/table_16S_${batch}.qza" \
    --o-representative-sequences "03_denoised_16S/repseq_16S_${batch}.qza" \
    --o-denoising-stats "03_denoised_16S/stats_16S_${batch}.qza"

  qiime metadata tabulate \
    --m-input-file "03_denoised_16S/stats_16S_${batch}.qza" \
    --o-visualization "03_denoised_16S/stats_16S_${batch}.qzv"
done

#### Fungi (ITS)
mkdir -p 03_denoised_ITS

for batch in GB1 GB2 GB3 GJ1 GJ2; do
  echo "=================================================="
  echo "Running DADA2 for Fungi_${batch}..."
  echo "=================================================="
  
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs "02_trimmed_ITS/trimmed_ITS_${batch}.qza" \
    --p-trunc-len-f 0 \
    --p-trunc-len-r 0 \
    --p-max-ee-f 2.0 \
    --p-max-ee-r 3.0 \
    --p-n-threads 8 \
    --o-table "03_denoised_ITS/table_ITS_${batch}.qza" \
    --o-representative-sequences "03_denoised_ITS/repseq_ITS_${batch}.qza" \
    --o-denoising-stats "03_denoised_ITS/stats_ITS_${batch}.qza"

  qiime metadata tabulate \
    --m-input-file "03_denoised_ITS/stats_ITS_${batch}.qza" \
    --o-visualization "03_denoised_ITS/stats_ITS_${batch}.qzv"
done



################# Merge Feature Tables & Representative Sequences (5 Batches) #################
# --- 1. 16S Bacteria Merging ---
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

# --- 2. ITS Fungi Merging ---
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



################# Downloading Classifiers #################
mkdir -p classifiers

# 1. SILVA 138 99% Naive Bayes Classifier (16S Bacteria)
if [[ ! -f "classifiers/silva-138-99-nb-classifier.qza" ]]; then
  wget -O classifiers/silva-138-99-nb-classifier.qza \
    "https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza"
fi

# 2. UNITE v10 99% Naive Bayes Classifier (ITS Fungi)
if [[ ! -f "classifiers/unite_ver10_99_04.04.2024-Q2-2024.5.qza" ]]; then
  wget -O classifiers/unite_ver10_99_04.04.2024-Q2-2024.5.qza \
    "https://github.com/colinbrislawn/unite-train/releases/download/v10.0-v04.04.2024-qiime2-2024.5/unite_ver10_99_04.04.2024-Q2-2024.5.qza"
fi



################# Taxonomic Classification #################
#### Bacteria (16S)
qiime feature-classifier classify-sklearn \
  --i-classifier classifiers/silva-138-99-nb-classifier.qza \
  --i-reads 04_merged_16S/repseq_16S_merged.qza \
  --p-n-jobs 8 \
  --o-classification 05_taxonomy_16S/taxonomy_16S.qza

#### Fungi (ITS)
qiime feature-classifier classify-sklearn \
  --i-classifier classifiers/unite_ver10_99_04.04.2024-Q2-2024.5.qza \
  --i-reads 04_merged_ITS/repseq_ITS_merged.qza \
  --p-n-jobs 4 \
  --o-classification 05_taxonomy_ITS/taxonomy_ITS.qza



################# Strict Taxonomic Contaminant Filtering #################
# Filter 16S Bacteria
qiime taxa filter-table \
  --i-table 04_merged_16S/table_16S_merged.qza \
  --i-taxonomy 05_taxonomy_16S/taxonomy_16S.qza \
  --p-exclude mitochondria,chloroplast,Archaea,Eukaryota,Cyanobacteria,Rickettsiales,Rickettsiellales,Unassigned \
  --o-filtered-table 05_taxonomy_16S/table_16S_clean.qza

# Filter ITS Fungi
qiime taxa filter-table \
  --i-table 04_merged_ITS/table_ITS_merged.qza \
  --i-taxonomy 05_taxonomy_ITS/taxonomy_ITS.qza \
  --p-include k__Fungi \
  --o-filtered-table 05_taxonomy_ITS/table_ITS_clean.qza



################# Generate Metadata Mapping File #################
python modules/make_metadata.py



################# Build 16S Tree & Generate Senior QC Deliverables #################
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences 04_merged_16S/repseq_16S_merged.qza \
  --p-n-threads 8 \
  --o-alignment 06_phylo_16S/aligned_repseq_16S.qza \
  --o-masked-alignment 06_phylo_16S/masked_aligned_repseq_16S.qza \
  --o-tree 06_phylo_16S/unrooted_tree_16S.qza \
  --o-rooted-tree 06_phylo_16S/rooted_tree_16S.qza



################# Generate Requested QC Reports #################
# 1. DADA2 Read Tracking (Raw -> Filtered -> Denoised -> Merged -> Non-Chimeric)
if [[ -f "merge_16S_dada2_stats.py" ]]; then
  python merge_16S_dada2_stats.py
  qiime metadata tabulate \
    --m-input-file 09_qc_reports/merged_16S_dada2_stats.tsv \
    --o-visualization 09_qc_reports/16S_dada2_filtering_stats.qzv
fi

if [[ -f "merge_ITS_dada2_stats.py" ]]; then
  python merge_ITS_dada2_stats.py
  qiime metadata tabulate \
    --m-input-file 09_qc_reports/merged_ITS_dada2_stats.tsv \
    --o-visualization 09_qc_reports/ITS_dada2_filtering_stats.qzv
fi

# 2. Sequencing Depth Summaries (Post-Filtering)
qiime feature-table summarize \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --o-visualization 09_qc_reports/16S_sequencing_depth_summary.qzv

qiime feature-table summarize \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --o-visualization 09_qc_reports/ITS_sequencing_depth_summary.qzv

# 3. Per-Sample Alpha Rarefaction Curves
qiime diversity alpha-rarefaction \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --p-max-depth 27000 \
  --m-metadata-file metadata_P_generation.tsv \
  --o-visualization 09_qc_reports/16S_alpha_rarefaction_curves.qzv

qiime diversity alpha-rarefaction \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --p-max-depth 20000 \
  --m-metadata-file metadata_P_generation.tsv \
  --o-visualization 09_qc_reports/ITS_alpha_rarefaction_curves.qzv



################# Core Diversity Calculations #################
# 16S Bacteria Diversity
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --p-sampling-depth 27000 \
  --m-metadata-file metadata_P_generation.tsv \
  --output-dir 07_diversity_16S/core_metrics

# ITS Fungi Diversity
qiime diversity core-metrics \
  --i-table 05_taxonomy_ITS/table_ITS_clean.qza \
  --p-sampling-depth 8600 \
  --m-metadata-file metadata_P_generation.tsv \
  --output-dir 07_diversity_ITS/core_metrics



################# PERMANOVA Comparison #################
#### Bacteria (16S)
mkdir -p 07_diversity_16S/permanova_audit

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_16S/core_metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_P_generation.tsv \
  --m-metadata-column Treatment_Group \
  --o-visualization 07_diversity_16S/permanova_audit/PERMANOVA_treatment.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_16S/core_metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_P_generation.tsv \
  --m-metadata-column Sequencing_Batch \
  --o-visualization 07_diversity_16S/permanova_audit/PERMANOVA_batch.qzv

#### Fungi (ITS)
mkdir -p 07_diversity_ITS/permanova_audit

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_ITS/core_metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_P_generation.tsv \
  --m-metadata-column Treatment_Group \
  --o-visualization 07_diversity_ITS/permanova_audit/PERMANOVA_treatment.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix 07_diversity_ITS/core_metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata_P_generation.tsv \
  --m-metadata-column Sequencing_Batch \
  --o-visualization 07_diversity_ITS/permanova_audit/PERMANOVA_batch.qzv



################# PICRUSt2 Functional Prediction #################
# Export artifacts from qiime2-amp environment
qiime tools export \
  --input-path 04_merged_16S/repseq_16S_merged.qza \
  --output-path picrust2_inputs

qiime tools export \
  --input-path 05_taxonomy_16S/table_16S_clean.qza \
  --output-path picrust2_inputs

################# Environment 2: PICRUSt2 (picrust2) #################
echo "Switching to picrust2 environment..."
conda activate picrust2
which picrust2_pipeline.py

picrust2_pipeline.py \
  -s picrust2_inputs/dna-sequences.fasta \
  -i picrust2_inputs/feature-table.biom \
  -o picrust2_out_16S \
  -p 10 \
  --verbose

echo "Pipeline execution finished successfully."