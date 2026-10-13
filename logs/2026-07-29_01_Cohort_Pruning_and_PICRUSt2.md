# Removal of Batch GJ3, Five-Batch Rebuild and PICRUSt2

**Date:** 2026-07-29 (PICRUSt2 finished on 2026-07-30)
**Stage:** removal of batch GJ3; five-batch merge, filtering and tree; QC reports; core diversity metrics; PICRUSt2
**Cohort:** five batches (GB1–GB3, GJ1–GJ2), original labels; N = 110 per assay (`GB_Control` 33, `GB_Treatment` 33, `GJ_Control` 22, `GJ_Treatment` 22)
**Script:** `modules/10_cohort_pruning_and_picrust2.sh`; the same steps are in `pipeline/run_pipeline.sh`
**Platform:** `qiime2-amp` and `picrust2` conda environments (`osx-64`)
**Input:** `03_denoised_16S/*`, `03_denoised_ITS/*` (GB1–GB3, GJ1–GJ2), `05_taxonomy_16S/taxonomy_16S.qza`, `05_taxonomy_ITS/taxonomy_ITS.qza`
**Output:** `04_merged_*`, `05_taxonomy_*/table_*_clean.qza`, `06_phylo_16S/*`, `07_diversity_16S/core_metrics/*`, `07_diversity_ITS/core_metrics/*`, `09_qc_reports/*`, `picrust2_out_16S/*`

---

## 1. Why batch GJ3 was removed

The reassignment in 2026-07-24_02 produced the expected group sizes (34 and 34), but it has two problems:

* **Circularity.** The labels were chosen from PCoA positions. Testing the same groups with PERMANOVA or differential abundance on the same data reduces within-group variation by construction and inflates the test statistics.
* **Unresolved records.** Without harvest records that settle the labels, the PCoA cannot rule out other handling problems in set 4.

All GJ3 libraries (`J4`, `jM4` and the seven x/y samples; 24 per assay) were therefore removed from both assays. The final design is unbalanced:

| Experiment | Batches | Control | Treated |
| :--- | :--- | ---: | ---: |
| GB | GB1, GB2, GB3 | 33 | 33 |
| GJ | GJ1, GJ2 | 22 | 22 |

All analyses from 2026-08-10 onwards use this set with the original labels in `metadata_P_generation.tsv`.

---

## 2. Rebuild of the five-batch tables

**Step 1. Removal of GJ3 files and old outputs.**

> **Warning:** these commands permanently delete the raw GJ3 FASTQ folders and earlier outputs. Keep a backup before running module 10.

```bash
rm -rf Bacteria_GJ3 Fungi_GJ3
rm -rf 04_merged_16S 04_merged_ITS 07_diversity_16S 07_diversity_ITS 08_exported
rm -f 01_imported/*GJ3* 02_trimmed_16S/*GJ3* 02_trimmed_ITS/*GJ3* 03_denoised_16S/*GJ3* 03_denoised_ITS/*GJ3*
rm -f 05_taxonomy_16S/table_* 05_taxonomy_ITS/table_*
```

**Step 2. Merge, filter and tree.**

* The per-batch DADA2 outputs of GB1–GB3 and GJ1–GJ2 were merged as in 2026-07-23_01.
* The same taxonomic filters were applied: the 16S exclusion list (mitochondria, chloroplast, Archaea, Eukaryota, Cyanobacteria, Rickettsiales, Rickettsiellales, Unassigned) and `k__Fungi` only for ITS2. The taxonomy files from the six-batch classification were reused; they cover all ASVs of the five batches.
* The 16S tree was rebuilt from the five-batch representative sequences (`align-to-tree-mafft-fasttree`, `--p-n-threads 8`).

---

## 3. QC reports (`09_qc_reports/`)

1. **DADA2 read tracking:** `16S_dada2_filtering_stats.qzv` and `ITS_dada2_filtering_stats.qzv`, made from tables written by `merge_16S_dada2_stats.py` and `merge_ITS_dada2_stats.py`. These two helper scripts are not in the repository; module 10 and `pipeline/run_pipeline.sh` skip this step if they are missing.
2. **Depth summaries** of the filtered tables: `16S_sequencing_depth_summary.qzv`, `ITS_sequencing_depth_summary.qzv`.
3. **Alpha rarefaction curves:** 16S up to 27,000 reads (with the tree); ITS2 up to 20,000 reads.

```bash
qiime diversity alpha-rarefaction \
  --i-table 05_taxonomy_16S/table_16S_clean.qza \
  --i-phylogeny 06_phylo_16S/rooted_tree_16S.qza \
  --p-max-depth 27000 \
  --m-metadata-file metadata_P_generation.tsv \
  --o-visualization 09_qc_reports/16S_alpha_rarefaction_curves.qzv
```

The original log reports that the curves level off below the chosen rarefaction depths (16S at 27,000 reads; ITS2 by about 6,000–8,000 reads). The QC reports are not included in the repository.

---

## 4. Core diversity metrics

| | 16S | ITS2 |
| :--- | :--- | :--- |
| Command | `core-metrics-phylogenetic` | `core-metrics` |
| Input table | `05_taxonomy_16S/table_16S_clean.qza` | `05_taxonomy_ITS/table_ITS_clean.qza` |
| Tree | `06_phylo_16S/rooted_tree_16S.qza` | none |
| Rarefaction depth | 27,000 | 8,600 |
| Alpha metrics | observed features, Shannon, Pielou's evenness, Faith's PD | observed features, Shannon, Pielou's evenness |
| Distances | Bray–Curtis, Jaccard, weighted and unweighted UniFrac | Bray–Curtis, Jaccard |
| Output | `07_diversity_16S/core_metrics/` | `07_diversity_ITS/core_metrics/` |

All 110 16S libraries are above 27,000 reads (all 134 were retained at this depth in the six-batch run).

---

## 5. PICRUSt2

```bash
qiime tools export --input-path 04_merged_16S/repseq_16S_merged.qza --output-path picrust2_inputs
qiime tools export --input-path 05_taxonomy_16S/table_16S_clean.qza --output-path picrust2_inputs

conda activate picrust2
picrust2_pipeline.py \
  -s picrust2_inputs/dna-sequences.fasta \
  -i picrust2_inputs/feature-table.biom \
  -o picrust2_out_16S \
  -p 10 \
  --verbose
```

From the run output (total run time 3,575 s):

* **Input.** The FASTA contained all 29,299 merged 16S ASVs (before taxonomic filtering). The BIOM table is the taxonomy-filtered, unrarefied table (27,767 ASVs, all present in the FASTA).
* **Placement.** Sequences were placed on the bacterial and archaeal reference trees with EPA-ng (v0.3.8) and GAPPA (v0.9.0), with `--min_align 0.8`. 98 sequences (bacterial reference) and 289 sequences (archaeal reference) aligned poorly and were not placed. After choosing the better domain for each sequence, 29,103 sequences were kept for bacteria and 103 for archaea.
* **Hidden-state prediction.** 16S copy number, EC numbers and KEGG orthologs were predicted with castor using maximum parsimony (`--hsp_method mp`), and NSTI was calculated.
* **NSTI filter.** 171 of 29,206 ASVs had NSTI > 2.0 and were excluded from the metagenome predictions (`--max_nsti 2.0`). Predicted gene abundances were normalized by predicted 16S copy number.
* **Pathways.** MetaCyc pathway abundances were inferred from the predicted EC numbers with `pathway_pipeline.py` (map file `metacyc_pathways_structured_filtered_v24.txt`). The unstratified table `picrust2_out_16S/pathways_out/path_abun_unstrat.tsv.gz` was used for differential abundance (2026-08-24_01).

PICRUSt2 output is a prediction from 16S sequences and reference genomes, not measured gene content.
