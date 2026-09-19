# Methodological Audit: Cohort Pruning, Five-Batch Synthesis, and Predictive Functional Metagenomics via PICRUSt2
**Date:** July 29, 2026  
**Project:** 04_MicroTom_P_Rhizosphere_Amplicon_Analysis  
**Target Assays:** Bacterial 16S rRNA (V3–V4) & Fungal ITS2  
**Platform:** Emulated `osx-64` (`qiime2-amp` and `picrust2` conda prefixes)  
**Input Artifacts:** Batches `GB1`, `GB2`, `GB3`, `GJ1`, `GJ2` (`03_denoised_*`), `classifiers/*`  
**Output Artifacts:** `04_merged_*`, `05_taxonomy_*`, `06_phylo_16S/*`, `07_diversity_*`, `09_qc_reports/*`, `picrust2_out_16S/*`  

---

## 1. Methodological Justification: The Decision to Prune Cohort GJ3

While ordination-guided thresholding previously resolved the biological clusters of the mislabeled samples and produced an exact numerical match ($N=34$ for both `GJ_Control` and `GJ_Treatment`), post-hoc sample reassignment introduces critical methodological risks:

* **Risk of Circular Inference ("Double Dipping"):** Relabeling sample identifiers based on their positions in PCoA space, followed by hypothesis testing (PERMANOVA, differential abundance) on those same groups, risks circular reasoning. Optimizing sample classifications to maximize ordination separation artificially deflates within-group variance ($SS_W$), inflates pseudo-$F$ statistics, and invalidates family-wise error rate controls.
* **Persistent Physical Ambiguity:** Without immutable physical chain-of-custody documentation from the harvest bench, algorithmic recovery cannot rule out cross-contamination, non-standard microenvironmental effects, or handling artifacts in set 4.
* **Integrity-First Remediation:** To preserve strict experimental reproducibility, all libraries derived from physical batch `GJ3` (comprising `J4`, `jM4`, and the ambiguous $x/y$ samples) were purged from both bacterial and fungal workflows.

### Final Pruned Cohort Architecture
The dataset transitions to an unbalanced but strictly verified five-batch design ($N = 110$ total biological samples):
* **Gijang B Lineage ($N = 66$):** 3 independent batches (`GB1`, `GB2`, `GB3`) $\to 33$ Control (`GB_Control`), 33 Treatment (`GB_Treatment`).
* **Gyeongju Lineage ($N = 44$):** 2 independent batches (`GJ1`, `GJ2`) $\to 22$ Control (`GJ_Control`), 22 Treatment (`GJ_Treatment`).

---

## 2. Pipeline Execution: Purge, Re-merge, and Taxonomic Decontamination

All downstream directories and intermediate files derived from the compromised six-batch dataset were purged to prevent file caching artifacts. Feature tables and representative sequences were consolidated across the five verified batches.

```
                      Verified Raw Denoised Batches (5 Batches)
                      [GB1, GB2, GB3] x [GJ1, GJ2]
                                    │
       ┌────────────────────────────┴────────────────────────────┐
       ▼                                                         ▼
16S Bacterial Pipeline (N=110)                            ITS Fungal Pipeline (N=110)
  │                                                         │
  ├─▶ Merge Tables & Rep-Seqs                               ├─▶ Merge Tables & Rep-Seqs
  ├─▶ Strict Exclusion Filtering                            ├─▶ Strict Inclusion Filtering
  │   (Mito, Chloro, Archaea, Euk, Cyanobacteria)           │   (Strictly k__Fungi)
  ├─▶ Alignment & Phylogeny (MAFFT + FastTree)              ├─▶ [Skip Tree Building]
  ├─▶ Core Diversity (27,000 Depth)                         ├─▶ Core Diversity (8,600 Depth)
  └─▶ Export to PICRUSt2 Pipeline                           └─▶ Export Matrices to R
```

### A. Taxonomic Decontamination Enforced
* **16S rRNA:** Negative exclusion filtering against SILVA 138 annotations purged off-target mitochondrial DNA, chloroplast sequences, non-bacterial domains (Archaea, Eukaryota), photosynthetic surface contaminants (Cyanobacteria), and degenerate organellar lineages (Rickettsiales, Rickettsiellales).
* **ITS2:** Positive inclusion filtering retained features classified strictly within `k__Fungi` from UNITE v10, excluding co-amplified host plant nuclear rDNA and soil microeukaryotes.

---

## 3. Comprehensive Quality Control Audits (`09_qc_reports/`)

To document library yields and technical filtering across all 110 retained libraries, three standard quality control reports were generated:

### 1. DADA2 Read Retention Tracing
Consolidated DADA2 tracking tables (`16S_dada2_filtering_stats.qzv` and `ITS_dada2_filtering_stats.qzv`) quantify read progression through quality filtering, error-model dereplication, paired-end merging, and bimera elimination.

### 2. Decontaminated Sequencing Depth Distributions
Summary artifacts (`16S_sequencing_depth_summary.qzv` and `ITS_sequencing_depth_summary.qzv`) verify the post-decontamination count distribution across all 110 libraries, ensuring that removing off-target host sequences did not cause library dropouts.

### 3. Alpha Rarefaction Saturation Curves
To evaluate whether the chosen sampling depths capture true biological diversity:
* **16S rRNA (`max-depth 27,000`):** Evaluates Observed Features, Shannon entropy, and Faith’s PD across incremental subsampling intervals up to 27,000 reads/sample. Curves confirm that richness estimators reach an asymptotic plateau across all five batches, validating $27,000$ reads as a robust depth that retains 100% of biological replicates ($N=110$).
* **ITS2 (`max-depth 20,000`):** Generates saturation curves up to 20,000 reads. Richness stabilizes by $\sim 6,000\text{--}8,000$ reads, confirming that the operational rarefaction threshold of $8,600$ reads captures community richness without unnecessary sample loss.

---

## 4. Core Diversity Recalculation

Phylogenetic (16S) and non-phylogenetic (ITS2) core diversity pipelines were executed on the decontaminated five-batch tables:

| Parameter / Metric | Bacterial 16S Pipeline | Fungal ITS2 Pipeline |
| :--- | :--- | :--- |
| **Input Feature Table** | `05_taxonomy_16S/table_16S_clean.qza` | `05_taxonomy_ITS/table_ITS_clean.qza` |
| **Phylogeny Artifact** | `06_phylo_16S/rooted_tree_16S.qza` | *None (Bypassed)* |
| **Rarefaction Depth** | **27,000 reads/sample** | **8,600 reads/sample** |
| **Alpha Diversity Indices** | Observed Features, Shannon, Faith's PD, Pielou | Observed Features, Shannon, Pielou |
| **Beta Diversity Metrics** | Bray-Curtis, Jaccard, Weighted & Unweighted UniFrac | Bray-Curtis, Jaccard |
| **Retained Libraries** | $N = 110$ (100% sample retention) | $N = 110$ (100% sample retention) |

---

## 5. Functional Metagenome Inference: PICRUSt2

To bridge taxonomic shifts with functional metabolic potential, de-novo Amplicon Sequence Variants were processed through the **PICRUSt2** (Phylogenetic Investigation of Communities by Reconstruction of Unobserved States) pipeline.

### Algorithmic Architecture

```
04_merged_16S/repseq_16S_merged.qza        05_taxonomy_16S/table_16S_clean.qza
                 │                                          │
       qiime tools export                         qiime tools export
                 │                                          │
       dna-sequences.fasta                        feature-table.biom
                 │                                          │
                 └────────────────────┬─────────────────────┘
                                      ▼
                           picrust2_pipeline.py
                                      │
     ┌────────────────────────────────┼────────────────────────────────┐
     ▼                                ▼                                ▼
1. Phylogenetic Placement      2. Hidden State Prediction      3. Pathway Inference
   (EPA-ng & GAPPA)               (castor R Package)              (MinPath Engine)
   - Place ASVs into reference    - Predict gene family copy      - Infer MetaCyc pathway
     tree without de novo           numbers (EC, KO) per ASV        abundances based on
     MSA errors                   - 16S copy-number normalization   functional gene presence
```

### Theoretical Mechanisms
1. **Phylogenetic Sequence Placement (`EPA-ng` & `GAPPA`):**
   Instead of constructing a global alignment from short amplicons, `EPA-ng` places query ASVs directly into an optimized reference tree derived from thousands of full-length prokaryotic genomes. `GAPPA` computes the most probable placement edge, outputting an annotated JPlace artifact.
2. **Hidden State Prediction (`castor`):**
   Using continuous-time Markov models of evolution, `castor` reconstructs ancestral character states to predict:
   * **16S rRNA Gene Copy Numbers:** Used to normalize raw feature counts, correcting for taxa with high ribosomal operon copy numbers (e.g., *Bacillus*) versus oligotrophic taxa with low copy numbers.
   * **Enzyme Commission (EC) Numbers & KEGG Orthologs (KOs):** Predicts gene family abundance vectors for each ASV.
3. **Nearest Sequenced Taxon Index (NSTI) Gating:**
   PICRUSt2 calculates the NSTI score for each ASV (branch-length distance to the nearest sequenced reference genome). ASVs with high NSTI are flagged, preventing speculative functional inference for novel or unrepresented soil lineages.
4. **Structured Pathway Reconstruction (`MinPath`):**
   Minimal Pathway Inference (`MinPath`) applies parsimony modeling to the predicted EC count matrix, evaluating whether complete biological pathways are present in the community rather than assuming pathway presence based on isolated enzymatic steps. Metabolic output is structured into **MetaCyc metabolic pathways**, providing functional profiles ready for downstream pathway-level differential testing in R via `ggpicrust2`.