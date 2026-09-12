# Micro-Tom Rhizosphere Amplicon Pipeline (16S & ITS2)
Dry-lab computational log and pipeline for profiling bacterial (16S rRNA) and fungal (ITS2) rhizosphere communities of Parent (P) generation *Solanum lycopersicum* cv. 'Micro-Tom' exposed to distinct regional soil inocula.

---

## 1. Project Background & Biological Context

This project evaluates the influence of distinct regional soil environments on plant phenotypes across generations, tracking potential microbiome-mediated and transgenerational memory effects. 

* **Gyeongju (GJ) Soil:** Associated with increased vegetative biomass (elevated Shoot Length [SL], Root Length [RL], Fresh Weight [FW], and Dry Weight [DW]).
* **Gijang B (GB) Soil:** Associated with accelerated flowering and fruit ripening kinetics.

```text
Multi-Generation Experimental Schema:
P Generation  : Soil Inoculation (GJ, GB, or Control) ──▶ [ CURRENT DATASET: P Rhizosphere (N=110) ]
      │
F1 Generation : Continuous Treatment (GJ, GB, or Control)
      │
F2 Generation : Continuous Treatment (GJ, GB, or Control)
      │
F3 Generation : Washout / No Treatment (Evaluation of inherited phenotypic memory)
```

---

## 2. Experimental Setup & Rhizosphere Sampling

* **Host Plant:** *Solanum lycopersicum* cv. 'Micro-Tom'
* **Base Substrate:** Peat plug
* **Treatment Regimen:**
  * **Control (CTRL):** Peat plug + MES buffer
  * **Gyeongju (GJ):** Peat plug + MES buffer + GJ soil suspension
  * **Gijang B (GB):** Peat plug + MES buffer + GB soil suspension
* **Sampled Compartment:** Rhizosphere soil (tightly root-adherent soil fraction) after ~96 days of cultivation.
* **Curated Sample Size:** $N = 110$ biological replicates across 5 validated sequencing batches (GB: $33\text{ Ctrl} + 33\text{ Treat}$; GJ: $22\text{ Ctrl} + 22\text{ Treat}$).

---

## 3. Sequencing Metadata & Sample Architecture

* **Sequencing Platform:** Illumina MiSeq (National Instrumentation Center for Environmental Management, NICEM)
* **Run Configuration:** Paired-End (PE), $2 \times 300$ bp
* **Sample Source:** Rhizosphere soil collected from **Parent (P) generation** *Solanum lycopersicum* cv. 'Micro-Tom'.

### Amplicon Targets & Primer Constructs

The library uses staggered heterogeneity spacers (0–3 bp 'N') inserted between the Illumina adapter overhang and the locus-specific primer to preserve base-calling diversity on the MiSeq flow cell.

#### 1. Bacterial 16S rRNA (V3–V4 Region)
* **Forward Primer Construct (341F):**  
  `5'- TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG [0-3N] CCTACGGGNGGCWGCAG -3'`
  * Overhang: `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`
  * Spacers: `""`, `N`, `NN`, `NNN`
  * Specific Primer: `CCTACGGGNGGCWGCAG`
* **Reverse Primer Construct (805R):**  
  `5'- GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG [0-3N] GACTACHVGGGTATCTAATCC -3'`
  * Overhang: `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`
  * Spacers: `""`, `N`, `NN`, `NNN`
  * Specific Primer: `GACTACHVGGGTATCTAATCC`

#### 2. Fungal ITS (ITS2 Region)
* **Forward Primer (ITS3):** `5'- GCATCGATGAAGAACGCAGC -3'` (+ N-spacers)
* **Reverse Primer (ITS4):** `5'- TCCTCCGCTTATTGATATGC -3'` (+ N-spacers)

---

### Filename Schema & Sample Demultiplexing

Raw sequencing files follow standard Illumina Casava naming conventions:

```text
[Domain]_[Treatment][Batch]-[Replicate]_S[SampleNo]_L001_R[1|2]_001.fastq.gz
```

* **`Domain` (Target Assay):**
  * `B_`: Bacterial 16S rRNA amplicon
  * `F_`: Fungal ITS2 amplicon
* **`Treatment` (Inoculation Condition):**
  * `B`: Gijang B (GB) soil treatment
  * `M`: GB Control (Peat plug + MES buffer)
  * `J`: Gyeongju (GJ) soil treatment
  * `jM`: GJ Control (Peat plug + MES buffer)
* **`Batch`:** Experimental cohort within the P generation (`B1`, `B2`, `B3`, `J1`, `J2`).
* **`Replicate`:** Plant individual replicate within a set (e.g., `-1` through `-11`).
* **`R[1|2]`:** Directional read index (`R1` = Forward, `R2` = Reverse).

---

## 4. Repository Structure

This repository follows the standardized dry-lab log architecture:

```text
.
├── assets/         # High-resolution output figures (Alpha/Beta diversity, Heatmaps, Volcano plots, Forest plots)
├── logs/           # Dated computational logs documenting parameter decisions and run outputs
├── scripts/        # Sequential, executable analysis scripts (Shell, R, QIIME 2, PICRUSt2)
└── README.md       # Project documentation and reproduction instructions
```

---

## 5. Computational Workflow

```text
                  miseq_data/ (Raw FastQ Batches)
                             │
     ┌───────────────────────┴───────────────────────┐
     ▼                                               ▼
[16S Bacteria (5 Batches)]               [ITS2 Fungi (5 Batches)]
     │                                               │
1. Import (Casava Format)                        1. Import (Casava Format)
     │                                               │
2. Cutadapt Primer Trimming                      2. Cutadapt Primer Trimming
   (341F / 805R + N-spacers)                        (ITS3 / ITS4 + N-spacers)
     │                                               │
3. Batch-Specific DADA2                          3. Batch-Specific DADA2
   (Trunc: F=260, R=200)                            (Trunc: F=0, R=0 for length variance)
     │                                               │
4. Merge Tables & Rep-Seqs                       4. Merge Tables & Rep-Seqs
     │                                               │
5. Taxonomy (SILVA 138)                          5. Taxonomy (UNITE v10)
   + Host Filtering (Mito/Chloro/Cyano)             + Non-Fungal Filtering (k__Fungi)
     │                                               │
6. MAFFT + FastTree                             6. [SKIP TREE BUILDING]
     │                                               │
7. Core Diversity (Rarefied: 27,000)             7. Core Diversity (Rarefied: 8,600)
     │                                               │
     ├───────────────────────┬───────────────────────┘
     │                       ▼
     │        8. Export to R (via qiime2R)
     │        9. Differential Abundance (LinDA)
     │        10. Taxonomic Agglomeration (Phylum to Genus)
     ▼
11. Functional Prediction (PICRUSt2)
     │
12. Pathway Differential Abundance (ggpicrust2)
```

---

### Step-by-Step Processing Strategy

#### 1. Ingestion & Demultiplexing
* Raw FastQ directories for the 5 verified experimental batches (`GB1`, `GB2`, `GB3`, `GJ1`, `GJ2`) are imported into QIIME 2 using the `CasavaOneEightSingleLanePerSampleDirFmt` protocol as `SampleData[PairedEndSequencesWithQuality]`.
* Ambiguous libraries from batch `GJ3` (Set 4 and mislabeled $x/y$ extras) were excluded to maintain analytical reproducibility and eliminate circular classification bias.

#### 2. Primer & Heterogeneity Spacer Cleavage
* `qiime cutadapt trim-paired` strips locus-specific primers and variable-length N-spacers (`341F/805R` for 16S, `ITS3/ITS4` for ITS2) under an error tolerance of $0.15$ with `--p-discard-untrimmed`.

#### 3. Batch-Specific DADA2 Denoising
Parametric error models are trained independently within each sequencing run prior to table consolidation:
* **16S V3–V4 Strategy:** Truncation set to `--p-trunc-len-f 260` and `--p-trunc-len-r 200` with $\text{maxEE} = (2.0, 3.0)$. This removes terminal Phred decay while ensuring $> 20$ bp overlap for paired contig assembly across the ~460 bp amplicon.
* **ITS2 Strategy:** Positional truncation is disabled (`--p-trunc-len-f 0 --p-trunc-len-r 0`) to prevent systemic exclusion of naturally short or long fungal ITS amplicons. Quality control is governed strictly by quality score and expected error filtering (`maxEE` $2.0 / 3.0$).

#### 4. Feature Consolidation & Decontamination
* Exact Amplicon Sequence Variants (ASVs) and representative sequences are merged across batches via `qiime feature-table merge` and `qiime feature-table merge-seqs`.
* **16S Filtering:** Classified against **SILVA 138** (99% OTUs). Features annotated as mitochondria, chloroplast, Cyanobacteria, Archaea, Eukaryota, Rickettsiales, Rickettsiellales, or Unassigned are purged.
* **ITS2 Filtering:** Classified against **UNITE v10** (dynamic 99% OTUs). Only features retained within `k__Fungi` are preserved.

#### 5. Phylogenetic Inference vs. Non-Phylogenetic Analysis
* **16S Pipeline:** Alignment via `MAFFT`, hypervariable positional masking, and midpoint-rooted tree construction via `FastTree` to support phylogenetic distance metrics (Faith's PD, Weighted/Unweighted UniFrac).
* **ITS2 Pipeline:** Global alignment and tree construction are bypassed due to non-homologous length variation across fungal phyla. Diversity is calculated strictly using non-phylogenetic distances (Bray-Curtis, Jaccard).

#### 6. Normalization & Ecological Modeling in R
* Feature tables are normalized via rarefaction thresholds established from saturation curves: **27,000 reads/sample** for 16S (100% sample retention) and **8,600 reads/sample** for ITS2.
* Downstream analysis in R (`aarch64`):
  * **Community Structure:** Multi-metric ordination (PCoA and NMDS) and non-parametric multivariate variance partitioning (PERMANOVA / `adonis2`).
  * **Differential Abundance:** Linear Models for Differential Abundance (**LinDA**) applied to synchronous pairwise contrasts (`GB_Treatment` vs. `GB_Control`; `GJ_Treatment` vs. `GJ_Control`), accounting for compositional effects without rank-tied $p$-value artifacts.
  * **Private Biomarker Mining:** Identification of group-exclusive rare taxa establishing uniquely under specific soil inocula.

#### 7. Functional Metagenome Inference (PICRUSt2)
* Decontaminated 16S ASVs and clean BIOM tables are processed through the `picrust2_pipeline.py` architecture (`EPA-ng`, `GAPPA`, `castor`, and `MinPath`).
* Predicted MetaCyc metabolic pathway abundances are analyzed in R using `ggpicrust2` to identify functional enrichments (carbon fixation, amino acid catabolism, secondary metabolism) differentiating regional rhizosphere states.

---

## 6. Software Dependencies & Computing Environments

Computational tasks are partitioned between an emulated x86_64 Conda architecture (for legacy binaries and QIIME 2 compatibility) and a native ARM64 Apple Silicon environment for high-performance statistical modeling in R.

---

### Architecture 1: Emulated x86_64 Environment (`osx-64` via Rosetta 2)

#### 1. Primary Amplicon Pipeline (`qiime2-amp` conda environment)
* **Core Framework:** QIIME 2 Amplicon Distribution (v2024.10.1)
* **Demultiplexing & Primer Trimming:** Cutadapt (v4.9 via `q2-cutadapt`)
* **Denoising & ASV Resolution:** DADA2 (v1.30.0 via `q2-dada2`)
* **Multiple Sequence Alignment:** MAFFT (v7.526)
* **Phylogenetic Reconstruction:** FastTree (v2.1.11)
* **Phylogenetic Distance Metrics:** UniFrac binaries (v1.4)

#### 2. Functional Metagenome Inference (`picrust2` conda environment)
* **Metagenomic Prediction:** PICRUSt2 (v2.6.3)
* **Phylogenetic Placement:** EPA-ng (v0.3.8) and GAPPA (v0.9.0)
* **Hidden State Prediction:** `castor` (v1.8.6)
* **Pathway Parsimony:** MinPath (v1.6)

---

### Architecture 2: Native ARM64 Environment (`aarch64-apple-darwin23`)

* **Base Platform:** R (v4.6.1)

#### Specialized R Packages
* **Artifact Ingestion & Containerization:**
  * `qiime2R` (v0.99.6) — Direct import of `.qza`/`.qzv` artifacts
  * `phyloseq` (v1.56.0) — Central data object orchestration
* **Community Ecology & Ordination:**
  * `vegan` (v2.7-6) — Alpha/Beta diversity, NMDS ordination, and PERMANOVA (`adonis2`)
  * `MicrobiomeStat` (v1.4) / `LinDA` — Linear regression modeling for differential abundance and compositional bias correction
  * `rstatix` (v1.0.0) — Non-parametric post-hoc testing
* **Functional Pathway Evaluation:**
  * `ggpicrust2` (v2.5.17) — MetaCyc pathway differential abundance, forest plotting, and Z-score heatmaps
* **Mixed-Effects & Compositional Modeling:**
  * `lme4` (v2.0-6) & `lmerTest` (v3.2-1) — Mixed-effects modeling accounting for batch structures
  * `zCompositions` (v1.6.2) — Zero-replacement strategies for compositional data
* **Data Visualization & Graphics:**
  * `ggplot2` (v4.0.3), `ggpubr` (v1.0.0), `patchwork` (v1.3.2)
  * `ComplexHeatmap` / `pheatmap` (v1.0.13) — Hierarchical clustering and effect-size heatmaps

---

### Reference Databases & Pre-Trained Classifiers

* **16S rRNA (Bacterial Community):**
  * Database: **SILVA 138** (99% OTUs, full-length 16S)
  * Classifier: Naive Bayes classifier trained via `scikit-learn` (v1.4.2)
  * Target artifact: `silva-138-99-nb-classifier.qza`
* **ITS2 (Fungal Community):**
  * Database: **UNITE v10** (ver10 99% threshold, release 04.04.2024)
  * Classifier: Naive Bayes classifier formatted for QIIME 2
  * Target artifact: `unite_ver10_99_04.04.2024-Q2-2024.5.qza`