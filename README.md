# Micro-Tom Rhizosphere Amplicon Sequencing (16S and ITS2)

Bacterial (16S rRNA, V3–V4) and fungal (ITS2) amplicon sequencing of rhizosphere soil from the parent (P) generation of tomato (*Solanum lycopersicum* cv. Micro-Tom). Plants were grown in peat plugs treated with a soil suspension from one of two sites, Gijang B (GB) or Gyeongju (GJ), or with buffer only. The repository contains the QIIME 2 and R code, dated analysis logs and figures. The final analysis uses 110 samples from five sequencing batches.

---

## 1. Overview

**Question.** Does treatment with soil suspensions from two regions change the bacterial and fungal communities of the Micro-Tom rhizosphere, and how do the communities from the two soils differ?

**Context.** This dataset is the P generation of a multi-generation experiment:

* P: soil treatment (GB, GJ or control) — this dataset;
* F1 and F2: treatment continued;
* F3: no treatment, to test for effects carried over from earlier generations.

In earlier phenotyping (not part of this repository), GJ soil was associated with larger plants (shoot length, root length, fresh and dry weight) and GB soil with earlier flowering and fruit ripening.

---

## 2. Experimental design

* **Plant:** *Solanum lycopersicum* cv. Micro-Tom, grown in peat plugs.
* **Treatments:** control, peat plug + MES buffer; GB, peat plug + MES buffer + GB soil suspension; GJ, peat plug + MES buffer + GJ soil suspension.
* **Two experiments:** the GB and GJ experiments were grown at different times, each with its own control group.
* **Sample:** rhizosphere soil (soil tightly attached to the roots), collected after about 96 days.

Each sequencing batch contains one treated set and one control set of 11 plants each. Sample codes: `B` = GB-treated, `M` = GB control, `J` = GJ-treated, `jM` = GJ control, followed by the set number.

| Batch | Treated set | Control set | Libraries per assay | Used in final analysis |
| :--- | :--- | :--- | ---: | :--- |
| GB1 | `B1` (11) | `M1` (11) | 22 | yes |
| GB2 | `B3` (11) | `M3` (11) | 22 | yes |
| GB3 | `B4` (11) | `M4` (11) | 22 | yes |
| GJ1 | `J2` (11) | `jM2` (11) | 22 | yes |
| GJ2 | `J3` (11) | `jM3` (11) | 22 | yes |
| GJ3 | `J4` (6) | `jM4` (11) | 24, including 7 samples with `x`/`y` codes | no (removed; see below) |

Labels in batch GJ3 were uncertain: seven tubes carried `x`/`y` codes that could not be assigned to a group, and some `J4`/`jM4` samples clustered with the opposite group. An exploratory reassignment of these labels from PCoA positions was circular, so the whole batch was removed (logs 2026-07-24_02 and 2026-07-29_01).

**Sample sets used in the logs**

| Set | Batches | N | `GB_Control` / `GB_Treatment` / `GJ_Control` / `GJ_Treatment` | Labels | Logs |
| :--- | :--- | ---: | :--- | :--- | :--- |
| Full | six | 134 | 33 / 33 / 33 / 28, plus 7 x/y samples | original | 2026-07-23_02, 2026-07-24_01 |
| Standard | six | 127 | 33 / 33 / 33 / 28 | original | 2026-07-23_02, 2026-07-24_01 |
| Reassigned (exploratory) | six | 134 | 33 / 33 / 34 / 34 | 17 GJ3 labels changed from PCoA | 2026-07-24_02, 2026-07-25_01 |
| **Final** | five (no GJ3) | **110** | 33 / 33 / 22 / 22 | original | 2026-07-29_01 to 2026-08-24_01 |

The final set uses the original labels in `metadata_P_generation.tsv`. None of its samples had their labels changed, because all reassigned samples were in GJ3.

---

## 3. Data

* **Sequencing:** Illumina MiSeq, 2 × 300 bp paired-end, National Instrumentation Center for Environmental Management (NICEM).
* **Primers:** each construct has the Illumina adapter overhang, a 0–3 nt N spacer and the locus-specific primer.

```text
16S V3–V4, 341F:  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG [0-3 N] CCTACGGGNGGCWGCAG
16S V3–V4, 805R:  GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG [0-3 N] GACTACHVGGGTATCTAATCC
ITS2, ITS3:       GCATCGATGAAGAACGCAGC  (+ N spacer)
ITS2, ITS4:       TCCTCCGCTTATTGATATGC  (+ N spacer)
```

* **File names:** `[B|F]_[group code][set]-[replicate]_S[sample number]_L001_R[1|2]_001.fastq.gz`, for example `B_jM3-8_S96_L001_R1_001.fastq.gz`. `B_` = 16S, `F_` = ITS2. The seven uncertain GJ3 samples have codes ending in `x` or `y` (for example `B_J2x`).
* **Re-sequenced libraries:** two 16S libraries (`B_J3-4`, `B_M4-1`) had more than one FASTQ pair; the runs were concatenated before import (2026-07-22_01).
* Raw reads, QIIME 2 artifacts and metadata files are not included in the repository.

---

## 4. Methods

Steps 1–6 were run in QIIME 2 2024.10; steps 7–10 in R 4.6.1 unless stated. Details and commands are in the logs.

1. **Import** of each batch separately (Casava 1.8 format).
2. **Primer removal** with cutadapt (`qiime cutadapt trim-paired`): primers given as 5′ adapters (this also removes the spacers), error rate 0.15, IUPAC wildcards allowed, read pairs without primers discarded.
3. **Denoising** with DADA2 for each batch separately. 16S: truncation at 260 (forward) and 200 (reverse), maximum expected errors 2.0 / 3.0. ITS2: no truncation, maximum expected errors 2.0 / 3.0.
4. **Merging** of the batch tables, **taxonomy** with naive Bayes classifiers (16S: SILVA 138, 99%; ITS2: UNITE v10, 99%) and **filtering**. 16S: features assigned to mitochondria, chloroplast, Archaea, Eukaryota, Cyanobacteria, Rickettsiales, Rickettsiellales or Unassigned were removed. ITS2: only `k__Fungi` features were kept. For the final five-batch tables, the taxonomy from the six-batch classification was reused (`pipeline/run_pipeline.sh` instead re-runs the classifiers on the five-batch sequences).
5. **16S tree** with MAFFT and FastTree, midpoint-rooted (`align-to-tree-mafft-fasttree`). No ITS2 tree.
6. **Rarefaction and core diversity metrics:** 27,000 reads per sample (16S) and 8,600 (ITS2).
7. **Alpha and beta diversity:** observed features, Shannon, Gini–Simpson and Faith's PD; Kruskal–Wallis tests with FDR-adjusted pairwise Wilcoxon tests. PCoA (`cmdscale`) and NMDS (`vegan::metaMDS`, k = 2). PERMANOVA (`vegan::adonis2`, `~ Treatment_Group`, 999 permutations, `set.seed(42)`) and dispersion tests (`betadisper`, `permutest`).
8. **Differential abundance of ASVs** with LinDA (`MicrobiomeStat::linda`) on the rarefied table, `~ Treatment_Group`, separately for `GB_Treatment` vs `GB_Control` and `GJ_Treatment` vs `GJ_Control`. ASVs in fewer than 10% of samples were removed (`prev.filter = 0.10`, `mean.abund.filter = 0.0001`). Significance: FDR < 0.05 and |log2FC| ≥ 1. The model has no batch or set term.
9. **Functional prediction** with PICRUSt2 from the unrarefied, taxonomy-filtered 16S table (27,767 ASVs): EPA-ng placement, castor maximum-parsimony prediction, NSTI cut-off 2.0 (171 ASVs excluded), MetaCyc pathways. Pathways were compared with LinDA (`~ Treatment_Group`, `prev.filter = 0.1`, `zero.handling = "pseudo-count"`). `ggpicrust2` was used only to add pathway names.
10. **Taxonomic composition:** mean relative abundance per group from phylum to genus (rarefied table).

---

## 5. Key results (final set, N = 110, 16S)

**Alpha diversity.** Treated rhizospheres had higher alpha diversity than controls (Kruskal–Wallis p < 10⁻¹⁵ for all four metrics). Median observed features were about 985 (`GB_Control`) and 900 (`GJ_Control`), compared with about 1,690 (`GB_Treatment`) and 1,340 (`GJ_Treatment`).

**Beta diversity (PERMANOVA, `~ Treatment_Group`, all p = 0.001).**

| Distance | Four groups, R² | Treated groups only, R² |
| :--- | ---: | ---: |
| Bray–Curtis | 0.342 | 0.290 |
| Weighted UniFrac | 0.303 | 0.196 |
| Unweighted UniFrac | 0.314 | 0.165 |
| Jaccard | 0.263 | 0.195 |

* Dispersion also differed between groups (betadisper p ≤ 0.002 in all eight tests), so these R² values include differences in spread as well as in group centres.
* The two control groups differed in composition (six-batch standard set: R² 0.08–0.11, p = 0.001). Each treatment was therefore compared with its own control.

**Differential abundance (LinDA).** In both contrasts, many more ASVs were enriched in the treated group than in the control. Most of the strongest enriched ASVs were absent or rare in the controls (for example *Acidibacter*, *Chryseolinea* and *Candidatus* Adlerbacteria in GB; WD2101_soil_group, Blastocatellaceae and *Porphyrobacter* in GJ).

**Taxonomic composition.** Both treated groups had less Firmicutes and more Acidobacteriota and Chloroflexi than their controls. In most plots from class to genus level, the most abundant taxa made up a smaller share of reads in treated samples.

**PICRUSt2.** Seven predicted MetaCyc pathways were among the 20 most significant in both contrasts, all higher in the treated groups (for example the 3-hydroxypropanoate cycle and glyoxylate assimilation). These are predictions from 16S data, not measured functions.

**ITS2** (six-batch full set only). Treated groups had more fungal ASVs than controls (Kruskal–Wallis p < 2 × 10⁻¹⁶). PERMANOVA R² was 0.258 (Bray–Curtis) and 0.187 (Jaccard). GB and GJ samples, including the two control groups, separated along PCoA 1, so this axis reflects the experiment as well as the soil.

---

## 6. Limitations and open points

* **Unbalanced and partly confounded design.** The GJ groups have 22 samples from two batches; the GB groups have 33 from three. Each batch contains one treated and one control set, and the GB and GJ experiments were grown at different times. GB vs GJ comparisons therefore mix soil origin with experiment (time), and set-to-set differences are visible in several analyses (2026-07-23_02, 2026-08-10_01, 2026-08-24_01). No model includes a set or batch term.
* **Metadata column name.** The column `Sequencing_Batch` holds the set code (`B1`, `M1`, …), not the batch name (`GB1`, …).
* **LinDA input.** LinDA was run on the rarefied table (Section 4, step 8).
* **Dispersion.** Group dispersions differ, which PERMANOVA does not separate from location effects.
* **PICRUSt2.** Several top pathways are named after archaeal, yeast or animal pathways, so the pathway results need checking before biological interpretation.
* **Exploratory logs.** Logs 2026-07-24_02 and 2026-07-25_01 use PCoA-derived labels. They are kept as a record and are not used for the final results.
* **Figure 088** is identical to Figure 045, the six-batch (N = 134) exclusive-ASV table.
* **Figure values.** Taxonomic percentages in logs 2026-07-25_01 and 2026-08-18_01 were measured from the bar heights in the figures (about ±0.2 percentage points), because the scripts do not export the group means.
* **ITS2.** The ITS2 results come from the six-batch set only (2026-07-24_01). The five-batch pipeline computes ITS2 diversity metrics (2026-07-29_01), but they were not analysed further. The number of ITS2 samples retained at 8,600 reads is not recorded.

---

## 7. How to run

```bash
# QIIME 2 and PICRUSt2, five-batch design (macOS; conda environments qiime2-amp and picrust2, osx-64)
bash pipeline/run_pipeline.sh

# R analyses of the final set (N = 110)
Rscript pipeline/run_analysis.R
```

* `pipeline/run_pipeline.sh` expects the raw folders `Bacteria_<batch>/` and `Fungi_<batch>/` (GB1–GB3, GJ1–GJ2) in the working directory.
* The DADA2 read-tracking step uses `merge_16S_dada2_stats.py` and `merge_ITS_dada2_stats.py`, which are not in the repository. The step is skipped if they are missing.
* `pipeline/run_analysis.R` sets the working directory with `setwd("/Users/heechan/Desktop/miseq_data_without_GJSet4")`; change this path before running. The script prints figures to the graphics device instead of saving them. The PICRUSt2 sections write CSV files to `figure/`, which must exist.
* `modules/` contains the same steps as numbered scripts (module-to-log mapping in Section 8). Modules 01–05 import and process only the five final batches. Modules 06–09 reproduce the six-batch analyses (N = 127 and 134) only if GJ3 is added back to modules 01–04.
* **Module 10 deletes the raw GJ3 folders** (`rm -rf Bacteria_GJ3 Fungi_GJ3`) and earlier outputs. Keep a backup before running it.
* Most figures in `assets/` were saved from interactive R sessions. The `ggsave()` calls in `modules/` were added later and use other image sizes; rerunning a module overwrites the corresponding figure. Modules 06–08 do not save every figure (for example 019, 020 and 028–040) and partly use other file names (module 08: `Bacteria_16S_Final_*.png`). Figures 001–014 are screenshots from QIIME 2 View; Figures 045 and 088 are screenshots of the exclusive-ASV table.

---

## 8. Repository layout

```text
.
├── README.md
├── assets/     # figures 001–088, numbered in log order
├── logs/       # analysis logs, 2026-07-22 to 2026-08-24
├── modules/    # step-by-step scripts (see table) and make_metadata.py
└── pipeline/   # run_pipeline.sh (QIIME 2, PICRUSt2), run_analysis.R (R, N = 110)
```

| Log | Content | Sample set | Modules | Figures |
| :--- | :--- | :--- | :--- | :--- |
| 2026-07-22_01 | environment, FASTQ merging, import | six batches | 00, 01 | — |
| 2026-07-22_02 | primer removal, DADA2 | six batches | 02, 03 | 001–005 |
| 2026-07-23_01 | merging, tree, taxonomy, metadata | six batches | 04, 05 | 006–008 |
| 2026-07-23_02 | core metrics, QIIME 2 PERMANOVA, x/y samples | full, standard | 05 | 009–014 |
| 2026-07-24_01 | alpha and beta diversity in R; ITS2 | full, standard | 06 | 015–022 |
| 2026-07-24_02 | label reassignment (exploratory) | reassigned | 07 | 023–030 |
| 2026-07-25_01 | diversity, taxonomy, LinDA (exploratory) | reassigned | 08, 09 | 031–045 |
| 2026-07-29_01 | removal of GJ3, five-batch rebuild, PICRUSt2 | final | 10 | — |
| 2026-08-10_01 | alpha and beta diversity | final | 11 | 046–066 |
| 2026-08-18_01 | taxonomic composition | final | 12 | 067–077 |
| 2026-08-24_01 | LinDA, PICRUSt2 pathways | final | 13 | 078–088 |

---

## 9. Software

* **QIIME 2** amplicon distribution 2024.10 (conda environment `qiime2-amp`, `osx-64` under Rosetta 2), with cutadapt 4.9, DADA2 1.30.0, MAFFT 7.526, FastTree 2.1.11 and UniFrac 1.4.
* **PICRUSt2** 2.6.3 (conda environment `picrust2`), with EPA-ng 0.3.8, GAPPA 0.9.0, castor 1.8.6 and MinPath 1.6.
* **R** 4.6.1 (`aarch64-apple-darwin23`): qiime2R 0.99.6, vegan 2.7-6, MicrobiomeStat 1.4 (LinDA), rstatix 1.0.0, ggpicrust2 2.5.17 (pathway names only), ggplot2 4.0.3, ggpubr 1.0.0, patchwork 1.3.2, pheatmap 1.0.13, ggrepel, RColorBrewer, tidyverse. phyloseq 1.56.0 and ComplexHeatmap are installed or loaded by some scripts but are not used for the reported results.
* **Python 3** (`modules/make_metadata.py`): pandas.

---

## 10. Reference databases

* **16S:** SILVA 138, 99% OTUs, full-length 16S; QIIME 2 naive Bayes classifier `silva-138-99-nb-classifier.qza` (scikit-learn 1.4.2), from data.qiime2.org.
* **ITS2:** UNITE v10, 99%, release 04.04.2024; classifier `unite_ver10_99_04.04.2024-Q2-2024.5.qza`, from github.com/colinbrislawn/unite-train.
* **PICRUSt2:** default reference files; MetaCyc pathway map `metacyc_pathways_structured_filtered_v24.txt`.
