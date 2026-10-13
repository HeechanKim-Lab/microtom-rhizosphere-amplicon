# Environment Setup, FASTQ Remediation and Import

**Date:** 2026-07-22
**Stage:** QIIME 2 environment; merging of re-sequenced libraries; import of raw reads
**Cohort:** six batches (GB1–GB3, GJ1–GJ3), 16S and ITS2; 134 libraries per assay
**Script:** `modules/01_remediation_and_import.sh` (as committed, it imports only five batches; see Section 3)
**Input:** `Bacteria_<batch>/`, `Fungi_<batch>/` (paired-end FASTQ, Casava 1.8 file names)
**Output:** `01_imported/demux_16S_<batch>.qza`, `01_imported/demux_ITS_<batch>.qza`

---

## 1. Environment

QIIME 2 was installed as an `osx-64` conda environment on an Apple Silicon Mac (macOS Sonoma, Darwin 23) and runs under Rosetta 2.

```bash
CONDA_NO_PLUGINS=true CONDA_SUBDIR=osx-64 conda env create --solver=classic -n qiime2-amp \
  --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml
conda activate qiime2-amp
conda config --env --set subdir osx-64
```

* QIIME 2 amplicon distribution 2024.10, Python 3.10.
* `CONDA_SUBDIR=osx-64` forces x86_64 packages. QIIME 2 depends on compiled tools (e.g. DADA2, UniFrac, FastTree) that are distributed for x86_64.
* The environment was built with the classic solver and with conda plugins disabled (`CONDA_NO_PLUGINS=true`).

---

## 2. Merging re-sequenced libraries

The Casava import format expects one R1/R2 file pair per sample. Two 16S libraries had more than one pair because they had been re-sequenced. Their runs were combined before import. Gzip files can be concatenated with `cat` and are read as one stream by gzip-compatible tools, so the files were not decompressed.

| Batch | Sample | File pairs | Action |
| :--- | :--- | :--- | :--- |
| `Bacteria_GJ2` | `B_J3-4` | `S7`, `S79` | Concatenated (R1 and R2 separately) and saved under the `S79` name. Originals moved to `Bacteria_GJ2/raw_runs/`. |
| `Bacteria_GB3` | `B_M4-1` | `S45`, `S78`, `S96` | `S45` (marked as a failed run in the script) moved to `raw_runs/`. `S78` and `S96` concatenated and saved under the `S96` name. Originals moved to `raw_runs/`. |

Commands for `B_M4-1` (the `B_J3-4` steps are analogous):

```bash
mkdir -p Bacteria_GB3/raw_runs
mv Bacteria_GB3/B_M4-1_S45_* Bacteria_GB3/raw_runs/
cat Bacteria_GB3/B_M4-1_S78_L001_R1_001.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz > Bacteria_GB3/B_M4-1_merged_R1.fastq.gz
cat Bacteria_GB3/B_M4-1_S78_L001_R2_001.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R2_001.fastq.gz > Bacteria_GB3/B_M4-1_merged_R2.fastq.gz
mv Bacteria_GB3/B_M4-1_S78_* Bacteria_GB3/B_M4-1_S96_* Bacteria_GB3/raw_runs/
mv Bacteria_GB3/B_M4-1_merged_R1.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R1_001.fastq.gz
mv Bacteria_GB3/B_M4-1_merged_R2.fastq.gz Bacteria_GB3/B_M4-1_S96_L001_R2_001.fastq.gz
```

---

## 3. Import

Each batch was imported into its own artifact, so that DADA2 could learn a separate error model for each batch (see 2026-07-22_02).

```bash
mkdir -p 01_imported
for batch in GB1 GB2 GB3 GJ1 GJ2 GJ3; do
  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path "Bacteria_${batch}" \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path "01_imported/demux_16S_${batch}.qza"
done
```

The ITS2 loop is the same, with `Fungi_${batch}` as input and `demux_ITS_${batch}.qza` as output. All 12 imports finished without errors.

Libraries per batch (from the FASTQ file listing in the run log):

| Batch | Sample sets (code prefix) | 16S | ITS2 |
| :--- | :--- | ---: | ---: |
| GB1 | `B1` (GB soil), `M1` (control) | 22 | 22 |
| GB2 | `B3`, `M3` | 22 | 22 |
| GB3 | `B4`, `M4` | 22 (after merging) | 22 |
| GJ1 | `J2` (GJ soil), `jM2` (control) | 22 | 22 |
| GJ2 | `J3`, `jM3` | 22 (after merging) | 22 |
| GJ3 | `J4` (6), `jM4` (11), 7 samples with `x`/`y` codes | 24 | 24 |
| **Total** | | **134** | **134** |

**Note.** `modules/01_remediation_and_import.sh` as committed imports only GB1–GB3 and GJ1–GJ2 (the final five-batch design; see 2026-07-29_01). The six-batch import above is what was run on 22 July and is the input for the logs dated 2026-07-22 to 2026-07-25.

---

## 4. Output folders

The commands write to these folders (not included in the repository):

```text
01_imported/      imported reads (.qza)
02_trimmed_16S/   02_trimmed_ITS/     primer-trimmed reads
03_denoised_16S/  03_denoised_ITS/    DADA2 tables, representative sequences, statistics
04_merged_16S/    04_merged_ITS/      merged tables and sequences
05_taxonomy_16S/  05_taxonomy_ITS/    taxonomy and filtered tables
06_phylo_16S/                         alignment and tree (16S only)
07_diversity_16S/ 07_diversity_ITS/   core diversity outputs
```
