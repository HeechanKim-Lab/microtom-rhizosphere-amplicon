# Taxonomic Composition from Phylum to Genus (N = 110)

**Date:** 2026-08-18
**Stage:** mean relative abundance per group at phylum, class, order, family and genus level
**Cohort:** five batches, original labels; N = 110 (`GB_Control` 33, `GB_Treatment` 33, `GJ_Control` 22, `GJ_Treatment` 22)
**Script:** `modules/12_taxonomic_profiling.R`; the same code is in `pipeline/run_analysis.R` (sections 3-1 to 3-5)
**Platform:** R 4.6.1 (`aarch64-apple-darwin23`)
**Input:** `07_diversity_16S/core_metrics/rarefied_table.qza`, `05_taxonomy_16S/taxonomy_16S.qza`, `metadata_P_generation.tsv`

---

## 1. Methods

* Relative abundances were computed per sample from the rarefied table (27,000 reads per sample). Taxonomy was parsed with `qiime2R::parse_taxonomy()`.
* Reads without a name at the rank being plotted (empty, `NA` or a bare prefix such as `g__`) were labelled "Unassigned".
* For each rank, taxa were ranked by their mean relative abundance over the samples in the plot, with "Unassigned" excluded from the ranking. The four-group plots show the top 10 phyla and the top 12 classes, orders, families and genera. The GB-only and GJ-only plots show the top 10 phyla, top 20 orders and top 30 genera, ranked within that experiment.
* All remaining taxa, including "Unassigned", were summed per sample into "Other" before the group means were taken, so each bar adds up to 100%.
* The GB-only and GJ-only plots use the same colour for the same taxon. The four-group plots use a separate palette.
* "uncultured" is a placeholder name that SILVA uses in many unrelated lineages. The scripts group taxa by name, so the "uncultured" genus bar pools ASVs from different families.
* No statistical test was applied at these ranks. Differential abundance was tested at ASV level (2026-08-24_01).
* The scripts do not export the group means. The percentages below were measured from the bar heights in the figures (calibrated on the axis ticks) and are accurate to about ±0.2 percentage points.

---

## 2. Phylum

![Phylum, four groups](../assets/067_16S_Relative_Abundance(Phylum).png)

![Phylum, GB](../assets/068_16S_GB_Relative_Abundance(Phylum).png)

![Phylum, GJ](../assets/069_16S_GJ_Relative_Abundance(Phylum).png)

| Phylum (%) | GB_Control | GB_Treatment | GJ_Control | GJ_Treatment |
| :--- | ---: | ---: | ---: | ---: |
| Proteobacteria | 26.0 | 29.6 | 29.5 | 28.9 |
| Bacteroidota | 15.5 | 15.3 | 17.9 | 16.0 |
| Actinobacteriota | 10.0 | 9.7 | 12.0 | 11.2 |
| Planctomycetota | 9.8 | 9.0 | 9.5 | 7.7 |
| Verrucomicrobiota | 8.9 | 10.0 | 8.2 | 6.9 |
| Patescibacteria | 7.7 | 6.7 | 6.2 | 8.5 |
| Firmicutes | 8.7 | 4.2 | 6.0 | 4.6 |
| Chloroflexi | 2.9 | 3.7 | 1.8 | 2.9 |
| Myxococcota | 2.9 | 3.0 | 1.8 | 2.8 |
| Acidobacteriota | 1.6 | 3.3 | 1.3 | 4.3 |
| Other | 6.1 | 5.5 | 5.8 | 6.1 |

The GB-only and GJ-only plots contain the same ten phyla; their group means are the same as in the four-group plot.

* Compared with their controls, both treated groups had less Firmicutes (GB 8.7 → 4.2%; GJ 6.0 → 4.6%) and more Acidobacteriota (GB 1.6 → 3.3%; GJ 1.3 → 4.3%) and Chloroflexi (GB 2.9 → 3.7%; GJ 1.8 → 2.9%).
* Proteobacteria were higher in `GB_Treatment` than in `GB_Control` (29.6 vs 26.0%) but similar in the two GJ groups (29.5 vs 28.9%).
* In GJ, Planctomycetota (9.5 → 7.7%) and Verrucomicrobiota (8.2 → 6.9%) were lower and Patescibacteria (6.2 → 8.5%) higher in the treated group. In GB, Patescibacteria were lower (7.7 → 6.7%).
* "Other" was 5.5–6.1% in all groups.

---

## 3. Class

![Class, four groups](../assets/070_16S_Relative_Abundance(Class).png)

| Class (%) | GB_Control | GB_Treatment | GJ_Control | GJ_Treatment |
| :--- | ---: | ---: | ---: | ---: |
| Alphaproteobacteria | 17.8 | 19.6 | 19.9 | 19.1 |
| Bacteroidia | 14.4 | 14.6 | 17.4 | 15.2 |
| Gammaproteobacteria | 8.1 | 9.9 | 9.6 | 9.8 |
| Actinobacteria | 7.6 | 7.3 | 10.1 | 8.7 |
| Verrucomicrobiae | 7.0 | 7.4 | 7.4 | 5.6 |
| Bacilli | 8.7 | 4.2 | 6.0 | 4.6 |
| Planctomycetes | 6.6 | 5.1 | 6.6 | 4.1 |
| Parcubacteria | 3.7 | 3.8 | 2.6 | 4.7 |
| Phycisphaerae | 3.2 | 3.7 | 2.9 | 3.6 |
| Saccharimonadia | 2.7 | 2.5 | 3.5 | 3.0 |
| Polyangia | 2.6 | 2.7 | 1.4 | 2.6 |
| Chlamydiae | 1.8 | 2.7 | 0.8 | 1.3 |
| Other | 15.7 | 16.6 | 11.7 | 17.7 |

* Bacilli were lower in both treated groups, in line with Firmicutes.
* Gammaproteobacteria were higher in `GB_Treatment` than in `GB_Control` (9.9 vs 8.1%) and similar in the two GJ groups.
* "Other" (classes outside the top 12) was larger in the treated groups (GB 15.7 → 16.6%; GJ 11.7 → 17.7%).

---

## 4. Order

![Order, four groups](../assets/071_16S_Relative_Abundance(Order).png)

| Order (%) | GB_Control | GB_Treatment | GJ_Control | GJ_Treatment |
| :--- | ---: | ---: | ---: | ---: |
| Rhizobiales | 6.5 | 8.1 | 7.3 | 7.3 |
| Chitinophagales | 5.8 | 6.4 | 5.6 | 5.9 |
| Burkholderiales | 5.1 | 5.5 | 6.4 | 6.5 |
| Bacillales | 6.6 | 3.3 | 4.6 | 3.6 |
| Sphingobacteriales | 4.0 | 3.0 | 7.0 | 4.7 |
| Sphingomonadales | 2.6 | 3.2 | 3.5 | 4.6 |
| Other | 50.5 | 52.1 | 44.9 | 50.5 |

![Order, GB (top 20)](../assets/072_16S_GB_Relative_Abundance(Order).png)

![Order, GJ (top 20)](../assets/073_16S_GJ_Relative_Abundance(Order).png)

* **GB (top 20 orders).** The top 20 covered 64.5% (control) and 64.6% (treated) of reads. Rhizobiales (6.5 → 8.1%), Cytophagales (1.7 → 2.6%), Chlamydiales (1.8 → 2.7%) and Xanthomonadales (1.4 → 2.1%) were higher in `GB_Treatment`; Bacillales (6.6 → 3.3%), Sphingobacteriales (4.0 → 3.0%) and Propionibacteriales (2.3 → 1.5%) were lower.
* **GJ (top 20 orders).** The top 20 covered 72.1% (control) and 62.1% (treated) of reads. Sphingobacteriales (7.0 → 4.7%), Micropepsales (4.1 → 2.6%), Planctomycetales (3.1 → 1.5%) and Gemmatales (2.0 → 0.8%) were lower in `GJ_Treatment`; Sphingomonadales (3.5 → 4.6%) and Polyangiales (1.0 → 1.9%) were higher.
* Gemmatales, *Candidatus* Kaiserbacteria, Streptosporangiales and Verrucomicrobiales are among the GJ top 20 but not the GB top 20. Chlamydiales, Cytophagales, Pirellulales and Xanthomonadales are among the GB top 20 only.

---

## 5. Family

![Family, four groups](../assets/074_16S_Relative_Abundance(Family).png)

| Family (%) | GB_Control | GB_Treatment | GJ_Control | GJ_Treatment |
| :--- | ---: | ---: | ---: | ---: |
| Chitinophagaceae | 5.6 | 6.1 | 5.2 | 5.6 |
| Sphingobacteriaceae | 3.5 | 2.2 | 6.4 | 3.7 |
| Bacillaceae | 5.3 | 2.2 | 3.3 | 2.6 |
| Opitutaceae | 3.9 | 4.0 | 3.3 | 2.3 |
| Sphingomonadaceae | 2.6 | 3.2 | 3.5 | 4.6 |
| Nocardioidaceae | 2.3 | 1.4 | 3.0 | 1.9 |
| Top 12 families combined | 37.4 | 32.7 | 42.9 | 35.6 |
| Other | 62.6 | 67.3 | 57.1 | 64.4 |

* The 12 most abundant families covered less of the community in the treated groups (GB 37.4 → 32.7%; GJ 42.9 → 35.6%).
* Bacillaceae, Sphingobacteriaceae and Nocardioidaceae were lower in both treated groups; Sphingomonadaceae were higher.

---

## 6. Genus

![Genus, four groups](../assets/075_16S_Relative_Abundance(Genus).png)

| Genus (%) | GB_Control | GB_Treatment | GJ_Control | GJ_Treatment |
| :--- | ---: | ---: | ---: | ---: |
| uncultured | 11.4 | 13.8 | 10.0 | 10.6 |
| *Bacillus* | 4.6 | 2.0 | 1.4 | 1.3 |
| *Mucilaginibacter* | 2.2 | 1.5 | 5.1 | 2.9 |
| *Flavobacterium* | 2.7 | 2.2 | 3.4 | 2.4 |
| *Nocardioides* | 2.3 | 1.4 | 3.0 | 1.8 |
| *Massilia* | 1.7 | 1.0 | 2.5 | 2.2 |
| WD2101_soil_group | 1.8 | 2.1 | 1.0 | 2.1 |
| Top 12 genera combined | 33.0 | 29.8 | 34.3 | 30.4 |
| Other | 67.0 | 70.2 | 65.7 | 69.6 |

* "uncultured" was the largest category in all groups (10.0–13.8%); see Methods for why it is not a single genus. No named genus exceeded 5.1%.
* *Bacillus* was lower in `GB_Treatment` than in `GB_Control` (2.0 vs 4.6%). It was low in both GJ groups (1.4 and 1.3%).
* *Mucilaginibacter*, *Flavobacterium* and *Nocardioides* were lower in both treated groups. WD2101_soil_group was higher in both.

![Genus, GB (top 30)](../assets/076_16S_GB_Relative_Abundance(Genus).png)

![Genus, GJ (top 30)](../assets/077_16S_GJ_Relative_Abundance(Genus).png)

* **GB (top 30 genera).** The top 30 covered 49.2% (control) and 45.3% (treated) of reads. *Bacillus* (4.6 → 2.0%), *Gemmata* (1.5 → 0.5%) and *Paenibacillus* (1.0 → 0.5%) were lower in `GB_Treatment`; *Sphingobium* (0.5 → 1.0%), *Candidatus* Adlerbacteria (0.7 → 1.0%), *Opitutus* (0.9 → 1.4%) and WD2101_soil_group (1.8 → 2.1%) were higher.
* **GJ (top 30 genera).** The top 30 covered 53.0% (control) and 46.1% (treated) of reads. *Mucilaginibacter* (5.1 → 2.9%), *Flavobacterium* (3.4 → 2.4%), *Nocardioides* (3.0 → 1.8%) and *Gemmata* (1.6 → 0.6%) were lower in `GJ_Treatment`; *Sphingobium* (1.1 → 1.9%), WD2101_soil_group (1.0 → 2.1%) and BIrii41 (0.8 → 1.5%) were higher.
* *Abditibacterium*, SH-PL14, *Schlesneria*, *Gemmatimonas*, *Taibaiella*, *Rugosimonospora* and *Micropepsis* are among the GJ top 30 only. *Caulobacter*, *Terrimonas*, *Paenibacillus*, *Edaphobaculum*, *Planctomicrobium*, *Burkholderia-Caballeronia-Paraburkholderia* and Parcubacteria are among the GB top 30 only. Being outside one experiment's top 30 does not mean that a genus was absent there.

---

## 7. Summary

* Firmicutes (Bacilli, Bacillales, Bacillaceae) were lower in both treated groups than in their controls. At genus level, the drop was seen for *Bacillus* in GB only.
* Acidobacteriota and Chloroflexi were higher in both treated groups.
* From class to genus level, the most abundant taxa made up a smaller share of reads in the treated groups, and "Other" was larger. The exceptions are the phylum plots and the GB top-20 order plot, where "Other" did not change. This is consistent with the higher number of observed ASVs in treated samples (2026-08-10_01). However, "Other" also contains reads without a name at that rank, and the scripts do not separate the two.
* These are descriptive group means of an unbalanced design (33 + 33 GB samples from three batches, 22 + 22 GJ samples from two batches).
