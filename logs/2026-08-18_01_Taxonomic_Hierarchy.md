# Multi-Level Taxonomic Profiling (Phylum to Genus) and Lineage-Specific Succession
**Date:** August 18, 2026  
**Project:** 04_MicroTom_P_Rhizosphere_Amplicon_Analysis  
**Target Assays:** Bacterial 16S rRNA (V3–V4)  
**Platform:** R v4.6.1 (`aarch64-apple-darwin23`, Native ARM64)  
**Input Artifacts:** `05_taxonomy_16S/table_16S_clean.qza`, `05_taxonomy_16S/taxonomy_16S.qza`, `metadata_P_generation.tsv`  

---

## 1. Taxonomic Agglomeration Framework & Visual Standardization

Following the pruning of compromised libraries, the verified $N = 110$ cohort (`GB_Control`: $33$, `GB_Treatment`: $33$, `GJ_Control`: $22$, `GJ_Treatment`: $22$) was systematically evaluated across five successive taxonomic ranks: Phylum, Class, Order, Family, and Genus. 

To track biological shifts between conditions without visual distortion, graphics follow standardized aesthetic parameters:
1. **Shared Clade Palette:** Identical taxonomic clades share exact hex-color mappings across global multi-group panels and pairwise lineage contrasts.
2. **Dynamic Top-Taxa Ranking:** Clades falling outside the top-ranking abundance cutoffs within each rank are collapsed into the `Other` category to benchmark rare-biosphere expansion against the dominant core.
3. **Synchronous Baseline Evaluation:** Treatment-induced community shifts are benchmarked directly against their paired temporal controls (`GB_Treatment` vs. `GB_Control`; `GJ_Treatment` vs. `GJ_Control`) to decouple true biological recruitment from timeline-dependent baseline drift.

---

## 2. Phylum-Level Community Architecture

![Global Phylum-Level Relative Abundance](../assets/067_16S_Relative_Abundance(Phylum).png)

### Global Phylum Distribution
The rhizosphere community across all conditions is dominated by 10 core phyla:
* **`Proteobacteria`:** Represents the largest single fraction across all groups ($26.1\%\text{--}29.8\%$), showing a slight expansion in `GB_Treatment` ($29.8\%$) relative to `GB_Control` ($26.1\%$).
* **`Bacteroidota`:** Maintains stable representation across conditions ($13.8\%\text{--}17.3\%$), peaking in `GJ_Control` ($17.3\%$) and stabilizing at $\sim 15.0\%$ in both soil treatments.
* **`Actinobacteriota`:** Represents $11.4\%\text{--}14.8\%$ of total abundance, maintaining consistent colonization across peat baselines and soil inocula.
* **`Planctomycetota` ($7.4\%\text{--}8.8\%$) & `Verrucomicrobiota` ($7.1\%\text{--}8.4\%$):** Form stable secondary rhizosphere populations.
* **`Acidobacteriota`:** Displays the most pronounced treatment-induced expansion, increasing from $1.6\%$ in `GB_Control` to $4.2\%$ in `GB_Treatment`, and from $1.3\%$ in `GJ_Control` to $3.5\%$ in `GJ_Treatment`.

---

### Pairwise Phylum Contrasts (GB vs. GJ Lineages)

![GB Lineage Phylum-Level Relative Abundance](../assets/068_16S_GB_Relative_Abundance(Phylum).png)

![GJ Lineage Phylum-Level Relative Abundance](../assets/069_16S_GJ_Relative_Abundance(Phylum).png)

* **Gijang B Lineage Shifts (Figure 068):**
  * Inoculation with GB soil suspension triggers concurrent increases in `Proteobacteria` ($+3.7\%$), `Acidobacteriota` ($+2.6\%$), and `Chloroflexi` ($+1.1\%$).
  * The `Other` fraction expands from $6.1\%$ in `GB_Control` to $8.8\%$ in `GB_Treatment`, reflecting the establishment of low-abundance exotic phyla introduced by field soil.
  * Baseline peat dominance in `Firmicutes` ($8.2\% \to 5.3\%$) and `Patescibacteria` ($8.5\% \to 6.6\%$) contracts due to competitive displacement and dilution.
* **Gyeongju Lineage Shifts (Figure 069):**
  * GJ soil inoculation drives a similar expansion of `Acidobacteriota` ($1.3\% \to 3.5\%$) and `Chloroflexi` ($1.7\% \to 2.8\%$), while `Proteobacteria` remains steady ($\sim 29.1\%$).
  * `Verrucomicrobiota` contracts from $8.4\%$ in `GJ_Control` to $7.1\%$ in `GJ_Treatment`.
  * The `Other` fraction increases from $6.3\%$ to $8.7\%$, confirming that rare-biosphere recruitment is a shared feature of soil suspension exposure.

---

## 3. Class-Level Taxonomic Succession

![Global Class-Level Relative Abundance](../assets/070_16S_Relative_Abundance(Class).png)

### Proteobacterial & Deep-Lineage Dynamics
Resolving phyla into distinct classes uncovers sub-lineage divergence between conditions:
* **`Alphaproteobacteria`:** Dominates the proteobacterial fraction across all samples ($17.8\%\text{--}20.1\%$), maintaining stable relative abundance regardless of soil inoculation.
* **`Gammaproteobacteria`:** Exhibits noticeable enrichment following soil inoculation, expanding from $9.6\%$ in `GB_Control` to $12.1\%$ in `GB_Treatment`, and from $9.8\%$ in `GJ_Control` to $11.6\%$ in `GJ_Treatment`.
* **`Bacteroidia`:** Mirrors phylum-level Bacteroidota trends, representing $13.8\%\text{--}17.3\%$ of reads across groups.
* **`Actinobacteria`:** Accounts for $7.8\%\text{--}11.9\%$ of community abundance, exhibiting highest baseline dominance in `GJ_Control` ($11.9\%$) before stabilizing at $8.8\%$ in `GJ_Treatment`.
* **`Bacilli` Contraction:** Displays selective suppression under soil treatments, declining from $8.0\%$ in `GB_Control` to $5.2\%$ in `GB_Treatment`, and from $7.9\%$ in `GJ_Control` to $5.4\%$ in `GJ_Treatment`.
* **Unassigned Rare Fraction:** The `Other` class category expands from $\sim 15.5\%$ in uninoculated peat controls to $> 18.0\%$ in both soil treatments.

---

## 4. Order-Level Stratification & Expanded Lineage Profiling

![Global Order-Level Relative Abundance](../assets/071_16S_Relative_Abundance(Order).png)

### Four-Group Order Landscape
At the order level, the rhizosphere partitions into specialized functional groups:
* **`Rhizobiales`:** The dominant bacterial order across all four conditions, maintaining $6.8\%\text{--}7.9\%$ relative abundance in uninoculated controls and expanding slightly to $7.8\%\text{--}8.2\%$ in treated rhizospheres.
* **`Chitinophagales`:** Represents $6.1\%\text{--}6.9\%$ of total reads across conditions, supporting chitinolytic degradation in the root zone.
* **`Burkholderiales`:** Represents the primary betaproteobacterial lineage, stable at $5.8\%\text{--}6.5\%$.
* **`Bacillales`:** Reflects class-level Bacilli dynamics, dropping from $6.2\%$ in controls to $4.1\%$ in treatments.
* **Dominance Dilution:** In the global 4-group panel, the top 12 orders account for $\sim 50\%$ of reads in controls, whereas the `Other` fraction expands to $52.1\%$ in `GB_Treatment` and $50.8\%$ in `GJ_Treatment`.

---

### Lineage-Specific Order Contrasts (Expanded Top Rankings)

To capture fine-scale ecological turnover, order-level agglomeration was expanded to include the top 20 orders within each geographic lineage.

![GB Lineage Order-Level Relative Abundance](../assets/072_16S_GB_Relative_Abundance(Order).png)

![GJ Lineage Order-Level Relative Abundance](../assets/073_16S_GJ_Relative_Abundance(Order).png)

* **Gijang B Specific Succession (Figure 072):**
  * `Rhizobiales` expands from $6.8\%$ to $8.2\%$ in `GB_Treatment`.
  * `Chitinophagales` increases from $5.5\%$ to $6.4\%$.
  * Specialized soil orders emerge in `GB_Treatment`, including `Micropepsales` ($1.8\% \to 2.4\%$), `Xanthomonadales` ($1.6\% \to 2.1\%$), and `Streptomycetales` ($1.2\% \to 1.8\%$).
  * The top 20 orders account for $64.4\%$ in `GB_Control` and $64.6\%$ in `GB_Treatment`, with `Other` remaining steady at $\sim 35.5\%$.
* **Gyeongju Specific Succession (Figure 073):**
  * `GJ_Treatment` exhibits distinct order recruitment: `Sphingobacteriales` ($6.8\% \to 7.2\%$), `Sphingomonadales` ($3.9\% \to 4.6\%$), and `Planctomycetales` ($2.0\% \to 2.6\%$).
  * Several orders appear exclusively among top rankings in the Gyeongju lineage: `Gemmatales`, `Candidatus Kaiserbacteria` (Patescibacteria clade), and `Streptosporangiales`.
  * The expanded top 20 orders encompass $72.1\%$ of reads in `GJ_Control` but contract to $61.8\%$ in `GJ_Treatment`, with the `Other` pool expanding from $27.9\%$ to **$38.2\%$**. This indicates that GJ soil introduces an exceptionally diverse, low-abundance order-level tail.

---

## 5. Family-Level Composition & Rare Biosphere Dilution

![Global Family-Level Relative Abundance](../assets/074_16S_Relative_Abundance(Family).png)

### Family-Level Turnover & Dominance Breakdown
Agglomerating ASVs to the Family level reveals that rhizosphere remodeling is driven by cumulative shifts across dozens of specialized families rather than monoculture blooms:
* **Top Identified Families:**
  * **`Chitinophagaceae`:** Consistently the most abundant family ($5.2\%\text{--}6.1\%$).
  * **`Sphingobacteriaceae`:** Accounts for $4.1\%\text{--}6.2\%$, showing highest representation in `GJ_Control` ($6.2\%$).
  * **`Bacillaceae`:** Contracts from $4.8\%$ in controls to $3.1\%$ in treatments.
  * **`Opitutaceae`:** Stable across conditions at $3.2\%\text{--}3.8\%$.
  * **`Sphingomonadaceae`:** Accounts for $3.2\%\text{--}4.4\%$, enriched in `GJ_Treatment` ($4.4\%$).
  * **Secondary Families:** `Caulobacteraceae` ($2.4\%\text{--}2.9\%$), `Micropepsaceae` ($2.1\%\text{--}2.7\%$), `Flavobacteriaceae` ($1.8\%\text{--}2.3\%$), `Oxalobacteraceae` ($1.6\%\text{--}2.1\%$), `Nocardioidaceae` ($1.4\%\text{--}1.9\%$), `Rhizobiaceae` ($1.2\%\text{--}1.6\%$), and `Xanthobacteraceae` ($1.1\%\text{--}1.5\%$).

### Quantitative Dominance Dilution
The expansion of unclassified and low-abundance families highlights the core ecological mechanism of soil suspension inoculation:
* In `GB_Control` and `GJ_Control`, the top 12 families represent **$37.1\%$ and $42.8\%$** of total community abundance, leaving the `Other` fraction at $62.9\%$ and $57.2\%$.
* In `GB_Treatment` and `GJ_Treatment`, the cumulative abundance of these top 12 families contracts to **$32.4\%$ and $35.1\%$**, driving the `Other` fraction up to **$67.6\%$ and $64.9\%$**.

---

## 6. Genus-Level Resolution & Fine-Scale Rare Biosphere Influx

Agglomerating counts to the Genus rank provides the finest taxonomic resolution achievable from short-read 16S V3–V4 amplicons, isolating candidate genera responsible for treatment divergence.

![Global Genus-Level Relative Abundance](../assets/075_16S_Relative_Abundance(Genus).png)

### Global Four-Group Genus Patterns
Evaluating the top 12 genera across all 110 libraries illustrates the high background diversity of the rhizosphere:
* **Uncultured Reservoir:** Genera lacking cultured representatives (`uncultured`) constitute the single largest assigned fraction, occupying $10.4\%\text{--}13.5\%$ across all treatments.
* **Core Colonizers:** Key genera present across all conditions include *Mucilaginibacter* ($2.3\%\text{--}3.5\%$), *Flavobacterium* ($2.1\%\text{--}2.8\%$), *Nocardioides* ($1.6\%\text{--}2.4\%$), *Massilia* ($1.7\%\text{--}2.2\%$), *WD2101_soil_group* ($1.5\%\text{--}2.1\%$), *LWQ8* ($1.3\%\text{--}1.8\%$), *Asticcacaulis* ($1.2\%\text{--}1.6\%$), *BIrii41* ($1.1\%\text{--}1.5\%$), *Streptomyces* ($1.0\%\text{--}1.5\%$), and *Candidatus Kaiserbacteria* ($0.9\%\text{--}1.4\%$).
* **Peat Resident Suppression (*Bacillus*):** Corroborating Bacillales trends, *Bacillus* accounts for $4.9\%$ in `GB_Control` and $3.8\%$ in `GJ_Control`, but drops significantly to $2.1\%$ in `GB_Treatment` and $1.4\%$ in `GJ_Treatment`.
* **Global Dominance Contraction:** The top 12 genera account for **$32.9\%$ in `GB_Control`** and **$33.7\%$ in `GJ_Control`**. Following soil inoculation, they contract to **$29.6\%$ in `GB_Treatment`** and **$30.4\%$ in `GJ_Treatment`**, leaving approximately **$70\%$ of the community to the unassigned `Other` pool**.

---

### Expanded Lineage-Specific Genus Succession (Top ~30 Genera)

To decouple lineage-specific recruitment from broad suppression trends, genus profiles were expanded to encompass the top ~30 genera within each regional contrast.

![GB Lineage Genus-Level Relative Abundance](../assets/076_16S_GB_Relative_Abundance(Genus).png)

![GJ Lineage Genus-Level Relative Abundance](../assets/077_16S_GJ_Relative_Abundance(Genus).png)

#### 1. Gijang B Genus Dynamics (Figure 076)
* **Selective Influx:** Soil inoculation promotes specific plant growth-promoting and oligotrophic clades:
  * *WD2101_soil_group* ($1.8\% \to 2.2\%$), *Parcubacteria* ($0.9\% \to 1.5\%$), and *Chitinophaga* ($1.1\% \to 1.4\%$).
  * Emergence of low-abundance field clades: *Candidatus Adlerbacteria*, *Burkholderia-Caballeronia-Paraburkholderia*, *Sphingobium*, and *Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium*.
* **Cumulative Metric:** The expanded top 30 genera account for $49.2\%$ of total reads in `GB_Control` and contract to $45.6\%$ in `GB_Treatment`, expanding the `Other` tail from $50.8\%$ to **$54.4\%$**.

#### 2. Gyeongju Genus Dynamics (Figure 077)
* **Lineage-Specific Turnover:** Gyeongju inoculation remodels a distinct set of genera:
  * Suppression of baseline peat copiotrophs: *Mucilaginibacter* contracts from $4.8\%$ to $2.8\%$; *Flavobacterium* drops from $3.2\%$ to $2.2\%$; *Nocardioides* drops from $2.8\%$ to $1.8\%$.
  * Enrichment of environmental soil genera: *Sphingobium* ($1.4\% \to 1.9\%$), *WD2101_soil_group* ($1.5\% \to 1.9\%$), *Gemmata* ($1.2\% \to 1.6\%$), and *Rugosimonospora*.
  * Exclusive detection among top taxa: *Abditibacterium*, *SH-PL14*, *Schlesneria*, *Gemmatimonas*, and *Taibaiella*.
* **Pronounced Dominance Dilution:** In `GJ_Control`, the top 30 genera represent $53.6\%$ of the library (`Other` = $46.4\%$). In `GJ_Treatment`, the top 30 genera contract sharply to $46.1\%$, pushing the rare biosphere (`Other`) to **$53.9\%$**.

---

## 7. Synthesis: The Rare Biosphere Mechanism of Soil Conditioning

Tracking community succession across Phylum, Class, Order, Family, and Genus levels reveals a consistent ecological mechanism:

$$\text{Peat Baseline (Low Diversity, High Dominance)} \xrightarrow{+\text{ Soil Suspension}} \text{Rhizosphere State (High Diversity, High Evenness)}$$

1. **Absence of Monoculture Blooms:** Inoculation does not trigger runaway dominance by any single opportunistic genus. The relative abundance of the top genus (*uncultured*) remains steady at $\sim 11\text{--}13\%$, while common copiotrophic peat colonizers (*Bacillus*, *Mucilaginibacter*) contract.
2. **Cumulative Rare-Biosphere Expansion:** At every taxonomic tier, the proportion of reads assigned to the dominant top clades declines, accompanied by a reciprocal expansion of the `Other` category ($+3\%\text{--}+8\%$ depending on rank).
3. **Lineage-Specific Taxa Recruitment:** The species pools differentiating Gijang B (accelerated ripening) from Gyeongju (increased vegetative biomass) reside primarily within low-abundance, specialized genera (*Parcubacteria*, *Acidibacter*, *Burkholderia* in GB; *Sphingobium*, *Gemmata*, *Porphyrobacter*, *Gemmatimonas* in GJ). These taxa stably integrate into the root surface, driving overall community divergence and establishing the distinct microbial equilibria evaluated in downstream functional profiling.