# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Script: R Environment Configuration & Dependency Installation
# Target Assays: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1, GB2, GB3, GJ1, GJ2; N=110)
# ==============================================================================

############################### R Studio Environment Set-up ###############################

# 1. CRAN Package Installation
cran_packages <- c(
  "tidyverse",
  "vegan",
  "ggpubr",
  "patchwork",
  "devtools",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "rstatix",
  "MicrobiomeStat",
  "ggpicrust2"
)

installed_cran <- rownames(installed.packages())
for (pkg in cran_packages) {
  if (!pkg %in% installed_cran) {
    install.packages(pkg)
  }
}

# 2. Bioconductor Manager & Core Packages Installation
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "phyloseq",
  "ANCOMBC",
  "ALDEx2",
  "microbiome",
  "TreeSummarizedExperiment",
  "mia",
  "ComplexHeatmap"
)

installed_bioc <- rownames(installed.packages())
for (pkg in bioc_packages) {
  if (!pkg %in% installed_bioc) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# 3. GitHub Package Installation (qiime2R for QIIME 2 Artifact Ingestion)
if (!requireNamespace("qiime2R", quietly = TRUE)) {
  devtools::install_github("jbisanz/qiime2R")
}

# 4. Working Directory Setup
setwd("/Users/heechan/Desktop/miseq_data_without_GJSet4")






############################### 1-1_16S_4Group_Alpha_Diversity ###############################
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggpubr)
library(rstatix)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter 4 Target Groups
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
    )
  )

# Color Palette for GB_Control, GB_Treatment, GJ_Control, GJ_Treatment
palette_4 <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4",
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 2. Import Artifacts from 07_diversity_16S/core_metrics/
# ==============================================================================
obs_features <- read_qza("07_diversity_16S/core_metrics/observed_features_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Observed_Features = observed_features)

shannon <- read_qza("07_diversity_16S/core_metrics/shannon_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Shannon = shannon_entropy)

faith_pd <- read_qza("07_diversity_16S/core_metrics/faith_pd_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Faith_PD = faith_pd)

# Calculate Simpson Index (1 - D) from Rarefied ASV Table
rarefied_table <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
asv_matrix <- t(rarefied_table)
simpson_vec <- diversity(asv_matrix, index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

# ==============================================================================
# 3. Merge & Format Alpha Diversity Metrics
# ==============================================================================
alpha_merged <- metadata %>%
  inner_join(obs_features, by = "sample-id") %>%
  inner_join(shannon, by = "sample-id") %>%
  inner_join(simpson_vec, by = "sample-id") %>%
  inner_join(faith_pd, by = "sample-id")

alpha_long <- alpha_merged %>%
  pivot_longer(
    cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"))
  )

# ==============================================================================
# 4. Statistical Testing (Kruskal-Wallis & Post-hoc Pairwise Wilcoxon FDR)
# ==============================================================================
kw_results <- alpha_long %>%
  group_by(Metric) %>%
  kruskal_test(Value ~ Treatment_Group)

pairwise_results <- alpha_long %>%
  group_by(Metric) %>%
  wilcox_test(Value ~ Treatment_Group, p.adjust.method = "fdr")

cat("=== Global Kruskal-Wallis Test Results ===\n")
print(kw_results)

cat("\n=== Pairwise Wilcoxon Post-Hoc (FDR Adjusted) ===\n")
print(pairwise_results)

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_alpha <- ggplot(alpha_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_4, name = "Treatment Group") +
  scale_color_manual(values = palette_4, name = "Treatment Group") +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment",
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  stat_compare_means(
    method = "kruskal.test",
    label.y.npc = "top",
    size = 3.8,
    family = font_family
  ) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = font_family, size = 10),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", color = "black", family = font_family, size = 10),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", size = 11, family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(
    title = "1. M(B), B, M(J), J - Alpha Diversity",
    y = "Calculated Diversity Value"
  )

# Display Plot
print(p_alpha)










############################### 1-2_16S_4Group_Beta_Diversity ###############################
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggpubr)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter 4 Target Groups
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
    )
  )

# Color Palette for 4 Groups
palette_4 <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4",
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 2. Function for PERMANOVA, Beta-Dispersion, and PCoA Ordination
# ==============================================================================
analyze_beta_diversity <- function(dist_qza_path, metric_name, metadata_df) {
  # A. Import Distance Matrix .qza
  raw_dist <- read_qza(dist_qza_path)$data
  dist_mat <- as.matrix(raw_dist)

  # B. Align Distance Matrix with Metadata
  common_samples <- intersect(metadata_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- metadata_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  # C. Statistical Testing: PERMANOVA (adonis2)
  set.seed(42)
  permanova_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  r2 <- permanova_res$R2[1]
  p_val <- permanova_res$`Pr(>F)`[1]

  # D. Statistical Testing: Multivariate Dispersion (betadisper)
  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  cat(sprintf("\n=== %s Distance Matrix Statistical Results ===\n", metric_name))
  cat(sprintf("PERMANOVA R²: %.4f | p-value: %.4f\n", r2, p_val))
  cat(sprintf("Beta-Dispersion p-value: %.4f\n", disp_p))

  # E. PCoA Ordination Computation
  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_explained <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  # F. Text Annotation Formatting
  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f",
    r2,
    ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    disp_p
  )

  # G. Subplot Generation
  p_ord <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(values = palette_4, name = "Treatment Group") +
    scale_fill_manual(values = palette_4, name = "Treatment Group") +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      text = element_text(family = font_family),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family, size = 10),
      axis.text = element_text(color = "black", family = font_family),
      axis.title = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate(
      "text", x = -Inf, y = Inf, label = annot_text,
      hjust = -0.05, vjust = 1.1, size = 3.5, fontface = "bold", family = font_family
    ) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2])
    )

  return(p_ord)
}

# ==============================================================================
# 3. Execute Beta Diversity Analysis (core_metrics directory)
# ==============================================================================
res_bc  <- analyze_beta_diversity("07_diversity_16S/core_metrics/bray_curtis_distance_matrix.qza", "Bray-Curtis", metadata)
res_wuf <- analyze_beta_diversity("07_diversity_16S/core_metrics/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", metadata)
res_uuf <- analyze_beta_diversity("07_diversity_16S/core_metrics/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", metadata)
res_jac <- analyze_beta_diversity("07_diversity_16S/core_metrics/jaccard_distance_matrix.qza", "Jaccard", metadata)

# ==============================================================================
# 4. Combine & Export Panel Plot
# ==============================================================================
combined_pcoa <- ggarrange(
  res_bc, res_wuf, res_uuf, res_jac,
  ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"
)

combined_pcoa_annotated <- annotate_figure(
  combined_pcoa,
  top = text_grob("1. M(B), B, M(J), J - Beta Diversity", face = "bold", size = 14, family = font_family)
)

# Display
print(combined_pcoa_annotated)









############################### 1-3_16S_4Group_Beta_Diversity(NMDS) ###############################
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggpubr)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter 4 Target Groups
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
    )
  )

# Color Palette for 4 Groups
palette_4 <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4",
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 2. Function for PERMANOVA, Beta-Dispersion, and NMDS Ordination
# ==============================================================================
analyze_beta_nmds <- function(dist_qza_path, metric_name, metadata_df) {
  # A. Import Distance Matrix .qza
  raw_dist <- read_qza(dist_qza_path)$data
  dist_mat <- as.matrix(raw_dist)

  # B. Align Distance Matrix with Metadata
  common_samples <- intersect(metadata_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- metadata_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  # C. Statistical Testing: PERMANOVA (adonis2)
  set.seed(42)
  permanova_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  r2 <- permanova_res$R2[1]
  p_val <- permanova_res$`Pr(>F)`[1]

  # D. Statistical Testing: Multivariate Dispersion (betadisper)
  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  # E. NMDS Computation (metaMDS)
  set.seed(42)
  nmds_res <- metaMDS(dist_sub, k = 2, trymax = 100, trace = FALSE)
  stress_val <- nmds_res$stress

  cat(sprintf("\n=== %s NMDS Statistical Results ===\n", metric_name))
  cat(sprintf("PERMANOVA R²: %.4f | p-value: %.4f\n", r2, p_val))
  cat(sprintf("Beta-Dispersion p-value: %.4f\n", disp_p))
  cat(sprintf("NMDS Stress: %.4f\n", stress_val))

  nmds_df <- as.data.frame(nmds_res$points) %>%
    rownames_to_column("sample-id") %>%
    rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
    inner_join(meta_sub, by = "sample-id")

  # F. Text Annotation Formatting
  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f\nStress = %.3f",
    r2,
    ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    disp_p,
    stress_val
  )

  # G. Subplot Generation
  p_ord <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(values = palette_4, name = "Treatment Group") +
    scale_fill_manual(values = palette_4, name = "Treatment Group") +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      text = element_text(family = font_family),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family, size = 10),
      axis.text = element_text(color = "black", family = font_family),
      axis.title = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate(
      "text", x = -Inf, y = Inf, label = annot_text,
      hjust = -0.05, vjust = 1.1, size = 3.3, fontface = "bold", family = font_family
    ) +
    labs(
      title = sprintf("NMDS - %s Distance", metric_name),
      x = "NMDS 1",
      y = "NMDS 2"
    )

  return(p_ord)
}

# ==============================================================================
# 3. Execute Beta Diversity Analysis (core_metrics directory)
# ==============================================================================
res_bc  <- analyze_beta_nmds("07_diversity_16S/core_metrics/bray_curtis_distance_matrix.qza", "Bray-Curtis", metadata)
res_wuf <- analyze_beta_nmds("07_diversity_16S/core_metrics/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", metadata)
res_uuf <- analyze_beta_nmds("07_diversity_16S/core_metrics/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", metadata)
res_jac <- analyze_beta_nmds("07_diversity_16S/core_metrics/jaccard_distance_matrix.qza", "Jaccard", metadata)

# ==============================================================================
# 4. Combine & Export Panel Plot
# ==============================================================================
combined_nmds <- ggarrange(
  res_bc, res_wuf, res_uuf, res_jac,
  ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"
)

combined_nmds_annotated <- annotate_figure(
  combined_nmds,
  top = text_grob("1. M(B), B, M(J), J - Beta Diversity", face = "bold", size = 14, family = font_family)
)

# Display
print(combined_nmds_annotated)










############################### 2-1_16S_Treatment_Beta_Diversity ###############################
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggpubr)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter Treatment Groups (GB_Treatment vs GJ_Treatment)
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Treatment", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Treatment", "GJ_Treatment")
    )
  )

# Color Palette for Treatment Groups
palette_2 <- c(
  "GB_Treatment" = "#41B6C4",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 2. Function for PERMANOVA, Beta-Dispersion, and PCoA Ordination
# ==============================================================================
analyze_beta_diversity <- function(dist_qza_path, metric_name, metadata_df) {
  # A. Import Distance Matrix .qza
  raw_dist <- read_qza(dist_qza_path)$data
  dist_mat <- as.matrix(raw_dist)

  # B. Align Distance Matrix with Filtered Metadata
  common_samples <- intersect(metadata_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- metadata_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  # C. Statistical Testing: PERMANOVA (adonis2)
  set.seed(42)
  permanova_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  r2 <- permanova_res$R2[1]
  p_val <- permanova_res$`Pr(>F)`[1]

  # D. Statistical Testing: Multivariate Dispersion (betadisper)
  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  cat(sprintf("\n=== %s Distance Matrix (B vs J) Statistical Results ===\n", metric_name))
  cat(sprintf("PERMANOVA R²: %.4f | p-value: %.4f\n", r2, p_val))
  cat(sprintf("Beta-Dispersion p-value: %.4f\n", disp_p))

  # E. PCoA Computation
  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_explained <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  # F. Formatting Statistical Text Box
  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f",
    r2,
    ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    disp_p
  )

  # G. Subplot Generation
  p_ord <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(values = palette_2, name = "Treatment Group") +
    scale_fill_manual(values = palette_2, name = "Treatment Group") +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      text = element_text(family = font_family),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family, size = 10),
      axis.text = element_text(color = "black", family = font_family),
      axis.title = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate(
      "text", x = -Inf, y = Inf, label = annot_text,
      hjust = -0.05, vjust = 1.1, size = 3.5, fontface = "bold", family = font_family
    ) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2])
    )

  return(p_ord)
}

# ==============================================================================
# 3. Execute Analysis across Core Distance Metrics (core_metrics directory)
# ==============================================================================
res_bc  <- analyze_beta_diversity("07_diversity_16S/core_metrics/bray_curtis_distance_matrix.qza", "Bray-Curtis", metadata)
res_wuf <- analyze_beta_diversity("07_diversity_16S/core_metrics/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", metadata)
res_uuf <- analyze_beta_diversity("07_diversity_16S/core_metrics/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", metadata)
res_jac <- analyze_beta_diversity("07_diversity_16S/core_metrics/jaccard_distance_matrix.qza", "Jaccard", metadata)

# ==============================================================================
# 4. Combine & Export Panel Plot
# ==============================================================================
combined_pcoa <- ggarrange(
  res_bc, res_wuf, res_uuf, res_jac,
  ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"
)

combined_pcoa_annotated <- annotate_figure(
  combined_pcoa,
  top = text_grob("2. B, J - Beta Diversity", face = "bold", size = 14, family = font_family)
)

# Display
print(combined_pcoa_annotated)










############################### 2-2_16S_Treatment_Beta_Diversity(NMDS) ###############################
library(qiime2R)
library(tidyverse)
library(vegan)
library(ggpubr)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter Treatment Groups (GB_Treatment vs GJ_Treatment)
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Treatment", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Treatment", "GJ_Treatment")
    )
  )

# Color Palette for Treatment Groups
palette_2 <- c(
  "GB_Treatment" = "#41B6C4",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 2. Function for PERMANOVA, Beta-Dispersion, and NMDS Ordination
# ==============================================================================
analyze_beta_nmds <- function(dist_qza_path, metric_name, metadata_df) {
  # A. Import Distance Matrix .qza
  raw_dist <- read_qza(dist_qza_path)$data
  dist_mat <- as.matrix(raw_dist)

  # B. Align Distance Matrix with Filtered Metadata
  common_samples <- intersect(metadata_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- metadata_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  # C. Statistical Testing: PERMANOVA (adonis2)
  set.seed(42)
  permanova_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  r2 <- permanova_res$R2[1]
  p_val <- permanova_res$`Pr(>F)`[1]

  # D. Statistical Testing: Multivariate Dispersion (betadisper)
  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  # E. NMDS Computation (metaMDS)
  set.seed(42)
  nmds_res <- metaMDS(dist_sub, k = 2, trymax = 100, trace = FALSE)
  stress_val <- nmds_res$stress

  cat(sprintf("\n=== %s NMDS (B vs J) Statistical Results ===\n", metric_name))
  cat(sprintf("PERMANOVA R²: %.4f | p-value: %.4f\n", r2, p_val))
  cat(sprintf("Beta-Dispersion p-value: %.4f\n", disp_p))
  cat(sprintf("NMDS Stress: %.4f\n", stress_val))

  nmds_df <- as.data.frame(nmds_res$points) %>%
    rownames_to_column("sample-id") %>%
    rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
    inner_join(meta_sub, by = "sample-id")

  # F. Formatting Statistical Text Box
  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f\nStress = %.3f",
    r2,
    ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    disp_p,
    stress_val
  )

  # G. Subplot Generation
  p_ord <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = 2.8, alpha = 0.85) +
    scale_color_manual(values = palette_2, name = "Treatment Group") +
    scale_fill_manual(values = palette_2, name = "Treatment Group") +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      text = element_text(family = font_family),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family, size = 10),
      axis.text = element_text(color = "black", family = font_family),
      axis.title = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate(
      "text", x = -Inf, y = Inf, label = annot_text,
      hjust = -0.05, vjust = 1.1, size = 3.3, fontface = "bold", family = font_family
    ) +
    labs(
      title = sprintf("NMDS - %s Distance", metric_name),
      x = "NMDS 1",
      y = "NMDS 2"
    )

  return(p_ord)
}

# ==============================================================================
# 3. Execute Analysis across Core Distance Metrics (core_metrics directory)
# ==============================================================================
res_bc  <- analyze_beta_nmds("07_diversity_16S/core_metrics/bray_curtis_distance_matrix.qza", "Bray-Curtis", metadata)
res_wuf <- analyze_beta_nmds("07_diversity_16S/core_metrics/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", metadata)
res_uuf <- analyze_beta_nmds("07_diversity_16S/core_metrics/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", metadata)
res_jac <- analyze_beta_nmds("07_diversity_16S/core_metrics/jaccard_distance_matrix.qza", "Jaccard", metadata)

# ==============================================================================
# 4. Combine & Export Panel Plot
# ==============================================================================
combined_nmds <- ggarrange(
  res_bc, res_wuf, res_uuf, res_jac,
  ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"
)

combined_nmds_annotated <- annotate_figure(
  combined_nmds,
  top = text_grob("2. B, J - Beta Diversity", face = "bold", size = 14, family = font_family)
)

# Display Plot
print(combined_nmds_annotated)



############################### Centralized Data Ingestion & Relative Abundance ###############################
library(qiime2R)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)

# Global font configuration
font_family <- "Times New Roman"

# Load Metadata & Filter 4 Target Groups (N = 110)
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
    )
  )

# ==============================================================================
# 1. Import Master Count Table & Taxonomy Artifacts
# ==============================================================================
counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw    <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

# Align samples between count matrix and metadata (N = 110)
common_samples <- intersect(metadata$`sample-id`, colnames(counts_raw))
counts_sub     <- counts_raw[, common_samples]

# ==============================================================================
# 2. Compute Master Long-Format Relative Abundance (%)
# ==============================================================================
rel_abund_long <- counts_sub %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  left_join(tax_clean, by = "Feature.ID") %>%
  inner_join(metadata, by = "sample-id")



############################### 3-1_16S_Relative_Abundance(Phylum) ###############################

# ==============================================================================
# 3. Phylum-Level Aggregation (Two-Step Summary for 100% Group Scaling)
# ==============================================================================
phylum_sample <- rel_abund_long %>%
  mutate(Phylum = if_else(is.na(Phylum) | Phylum == "" | Phylum == "p__", "Unassigned", Phylum)) %>%
  group_by(`sample-id`, Treatment_Group, Phylum) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# Identify Top 10 Phyla by global mean relative abundance
top10_phyla <- phylum_sample %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Phylum != "Unassigned") %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# STEP 1: Sum Top 10 + "Other" PER SAMPLE
phylum_sample_grouped <- phylum_sample %>%
  mutate(Phylum_Group = if_else(Phylum %in% top10_phyla, Phylum, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Phylum_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across Treatment Groups
phylum_group_summary <- phylum_sample_grouped %>%
  group_by(Treatment_Group, Phylum_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Set factor levels for stacking order (Top taxa on top, "Other" at the base)
phylum_levels <- c(top10_phyla, "Other")
phylum_group_summary <- phylum_group_summary %>%
  mutate(Phylum_Group = factor(Phylum_Group, levels = rev(phylum_levels)))

# Color palette for 10 taxa + 1 grey for Other
taxa_colors_phylum <- c(brewer.pal(10, "Paired"), "#B0B0B0")
names(taxa_colors_phylum) <- phylum_levels

# ==============================================================================
# 4. Visualization
# ==============================================================================
p_phylum <- ggplot(phylum_group_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Phylum_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_phylum) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment",
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. Relative Abundance (Phylum Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_phylum)










############################### 3-1_16S_GB_Relative_Abundance(Phylum) ###############################

# ==============================================================================
# 3. Cohort-Specific Top 10 Identification & Master Color Dictionary
# ==============================================================================
# Top 10 Phyla in GB cohort
top10_gb <- phylum_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Phylum != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Top 10 Phyla in GJ cohort (for cross-cohort color consistency)
top10_gj <- phylum_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Phylum != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Union of all top phyla across GB and GJ cohorts
all_top_phyla <- unique(c(top10_gb, top10_gj))

# Build Master Named Palette (Ensures shared phyla have identical hex codes)
palette_pool <- colorRampPalette(brewer.pal(12, "Paired"))(length(all_top_phyla))
master_taxa_colors <- c(setNames(palette_pool, all_top_phyla), "Other" = "#B0B0B0")

# ==============================================================================
# 4. GB Phylum-Level Aggregation (Two-Step Summary)
# ==============================================================================
# Filter dataset down to GB samples
phylum_sample_gb <- phylum_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment"))

# STEP 1: Collapse non-top 10 GB phyla into "Other"
phylum_sample_gb_grouped <- phylum_sample_gb %>%
  mutate(Phylum_Group = if_else(Phylum %in% top10_gb, Phylum, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Phylum_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across GB Groups
phylum_gb_summary <- phylum_sample_gb_grouped %>%
  group_by(Treatment_Group, Phylum_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Stacking order (Top GB taxa on top, "Other" at base)
gb_phylum_levels <- c(top10_gb, "Other")
phylum_gb_summary <- phylum_gb_summary %>%
  mutate(Phylum_Group = factor(Phylum_Group, levels = rev(gb_phylum_levels)))

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_phylum_gb <- ggplot(phylum_gb_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Phylum_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. GB Relative Abundance (Phylum Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_phylum_gb)










############################### 3-1_16S_GJ_Relative_Abundance(Phylum) ###############################

# ==============================================================================
# 3. Cohort-Specific Top 10 Identification & Master Color Dictionary
# ==============================================================================
# Top 10 Phyla in GB cohort
top10_gb <- phylum_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Phylum != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Top 10 Phyla in GJ cohort
top10_gj <- phylum_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Phylum != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

# Union of all top phyla across GB and GJ cohorts
all_top_phyla <- unique(c(top10_gb, top10_gj))

# Build Master Named Palette (Shared phyla retain identical hex codes)
palette_pool <- colorRampPalette(brewer.pal(12, "Paired"))(length(all_top_phyla))
master_taxa_colors <- c(setNames(palette_pool, all_top_phyla), "Other" = "#B0B0B0")

# ==============================================================================
# 4. GJ Phylum-Level Aggregation (Two-Step Summary)
# ==============================================================================
# Filter dataset down to GJ samples
phylum_sample_gj <- phylum_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment"))

# STEP 1: Collapse non-top 10 GJ phyla into "Other"
phylum_sample_gj_grouped <- phylum_sample_gj %>%
  mutate(Phylum_Group = if_else(Phylum %in% top10_gj, Phylum, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Phylum_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across GJ Groups
phylum_gj_summary <- phylum_sample_gj_grouped %>%
  group_by(Treatment_Group, Phylum_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Stacking order (Top GJ taxa on top, "Other" at base)
gj_phylum_levels <- c(top10_gj, "Other")
phylum_gj_summary <- phylum_gj_summary %>%
  mutate(Phylum_Group = factor(Phylum_Group, levels = rev(gj_phylum_levels)))

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_phylum_gj <- ggplot(phylum_gj_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Phylum_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors) +
  scale_x_discrete(labels = c(
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. GJ Relative Abundance (Phylum Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_phylum_gj)










############################### 3-2_16S_Relative_Abundance(Class) ###############################

# ==============================================================================
# 3. Class-Level Aggregation (Two-Step Summary for 100% Group Scaling)
# ==============================================================================
class_sample <- rel_abund_long %>%
  mutate(Class = if_else(is.na(Class) | Class == "" | Class == "c__", "Unassigned", Class)) %>%
  group_by(`sample-id`, Treatment_Group, Class) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# Identify Top 12 Classes by global mean relative abundance
top12_classes <- class_sample %>%
  group_by(Class) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Class != "Unassigned") %>%
  slice_head(n = 12) %>%
  pull(Class)

# STEP 1: Sum Top 12 + "Other" PER SAMPLE
class_sample_grouped <- class_sample %>%
  mutate(Class_Group = if_else(Class %in% top12_classes, Class, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Class_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across Treatment Groups
class_group_summary <- class_sample_grouped %>%
  group_by(Treatment_Group, Class_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Set factor levels for stacking order (Top taxa on top, "Other" at the base)
class_levels <- c(top12_classes, "Other")
class_group_summary <- class_group_summary %>%
  mutate(Class_Group = factor(Class_Group, levels = rev(class_levels)))

# Color palette for 12 taxa + 1 grey for Other
taxa_colors_class <- c(brewer.pal(12, "Set3"), "#B0B0B0")
names(taxa_colors_class) <- class_levels

# ==============================================================================
# 4. Visualization
# ==============================================================================
p_class <- ggplot(class_group_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Class_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_class) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment",
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. Relative Abundance (Class Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_class)










############################### 3-3_16S_Relative_Abundance(Order) ###############################

# ==============================================================================
# 3. Order-Level Aggregation (Two-Step Summary for 100% Group Scaling)
# ==============================================================================
order_sample <- rel_abund_long %>%
  mutate(Order = if_else(is.na(Order) | Order == "" | Order == "o__", "Unassigned", Order)) %>%
  group_by(`sample-id`, Treatment_Group, Order) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# Identify Top 12 Orders by global mean relative abundance
top12_orders <- order_sample %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Order != "Unassigned") %>%
  slice_head(n = 12) %>%
  pull(Order)

# STEP 1: Sum Top 12 + "Other" PER SAMPLE
order_sample_grouped <- order_sample %>%
  mutate(Order_Group = if_else(Order %in% top12_orders, Order, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Order_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across Treatment Groups
order_group_summary <- order_sample_grouped %>%
  group_by(Treatment_Group, Order_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Set factor levels for stacking order (Top taxa on top, "Other" at the base)
order_levels <- c(top12_orders, "Other")
order_group_summary <- order_group_summary %>%
  mutate(Order_Group = factor(Order_Group, levels = rev(order_levels)))

# Color palette for 12 taxa + 1 grey for Other
taxa_colors_order <- c(brewer.pal(12, "Set3"), "#B0B0B0")
names(taxa_colors_order) <- order_levels

# ==============================================================================
# 4. Visualization
# ==============================================================================
p_order <- ggplot(order_group_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Order_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_order) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment",
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. Relative Abundance (Order Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_order)










############################### 3-3_16S_GB_Relative_Abundance(Order) ###############################

# ==============================================================================
# 3. Cohort-Specific Top 12 Identification & Master Color Dictionary
# ==============================================================================
# Top 12 Orders in GB cohort
top12_gb <- order_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Order != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Order)

# Top 12 Orders in GJ cohort (for cross-cohort color consistency)
top12_gj <- order_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Order != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Order)

# Union of all top orders across GB and GJ cohorts
all_top_orders <- unique(c(top12_gb, top12_gj))

# Build Master Named Palette (Shared orders retain identical hex codes)
palette_pool <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_top_orders))
master_taxa_colors <- c(setNames(palette_pool, all_top_orders), "Other" = "#B0B0B0")

# ==============================================================================
# 4. GB Order-Level Aggregation (Two-Step Summary)
# ==============================================================================
# Filter dataset down to GB samples
order_sample_gb <- order_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment"))

# STEP 1: Collapse non-top 12 GB orders into "Other"
order_sample_gb_grouped <- order_sample_gb %>%
  mutate(Order_Group = if_else(Order %in% top12_gb, Order, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Order_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across GB Groups
order_gb_summary <- order_sample_gb_grouped %>%
  group_by(Treatment_Group, Order_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Stacking order (Top GB taxa on top, "Other" at base)
gb_order_levels <- c(top12_gb, "Other")
order_gb_summary <- order_gb_summary %>%
  mutate(Order_Group = factor(Order_Group, levels = rev(gb_order_levels)))

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_order_gb <- ggplot(order_gb_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Order_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. GB Relative Abundance (Order Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_order_gb)









############################### 3-3_16S_GJ_Relative_Abundance(Order) ###############################

# ==============================================================================
# 3. Cohort-Specific Top 12 Identification & Master Color Dictionary
# ==============================================================================
# Top 12 Orders in GB cohort
top12_gb <- order_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Order != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Order)

# Top 12 Orders in GJ cohort
top12_gj <- order_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Order != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Order)

# Union of all top orders across GB and GJ cohorts
all_top_orders <- unique(c(top12_gb, top12_gj))

# Build Master Named Palette (Shared orders retain identical hex codes)
palette_pool <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_top_orders))
master_taxa_colors <- c(setNames(palette_pool, all_top_orders), "Other" = "#B0B0B0")

# ==============================================================================
# 4. GJ Order-Level Aggregation (Two-Step Summary)
# ==============================================================================
# Filter dataset down to GJ samples
order_sample_gj <- order_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment"))

# STEP 1: Collapse non-top 12 GJ orders into "Other"
order_sample_gj_grouped <- order_sample_gj %>%
  mutate(Order_Group = if_else(Order %in% top12_gj, Order, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Order_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across GJ Groups
order_gj_summary <- order_sample_gj_grouped %>%
  group_by(Treatment_Group, Order_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Stacking order (Top GJ taxa on top, "Other" at base)
gj_order_levels <- c(top12_gj, "Other")
order_gj_summary <- order_gj_summary %>%
  mutate(Order_Group = factor(Order_Group, levels = rev(gj_order_levels)))

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_order_gj <- ggplot(order_gj_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Order_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors) +
  scale_x_discrete(labels = c(
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. GJ Relative Abundance (Order Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_order_gj)











############################### 3-4_16S_Relative_Abundance(Family) ###############################

# ==============================================================================
# 3. Family-Level Aggregation (Two-Step Summary for 100% Group Scaling)
# ==============================================================================
family_sample <- rel_abund_long %>%
  mutate(Family = if_else(is.na(Family) | Family == "" | Family == "f__", "Unassigned", Family)) %>%
  group_by(`sample-id`, Treatment_Group, Family) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# Identify Top 12 Families by global mean relative abundance
top12_families <- family_sample %>%
  group_by(Family) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Family != "Unassigned") %>%
  slice_head(n = 12) %>%
  pull(Family)

# STEP 1: Sum Top 12 + "Other" PER SAMPLE
family_sample_grouped <- family_sample %>%
  mutate(Family_Group = if_else(Family %in% top12_families, Family, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Family_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across Treatment Groups
family_group_summary <- family_sample_grouped %>%
  group_by(Treatment_Group, Family_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Set factor levels for stacking order (Top taxa on top, "Other" at the base)
family_levels <- c(top12_families, "Other")
family_group_summary <- family_group_summary %>%
  mutate(Family_Group = factor(Family_Group, levels = rev(family_levels)))

# Color palette for 12 taxa + 1 grey for Other
taxa_colors_family <- c(brewer.pal(12, "Set3"), "#B0B0B0")
names(taxa_colors_family) <- family_levels

# ==============================================================================
# 4. Visualization
# ==============================================================================
p_family <- ggplot(family_group_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Family_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_family) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment",
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. Relative Abundance (Family Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_family)










############################### 3-5_16S_Relative_Abundance(Genus) ###############################

# ==============================================================================
# 3. Genus-Level Aggregation (Two-Step Summary for 100% Group Scaling)
# ==============================================================================
genus_sample <- rel_abund_long %>%
  mutate(Genus = if_else(is.na(Genus) | Genus == "" | Genus == "g__", "Unassigned", Genus)) %>%
  group_by(`sample-id`, Treatment_Group, Genus) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# Identify Top 12 Genera by global mean relative abundance
top12_genera <- genus_sample %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Genus != "Unassigned") %>%
  slice_head(n = 12) %>%
  pull(Genus)

# STEP 1: Sum Top 12 + "Other" PER SAMPLE
genus_sample_grouped <- genus_sample %>%
  mutate(Genus_Group = if_else(Genus %in% top12_genera, Genus, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Genus_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across Treatment Groups
genus_group_summary <- genus_sample_grouped %>%
  group_by(Treatment_Group, Genus_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Set factor levels for stacking order (Top taxa on top, "Other" at the base)
genus_levels <- c(top12_genera, "Other")
genus_group_summary <- genus_group_summary %>%
  mutate(Genus_Group = factor(Genus_Group, levels = rev(genus_levels)))

# Color palette for 12 taxa + 1 grey for Other
taxa_colors_genus <- c(brewer.pal(12, "Set3"), "#B0B0B0")
names(taxa_colors_genus) <- genus_levels

# ==============================================================================
# 4. Visualization
# ==============================================================================
p_genus <- ggplot(genus_group_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Genus_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_genus) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment",
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. Relative Abundance (Genus Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_genus)











############################### 3-5_16S_GB_Relative_Abundance(Genus) ###############################

# ==============================================================================
# 3. Cohort-Specific Top 12 Identification & Master Color Dictionary
# ==============================================================================
# Top 12 Genera in GB cohort
top12_gb <- genus_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Genus != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Genus)

# Top 12 Genera in GJ cohort (for cross-cohort color consistency)
top12_gj <- genus_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Genus != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Genus)

# Union of all top genera across GB and GJ cohorts
all_top_genera <- unique(c(top12_gb, top12_gj))

# Build Master Named Palette (Shared genera retain identical hex codes)
palette_pool <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_top_genera))
master_taxa_colors <- c(setNames(palette_pool, all_top_genera), "Other" = "#B0B0B0")

# ==============================================================================
# 4. GB Genus-Level Aggregation (Two-Step Summary)
# ==============================================================================
# Filter dataset down to GB samples
genus_sample_gb <- genus_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment"))

# STEP 1: Collapse non-top 12 GB genera into "Other"
genus_sample_gb_grouped <- genus_sample_gb %>%
  mutate(Genus_Group = if_else(Genus %in% top12_gb, Genus, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Genus_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across GB Groups
genus_gb_summary <- genus_sample_gb_grouped %>%
  group_by(Treatment_Group, Genus_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Stacking order (Top GB taxa on top, "Other" at base)
gb_genus_levels <- c(top12_gb, "Other")
genus_gb_summary <- genus_gb_summary %>%
  mutate(Genus_Group = factor(Genus_Group, levels = rev(gb_genus_levels)))

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_genus_gb <- ggplot(genus_gb_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Genus_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. GB Relative Abundance (Genus Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_genus_gb)









############################### 3-5_16S_GJ_Relative_Abundance(Genus) ###############################

# ==============================================================================
# 3. Cohort-Specific Top 12 Identification & Master Color Dictionary
# ==============================================================================
# Top 12 Genera in GB cohort
top12_gb <- genus_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Genus != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Genus)

# Top 12 Genera in GJ cohort
top12_gj <- genus_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Genus != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 12) %>%
  pull(Genus)

# Union of all top genera across GB and GJ cohorts
all_top_genera <- unique(c(top12_gb, top12_gj))

# Build Master Named Palette (Shared genera retain identical hex codes)
palette_pool <- colorRampPalette(brewer.pal(12, "Set3"))(length(all_top_genera))
master_taxa_colors <- c(setNames(palette_pool, all_top_genera), "Other" = "#B0B0B0")

# ==============================================================================
# 4. GJ Genus-Level Aggregation (Two-Step Summary)
# ==============================================================================
# Filter dataset down to GJ samples
genus_sample_gj <- genus_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment"))

# STEP 1: Collapse non-top 12 GJ genera into "Other"
genus_sample_gj_grouped <- genus_sample_gj %>%
  mutate(Genus_Group = if_else(Genus %in% top12_gj, Genus, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Genus_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

# STEP 2: Calculate Mean Relative Abundance across GJ Groups
genus_gj_summary <- genus_sample_gj_grouped %>%
  group_by(Treatment_Group, Genus_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop")

# Stacking order (Top GJ taxa on top, "Other" at base)
gj_genus_levels <- c(top12_gj, "Other")
genus_gj_summary <- genus_gj_summary %>%
  mutate(Genus_Group = factor(Genus_Group, levels = rev(gj_genus_levels)))

# ==============================================================================
# 5. Visualization
# ==============================================================================
p_genus_gj <- ggplot(genus_gj_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Genus_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors) +
  scale_x_discrete(labels = c(
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(
    title = "3. GJ Relative Abundance (Genus Level)",
    y = "Mean Relative Abundance (%)"
  )

# Display Plot
print(p_genus_gj)









############################### 4-1_16S_GB_Volcano_Plot ###############################
library(qiime2R)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(MicrobiomeStat)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter for GB Cohort (GB_Control = M(B), GB_Treatment = B)
# ==============================================================================
meta_gb <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment"))
  )

# ==============================================================================
# 2. Import Feature Table & Taxonomy Artifacts
# ==============================================================================
counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw    <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

# Subset feature table to GB samples
common_samples <- intersect(meta_gb$`sample-id`, colnames(counts_raw))
counts_gb <- counts_raw[, common_samples]

# Calculate Mean Relative Abundance (%) per feature for dot sizing
rel_abund_mat <- apply(counts_gb, 2, function(x) (x / sum(x)) * 100)
mean_rel_abund <- rowMeans(rel_abund_mat) %>%
  enframe(name = "Feature.ID", value = "mean_rel_abund")

# ==============================================================================
# 3. Differential Abundance Analysis with LinDA
# ==============================================================================
meta_linda <- meta_gb %>%
  column_to_rownames("sample-id")

set.seed(42)
linda_res <- linda(
  feature.dat = counts_gb,
  meta.dat = meta_linda,
  formula = "~ Treatment_Group",
  alpha = 0.05,
  prev.filter = 0.10,
  mean.abund.filter = 0.0001
)

# Extract comparison results (GB_Treatment vs baseline GB_Control)
linda_df <- linda_res$output$Treatment_GroupGB_Treatment %>%
  rownames_to_column("Feature.ID") %>%
  mutate(
    lfc            = log2FoldChange,
    neg_log10_padj = -log10(padj)
  )

# Merge taxonomy, mean relative abundance, and define significance groups
volcano_df <- linda_df %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  inner_join(mean_rel_abund, by = "Feature.ID") %>%
  mutate(
    Significance = case_when(
      lfc >= 1.0 & padj < 0.05 ~ "Enriched in B",
      lfc <= -1.0 & padj < 0.05 ~ "Enriched in M(B)",
      TRUE ~ "Not Significant"
    ),
    Significance = factor(Significance, levels = c("Enriched in B", "Enriched in M(B)", "Not Significant")),
    Label_Name = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != "" ~ Genus,
      !is.na(Family) & Family != "f__" & Family != "" ~ paste("Fam.", Family),
      TRUE ~ paste("Phy.", Phylum)
    )
  )

# Select top significant features per group for text repelling
top_labels <- volcano_df %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_max(order_by = neg_log10_padj, n = 8, with_ties = FALSE) %>%
  ungroup()

# ==============================================================================
# 4. Volcano Plot Generation
# ==============================================================================
volcano_palette <- c(
  "Enriched in B"    = "#41B6C4",
  "Enriched in M(B)" = "#2B5C8F",
  "Not Significant"  = "#C0C0C0"
)

p_volcano <- ggplot(volcano_df, aes(x = lfc, y = neg_log10_padj, color = Significance)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(aes(size = mean_rel_abund), alpha = 0.75) +
  scale_size_continuous(range = c(1.5, 5.5), name = "Mean Rel. Abund (%)") +
  scale_color_manual(values = volcano_palette, name = "Significance") +
  scale_y_continuous(
    name = expression(bold("-Log"[10] * " (FDR-adjusted " * italic("p") * "-value)")),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  geom_text_repel(
    data = top_labels,
    aes(label = Label_Name, color = Significance),
    size = 3.5,
    fontface = "italic",
    family = font_family,
    max.overlaps = 15,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.4,
    show.legend = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = font_family),
    legend.text = element_text(size = 9, family = font_family),
    axis.text = element_text(color = "black", family = font_family, size = 11),
    axis.title = element_text(face = "bold", family = font_family, size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "4. B vs. M(B) Volcano Plot (LinDA)",
    x = expression(bold("Log"[2] * " Fold Change (B / M(B))"))
  )

# Display Plot
print(p_volcano)










############################### 4-1_16S_GB_Top20_Shift_Taxa_Boxplots ###############################
library(tidyverse)
library(ggpubr)

# Global font & palette configuration
font_family <- "Times New Roman"
palette_gb  <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4"
)

# ==============================================================================
# 1. Extract Top 20 Large-Shift Taxa Enriched in GB Treatment
# ==============================================================================
# Filter for significant features enriched in GB_Treatment ordered by highest Log2FC
top20_shift_taxa <- volcano_df %>%
  filter(Significance == "Enriched in B") %>%
  slice_max(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(desc(lfc)) %>%
  mutate(
    # Append short Feature ID to prevent facet merging if multiple ASVs share a Genus
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Name, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name) # Preserves LFC rank order
  )

# ==============================================================================
# 2. Subset Abundance Table & Format Metadata
# ==============================================================================
target_abund_df <- rel_abund_long %>%
  filter(
    Treatment_Group %in% c("GB_Control", "GB_Treatment"),
    Feature.ID %in% top20_shift_taxa$Feature.ID
  ) %>%
  inner_join(
    top20_shift_taxa %>% select(Feature.ID, Facet_Label, lfc),
    by = "Feature.ID"
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment"))
  )

# ==============================================================================
# 3. Generate Faceted Boxplot with Jitter
# ==============================================================================
p_top20_box <- ggplot(target_abund_df, aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(
    alpha = 0.65,
    outlier.shape = NA,
    width = 0.45,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = Treatment_Group),
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gb, name = "Treatment Group") +
  scale_color_manual(values = palette_gb, name = "Treatment Group") +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18)),
    limits = c(0, NA)
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = font_family, size = 11),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 8.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12, family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold.italic", size = 9, family = font_family),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(
    title = "Relative Abundance of Top 20 Large-Shift Taxa (GB Treatment)",
    y = "Relative Abundance (%)"
  )

# Display Plot
print(p_top20_box)












############################### 4-1_16S_GB_Top20_Depleted_Taxa_Boxplots ###############################
library(tidyverse)
library(ggpubr)

# Global font & palette configuration for GB cohort
font_family <- "Times New Roman"
palette_gb  <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4"
)

# ==============================================================================
# 1. Extract Top 20 Left-Side Taxa (Enriched in GB_Control / Depleted in GB_Treatment)
# ==============================================================================
top20_depleted_taxa_gb <- volcano_df %>%
  filter(Significance == "Enriched in M(B)") %>%
  slice_min(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(lfc) %>% # Sorts from most negative LFC upward
  mutate(
    # Clean up long genus names to avoid strip title truncation
    Label_Name_Clean = case_when(
      str_detect(Label_Name, "Burkholderia") ~ "Burkholderia group",
      nchar(Label_Name) > 22 ~ paste0(substr(Label_Name, 1, 20), ".."),
      TRUE ~ Label_Name
    ),
    # Append short Feature ID to guarantee unique facet keys
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Name_Clean, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name) # Preserves negative LFC rank order
  )

# ==============================================================================
# 2. Compute Long-Format Relative Abundance Matrix for GB Cohort
# ==============================================================================
rel_abund_gb <- counts_gb %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  inner_join(meta_gb, by = "sample-id")

target_abund_gb_depleted_df <- rel_abund_gb %>%
  filter(Feature.ID %in% top20_depleted_taxa_gb$Feature.ID) %>%
  inner_join(
    top20_depleted_taxa_gb %>% select(Feature.ID, Facet_Label, lfc),
    by = "Feature.ID"
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment"))
  )

# ==============================================================================
# 3. Generate Faceted Boxplot with Jitter
# ==============================================================================
p_top20_depleted_box_gb <- ggplot(target_abund_gb_depleted_df, aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(
    alpha = 0.65,
    outlier.shape = NA,
    width = 0.45,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = Treatment_Group),
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gb, name = "Treatment Group") +
  scale_color_manual(values = palette_gb, name = "Treatment Group") +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control",
    "GB_Treatment" = "GB_Treatment"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18)),
    limits = c(0, NA)
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = font_family, size = 11),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 8.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12, family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold.italic", size = 8.5, family = font_family),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(
    title = "Relative Abundance of Top 20 Depleted Taxa (Enriched in GB Control)",
    y = "Relative Abundance (%)"
  )

# Display Plot
print(p_top20_depleted_box_gb)












############################### 4-2_16S_GJ_Volcano_Plot ###############################
library(qiime2R)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(MicrobiomeStat)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter for GJ Cohort (GJ_Control = M(J), GJ_Treatment = J)
# ==============================================================================
meta_gj <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GJ_Control", "GJ_Treatment"))
  )

# ==============================================================================
# 2. Import Feature Table & Taxonomy Artifacts
# ==============================================================================
counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw    <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

# Subset feature table to GJ samples
common_samples <- intersect(meta_gj$`sample-id`, colnames(counts_raw))
counts_gj <- counts_raw[, common_samples]

# Calculate Mean Relative Abundance (%) per feature for dot sizing
rel_abund_mat <- apply(counts_gj, 2, function(x) (x / sum(x)) * 100)
mean_rel_abund <- rowMeans(rel_abund_mat) %>%
  enframe(name = "Feature.ID", value = "mean_rel_abund")

# ==============================================================================
# 3. Differential Abundance Analysis with LinDA
# ==============================================================================
meta_linda <- meta_gj %>%
  column_to_rownames("sample-id")

set.seed(42)
linda_res <- linda(
  feature.dat = counts_gj,
  meta.dat = meta_linda,
  formula = "~ Treatment_Group",
  alpha = 0.05,
  prev.filter = 0.10,
  mean.abund.filter = 0.0001
)

# Extract comparison results (GJ_Treatment vs baseline GJ_Control)
linda_df <- linda_res$output$Treatment_GroupGJ_Treatment %>%
  rownames_to_column("Feature.ID") %>%
  mutate(
    lfc            = log2FoldChange,
    neg_log10_padj = -log10(padj)
  )

# Merge taxonomy, mean relative abundance, and define significance groups
volcano_df <- linda_df %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  inner_join(mean_rel_abund, by = "Feature.ID") %>%
  mutate(
    Significance = case_when(
      lfc >= 1.0 & padj < 0.05 ~ "Enriched in J",
      lfc <= -1.0 & padj < 0.05 ~ "Enriched in M(J)",
      TRUE ~ "Not Significant"
    ),
    Significance = factor(Significance, levels = c("Enriched in J", "Enriched in M(J)", "Not Significant")),
    Label_Name = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != "" ~ Genus,
      !is.na(Family) & Family != "f__" & Family != "" ~ paste("Fam.", Family),
      TRUE ~ paste("Phy.", Phylum)
    )
  )

# Select top significant features per group for text repelling
top_labels <- volcano_df %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_max(order_by = neg_log10_padj, n = 8, with_ties = FALSE) %>%
  ungroup()

# ==============================================================================
# 4. Volcano Plot Generation
# ==============================================================================
volcano_palette <- c(
  "Enriched in J"    = "#FE9929",
  "Enriched in M(J)" = "#D95F02",
  "Not Significant"  = "#C0C0C0"
)

p_volcano <- ggplot(volcano_df, aes(x = lfc, y = neg_log10_padj, color = Significance)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(aes(size = mean_rel_abund), alpha = 0.75) +
  scale_size_continuous(range = c(1.5, 5.5), name = "Mean Rel. Abund (%)") +
  scale_color_manual(values = volcano_palette, name = "Significance") +
  scale_y_continuous(
    name = expression(bold("-Log"[10] * " (FDR-adjusted " * italic("p") * "-value)")),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  geom_text_repel(
    data = top_labels,
    aes(label = Label_Name, color = Significance),
    size = 3.5,
    fontface = "italic",
    family = font_family,
    max.overlaps = 15,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.4,
    show.legend = FALSE
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = font_family),
    legend.text = element_text(size = 9, family = font_family),
    axis.text = element_text(color = "black", family = font_family, size = 11),
    axis.title = element_text(face = "bold", family = font_family, size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "4. J vs. M(J) Volcano Plot (LinDA)",
    x = expression(bold("Log"[2] * " Fold Change (J / M(J))"))
  )

# Display Plot
print(p_volcano)










############################### 4-2_16S_GJ_Top20_Shift_Taxa_Boxplots ###############################
library(tidyverse)
library(ggpubr)

# Global font & palette configuration for GJ cohort
font_family <- "Times New Roman"
palette_gj  <- c(
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 1. Extract Top 20 Large-Shift Taxa Enriched in GJ Treatment
# ==============================================================================
top20_shift_taxa_gj <- volcano_df %>%
  filter(Significance == "Enriched in J") %>%
  slice_max(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(desc(lfc)) %>%
  mutate(
    # Clean up long genus names to avoid strip title truncation
    Label_Name_Clean = case_when(
      str_detect(Label_Name, "Burkholderia") ~ "Burkholderia group",
      nchar(Label_Name) > 22 ~ paste0(substr(Label_Name, 1, 20), ".."),
      TRUE ~ Label_Name
    ),
    # Append short Feature ID to guarantee unique facet keys
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Name_Clean, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name) # Preserves LFC rank order
  )

# ==============================================================================
# 2. Compute Long-Format Relative Abundance Matrix for GJ Cohort
# ==============================================================================
rel_abund_gj <- counts_gj %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  inner_join(meta_gj, by = "sample-id")

target_abund_gj_df <- rel_abund_gj %>%
  filter(Feature.ID %in% top20_shift_taxa_gj$Feature.ID) %>%
  inner_join(
    top20_shift_taxa_gj %>% select(Feature.ID, Facet_Label, lfc),
    by = "Feature.ID"
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GJ_Control", "GJ_Treatment"))
  )

# ==============================================================================
# 3. Generate Faceted Boxplot with Jitter
# ==============================================================================
p_top20_box_gj <- ggplot(target_abund_gj_df, aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(
    alpha = 0.65,
    outlier.shape = NA,
    width = 0.45,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = Treatment_Group),
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gj, name = "Treatment Group") +
  scale_color_manual(values = palette_gj, name = "Treatment Group") +
  scale_x_discrete(labels = c(
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18)),
    limits = c(0, NA)
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = font_family, size = 11),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 8.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12, family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold.italic", size = 8.5, family = font_family),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(
    title = "Relative Abundance of Top 20 Large-Shift Taxa (GJ Treatment)",
    y = "Relative Abundance (%)"
  )

# Display Plot
print(p_top20_box_gj)











############################### 4-2_16S_GJ_Top20_Depleted_Taxa_Boxplots ###############################
library(tidyverse)
library(ggpubr)

# Global font & palette configuration for GJ cohort
font_family <- "Times New Roman"
palette_gj  <- c(
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929"
)

# ==============================================================================
# 1. Extract Top 20 Left-Side Taxa (Enriched in GJ_Control / Depleted in GJ_Treatment)
# ==============================================================================
top20_depleted_taxa_gj <- volcano_df %>%
  filter(Significance == "Enriched in M(J)") %>%
  slice_min(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(lfc) %>% # Sorts from most negative LFC upward
  mutate(
    # Clean up long genus names to avoid strip title truncation
    Label_Name_Clean = case_when(
      str_detect(Label_Name, "Burkholderia") ~ "Burkholderia group",
      nchar(Label_Name) > 22 ~ paste0(substr(Label_Name, 1, 20), ".."),
      TRUE ~ Label_Name
    ),
    # Append short Feature ID to guarantee unique facet keys
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Name_Clean, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name) # Preserves negative LFC rank order
  )

# ==============================================================================
# 2. Compute Long-Format Relative Abundance Matrix for GJ Cohort
# ==============================================================================
rel_abund_gj <- counts_gj %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  inner_join(meta_gj, by = "sample-id")

target_abund_gj_depleted_df <- rel_abund_gj %>%
  filter(Feature.ID %in% top20_depleted_taxa_gj$Feature.ID) %>%
  inner_join(
    top20_depleted_taxa_gj %>% select(Feature.ID, Facet_Label, lfc),
    by = "Feature.ID"
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GJ_Control", "GJ_Treatment"))
  )

# ==============================================================================
# 3. Generate Faceted Boxplot with Jitter
# ==============================================================================
p_top20_depleted_box_gj <- ggplot(target_abund_gj_depleted_df, aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(
    alpha = 0.65,
    outlier.shape = NA,
    width = 0.45,
    color = "black",
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = Treatment_Group),
    width = 0.15,
    size = 1.8,
    alpha = 0.8
  ) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gj, name = "Treatment Group") +
  scale_color_manual(values = palette_gj, name = "Treatment Group") +
  scale_x_discrete(labels = c(
    "GJ_Control"   = "GJ_Control",
    "GJ_Treatment" = "GJ_Treatment"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.18)),
    limits = c(0, NA)
  ) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = font_family, size = 11),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", color = "black", size = 9),
    axis.text.y = element_text(color = "black", size = 8.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", size = 12, family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold.italic", size = 8.5, family = font_family),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(
    title = "Relative Abundance of Top 20 Depleted Taxa (Enriched in GJ Control)",
    y = "Relative Abundance (%)"
  )

# Display Plot
print(p_top20_depleted_box_gj)












############################### 5_16S_Heatmap ###############################
library(qiime2R)
library(tidyverse)
library(MicrobiomeStat)
library(pheatmap)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Artifacts (core_metrics directory)
# ==============================================================================
meta_raw <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(Domain == "Bacteria")

counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data

# ==============================================================================
# 2. Run LinDA Differential Abundance Analysis
# ==============================================================================
# --- Cohort 1: GB_Treatment vs. GB_Control ---
meta_gb <- meta_raw %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment"))) %>%
  column_to_rownames("sample-id")

common_samples_gb <- intersect(rownames(meta_gb), colnames(counts_raw))
counts_gb <- counts_raw[, common_samples_gb]
meta_gb   <- meta_gb[common_samples_gb, , drop = FALSE]

set.seed(42)
linda_gb <- linda(
  feature.dat       = counts_gb,
  meta.dat          = meta_gb,
  formula           = "~ Treatment_Group",
  alpha             = 0.05,
  prev.filter       = 0.10,
  mean.abund.filter = 0.0001
)

df_gb <- linda_gb$output$Treatment_GroupGB_Treatment %>%
  rownames_to_column("Feature.ID") %>%
  transmute(
    Feature.ID,
    Comparison = "GB_Treatment vs. GB_Control",
    lfc        = log2FoldChange,
    padj       = padj
  )

# --- Cohort 2: GJ_Treatment vs. GJ_Control ---
meta_gj <- meta_raw %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GJ_Control", "GJ_Treatment"))) %>%
  column_to_rownames("sample-id")

common_samples_gj <- intersect(rownames(meta_gj), colnames(counts_raw))
counts_gj <- counts_raw[, common_samples_gj]
meta_gj   <- meta_gj[common_samples_gj, , drop = FALSE]

set.seed(42)
linda_gj <- linda(
  feature.dat       = counts_gj,
  meta.dat          = meta_gj,
  formula           = "~ Treatment_Group",
  alpha             = 0.05,
  prev.filter       = 0.10,
  mean.abund.filter = 0.0001
)

df_gj <- linda_gj$output$Treatment_GroupGJ_Treatment %>%
  rownames_to_column("Feature.ID") %>%
  transmute(
    Feature.ID,
    Comparison = "GJ_Treatment vs. GJ_Control",
    lfc        = log2FoldChange,
    padj       = padj
  )

# ==============================================================================
# 3. Combine & Reshape Matrix for pheatmap
# ==============================================================================
combined_linda <- bind_rows(df_gb, df_gj)

# Identify features significant (padj < 0.05 & |lfc| >= 1.0) in AT LEAST ONE comparison
sig_features <- combined_linda %>%
  filter(padj < 0.05 & abs(lfc) >= 1.0) %>%
  pull(Feature.ID) %>%
  unique()

# Pivot wide into numeric matrix (Feature.ID x Comparison)
lfc_matrix <- combined_linda %>%
  filter(Feature.ID %in% sig_features) %>%
  select(Feature.ID, Comparison, lfc) %>%
  pivot_wider(names_from = Comparison, values_from = lfc, values_fill = 0) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

# Ensure explicit column ordering
lfc_matrix <- lfc_matrix[, c("GB_Treatment vs. GB_Control", "GJ_Treatment vs. GJ_Control")]

# Clamp extreme LFC values to [-8, 8] for balanced visual range
lfc_matrix[lfc_matrix > 8]  <- 8
lfc_matrix[lfc_matrix < -8] <- -8

# ==============================================================================
# 4. Generate Heatmap via pheatmap (Without Taxon Names)
# ==============================================================================
color_palette <- colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100)
breaks_seq    <- seq(-8, 8, length.out = 101)

p_heat <- pheatmap(
  mat                      = lfc_matrix,
  cluster_rows             = TRUE,
  cluster_cols             = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method        = "complete",
  show_rownames            = FALSE, # Hides individual taxon/ASV names
  show_colnames            = TRUE,
  color                    = color_palette,
  breaks                   = breaks_seq,
  treeheight_row           = 45,
  border_color             = NA,
  fontsize_col             = 11,
  angle_col                = "0",
  fontfamily               = font_family,
  main                     = "5. Heatmap (LinDA)"
)

# Display Plot
print(p_heat)








################################# Exclusive Rare Taxa Analysis #################################
library(qiime2R)
library(tidyverse)

# ==============================================================================
# 1. Load Final Corrected Metadata & Filter Groups
# ==============================================================================
meta_final <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Group_Code = case_when(
      Treatment_Group == "GB_Control" ~ "M(B)",
      Treatment_Group == "GB_Treatment" ~ "B",
      Treatment_Group == "GJ_Control" ~ "M(J)",
      Treatment_Group == "GJ_Treatment" ~ "J"
    )
  )

# ==============================================================================
# 2. Import Feature Table & Taxonomy Artifacts
# ==============================================================================
counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

# Subset feature table to valid samples (N=134)
common_samples <- intersect(meta_final$`sample-id`, colnames(counts_raw))
counts_sub <- counts_raw[, common_samples]

# Convert counts to relative abundance (%)
rel_abund_mat <- apply(counts_sub, 2, function(x) (x / sum(x)) * 100)

# ==============================================================================
# 3. Calculate Group-wise Prevalence & Relative Abundance Metrics
# ==============================================================================
df_long <- rel_abund_mat %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "RelAbund") %>%
  inner_join(meta_final %>% select(`sample-id`, Group_Code), by = "sample-id")

group_stats <- df_long %>%
  group_by(Feature.ID, Group_Code) %>%
  summarise(
    Sample_Count      = n(),
    Present_Count     = sum(RelAbund > 0),
    Prevalence_Pct    = (sum(RelAbund > 0) / n()) * 100,
    Mean_RelAbund_Pct = mean(RelAbund),
    .groups           = "drop"
  )

# Pivot prevalence and mean abundance for exclusive masking
prev_wide <- group_stats %>%
  select(Feature.ID, Group_Code, Prevalence_Pct) %>%
  pivot_wider(names_from = Group_Code, values_from = Prevalence_Pct, values_fill = 0)

abund_wide <- group_stats %>%
  select(Feature.ID, Group_Code, Mean_RelAbund_Pct) %>%
  pivot_wider(names_from = Group_Code, values_from = Mean_RelAbund_Pct, values_fill = 0, names_prefix = "MeanAbund_")

# ==============================================================================
# 4. Filter Exclusive Taxa (Uniquely Present in Exactly 1 Group)
# ==============================================================================
# Set minimum detection threshold to avoid single-read artifact false positives
min_samples_detected <- 2 # Must be present in >= 2 samples within the group

sample_counts_wide <- group_stats %>%
  select(Feature.ID, Group_Code, Present_Count) %>%
  pivot_wider(names_from = Group_Code, values_from = Present_Count, values_fill = 0)

exclusive_taxa_df <- prev_wide %>%
  inner_join(sample_counts_wide, by = "Feature.ID", suffix = c("_prev", "_count")) %>%
  mutate(
    Exclusive_Group = case_when(
      `M(B)_count` >= min_samples_detected & `B_count` == 0 & `M(J)_count` == 0 & `J_count` == 0 ~ "M(B)",
      `B_count` >= min_samples_detected & `M(B)_count` == 0 & `M(J)_count` == 0 & `J_count` == 0 ~ "B",
      `M(J)_count` >= min_samples_detected & `M(B)_count` == 0 & `B_count` == 0 & `J_count` == 0 ~ "M(J)",
      `J_count` >= min_samples_detected & `M(B)_count` == 0 & `B_count` == 0 & `M(J)_count` == 0 ~ "J",
      TRUE ~ "Shared / Below Threshold"
    )
  ) %>%
  filter(Exclusive_Group != "Shared / Below Threshold") %>%
  inner_join(abund_wide, by = "Feature.ID") %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  mutate(
    Detected_Sample_Count = case_when(
      Exclusive_Group == "M(B)" ~ `M(B)_count`,
      Exclusive_Group == "B" ~ `B_count`,
      Exclusive_Group == "M(J)" ~ `M(J)_count`,
      Exclusive_Group == "J" ~ `J_count`
    ),
    Prevalence_Pct = case_when(
      Exclusive_Group == "M(B)" ~ `M(B)_prev`,
      Exclusive_Group == "B" ~ `B_prev`,
      Exclusive_Group == "M(J)" ~ `M(J)_prev`,
      Exclusive_Group == "J" ~ `J_prev`
    ),
    Mean_RelAbund_Pct = case_when(
      Exclusive_Group == "M(B)" ~ `MeanAbund_M(B)`,
      Exclusive_Group == "B" ~ `MeanAbund_B`,
      Exclusive_Group == "M(J)" ~ `MeanAbund_M(J)`,
      Exclusive_Group == "J" ~ `MeanAbund_J`
    )
  ) %>%
  select(
    Feature.ID, Exclusive_Group, Detected_Sample_Count, Prevalence_Pct, Mean_RelAbund_Pct,
    Phylum, Class, Order, Family, Genus
  ) %>%
  arrange(Exclusive_Group, desc(Mean_RelAbund_Pct))

# ==============================================================================
# 5. Output Console Report & Export TSV
# ==============================================================================
cat("======================================================================\n")
cat("          EXCLUSIVE RARE TAXA COUNTS PER GROUP (Detected >= 2)        \n")
cat("======================================================================\n\n")
print(table(exclusive_taxa_df$Exclusive_Group))

write_tsv(exclusive_taxa_df, "Bacteria_16S_Exclusive_Rare_Taxa_Table.tsv")
cat("\nExclusive rare taxa table saved to 'Bacteria_16S_Exclusive_Rare_Taxa_Table.tsv'.\n")










############################### LinDA Differential Pathway Analysis (GB Group) ###############################
library(tidyverse)
library(MicrobiomeStat) # Provides linda()
library(ggpicrust2)     # For MetaCyc ID annotation
library(ggpubr)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter Samples (GB_Control vs GB_Treatment)
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment")
  ) %>%
  column_to_rownames("sample-id") %>%
  mutate(
    # GB_Control set as reference level (level 1)
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Control", "GB_Treatment")
    )
  )

# ==============================================================================
# 2. Load PICRUSt2 Pathway Abundance Matrix
# ==============================================================================
pathway_file <- "picrust2_out_16S/pathways_out/path_abun_unstrat.tsv.gz"

pathway_df <- read_delim(pathway_file, delim = "\t", comment = "#") %>%
  rename(pathway = 1)

path_matrix <- pathway_df %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Align sample IDs between metadata and abundance matrix
common_samples <- intersect(rownames(metadata), colnames(path_matrix))
path_matrix <- path_matrix[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# ==============================================================================
# 3. MetaCyc Pathway ID Mapping Function
# ==============================================================================
get_metacyc_annotations <- function(pathway_ids) {
  if (requireNamespace("ggpicrust2", quietly = TRUE)) {
    annot <- ggpicrust2::pathway_annotation(
      pathway = "MetaCyc",
      daa_results_df = data.frame(feature = pathway_ids),
      ko_to_kegg = FALSE
    )
    return(annot %>% select(feature, description))
  } else {
    map_url <- "https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/pathway_mapfiles/metacyc_pathways.tsv"
    map_df <- tryCatch(
      read_tsv(map_url, col_names = c("feature", "description"), show_col_types = FALSE),
      error = function(e) data.frame(feature = pathway_ids, description = pathway_ids)
    )
    return(map_df)
  }
}

# Annotate full pathway list
annot_lookup <- get_metacyc_annotations(rownames(path_matrix))

# ==============================================================================
# 4. Run LinDA Differential Abundance Analysis
# ==============================================================================
set.seed(42)
linda_res <- linda(
  feature.dat = path_matrix,
  meta.dat = metadata,
  formula = "~ Treatment_Group",
  alpha = 0.05,
  prev.filter = 0.1,  # Retain features present in at least 10% of samples
  zero.handling = "pseudo-count"
)

# Extract single contrast result table (GB_Treatment vs GB_Control)
linda_df <- linda_res$output[[1]] %>%
  rownames_to_column("feature") %>%
  rename(
    log2FC = log2FoldChange,
    p_value = pvalue,
    se = lfcSE
  ) %>%
  left_join(annot_lookup, by = "feature") %>%
  mutate(
    description = coalesce(description, feature),
    
    # Manual patch for missing/unmapped MetaCyc IDs
    description = case_when(
      feature == "PWY-8190" ~ "L-glutamate degradation XI (reductive Stickland reaction)",
      TRUE ~ description
    ),
    
    neg_log10_padj = -log10(padj),
    ci_lower = log2FC - 1.96 * se,
    ci_upper = log2FC + 1.96 * se,
    Enrichment = ifelse(log2FC > 0, "GB_Treatment", "GB_Control")
  )

# Filter significant pathways (FDR < 0.05)
sig_pathways <- linda_df %>%
  filter(padj < 0.05) %>%
  arrange(padj)

print(head(sig_pathways %>% select(feature, description, log2FC, se, padj), 10))

# Save full LinDA results table
write_csv(linda_df, "figure/PICRUSt2_LinDA_GB_Control_vs_GB_Treatment.csv")

# ==============================================================================
# 5. Forest Plot with FDR Size Mapping & 95% CIs
# ==============================================================================
# Take top 20 most significant pathways (or all significant if < 20)
top_sig <- sig_pathways %>%
  slice_min(padj, n = 20)

# Colors matching your project palette
gb_colors <- c(
  "GB_Treatment" = "#7BCCC4",
  "GB_Control"   = "#2B8CBE"
)

p_linda_forest <- ggplot(top_sig, aes(x = log2FC, y = reorder(description, log2FC))) +
  # Reference zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  
  # 95% Confidence Interval Error Bars
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = Enrichment), height = 0.25, linewidth = 0.7) +
  
  # Point size scaled by -log10(FDR)
  geom_point(aes(size = neg_log10_padj, fill = Enrichment), shape = 21, color = "black", stroke = 0.5) +
  
  # Color & Fill scale
  scale_color_manual(values = gb_colors) +
  scale_fill_manual(values = gb_colors) +
  
  # Continuous size scale for FDR visualization
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(2.5, 6.5)
  ) +
  
  # Publication Theme
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(color = "black", family = font_family),
    axis.title = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = font_family),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Differentially Abundant Pathways (GB_Treatment vs GB_Control)",
    x = expression(Log[2]~Fold~Change),
    y = "MetaCyc Pathway Description",
    color = "Enriched Group",
    fill = "Enriched Group"
  )

# Display Output
print(p_linda_forest)










############################### LinDA Differential Pathway Analysis (GJ Group) ###############################
library(tidyverse)
library(MicrobiomeStat) # Provides linda()
library(ggpicrust2)     # For MetaCyc ID annotation
library(ggpubr)

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter Samples (GJ_Control vs GJ_Treatment)
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GJ_Control", "GJ_Treatment")
  ) %>%
  column_to_rownames("sample-id") %>%
  mutate(
    # GJ_Control set as reference level (level 1)
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GJ_Control", "GJ_Treatment")
    )
  )

# ==============================================================================
# 2. Load PICRUSt2 Pathway Abundance Matrix
# ==============================================================================
pathway_file <- "picrust2_out_16S/pathways_out/path_abun_unstrat.tsv.gz"

pathway_df <- read_delim(pathway_file, delim = "\t", comment = "#") %>%
  rename(pathway = 1)

path_matrix <- pathway_df %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Align sample IDs between metadata and abundance matrix
common_samples <- intersect(rownames(metadata), colnames(path_matrix))
path_matrix <- path_matrix[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# ==============================================================================
# 3. MetaCyc Pathway ID Mapping Function
# ==============================================================================
get_metacyc_annotations <- function(pathway_ids) {
  if (requireNamespace("ggpicrust2", quietly = TRUE)) {
    annot <- ggpicrust2::pathway_annotation(
      pathway = "MetaCyc",
      daa_results_df = data.frame(feature = pathway_ids),
      ko_to_kegg = FALSE
    )
    return(annot %>% select(feature, description))
  } else {
    map_url <- "https://raw.githubusercontent.com/picrust/picrust2/master/picrust2/default_files/pathway_mapfiles/metacyc_pathways.tsv"
    map_df <- tryCatch(
      read_tsv(map_url, col_names = c("feature", "description"), show_col_types = FALSE),
      error = function(e) data.frame(feature = pathway_ids, description = pathway_ids)
    )
    return(map_df)
  }
}

# Annotate full pathway list
annot_lookup <- get_metacyc_annotations(rownames(path_matrix))

# ==============================================================================
# 4. Run LinDA Differential Abundance Analysis
# ==============================================================================
set.seed(42)
linda_res <- linda(
  feature.dat = path_matrix,
  meta.dat = metadata,
  formula = "~ Treatment_Group",
  alpha = 0.05,
  prev.filter = 0.1,  # Retain features present in at least 10% of samples
  zero.handling = "pseudo-count"
)

# Extract single contrast result table (GJ_Treatment vs GJ_Control)
linda_df <- linda_res$output[[1]] %>%
  rownames_to_column("feature") %>%
  rename(
    log2FC = log2FoldChange,
    p_value = pvalue,
    se = lfcSE
  ) %>%
  left_join(annot_lookup, by = "feature") %>%
  mutate(
    description = coalesce(description, feature),
    
    # Manual patch for missing/unmapped MetaCyc IDs
    description = case_when(
      feature == "PWY-8190" ~ "L-glutamate degradation XI (reductive Stickland reaction)",
      TRUE ~ description
    ),
    
    neg_log10_padj = -log10(padj),
    ci_lower = log2FC - 1.96 * se,
    ci_upper = log2FC + 1.96 * se,
    Enrichment = ifelse(log2FC > 0, "GJ_Treatment", "GJ_Control")
  )

# Filter significant pathways (FDR < 0.05)
sig_pathways <- linda_df %>%
  filter(padj < 0.05) %>%
  arrange(padj)

print(head(sig_pathways %>% select(feature, description, log2FC, se, padj), 10))

# Save full LinDA results table
write_csv(linda_df, "figure/PICRUSt2_LinDA_GJ_Control_vs_GJ_Treatment.csv")

# ==============================================================================
# 5. Forest Plot with FDR Size Mapping & 95% CIs
# ==============================================================================
# Take top 20 most significant pathways (or all significant if < 20)
top_sig <- sig_pathways %>%
  slice_min(padj, n = 20)

# Colors matching your project palette for GJ group
gj_colors <- c(
  "GJ_Treatment" = "#FDAE6B",
  "GJ_Control"   = "#E6550D"
)

p_linda_forest <- ggplot(top_sig, aes(x = log2FC, y = reorder(description, log2FC))) +
  # Reference zero line
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  
  # 95% Confidence Interval Error Bars
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = Enrichment), height = 0.25, linewidth = 0.7) +
  
  # Point size scaled by -log10(FDR)
  geom_point(aes(size = neg_log10_padj, fill = Enrichment), shape = 21, color = "black", stroke = 0.5) +
  
  # Color & Fill scale
  scale_color_manual(values = gj_colors) +
  scale_fill_manual(values = gj_colors) +
  
  # Continuous size scale for FDR visualization
  scale_size_continuous(
    name = expression(-log[10](FDR)),
    range = c(2.5, 6.5)
  ) +
  
  # Publication Theme
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(color = "black", family = font_family),
    axis.title = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = font_family),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Differentially Abundant Pathways (GJ_Treatment vs GJ_Control)",
    x = expression(Log[2]~Fold~Change),
    y = "MetaCyc Pathway Description",
    color = "Enriched Group",
    fill = "Enriched Group"
  )

# Display Output
print(p_linda_forest)










############################### PICRUSt2 Sample-Level Heatmap (Z-Score) ###############################
library(tidyverse)
library(pheatmap)
library(ggpicrust2) # For MetaCyc ID annotation

# Global font configuration
font_family <- "Times New Roman"

# ==============================================================================
# 1. Load Metadata & Filter Samples
# ==============================================================================
metadata <- read_tsv("metadata_P_generation.tsv") %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(
      Treatment_Group,
      levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
    )
  ) %>%
  arrange(Treatment_Group) %>% # Order samples sequentially by group
  column_to_rownames("sample-id")

# ==============================================================================
# 2. Load PICRUSt2 Pathway Abundance Matrix
# ==============================================================================
pathway_file <- "picrust2_out_16S/pathways_out/path_abun_unstrat.tsv.gz"

pathway_df <- read_delim(pathway_file, delim = "\t", comment = "#") %>%
  rename(pathway = 1)

path_matrix <- pathway_df %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Align samples between abundance matrix and ordered metadata
common_samples <- intersect(rownames(metadata), colnames(path_matrix))
path_matrix <- path_matrix[, common_samples]
metadata <- metadata[common_samples, , drop = FALSE]

# Total Sum Scaling (TSS) relative abundance normalization
path_rel <- sweep(path_matrix, 2, colSums(path_matrix), FUN = "/")

# ==============================================================================
# 3. Select Top Variable Pathways for Heatmap Visualization
# ==============================================================================
# Select top 35 most variable pathways across samples
top_n <- 35
pathway_variances <- apply(path_rel, 1, sd)
top_pathways <- names(sort(pathway_variances, decreasing = TRUE))[1:top_n]

path_sub <- path_rel[top_pathways, ]

# ==============================================================================
# 4. MetaCyc Pathway ID Mapping & Label Cleaning
# ==============================================================================
annot_df <- get_metacyc_annotations(rownames(path_sub))
mapped_names <- setNames(annot_df$description, annot_df$feature)

row_labels <- ifelse(!is.na(mapped_names[rownames(path_sub)]), 
                     mapped_names[rownames(path_sub)], 
                     rownames(path_sub))

# Clean HTML formatting entities
row_labels <- row_labels %>%
  stringr::str_replace_all("&beta;", "beta") %>%
  stringr::str_replace_all("&alpha;", "alpha")

# Apply manual patches for missing IDs
row_labels <- dplyr::case_when(
  row_labels == "PWY-7883" ~ "4-hydroxy-2-nonenal degradation",
  row_labels == "PWY-8004" ~ "ent-kaurene biosynthesis I",
  TRUE ~ row_labels
) %>%
  stringr::str_trunc(width = 75, side = "right")

rownames(path_sub) <- make.unique(row_labels)

# ==============================================================================
# 5. Configure Heatmap Annotations & Colors
# ==============================================================================
# Column annotation frame for pheatmap
annotation_col <- metadata %>%
  select(Treatment_Group)

# Custom color scheme for treatment groups
annotation_colors <- list(
  Treatment_Group = c(
    "GB_Control"   = "#2B8CBE",
    "GB_Treatment" = "#7BCCC4",
    "GJ_Control"   = "#E6550D",
    "GJ_Treatment" = "#FDAE6B"
  )
)

# Heatmap gradient (Blue -> White -> Red for Z-scores)
heatmap_colors <- colorRampPalette(c("#1F618D", "#FFFFFF", "#B03A2E"))(100)

# Calculate group boundaries for column separators (gaps)
group_counts <- table(metadata$Treatment_Group)
gaps_col <- cumsum(group_counts)[1:(length(group_counts) - 1)]

# ==============================================================================
# 6. Generate & Print Heatmap
# ==============================================================================
p_heatmap <- pheatmap(
  mat = path_sub,
  scale = "row",
  color = heatmap_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  gaps_col = gaps_col,
  show_colnames = FALSE,            # Hide overlapping sample IDs on X-axis
  show_rownames = TRUE,
  fontsize = 9,
  fontsize_row = 8,
  main = "Sample-level Predicted Pathway Profiles (Z-score)"
)

# Render to screen
grid::grid.newpage()
grid::grid.draw(p_heatmap$gtable)
