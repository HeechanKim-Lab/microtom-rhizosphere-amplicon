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
  "ggpicrust2",
  "randomForest"
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
counts_raw <- read_qza("05_taxonomy_16S/table_16S_clean.qza")$data  # unrarefied counts; LinDA corrects for read depth itself
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
  mean.abund.filter = 0.0001,
  adaptive          = FALSE,  # zeros imputed; set explicitly (log 2026-08-24_01, Section 1)
  zero.handling     = "imputation"
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
write_csv(volcano_df, "figure/LinDA_ASV_GB_Treatment_vs_GB_Control.csv")










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
# Figure 080: relative abundance from counts_gb (unrarefied), the same way the Depleted section builds rel_abund_gb
rel_abund_long <- counts_gb %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  inner_join(meta_gb, by = "sample-id")

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
counts_raw <- read_qza("05_taxonomy_16S/table_16S_clean.qza")$data  # unrarefied counts; LinDA corrects for read depth itself
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
  mean.abund.filter = 0.0001,
  adaptive          = FALSE,  # zeros imputed; set explicitly (log 2026-08-24_01, Section 1)
  zero.handling     = "imputation"
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
write_csv(volcano_df, "figure/LinDA_ASV_GJ_Treatment_vs_GJ_Control.csv")










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

counts_raw <- read_qza("05_taxonomy_16S/table_16S_clean.qza")$data  # unrarefied counts; LinDA corrects for read depth itself

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
  mean.abund.filter = 0.0001,
  adaptive          = FALSE,  # zeros imputed; set explicitly (log 2026-08-24_01, Section 1)
  zero.handling     = "imputation"
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
  mean.abund.filter = 0.0001,
  adaptive          = FALSE,  # zeros imputed; set explicitly (log 2026-08-24_01, Section 1)
  zero.handling     = "imputation"
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










############################### 9-1_16S_RandomForest_Treatment ###############################
# Ensure randomForest is installed
if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}

library(qiime2R)
library(tidyverse)
library(randomForest)
library(pheatmap)
library(RColorBrewer)

# Global font & aesthetic configuration
font_family <- "Times New Roman"
palette_rf <- c(
  "GB_Treatment" = "#1A365D",
  "GJ_Treatment" = "#881337"
)

# ==============================================================================
# 1. Load Metadata & Rarefied Count Table (Isolate Model 1: GB vs. GJ Treatments)
# ==============================================================================
# Filter strictly to the 55 inoculated plant samples
meta_rf <- read_tsv("metadata_P_generation.tsv", show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Treatment", "GJ_Treatment")
  ) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Treatment", "GJ_Treatment"))) %>%
  column_to_rownames("sample-id")

# Import ASV table and taxonomy
counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw    <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

# Calculate global relative abundance (%)
rel_abund_mat <- sweep(counts_raw, 2, colSums(counts_raw), FUN = "/") * 100

# Subset abundance matrix to the 55 target samples
common_samples <- intersect(rownames(meta_rf), colnames(rel_abund_mat))
meta_rf <- meta_rf[common_samples, , drop = FALSE]
rel_rf  <- rel_abund_mat[, common_samples]

# ==============================================================================
# 2. Unsupervised Pre-Filtering (Noise Reduction)
# ==============================================================================
# Retain ASVs present in >= 15% of samples (>= 9 samples) AND mean relative abundance >= 0.005%
# This is unsupervised (label-agnostic) and preserves FDR / model integrity
prev_threshold <- 0.15 * ncol(rel_rf)
keep_features <- rownames(rel_rf)[
  rowSums(rel_rf > 0) >= prev_threshold & rowMeans(rel_rf) >= 0.005
]

cat(sprintf("\nFiltering Complete: Retained %d ASVs out of %d for Random Forest classification.\n", 
            length(keep_features), nrow(rel_rf)))

# Transpose for modeling (Rows = Samples, Columns = ASVs)
rf_input <- as.data.frame(t(rel_rf[keep_features, ]))
# Explicitly separate X (features) and Y (target) to avoid R formula syntax breaking on ASV hash names
x_data <- rf_input
y_data <- meta_rf$Treatment_Group

# ==============================================================================
# 3. Model Training (Direct Parameterization)
# ==============================================================================
set.seed(42)

# Set mtry to standard theoretical optimum: floor(sqrt(p))
best_mtry <- floor(sqrt(ncol(x_data))) # 47 for 2,261 ASVs

cat(sprintf("Proceeding with mtry = %d across 2,000 trees.\n", best_mtry))

# Train the production Random Forest model
rf_model <- randomForest(
  x = x_data,
  y = y_data,
  ntree = 2000,               # Ensures Mean Decrease Accuracy (MDA) scores stabilize
  mtry = best_mtry,           # Evaluates 47 features per split
  importance = TRUE,          # Enables MDA calculation via OOB permutation
  proximity = TRUE,
  sampsize = c(22, 22),       # Stratified 1:1 sampling per bootstrap to balance classes
  nodesize = 2                # Structural regularization against singleton noise
)

cat("\n=== Random Forest OOB Classification Report ===\n")
print(rf_model)

# ==============================================================================
# 4. Extract Permutation Importance (MDA) & Assign Taxonomy
# ==============================================================================
imp_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column("Feature.ID") %>%
  rename(
    MeanDecreaseAccuracy = MeanDecreaseAccuracy,
    MeanDecreaseGini     = MeanDecreaseGini
  ) %>%
  left_join(tax_clean, by = "Feature.ID") %>%
  mutate(
    Clean_Name = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != "" ~ Genus,
      !is.na(Family) & Family != "f__" & Family != "" ~ paste("Fam.", Family),
      !is.na(Order) & Order != "o__" & Order != "" ~ paste("Ord.", Order),
      TRUE ~ paste("Phy.", Phylum)
    ),
    Short_ID = substr(Feature.ID, 1, 5),
    DisplayName = paste0(Clean_Name, " (", Short_ID, ")")
  ) %>%
  arrange(desc(MeanDecreaseAccuracy))

# Extract Top 30 Predictive Biomarkers
top30_rf <- imp_df %>%
  slice_head(n = 30)

# Determine enrichment direction (GB_Treatment vs. GJ_Treatment)
group_means <- rf_input %>%
  mutate(Target_Class = y_data) %>%
  group_by(Target_Class) %>%
  summarise(across(all_of(top30_rf$Feature.ID), mean)) %>%
  pivot_longer(-Target_Class, names_to = "Feature.ID", values_to = "Mean_Abund")

enrich_dir <- group_means %>%
  group_by(Feature.ID) %>%
  slice_max(Mean_Abund, n = 1) %>%
  select(Feature.ID, Enriched_In = Target_Class)

# Merge enrichment direction and order factor levels for plotting
top30_rf <- top30_rf %>%
  left_join(enrich_dir, by = "Feature.ID") %>%
  mutate(DisplayName = factor(DisplayName, levels = rev(DisplayName)))

# Export Diagnostic Biomarkers
dir.create("figure", showWarnings = FALSE)
write_csv(top30_rf, "figure/RF_16S_Top30_Diagnostic_Biomarkers.csv")
cat("\nTop 30 diagnostic biomarkers exported to 'figure/RF_16S_Top30_Diagnostic_Biomarkers.csv'.\n")

# ==============================================================================
# 5. Visualization: Variable Importance Lollipop Chart
# ==============================================================================
p_importance <- ggplot(top30_rf, aes(x = MeanDecreaseAccuracy, y = DisplayName)) +
  geom_segment(
    aes(x = 0, xend = MeanDecreaseAccuracy, y = DisplayName, yend = DisplayName), 
    color = "gray60", 
    linewidth = 0.6
  ) +
  geom_point(
    aes(size = MeanDecreaseGini, fill = Enriched_In), 
    shape = 21, 
    stroke = 0.5, 
    color = "black"
  ) +
  scale_fill_manual(values = palette_rf, name = "Enriched Group") +
  scale_size_continuous(range = c(2.5, 6.0), name = "Mean Decrease Gini") +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    axis.text.y = element_text(face = "italic", size = 9.5, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "6. Top 30 Diagnostic Biomarkers (GB vs. GJ Rhizosphere)",
    x = "Variable Importance: Mean Decrease Accuracy (OOB Permutation)",
    y = "Taxonomic Assignment (ASV ID)"
  )

print(p_importance)

# ==============================================================================
# 6. Sample-Level Abundance Profiling (Z-Score Heatmap)
# ==============================================================================
top20_ids <- top30_rf$Feature.ID[1:20]
rf_heatmap_mat <- rel_rf[top20_ids, ]
rownames(rf_heatmap_mat) <- top30_rf$DisplayName[match(rownames(rf_heatmap_mat), top30_rf$Feature.ID)]

# Order columns by treatment group
col_annot <- meta_rf %>%
  select(Treatment_Group) %>%
  arrange(Treatment_Group)

rf_heatmap_mat <- rf_heatmap_mat[, rownames(col_annot)]

heatmap_colors <- colorRampPalette(c("#1F618D", "#FFFFFF", "#B03A2E"))(100)

p_rf_heat <- pheatmap(
  mat = rf_heatmap_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_annot,
  annotation_colors = list(Treatment_Group = palette_rf),
  color = heatmap_colors,
  show_colnames = FALSE,
  fontsize_row = 10,
  fontfamily = font_family,
  border_color = NA,
  main = "6. Relative Abundance of Top Diagnostic RF Biomarkers (Z-Score)"
)

grid::grid.draw(p_rf_heat$gtable)










############################### 9-3_16S_RandomForest_4Groups ###############################
library(qiime2R)
library(tidyverse)
library(randomForest)
library(pheatmap)
library(RColorBrewer)

# Global font & aesthetic configuration matching previous logs
font_family <- "Times New Roman"
palette_4 <- c(
  "GB_Control"   = "#4A90E2",
  "GB_Treatment" = "#1A365D",
  "GJ_Control"   = "#F87171",
  "GJ_Treatment" = "#881337"
)

# ==============================================================================
# 1. Load Metadata & Rarefied Count Table (Complete N = 110 Cohort)
# ==============================================================================
meta_4g <- read_tsv("metadata_P_generation.tsv", show_col_types = FALSE) %>%
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
  column_to_rownames("sample-id")

# Import ASV table and taxonomy
counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw    <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

# Calculate relative abundance (%)
rel_abund_mat <- sweep(counts_raw, 2, colSums(counts_raw), FUN = "/") * 100

# Align samples to locked N = 110 cohort
common_samples <- intersect(rownames(meta_4g), colnames(rel_abund_mat))
meta_4g <- meta_4g[common_samples, , drop = FALSE]
rel_4g  <- rel_abund_mat[, common_samples]

# ==============================================================================
# 2. Unsupervised Pre-Filtering (N = 110)
# ==============================================================================
# Retain ASVs present in >= 10% of samples (>= 11 samples) with mean relative abundance >= 0.005%
prev_threshold <- 0.10 * ncol(rel_4g)
keep_features <- rownames(rel_4g)[
  rowSums(rel_4g > 0) >= prev_threshold & rowMeans(rel_4g) >= 0.005
]

cat(sprintf("\nFiltering Complete: Retained %d ASVs out of %d across 110 samples.\n", 
            length(keep_features), nrow(rel_4g)))

x_data <- as.data.frame(t(rel_4g[keep_features, ]))
y_data <- meta_4g$Treatment_Group

# ==============================================================================
# 3. Multiclass Random Forest Training
# ==============================================================================
set.seed(42)

# Theoretical classification mtry: floor(sqrt(p))
best_mtry <- floor(sqrt(ncol(x_data)))

cat(sprintf("Training Multiclass Random Forest (mtry = %d, ntree = 2000, 4 balanced classes).\n", best_mtry))

rf_4g_model <- randomForest(
  x = x_data,
  y = y_data,
  ntree = 2000,
  mtry = best_mtry,
  importance = TRUE,
  proximity = TRUE,
  sampsize = c(22, 22, 22, 22), # Balances 33:33:22:22 representation across trees
  nodesize = 2
)

cat("\n=== Multiclass Random Forest Confusion Matrix & OOB Report ===\n")
print(rf_4g_model)

# ==============================================================================
# 4. Extract Multiclass Permutation Importance & Assign Taxonomy
# ==============================================================================
imp_matrix <- importance(rf_4g_model)

# Extract global MeanDecreaseAccuracy alongside per-class accuracy drops
imp_df <- as.data.frame(imp_matrix) %>%
  rownames_to_column("Feature.ID") %>%
  left_join(tax_clean, by = "Feature.ID") %>%
  mutate(
    Clean_Name = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != "" ~ Genus,
      !is.na(Family) & Family != "f__" & Family != "" ~ paste("Fam.", Family),
      !is.na(Order) & Order != "o__" & Order != "" ~ paste("Ord.", Order),
      TRUE ~ paste("Phy.", Phylum)
    ),
    Short_ID = substr(Feature.ID, 1, 5),
    DisplayName = paste0(Clean_Name, " (", Short_ID, ")")
  ) %>%
  arrange(desc(MeanDecreaseAccuracy))

# Extract Top 30 Predictive Multiclass Biomarkers
top30_rf_4g <- imp_df %>%
  slice_head(n = 30)

# Determine Group of Maximum Mean Abundance (Enrichment Core)
group_means <- x_data %>%
  mutate(Target_Class = y_data) %>%
  group_by(Target_Class) %>%
  summarise(across(all_of(top30_rf_4g$Feature.ID), mean)) %>%
  pivot_longer(-Target_Class, names_to = "Feature.ID", values_to = "Mean_Abund")

enrich_dir <- group_means %>%
  group_by(Feature.ID) %>%
  slice_max(Mean_Abund, n = 1) %>%
  select(Feature.ID, Enriched_In = Target_Class)

top30_rf_4g <- top30_rf_4g %>%
  left_join(enrich_dir, by = "Feature.ID") %>%
  mutate(DisplayName = factor(DisplayName, levels = rev(DisplayName)))

dir.create("figure", showWarnings = FALSE)
write_csv(top30_rf_4g, "figure/RF_16S_4Group_Top30_Biomarkers.csv")
cat("\nMulticlass biomarkers exported to 'figure/RF_16S_4Group_Top30_Biomarkers.csv'.\n")

# ==============================================================================
# 5. Visualization: Multiclass Variable Importance Lollipop Chart
# ==============================================================================
p_importance_4g <- ggplot(top30_rf_4g, aes(x = MeanDecreaseAccuracy, y = DisplayName)) +
  geom_segment(
    aes(x = 0, xend = MeanDecreaseAccuracy, y = DisplayName, yend = DisplayName), 
    color = "gray60", 
    linewidth = 0.6
  ) +
  geom_point(
    aes(size = MeanDecreaseGini, fill = Enriched_In), 
    shape = 21, 
    stroke = 0.5, 
    color = "black"
  ) +
  scale_fill_manual(values = palette_4, name = "Primary Enriched Group") +
  scale_size_continuous(range = c(2.5, 6.0), name = "Mean Decrease Gini") +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    axis.text.y = element_text(face = "italic", size = 9.5, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "7. Top 30 Multiclass Biomarkers (4-Group Model)",
    x = "Variable Importance: Mean Decrease Accuracy (OOB Permutation)",
    y = "Taxonomic Assignment (ASV ID)"
  )

print(p_importance_4g)

# ==============================================================================
# 6. Sample-Level Abundance Profiling (N = 110 Heatmap)
# ==============================================================================
top20_ids <- top30_rf_4g$Feature.ID[1:20]
rf_heat_mat <- rel_4g[top20_ids, ]
rownames(rf_heat_mat) <- top30_rf_4g$DisplayName[match(rownames(rf_heat_mat), top30_rf_4g$Feature.ID)]

# Order columns sequentially across all 4 groups
col_annot <- meta_4g %>%
  select(Treatment_Group) %>%
  arrange(Treatment_Group)

rf_heat_mat <- rf_heat_mat[, rownames(col_annot)]

# Group boundary gaps for pheatmap visual clarity
group_counts <- table(col_annot$Treatment_Group)
gaps_col <- cumsum(group_counts)[1:(length(group_counts) - 1)]

heatmap_colors <- colorRampPalette(c("#1F618D", "#FFFFFF", "#B03A2E"))(100)

p_rf_4g_heat <- pheatmap(
  mat = rf_heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_annot,
  annotation_colors = list(Treatment_Group = palette_4),
  gaps_col = gaps_col,
  color = heatmap_colors,
  show_colnames = FALSE,
  fontsize_row = 9.5,
  fontfamily = font_family,
  border_color = NA,
  main = "7. Relative Abundance of Top 20 Multiclass RF Biomarkers (N = 110, Z-Score)"
)

grid::grid.draw(p_rf_4g_heat$gtable)










############################### 10-1_16S_RF_Marker_Selection ###############################
# Input:  figure/RF_16S_Top30_Diagnostic_Biomarkers.csv                         (section 9-1)
#         figure/LinDA_ASV_GB_Treatment_vs_GB_Control.csv                        (section 4-1)
#         figure/LinDA_ASV_GJ_Treatment_vs_GJ_Control.csv                        (section 4-2)
#         picrust2_out_16S/EC_metagenome_out/seqtab_norm.tsv.gz
#         picrust2_out_16S/combined_EC_predicted.tsv.gz
# Output: rf_picrust2/rf_markers.tsv
#         rf_picrust2/seqtab_norm_markers.tsv.gz, rf_picrust2/EC_predicted_markers.tsv.gz (for the Terminal step below)
library(tidyverse)

out_dir <- "rf_picrust2"
dir.create(out_dir, showWarnings = FALSE)

# ==============================================================================
# 1. RF top 30 (section 9-1), checked against LinDA (treatment vs its own control)
# ==============================================================================
rf_top30 <- read_csv("figure/RF_16S_Top30_Diagnostic_Biomarkers.csv", show_col_types = FALSE)

read_linda_calls <- function(file, soil) {
  read_csv(file, show_col_types = FALSE) %>%
    transmute(
      Feature.ID,
      Soil = soil,
      Call = case_when(
        Significance %in% c("Enriched in B", "Enriched in J")       ~ "up",
        Significance %in% c("Enriched in M(B)", "Enriched in M(J)") ~ "down",
        TRUE                                                        ~ "ns"
      )
    )
}

linda_calls <- bind_rows(
  read_linda_calls("figure/LinDA_ASV_GB_Treatment_vs_GB_Control.csv", "LinDA_GB"),
  read_linda_calls("figure/LinDA_ASV_GJ_Treatment_vs_GJ_Control.csv", "LinDA_GJ")
) %>%
  pivot_wider(names_from = Soil, values_from = Call)

# GB-specific: RF marker of GB_Treatment, enriched vs GB_Control, and not enriched in GJ
# GJ-specific: RF marker of GJ_Treatment, enriched vs GJ_Control, and not enriched in GB
# "not tested" = removed by LinDA's 10% prevalence filter (absent or rare in that soil)
markers <- rf_top30 %>%
  select(Feature.ID, DisplayName, Enriched_In, MeanDecreaseAccuracy, Phylum, Family, Genus) %>%
  left_join(linda_calls, by = "Feature.ID") %>%
  mutate(
    across(c(LinDA_GB, LinDA_GJ), ~ replace_na(.x, "not tested")),
    Marker_Set = case_when(
      Enriched_In == "GB_Treatment" & LinDA_GB == "up" & LinDA_GJ != "up" ~ "GB-specific",
      Enriched_In == "GJ_Treatment" & LinDA_GJ == "up" & LinDA_GB != "up" ~ "GJ-specific",
      TRUE                                                                ~ "excluded"
    )
  ) %>%
  arrange(Marker_Set, desc(MeanDecreaseAccuracy))

print(count(markers, Enriched_In, Marker_Set))
print(markers %>% filter(Marker_Set == "excluded") %>% select(DisplayName, Enriched_In, LinDA_GB, LinDA_GJ))

write_tsv(markers, file.path(out_dir, "rf_markers.tsv"))

# ==============================================================================
# 2. Marker-only copies of the two per-ASV PICRUSt2 tables (rows copied unchanged)
# ==============================================================================
marker_ids <- markers %>% filter(Marker_Set != "excluded") %>% pull(Feature.ID)

subset_rows <- function(infile, outfile, ids) {
  lines <- readLines(gzfile(infile))
  row_ids <- sub("\t.*$", "", lines[-1])
  keep <- c(lines[1], lines[-1][row_ids %in% ids])
  con <- gzfile(outfile, "w")
  writeLines(keep, con)
  close(con)
  missing <- setdiff(ids, row_ids)
  if (length(missing) > 0) warning(infile, " has no row for: ", paste(missing, collapse = ", "))
  cat(sprintf("%s: %d of %d marker ASVs\n", outfile, length(keep) - 1, length(ids)))
}

subset_rows("picrust2_out_16S/EC_metagenome_out/seqtab_norm.tsv.gz",
            file.path(out_dir, "seqtab_norm_markers.tsv.gz"), marker_ids)
subset_rows("picrust2_out_16S/combined_EC_predicted.tsv.gz",
            file.path(out_dir, "EC_predicted_markers.tsv.gz"), marker_ids)










############################### PICRUSt2 Pathways of RF Marker ASVs (Terminal) ###############################
# Run in Terminal (not in R) after 10-1 and before 10-2, from the same working directory.
# --per_sequence_contrib runs MinPath on each marker ASV's predicted genome separately.
#
#   conda activate picrust2
#
#   pathway_pipeline.py \
#     -i picrust2_out_16S/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
#     -o rf_picrust2/pathways_markers \
#     --per_sequence_contrib \
#     --per_sequence_abun rf_picrust2/seqtab_norm_markers.tsv.gz \
#     --per_sequence_function rf_picrust2/EC_predicted_markers.tsv.gz \
#     --coverage \
#     -p 4
#
#   add_descriptions.py \
#     -i rf_picrust2/pathways_markers/path_abun_unstrat_per_seq.tsv.gz \
#     -m METACYC \
#     -o rf_picrust2/pathways_markers/path_abun_unstrat_per_seq_descrip.tsv.gz
#
# The same commands are in modules/15_marker_pathways_picrust2.sh.










############################### 10-2_16S_RF_Marker_Pathways ###############################
# Input:  rf_picrust2/rf_markers.tsv                                             (section 10-1)
#         rf_picrust2/pathways_markers/path_abun_predictions.tsv.gz               (Terminal step)
#         rf_picrust2/pathways_markers/path_abun_unstrat_per_seq_descrip.tsv.gz   (Terminal step)
#         figure/PICRUSt2_LinDA_GB_Control_vs_GB_Treatment.csv                    (pathway LinDA, GB)
#         figure/PICRUSt2_LinDA_GJ_Control_vs_GJ_Treatment.csv                    (pathway LinDA, GJ)
# Output: rf_picrust2/pathway_comparison_GB_vs_GJ.tsv, marker_pathway_presence.tsv, marker_pathway_counts.tsv
#         rf_picrust2/RF_markers_pathway_share.png, rf_picrust2/RF_markers_pathway_heatmap.png
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(RColorBrewer)

font_family <- "Times New Roman"
out_dir <- "rf_picrust2"
pwy_dir <- file.path(out_dir, "pathways_markers")

palette_set <- c("GB-specific" = "#1A365D", "GJ-specific" = "#881337")

# ==============================================================================
# 1. Load marker list, per-ASV pathway predictions and pathway descriptions
# ==============================================================================
markers <- read_tsv(file.path(out_dir, "rf_markers.tsv"), show_col_types = FALSE) %>%
  filter(Marker_Set != "excluded")
n_set <- count(markers, Marker_Set) %>% deframe()

# Predicted abundance of each pathway in each marker ASV's genome (MinPath run per genome).
# A value > 0 means the pathway is predicted to be present in that ASV.
pwy_by_asv <- read_tsv(file.path(pwy_dir, "path_abun_predictions.tsv.gz"), show_col_types = FALSE) %>%
  pivot_longer(-sequence, names_to = "pathway", values_to = "genome_abun") %>%
  rename(Feature.ID = sequence) %>%
  inner_join(markers %>% select(Feature.ID, DisplayName, Marker_Set, Phylum), by = "Feature.ID") %>%
  mutate(Present = genome_abun > 0)

# MetaCyc descriptions added by add_descriptions.py in the Terminal step (HTML entities such as &beta; -> beta)
pwy_desc <- read_tsv(file.path(pwy_dir, "path_abun_unstrat_per_seq_descrip.tsv.gz"), show_col_types = FALSE) %>%
  transmute(pathway, description = str_replace_all(description, "&([A-Za-z]+);", "\\1"))

# Community-level result for each pathway (pathway LinDA, treated vs its own control)
pwy_linda <- bind_rows(
  read_csv("figure/PICRUSt2_LinDA_GB_Control_vs_GB_Treatment.csv", show_col_types = FALSE) %>% mutate(Soil = "GB"),
  read_csv("figure/PICRUSt2_LinDA_GJ_Control_vs_GJ_Treatment.csv", show_col_types = FALSE) %>% mutate(Soil = "GJ")
) %>%
  transmute(
    pathway = feature,
    Soil,
    Community = case_when(
      padj < 0.05 & log2FC > 0 ~ "higher in treated",
      padj < 0.05 & log2FC < 0 ~ "lower in treated",
      TRUE                     ~ "no change"
    )
  ) %>%
  pivot_wider(names_from = Soil, values_from = Community, names_prefix = "Community_")

# ==============================================================================
# 2. Compare GB-specific vs GJ-specific markers, pathway by pathway
# ==============================================================================
asv_summary <- pwy_by_asv %>%
  group_by(Marker_Set, DisplayName, Phylum) %>%
  summarise(n_pathways = sum(Present), .groups = "drop")

cat("\nPredicted pathways per marker ASV:\n")
print(asv_summary %>% group_by(Marker_Set) %>%
        summarise(n_ASVs = n(), median = median(n_pathways), min = min(n_pathways), max = max(n_pathways)))

comparison <- pwy_by_asv %>%
  group_by(pathway) %>%
  summarise(
    n_GB = sum(Present & Marker_Set == "GB-specific"),
    n_GJ = sum(Present & Marker_Set == "GJ-specific"),
    .groups = "drop"
  ) %>%
  filter(n_GB + n_GJ > 0) %>%
  mutate(
    pct_GB = 100 * n_GB / n_set[["GB-specific"]],
    pct_GJ = 100 * n_GJ / n_set[["GJ-specific"]],
    p_fisher = map2_dbl(n_GB, n_GJ, ~ fisher.test(matrix(
      c(.x, n_set[["GB-specific"]] - .x, .y, n_set[["GJ-specific"]] - .y), nrow = 2))$p.value),
    fdr = p.adjust(p_fisher, method = "BH"),
    Higher_In = case_when(
      p_fisher < 0.05 & pct_GB > pct_GJ ~ "GB-specific",
      p_fisher < 0.05 & pct_GJ > pct_GB ~ "GJ-specific",
      TRUE                              ~ "neither"
    )
  ) %>%
  left_join(pwy_desc, by = "pathway") %>%
  left_join(pwy_linda, by = "pathway") %>%
  mutate(
    across(c(Community_GB, Community_GJ), ~ replace_na(.x, "not tested")),
    label = paste0(str_trunc(coalesce(description, pathway), 50), " [", pathway, "]")
  ) %>%
  arrange(p_fisher, desc(abs(pct_GB - pct_GJ)))

n_p05 <- sum(comparison$p_fisher < 0.05)
n_fdr <- sum(comparison$fdr < 0.05)
cat(sprintf("\n%d pathways are predicted in at least one marker ASV. %d differ between the sets at p < 0.05 (Fisher, unadjusted); %d at FDR < 0.05.\n",
            nrow(comparison), n_p05, n_fdr))
print(comparison %>% filter(p_fisher < 0.05) %>%
        select(pathway, description, n_GB, n_GJ, p_fisher, fdr, Community_GB, Community_GJ), width = Inf)

write_tsv(comparison %>% select(-label), file.path(out_dir, "pathway_comparison_GB_vs_GJ.tsv"))

presence_wide <- pwy_by_asv %>%
  transmute(pathway, DisplayName, Present = as.integer(Present)) %>%
  pivot_wider(names_from = DisplayName, values_from = Present) %>%
  left_join(pwy_desc, by = "pathway") %>%
  relocate(description, .after = pathway)
write_tsv(presence_wide, file.path(out_dir, "marker_pathway_presence.tsv"))
write_tsv(asv_summary, file.path(out_dir, "marker_pathway_counts.tsv"))

# ==============================================================================
# 3. Figure 1: share of GB- vs GJ-specific markers carrying each pathway
# ==============================================================================
bubble_df <- comparison %>% count(pct_GB, pct_GJ, name = "n_pathways")
sig_df <- comparison %>% filter(p_fisher < 0.05)

p_share <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_point(data = bubble_df, aes(x = pct_GJ, y = pct_GB, size = n_pathways),
             shape = 21, fill = "gray80", color = "gray45", stroke = 0.3, alpha = 0.7) +
  geom_point(data = sig_df, aes(x = pct_GJ, y = pct_GB, fill = Higher_In),
             shape = 21, size = 3, color = "black", stroke = 0.4) +
  geom_text_repel(data = sig_df, aes(x = pct_GJ, y = pct_GB, label = str_trunc(coalesce(description, pathway), 40)),
                  size = 3, family = font_family, max.overlaps = 30, box.padding = 0.4,
                  segment.size = 0.3, segment.color = "gray40", seed = 42) +
  scale_fill_manual(values = palette_set, name = "p < 0.05, higher in") +
  scale_size_continuous(range = c(1.5, 8), name = "Pathways") +
  scale_x_continuous(limits = c(-2, 102), breaks = seq(0, 100, 25)) +
  scale_y_continuous(limits = c(-2, 102), breaks = seq(0, 100, 25)) +
  coord_equal() +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold", size = 9.5)
  ) +
  labs(
    title = "Predicted MetaCyc Pathways of the RF Marker ASVs",
    subtitle = sprintf("%d pathways; %d differ at p < 0.05 (Fisher, unadjusted), %d at FDR < 0.05",
                       nrow(comparison), n_p05, n_fdr),
    x = sprintf("GJ-specific markers carrying the pathway (%%, n = %d)", n_set[["GJ-specific"]]),
    y = sprintf("GB-specific markers carrying the pathway (%%, n = %d)", n_set[["GB-specific"]])
  )

print(p_share)
ggsave(file.path(out_dir, "RF_markers_pathway_share.png"), p_share,
       width = 1400 / 192, height = 1400 / 192, dpi = 192)

# ==============================================================================
# 4. Figure 2: presence of the differing pathways in each marker ASV
# ==============================================================================
if (n_p05 >= 2) {
  col_info <- markers %>% arrange(Marker_Set, Phylum, desc(MeanDecreaseAccuracy))

  mat <- pwy_by_asv %>%
    filter(pathway %in% sig_df$pathway) %>%
    transmute(pathway, Feature.ID, Present = as.integer(Present)) %>%
    pivot_wider(names_from = Feature.ID, values_from = Present) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  mat <- mat[, col_info$Feature.ID, drop = FALSE]
  colnames(mat) <- col_info$DisplayName
  rownames(mat) <- sig_df$label[match(rownames(mat), sig_df$pathway)]

  annotation_col <- data.frame(Marker = col_info$Marker_Set, Phylum = col_info$Phylum,
                               row.names = col_info$DisplayName)
  annotation_row <- sig_df %>%
    transmute(label, `Community GB` = Community_GB, `Community GJ` = Community_GJ) %>%
    column_to_rownames("label")

  phyla <- sort(unique(col_info$Phylum))
  community_colors <- c("higher in treated" = "#B03A2E", "lower in treated" = "#1F618D",
                        "no change" = "gray85", "not tested" = "white")
  annotation_colors <- list(
    Marker = palette_set,
    Phylum = setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(phyla)), phyla),
    `Community GB` = community_colors,
    `Community GJ` = community_colors
  )

  p_heat <- pheatmap(
    mat,
    color             = c("gray95", "gray20"),
    breaks            = c(-0.5, 0.5, 1.5),
    legend_breaks     = c(0, 1),
    legend_labels     = c("absent", "present"),
    cluster_rows      = nrow(mat) > 2,
    cluster_cols      = FALSE,
    gaps_col          = n_set[["GB-specific"]],
    annotation_col    = annotation_col,
    annotation_row    = annotation_row,
    annotation_colors = annotation_colors,
    border_color      = "white",
    fontsize          = 9,
    fontsize_row      = 8.5,
    fontsize_col      = 8.5,
    angle_col         = "90",
    fontfamily        = font_family,
    main              = "Pathways that differ between GB- and GJ-specific marker ASVs (p < 0.05, unadjusted)",
    silent            = TRUE
  )

  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  png(file.path(out_dir, "RF_markers_pathway_heatmap.png"), width = 2800, height = 1400, res = 192)
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  dev.off()
}
