#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 11_publication_diversity_figures.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1-3, GJ1-2; N=110 Locked Cohort)
# Platform: R (Native ARM64 / macOS)
# Description: Generates publication-ready Alpha and Beta diversity figures
#              (PCoA and NMDS across Bray-Curtis, Weighted/Unweighted UniFrac,
#              and Jaccard), exporting both composite panels and individual
#              standalone plots directly to assets/ (Figures 046-066).
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Dependencies, Global Typography & Directory Configuration
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(rstatix)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing Publication Diversity Figures Module ===")

# Verify metadata presence
meta_file <- "metadata_P_generation.tsv"
if (!file.exists(meta_file)) {
  stop("[-] ERROR: Required metadata '", meta_file, "' not found.")
}

# ------------------------------------------------------------------------------
# 1. Global Aesthetic Mappings (Pruned 5-Batch Cohort)
# ------------------------------------------------------------------------------
# High-contrast palette: Brighter tints for Controls, saturated dark tones for Treatments
palette_4 <- c(
  "GB_Control"   = "#4A90E2",
  "GB_Treatment" = "#1A365D",
  "GJ_Control"   = "#F87171",
  "GJ_Treatment" = "#881337"
)

# Geographic Lineage Shapes: Circle (16) for Gijang B, Triangle (17) for Gyeongju
shape_4 <- c(
  "GB_Control"   = 16,
  "GB_Treatment" = 16,
  "GJ_Control"   = 17,
  "GJ_Treatment" = 17
)

palette_2 <- c(
  "GB_Treatment" = "#1A365D",
  "GJ_Treatment" = "#881337"
)

shape_2 <- c(
  "GB_Treatment" = 16,
  "GJ_Treatment" = 17
)

# Core Metrics Input Directory
cm_dir <- "07_diversity_16S/core_metrics"

# ==============================================================================
# 2. 4-Group Alpha Diversity Profiling (Figure 046)
# ==============================================================================
message("--- Step 1: Generating 4-Group Alpha Diversity Figures ---")

meta_4g <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment"))
  )

obs_features <- read_qza(file.path(cm_dir, "observed_features_vector.qza"))$data %>%
  rownames_to_column("sample-id") %>%
  rename(Observed_Features = observed_features)

shannon <- read_qza(file.path(cm_dir, "shannon_vector.qza"))$data %>%
  rownames_to_column("sample-id") %>%
  rename(Shannon = shannon_entropy)

faith_pd <- read_qza(file.path(cm_dir, "faith_pd_vector.qza"))$data %>%
  rownames_to_column("sample-id") %>%
  rename(Faith_PD = faith_pd)

rarefied_table <- read_qza(file.path(cm_dir, "rarefied_table.qza"))$data
simpson_vec <- diversity(t(rarefied_table), index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

alpha_merged <- meta_4g %>%
  inner_join(obs_features, by = "sample-id") %>%
  inner_join(shannon, by = "sample-id") %>%
  inner_join(simpson_vec, by = "sample-id") %>%
  inner_join(faith_pd, by = "sample-id")

alpha_long <- alpha_merged %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

kw_results <- alpha_long %>% group_by(Metric) %>% kruskal_test(Value ~ Treatment_Group)
pairwise_results <- alpha_long %>% group_by(Metric) %>% wilcox_test(Value ~ Treatment_Group, p.adjust.method = "fdr")

cat("=== Global Kruskal-Wallis Test Results (N=110) ===\n")
print(kw_results)

p_046 <- ggplot(alpha_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_4, name = "Treatment Group") +
  scale_color_manual(values = palette_4, name = "Treatment Group") +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 3.8, family = font_family) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
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
  labs(title = "1. M(B), B, M(J), J - Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/046_16S_4Group_Alpha_Diversity.png", p_046, width = 9.5, height = 7.5, dpi = 300)

# ==============================================================================
# 3. 4-Group Beta Diversity: PCoA (Figures 047 & 048-051)
# ==============================================================================
message("--- Step 2: Generating 4-Group PCoA Ordinations ---")

generate_pcoa <- function(dist_path, metric_name, meta_df, pal, shp, standalone = FALSE) {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  common <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common, common])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common) %>%
    arrange(match(`sample-id`, common))

  set.seed(42)
  perm_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  r2 <- perm_res$R2[1]
  p_val <- perm_res$`Pr(>F)`[1]

  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_exp <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f",
    r2,
    ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    disp_p
  )

  pt_size <- if (standalone) 3.2 else 2.8

  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group, shape = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = pt_size, alpha = 0.85) +
    scale_color_manual(values = pal, name = "Group") +
    scale_fill_manual(values = pal, name = "Group") +
    scale_shape_manual(values = shp, name = "Group") +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      text = element_text(family = font_family),
      legend.position = if (standalone) "right" else "bottom",
      legend.title = if (standalone) element_text(face = "bold", size = 11) else element_blank(),
      legend.text = element_text(size = 10),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = if (standalone) 14 else 12),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
    )

  if (!standalone) {
    p <- p + annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
                      size = 3.5, fontface = "bold", family = font_family)
  }

  return(p)
}

# Composite subplots (with stats)
p_4g_bc  <- generate_pcoa(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_4g, palette_4, shape_4, FALSE)
p_4g_wuf <- generate_pcoa(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_4g, palette_4, shape_4, FALSE)
p_4g_uuf <- generate_pcoa(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_4g, palette_4, shape_4, FALSE)
p_4g_jac <- generate_pcoa(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_4g, palette_4, shape_4, FALSE)

p_047 <- annotate_figure(
  ggarrange(p_4g_bc, p_4g_wuf, p_4g_uuf, p_4g_jac, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("1. M(B), B, M(J), J - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/047_16S_4Group_Beta_Diversity.png", p_047, width = 9.5, height = 9, dpi = 300)

# Standalone clean plots (no stat overlay, legend on right)
ggsave("assets/048_16S_4Group_Beta_Diversity_BC.png",
       generate_pcoa(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/049_16S_4Group_Beta_Diversity_WUD.png",
       generate_pcoa(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/050_16S_4Group_Beta_Diversity_UUF.png",
       generate_pcoa(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/051_16S_4Group_Beta_Diversity_JD.png",
       generate_pcoa(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

# ==============================================================================
# 4. 4-Group Beta Diversity: NMDS (Figures 052 & 053-056)
# ==============================================================================
message("--- Step 3: Generating 4-Group NMDS Ordinations ---")

generate_nmds <- function(dist_path, metric_name, meta_df, pal, shp, standalone = FALSE) {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  common <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common, common])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common) %>%
    arrange(match(`sample-id`, common))

  set.seed(42)
  perm_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  r2 <- perm_res$R2[1]
  p_val <- perm_res$`Pr(>F)`[1]

  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  set.seed(42)
  nmds_res <- metaMDS(dist_sub, k = 2, trymax = 100, trace = FALSE)
  stress_val <- nmds_res$stress

  nmds_df <- as.data.frame(nmds_res$points) %>%
    rownames_to_column("sample-id") %>%
    rename(NMDS1 = MDS1, NMDS2 = MDS2) %>%
    inner_join(meta_sub, by = "sample-id")

  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f\nStress = %.3f",
    r2,
    ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val)),
    disp_p,
    stress_val
  )

  pt_size <- if (standalone) 3.2 else 2.8

  p <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Treatment_Group, fill = Treatment_Group, shape = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7, show.legend = FALSE) +
    geom_point(size = pt_size, alpha = 0.85) +
    scale_color_manual(values = pal, name = "Group") +
    scale_fill_manual(values = pal, name = "Group") +
    scale_shape_manual(values = shp, name = "Group") +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      text = element_text(family = font_family),
      legend.position = if (standalone) "right" else "bottom",
      legend.title = if (standalone) element_text(face = "bold", size = 11) else element_blank(),
      legend.text = element_text(size = 10),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = if (standalone) 14 else 12),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = sprintf("NMDS - %s Distance", metric_name),
      x = "NMDS 1",
      y = "NMDS 2"
    )

  if (!standalone) {
    p <- p + annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
                      size = 3.3, fontface = "bold", family = font_family)
  }

  return(p)
}

# Composite subplots (with stats)
p_4g_nmds_bc  <- generate_nmds(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_4g, palette_4, shape_4, FALSE)
p_4g_nmds_wuf <- generate_nmds(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_4g, palette_4, shape_4, FALSE)
p_4g_nmds_uuf <- generate_nmds(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_4g, palette_4, shape_4, FALSE)
p_4g_nmds_jac <- generate_nmds(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_4g, palette_4, shape_4, FALSE)

p_052 <- annotate_figure(
  ggarrange(p_4g_nmds_bc, p_4g_nmds_wuf, p_4g_nmds_uuf, p_4g_nmds_jac, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("1. M(B), B, M(J), J - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/052_16S_4Group_Beta_Diversity(NMDS).png", p_052, width = 9.5, height = 9, dpi = 300)

# Standalone clean plots
ggsave("assets/053_16S_4Group_Beta_Diversity(NMDS)_BC.png",
       generate_nmds(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/054_16S_4Group_Beta_Diversity(NMDS)_WUF.png",
       generate_nmds(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/055_16S_4Group_Beta_Diversity(NMDS)_UUF.png",
       generate_nmds(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/056_16S_4Group_Beta_Diversity(NMDS)_JD.png",
       generate_nmds(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_4g, palette_4, shape_4, TRUE),
       width = 7.5, height = 6, dpi = 300)

# ==============================================================================
# 5. Treatment-Only Beta Diversity: PCoA (Figures 057 & 058-061)
# ==============================================================================
message("--- Step 4: Generating Treatment-Only PCoA Ordinations ---")

meta_treat <- meta_4g %>%
  filter(Treatment_Group %in% c("GB_Treatment", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Treatment", "GJ_Treatment")))

p_tr_bc  <- generate_pcoa(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_treat, palette_2, shape_2, FALSE)
p_tr_wuf <- generate_pcoa(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_treat, palette_2, shape_2, FALSE)
p_tr_uuf <- generate_pcoa(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_treat, palette_2, shape_2, FALSE)
p_tr_jac <- generate_pcoa(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_treat, palette_2, shape_2, FALSE)

p_057 <- annotate_figure(
  ggarrange(p_tr_bc, p_tr_wuf, p_tr_uuf, p_tr_jac, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("2. B, J - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/057_16S_Treatment_Beta_Diversity.png", p_057, width = 9.5, height = 9, dpi = 300)

# Standalone clean treatment plots
ggsave("assets/058_16S_Treatment_Beta_Diversity_BC.png",
       generate_pcoa(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/059_16S_Treatment_Beta_Diversity_WUF.png",
       generate_pcoa(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/060_16S_Treatment_Beta_Diversity_UUF.png",
       generate_pcoa(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/061_16S_Treatment_Beta_Diversity_JD.png",
       generate_pcoa(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

# ==============================================================================
# 6. Treatment-Only Beta Diversity: NMDS (Figures 062 & 063-066)
# ==============================================================================
message("--- Step 5: Generating Treatment-Only NMDS Ordinations ---")

p_tr_nmds_bc  <- generate_nmds(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_treat, palette_2, shape_2, FALSE)
p_tr_nmds_wuf <- generate_nmds(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_treat, palette_2, shape_2, FALSE)
p_tr_nmds_uuf <- generate_nmds(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_treat, palette_2, shape_2, FALSE)
p_tr_nmds_jac <- generate_nmds(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_treat, palette_2, shape_2, FALSE)

p_062 <- annotate_figure(
  ggarrange(p_tr_nmds_bc, p_tr_nmds_wuf, p_tr_nmds_uuf, p_tr_nmds_jac, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("2. B, J - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/062_16S_Treatment_Beta_Diversity(NMDS).png", p_062, width = 9.5, height = 9, dpi = 300)

# Standalone clean treatment NMDS plots
ggsave("assets/063_16S_Treatment_Beta_Diversity(NMDS)_BC.png",
       generate_nmds(file.path(cm_dir, "bray_curtis_distance_matrix.qza"), "Bray-Curtis", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/064_16S_Treatment_Beta_Diversity(NMDS)_WUF.png",
       generate_nmds(file.path(cm_dir, "weighted_unifrac_distance_matrix.qza"), "Weighted UniFrac", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/065_16S_Treatment_Beta_Diversity(NMDS)_UUF.png",
       generate_nmds(file.path(cm_dir, "unweighted_unifrac_distance_matrix.qza"), "Unweighted UniFrac", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

ggsave("assets/066_16S_Treatment_Beta_Diversity(NMDS)_JD.png",
       generate_nmds(file.path(cm_dir, "jaccard_distance_matrix.qza"), "Jaccard", meta_treat, palette_2, shape_2, TRUE),
       width = 7.5, height = 6, dpi = 300)

message("=== [", format(Sys.time(), "%T"), "] Publication Diversity Figures Successfully Exported to assets/ ===")