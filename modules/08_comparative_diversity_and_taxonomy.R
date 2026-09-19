#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 08_comparative_diversity_and_taxonomy.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Reclassified 4-Group Dataset (N=134; Sensitivity Model)
# Platform: R (Native ARM64 / macOS)
# Description: Evaluates pairwise and 3-group alpha/beta diversity contrasts
#              (Control-only, Treatment-only, GB-Control baseline, GJ-Control baseline)
#              and computes stacked relative abundance profiles at Phylum and Genus levels.
# ==============================================================================

suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(rstatix)
  library(RColorBrewer)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing Comparative Diversity and Taxonomy Module ===")

# ==============================================================================
# 0. Load Curated Metadata & Core Vectors
# ==============================================================================
meta_file <- "metadata_P_generation_corrected_final.tsv"
if (!file.exists(meta_file)) {
  stop("[-] ERROR: Required metadata '", meta_file, "' not found.")
}

meta_base <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types", Domain == "Bacteria")

obs_features <- read_qza("07_diversity_16S/full/observed_features_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Observed_Features = observed_features)

shannon <- read_qza("07_diversity_16S/full/shannon_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Shannon = shannon_entropy)

faith_pd <- read_qza("07_diversity_16S/full/faith_pd_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Faith_PD = faith_pd)

rarefied_table <- read_qza("07_diversity_16S/full/rarefied_table.qza")$data
simpson_vec <- diversity(t(rarefied_table), index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

# Generic PCoA Runner
run_pcoa_contrast <- function(dist_path, metric_name, meta_sub, pal, title_str) {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  common <- intersect(meta_sub$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common, common])

  meta_sub <- meta_sub %>%
    filter(`sample-id` %in% common) %>%
    arrange(match(`sample-id`, common))

  set.seed(42)
  perm_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  disp_res <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(disp_res, permutations = 999)$tab$`Pr(>F)`[1]

  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_exp <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f",
    perm_res$R2[1],
    ifelse(perm_res$`Pr(>F)`[1] < 0.001, "< 0.001", sprintf("%.3f", perm_res$`Pr(>F)`[1])),
    disp_p
  )

  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, color = NA, show.legend = FALSE) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.8, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
             size = 3.3, fontface = "bold", family = font_family) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
    )
}

# ==============================================================================
# 1. Final Control-Only Diversity (GB_Control vs. GJ_Control)
# ==============================================================================
message("--- Step 1: Control-Only Diversity Profiling ---")

meta_ctrl <- meta_base %>%
  filter(Treatment_Group %in% c("GB_Control", "GJ_Control")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GJ_Control")))

ctrl_pal <- c("GB_Control" = "#2B5C8F", "GJ_Control" = "#D95F02")

alpha_ctrl <- meta_ctrl %>%
  inner_join(obs_features, by = "sample-id") %>%
  inner_join(shannon, by = "sample-id") %>%
  inner_join(simpson_vec, by = "sample-id") %>%
  inner_join(faith_pd, by = "sample-id") %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

p_alpha_ctrl <- ggplot(alpha_ctrl, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = ctrl_pal) +
  scale_color_manual(values = c("GB_Control" = "#152E48", "GJ_Control" = "#7A3501")) +
  scale_x_discrete(labels = c("GB_Control" = "GB_Control (N=33)", "GJ_Control" = "GJ_Control (N=34)")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", size = 4, family = font_family) +
  theme_bw(base_size = 13, base_family = font_family) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "6. M(B), M(J) - Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/Bacteria_16S_Final_Control_Alpha_Diversity.png", p_alpha_ctrl, width = 8.5, height = 7, dpi = 300)

p_bc_ctrl  <- run_pcoa_contrast("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_ctrl, ctrl_pal, "Controls")
p_wuf_ctrl <- run_pcoa_contrast("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_ctrl, ctrl_pal, "Controls")
p_uuf_ctrl <- run_pcoa_contrast("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_ctrl, ctrl_pal, "Controls")
p_jac_ctrl <- run_pcoa_contrast("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_ctrl, ctrl_pal, "Controls")

p_beta_ctrl <- annotate_figure(
  ggarrange(p_bc_ctrl, p_wuf_ctrl, p_uuf_ctrl, p_jac_ctrl, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("6. M(B), M(J) - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/Bacteria_16S_Final_Control_Beta_Diversity.png", p_beta_ctrl, width = 9.5, height = 9, dpi = 300)

# ==============================================================================
# 2. Final Treatment-Only Diversity (GB_Treatment vs. GJ_Treatment)
# ==============================================================================
message("--- Step 2: Treatment-Only Diversity Profiling ---")

meta_treat <- meta_base %>%
  filter(Treatment_Group %in% c("GB_Treatment", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Treatment", "GJ_Treatment")))

treat_pal <- c("GB_Treatment" = "#41B6C4", "GJ_Treatment" = "#FE9929")

alpha_treat <- meta_treat %>%
  inner_join(obs_features, by = "sample-id") %>%
  inner_join(shannon, by = "sample-id") %>%
  inner_join(simpson_vec, by = "sample-id") %>%
  inner_join(faith_pd, by = "sample-id") %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

p_alpha_treat <- ggplot(alpha_treat, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = treat_pal) +
  scale_color_manual(values = c("GB_Treatment" = "#1F6872", "GJ_Treatment" = "#9A5B10")) +
  scale_x_discrete(labels = c("GB_Treatment" = "GB_Treatment (N=33)", "GJ_Treatment" = "GJ_Treatment (N=34)")) +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", size = 4, family = font_family) +
  theme_bw(base_size = 13, base_family = font_family) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(face = "bold", color = "black", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "7. B, J - Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/Bacteria_16S_Final_Treatment_Alpha_Diversity.png", p_alpha_treat, width = 8.5, height = 7, dpi = 300)

p_bc_treat  <- run_pcoa_contrast("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_treat, treat_pal, "Treatments")
p_wuf_treat <- run_pcoa_contrast("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_treat, treat_pal, "Treatments")
p_uuf_treat <- run_pcoa_contrast("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_treat, treat_pal, "Treatments")
p_jac_treat <- run_pcoa_contrast("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_treat, treat_pal, "Treatments")

p_beta_treat <- annotate_figure(
  ggarrange(p_bc_treat, p_wuf_treat, p_uuf_treat, p_jac_treat, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("7. B, J - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/Bacteria_16S_Final_Treatment_Beta_Diversity.png", p_beta_treat, width = 9.5, height = 9, dpi = 300)

# ==============================================================================
# 3. Final 3-Group Diversity Models (GB_Control vs. Inocula)
# ==============================================================================
message("--- Step 3: Tripartite Diversity Models ---")

meta_3g_gb <- meta_base %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Treatment")))

pal_3g_gb <- c("GB_Control" = "#2B5C8F", "GB_Treatment" = "#41B6C4", "GJ_Treatment" = "#FE9929")

alpha_3g_gb <- meta_3g_gb %>%
  inner_join(obs_features, by = "sample-id") %>%
  inner_join(shannon, by = "sample-id") %>%
  inner_join(simpson_vec, by = "sample-id") %>%
  inner_join(faith_pd, by = "sample-id") %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

p_alpha_3g_gb <- ggplot(alpha_3g_gb, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = pal_3g_gb) +
  scale_color_manual(values = pal_3g_gb) +
  scale_x_discrete(labels = c("GB_Control" = "GB_Control (N=33)", "GB_Treatment" = "GB_Treatment (N=33)", "GJ_Treatment" = "GJ_Treatment (N=34)")) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 3.8, family = font_family) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", color = "black", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "7. M(B), B, J - Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/Bacteria_16S_Final_GBctrl_GBtreat_GJtreat_Alpha_Diversity.png", p_alpha_3g_gb, width = 9.5, height = 7.5, dpi = 300)

p_bc_3g_gb  <- run_pcoa_contrast("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_3g_gb, pal_3g_gb, "3-Group")
p_wuf_3g_gb <- run_pcoa_contrast("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_3g_gb, pal_3g_gb, "3-Group")
p_uuf_3g_gb <- run_pcoa_contrast("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_3g_gb, pal_3g_gb, "3-Group")
p_jac_3g_gb <- run_pcoa_contrast("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_3g_gb, pal_3g_gb, "3-Group")

p_beta_3g_gb <- annotate_figure(
  ggarrange(p_bc_3g_gb, p_wuf_3g_gb, p_uuf_3g_gb, p_jac_3g_gb, ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom"),
  top = text_grob("7. M(B), B, J - Beta Diversity", face = "bold", size = 14, family = font_family)
)
ggsave("assets/Bacteria_16S_Final_GBctrl_GBtreat_GJtreat_Beta_Diversity.png", p_beta_3g_gb, width = 9.5, height = 9, dpi = 300)

# ==============================================================================
# 4. Taxonomic Composition Profiling (Phylum & Genus Levels)
# ==============================================================================
message("--- Step 4: Taxonomic Aggregation and Stacked Composition ---")

counts_raw <- read_qza("07_diversity_16S/full/rarefied_table.qza")$data
tax_clean <- parse_taxonomy(read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data) %>%
  rownames_to_column("Feature.ID")

common_samples <- intersect(meta_base$`sample-id`, colnames(counts_raw))
counts_sub <- counts_raw[, common_samples]

rel_abund_long <- counts_sub %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  left_join(tax_clean, by = "Feature.ID") %>%
  inner_join(meta_base, by = "sample-id")

# 4.1 Phylum Aggregation
phylum_sample <- rel_abund_long %>%
  mutate(Phylum = if_else(is.na(Phylum) | Phylum == "", "Unassigned", Phylum)) %>%
  group_by(`sample-id`, Treatment_Group, Phylum) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

top10_phyla <- phylum_sample %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Phylum != "Unassigned") %>%
  slice_head(n = 10) %>%
  pull(Phylum)

phylum_summary <- phylum_sample %>%
  mutate(Phylum_Group = if_else(Phylum %in% top10_phyla, Phylum, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Phylum_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Phylum_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Phylum_Group = factor(Phylum_Group, levels = rev(c(top10_phyla, "Other"))))

taxa_colors_phylum <- setNames(c(brewer.pal(10, "Paired"), "#B0B0B0"), c(top10_phyla, "Other"))

p_phylum <- ggplot(phylum_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Phylum_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_phylum) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control\n(N=33)",
    "GB_Treatment" = "GB_Treatment\n(N=33)",
    "GJ_Control"   = "GJ_Control\n(N=34)",
    "GJ_Treatment" = "GJ_Treatment\n(N=34)"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    panel.grid = element_blank()
  ) +
  labs(title = "8. Relative Abundance (Phylum Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/Bacteria_16S_Final_Phylum_Relative_Abundance.png", p_phylum, width = 8.5, height = 6.5, dpi = 300)

# 4.2 Genus Aggregation
genus_sample <- rel_abund_long %>%
  mutate(Genus = if_else(is.na(Genus) | Genus == "" | Genus == "g__", "Unassigned", Genus)) %>%
  group_by(`sample-id`, Treatment_Group, Genus) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

top12_genera <- genus_sample %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  arrange(desc(MeanAbund)) %>%
  filter(Genus != "Unassigned") %>%
  slice_head(n = 12) %>%
  pull(Genus)

genus_summary <- genus_sample %>%
  mutate(Genus_Group = if_else(Genus %in% top12_genera, Genus, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Genus_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Genus_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Genus_Group = factor(Genus_Group, levels = rev(c(top12_genera, "Other"))))

taxa_colors_genus <- setNames(c(brewer.pal(12, "Set3"), "#B0B0B0"), c(top12_genera, "Other"))

p_genus <- ggplot(genus_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Genus_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = taxa_colors_genus) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control\n(N=33)",
    "GB_Treatment" = "GB_Treatment\n(N=33)",
    "GJ_Control"   = "GJ_Control\n(N=34)",
    "GJ_Treatment" = "GJ_Treatment\n(N=34)"
  )) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 9, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 11),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    panel.grid = element_blank()
  ) +
  labs(title = "8. Relative Abundance (Genus Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/Bacteria_16S_Final_Genus_Relative_Abundance.png", p_genus, width = 8.5, height = 6.5, dpi = 300)

message("=== [", format(Sys.time(), "%T"), "] Comparative Diversity and Taxonomy Complete ===")