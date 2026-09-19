#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 06_alpha_beta_diversity.R
# Target: 16S rRNA & ITS2 Community Diversity Profiling
# Platform: R (Native ARM64 / macOS)
# Description: Evaluates alpha diversity (Observed Features, Shannon, Simpson,
#              Faith's PD, Pielou) and beta diversity (Bray-Curtis, Jaccard,
#              Weighted/Unweighted UniFrac) with PERMANOVA and betadisper models.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment & Package Prerequisites
# ------------------------------------------------------------------------------
required_bioc <- c("phyloseq", "microbiome")
required_cran <- c("tidyverse", "vegan", "ggpubr", "rstatix", "patchwork", "devtools")

for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in required_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, update = FALSE)
}

if (!requireNamespace("qiime2R", quietly = TRUE)) {
  devtools::install_github("jbisanz/qiime2R")
}

suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(rstatix)
})

# Typography setup
font_family <- "Times New Roman"

# Ensure assets directory exists for export
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing Diversity Modeling Module ===")


# ==============================================================================
# 1. Validation of Control Invariance: GB_Control vs. GJ_Control (16S)
# ==============================================================================
message("--- Step 1: Evaluating Control Invariance (GB_Control vs. GJ_Control) ---")

meta_ctrl <- read_tsv("metadata_P_generation.tsv", show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Sample_Type == "Standard",
    Treatment_Group %in% c("GB_Control", "GJ_Control")
  ) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GJ_Control")))

# 1.1 Alpha Diversity Data Ingestion & Calculation
obs_ctrl <- read_qza("07_diversity_16S/standard/observed_features_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Observed_Features = observed_features)

shannon_ctrl <- read_qza("07_diversity_16S/standard/shannon_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Shannon = shannon_entropy)

faith_ctrl <- read_qza("07_diversity_16S/standard/faith_pd_vector.qza")$data %>%
  rownames_to_column("sample-id") %>%
  rename(Faith_PD = faith_pd)

rarefied_tbl_ctrl <- read_qza("07_diversity_16S/standard/rarefied_table.qza")$data
simpson_ctrl <- diversity(t(rarefied_tbl_ctrl), index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

alpha_ctrl_long <- meta_ctrl %>%
  inner_join(obs_ctrl, by = "sample-id") %>%
  inner_join(shannon_ctrl, by = "sample-id") %>%
  inner_join(simpson_ctrl, by = "sample-id") %>%
  inner_join(faith_ctrl, by = "sample-id") %>%
  pivot_longer(
    cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

p_alpha_ctrl <- ggplot(alpha_ctrl_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black") +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("GB_Control" = "#2B5C8F", "GJ_Control" = "#D95F02")) +
  scale_color_manual(values = c("GB_Control" = "#152E48", "GJ_Control" = "#7A3501")) +
  scale_x_discrete(labels = c("GB_Control" = "Gijang_B Control", "GJ_Control" = "Gyeongju Control")) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.x.npc = "center",
    size = 4,
    family = font_family
  ) +
  theme_bw(base_size = 13, base_family = font_family) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", size = 11, family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = font_family),
    axis.text.x = element_text(face = "bold", color = "black", family = font_family),
    axis.text.y = element_text(family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "1. M(B), M(J) - Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/015_16S_Control_Alpha_Diversity.png", p_alpha_ctrl, width = 8, height = 8, dpi = 300)

# 1.2 Beta Diversity Ordination Function for Controls
analyze_ctrl_beta <- function(dist_path, metric_name, meta_df) {
  dist_mat <- as.matrix(read_qza(dist_path)$data)
  common_samples <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

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

  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, color = NA) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.8) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_manual(values = c("GB_Control" = "#2B5C8F", "GJ_Control" = "#D95F02"),
                       labels = c("GB_Control" = "Gijang_B Control", "GJ_Control" = "Gyeongju Control")) +
    scale_fill_manual(values = c("GB_Control" = "#2B5C8F", "GJ_Control" = "#D95F02"),
                      labels = c("GB_Control" = "Gijang_B Control", "GJ_Control" = "Gyeongju Control")) +
    theme_bw(base_size = 13, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
             size = 3.8, fontface = "bold", family = font_family) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
    )

  return(p)
}

p_ctrl_bc  <- analyze_ctrl_beta("07_diversity_16S/standard/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_ctrl)
p_ctrl_wuf <- analyze_ctrl_beta("07_diversity_16S/standard/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_ctrl)
p_ctrl_uuf <- analyze_ctrl_beta("07_diversity_16S/standard/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_ctrl)
p_ctrl_jac <- analyze_ctrl_beta("07_diversity_16S/standard/jaccard_distance_matrix.qza", "Jaccard", meta_ctrl)

combined_ctrl_pcoa <- ggarrange(p_ctrl_bc, p_ctrl_wuf, p_ctrl_uuf, p_ctrl_jac, ncol = 2, nrow = 2,
                                common.legend = TRUE, legend = "bottom")
final_ctrl_pcoa <- annotate_figure(
  combined_ctrl_pcoa,
  top = text_grob("1. M(B), M(J) - Beta Diversity", face = "bold", size = 15, family = font_family)
)

ggsave("assets/016_16S_Control_Beta_Diversity.png", final_ctrl_pcoa, width = 9, height = 9, dpi = 300)


# ==============================================================================
# 2. 16S Full Cohort Diversity Profiling (Including Extra x/y Samples)
# ==============================================================================
message("--- Step 2: Evaluating 16S Full Cohort (5-Group Sensitivity Model) ---")

meta_full <- read_tsv("metadata_P_generation.tsv", show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(Domain == "Bacteria") %>%
  mutate(
    Treatment_Group = if_else(Treatment_Group == "GJ_Treatment_Extra", "GJ_Extra", Treatment_Group),
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment", "GJ_Extra"))
  )

palette_5 <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4",
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929",
  "GJ_Extra"     = "#7A0177"
)

group_counts_full <- table(meta_full$Treatment_Group)
labels_5 <- setNames(paste0(names(group_counts_full), " (N=", as.numeric(group_counts_full), ")"), names(group_counts_full))

# 2.1 Full Alpha Diversity
obs_full <- read_qza("07_diversity_16S/full/observed_features_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Observed_Features = observed_features)
shannon_full <- read_qza("07_diversity_16S/full/shannon_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Shannon = shannon_entropy)
faith_full <- read_qza("07_diversity_16S/full/faith_pd_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Faith_PD = faith_pd)

rarefied_tbl_full <- read_qza("07_diversity_16S/full/rarefied_table.qza")$data
simpson_full <- diversity(t(rarefied_tbl_full), index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

alpha_full_long <- meta_full %>%
  inner_join(obs_full, by = "sample-id") %>%
  inner_join(shannon_full, by = "sample-id") %>%
  inner_join(simpson_full, by = "sample-id") %>%
  inner_join(faith_full, by = "sample-id") %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

p_alpha_full <- ggplot(alpha_full_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_5, labels = labels_5) +
  scale_color_manual(values = palette_5, labels = labels_5) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 3.8, family = font_family) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", size = 11, family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "1. M(B), M(J), B, J, J(Extra) - Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/017_16S_Extra_Alpha_Diversity.png", p_alpha_full, width = 8, height = 8, dpi = 300)

# 2.2 Full Beta Diversity
analyze_full_beta <- function(dist_path, metric_name, meta_df) {
  dist_mat <- as.matrix(read_qza(dist_path)$data)
  common_samples <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

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

  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = palette_5) +
    scale_fill_manual(values = palette_5) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
             size = 3.5, fontface = "bold", family = font_family) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
    )

  return(p)
}

p_full_bc  <- analyze_full_beta("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_full)
p_full_wuf <- analyze_full_beta("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_full)
p_full_uuf <- analyze_full_beta("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_full)
p_full_jac <- analyze_full_beta("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_full)

combined_full_pcoa <- ggarrange(p_full_bc, p_full_wuf, p_full_uuf, p_full_jac, ncol = 2, nrow = 2,
                                common.legend = TRUE, legend = "bottom")
final_full_pcoa <- annotate_figure(
  combined_full_pcoa,
  top = text_grob("1. M(B), M(J), B, J, J(Extra) - Beta Diversity", face = "bold", size = 14, family = font_family)
)

ggsave("assets/018_16S_Extra_Beta_Diversity.png", final_full_pcoa, width = 9, height = 9, dpi = 300)


# ==============================================================================
# 3. Fungal ITS2 Diversity Profiling (Full Cohort)
# ==============================================================================
message("--- Step 3: Evaluating Fungal ITS2 Diversity ---")

meta_its <- read_tsv("metadata_P_generation.tsv", show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(Domain == "Fungi") %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment", "GJ_Treatment_Extra"))
  )

palette_its <- c(
  "GB_Control"         = "#2B5C8F",
  "GB_Treatment"       = "#41B6C4",
  "GJ_Control"         = "#D95F02",
  "GJ_Treatment"       = "#FE9929",
  "GJ_Treatment_Extra" = "#7A0177"
)

# 3.1 ITS Alpha Diversity
obs_its <- read_qza("07_diversity_ITS/full/observed_features_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Observed_Features = observed_features)
shannon_its <- read_qza("07_diversity_ITS/full/shannon_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Shannon = shannon_entropy)
evenness_its <- read_qza("07_diversity_ITS/full/evenness_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Pielou_Evenness = pielou_evenness)

rarefied_tbl_its <- read_qza("07_diversity_ITS/full/rarefied_table.qza")$data
simpson_its <- diversity(t(rarefied_tbl_its), index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

alpha_its_long <- meta_its %>%
  inner_join(obs_its, by = "sample-id") %>%
  inner_join(shannon_its, by = "sample-id") %>%
  inner_join(simpson_its, by = "sample-id") %>%
  inner_join(evenness_its, by = "sample-id") %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Pielou_Evenness"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Pielou_Evenness")))

p_alpha_its <- ggplot(alpha_its_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_its) +
  scale_color_manual(values = palette_its) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 3.8, family = font_family) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", size = 11, family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "1. M(B), M(J), B, J - ITS Alpha Diversity", y = "Calculated Diversity Value")

ggsave("assets/021_ITS_Extra_Alpha_Diversity.png", p_alpha_its, width = 8, height = 8, dpi = 300)

# 3.2 ITS Beta Diversity
analyze_its_beta <- function(dist_path, metric_name, meta_df) {
  dist_mat <- as.matrix(read_qza(dist_path)$data)
  common_samples <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

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

  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.12, level = 0.95, color = NA) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.7) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = palette_its) +
    scale_fill_manual(values = palette_its) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
             size = 3.5, fontface = "bold", family = font_family) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_exp[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_exp[2])
    )

  return(p)
}

p_its_bc  <- analyze_its_beta("07_diversity_ITS/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_its)
p_its_jac <- analyze_its_beta("07_diversity_ITS/full/jaccard_distance_matrix.qza", "Jaccard", meta_its)

combined_its_pcoa <- ggarrange(p_its_bc, p_its_jac, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
final_its_pcoa <- annotate_figure(
  combined_its_pcoa,
  top = text_grob("1. M(B), M(J), B, J - ITS Beta Diversity", face = "bold", size = 15, family = font_family)
)

ggsave("assets/022_ITS_Extra_Beta_Diversity.png", final_its_pcoa, width = 9, height = 4.5, dpi = 300)

message("=== [", format(Sys.time(), "%T"), "] Diversity Modeling Completed Successfully ===")