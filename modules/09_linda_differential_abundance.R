#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 09_linda_differential_abundance.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Reclassified 4-Group Dataset (N=134; Sensitivity Model)
# Platform: R (Native ARM64 / macOS)
# Description: Executes LinDA differential abundance modeling across synchronous
#              cohorts (B vs. M(B) and J vs. M(J)), plots volcano distributions,
#              generates Global and Top 30 heatmaps, and filters exclusive rare taxa.
# ==============================================================================

suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(MicrobiomeStat)
  library(pheatmap)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing LinDA Differential Abundance Module ===")

# ==============================================================================
# 0. Load Artifacts & Metadata
# ==============================================================================
meta_file <- "metadata_P_generation_corrected_final.tsv"
if (!file.exists(meta_file)) {
  stop("[-] ERROR: Required metadata '", meta_file, "' not found.")
}

meta_final <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types", Domain == "Bacteria")

counts_raw <- read_qza("07_diversity_16S/full/rarefied_table.qza")$data
tax_clean <- parse_taxonomy(read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data) %>%
  rownames_to_column("Feature.ID")

# ==============================================================================
# 1. LinDA Differential Modeling: GB Cohort [B vs. M(B)]
# ==============================================================================
message("--- Step 1: LinDA Modeling for GB Cohort [B vs. M(B)] ---")

meta_gb <- meta_final %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment"))) %>%
  column_to_rownames("sample-id")

counts_gb <- counts_raw[, rownames(meta_gb)]
rel_abund_gb <- apply(counts_gb, 2, function(x) (x / sum(x)) * 100)

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
  mutate(
    lfc            = log2FoldChange,
    neg_log10_padj = -log10(padj)
  )

volcano_df_gb <- df_gb %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  inner_join(enframe(rowMeans(rel_abund_gb), name = "Feature.ID", value = "mean_rel_abund"), by = "Feature.ID") %>%
  mutate(
    Significance = case_when(
      lfc >= 1.0 & padj < 0.05  ~ "Enriched in B",
      lfc <= -1.0 & padj < 0.05 ~ "Enriched in M(B)",
      TRUE                      ~ "Not Significant"
    ),
    Significance = factor(Significance, levels = c("Enriched in B", "Enriched in M(B)", "Not Significant")),
    Label_Name = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != ""    ~ Genus,
      !is.na(Family) & Family != "f__" & Family != "" ~ paste("Fam.", Family),
      TRUE                                            ~ paste("Phy.", Phylum)
    )
  )

top_labels_gb <- volcano_df_gb %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_max(order_by = neg_log10_padj, n = 8, with_ties = FALSE) %>%
  ungroup()

p_volcano_gb <- ggplot(volcano_df_gb, aes(x = lfc, y = neg_log10_padj, color = Significance)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(aes(size = mean_rel_abund), alpha = 0.75) +
  scale_size_continuous(range = c(1.0, 5.0), name = "Mean Rel. Abund (%)") +
  scale_color_manual(values = c("Enriched in B" = "#008B8B", "Enriched in M(B)" = "#2B5C8F", "Not Significant" = "#C0C0C0")) +
  scale_y_continuous(name = expression(bold("-Log"[10] * " (FDR-adjusted " * italic("p") * "-value)")), expand = expansion(mult = c(0, 0.05))) +
  geom_text_repel(
    data = top_labels_gb, aes(label = Label_Name, color = Significance),
    size = 3.5, fontface = "italic", family = font_family, max.overlaps = 15, show.legend = FALSE
  ) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "9. B vs. M(B) - Volcano Plot (LinDA)", x = expression(bold("Log"[2] * " Fold Change (B / M(B))")))

ggsave("assets/041_16S_B_Volcano_Plot.png", p_volcano_gb, width = 8.5, height = 7.5, dpi = 300)

# ==============================================================================
# 2. LinDA Differential Modeling: GJ Cohort [J vs. M(J)]
# ==============================================================================
message("--- Step 2: LinDA Modeling for GJ Cohort [J vs. M(J)] ---")

meta_gj <- meta_final %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GJ_Control", "GJ_Treatment"))) %>%
  column_to_rownames("sample-id")

counts_gj <- counts_raw[, rownames(meta_gj)]
rel_abund_gj <- apply(counts_gj, 2, function(x) (x / sum(x)) * 100)

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
  mutate(
    lfc            = log2FoldChange,
    neg_log10_padj = -log10(padj)
  )

volcano_df_gj <- df_gj %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  inner_join(enframe(rowMeans(rel_abund_gj), name = "Feature.ID", value = "mean_rel_abund"), by = "Feature.ID") %>%
  mutate(
    Significance = case_when(
      lfc >= 1.0 & padj < 0.05  ~ "Enriched in J",
      lfc <= -1.0 & padj < 0.05 ~ "Enriched in M(J)",
      TRUE                      ~ "Not Significant"
    ),
    Significance = factor(Significance, levels = c("Enriched in J", "Enriched in M(J)", "Not Significant")),
    Label_Name = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != ""    ~ Genus,
      !is.na(Family) & Family != "f__" & Family != "" ~ paste("Fam.", Family),
      TRUE                                            ~ paste("Phy.", Phylum)
    )
  )

top_labels_gj <- volcano_df_gj %>%
  filter(Significance != "Not Significant") %>%
  group_by(Significance) %>%
  slice_max(order_by = neg_log10_padj, n = 8, with_ties = FALSE) %>%
  ungroup()

p_volcano_gj <- ggplot(volcano_df_gj, aes(x = lfc, y = neg_log10_padj, color = Significance)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(aes(size = mean_rel_abund), alpha = 0.75) +
  scale_size_continuous(range = c(1.0, 5.0), name = "Mean Rel. Abund (%)") +
  scale_color_manual(values = c("Enriched in J" = "#FE9929", "Enriched in M(J)" = "#D95F02", "Not Significant" = "#C0C0C0")) +
  scale_y_continuous(name = expression(bold("-Log"[10] * " (FDR-adjusted " * italic("p") * "-value)")), expand = expansion(mult = c(0, 0.05))) +
  geom_text_repel(
    data = top_labels_gj, aes(label = Label_Name, color = Significance),
    size = 3.5, fontface = "italic", family = font_family, max.overlaps = 15, show.legend = FALSE
  ) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10, family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "9. J vs. M(J) - Volcano Plot (LinDA)", x = expression(bold("Log"[2] * " Fold Change (J / M(J))")))

ggsave("assets/042_16S_J_Volcano_Plot.png", p_volcano_gj, width = 8.5, height = 7.5, dpi = 300)

# ==============================================================================
# 3. Global Comparative Heatmap (LinDA Clamped Effect Sizes)
# ==============================================================================
message("--- Step 3: Generating Global Differential Heatmap ---")

combined_linda <- bind_rows(
  df_gb %>% transmute(Feature.ID, Comparison = "B vs. M(B)", lfc, padj),
  df_gj %>% transmute(Feature.ID, Comparison = "J vs. M(J)", lfc, padj)
)

sig_features <- combined_linda %>%
  filter(padj < 0.05, abs(lfc) >= 1.0) %>%
  pull(Feature.ID) %>%
  unique()

lfc_matrix_global <- combined_linda %>%
  filter(Feature.ID %in% sig_features) %>%
  select(Feature.ID, Comparison, lfc) %>%
  pivot_wider(names_from = Comparison, values_from = lfc, values_fill = 0) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

lfc_matrix_global <- lfc_matrix_global[, c("B vs. M(B)", "J vs. M(J)")]
lfc_matrix_global[lfc_matrix_global > 8]  <- 8
lfc_matrix_global[lfc_matrix_global < -8] <- -8

pheatmap(
  mat                      = lfc_matrix_global,
  cluster_rows             = TRUE,
  cluster_cols             = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method        = "complete",
  show_rownames            = FALSE,
  show_colnames            = TRUE,
  color                    = colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100),
  breaks                   = seq(-8, 8, length.out = 101),
  treeheight_row           = 45,
  border_color             = NA,
  fontsize_col             = 11,
  angle_col                = "0",
  fontfamily               = font_family,
  main                     = "10. B vs. M(B), J vs. M(J) - Global Heatmap (LinDA)",
  filename                 = "assets/043_16S_Global_Heatmap.png",
  width                    = 8.5,
  height                   = 8.5
)

# ==============================================================================
# 4. Top 30 Deterministic Biomarker Heatmap
# ==============================================================================
message("--- Step 4: Generating Top 30 Deterministic Heatmap ---")

rel_abund_all <- apply(counts_raw[, meta_final$`sample-id`], 2, function(x) (x / sum(x)) * 100)
mean_rel_abund <- rowMeans(rel_abund_all) %>%
  enframe(name = "Feature.ID", value = "mean_rel_abund")

top30_features <- combined_linda %>%
  filter(padj < 0.05) %>%
  group_by(Feature.ID) %>%
  summarise(min_padj = min(padj, na.rm = TRUE), .groups = "drop") %>%
  inner_join(mean_rel_abund, by = "Feature.ID") %>%
  arrange(min_padj, desc(mean_rel_abund)) %>%
  slice_head(n = 30) %>%
  pull(Feature.ID)

lfc_wide_top30 <- combined_linda %>%
  filter(Feature.ID %in% top30_features) %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  mutate(
    Taxa_Label = case_when(
      !is.na(Genus) & Genus != "g__" & Genus != ""    ~ paste0(Genus, " (", str_sub(Feature.ID, 1, 6), ")"),
      !is.na(Family) & Family != "f__" & Family != "" ~ paste0("Fam. ", Family, " (", str_sub(Feature.ID, 1, 6), ")"),
      TRUE                                            ~ paste0("Phy. ", Phylum, " (", str_sub(Feature.ID, 1, 6), ")")
    )
  ) %>%
  select(Taxa_Label, Comparison, lfc) %>%
  pivot_wider(names_from = Comparison, values_from = lfc, values_fill = 0) %>%
  arrange(desc(`B vs. M(B)`))

lfc_matrix_top30 <- lfc_wide_top30 %>%
  column_to_rownames("Taxa_Label") %>%
  as.matrix()

lfc_matrix_top30 <- lfc_matrix_top30[, c("B vs. M(B)", "J vs. M(J)")]
lfc_matrix_top30[lfc_matrix_top30 > 8]  <- 8
lfc_matrix_top30[lfc_matrix_top30 < -8] <- -8

pheatmap(
  mat           = lfc_matrix_top30,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color         = colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100),
  breaks        = seq(-8, 8, length.out = 101),
  border_color  = "gray80",
  fontsize_row  = 9,
  fontsize_col  = 11,
  angle_col     = "0",
  fontfamily    = font_family,
  main          = "10. B vs. M(B), J vs. M(J) - Top 30 Heatmap (LinDA)",
  filename      = "assets/044_16S_Top30_Heatmap.png",
  width         = 8.5,
  height        = 8.5
)

# ==============================================================================
# 5. Group-Exclusive Rare Taxa Mining
# ==============================================================================
message("--- Step 5: Mining Group-Exclusive Rare Biosphere Taxa ---")

meta_rare <- meta_final %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")) %>%
  mutate(
    Group_Code = case_when(
      Treatment_Group == "GB_Control"   ~ "M(B)",
      Treatment_Group == "GB_Treatment" ~ "B",
      Treatment_Group == "GJ_Control"   ~ "M(J)",
      Treatment_Group == "GJ_Treatment" ~ "J"
    )
  )

counts_sub_rare <- counts_raw[, meta_rare$`sample-id`]
rel_abund_rare <- apply(counts_sub_rare, 2, function(x) (x / sum(x)) * 100)

df_long_rare <- rel_abund_rare %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "RelAbund") %>%
  inner_join(meta_rare %>% select(`sample-id`, Group_Code), by = "sample-id")

group_stats <- df_long_rare %>%
  group_by(Feature.ID, Group_Code) %>%
  summarise(
    Sample_Count      = n(),
    Present_Count     = sum(RelAbund > 0),
    Prevalence_Pct    = (sum(RelAbund > 0) / n()) * 100,
    Mean_RelAbund_Pct = mean(RelAbund),
    .groups           = "drop"
  )

prev_wide <- group_stats %>%
  select(Feature.ID, Group_Code, Prevalence_Pct) %>%
  pivot_wider(names_from = Group_Code, values_from = Prevalence_Pct, values_fill = 0)

sample_counts_wide <- group_stats %>%
  select(Feature.ID, Group_Code, Present_Count) %>%
  pivot_wider(names_from = Group_Code, values_from = Present_Count, values_fill = 0)

abund_wide <- group_stats %>%
  select(Feature.ID, Group_Code, Mean_RelAbund_Pct) %>%
  pivot_wider(names_from = Group_Code, values_from = Mean_RelAbund_Pct, values_fill = 0, names_prefix = "MeanAbund_")

min_samples_detected <- 2

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
      Exclusive_Group == "B"    ~ `B_count`,
      Exclusive_Group == "M(J)" ~ `M(J)_count`,
      Exclusive_Group == "J"    ~ `J_count`
    ),
    Prevalence_Pct = case_when(
      Exclusive_Group == "M(B)" ~ `M(B)_prev`,
      Exclusive_Group == "B"    ~ `B_prev`,
      Exclusive_Group == "M(J)" ~ `M(J)_prev`,
      Exclusive_Group == "J"    ~ `J_prev`
    ),
    Mean_RelAbund_Pct = case_when(
      Exclusive_Group == "M(B)" ~ `MeanAbund_M(B)`,
      Exclusive_Group == "B"    ~ `MeanAbund_B`,
      Exclusive_Group == "M(J)" ~ `MeanAbund_M(J)`,
      Exclusive_Group == "J"    ~ `MeanAbund_J`
    )
  ) %>%
  select(
    Feature.ID, Exclusive_Group, Detected_Sample_Count, Prevalence_Pct, Mean_RelAbund_Pct,
    Phylum, Class, Order, Family, Genus
  ) %>%
  arrange(Exclusive_Group, desc(Mean_RelAbund_Pct))

out_tsv <- "Bacteria_16S_Exclusive_Rare_Taxa_Table.tsv"
write_tsv(exclusive_taxa_df, out_tsv)
cat("[+] Saved:", out_tsv, "\n")
print(table(exclusive_taxa_df$Exclusive_Group))

message("=== [", format(Sys.time(), "%T"), "] LinDA Differential Testing & Rare Taxa Analysis Complete ===")