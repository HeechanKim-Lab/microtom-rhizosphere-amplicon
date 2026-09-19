#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 13_linda_and_picrust2_functional.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1-3, GJ1-2; N=110 Locked Cohort)
# Platform: R (Native ARM64 / macOS)
# Description: Executes LinDA differential abundance modeling on ASVs, generates
#              volcano plots, validates immigrant colonization vs. competitive
#              depletion via top-20 boxplots, renders global effect heatmaps,
#              and performs PICRUSt2 functional pathway enrichment modeling.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Dependencies, Global Typography & Directory Setup
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(ggrepel)
  library(ggpubr)
  library(MicrobiomeStat)
  library(pheatmap)
  library(ggpicrust2)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)
if (!dir.exists("tables")) dir.create("tables", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing LinDA and PICRUSt2 Functional Module ===")

# Verify core inputs
meta_file <- "metadata_P_generation.tsv"
table_file <- "07_diversity_16S/core_metrics/rarefied_table.qza"
tax_file   <- "05_taxonomy_16S/taxonomy_16S.qza"
path_file  <- "picrust2_out_16S/pathways_out/path_abun_unstrat.tsv.gz"

for (f in c(meta_file, table_file, tax_file)) {
  if (!file.exists(f)) stop("[-] ERROR: Required file '", f, "' not found.")
}

# ==============================================================================
# 0.1 Helper: MetaCyc Pathway Annotation Function
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

# ------------------------------------------------------------------------------
# 0.2 Ingest Metadata, Counts, and Taxonomy
# ------------------------------------------------------------------------------
metadata_all <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types", Domain == "Bacteria")

counts_raw <- read_qza(table_file)$data
tax_clean  <- parse_taxonomy(read_qza(tax_file)$data) %>%
  rownames_to_column("Feature.ID")

# ==============================================================================
# 1. Gijang B (GB) LinDA Differential Abundance & Boxplot Audits
# ==============================================================================
message("--- Step 1: LinDA Modeling & Diagnostics for GB Cohort ---")

meta_gb <- metadata_all %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment")))

common_gb <- intersect(meta_gb$`sample-id`, colnames(counts_raw))
counts_gb <- counts_raw[, common_gb]
meta_gb   <- meta_gb %>% filter(`sample-id` %in% common_gb)

rel_abund_gb <- apply(counts_gb, 2, function(x) (x / sum(x)) * 100)
mean_rel_abund_gb <- rowMeans(rel_abund_gb) %>%
  enframe(name = "Feature.ID", value = "mean_rel_abund")

set.seed(42)
linda_res_gb <- linda(
  feature.dat       = counts_gb,
  meta.dat          = meta_gb %>% column_to_rownames("sample-id"),
  formula           = "~ Treatment_Group",
  alpha             = 0.05,
  prev.filter       = 0.10,
  mean.abund.filter = 0.0001
)

df_gb <- linda_res_gb$output$Treatment_GroupGB_Treatment %>%
  rownames_to_column("Feature.ID") %>%
  mutate(lfc = log2FoldChange, neg_log10_padj = -log10(padj))

volcano_df_gb <- df_gb %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  inner_join(mean_rel_abund_gb, by = "Feature.ID") %>%
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

# Figure 078: GB Volcano Plot
p_078 <- ggplot(volcano_df_gb, aes(x = lfc, y = neg_log10_padj, color = Significance)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(aes(size = mean_rel_abund), alpha = 0.75) +
  scale_size_continuous(range = c(1.5, 5.5), name = "Mean Rel. Abund (%)") +
  scale_color_manual(values = c("Enriched in B" = "#41B6C4", "Enriched in M(B)" = "#2B5C8F", "Not Significant" = "#C0C0C0")) +
  scale_y_continuous(name = expression(bold("-Log"[10] * " (FDR-adjusted " * italic("p") * "-value)")), expand = expansion(mult = c(0, 0.05))) +
  geom_text_repel(
    data = top_labels_gb, aes(label = Label_Name, color = Significance),
    size = 3.5, fontface = "italic", family = font_family, max.overlaps = 15, show.legend = FALSE
  ) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "4. B vs. M(B) Volcano Plot (LinDA)", x = expression(bold("Log"[2] * " Fold Change (B / M(B))")))

ggsave("assets/078_16S_GB_Volcano_Plot.png", p_078, width = 8.5, height = 7.5, dpi = 300)

# Long relative abundance table for GB
rel_abund_gb_long <- counts_gb %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  inner_join(meta_gb, by = "sample-id")

palette_gb <- c("GB_Control" = "#2B5C8F", "GB_Treatment" = "#41B6C4")

# Figure 080: GB Top 20 Positive Shift Taxa (Immigrant Confirmation)
top20_shift_gb <- volcano_df_gb %>%
  filter(Significance == "Enriched in B") %>%
  slice_max(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(desc(lfc)) %>%
  mutate(
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Name, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name)
  )

p_080 <- rel_abund_gb_long %>%
  filter(Feature.ID %in% top20_shift_gb$Feature.ID) %>%
  inner_join(top20_shift_gb %>% select(Feature.ID, Facet_Label), by = "Feature.ID") %>%
  ggplot(aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.35) +
  geom_jitter(aes(color = Treatment_Group), width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gb) +
  scale_color_manual(values = palette_gb) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)), limits = c(0, NA)) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", size = 9),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold.italic", size = 8.5),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
  ) +
  labs(title = "Relative Abundance of Top 20 Large-Shift Taxa (GB Treatment)", y = "Relative Abundance (%)", x = NULL)

ggsave("assets/080_16S_GB_Top20_LargeShift_Taxa_Boxplot.png", p_080, width = 12, height = 7.5, dpi = 300)

# Figure 082: GB Top 20 Depleted Taxa (Competitive Dilution Confirmation)
top20_depleted_gb <- volcano_df_gb %>%
  filter(Significance == "Enriched in M(B)") %>%
  slice_min(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(lfc) %>%
  mutate(
    Label_Clean  = if_else(nchar(Label_Name) > 22, paste0(substr(Label_Name, 1, 20), ".."), Label_Name),
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Clean, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name)
  )

p_082 <- rel_abund_gb_long %>%
  filter(Feature.ID %in% top20_depleted_gb$Feature.ID) %>%
  inner_join(top20_depleted_gb %>% select(Feature.ID, Facet_Label), by = "Feature.ID") %>%
  ggplot(aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.35) +
  geom_jitter(aes(color = Treatment_Group), width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gb) +
  scale_color_manual(values = palette_gb) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)), limits = c(0, NA)) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", size = 9),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold.italic", size = 8.5),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
  ) +
  labs(title = "Relative Abundance of Top 20 Depleted Taxa (Enriched in GB Control)", y = "Relative Abundance (%)", x = NULL)

ggsave("assets/082_16S_GB_Top20_N_LargeShift_Taxa_Boxplot.png", p_082, width = 12, height = 7.5, dpi = 300)

# ==============================================================================
# 2. Gyeongju (GJ) LinDA Differential Abundance & Boxplot Audits
# ==============================================================================
message("--- Step 2: LinDA Modeling & Diagnostics for GJ Cohort ---")

meta_gj <- metadata_all %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GJ_Control", "GJ_Treatment")))

common_gj <- intersect(meta_gj$`sample-id`, colnames(counts_raw))
counts_gj <- counts_raw[, common_gj]
meta_gj   <- meta_gj %>% filter(`sample-id` %in% common_gj)

rel_abund_gj <- apply(counts_gj, 2, function(x) (x / sum(x)) * 100)
mean_rel_abund_gj <- rowMeans(rel_abund_gj) %>%
  enframe(name = "Feature.ID", value = "mean_rel_abund")

set.seed(42)
linda_res_gj <- linda(
  feature.dat       = counts_gj,
  meta.dat          = meta_gj %>% column_to_rownames("sample-id"),
  formula           = "~ Treatment_Group",
  alpha             = 0.05,
  prev.filter       = 0.10,
  mean.abund.filter = 0.0001
)

df_gj <- linda_res_gj$output$Treatment_GroupGJ_Treatment %>%
  rownames_to_column("Feature.ID") %>%
  mutate(lfc = log2FoldChange, neg_log10_padj = -log10(padj))

volcano_df_gj <- df_gj %>%
  inner_join(tax_clean, by = "Feature.ID") %>%
  inner_join(mean_rel_abund_gj, by = "Feature.ID") %>%
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

# Figure 079: GJ Volcano Plot
p_079 <- ggplot(volcano_df_gj, aes(x = lfc, y = neg_log10_padj, color = Significance)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_point(aes(size = mean_rel_abund), alpha = 0.75) +
  scale_size_continuous(range = c(1.5, 5.5), name = "Mean Rel. Abund (%)") +
  scale_color_manual(values = c("Enriched in J" = "#FE9929", "Enriched in M(J)" = "#D95F02", "Not Significant" = "#C0C0C0")) +
  scale_y_continuous(name = expression(bold("-Log"[10] * " (FDR-adjusted " * italic("p") * "-value)")), expand = expansion(mult = c(0, 0.05))) +
  geom_text_repel(
    data = top_labels_gj, aes(label = Label_Name, color = Significance),
    size = 3.5, fontface = "italic", family = font_family, max.overlaps = 15, show.legend = FALSE
  ) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "4. J vs. M(J) Volcano Plot (LinDA)", x = expression(bold("Log"[2] * " Fold Change (J / M(J))")))

ggsave("assets/079_16S_GJ_Volcano_Plot.png", p_079, width = 8.5, height = 7.5, dpi = 300)

# Long relative abundance table for GJ
rel_abund_gj_long <- counts_gj %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  inner_join(meta_gj, by = "sample-id")

palette_gj <- c("GJ_Control" = "#D95F02", "GJ_Treatment" = "#FE9929")

# Figure 081: GJ Top 20 Positive Shift Taxa (Immigrant Confirmation)
top20_shift_gj <- volcano_df_gj %>%
  filter(Significance == "Enriched in J") %>%
  slice_max(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(desc(lfc)) %>%
  mutate(
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Name, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name)
  )

p_081 <- rel_abund_gj_long %>%
  filter(Feature.ID %in% top20_shift_gj$Feature.ID) %>%
  inner_join(top20_shift_gj %>% select(Feature.ID, Facet_Label), by = "Feature.ID") %>%
  ggplot(aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.35) +
  geom_jitter(aes(color = Treatment_Group), width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gj) +
  scale_color_manual(values = palette_gj) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)), limits = c(0, NA)) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", size = 9),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold.italic", size = 8.5),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
  ) +
  labs(title = "Relative Abundance of Top 20 Large-Shift Taxa (GJ Treatment)", y = "Relative Abundance (%)", x = NULL)

ggsave("assets/081_16S_GJ_Top20_LargeShift_Taxa_Boxplot.png", p_081, width = 12, height = 7.5, dpi = 300)

# Figure 083: GJ Top 20 Depleted Taxa (Competitive Dilution Confirmation)
top20_depleted_gj <- volcano_df_gj %>%
  filter(Significance == "Enriched in M(J)") %>%
  slice_min(order_by = lfc, n = 20, with_ties = FALSE) %>%
  arrange(lfc) %>%
  mutate(
    Label_Clean  = if_else(nchar(Label_Name) > 22, paste0(substr(Label_Name, 1, 20), ".."), Label_Name),
    Short_ID     = substr(Feature.ID, 1, 5),
    Display_Name = paste0(Label_Clean, " (", Short_ID, ")"),
    Facet_Label  = factor(Display_Name, levels = Display_Name)
  )

p_083 <- rel_abund_gj_long %>%
  filter(Feature.ID %in% top20_depleted_gj$Feature.ID) %>%
  inner_join(top20_depleted_gj %>% select(Feature.ID, Facet_Label), by = "Feature.ID") %>%
  ggplot(aes(x = Treatment_Group, y = RelAbund, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.45, color = "black", linewidth = 0.35) +
  geom_jitter(aes(color = Treatment_Group), width = 0.15, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Facet_Label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = palette_gj) +
  scale_color_manual(values = palette_gj) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18)), limits = c(0, NA)) +
  theme_bw(base_size = 11, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold", size = 9),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold.italic", size = 8.5),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
  ) +
  labs(title = "Relative Abundance of Top 20 Depleted Taxa (Enriched in GJ Control)", y = "Relative Abundance (%)", x = NULL)

ggsave("assets/083_16S_GJ_Top20_N_LargeShift_Taxa_Boxplot.png", p_083, width = 12, height = 7.5, dpi = 300)

# ==============================================================================
# 3. Global Effect Heatmap (Figure 084)
# ==============================================================================
message("--- Step 3: Generating Global Differential Heatmap ---")

combined_linda <- bind_rows(
  df_gb %>% transmute(Feature.ID, Comparison = "GB_Treatment vs. GB_Control", lfc, padj),
  df_gj %>% transmute(Feature.ID, Comparison = "GJ_Treatment vs. GJ_Control", lfc, padj)
)

sig_features <- combined_linda %>%
  filter(padj < 0.05, abs(lfc) >= 1.0) %>%
  pull(Feature.ID) %>%
  unique()

lfc_mat_global <- combined_linda %>%
  filter(Feature.ID %in% sig_features) %>%
  select(Feature.ID, Comparison, lfc) %>%
  pivot_wider(names_from = Comparison, values_from = lfc, values_fill = 0) %>%
  column_to_rownames("Feature.ID") %>%
  as.matrix()

lfc_mat_global <- lfc_mat_global[, c("GB_Treatment vs. GB_Control", "GJ_Treatment vs. GJ_Control")]
lfc_mat_global[lfc_mat_global > 8]  <- 8
lfc_mat_global[lfc_mat_global < -8] <- -8

pheatmap(
  mat                      = lfc_mat_global,
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
  main                     = "5. Heatmap (LinDA)",
  filename                 = "assets/084_16S_Heatmap.png",
  width                    = 8.5,
  height                   = 8.5
)

# ==============================================================================
# 4. PICRUSt2 Functional Pathway Enrichment (Figures 085, 086)
# ==============================================================================
message("--- Step 4: PICRUSt2 Pathway Differential Modeling & Forest Plots ---")

if (file.exists(path_file)) {
  path_raw <- read_delim(path_file, delim = "\t", comment = "#", show_col_types = FALSE) %>%
    rename(pathway = 1) %>%
    column_to_rownames("pathway") %>%
    as.matrix()

  # Common annotator lookup
  annot_lookup <- get_metacyc_annotations(rownames(path_raw))

  # 4.1 GB Pathway Differential Modeling
  common_pwy_gb <- intersect(meta_gb$`sample-id`, colnames(path_raw))
  path_gb <- path_raw[, common_pwy_gb]
  meta_pwy_gb <- meta_gb %>% filter(`sample-id` %in% common_pwy_gb) %>% column_to_rownames("sample-id")

  set.seed(42)
  linda_pwy_gb <- linda(
    feature.dat   = path_gb,
    meta.dat      = meta_pwy_gb,
    formula       = "~ Treatment_Group",
    alpha         = 0.05,
    prev.filter   = 0.10,
    zero.handling = "pseudo-count"
  )

  df_pwy_gb <- linda_pwy_gb$output[[1]] %>%
    rownames_to_column("feature") %>%
    rename(log2FC = log2FoldChange, se = lfcSE) %>%
    left_join(annot_lookup, by = "feature") %>%
    mutate(
      description = coalesce(description, feature),
      description = case_when(feature == "PWY-8190" ~ "L-glutamate degradation XI (reductive Stickland reaction)", TRUE ~ description),
      neg_log10_padj = -log10(padj),
      ci_lower = log2FC - 1.96 * se,
      ci_upper = log2FC + 1.96 * se,
      Enrichment = ifelse(log2FC > 0, "GB_Treatment", "GB_Control")
    )

  write_csv(df_pwy_gb, "tables/PICRUSt2_LinDA_GB_Control_vs_GB_Treatment.csv")

  top_sig_gb <- df_pwy_gb %>% filter(padj < 0.05) %>% slice_min(padj, n = 20)

  p_085 <- ggplot(top_sig_gb, aes(x = log2FC, y = reorder(description, log2FC))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = Enrichment), height = 0.25, linewidth = 0.7) +
    geom_point(aes(size = neg_log10_padj, fill = Enrichment), shape = 21, color = "black", stroke = 0.5) +
    scale_color_manual(values = c("GB_Treatment" = "#7BCCC4", "GB_Control" = "#2B8CBE")) +
    scale_fill_manual(values = c("GB_Treatment" = "#7BCCC4", "GB_Control" = "#2B8CBE")) +
    scale_size_continuous(name = expression(-log[10](FDR)), range = c(2.5, 6.5)) +
    theme_bw(base_size = 11, base_family = font_family) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Differentially Abundant Pathways (GB_Treatment vs GB_Control)", x = expression(Log[2]~Fold~Change), y = "MetaCyc Pathway Description")

  ggsave("assets/085_16S_GB_PICRUSt2_Forestplot.png", p_085, width = 9.5, height = 7.5, dpi = 300)

  # 4.2 GJ Pathway Differential Modeling
  common_pwy_gj <- intersect(meta_gj$`sample-id`, colnames(path_raw))
  path_gj <- path_raw[, common_pwy_gj]
  meta_pwy_gj <- meta_gj %>% filter(`sample-id` %in% common_pwy_gj) %>% column_to_rownames("sample-id")

  set.seed(42)
  linda_pwy_gj <- linda(
    feature.dat   = path_gj,
    meta.dat      = meta_pwy_gj,
    formula       = "~ Treatment_Group",
    alpha         = 0.05,
    prev.filter   = 0.10,
    zero.handling = "pseudo-count"
  )

  df_pwy_gj <- linda_pwy_gj$output[[1]] %>%
    rownames_to_column("feature") %>%
    rename(log2FC = log2FoldChange, se = lfcSE) %>%
    left_join(annot_lookup, by = "feature") %>%
    mutate(
      description = coalesce(description, feature),
      description = case_when(feature == "PWY-8154" ~ "superpathway of cholesterol degradation III (oxidase)", TRUE ~ description),
      neg_log10_padj = -log10(padj),
      ci_lower = log2FC - 1.96 * se,
      ci_upper = log2FC + 1.96 * se,
      Enrichment = ifelse(log2FC > 0, "GJ_Treatment", "GJ_Control")
    )

  write_csv(df_pwy_gj, "tables/PICRUSt2_LinDA_GJ_Control_vs_GJ_Treatment.csv")

  top_sig_gj <- df_pwy_gj %>% filter(padj < 0.05) %>% slice_min(padj, n = 20)

  p_086 <- ggplot(top_sig_gj, aes(x = log2FC, y = reorder(description, log2FC))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = Enrichment), height = 0.25, linewidth = 0.7) +
    geom_point(aes(size = neg_log10_padj, fill = Enrichment), shape = 21, color = "black", stroke = 0.5) +
    scale_color_manual(values = c("GJ_Treatment" = "#FDAE6B", "GJ_Control" = "#E6550D")) +
    scale_fill_manual(values = c("GJ_Treatment" = "#FDAE6B", "GJ_Control" = "#E6550D")) +
    scale_size_continuous(name = expression(-log[10](FDR)), range = c(2.5, 6.5)) +
    theme_bw(base_size = 11, base_family = font_family) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Differentially Abundant Pathways (GJ_Treatment vs GJ_Control)", x = expression(Log[2]~Fold~Change), y = "MetaCyc Pathway Description")

  ggsave("assets/086_16S_GJ_PICRUSt2_Forestplot.png", p_086, width = 9.5, height = 7.5, dpi = 300)

  # ==============================================================================
  # 5. Sample-Level PICRUSt2 Z-Score Heatmap (Figure 087)
  # ==============================================================================
  message("--- Step 5: Generating Sample-Level PICRUSt2 Z-Score Heatmap ---")

  meta_ordered <- metadata_all %>%
    filter(Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")) %>%
    mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment"))) %>%
    arrange(Treatment_Group) %>%
    column_to_rownames("sample-id")

  common_samples_all <- intersect(rownames(meta_ordered), colnames(path_raw))
  path_rel <- sweep(path_raw[, common_samples_all], 2, colSums(path_raw[, common_samples_all]), FUN = "/")

  # Top 35 variable pathways
  top35_pwy <- names(sort(apply(path_rel, 1, sd), decreasing = TRUE))[1:35]
  path_sub  <- path_rel[top35_pwy, ]

  annot_sub <- get_metacyc_annotations(rownames(path_sub))
  mapped_sub <- setNames(annot_sub$description, annot_sub$feature)

  row_labs <- ifelse(!is.na(mapped_sub[rownames(path_sub)]), mapped_sub[rownames(path_sub)], rownames(path_sub)) %>%
    stringr::str_replace_all("&beta;", "beta") %>%
    stringr::str_replace_all("&alpha;", "alpha") %>%
    stringr::str_trunc(width = 75, side = "right")

  rownames(path_sub) <- make.unique(row_labs)

  annotation_col <- meta_ordered[common_samples_all, "Treatment_Group", drop = FALSE]
  annotation_colors <- list(
    Treatment_Group = c(
      "GB_Control"   = "#2B8CBE",
      "GB_Treatment" = "#7BCCC4",
      "GJ_Control"   = "#E6550D",
      "GJ_Treatment" = "#FDAE6B"
    )
  )

  group_counts <- table(annotation_col$Treatment_Group)
  gaps_col <- cumsum(group_counts)[1:(length(group_counts) - 1)]

  pheatmap(
    mat               = path_sub,
    scale             = "row",
    color             = colorRampPalette(c("#1F618D", "#FFFFFF", "#B03A2E"))(100),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    annotation_col    = annotation_col,
    annotation_colors = annotation_colors,
    gaps_col          = gaps_col,
    show_colnames     = FALSE,
    show_rownames     = TRUE,
    fontsize          = 9,
    fontsize_row      = 8,
    fontfamily        = font_family,
    main              = "Sample-level Predicted Pathway Profiles (Z-score)",
    filename          = "assets/087_16S_PICRUSt2_Heatmap(zscore).png",
    width             = 11.5,
    height            = 8.5
  )
} else {
  message("[!] Notice: PICRUSt2 abundance output '", path_file, "' not found. Skipping Steps 4 & 5.")
}

message("=== [", format(Sys.time(), "%T"), "] LinDA Modeling & PICRUSt2 Enrichment Complete ===")