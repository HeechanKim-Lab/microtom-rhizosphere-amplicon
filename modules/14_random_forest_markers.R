#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 14_random_forest_markers.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1-3, GJ1-2; N=110 Locked Cohort)
# Platform: R (Native ARM64 / macOS)
# Description: Random Forest on the rarefied 16S table (GB_Treatment vs
#              GJ_Treatment; all four groups), top 30 ASVs by permutation
#              importance, selection of GB- and GJ-specific marker ASVs with the
#              LinDA results of module 13, and marker-only copies of two
#              PICRUSt2 tables for module 15.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Dependencies, Global Typography & Directory Setup
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(randomForest)
  library(pheatmap)
  library(RColorBrewer)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)
if (!dir.exists("tables")) dir.create("tables", recursive = TRUE)
out_dir <- "rf_picrust2"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing Random Forest and Marker Selection Module ===")

# Verify core inputs
meta_file  <- "metadata_P_generation.tsv"
table_file <- "07_diversity_16S/core_metrics/rarefied_table.qza"  # rarefied, 27,000 reads per sample
tax_file   <- "05_taxonomy_16S/taxonomy_16S.qza"
linda_gb_file <- "tables/LinDA_ASV_GB_Treatment_vs_GB_Control.csv"  # module 13
linda_gj_file <- "tables/LinDA_ASV_GJ_Treatment_vs_GJ_Control.csv"  # module 13
seqtab_file   <- "picrust2_out_16S/EC_metagenome_out/seqtab_norm.tsv.gz"
ec_pred_file  <- "picrust2_out_16S/combined_EC_predicted.tsv.gz"

for (f in c(meta_file, table_file, tax_file, linda_gb_file, linda_gj_file, seqtab_file, ec_pred_file)) {
  if (!file.exists(f)) stop("[-] ERROR: Required file '", f, "' not found.")
}

# Rarefied counts -> relative abundance (%) per sample; taxonomy
counts_raw <- read_qza(table_file)$data
tax_raw    <- read_qza(tax_file)$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

rel_abund_mat <- sweep(counts_raw, 2, colSums(counts_raw), FUN = "/") * 100

heatmap_colors <- colorRampPalette(c("#1F618D", "#FFFFFF", "#B03A2E"))(100)

# Taxon label used in the figures and tables: genus, else family, order or phylum
add_display_name <- function(df) {
  df %>%
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
    )
}

# ==============================================================================
# 1. Random Forest: GB_Treatment vs. GJ_Treatment (55 samples; Figures 089, 090)
# ==============================================================================
message("--- Step 1: Random Forest, GB_Treatment vs. GJ_Treatment ---")

palette_rf <- c(
  "GB_Treatment" = "#1A365D",
  "GJ_Treatment" = "#881337"
)

meta_rf <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Treatment", "GJ_Treatment")
  ) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Treatment", "GJ_Treatment"))) %>%
  column_to_rownames("sample-id")

common_samples <- intersect(rownames(meta_rf), colnames(rel_abund_mat))
meta_rf <- meta_rf[common_samples, , drop = FALSE]
rel_rf  <- rel_abund_mat[, common_samples]

# 1.1 Unsupervised pre-filter: present in >= 15% of samples (>= 9 samples) and mean relative abundance >= 0.005%
prev_threshold <- 0.15 * ncol(rel_rf)
keep_features <- rownames(rel_rf)[
  rowSums(rel_rf > 0) >= prev_threshold & rowMeans(rel_rf) >= 0.005
]
message(sprintf("    Retained %d ASVs out of %d.", length(keep_features), nrow(rel_rf)))

rf_input <- as.data.frame(t(rel_rf[keep_features, ]))
x_data <- rf_input
y_data <- meta_rf$Treatment_Group

# 1.2 Model training
set.seed(42)
best_mtry <- floor(sqrt(ncol(x_data)))  # 47 for 2,261 ASVs

rf_model <- randomForest(
  x = x_data,
  y = y_data,
  ntree = 2000,
  mtry = best_mtry,
  importance = TRUE,
  proximity = TRUE,
  sampsize = c(22, 22),  # 22 samples per class in each bootstrap sample
  nodesize = 2
)
print(rf_model)

# 1.3 Permutation importance (MDA) and taxonomy; top 30 ASVs
imp_df <- as.data.frame(importance(rf_model)) %>%
  rownames_to_column("Feature.ID") %>%
  add_display_name() %>%
  arrange(desc(MeanDecreaseAccuracy))

top30_rf <- imp_df %>%
  slice_head(n = 30)

# Enriched_In = group with the higher mean relative abundance
group_means <- rf_input %>%
  mutate(Target_Class = y_data) %>%
  group_by(Target_Class) %>%
  summarise(across(all_of(top30_rf$Feature.ID), mean)) %>%
  pivot_longer(-Target_Class, names_to = "Feature.ID", values_to = "Mean_Abund")

enrich_dir <- group_means %>%
  group_by(Feature.ID) %>%
  slice_max(Mean_Abund, n = 1) %>%
  select(Feature.ID, Enriched_In = Target_Class)

top30_rf <- top30_rf %>%
  left_join(enrich_dir, by = "Feature.ID") %>%
  mutate(DisplayName = factor(DisplayName, levels = rev(DisplayName)))

write_csv(top30_rf, "tables/RF_16S_Top30_Diagnostic_Biomarkers.csv")

# 1.4 Figure 089: variable importance lollipop chart
p_089 <- ggplot(top30_rf, aes(x = MeanDecreaseAccuracy, y = DisplayName)) +
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

ggsave("assets/089_16S_RandomForest_Treatment.png", p_089, width = 1600 / 192, height = 1200 / 192, dpi = 192)

# 1.5 Figure 090: z-score heatmap of the top 20 ASVs
top20_ids <- top30_rf$Feature.ID[1:20]
rf_heatmap_mat <- rel_rf[top20_ids, ]
rownames(rf_heatmap_mat) <- top30_rf$DisplayName[match(rownames(rf_heatmap_mat), top30_rf$Feature.ID)]

col_annot <- meta_rf %>%
  select(Treatment_Group) %>%
  arrange(Treatment_Group)

rf_heatmap_mat <- rf_heatmap_mat[, rownames(col_annot)]

p_090 <- pheatmap(
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
  main = "6. Relative Abundance of Top Diagnostic RF Biomarkers (Z-Score)",
  silent = TRUE
)

png("assets/090_16S_RandomForest_Treatment_Z.png", width = 2400, height = 1400, res = 192)
grid::grid.newpage()
grid::grid.draw(p_090$gtable)
invisible(dev.off())

# ==============================================================================
# 2. Random Forest: four groups (110 samples; Figures 091, 092)
# ==============================================================================
message("--- Step 2: Random Forest, four groups ---")

palette_4 <- c(
  "GB_Control"   = "#4A90E2",
  "GB_Treatment" = "#1A365D",
  "GJ_Control"   = "#F87171",
  "GJ_Treatment" = "#881337"
)

meta_4g <- read_tsv(meta_file, show_col_types = FALSE) %>%
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

common_samples <- intersect(rownames(meta_4g), colnames(rel_abund_mat))
meta_4g <- meta_4g[common_samples, , drop = FALSE]
rel_4g  <- rel_abund_mat[, common_samples]

# 2.1 Unsupervised pre-filter: present in >= 10% of samples (>= 11 samples) and mean relative abundance >= 0.005%
prev_threshold <- 0.10 * ncol(rel_4g)
keep_features_4g <- rownames(rel_4g)[
  rowSums(rel_4g > 0) >= prev_threshold & rowMeans(rel_4g) >= 0.005
]
message(sprintf("    Retained %d ASVs out of %d.", length(keep_features_4g), nrow(rel_4g)))

x_4g <- as.data.frame(t(rel_4g[keep_features_4g, ]))
y_4g <- meta_4g$Treatment_Group

# 2.2 Model training
set.seed(42)
best_mtry_4g <- floor(sqrt(ncol(x_4g)))

rf_4g_model <- randomForest(
  x = x_4g,
  y = y_4g,
  ntree = 2000,
  mtry = best_mtry_4g,
  importance = TRUE,
  proximity = TRUE,
  sampsize = c(22, 22, 22, 22),  # 22 samples per class in each bootstrap sample
  nodesize = 2
)
print(rf_4g_model)

# 2.3 Permutation importance and taxonomy; top 30 ASVs
imp_df_4g <- as.data.frame(importance(rf_4g_model)) %>%
  rownames_to_column("Feature.ID") %>%
  add_display_name() %>%
  arrange(desc(MeanDecreaseAccuracy))

top30_rf_4g <- imp_df_4g %>%
  slice_head(n = 30)

# Enriched_In = group with the highest mean relative abundance
group_means_4g <- x_4g %>%
  mutate(Target_Class = y_4g) %>%
  group_by(Target_Class) %>%
  summarise(across(all_of(top30_rf_4g$Feature.ID), mean)) %>%
  pivot_longer(-Target_Class, names_to = "Feature.ID", values_to = "Mean_Abund")

enrich_dir_4g <- group_means_4g %>%
  group_by(Feature.ID) %>%
  slice_max(Mean_Abund, n = 1) %>%
  select(Feature.ID, Enriched_In = Target_Class)

top30_rf_4g <- top30_rf_4g %>%
  left_join(enrich_dir_4g, by = "Feature.ID") %>%
  mutate(DisplayName = factor(DisplayName, levels = rev(DisplayName)))

write_csv(top30_rf_4g, "tables/RF_16S_4Group_Top30_Biomarkers.csv")

# 2.4 Figure 091: variable importance lollipop chart
p_091 <- ggplot(top30_rf_4g, aes(x = MeanDecreaseAccuracy, y = DisplayName)) +
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

ggsave("assets/091_16S_RandomForest_4Groups.png", p_091, width = 1600 / 192, height = 1200 / 192, dpi = 192)

# 2.5 Figure 092: z-score heatmap of the top 20 ASVs, samples ordered by group
top20_ids_4g <- top30_rf_4g$Feature.ID[1:20]
rf_heat_mat <- rel_4g[top20_ids_4g, ]
rownames(rf_heat_mat) <- top30_rf_4g$DisplayName[match(rownames(rf_heat_mat), top30_rf_4g$Feature.ID)]

col_annot_4g <- meta_4g %>%
  select(Treatment_Group) %>%
  arrange(Treatment_Group)

rf_heat_mat <- rf_heat_mat[, rownames(col_annot_4g)]

group_counts <- table(col_annot_4g$Treatment_Group)
gaps_col <- cumsum(group_counts)[1:(length(group_counts) - 1)]

p_092 <- pheatmap(
  mat = rf_heat_mat,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = col_annot_4g,
  annotation_colors = list(Treatment_Group = palette_4),
  gaps_col = gaps_col,
  color = heatmap_colors,
  show_colnames = FALSE,
  fontsize_row = 9.5,
  fontfamily = font_family,
  border_color = NA,
  main = "7. Relative Abundance of Top 20 Multiclass RF Biomarkers (N = 110, Z-Score)",
  silent = TRUE
)

png("assets/092_16S_RandomForest_4Groups_Z.png", width = 2400, height = 1400, res = 192)
grid::grid.newpage()
grid::grid.draw(p_092$gtable)
invisible(dev.off())

# ==============================================================================
# 3. GB- and GJ-specific marker ASVs (RF top 30 of Step 1, checked against LinDA)
# ==============================================================================
message("--- Step 3: Marker selection with LinDA (module 13) ---")

rf_top30 <- read_csv("tables/RF_16S_Top30_Diagnostic_Biomarkers.csv", show_col_types = FALSE)

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
  read_linda_calls(linda_gb_file, "LinDA_GB"),
  read_linda_calls(linda_gj_file, "LinDA_GJ")
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
# 4. Marker-only copies of the two per-ASV PICRUSt2 tables (input of module 15)
# ==============================================================================
message("--- Step 4: Marker-only PICRUSt2 inputs ---")

marker_ids <- markers %>% filter(Marker_Set != "excluded") %>% pull(Feature.ID)

# Rows are copied unchanged
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

subset_rows(seqtab_file, file.path(out_dir, "seqtab_norm_markers.tsv.gz"), marker_ids)
subset_rows(ec_pred_file, file.path(out_dir, "EC_predicted_markers.tsv.gz"), marker_ids)

message("=== [", format(Sys.time(), "%T"), "] Random Forest & Marker Selection Complete ===")
