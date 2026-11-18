#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 16_marker_pathway_comparison.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1-3, GJ1-2; N=110 Locked Cohort)
# Platform: R (Native ARM64 / macOS)
# Description: Compares the predicted MetaCyc pathways (module 15) of the
#              GB-specific and GJ-specific RF marker ASVs (module 14), pathway
#              by pathway (Fisher's exact test, BH correction), adds the
#              community-level pathway LinDA results of module 13, and draws
#              the share plot (093) and presence heatmap (094).
#              PICRUSt2 predicts genomic potential, not activity.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Dependencies, Global Typography & Directory Setup
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(RColorBrewer)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)
out_dir <- "rf_picrust2"
pwy_dir <- file.path(out_dir, "pathways_markers")

message("=== [", format(Sys.time(), "%T"), "] Initializing Marker Pathway Comparison Module ===")

# Verify core inputs
pwy_linda_gb_file <- "tables/PICRUSt2_LinDA_GB_Control_vs_GB_Treatment.csv"  # module 13
pwy_linda_gj_file <- "tables/PICRUSt2_LinDA_GJ_Control_vs_GJ_Treatment.csv"  # module 13

for (f in c(file.path(out_dir, "rf_markers.tsv"),                           # module 14
            file.path(pwy_dir, "path_abun_predictions.tsv.gz"),             # module 15
            file.path(pwy_dir, "path_abun_unstrat_per_seq_descrip.tsv.gz"), # module 15
            pwy_linda_gb_file, pwy_linda_gj_file)) {
  if (!file.exists(f)) stop("[-] ERROR: Required file '", f, "' not found.")
}

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

# MetaCyc descriptions added by add_descriptions.py in module 15 (HTML entities such as &beta; -> beta)
pwy_desc <- read_tsv(file.path(pwy_dir, "path_abun_unstrat_per_seq_descrip.tsv.gz"), show_col_types = FALSE) %>%
  transmute(pathway, description = str_replace_all(description, "&([A-Za-z]+);", "\\1"))

# Community-level result for each pathway (pathway LinDA, treated vs its own control)
pwy_linda <- bind_rows(
  read_csv(pwy_linda_gb_file, show_col_types = FALSE) %>% mutate(Soil = "GB"),
  read_csv(pwy_linda_gj_file, show_col_types = FALSE) %>% mutate(Soil = "GJ")
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

message("    Predicted pathways per marker ASV:")
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
# 3. Figure 093: share of GB- vs GJ-specific markers carrying each pathway
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

ggsave("assets/093_16S_RF_Markers_Pathway_Share.png", p_share,
       width = 1400 / 192, height = 1400 / 192, dpi = 192)

# ==============================================================================
# 4. Figure 094: presence of the differing pathways in each marker ASV
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

  png("assets/094_16S_RF_Markers_Pathway_Heatmap.png", width = 2800, height = 1400, res = 192)
  grid::grid.newpage()
  grid::grid.draw(p_heat$gtable)
  invisible(dev.off())
}

message("=== [", format(Sys.time(), "%T"), "] Marker Pathway Comparison Complete ===")
