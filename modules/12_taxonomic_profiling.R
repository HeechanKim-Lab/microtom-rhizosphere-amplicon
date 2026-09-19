#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 12_taxonomic_profiling.R
# Target: 16S rRNA (V3-V4) Rhizosphere Amplicon Sequencing
# Cohort: Verified 5-Batch Design (GB1-3, GJ1-2; N=110 Locked Cohort)
# Platform: R (Native ARM64 / macOS)
# Description: Generates multi-level taxonomic composition bar plots across
#              Phylum, Class, Order, Family, and Genus ranks (Figures 067-077),
#              enforcing exact two-step 100% group scaling and a shared,
#              harmonized color palette for regional contrasts.
# ==============================================================================

suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(ggpubr)
  library(RColorBrewer)
})

font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing Taxonomic Profiling Module ===")

# ==============================================================================
# 0. Load Metadata & Decontaminated Feature Tables
# ==============================================================================
meta_file <- "metadata_P_generation.tsv"
if (!file.exists(meta_file)) {
  stop("[-] ERROR: Required metadata '", meta_file, "' not found.")
}

metadata_all <- read_tsv(meta_file, show_col_types = FALSE) %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(
    Domain == "Bacteria",
    Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")
  ) %>%
  mutate(
    Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment"))
  )

counts_raw <- read_qza("07_diversity_16S/core_metrics/rarefied_table.qza")$data
tax_raw    <- read_qza("05_taxonomy_16S/taxonomy_16S.qza")$data

tax_clean <- parse_taxonomy(tax_raw) %>%
  rownames_to_column("Feature.ID")

common_samples_all <- intersect(metadata_all$`sample-id`, colnames(counts_raw))
counts_sub_all     <- counts_raw[, common_samples_all]

# Transform counts to sample-level relative abundance (%)
rel_abund_long <- counts_sub_all %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  pivot_longer(-Feature.ID, names_to = "sample-id", values_to = "Count") %>%
  group_by(`sample-id`) %>%
  mutate(RelAbund = (Count / sum(Count)) * 100) %>%
  ungroup() %>%
  left_join(tax_clean, by = "Feature.ID") %>%
  inner_join(metadata_all, by = "sample-id")

# Helper function for global 4-group summary and plotting
plot_global_rank <- function(df_long, rank_col, top_n, pal, title_str, out_png, is_italic = FALSE) {
  rank_sym <- sym(rank_col)
  prefix   <- paste0(tolower(substr(rank_col, 1, 1)), "__")
  
  sample_df <- df_long %>%
    mutate(!!rank_sym := if_else(is.na(!!rank_sym) | !!rank_sym == "" | !!rank_sym == prefix, "Unassigned", !!rank_sym)) %>%
    group_by(`sample-id`, Treatment_Group, !!rank_sym) %>%
    summarise(RelAbund = sum(RelAbund), .groups = "drop")
  
  top_taxa <- sample_df %>%
    group_by(!!rank_sym) %>%
    summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
    filter(!!rank_sym != "Unassigned") %>%
    arrange(desc(MeanAbund)) %>%
    slice_head(n = top_n) %>%
    pull(!!rank_sym)
  
  group_summary <- sample_df %>%
    mutate(Taxa_Group = if_else(!!rank_sym %in% top_taxa, !!rank_sym, "Other")) %>%
    group_by(`sample-id`, Treatment_Group, Taxa_Group) %>%
    summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
    group_by(Treatment_Group, Taxa_Group) %>%
    summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
    mutate(Taxa_Group = factor(Taxa_Group, levels = rev(c(top_taxa, "Other"))))
  
  taxa_colors <- setNames(c(pal, "#B0B0B0"), c(top_taxa, "Other"))
  
  p <- ggplot(group_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Taxa_Group)) +
    geom_bar(stat = "identity", position = "stack", width = 0.65, color = "black", linewidth = 0.25) +
    scale_fill_manual(values = taxa_colors) +
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
      legend.text = element_text(size = 9, family = font_family, face = if (is_italic) "italic" else "plain"),
      axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
      axis.text.y = element_text(color = "black", family = font_family),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
      panel.grid = element_blank()
    ) +
    labs(title = title_str, y = "Mean Relative Abundance (%)")
  
  ggsave(out_png, p, width = 8.5, height = 7.5, dpi = 300)
  return(p)
}

# ==============================================================================
# 1. Phylum-Level Composition (Figures 067, 068, 069)
# ==============================================================================
message("--- Step 1: Processing Phylum-Level Relative Abundance ---")

# Figure 067: Global Phylum
p_phylum <- plot_global_rank(
  df_long   = rel_abund_long,
  rank_col  = "Phylum",
  top_n     = 10,
  pal       = brewer.pal(10, "Paired"),
  title_str = "3. Relative Abundance (Phylum Level)",
  out_png   = "assets/067_16S_Relative_Abundance(Phylum).png"
)

# Shared Palette Preparation for Regional GB and GJ Contrasts
phylum_sample <- rel_abund_long %>%
  mutate(Phylum = if_else(is.na(Phylum) | Phylum == "" | Phylum == "p__", "Unassigned", Phylum)) %>%
  group_by(`sample-id`, Treatment_Group, Phylum) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

top10_gb_phy <- phylum_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Phylum != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

top10_gj_phy <- phylum_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Phylum) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Phylum != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

all_top_phyla <- unique(c(top10_gb_phy, top10_gj_phy))
palette_pool_phy <- colorRampPalette(brewer.pal(12, "Paired"))(length(all_top_phyla))
master_taxa_colors_phy <- c(setNames(palette_pool_phy, all_top_phyla), "Other" = "#B0B0B0")

# Figure 068: GB Lineage Phylum
gb_phylum_summary <- phylum_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  mutate(Phylum_Group = if_else(Phylum %in% top10_gb_phy, Phylum, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Phylum_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Phylum_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Phylum_Group = factor(Phylum_Group, levels = rev(c(top10_gb_phy, "Other"))))

p_phylum_gb <- ggplot(gb_phylum_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Phylum_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors_phy) +
  scale_x_discrete(labels = c("GB_Control" = "GB_Control", "GB_Treatment" = "GB_Treatment")) +
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
  labs(title = "3. GB Relative Abundance (Phylum Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/068_16S_GB_Relative_Abundance(Phylum).png", p_phylum_gb, width = 7.5, height = 7.5, dpi = 300)

# Figure 069: GJ Lineage Phylum
gj_phylum_summary <- phylum_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  mutate(Phylum_Group = if_else(Phylum %in% top10_gj_phy, Phylum, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Phylum_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Phylum_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Phylum_Group = factor(Phylum_Group, levels = rev(c(top10_gj_phy, "Other"))))

p_phylum_gj <- ggplot(gj_phylum_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Phylum_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors_phy) +
  scale_x_discrete(labels = c("GJ_Control" = "GJ_Control", "GJ_Treatment" = "GJ_Treatment")) +
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
  labs(title = "3. GJ Relative Abundance (Phylum Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/069_16S_GJ_Relative_Abundance(Phylum).png", p_phylum_gj, width = 7.5, height = 7.5, dpi = 300)

# ==============================================================================
# 2. Class-Level Composition (Figure 070)
# ==============================================================================
message("--- Step 2: Processing Class-Level Relative Abundance ---")

p_class <- plot_global_rank(
  df_long   = rel_abund_long,
  rank_col  = "Class",
  top_n     = 12,
  pal       = brewer.pal(12, "Set3"),
  title_str = "3. Relative Abundance (Class Level)",
  out_png   = "assets/070_16S_Relative_Abundance(Class).png",
  is_italic = TRUE
)

# ==============================================================================
# 3. Order-Level Composition (Figures 071, 072, 073)
# ==============================================================================
message("--- Step 3: Processing Order-Level Relative Abundance ---")

# Figure 071: Global Order
p_order <- plot_global_rank(
  df_long   = rel_abund_long,
  rank_col  = "Order",
  top_n     = 12,
  pal       = brewer.pal(12, "Set3"),
  title_str = "3. Relative Abundance (Order Level)",
  out_png   = "assets/071_16S_Relative_Abundance(Order).png",
  is_italic = TRUE
)

# Shared Palette Preparation for Regional GB and GJ Order Contrasts (Top 20)
order_sample <- rel_abund_long %>%
  mutate(Order = if_else(is.na(Order) | Order == "" | Order == "o__", "Unassigned", Order)) %>%
  group_by(`sample-id`, Treatment_Group, Order) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

top20_gb_ord <- order_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Order != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 20) %>%
  pull(Order)

top20_gj_ord <- order_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Order) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Order != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 20) %>%
  pull(Order)

all_top_orders <- unique(c(top20_gb_ord, top20_gj_ord))
base_colors_ord <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"), brewer.pal(8, "Set2"))
palette_pool_ord <- colorRampPalette(base_colors_ord)(length(all_top_orders))
master_taxa_colors_ord <- c(setNames(palette_pool_ord, all_top_orders), "Other" = "#D3D3D3")

# Figure 072: GB Lineage Order
gb_order_summary <- order_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  mutate(Order_Group = if_else(Order %in% top20_gb_ord, Order, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Order_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Order_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Order_Group = factor(Order_Group, levels = rev(c(top20_gb_ord, "Other"))))

p_order_gb <- ggplot(gb_order_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Order_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors_ord, guide = guide_legend(ncol = 1)) +
  scale_x_discrete(labels = c("GB_Control" = "GB_Control", "GB_Treatment" = "GB_Treatment")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8.5, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(title = "3. GB Relative Abundance (Order Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/072_16S_GB_Relative_Abundance(Order).png", p_order_gb, width = 8.5, height = 7.5, dpi = 300)

# Figure 073: GJ Lineage Order
gj_order_summary <- order_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  mutate(Order_Group = if_else(Order %in% top20_gj_ord, Order, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Order_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Order_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Order_Group = factor(Order_Group, levels = rev(c(top20_gj_ord, "Other"))))

p_order_gj <- ggplot(gj_order_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Order_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors_ord, guide = guide_legend(ncol = 1)) +
  scale_x_discrete(labels = c("GJ_Control" = "GJ_Control", "GJ_Treatment" = "GJ_Treatment")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8.5, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(title = "3. GJ Relative Abundance (Order Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/073_16S_GJ_Relative_Abundance(Order).png", p_order_gj, width = 8.5, height = 7.5, dpi = 300)

# ==============================================================================
# 4. Family-Level Composition (Figure 074)
# ==============================================================================
message("--- Step 4: Processing Family-Level Relative Abundance ---")

p_family <- plot_global_rank(
  df_long   = rel_abund_long,
  rank_col  = "Family",
  top_n     = 12,
  pal       = brewer.pal(12, "Set3"),
  title_str = "3. Relative Abundance (Family Level)",
  out_png   = "assets/074_16S_Relative_Abundance(Family).png",
  is_italic = TRUE
)

# ==============================================================================
# 5. Genus-Level Composition (Figures 075, 076, 077)
# ==============================================================================
message("--- Step 5: Processing Genus-Level Relative Abundance ---")

# Figure 075: Global Genus
p_genus <- plot_global_rank(
  df_long   = rel_abund_long,
  rank_col  = "Genus",
  top_n     = 12,
  pal       = brewer.pal(12, "Set3"),
  title_str = "3. Relative Abundance (Genus Level)",
  out_png   = "assets/075_16S_Relative_Abundance(Genus).png",
  is_italic = TRUE
)

# Shared Palette Preparation for Regional GB and GJ Genus Contrasts (Top 30)
genus_sample <- rel_abund_long %>%
  mutate(Genus = if_else(is.na(Genus) | Genus == "" | Genus == "g__", "Unassigned", Genus)) %>%
  group_by(`sample-id`, Treatment_Group, Genus) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop")

top30_gb_gen <- genus_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Genus != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 30) %>%
  pull(Genus)

top30_gj_gen <- genus_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  group_by(Genus) %>%
  summarise(MeanAbund = mean(RelAbund), .groups = "drop") %>%
  filter(Genus != "Unassigned") %>%
  arrange(desc(MeanAbund)) %>%
  slice_head(n = 30) %>%
  pull(Genus)

all_top_genera <- unique(c(top30_gb_gen, top30_gj_gen))
base_colors_gen <- c(
  brewer.pal(12, "Paired"),
  brewer.pal(8, "Dark2"),
  brewer.pal(8, "Set2"),
  brewer.pal(12, "Set3")
)
palette_pool_gen <- colorRampPalette(base_colors_gen)(length(all_top_genera))
master_taxa_colors_gen <- c(setNames(palette_pool_gen, all_top_genera), "Other" = "#D3D3D3")

# Figure 076: GB Lineage Genus
gb_genus_summary <- genus_sample %>%
  filter(Treatment_Group %in% c("GB_Control", "GB_Treatment")) %>%
  mutate(Genus_Group = if_else(Genus %in% top30_gb_gen, Genus, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Genus_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Genus_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Genus_Group = factor(Genus_Group, levels = rev(c(top30_gb_gen, "Other"))))

p_genus_gb <- ggplot(gb_genus_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Genus_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors_gen, guide = guide_legend(ncol = 1)) +
  scale_x_discrete(labels = c("GB_Control" = "GB_Control", "GB_Treatment" = "GB_Treatment")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(title = "3. GB Relative Abundance (Genus Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/076_16S_GB_Relative_Abundance(Genus).png", p_genus_gb, width = 8.5, height = 7.5, dpi = 300)

# Figure 077: GJ Lineage Genus
gj_genus_summary <- genus_sample %>%
  filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment")) %>%
  mutate(Genus_Group = if_else(Genus %in% top30_gj_gen, Genus, "Other")) %>%
  group_by(`sample-id`, Treatment_Group, Genus_Group) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment_Group, Genus_Group) %>%
  summarise(MeanRelAbund = mean(RelAbund), .groups = "drop") %>%
  mutate(Genus_Group = factor(Genus_Group, levels = rev(c(top30_gj_gen, "Other"))))

p_genus_gj <- ggplot(gj_genus_summary, aes(x = Treatment_Group, y = MeanRelAbund, fill = Genus_Group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8, color = "black", linewidth = 0.25) +
  scale_fill_manual(values = master_taxa_colors_gen, guide = guide_legend(ncol = 1)) +
  scale_x_discrete(labels = c("GJ_Control" = "GJ_Control", "GJ_Treatment" = "GJ_Treatment")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.01)) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    text = element_text(family = font_family),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 8, family = font_family, face = "italic"),
    axis.text.x = element_text(face = "bold", color = "black", size = 10, family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold", family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family),
    panel.grid = element_blank()
  ) +
  labs(title = "3. GJ Relative Abundance (Genus Level)", y = "Mean Relative Abundance (%)")

ggsave("assets/077_16S_GJ_Relative_Abundance(Genus).png", p_genus_gj, width = 8.5, height = 7.5, dpi = 300)

message("=== [", format(Sys.time(), "%T"), "] Taxonomic Profiling Completed Successfully (Figures 067-077 Exported) ===")