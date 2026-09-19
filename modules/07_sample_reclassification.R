#!/usr/bin/env Rscript
# ==============================================================================
# Pipeline: 04_MicroTom_P_Rhizosphere_Amplicon_Analysis
# Module: 07_sample_reclassification.R
# Target: 16S rRNA (V3-V4) & ITS2 Rhizosphere Amplicon Sequencing
# Cohort: Sensitivity Diagnostic (Extra Sample Rescue & GJ Set 4 Label-Swap Audit)
# Platform: R (Native ARM64 / macOS)
# Description: Automates 1D PCoA Axis 1 centroid thresholding (T1) to evaluate
#              reclassification of Extra technical replicates and Set 4 batch
#              swaps, audits cross-metric concordance, and exports corrected
#              diagnostic metadata artifacts.
# ==============================================================================

suppressPackageStartupMessages({
  library(qiime2R)
  library(tidyverse)
  library(vegan)
  library(ggpubr)
  library(rstatix)
})

# Global typography & directory configuration
font_family <- "Times New Roman"
if (!dir.exists("assets")) dir.create("assets", recursive = TRUE)

message("=== [", format(Sys.time(), "%T"), "] Initializing Sample Reclassification and Audit Module ===")

# ==============================================================================
# 1. Load Raw Metadata & Export Rescued Extra Metadata
# ==============================================================================
message("--- Step 1: Processing Extra Sample Reclassification Metadata ---")

if (!file.exists("metadata_P_generation.tsv")) {
  stop("[-] ERROR: Base metadata 'metadata_P_generation.tsv' not found.")
}

meta_raw <- read_tsv("metadata_P_generation.tsv", show_col_types = FALSE)

# Extra sample target vectors derived from PCoA Axis 1 decision boundaries
extra_to_treatment <- c("J7x", "J8x", "J2y", "J5y")
extra_to_control   <- c("J13x", "J5x", "J2x")

meta_corrected_extra <- meta_raw %>%
  mutate(
    Treatment_Group_Original = Treatment_Group,
    Sample_Type_Original     = Sample_Type,
    Sample_Suffix            = str_remove(`sample-id`, "^[BF]_"),
    Treatment_Group_Updated  = case_when(
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") & Sample_Suffix %in% extra_to_treatment ~ "GJ_Treatment",
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") & Sample_Suffix %in% extra_to_control   ~ "GJ_Control",
      TRUE ~ Treatment_Group
    ),
    Sample_Type_Updated      = case_when(
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") ~ "Standard",
      TRUE ~ Sample_Type
    ),
    Reclassification_Note    = case_when(
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") ~ "Rescued Extra Sample via PCoA T1 Threshold",
      TRUE ~ "Original"
    )
  ) %>%
  mutate(
    Treatment_Group = Treatment_Group_Updated,
    Sample_Type     = Sample_Type_Updated
  ) %>%
  select(-Sample_Suffix, -Treatment_Group_Updated, -Sample_Type_Updated)

write_tsv(meta_corrected_extra, "metadata_P_generation_corrected_extra.tsv")
cat("[+] Saved: metadata_P_generation_corrected_extra.tsv\n")

# ==============================================================================
# 2. Extra Sample Classification Diagnostics (1D & 2D PCoA)
# ==============================================================================
message("--- Step 2: Evaluating Extra Technical Samples via Decision Boundaries ---")

classify_extra_samples <- function(dist_path, metric_name, meta_df) {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  common_samples <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_explained <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  # Centroids & Decision Threshold (T)
  mean_ctrl <- mean(pcoa_df$PCoA1[pcoa_df$Treatment_Group == "GJ_Control"])
  mean_treat <- mean(pcoa_df$PCoA1[pcoa_df$Treatment_Group == "GJ_Treatment"])
  threshold <- (mean_ctrl + mean_treat) / 2

  extra_df <- pcoa_df %>%
    filter(Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra"))

  # Subplot 1: 1D Threshold Verification
  p1_threshold <- ggplot(
    pcoa_df %>% filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment", "GJ_Treatment_Extra", "GJ_Extra")),
    aes(x = PCoA1, y = Treatment_Group, color = Treatment_Group)
  ) +
    geom_jitter(height = 0.15, size = 2.5, alpha = 0.7) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 0.9) +
    geom_text(
      data = extra_df,
      aes(x = PCoA1, y = Treatment_Group, label = str_remove(`sample-id`, "^B_")),
      vjust = -0.9, size = 2.8, color = "black", fontface = "bold", family = font_family
    ) +
    scale_color_manual(values = c("GJ_Control" = "#D95F02", "GJ_Treatment" = "#FE9929",
                                  "GJ_Treatment_Extra" = "#7A0177", "GJ_Extra" = "#7A0177")) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.position = "none",
      axis.text = element_text(color = "black", family = font_family),
      axis.title = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11, family = font_family)
    ) +
    labs(
      title = sprintf("%s: PCoA Threshold", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = "Group"
    )

  # Subplot 2: 2D Projection with Reassigned Symbols
  pcoa_reassigned <- pcoa_df %>%
    filter(Treatment_Group %in% c("GJ_Control", "GJ_Treatment", "GJ_Treatment_Extra", "GJ_Extra")) %>%
    mutate(
      Reclassified_Group = if_else(
        Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra"),
        if_else(PCoA1 < threshold, "GJ_Treatment (Reclassified)", "GJ_Control (Reclassified)"),
        as.character(Treatment_Group)
      )
    )

  p2_reclassified <- ggplot(pcoa_reassigned, aes(x = PCoA1, y = PCoA2)) +
    stat_ellipse(
      data = filter(pcoa_reassigned, !Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra")),
      aes(color = Treatment_Group, fill = Treatment_Group),
      geom = "polygon", alpha = 0.15, level = 0.95, color = NA
    ) +
    stat_ellipse(
      data = filter(pcoa_reassigned, !Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra")),
      aes(color = Treatment_Group),
      geom = "path", level = 0.95, linewidth = 0.7
    ) +
    geom_point(
      data = filter(pcoa_reassigned, !Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra")),
      aes(color = Treatment_Group), size = 2.5, alpha = 0.8
    ) +
    geom_point(
      data = filter(pcoa_reassigned, Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra")),
      aes(shape = Reclassified_Group), color = "#7A0177", size = 3.5, stroke = 1.1
    ) +
    scale_color_manual(values = c("GJ_Control" = "#D95F02", "GJ_Treatment" = "#FE9929")) +
    scale_fill_manual(values = c("GJ_Control" = "#D95F02", "GJ_Treatment" = "#FE9929")) +
    scale_shape_manual(values = c("GJ_Control (Reclassified)" = 17, "GJ_Treatment (Reclassified)" = 15)) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family),
      axis.text = element_text(color = "black", family = font_family),
      axis.title = element_text(face = "bold", family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11, family = font_family)
    ) +
    labs(
      title = sprintf("%s: Reclassified Samples", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2])
    )

  return(list(p_threshold = p1_threshold, p_reclassified = p2_reclassified))
}

meta_gj <- meta_raw %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(Domain == "Bacteria", Treatment_Group %in% c("GJ_Control", "GJ_Treatment", "GJ_Treatment_Extra", "GJ_Extra"))

res_uuf_extra <- classify_extra_samples("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_gj)
res_jac_extra <- classify_extra_samples("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_gj)

comb_extra_verification <- ggarrange(
  res_uuf_extra$p_threshold, res_uuf_extra$p_reclassified,
  res_jac_extra$p_threshold, res_jac_extra$p_reclassified,
  ncol = 2, nrow = 2
)

final_extra_verification <- annotate_figure(
  comb_extra_verification,
  top = text_grob("2. J(Extra) Reclassification", face = "bold", size = 14, family = font_family)
)

ggsave("assets/023_16S_Extra_Reclassification_Verification.png", final_extra_verification, width = 11, height = 9, dpi = 300)

# ==============================================================================
# 3. Corrected Samples Beta Diversity Visualization (N=134 Sensitivity Model)
# ==============================================================================
message("--- Step 3: Visualizing Corrected Beta Diversity (Extra Rescued) ---")

# Use in-memory meta_corrected_extra if physical metadata_P_generation_corrected.tsv is not present
meta_corrected_cohort <- if (file.exists("metadata_P_generation_corrected.tsv")) {
  read_tsv("metadata_P_generation_corrected.tsv", show_col_types = FALSE)
} else {
  meta_corrected_extra
} %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(Domain == "Bacteria") %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")))

palette_4 <- c(
  "GB_Control"   = "#2B5C8F",
  "GB_Treatment" = "#41B6C4",
  "GJ_Control"   = "#D95F02",
  "GJ_Treatment" = "#FE9929"
)

analyze_beta_pcoa <- function(dist_path, metric_name, meta_df, title_prefix = "PCoA") {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  common_samples <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  set.seed(42)
  permanova_res <- adonis2(dist_sub ~ Treatment_Group, data = meta_sub, permutations = 999)
  dispersion <- betadisper(dist_sub, meta_sub$Treatment_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_explained <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f",
    permanova_res$R2[1],
    ifelse(permanova_res$`Pr(>F)`[1] < 0.001, "< 0.001", sprintf("%.3f", permanova_res$`Pr(>F)`[1])),
    disp_p
  )

  p_ord <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Treatment_Group, fill = Treatment_Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.95, color = NA) +
    stat_ellipse(geom = "path", level = 0.95, linewidth = 0.8) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = palette_4) +
    scale_fill_manual(values = palette_4) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
             size = 3.5, fontface = "bold", family = font_family) +
    labs(
      title = sprintf("%s - %s Distance", title_prefix, metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2])
    )

  return(p_ord)
}

p_corr_bc  <- analyze_beta_pcoa("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_corrected_cohort)
p_corr_wuf <- analyze_beta_pcoa("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_corrected_cohort)
p_corr_uuf <- analyze_beta_pcoa("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_corrected_cohort)
p_corr_jac <- analyze_beta_pcoa("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_corrected_cohort)

combined_corr_beta <- ggarrange(p_corr_bc, p_corr_wuf, p_corr_uuf, p_corr_jac,
                                ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
final_corr_beta <- annotate_figure(
  combined_corr_beta,
  top = text_grob("3. M(B), M(J), B, J - Beta Diversity (Rescued)", face = "bold", size = 15, family = font_family)
)

ggsave("assets/024_16S_Corrected_Extra_Beta_Diversity.png", final_corr_beta, width = 9.5, height = 9, dpi = 300)

# ==============================================================================
# 4. Outlier Detection and Label Swap Extraction: GJ Set 4 (jM4 & J4)
# ==============================================================================
message("--- Step 4: Extracting Potential Set 4 Batch Label-Swaps ---")

extract_all_gj_outliers <- function(dist_path, metric_name) {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  meta_filtered <- meta_corrected_extra %>%
    filter(Domain == "Bacteria", Treatment_Group %in% c("GJ_Control", "GJ_Treatment"))

  common_samples <- intersect(meta_filtered$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_filtered %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  mean_ctrl <- mean(pcoa_df$PCoA1[pcoa_df$Treatment_Group == "GJ_Control"])
  mean_treat <- mean(pcoa_df$PCoA1[pcoa_df$Treatment_Group == "GJ_Treatment"])
  threshold <- (mean_ctrl + mean_treat) / 2

  if (mean_ctrl > mean_treat) {
    ctrl_outliers <- pcoa_df %>% filter(Treatment_Group == "GJ_Control", PCoA1 < threshold)
    treat_outliers <- pcoa_df %>% filter(Treatment_Group == "GJ_Treatment", PCoA1 >= threshold)
  } else {
    ctrl_outliers <- pcoa_df %>% filter(Treatment_Group == "GJ_Control", PCoA1 >= threshold)
    treat_outliers <- pcoa_df %>% filter(Treatment_Group == "GJ_Treatment", PCoA1 < threshold)
  }

  return(list(
    ctrl_outliers  = ctrl_outliers %>% select(`sample-id`, PCoA1, Sequencing_Batch, Generation) %>% arrange(`sample-id`),
    treat_outliers = treat_outliers %>% select(`sample-id`, PCoA1, Sequencing_Batch, Generation) %>% arrange(`sample-id`),
    threshold      = threshold
  ))
}

uuf_audit <- extract_all_gj_outliers("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac")
jac_audit <- extract_all_gj_outliers("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard")

cat("\n======================================================================\n")
cat("          GJ SET 4 OUTLIER SAMPLE NAME VERIFICATION REPORT            \n")
cat("======================================================================\n")
cat("Unweighted UniFrac Control Outliers (jM4): ", paste(uuf_audit$ctrl_outliers$`sample-id`, collapse = ", "), "\n")
cat("Unweighted UniFrac Treatment Outliers (J4): ", paste(uuf_audit$treat_outliers$`sample-id`, collapse = ", "), "\n")

# ==============================================================================
# 5. Metadata Final Correction Generation
# ==============================================================================
message("--- Step 5: Generating metadata_P_generation_corrected_final.tsv ---")

jm4_swap_to_treatment <- c("jM4-2", "jM4-3", "jM4-4", "jM4-5", "jM4-6", "jM4-7")
j4_swap_to_control     <- c("J4-3", "J4-4", "J4-6", "J4-10")

meta_final <- meta_raw %>%
  mutate(
    Treatment_Group_Original = Treatment_Group,
    Sample_Type_Original     = Sample_Type,
    Sample_Suffix            = str_remove(`sample-id`, "^[BF]_"),
    Treatment_Group_Updated  = case_when(
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") & Sample_Suffix %in% extra_to_treatment ~ "GJ_Treatment",
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") & Sample_Suffix %in% extra_to_control   ~ "GJ_Control",
      Treatment_Group == "GJ_Control" & Sample_Suffix %in% jm4_swap_to_treatment                        ~ "GJ_Treatment",
      Treatment_Group == "GJ_Treatment" & Sample_Suffix %in% j4_swap_to_control                        ~ "GJ_Control",
      TRUE ~ Treatment_Group
    ),
    Sample_Type_Updated      = case_when(
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") ~ "Standard",
      TRUE ~ Sample_Type
    ),
    Reclassification_Note    = case_when(
      Treatment_Group %in% c("GJ_Treatment_Extra", "GJ_Extra") ~ "Rescued Extra Sample via PCoA T1 Threshold",
      Sample_Suffix %in% jm4_swap_to_treatment                 ~ "Corrected Set 4 Batch Swap Outlier (jM4 -> GJ_Treatment)",
      Sample_Suffix %in% j4_swap_to_control                    ~ "Corrected Set 4 Batch Swap Outlier (J4 -> GJ_Control)",
      TRUE ~ "Original"
    )
  ) %>%
  mutate(
    Treatment_Group = Treatment_Group_Updated,
    Sample_Type     = Sample_Type_Updated
  ) %>%
  select(-Sample_Suffix, -Treatment_Group_Updated, -Sample_Type_Updated)

write_tsv(meta_final, "metadata_P_generation_corrected_final.tsv")
cat("[+] Saved: metadata_P_generation_corrected_final.tsv\n")

# ==============================================================================
# 6. Sensitivity Beta Diversity: Highlighting Set 4 (with B_J4-1)
# ==============================================================================
message("--- Step 6: Plotting PCoA Highlighting Set 4 Swaps ---")

jm4_outliers_full <- c("B_jM4-2", "B_jM4-3", "B_jM4-4", "B_jM4-5", "B_jM4-6", "B_jM4-7")
j4_outliers_full  <- c("B_J4-3", "B_J4-4", "B_J4-6", "B_J4-10")
j4_non_swapped    <- c("B_J4-1")

meta_set4_marked <- meta_corrected_extra %>%
  filter(Domain == "Bacteria") %>%
  mutate(
    Display_Group = case_when(
      `sample-id` %in% jm4_outliers_full ~ "jM4 Outliers (Control -> Treatment)",
      `sample-id` %in% j4_outliers_full  ~ "J4 Outliers (Treatment -> Control)",
      `sample-id` %in% j4_non_swapped    ~ "J4-1 (Confirmed Treatment)",
      TRUE ~ as.character(Treatment_Group)
    ),
    Parent_Group = case_when(
      Display_Group == "jM4 Outliers (Control -> Treatment)" ~ "GJ_Treatment",
      Display_Group == "J4 Outliers (Treatment -> Control)"  ~ "GJ_Control",
      Display_Group == "J4-1 (Confirmed Treatment)"          ~ "GJ_Treatment",
      TRUE ~ as.character(Treatment_Group)
    )
  ) %>%
  mutate(
    Parent_Group = factor(Parent_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")),
    Display_Group = factor(Display_Group, levels = c(
      "GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment",
      "jM4 Outliers (Control -> Treatment)", "J4 Outliers (Treatment -> Control)", "J4-1 (Confirmed Treatment)"
    ))
  )

color_pal_set4 <- c(
  "GB_Control"                          = "#2B5C8F",
  "GB_Treatment"                        = "#41B6C4",
  "GJ_Control"                          = "#D95F02",
  "GJ_Treatment"                        = "#FE9929",
  "jM4 Outliers (Control -> Treatment)" = "#7A0177",
  "J4 Outliers (Treatment -> Control)"  = "#E7298A",
  "J4-1 (Confirmed Treatment)"          = "#33A02C"
)

shape_pal_set4 <- c(
  "GB_Control"                          = 16,
  "GB_Treatment"                        = 16,
  "GJ_Control"                          = 16,
  "GJ_Treatment"                        = 16,
  "jM4 Outliers (Control -> Treatment)" = 17,
  "J4 Outliers (Treatment -> Control)"  = 15,
  "J4-1 (Confirmed Treatment)"          = 18
)

plot_pcoa_set4 <- function(dist_path, metric_name, meta_df) {
  raw_dist <- read_qza(dist_path)$data
  dist_mat <- as.matrix(raw_dist)

  common_samples <- intersect(meta_df$`sample-id`, rownames(dist_mat))
  dist_sub <- as.dist(dist_mat[common_samples, common_samples])

  meta_sub <- meta_df %>%
    filter(`sample-id` %in% common_samples) %>%
    arrange(match(`sample-id`, common_samples))

  set.seed(42)
  permanova_res <- adonis2(dist_sub ~ Parent_Group, data = meta_sub, permutations = 999)
  dispersion <- betadisper(dist_sub, meta_sub$Parent_Group)
  disp_p <- permutest(dispersion, permutations = 999)$tab$`Pr(>F)`[1]

  pcoa_cmd <- cmdscale(dist_sub, k = 2, eig = TRUE)
  var_explained <- round(100 * (pcoa_cmd$eig[1:2] / sum(pcoa_cmd$eig[pcoa_cmd$eig > 0])), 2)

  pcoa_df <- as.data.frame(pcoa_cmd$points) %>%
    rownames_to_column("sample-id") %>%
    rename(PCoA1 = V1, PCoA2 = V2) %>%
    inner_join(meta_sub, by = "sample-id")

  annot_text <- sprintf(
    "PERMANOVA R² = %.3f\nPERMANOVA p = %s\nDispersion p = %.3f",
    permanova_res$R2[1],
    ifelse(permanova_res$`Pr(>F)`[1] < 0.001, "< 0.001", sprintf("%.3f", permanova_res$`Pr(>F)`[1])),
    disp_p
  )

  p_ord <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2)) +
    stat_ellipse(
      aes(group = Parent_Group, fill = Parent_Group, color = after_scale(fill)),
      geom = "polygon", alpha = 0.12, level = 0.95, linewidth = 0.7, show.legend = FALSE
    ) +
    geom_point(aes(color = Display_Group, shape = Display_Group), size = 2.8, stroke = 1.0, alpha = 0.85) +
    scale_color_manual(values = color_pal_set4, name = NULL) +
    scale_fill_manual(values = palette_4, guide = "none") +
    scale_shape_manual(values = shape_pal_set4, name = NULL) +
    theme_bw(base_size = 12, base_family = font_family) +
    theme(
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(family = font_family),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, family = font_family),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", x = -Inf, y = Inf, label = annot_text, hjust = -0.05, vjust = 1.1,
             size = 3.4, fontface = "bold", family = font_family) +
    labs(
      title = sprintf("PCoA - %s Distance", metric_name),
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2])
    )

  return(p_ord)
}

p_set4_bc  <- plot_pcoa_set4("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_set4_marked)
p_set4_wuf <- plot_pcoa_set4("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_set4_marked)
p_set4_uuf <- plot_pcoa_set4("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_set4_marked)
p_set4_jac <- plot_pcoa_set4("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_set4_marked)

combined_set4_beta <- ggarrange(p_set4_bc, p_set4_wuf, p_set4_uuf, p_set4_jac,
                                ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
final_set4_beta <- annotate_figure(
  combined_set4_beta,
  top = text_grob("5. M(B), M(J), B, J - Beta Diversity (GJ_Set4 Marked)", face = "bold", size = 14, family = font_family)
)

ggsave("assets/025_16S_GJ_Set4_Marked_Beta_Diversity.png", final_set4_beta, width = 12, height = 10.5, dpi = 300)

# ==============================================================================
# 7. Final Corrected Diversity Visualization (Hypothesis Testing Climax)
# ==============================================================================
message("--- Step 7: Visualizing Final Corrected Diversity Models ---")

meta_final_cohort <- meta_final %>%
  filter(`sample-id` != "#q2:types") %>%
  filter(Domain == "Bacteria", Treatment_Group %in% c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")) %>%
  mutate(Treatment_Group = factor(Treatment_Group, levels = c("GB_Control", "GB_Treatment", "GJ_Control", "GJ_Treatment")))

# 7.1 Final Alpha Diversity
obs_final <- read_qza("07_diversity_16S/full/observed_features_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Observed_Features = observed_features)
shannon_final <- read_qza("07_diversity_16S/full/shannon_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Shannon = shannon_entropy)
faith_final <- read_qza("07_diversity_16S/full/faith_pd_vector.qza")$data %>%
  rownames_to_column("sample-id") %>% rename(Faith_PD = faith_pd)

rarefied_tbl_final <- read_qza("07_diversity_16S/full/rarefied_table.qza")$data
simpson_final <- diversity(t(rarefied_tbl_final), index = "simpson") %>%
  enframe(name = "sample-id", value = "Simpson")

alpha_final_long <- meta_final_cohort %>%
  inner_join(obs_final, by = "sample-id") %>%
  inner_join(shannon_final, by = "sample-id") %>%
  inner_join(simpson_final, by = "sample-id") %>%
  inner_join(faith_final, by = "sample-id") %>%
  pivot_longer(cols = c("Observed_Features", "Shannon", "Simpson", "Faith_PD"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Observed_Features", "Shannon", "Simpson", "Faith_PD")))

p_alpha_final <- ggplot(alpha_final_long, aes(x = Treatment_Group, y = Value, fill = Treatment_Group)) +
  geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.8, aes(color = Treatment_Group)) +
  facet_wrap(~Metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = palette_4) +
  scale_color_manual(values = palette_4) +
  scale_x_discrete(labels = c(
    "GB_Control"   = "GB_Control (N=33)",
    "GB_Treatment" = "GB_Treatment (N=33)",
    "GJ_Control"   = "GJ_Control (N=34)",
    "GJ_Treatment" = "GJ_Treatment (N=34)"
  )) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "top", size = 3.8, family = font_family) +
  theme_bw(base_size = 12, base_family = font_family) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold", color = "black", family = font_family),
    axis.text.y = element_text(color = "black", family = font_family),
    strip.background = element_rect(fill = "#EFEFEF", color = "black"),
    strip.text = element_text(face = "bold", size = 11, family = font_family),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14, family = font_family)
  ) +
  labs(title = "5. M(B), M(J), B, J - Alpha Diversity (Reclassified)", y = "Calculated Diversity Value")

ggsave("assets/026_16S_Final_Corrected_Alpha_Diversity.png", p_alpha_final, width = 9.5, height = 7.5, dpi = 300)

# 7.2 Final Beta Diversity
p_fin_bc  <- analyze_beta_pcoa("07_diversity_16S/full/bray_curtis_distance_matrix.qza", "Bray-Curtis", meta_final_cohort)
p_fin_wuf <- analyze_beta_pcoa("07_diversity_16S/full/weighted_unifrac_distance_matrix.qza", "Weighted UniFrac", meta_final_cohort)
p_fin_uuf <- analyze_beta_pcoa("07_diversity_16S/full/unweighted_unifrac_distance_matrix.qza", "Unweighted UniFrac", meta_final_cohort)
p_fin_jac <- analyze_beta_pcoa("07_diversity_16S/full/jaccard_distance_matrix.qza", "Jaccard", meta_final_cohort)

combined_final_beta <- ggarrange(p_fin_bc, p_fin_wuf, p_fin_uuf, p_fin_jac,
                                 ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")
final_annotated_beta <- annotate_figure(
  combined_final_beta,
  top = text_grob("5. M(B), M(J), B, J - Beta Diversity (Reclassified)", face = "bold", size = 15, family = font_family)
)

ggsave("assets/027_16S_Final_Corrected_Beta_Diversity.png", final_annotated_beta, width = 9.5, height = 9, dpi = 300)

message("=== [", format(Sys.time(), "%T"), "] Sample Reclassification and Audit Completed Successfully ===")