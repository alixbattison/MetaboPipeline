suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# ── Main heatmap entry point ──────────────────────────────────────────────────

run_heatmap <- function(data_obj, diff_results, config, organ_dir, mapping_tbl = NULL) {
  log_info("Generating heatmaps for: ", data_obj$organ_name)
  hm_dir <- file.path(organ_dir, "heatmaps")
  dir.create(hm_dir, recursive = TRUE, showWarnings = FALSE)

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    log_warn("ComplexHeatmap or circlize not installed. Falling back to pheatmap.")
    use_complex <- FALSE
  } else {
    use_complex <- TRUE
  }

  # ── Class A: top significant features across all pairwise comparisons
  # Use adj_p threshold only (consistent with volcano plot colouring)
  all_sig <- purrr::map_dfr(diff_results$pairwise, function(df) {
    df %>%
      filter(adj_p < config$fdr_threshold) %>%
      select(feature, log2FC, adj_p)
  }) %>%
    group_by(feature) %>%
    summarise(min_adj_p = min(adj_p), .groups = "drop") %>%
    arrange(min_adj_p) %>%
    slice_head(n = config$heatmap_top_n) %>%
    pull(feature) %>%
    # Guard against features not present in the intensity matrix
    intersect(rownames(data_obj$intensity))

  if (length(all_sig) >= 2) {
    draw_quantitative_heatmap(
      data_obj, all_sig, config, organ_dir = hm_dir,
      filename  = "quantitative_top_features.pdf",
      mapping_tbl = mapping_tbl,
      use_complex = use_complex
    )
  } else {
    log_warn("Fewer than 2 significant Class A features — skipping quantitative heatmap.")
  }

  # ── Class B: presence/absence heatmap
  all_class_b <- purrr::map_dfr(diff_results$class_b, identity)
  if (nrow(all_class_b) > 0) {
    draw_presence_heatmap(data_obj, all_class_b, config, hm_dir, mapping_tbl, use_complex)
  }

  log_info("  Heatmaps saved to ", hm_dir)
}

# ── Quantitative heatmap (Class A) ────────────────────────────────────────────

draw_quantitative_heatmap <- function(data_obj, features, config, organ_dir,
                                      filename, mapping_tbl, use_complex) {
  mat <- data_obj$intensity[features, , drop = FALSE]

  # Z-score per feature (row)
  mat_z <- t(scale(t(mat)))
  mat_z[is.nan(mat_z)] <- 0

  # Row labels: use compound name if mapped, otherwise feature ID
  row_labels <- rownames(mat_z)
  if (!is.null(mapping_tbl)) {
    name_map <- setNames(mapping_tbl$compound_name, mapping_tbl$feature)
    row_labels <- ifelse(!is.na(name_map[row_labels]) & name_map[row_labels] != "",
                         name_map[row_labels], row_labels)
  }

  group_colours <- config$group_colours

  sample_groups <- data_obj$sample_meta$group[
    match(colnames(mat_z), data_obj$sample_meta$sample)
  ]

  pdf(file.path(organ_dir, filename), width = 14, height = max(8, length(features) * 0.25))

  if (use_complex) {
    col_fun <- circlize::colorRamp2(
      c(-2, 0, 2),
      c("#4575b4", "white", "#d73027")
    )
    col_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Group = sample_groups,
      col   = list(Group = group_colours[unique(sample_groups)]),
      annotation_name_side = "left"
    )
    hm <- ComplexHeatmap::Heatmap(
      mat_z,
      name              = "Z-score",
      col               = col_fun,
      top_annotation    = col_annotation,
      row_labels        = row_labels,
      row_names_gp      = grid::gpar(fontsize = 8),
      column_names_gp   = grid::gpar(fontsize = 8),
      cluster_rows      = TRUE,
      cluster_columns   = TRUE,
      show_column_names = TRUE,
      column_title      = paste("Top", length(features), "Differential Features (Z-score) —",
                                data_obj$organ_name)
    )
    ComplexHeatmap::draw(hm)
  } else {
    ann_df <- data.frame(Group = sample_groups,
                         row.names = colnames(mat_z))
    ann_colors <- list(Group = group_colours[unique(sample_groups)])
    rownames(mat_z) <- row_labels
    pheatmap::pheatmap(
      mat_z,
      annotation_col  = ann_df,
      annotation_colors = ann_colors,
      fontsize_row    = 7,
      fontsize_col    = 7,
      main            = paste("Top", length(features), "Differential Features (Z-score) —",
                               data_obj$organ_name)
    )
  }
  dev.off()
}

# ── Presence/absence heatmap (Class B) ───────────────────────────────────────

draw_presence_heatmap <- function(data_obj, class_b_df, config, organ_dir,
                                   mapping_tbl, use_complex) {
  b_features <- unique(class_b_df$feature)
  if (length(b_features) < 2) return(invisible(NULL))
  b_features <- b_features[seq_len(min(config$heatmap_top_n, length(b_features)))]

  # Build presence matrix per group (1 = present, 0 = absent)
  pres_wide <- data_obj$presence %>%
    filter(feature %in% b_features) %>%
    select(feature, group, is_present) %>%
    pivot_wider(names_from = group, values_from = is_present) %>%
    column_to_rownames("feature")
  pres_mat <- as.matrix(pres_wide) * 1L

  row_labels <- rownames(pres_mat)
  if (!is.null(mapping_tbl)) {
    name_map <- setNames(mapping_tbl$compound_name, mapping_tbl$feature)
    row_labels <- ifelse(!is.na(name_map[row_labels]) & name_map[row_labels] != "",
                         name_map[row_labels], row_labels)
  }

  pdf(file.path(organ_dir, "presence_absence_classB.pdf"),
      width = 10, height = max(6, length(b_features) * 0.22))

  if (use_complex) {
    col_fun <- circlize::colorRamp2(c(0, 1), c("#f7f7f7", "#2166ac"))
    hm <- ComplexHeatmap::Heatmap(
      pres_mat,
      name            = "Present",
      col             = col_fun,
      row_labels      = row_labels,
      row_names_gp    = grid::gpar(fontsize = 8),
      column_names_gp = grid::gpar(fontsize = 9),
      cluster_rows    = TRUE,
      cluster_columns = FALSE,
      column_title    = paste("Presence/Absence (Class B) —", data_obj$organ_name),
      heatmap_legend_param = list(at = c(0, 1), labels = c("Absent", "Present"))
    )
    ComplexHeatmap::draw(hm)
  } else {
    rownames(pres_mat) <- row_labels
    pheatmap::pheatmap(
      pres_mat,
      color        = c("#f7f7f7", "#2166ac"),
      cluster_cols = FALSE,
      fontsize_row = 7,
      main         = paste("Presence/Absence (Class B) —", data_obj$organ_name)
    )
  }
  dev.off()
}
