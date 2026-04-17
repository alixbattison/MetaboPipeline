suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ── Main QC entry point ───────────────────────────────────────────────────────

run_qc <- function(data_obj, config, organ_dir) {
  log_info("Running QC for: ", data_obj$organ_name)
  qc_dir <- file.path(organ_dir, "qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  p_box     <- plot_sample_boxplots(data_obj, config)
  p_cv      <- plot_cv_per_group(data_obj, config)
  p_missing <- plot_missing_value_heatmap(data_obj, config)

  save_plot(p_box,     file.path(qc_dir, "sample_boxplots.pdf"),     width = 14, height = 6)
  save_plot(p_cv,      file.path(qc_dir, "cv_per_group.pdf"),         width = 10, height = 6)
  save_plot(p_missing, file.path(qc_dir, "missing_value_heatmap.pdf"), width = 12, height = 8)

  cv_tbl <- compute_cv_table(data_obj, config)
  utils::write.csv(cv_tbl, file.path(qc_dir, "cv_summary.csv"), row.names = FALSE)

  log_info("  QC plots saved to ", qc_dir)
  invisible(list(cv_table = cv_tbl))
}

# ── Sample intensity boxplots ─────────────────────────────────────────────────

plot_sample_boxplots <- function(data_obj, config) {
  long_df <- as.data.frame(data_obj$intensity) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "intensity") %>%
    left_join(data_obj$sample_meta, by = "sample") %>%
    filter(!is.na(intensity))

  group_colours <- config$group_colours[unique(long_df$group)]

  ggplot(long_df, aes(x = sample, y = intensity, fill = group)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
    scale_fill_manual(values = group_colours,
                      labels = config$group_labels[names(group_colours)]) +
    labs(title = paste("Sample Intensity Distributions —", data_obj$organ_name),
         x = NULL, y = "log2 Intensity", fill = "Group") +
    theme_metabo() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
}

# ── CV per group ──────────────────────────────────────────────────────────────

compute_cv_table <- function(data_obj, config) {
  groups <- unique(data_obj$sample_meta$group)
  purrr::map_dfr(groups, function(g) {
    g_cols <- data_obj$sample_meta$sample[data_obj$sample_meta$group == g]
    g_mat  <- data_obj$intensity[, g_cols, drop = FALSE]
    # CV on log2 scale: SD / mean × 100 (approximate)
    m  <- rowMeans(g_mat, na.rm = TRUE)
    s  <- apply(g_mat, 1, sd, na.rm = TRUE)
    cv <- s / m * 100
    tibble(feature = rownames(data_obj$intensity), group = g, mean = m, sd = s, cv_pct = cv)
  }) %>%
    filter(!is.na(cv_pct), is.finite(cv_pct))
}

plot_cv_per_group <- function(data_obj, config) {
  cv_df <- compute_cv_table(data_obj, config)
  group_colours <- config$group_colours[unique(cv_df$group)]

  ggplot(cv_df, aes(x = group, y = cv_pct, fill = group)) +
    geom_violin(alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.size = 0.4, fill = "white") +
    scale_fill_manual(values = group_colours,
                      labels = config$group_labels[names(group_colours)]) +
    labs(title = paste("Coefficient of Variation —", data_obj$organ_name),
         x = "Group", y = "CV (%) on log2 scale", fill = "Group") +
    theme_metabo() +
    theme(legend.position = "none")
}

# ── Missing value heatmap ─────────────────────────────────────────────────────

plot_missing_value_heatmap <- function(data_obj, config) {
  # Show presence/absence per feature per group (summarised)
  pres_wide <- data_obj$presence %>%
    select(feature, group, prop_present) %>%
    tidyr::pivot_wider(names_from = group, values_from = prop_present) %>%
    tibble::column_to_rownames("feature")

  # Only plot a manageable number of features (top 60 by total missingness)
  total_present <- rowSums(pres_wide, na.rm = TRUE)
  top_features  <- names(sort(total_present))[seq_len(min(60, nrow(pres_wide)))]
  plot_df <- pres_wide[top_features, , drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "group", values_to = "prop_present")

  ggplot(plot_df, aes(x = group, y = feature, fill = prop_present)) +
    geom_tile(colour = "white") +
    scale_fill_gradient2(low = "#d73027", mid = "#fee08b", high = "#1a9850",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Proportion\npresent") +
    labs(title = paste("Feature Presence per Group —", data_obj$organ_name),
         x = "Group", y = "Feature") +
    theme_metabo() +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(angle = 30, hjust = 1))
}
