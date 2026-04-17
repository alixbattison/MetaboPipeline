suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggrepel)
})

# ── Main PCA entry point ──────────────────────────────────────────────────────

run_pca <- function(data_obj, config, organ_dir) {
  log_info("Running PCA for: ", data_obj$organ_name)
  pca_dir <- file.path(organ_dir, "pca")
  dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(config$seed)

  # Use features present in >50 % of replicates in at least one group
  keep_features <- data_obj$presence %>%
    group_by(feature) %>%
    summarise(any_present = any(is_present), .groups = "drop") %>%
    filter(any_present) %>%
    pull(feature)

  mat <- data_obj$intensity[keep_features, , drop = FALSE]

  # Mean-impute remaining NAs within each sample for PCA only
  mat_imp <- apply(mat, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })

  if (any(is.nan(mat_imp) | !is.finite(mat_imp))) {
    mat_imp[!is.finite(mat_imp)] <- 0
  }

  # PCA on transposed matrix (samples as rows)
  pca_res <- prcomp(t(mat_imp), center = TRUE, scale. = TRUE)

  var_exp <- summary(pca_res)$importance["Proportion of Variance", ] * 100

  scores <- as.data.frame(pca_res$x) %>%
    rownames_to_column("sample") %>%
    left_join(data_obj$sample_meta, by = "sample")

  loadings <- as.data.frame(pca_res$rotation) %>%
    rownames_to_column("feature")

  # Plots
  p12     <- plot_pca_scores(scores, var_exp, 1, 2, config, data_obj$organ_name)
  p13     <- plot_pca_scores(scores, var_exp, 1, 3, config, data_obj$organ_name)
  p_scree <- plot_scree(var_exp, data_obj$organ_name)
  p_load  <- plot_loadings(loadings, var_exp, data_obj$organ_name)

  save_plot(p12,     file.path(pca_dir, "pca_PC1_PC2.pdf"),   width = 9,  height = 7)
  save_plot(p13,     file.path(pca_dir, "pca_PC1_PC3.pdf"),   width = 9,  height = 7)
  save_plot(p_scree, file.path(pca_dir, "scree_plot.pdf"),     width = 7,  height = 5)
  save_plot(p_load,  file.path(pca_dir, "loadings_plot.pdf"), width = 9,  height = 7)

  utils::write.csv(scores,   file.path(pca_dir, "pca_scores.csv"),   row.names = FALSE)
  utils::write.csv(loadings, file.path(pca_dir, "pca_loadings.csv"), row.names = FALSE)

  log_info("  PCA plots saved to ", pca_dir)
  invisible(list(pca = pca_res, scores = scores, loadings = loadings, var_exp = var_exp))
}

# ── Score biplot ──────────────────────────────────────────────────────────────

plot_pca_scores <- function(scores, var_exp, pc_x, pc_y, config, organ_name) {
  xcol <- paste0("PC", pc_x)
  ycol <- paste0("PC", pc_y)
  xl   <- sprintf("PC%d (%.1f%%)", pc_x, var_exp[pc_x])
  yl   <- sprintf("PC%d (%.1f%%)", pc_y, var_exp[pc_y])

  group_colours <- config$group_colours[unique(scores$group)]

  ggplot(scores, aes(x = .data[[xcol]], y = .data[[ycol]],
                     colour = group, shape = group, label = sample)) +
    stat_ellipse(aes(group = group), type = "norm", level = 0.95,
                 linetype = "dashed", alpha = 0.5, show.legend = FALSE) +
    geom_point(size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 15,
                             show.legend = FALSE, segment.alpha = 0.4) +
    scale_colour_manual(values = group_colours,
                        labels = config$group_labels[names(group_colours)]) +
    scale_shape_manual(values = setNames(seq_along(group_colours),
                                         names(group_colours)),
                       labels = config$group_labels[names(group_colours)]) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dotted") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dotted") +
    labs(title = paste("PCA Score Plot —", organ_name),
         x = xl, y = yl, colour = "Group", shape = "Group") +
    theme_metabo()
}

# ── Scree plot ────────────────────────────────────────────────────────────────

plot_scree <- function(var_exp, organ_name, n_pcs = 10) {
  n_show <- min(n_pcs, length(var_exp))
  df <- tibble(
    PC       = factor(paste0("PC", seq_len(n_show)), levels = paste0("PC", seq_len(n_show))),
    Variance = var_exp[seq_len(n_show)],
    Cumulative = cumsum(var_exp[seq_len(n_show)])
  )

  ggplot(df, aes(x = PC)) +
    geom_col(aes(y = Variance), fill = "#4292c6", alpha = 0.85) +
    geom_line(aes(y = Cumulative, group = 1), colour = "#cb181d", linewidth = 1) +
    geom_point(aes(y = Cumulative), colour = "#cb181d", size = 2.5) +
    scale_y_continuous(
      name     = "Variance explained (%)",
      sec.axis = sec_axis(~., name = "Cumulative variance (%)")
    ) +
    labs(title = paste("Scree Plot —", organ_name), x = NULL) +
    theme_metabo()
}

# ── Loadings plot (top contributors) ─────────────────────────────────────────

plot_loadings <- function(loadings, var_exp, organ_name, top_n = 20) {
  contrib <- loadings %>%
    mutate(contrib = PC1^2 + PC2^2) %>%
    arrange(desc(contrib)) %>%
    slice_head(n = top_n) %>%
    mutate(feature = factor(feature, levels = rev(feature)))

  ggplot(contrib, aes(y = feature)) +
    geom_segment(aes(x = 0, xend = PC1, yend = feature),
                 colour = "#4292c6", linewidth = 0.8) +
    geom_point(aes(x = PC1), colour = "#4292c6", size = 2.5) +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    labs(title = paste("Top", top_n, "PC1 Loadings —", organ_name),
         x = sprintf("PC1 loading (%.1f%% variance)", var_exp[1]),
         y = "Feature") +
    theme_metabo()
}
