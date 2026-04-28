suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(limma)
}) 

# ── Main differential analysis entry point ───────────────────────────────────

run_differential <- function(data_obj, config, organ_dir, mapping_tbl = NULL) {
  log_info("Running differential analysis for: ", data_obj$organ_name)
  diff_dir <- file.path(organ_dir, "differential")
  dir.create(file.path(diff_dir, "anova"),            recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(diff_dir, "pairwise"),         recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(diff_dir, "volcano_plots"),    recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(diff_dir, "presence_absence"), recursive = TRUE, showWarnings = FALSE)

  groups        <- unique(data_obj$sample_meta$group)
  pairs         <- all_pairs(groups)
  all_pairwise  <- list()
  all_class_b   <- list()

  # ── ANOVA (Class A features present in ≥2 groups) ─────────────────────────
  log_info("  Running one-way ANOVA + Kruskal-Wallis")
  anova_res <- run_anova(data_obj, config)
  utils::write.csv(anova_res,
                   file.path(diff_dir, "anova", "anova_results.csv"),
                   row.names = FALSE)

  # ── Pairwise comparisons ──────────────────────────────────────────────────
  log_info("  Running ", length(pairs), " pairwise comparisons (limma + Kruskal-Wallis)")

  for (pair in pairs) {
    lbl     <- pair_label(pair)
    pw_dir  <- file.path(diff_dir, "pairwise", lbl)
    dir.create(pw_dir, recursive = TRUE, showWarnings = FALSE)

    feat_class <- classify_features(data_obj$presence, pair)

    # Class A: limma
    class_a_res <- run_limma_pairwise(data_obj, pair, feat_class, config)
    if (!is.null(class_a_res) && nrow(class_a_res) > 0) {
      utils::write.csv(class_a_res,
                       file.path(pw_dir, paste0(lbl, "_classA_limma.csv")),
                       row.names = FALSE)

      kw_res <- run_kruskal_pairwise(data_obj, pair, feat_class)
      utils::write.csv(kw_res,
                       file.path(pw_dir, paste0(lbl, "_classA_kruskal.csv")),
                       row.names = FALSE)

      p_volc <- plot_volcano(class_a_res, config, lbl, mapping_tbl)
      save_plot(p_volc,
                file.path(diff_dir, "volcano_plots", paste0(lbl, "_volcano.pdf")),
                width = 9, height = 7)

      all_pairwise[[lbl]] <- class_a_res
    }

    # Class B: presence/absence Fisher test
    class_b_res <- run_fisher_classB(data_obj, pair, feat_class)
    if (!is.null(class_b_res) && nrow(class_b_res) > 0) {
      utils::write.csv(class_b_res,
                       file.path(diff_dir, "presence_absence",
                                 paste0(lbl, "_classB_fisher.csv")),
                       row.names = FALSE)
      all_class_b[[lbl]] <- class_b_res
    }
  }

  # ── UpSet plot across all pairwise comparisons ────────────────────────────
  if (length(all_pairwise) > 0 &&
      requireNamespace("UpSetR", quietly = TRUE)) {
    p_upset <- make_upset_plot(all_pairwise, config)
    if (!is.null(p_upset)) {
      pdf(file.path(diff_dir, "upset_significant_features.pdf"),
          width = 12, height = 7)
      print(p_upset)
      dev.off()
    }
  }

  log_info("  Differential analysis complete")
  invisible(list(anova = anova_res,
                 pairwise = all_pairwise,
                 class_b  = all_class_b))
}

# ── One-way ANOVA ─────────────────────────────────────────────────────────────

run_anova <- function(data_obj, config) {
  int_mat  <- data_obj$intensity
  groups   <- data_obj$sample_meta$group
  features <- rownames(int_mat)

  # Only test features present in ≥2 groups
  n_groups_present <- data_obj$presence %>%
    group_by(feature) %>%
    summarise(n_present_groups = sum(is_present), .groups = "drop")

  test_features <- n_groups_present %>%
    filter(n_present_groups >= 2) %>%
    pull(feature)

  if (length(test_features) == 0) {
    log_warn("No features present in ≥2 groups for ANOVA.")
    return(tibble())
  }

  results <- purrr::map_dfr(test_features, function(feat) {
    vals   <- int_mat[feat, ]
    grp    <- groups[!is.na(vals)]
    vals   <- vals[!is.na(vals)]
    if (length(unique(grp)) < 2) return(NULL)

    aov_p  <- tryCatch({
      fit <- aov(vals ~ grp)
      summary(fit)[[1]][["Pr(>F)"]][1]
    }, error = function(e) NA_real_)

    kw_p   <- tryCatch(
      kruskal.test(vals ~ factor(grp))$p.value,
      error = function(e) NA_real_
    )

    tibble(feature = feat, anova_p = aov_p, kruskal_p = kw_p)
  })

  results %>%
    mutate(
      anova_fdr   = p.adjust(anova_p,   method = "BH"),
      kruskal_fdr = p.adjust(kruskal_p, method = "BH"),
      significant_anova   = anova_fdr   < config$fdr_threshold,
      significant_kruskal = kruskal_fdr < config$fdr_threshold
    ) %>%
    arrange(anova_fdr)
}

# ── limma pairwise ────────────────────────────────────────────────────────────

run_limma_pairwise <- function(data_obj, pair, feat_class, config) {
  g1 <- pair[1]; g2 <- pair[2]

  class_a_features <- feat_class %>%
    filter(class == "A") %>%
    pull(feature)

  if (length(class_a_features) < 2) return(NULL)

  keep_samples <- data_obj$sample_meta$group %in% pair
  sub_meta     <- data_obj$sample_meta[keep_samples, ]
  sub_mat      <- data_obj$intensity[class_a_features, sub_meta$sample, drop = FALSE]

  # Mean-impute NAs within Class A for limma
  sub_mat <- apply(sub_mat, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })
  sub_mat[!is.finite(sub_mat)] <- 0

  group_factor <- factor(sub_meta$group, levels = pair)
  design       <- model.matrix(~0 + group_factor)
  colnames(design) <- pair
  contrast_mat <- limma::makeContrasts(
    contrasts = paste0(g1, "-", g2),
    levels    = design
  )

  fit  <- limma::lmFit(sub_mat, design)
  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2)

  tt <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "P") %>%
    rownames_to_column("feature") %>%
    as_tibble() %>%
    rename(log2FC = logFC, p_value = P.Value, adj_p = adj.P.Val,
           avg_expr = AveExpr) %>%
    mutate(
      group1     = g1,
      group2     = g2,
      class      = "A",
      significant = adj_p < config$fdr_threshold &
                    abs(log2FC) >= config$fc_threshold,
      direction  = case_when(
        significant & log2FC > 0 ~ paste0("Up in ", g1),
        significant & log2FC < 0 ~ paste0("Up in ", g2),
        TRUE ~ "NS"
      )
    )

  tt
}

# ── Kruskal-Wallis pairwise (non-parametric alternative) ─────────────────────

run_kruskal_pairwise <- function(data_obj, pair, feat_class) {
  g1 <- pair[1]; g2 <- pair[2]

  class_a_features <- feat_class %>% filter(class == "A") %>% pull(feature)
  if (length(class_a_features) == 0) return(tibble())

  keep_samples <- data_obj$sample_meta$group %in% pair
  sub_meta     <- data_obj$sample_meta[keep_samples, ]

  purrr::map_dfr(class_a_features, function(feat) {
    vals <- data_obj$intensity[feat, sub_meta$sample]
    grp  <- sub_meta$group[!is.na(vals)]
    vals <- vals[!is.na(vals)]
    if (length(unique(grp)) < 2 || length(vals) < 4) return(NULL)
    kw_p <- tryCatch(kruskal.test(vals ~ factor(grp))$p.value, error = function(e) NA_real_)
    tibble(feature = feat, group1 = g1, group2 = g2, kruskal_p = kw_p)
  }) %>%
    mutate(kruskal_fdr = p.adjust(kruskal_p, method = "BH")) %>%
    arrange(kruskal_fdr)
}

# ── Fisher exact test for Class B ─────────────────────────────────────────────

run_fisher_classB <- function(data_obj, pair, feat_class) {
  g1 <- pair[1]; g2 <- pair[2]

  class_b <- feat_class %>% filter(class == "B")
  if (nrow(class_b) == 0) return(tibble())

  g1_samples <- data_obj$sample_meta$sample[data_obj$sample_meta$group == g1]
  g2_samples <- data_obj$sample_meta$sample[data_obj$sample_meta$group == g2]

  purrr::map_dfr(class_b$feature, function(feat) {
    g1_vals  <- data_obj$intensity[feat, g1_samples]
    g2_vals  <- data_obj$intensity[feat, g2_samples]
    g1_pres  <- sum(!is.na(g1_vals)); g1_abs  <- sum(is.na(g1_vals))
    g2_pres  <- sum(!is.na(g2_vals)); g2_abs  <- sum(is.na(g2_vals))
    cont_tbl <- matrix(c(g1_pres, g1_abs, g2_pres, g2_abs), nrow = 2)
    fp       <- tryCatch(fisher.test(cont_tbl)$p.value, error = function(e) NA_real_)
    tibble(
      feature   = feat,
      group1    = g1, group2 = g2,
      n_present_g1 = g1_pres, n_total_g1 = length(g1_vals),
      n_present_g2 = g2_pres, n_total_g2 = length(g2_vals),
      direction = class_b$direction[class_b$feature == feat],
      fisher_p  = fp,
      class     = "B"
    )
  }) %>%
    mutate(fisher_fdr = p.adjust(fisher_p, method = "BH"),
           significant = fisher_fdr < 0.05) %>%
    arrange(fisher_fdr)
}

# ── Volcano plot ──────────────────────────────────────────────────────────────

plot_volcano <- function(results, config, comparison_label, mapping_tbl = NULL) {
  df <- results

  if (!is.null(mapping_tbl)) {
    df <- df %>%
      left_join(mapping_tbl %>% select(feature, compound_name), by = "feature") %>%
      mutate(label = if_else(!is.na(compound_name), compound_name, feature))
  } else {
    df <- df %>% mutate(label = feature)
  }

  df <- df %>%
    mutate(
      # Use adjusted p-value on y-axis so the threshold line perfectly
      # separates coloured from grey dots
      neg_log10_p  = -log10(pmax(adj_p, 1e-300)),
      sig_by_p     = adj_p < config$fdr_threshold,
      colour_group = case_when(
        sig_by_p & log2FC > 0 ~ "Up",
        sig_by_p & log2FC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )

  colour_map <- c(Up = "#e31a1c", Down = "#1f78b4", NS = "grey75")
  size_map   <- c(Up = 2.2, Down = 2.2, NS = 1.4)
  alpha_map  <- c(Up = 0.9, Down = 0.9, NS = 0.4)

  # Label ALL features that pass the p-value threshold
  label_df <- df %>% filter(sig_by_p)

  ggplot(df, aes(x = log2FC, y = neg_log10_p, colour = colour_group,
                 size = colour_group, alpha = colour_group)) +
    geom_point() +
    ggrepel::geom_text_repel(
      data          = label_df,
      aes(label     = label),
      colour        = "black",
      size          = 2.8,
      fontface      = "italic",
      box.padding   = 0.4,
      point.padding = 0.2,
      max.overlaps  = Inf,
      segment.colour = "grey50",
      segment.alpha = 0.6,
      show.legend   = FALSE
    ) +
    geom_hline(yintercept = -log10(config$fdr_threshold),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = c(-config$fc_threshold, config$fc_threshold),
               linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    scale_colour_manual(values = colour_map,
                        labels = c(Up   = paste0("Up in ", unique(results$group1)),
                                   Down = paste0("Up in ", unique(results$group2)),
                                   NS   = "Not significant")) +
    scale_size_manual(values  = size_map,  guide = "none") +
    scale_alpha_manual(values = alpha_map, guide = "none") +
    labs(
      title    = paste("Volcano Plot:", gsub("_", " ", comparison_label)),
      subtitle = paste0("FDR < ", config$fdr_threshold,
                        "  |  |log2FC| \u2265 ", config$fc_threshold,
                        "  |  Red = up, Blue = down"),
      x        = "log2 Fold Change",
      y        = expression(-log[10](adjusted~p-value)),
      colour   = NULL
    ) +
    theme_metabo()
}

# ── UpSet plot across comparisons ─────────────────────────────────────────────

make_upset_plot <- function(pairwise_list, config) {
  sig_lists <- purrr::map(pairwise_list, function(df) {
    df %>% filter(significant) %>% pull(feature)
  })
  sig_lists <- sig_lists[lengths(sig_lists) > 0]
  if (length(sig_lists) < 2) return(NULL)

  tryCatch(
    UpSetR::upset(
      UpSetR::fromList(sig_lists),
      order.by  = "freq",
      nsets     = length(sig_lists),
      mainbar.y.label = "Shared significant features",
      sets.x.label    = "Significant features per comparison"
    ),
    error = function(e) {
      log_warn("UpSet plot failed: ", e$message)
      NULL
    }
  )
}
