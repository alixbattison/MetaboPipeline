suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

# ── Logging ──────────────────────────────────────────────────────────────────

log_info  <- function(...) message("[INFO]  ", ..., appendLF = TRUE)
log_warn  <- function(...) message("[WARN]  ", ..., appendLF = TRUE)
log_error <- function(...) message("[ERROR] ", ..., appendLF = TRUE)

# ── Directory helpers ─────────────────────────────────────────────────────────

make_output_dirs <- function(base_dir, organ_name) {
  dirs <- file.path(base_dir, organ_name, c(
    "qc",
    "pca",
    "id_mapping",
    "differential/anova",
    "differential/pairwise",
    "differential/volcano_plots",
    "differential/presence_absence",
    "heatmaps",
    "pathway_enrichment/kegg",
    "pathway_enrichment/hmdb",
    "pathway_enrichment/reactome",
    "pathway_enrichment/wikipathways"
  ))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  file.path(base_dir, organ_name)
}

out_path <- function(organ_dir, ...) file.path(organ_dir, ...)

# ── Group parsing ─────────────────────────────────────────────────────────────

parse_group <- function(sample_names, groups) {
  # Build alternation in priority order (longest first prevents partial matches)
  pattern <- paste0("(", paste(groups, collapse = "|"), ")(?=\\d)")
  str_extract(sample_names, pattern)
}

# ── Pairwise combination helper ───────────────────────────────────────────────

all_pairs <- function(groups) {
  combn(groups, 2, simplify = FALSE)
}

pair_label <- function(pair) paste(pair[1], "vs", pair[2], sep = "_")

# ── Safe save wrappers ────────────────────────────────────────────────────────

save_plot <- function(plot, path, width = 10, height = 8, dpi = 300) {
  tryCatch(
    ggplot2::ggsave(path, plot = plot, width = width, height = height,
                    dpi = dpi, bg = "white"),
    error = function(e) log_warn("Could not save plot to ", path, ": ", e$message)
  )
}

save_table <- function(df, path) {
  tryCatch(
    readr_write(df, path),
    error = function(e) {
      tryCatch(
        utils::write.csv(df, path, row.names = FALSE),
        error = function(e2) log_warn("Could not save table to ", path)
      )
    }
  )
}

readr_write <- function(df, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(df, path)
  } else {
    utils::write.csv(df, path, row.names = FALSE)
  }
}

# ── Theme ─────────────────────────────────────────────────────────────────────

theme_metabo <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.background  = ggplot2::element_rect(fill = "grey92", colour = NA),
      legend.position   = "right",
      plot.title        = ggplot2::element_text(face = "bold", size = 13),
      axis.title        = ggplot2::element_text(size = 11)
    )
}
