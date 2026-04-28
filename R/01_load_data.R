suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
})

# ── Main loader ───────────────────────────────────────────────────────────────

load_metabolomics_data <- function(file_path, config) {
  log_info("Loading: ", basename(file_path))

  # Peek at column count so we can force the Name column to text.
  # readxl infers mixed (numeric + string) Name columns as numeric and
  # silently coerces compound name strings to NA.
  n_cols    <- ncol(readxl::read_excel(file_path, sheet = 1, n_max = 0))
  col_types <- c("text", rep("guess", n_cols - 1))

  raw <- tryCatch(
    readxl::read_excel(file_path, sheet = 1, col_types = col_types),
    error = function(e) stop("Cannot read '", basename(file_path), "': ", e$message)
  )

  name_col    <- names(raw)[1]
  sample_cols <- names(raw)[-1]

  # Parse groups from sample names
  groups_vec <- parse_group(sample_cols, config$groups)
  unrecognised <- sample_cols[is.na(groups_vec)]
  if (length(unrecognised) > 0) {
    log_warn("Could not assign group to: ", paste(unrecognised, collapse = ", "),
             " — these columns will be dropped.")
    keep        <- !is.na(groups_vec)
    sample_cols <- sample_cols[keep]
    groups_vec  <- groups_vec[keep]
  }

  sample_meta <- tibble(
    sample = sample_cols,
    group  = groups_vec
  )

  log_info("  Groups detected: ",
           paste(names(table(groups_vec)), table(groups_vec),
                 sep = "=", collapse = ", "))

  # Build intensity matrix (rows = features, cols = samples)
  int_mat <- as.matrix(raw[, sample_cols])
  mode(int_mat) <- "numeric"
  rownames(int_mat) <- as.character(raw[[name_col]])

  # Zeros → NA  (absent, not a valid log2 intensity)
  int_mat[int_mat == 0 | is.nan(int_mat)] <- NA

  # Per-feature, per-group presence
  presence_df <- compute_presence(int_mat, sample_meta, config$presence_threshold)

  # Group means (only where present in >50 % of replicates)
  group_means <- compute_group_means(int_mat, sample_meta, presence_df)

  list(
    intensity   = int_mat,
    sample_meta = sample_meta,
    feature_ids = rownames(int_mat),
    presence    = presence_df,
    group_means = group_means,
    file_name   = basename(file_path),
    organ_name  = tools::file_path_sans_ext(basename(file_path))
  )
}

# ── Presence calculation ──────────────────────────────────────────────────────

compute_presence <- function(int_mat, sample_meta, threshold) {
  groups <- unique(sample_meta$group)

  purrr::map_dfr(groups, function(g) {
    g_cols  <- sample_meta$sample[sample_meta$group == g]
    g_mat   <- int_mat[, g_cols, drop = FALSE]
    n_total <- ncol(g_mat)

    tibble(
      feature      = rownames(int_mat),
      group        = g,
      n_present    = rowSums(!is.na(g_mat)),
      n_total      = n_total,
      prop_present = rowSums(!is.na(g_mat)) / n_total,
      is_present   = rowSums(!is.na(g_mat)) / n_total > threshold
    )
  })
}

# ── Group means (respecting presence threshold) ───────────────────────────────

compute_group_means <- function(int_mat, sample_meta, presence_df) {
  groups <- unique(sample_meta$group)

  purrr::map_dfr(groups, function(g) {
    g_cols    <- sample_meta$sample[sample_meta$group == g]
    g_mat     <- int_mat[, g_cols, drop = FALSE]
    raw_means <- rowMeans(g_mat, na.rm = TRUE)

    present_features <- presence_df %>%
      filter(group == g, is_present) %>%
      pull(feature)

    raw_means[!(names(raw_means) %in% present_features)] <- NA

    tibble(
      feature        = names(raw_means),
      group          = g,
      mean_intensity = raw_means
    )
  })
}

# ── Classify features for each pairwise comparison ───────────────────────────

classify_features <- function(presence_df, pair) {
  g1 <- pair[1]; g2 <- pair[2]

  p1 <- presence_df %>% filter(group == g1) %>% select(feature, is_present) %>%
    rename(present_g1 = is_present)
  p2 <- presence_df %>% filter(group == g2) %>% select(feature, is_present) %>%
    rename(present_g2 = is_present)

  full_join(p1, p2, by = "feature") %>%
    mutate(
      class = case_when(
        present_g1 & present_g2  ~ "A",   # quantitative
        present_g1 & !present_g2 ~ "B",   # gained in g1 / lost in g2
        !present_g1 & present_g2 ~ "B",   # gained in g2 / lost in g1
        TRUE                     ~ "absent_both"
      ),
      direction = case_when(
        class == "B" & present_g1  ~ paste0("present_in_", g1),
        class == "B" & present_g2  ~ paste0("present_in_", g2),
        TRUE ~ NA_character_
      )
    )
}

# ── QC summary table ──────────────────────────────────────────────────────────

data_qc_summary <- function(data_obj) {
  tibble(
    file          = data_obj$file_name,
    n_features    = nrow(data_obj$intensity),
    n_samples     = ncol(data_obj$intensity),
    n_groups      = length(unique(data_obj$sample_meta$group)),
    pct_missing   = round(100 * mean(is.na(data_obj$intensity)), 1)
  )
}
