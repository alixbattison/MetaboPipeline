suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(httr)
  library(jsonlite)
})

# ── Main mapping entry point ──────────────────────────────────────────────────

run_id_mapping <- function(data_obj, config, organ_dir, mz_rt_df = NULL) {
  log_info("Running ID mapping for: ", data_obj$organ_name)

  feature_ids <- data_obj$feature_ids

  # Separate named compounds from numeric mystery IDs
  is_named <- !grepl("^[0-9]+\\.?[0-9]*$", feature_ids)
  named_ids  <- feature_ids[is_named]
  numeric_ids <- feature_ids[!is_named]

  log_info("  Named features: ", length(named_ids),
           " | Numeric/unknown IDs: ", length(numeric_ids))

  results <- tibble(
    feature         = feature_ids,
    compound_name   = NA_character_,
    hmdb_id         = NA_character_,
    kegg_id         = NA_character_,
    mapping_method  = NA_character_,
    confidence      = NA_character_
  )

  # Strategy 1: name-based mapping
  if (length(named_ids) > 0) {
    log_info("  Strategy 1: name-based HMDB/KEGG mapping (", length(named_ids), " features)")
    name_map <- map_by_name(named_ids)
    idx <- match(name_map$feature, results$feature)
    results$compound_name[idx]  <- name_map$compound_name
    results$hmdb_id[idx]        <- name_map$hmdb_id
    results$kegg_id[idx]        <- name_map$kegg_id
    results$mapping_method[idx] <- name_map$mapping_method
    results$confidence[idx]     <- name_map$confidence
  }

  # Strategy 2: m/z-based mapping (if m/z + RT supplied)
  if (!is.null(mz_rt_df)) {
    log_info("  Strategy 2: m/z-based mass matching (ppm tol = ",
             config$ppm_tolerance, ")")
    mz_map <- map_by_mz(mz_rt_df, config$ppm_tolerance, config$adducts)
    unmapped <- is.na(results$hmdb_id)
    idx <- match(mz_map$feature[unmapped[match(mz_map$feature, results$feature)]], results$feature)
    idx <- idx[!is.na(idx)]
    for (i in idx) {
      feat <- results$feature[i]
      hit  <- mz_map[mz_map$feature == feat, ]
      if (nrow(hit) > 0 && is.na(results$hmdb_id[i])) {
        results$compound_name[i]  <- hit$compound_name[1]
        results$hmdb_id[i]        <- hit$hmdb_id[1]
        results$kegg_id[i]        <- hit$kegg_id[1]
        results$mapping_method[i] <- "mz_match"
        results$confidence[i]     <- hit$confidence[1]
      }
    }
  }

  # Numeric IDs with no m/z: flag as unmapped
  results <- results %>%
    mutate(
      compound_name = if_else(is.na(compound_name) & !is_named[match(feature, feature_ids)],
                              paste0("Unknown_", feature), compound_name),
      mapping_method = if_else(is.na(mapping_method), "none", mapping_method),
      confidence     = if_else(is.na(confidence),     "none", confidence)
    )

  mapped   <- results %>% filter(confidence != "none")
  unmapped <- results %>% filter(confidence == "none")

  log_info("  Mapped: ", nrow(mapped), " | Unmapped: ", nrow(unmapped))

  dir.create(file.path(organ_dir, "id_mapping"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, file.path(organ_dir, "id_mapping", "mapping_table.csv"),
                   row.names = FALSE)
  utils::write.csv(unmapped, file.path(organ_dir, "id_mapping", "unmapped_features.csv"),
                   row.names = FALSE)

  results
}

# ── Strategy 1: name-based mapping ───────────────────────────────────────────

map_by_name <- function(compound_names) {
  hmdb_results <- purrr::map_dfr(compound_names, query_hmdb_by_name)
  kegg_results <- purrr::map_dfr(compound_names, query_kegg_by_name)

  dplyr::left_join(hmdb_results, kegg_results, by = "feature") %>%
    mutate(
      compound_name  = coalesce(hmdb_name, feature),
      mapping_method = case_when(
        !is.na(hmdb_id) & !is.na(kegg_id) ~ "name_hmdb_kegg",
        !is.na(hmdb_id)                    ~ "name_hmdb",
        !is.na(kegg_id)                    ~ "name_kegg",
        TRUE                               ~ "none"
      ),
      confidence = case_when(
        mapping_method == "name_hmdb_kegg" ~ "high",
        mapping_method %in% c("name_hmdb", "name_kegg") ~ "medium",
        TRUE ~ "none"
      )
    ) %>%
    select(feature, compound_name, hmdb_id, kegg_id, mapping_method, confidence)
}

query_hmdb_by_name <- function(name) {
  base_url <- "https://hmdb.ca/metabolites.json"
  result   <- tibble(feature = name, hmdb_id = NA_character_, hmdb_name = NA_character_)

  resp <- tryCatch(
    httr::GET(base_url,
              query = list(search_query = name, search_field = "name"),
              httr::timeout(10)),
    error = function(e) NULL
  )

  if (is.null(resp) || httr::status_code(resp) != 200) return(result)

  parsed <- tryCatch(jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                        simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed) || length(parsed) == 0) return(result)

  # Take the top hit if name matches closely
  top <- parsed[[1]]
  if (!is.null(top$accession)) {
    result$hmdb_id   <- top$accession
    result$hmdb_name <- top$name %||% name
  }
  result
}

query_kegg_by_name <- function(name) {
  result <- tibble(feature = name, kegg_id = NA_character_)

  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(result)

  hits <- tryCatch(
    KEGGREST::keggFind("compound", name),
    error = function(e) NULL
  )
  if (is.null(hits) || length(hits) == 0) return(result)

  result$kegg_id <- names(hits)[1]
  result
}

# ── Strategy 2: m/z-based mapping ────────────────────────────────────────────

map_by_mz <- function(mz_rt_df, ppm_tol, adducts) {
  # mz_rt_df must have columns: feature, mz, rt
  if (!all(c("feature", "mz") %in% names(mz_rt_df))) {
    log_warn("mz_rt_df must have columns 'feature' and 'mz'. Skipping m/z mapping.")
    return(tibble(feature = character()))
  }

  ref_db <- get_hmdb_mass_reference()
  if (is.null(ref_db) || nrow(ref_db) == 0) {
    log_warn("Could not retrieve HMDB mass reference. Skipping m/z mapping.")
    return(tibble(feature = character()))
  }

  all_adducts <- c(adducts$pos, adducts$neg)

  purrr::map_dfr(seq_len(nrow(mz_rt_df)), function(i) {
    obs_mz   <- mz_rt_df$mz[i]
    feat     <- mz_rt_df$feature[i]
    best_hit <- NULL

    for (adduct_mass in all_adducts) {
      neutral_mass <- obs_mz - adduct_mass
      ppm_diff     <- abs((ref_db$monisotopic_molecular_weight - neutral_mass) /
                            neutral_mass * 1e6)
      hits <- ref_db[ppm_diff <= ppm_tol, ]
      if (nrow(hits) > 0) {
        hits$ppm_error <- ppm_diff[ppm_diff <= ppm_tol]
        hits$adduct    <- names(all_adducts)[all_adducts == adduct_mass][1]
        if (is.null(best_hit) || hits$ppm_error[1] < best_hit$ppm_error[1]) {
          best_hit <- hits[1, ]
        }
      }
    }

    if (!is.null(best_hit)) {
      tibble(
        feature       = feat,
        compound_name = best_hit$name,
        hmdb_id       = best_hit$accession,
        kegg_id       = NA_character_,
        confidence    = if_else(best_hit$ppm_error < 2, "high", "medium"),
        ppm_error     = best_hit$ppm_error
      )
    } else {
      tibble(feature = feat, compound_name = NA_character_,
             hmdb_id = NA_character_, kegg_id = NA_character_,
             confidence = NA_character_, ppm_error = NA_real_)
    }
  })
}

get_hmdb_mass_reference <- function() {
  # Download a small HMDB metabolite reference table via REST API
  # Returns a data frame with accession, name, monisotopic_molecular_weight
  url <- "https://hmdb.ca/metabolites.json?page=1"
  resp <- tryCatch(
    httr::GET(url, httr::timeout(15)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) return(NULL)

  parsed <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(parsed)) return(NULL)
  as_tibble(parsed)
}

# ── Null coalescing operator ──────────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b
