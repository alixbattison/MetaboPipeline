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

  is_named    <- !grepl("^[0-9]+\\.?[0-9]*$", feature_ids)
  named_ids   <- feature_ids[is_named]
  numeric_ids <- feature_ids[!is_named]

  log_info("  Named features: ", length(named_ids),
           " | Numeric/unknown IDs: ", length(numeric_ids))

  results <- tibble(
    feature        = feature_ids,
    compound_name  = NA_character_,
    hmdb_id        = NA_character_,
    kegg_id        = NA_character_,
    mapping_method = NA_character_,
    confidence     = NA_character_
  )

  # Strategy 1: bulk name-based mapping (single API calls, then local matching)
  if (length(named_ids) > 0) {
    log_info("  Strategy 1: bulk name-based mapping (", length(named_ids), " features)")
    name_map <- map_by_name_bulk(named_ids, config)
    idx <- match(name_map$feature, results$feature)
    results$compound_name[idx]  <- name_map$compound_name
    results$hmdb_id[idx]        <- name_map$hmdb_id
    results$kegg_id[idx]        <- name_map$kegg_id
    results$mapping_method[idx] <- name_map$mapping_method
    results$confidence[idx]     <- name_map$confidence
  }

  # Strategy 2: m/z-based mapping
  if (!is.null(mz_rt_df)) {
    log_info("  Strategy 2: m/z-based mass matching (ppm tol = ",
             config$ppm_tolerance, ")")
    mz_map <- map_by_mz(mz_rt_df, config$ppm_tolerance, config$adducts)
    for (i in seq_len(nrow(mz_map))) {
      feat <- mz_map$feature[i]
      ri   <- match(feat, results$feature)
      if (!is.na(ri) && is.na(results$hmdb_id[ri])) {
        results$compound_name[ri]  <- mz_map$compound_name[i]
        results$hmdb_id[ri]        <- mz_map$hmdb_id[i]
        results$kegg_id[ri]        <- mz_map$kegg_id[i]
        results$mapping_method[ri] <- "mz_match"
        results$confidence[ri]     <- mz_map$confidence[i]
      }
    }
  }

  results <- results %>%
    mutate(
      compound_name  = if_else(is.na(compound_name),
                               if_else(!is_named[match(feature, feature_ids)],
                                       paste0("Unknown_", feature), feature),
                               compound_name),
      mapping_method = if_else(is.na(mapping_method), "none", mapping_method),
      confidence     = if_else(is.na(confidence),     "none", confidence)
    )

  mapped   <- results %>% filter(confidence != "none")
  unmapped <- results %>% filter(confidence == "none")
  log_info("  Mapped: ", nrow(mapped), " | Unmapped: ", nrow(unmapped))

  dir.create(file.path(organ_dir, "id_mapping"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results,  file.path(organ_dir, "id_mapping", "mapping_table.csv"),
                   row.names = FALSE)
  utils::write.csv(unmapped, file.path(organ_dir, "id_mapping", "unmapped_features.csv"),
                   row.names = FALSE)
  results
}

# ── Strategy 1: bulk name-based mapping ──────────────────────────────────────
# Fetches the entire KEGG compound list in ONE call and the entire cached HMDB
# lookup, then matches locally — avoids 1000s of individual API requests.

map_by_name_bulk <- function(compound_names, config) {
  cache_dir <- here::here("data", ".mapping_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  kegg_lookup <- get_kegg_lookup(cache_dir)
  hmdb_lookup <- get_hmdb_lookup(cache_dir)

  norm <- function(x) tolower(trimws(gsub("[^a-zA-Z0-9 ]", "", x)))

  names_norm <- norm(compound_names)

  results <- purrr::map_dfr(seq_along(compound_names), function(i) {
    if (i %% 200 == 0)
      log_info("    Mapping progress: ", i, "/", length(compound_names))

    name <- compound_names[i]
    nn   <- names_norm[i]

    kegg_id <- NA_character_
    hmdb_id <- NA_character_
    hmdb_nm <- NA_character_

    # KEGG match
    if (!is.null(kegg_lookup)) {
      ki <- match(nn, kegg_lookup$name_norm)
      if (!is.na(ki)) kegg_id <- kegg_lookup$kegg_id[ki]
    }

    # HMDB match
    if (!is.null(hmdb_lookup)) {
      hi <- match(nn, hmdb_lookup$name_norm)
      if (!is.na(hi)) {
        hmdb_id <- hmdb_lookup$hmdb_id[hi]
        hmdb_nm <- hmdb_lookup$name[hi]
      }
    }

    tibble(
      feature       = name,
      compound_name = coalesce(hmdb_nm, name),
      hmdb_id       = hmdb_id,
      kegg_id       = kegg_id,
      mapping_method = case_when(
        !is.na(hmdb_id) & !is.na(kegg_id) ~ "name_hmdb_kegg",
        !is.na(hmdb_id)                   ~ "name_hmdb",
        !is.na(kegg_id)                   ~ "name_kegg",
        TRUE                              ~ "none"
      ),
      confidence = case_when(
        mapping_method == "name_hmdb_kegg"              ~ "high",
        mapping_method %in% c("name_hmdb", "name_kegg") ~ "medium",
        TRUE                                             ~ "none"
      )
    )
  })

  results
}

# ── KEGG bulk lookup (one API call for ~18 000 compounds) ────────────────────

get_kegg_lookup <- function(cache_dir) {
  cache_file <- file.path(cache_dir, "kegg_compounds.rds")

  if (file.exists(cache_file)) {
    log_info("  Loading KEGG compound list from cache")
    return(readRDS(cache_file))
  }

  if (!requireNamespace("KEGGREST", quietly = TRUE)) {
    log_warn("KEGGREST not installed — skipping KEGG name mapping.")
    return(NULL)
  }

  log_info("  Downloading full KEGG compound list (one-time, ~30 s) ...")
  all_cpds <- tryCatch(
    KEGGREST::keggList("compound"),
    error = function(e) {
      log_warn("KEGG compound list download failed: ", e$message)
      NULL
    }
  )
  if (is.null(all_cpds)) return(NULL)

  # keggList returns named vector: names = "cpd:C00001", values = "Water; H2O"
  lookup <- tibble(
    kegg_id   = sub("^cpd:", "", names(all_cpds)),
    name_raw  = as.character(all_cpds)
  ) %>%
    mutate(
      # KEGG names are semicolon-separated; take the first as primary
      name_primary = trimws(sub(";.*", "", name_raw)),
      name_norm    = tolower(trimws(gsub("[^a-zA-Z0-9 ]", "", name_primary)))
    )

  saveRDS(lookup, cache_file)
  log_info("  KEGG lookup cached (", nrow(lookup), " compounds)")
  lookup
}

# ── HMDB bulk lookup (paginated download, cached locally) ────────────────────

get_hmdb_lookup <- function(cache_dir) {
  cache_file <- file.path(cache_dir, "hmdb_compounds.rds")

  if (file.exists(cache_file)) {
    log_info("  Loading HMDB compound list from cache")
    return(readRDS(cache_file))
  }

  log_info("  Downloading HMDB compound list (one-time, may take 1-2 min) ...")

  # HMDB REST API — paginated, 10 metabolites per page
  # Fetch first 5000 entries (covers the most common plasma/tissue metabolites)
  max_pages <- 500
  all_rows  <- vector("list", max_pages)
  fetched   <- 0

  for (page in seq_len(max_pages)) {
    resp <- tryCatch(
      httr::GET("https://hmdb.ca/metabolites.json",
                query   = list(page = page),
                httr::timeout(20)),
      error = function(e) NULL
    )
    if (is.null(resp) || httr::status_code(resp) != 200) break

    parsed <- tryCatch(
      jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                         simplifyDataFrame = TRUE),
      error = function(e) NULL
    )
    if (is.null(parsed) || length(parsed) == 0 ||
        (is.data.frame(parsed) && nrow(parsed) == 0)) break

    if (is.data.frame(parsed)) {
      all_rows[[page]] <- parsed %>%
        select(any_of(c("accession", "name"))) %>%
        rename(hmdb_id = accession)
      fetched <- fetched + nrow(parsed)
    }

    if (page %% 50 == 0)
      log_info("    HMDB pages fetched: ", page, " (", fetched, " compounds)")
  }

  all_rows <- Filter(Negate(is.null), all_rows)
  if (length(all_rows) == 0) {
    log_warn("Could not download HMDB compound list — HMDB mapping disabled.")
    return(NULL)
  }

  lookup <- bind_rows(all_rows) %>%
    distinct(hmdb_id, .keep_all = TRUE) %>%
    mutate(name_norm = tolower(trimws(gsub("[^a-zA-Z0-9 ]", "", name))))

  saveRDS(lookup, cache_file)
  log_info("  HMDB lookup cached (", nrow(lookup), " compounds)")
  lookup
}

# ── Strategy 2: m/z-based mapping ────────────────────────────────────────────

map_by_mz <- function(mz_rt_df, ppm_tol, adducts) {
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
      neutral  <- obs_mz - adduct_mass
      ppm_err  <- abs((ref_db$monisotopic_molecular_weight - neutral) / neutral * 1e6)
      hits     <- ref_db[ppm_err <= ppm_tol, ]
      if (nrow(hits) > 0) {
        hits$ppm_error <- ppm_err[ppm_err <= ppm_tol]
        if (is.null(best_hit) || hits$ppm_error[1] < best_hit$ppm_error[1])
          best_hit <- hits[1, ]
      }
    }

    if (!is.null(best_hit)) {
      tibble(feature = feat, compound_name = best_hit$name,
             hmdb_id = best_hit$accession, kegg_id = NA_character_,
             confidence = if_else(best_hit$ppm_error < 2, "high", "medium"),
             ppm_error  = best_hit$ppm_error)
    } else {
      tibble(feature = feat, compound_name = NA_character_,
             hmdb_id = NA_character_, kegg_id = NA_character_,
             confidence = NA_character_, ppm_error = NA_real_)
    }
  })
}

get_hmdb_mass_reference <- function() {
  resp <- tryCatch(
    httr::GET("https://hmdb.ca/metabolites.json?page=1", httr::timeout(15)),
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

`%||%` <- function(a, b) if (!is.null(a)) a else b
