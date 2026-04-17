suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(stringr)
  library(httr)
  library(jsonlite)
  library(purrr)
})

# Bump this string any time the mapping logic changes — old caches are
# automatically discarded and rebuilt from the APIs.
.CACHE_VERSION <- "v3"

# ── Cache helpers ─────────────────────────────────────────────────────────────

load_cache <- function(cache_file) {
  if (!file.exists(cache_file)) return(list())
  cache <- tryCatch(readRDS(cache_file), error = function(e) list())
  # Discard if version doesn't match
  if (!identical(cache[["__version__"]], .CACHE_VERSION)) {
    log_info("  Cache version mismatch — rebuilding: ", basename(cache_file))
    return(list())
  }
  cache
}

save_cache <- function(cache, cache_file) {
  cache[["__version__"]] <- .CACHE_VERSION
  saveRDS(cache, cache_file)
}

# ── Main mapping entry point ──────────────────────────────────────────────────

run_id_mapping <- function(data_obj, config, organ_dir, mz_rt_df = NULL) {
  log_info("Running ID mapping for: ", data_obj$organ_name)

  feature_ids <- data_obj$feature_ids
  cache_dir   <- here::here("data", ".mapping_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  # Clean names: strip .mol, .sdf artifacts and normalise whitespace
  clean_ids <- clean_compound_names(feature_ids)

  is_named    <- !grepl("^[0-9]+\\.?[0-9]*$", clean_ids)
  named_clean <- clean_ids[is_named]
  named_raw   <- feature_ids[is_named]

  log_info("  Named features: ", length(named_clean),
           " | Numeric/unknown IDs: ", sum(!is_named))

  results <- tibble(
    feature        = feature_ids,
    clean_name     = clean_ids,
    compound_name  = NA_character_,
    hmdb_id        = NA_character_,
    kegg_id        = NA_character_,
    lmid           = NA_character_,
    lipid_class    = NA_character_,
    lipid_subclass = NA_character_,
    inchikey       = NA_character_,
    mapping_method = NA_character_,
    confidence     = NA_character_
  )

  if (length(named_clean) > 0) {
    # Route lipid-style names to LIPID MAPS; everything else to PubChem
    is_lipid  <- detect_lipid_names(named_clean)
    lipid_names   <- named_clean[is_lipid]
    general_names <- named_clean[!is_lipid]

    # Strategy 1a: LIPID MAPS for lipid shorthand names
    if (length(lipid_names) > 0) {
      log_info("  Strategy 1a: LIPID MAPS mapping (", length(lipid_names), " lipid features)")
      lm_map <- map_via_lipidmaps(lipid_names, cache_dir)
      for (i in seq_len(nrow(lm_map))) {
        ri <- match(lm_map$feature[i], results$clean_name)
        if (!is.na(ri)) {
          results$compound_name[ri]  <- lm_map$compound_name[i]
          results$hmdb_id[ri]        <- lm_map$hmdb_id[i]
          results$kegg_id[ri]        <- lm_map$kegg_id[i]
          results$lmid[ri]           <- lm_map$lmid[i]
          results$lipid_class[ri]    <- lm_map$lipid_class[i]
          results$lipid_subclass[ri] <- lm_map$lipid_subclass[i]
          results$mapping_method[ri] <- lm_map$mapping_method[i]
          results$confidence[ri]     <- lm_map$confidence[i]
        }
      }
    }

    # Strategy 1b: PubChem for all other named compounds
    if (length(general_names) > 0) {
      log_info("  Strategy 1b: PubChem mapping (", length(general_names), " features)")
      pc_map <- map_via_pubchem(general_names, cache_dir)
      for (i in seq_len(nrow(pc_map))) {
        ri <- match(pc_map$feature[i], results$clean_name)
        if (!is.na(ri) && is.na(results$hmdb_id[ri])) {
          results$compound_name[ri]  <- pc_map$compound_name[i]
          results$hmdb_id[ri]        <- pc_map$hmdb_id[i]
          results$kegg_id[ri]        <- pc_map$kegg_id[i]
          results$inchikey[ri]       <- pc_map$inchikey[i]
          results$mapping_method[ri] <- pc_map$mapping_method[i]
          results$confidence[ri]     <- pc_map$confidence[i]
        }
      }
    }

    # Strategy 1c: PubChem fallback for any lipids that LIPID MAPS missed
    lm_missed <- named_clean[is_lipid][is.na(results$hmdb_id[match(
      named_clean[is_lipid], results$clean_name)])]
    if (length(lm_missed) > 0) {
      log_info("  Strategy 1c: PubChem fallback for ", length(lm_missed),
               " unmatched lipids")
      pc_fallback <- map_via_pubchem(lm_missed, cache_dir)
      for (i in seq_len(nrow(pc_fallback))) {
        ri <- match(pc_fallback$feature[i], results$clean_name)
        if (!is.na(ri) && is.na(results$hmdb_id[ri])) {
          results$compound_name[ri]  <- pc_fallback$compound_name[i]
          results$hmdb_id[ri]        <- pc_fallback$hmdb_id[i]
          results$kegg_id[ri]        <- pc_fallback$kegg_id[i]
          results$inchikey[ri]       <- pc_fallback$inchikey[i]
          results$mapping_method[ri] <- pc_fallback$mapping_method[i]
          results$confidence[ri]     <- pc_fallback$confidence[i]
        }
      }
    }
  }

  # Strategy 2: m/z-based mapping for anything still unmapped
  if (!is.null(mz_rt_df)) {
    log_info("  Strategy 2: m/z-based mass matching")
    mz_map <- map_by_mz(mz_rt_df, config$ppm_tolerance, config$adducts)
    for (i in seq_len(nrow(mz_map))) {
      ri <- match(mz_map$feature[i], results$feature)
      if (!is.na(ri) && is.na(results$hmdb_id[ri])) {
        results$compound_name[ri]  <- mz_map$compound_name[i]
        results$hmdb_id[ri]        <- mz_map$hmdb_id[i]
        results$mapping_method[ri] <- "mz_match"
        results$confidence[ri]     <- mz_map$confidence[i]
      }
    }
  }

  results <- results %>%
    mutate(
      compound_name  = coalesce(compound_name, clean_name),
      mapping_method = replace_na(mapping_method, "none"),
      confidence     = replace_na(confidence, "none")
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

# ── Name cleaning ─────────────────────────────────────────────────────────────

clean_compound_names <- function(names) {
  names %>%
    str_remove_all("\\.mol$|\\.sdf$|\\.txt$") %>%  # strip file extensions
    str_squish() %>%                                # normalise whitespace
    str_trim()
}

# ── Lipid name detection ──────────────────────────────────────────────────────

detect_lipid_names <- function(names) {
  lipid_prefix <- paste0(
    "^(PC|PE|PI|PG|PS|PA|PT|PIP|",
    "LPC|LPE|LPI|LPG|LPS|LPA|",
    "TG|DG|MG|TAG|DAG|MAG|",
    "Cer|SM|HexCer|LacCer|GlcCer|GalCer|SHexCer|",
    "FA|CAR|AC|FAHFA|",
    "CE|FC|ChE|",
    "WE|SE|ST|",
    "CoA|AcCa|",
    "Hex[0-9]?Cer|",
    "GM[0-9]|GD[0-9]|GT[0-9])",
    "\\s*[\\(\\d\\s]"
  )
  # Also catch anything with fatty acid chain notation like (16:0/18:2)
  grepl(lipid_prefix, names, perl = TRUE) |
    grepl("\\(\\d{1,2}:\\d{1,2}[/\\\\,]\\d{1,2}:\\d{1,2}", names, perl = TRUE)
}

# ── Strategy 1a: LIPID MAPS API ───────────────────────────────────────────────

map_via_lipidmaps <- function(lipid_names, cache_dir) {
  cache_file <- file.path(cache_dir, "lipidmaps_cache.rds")
  cache      <- load_cache(cache_file)

  results <- map_dfr(seq_along(lipid_names), function(i) {
    name <- lipid_names[i]

    if (i %% 100 == 0)
      log_info("    LIPID MAPS progress: ", i, "/", length(lipid_names))

    if (!is.null(cache[[name]])) return(cache[[name]])

    result <- query_lipidmaps(name)
    cache[[name]] <<- result
    Sys.sleep(0.1)
    result
  })

  save_cache(cache, cache_file)
  results
}

query_lipidmaps <- function(name) {
  blank <- tibble(
    feature = name, compound_name = NA_character_,
    hmdb_id = NA_character_, kegg_id = NA_character_,
    lmid = NA_character_, lipid_class = NA_character_,
    lipid_subclass = NA_character_,
    mapping_method = "none", confidence = "none"
  )

  encoded <- utils::URLencode(name, repeated = TRUE)
  url     <- paste0("https://www.lipidmaps.org/rest/compound/name/",
                    encoded, "/all/json")

  resp <- tryCatch(
    httr::GET(url, httr::timeout(15)),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) return(blank)

  parsed <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                       simplifyVector = TRUE),
    error = function(e) NULL
  )
  if (is.null(parsed) || length(parsed) == 0 ||
      (is.data.frame(parsed) && nrow(parsed) == 0)) return(blank)

  # fromJSON may return a list or data frame depending on result count
  if (is.data.frame(parsed)) {
    hit <- parsed[1, ]
  } else {
    hit <- as.data.frame(as.list(parsed), stringsAsFactors = FALSE)
  }

  lmid    <- hit$lm_id    %||% NA_character_
  kegg_id <- hit$kegg_id  %||% NA_character_
  hmdb_id <- hit$hmdb_id  %||% NA_character_
  cname   <- hit$name     %||% name
  mclass  <- hit$main_class %||% NA_character_
  sclass  <- hit$sub_class  %||% NA_character_

  # Normalise HMDB ID format
  if (!is.na(hmdb_id) && nzchar(hmdb_id) && !grepl("^HMDB", hmdb_id))
    hmdb_id <- paste0("HMDB", str_pad(hmdb_id, 7, pad = "0"))

  has_id <- (!is.na(lmid) & nzchar(lmid)) |
            (!is.na(kegg_id) & nzchar(kegg_id)) |
            (!is.na(hmdb_id) & nzchar(hmdb_id))

  tibble(
    feature        = name,
    compound_name  = if_else(nzchar(cname %||% ""), cname, name),
    hmdb_id        = if_else(nzchar(hmdb_id %||% ""), hmdb_id, NA_character_),
    kegg_id        = if_else(nzchar(kegg_id %||% ""), kegg_id, NA_character_),
    lmid           = if_else(nzchar(lmid %||% ""), lmid, NA_character_),
    lipid_class    = mclass,
    lipid_subclass = sclass,
    mapping_method = if_else(has_id, "lipidmaps", "none"),
    confidence     = if_else(has_id, "high", "none")
  )
}

# ── Strategy 1b: PubChem API ─────────────────────────────────────────────────
# Batch POST requests → CID → synonyms (contains HMDB + KEGG IDs)

map_via_pubchem <- function(compound_names, cache_dir) {
  cache_file <- file.path(cache_dir, "pubchem_cache.rds")
  cache      <- load_cache(cache_file)

  uncached <- compound_names[!compound_names %in% names(cache)]

  if (length(uncached) > 0) {
    log_info("    Querying PubChem for ", length(uncached), " compounds ...")
    # Process in batches of 50
    batches <- split(uncached, ceiling(seq_along(uncached) / 50))
    for (b_idx in seq_along(batches)) {
      batch   <- batches[[b_idx]]
      log_info("    PubChem batch ", b_idx, "/", length(batches))
      cid_map <- pubchem_names_to_cids(batch)

      for (name in batch) {
        cid <- cid_map[[name]]
        if (is.null(cid) || is.na(cid)) {
          cache[[name]] <- pubchem_blank(name)
        } else {
          cache[[name]] <- pubchem_cid_to_ids(name, cid)
          Sys.sleep(0.15)
        }
      }
      Sys.sleep(0.5)
    }
    save_cache(cache, cache_file)
  }

  map_dfr(compound_names, function(n) cache[[n]] %||% pubchem_blank(n))
}

pubchem_names_to_cids <- function(names) {
  body_str <- paste(paste0("name=", utils::URLencode(names, repeated = TRUE)),
                    collapse = "&")
  resp <- tryCatch(
    httr::POST(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/cids/JSON",
      body    = body_str,
      encode  = "raw",
      httr::content_type("application/x-www-form-urlencoded"),
      httr::timeout(30)
    ),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200) {
    return(setNames(rep(list(NA), length(names)), names))
  }
  parsed <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(parsed)) return(setNames(rep(list(NA), length(names)), names))

  # Response: IdentifierList$CID (one CID per name, in order submitted)
  cids <- parsed$IdentifierList$CID
  if (is.null(cids)) return(setNames(rep(list(NA), length(names)), names))
  setNames(as.list(cids), names)
}

pubchem_cid_to_ids <- function(name, cid) {
  url  <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                 cid, "/synonyms/JSON")
  resp <- tryCatch(httr::GET(url, httr::timeout(15)), error = function(e) NULL)

  blank <- pubchem_blank(name)
  if (is.null(resp) || httr::status_code(resp) != 200) return(blank)

  parsed <- tryCatch(
    jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(parsed)) return(blank)

  syns <- parsed$InformationList$Information[[1]]$Synonym
  if (is.null(syns)) return(blank)

  # Extract HMDB ID (format: HMDB0000001 or HMDB00001)
  hmdb_hits <- syns[grepl("^HMDB\\d+$", syns, ignore.case = TRUE)]
  hmdb_id   <- if (length(hmdb_hits) > 0) hmdb_hits[1] else NA_character_

  # Extract KEGG ID (format: C##### or G#####)
  kegg_hits <- syns[grepl("^[CG]\\d{5}$", syns)]
  kegg_id   <- if (length(kegg_hits) > 0) kegg_hits[1] else NA_character_

  # Get InChIKey
  ikey_url  <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
                      cid, "/property/InChIKey,IUPACName/JSON")
  ikey_resp <- tryCatch(httr::GET(ikey_url, httr::timeout(10)), error = function(e) NULL)
  inchikey  <- NA_character_
  iupac     <- NA_character_
  if (!is.null(ikey_resp) && httr::status_code(ikey_resp) == 200) {
    ikey_parsed <- tryCatch(
      jsonlite::fromJSON(httr::content(ikey_resp, "text", encoding = "UTF-8")),
      error = function(e) NULL
    )
    if (!is.null(ikey_parsed)) {
      props    <- ikey_parsed$PropertyTable$Properties
      inchikey <- props$InChIKey[1] %||% NA_character_
      iupac    <- props$IUPACName[1] %||% NA_character_
    }
  }

  has_id <- !is.na(hmdb_id) | !is.na(kegg_id)

  tibble(
    feature        = name,
    compound_name  = syns[1] %||% name,
    hmdb_id        = hmdb_id,
    kegg_id        = kegg_id,
    lmid           = NA_character_,
    lipid_class    = NA_character_,
    lipid_subclass = NA_character_,
    inchikey       = inchikey,
    mapping_method = if_else(has_id, "pubchem", "pubchem_no_xref"),
    confidence     = case_when(
      !is.na(hmdb_id) & !is.na(kegg_id) ~ "high",
      has_id                              ~ "medium",
      !is.na(inchikey)                    ~ "low",
      TRUE                                ~ "none"
    )
  )
}

pubchem_blank <- function(name) {
  tibble(
    feature = name, compound_name = name,
    hmdb_id = NA_character_, kegg_id = NA_character_,
    lmid = NA_character_, lipid_class = NA_character_,
    lipid_subclass = NA_character_, inchikey = NA_character_,
    mapping_method = "none", confidence = "none"
  )
}

# ── Strategy 2: m/z-based mapping ────────────────────────────────────────────

map_by_mz <- function(mz_rt_df, ppm_tol, adducts) {
  if (!all(c("feature", "mz") %in% names(mz_rt_df))) {
    log_warn("mz_rt_df missing 'feature'/'mz' columns — skipping m/z mapping.")
    return(tibble(feature = character()))
  }

  ref_db <- get_hmdb_mass_reference()
  if (is.null(ref_db) || nrow(ref_db) == 0) {
    log_warn("Could not retrieve HMDB mass reference — skipping m/z mapping.")
    return(tibble(feature = character()))
  }

  all_adducts <- c(adducts$pos, adducts$neg)

  map_dfr(seq_len(nrow(mz_rt_df)), function(i) {
    obs_mz   <- mz_rt_df$mz[i]
    feat     <- mz_rt_df$feature[i]
    best_hit <- NULL

    for (adduct_mass in all_adducts) {
      neutral <- obs_mz - adduct_mass
      ppm_err <- abs((ref_db$monisotopic_molecular_weight - neutral) / neutral * 1e6)
      hits    <- ref_db[ppm_err <= ppm_tol, ]
      if (nrow(hits) > 0) {
        hits$ppm_error <- ppm_err[ppm_err <= ppm_tol]
        if (is.null(best_hit) || hits$ppm_error[1] < best_hit$ppm_error[1])
          best_hit <- hits[1, ]
      }
    }

    if (!is.null(best_hit)) {
      tibble(feature = feat, compound_name = best_hit$name,
             hmdb_id = best_hit$accession, kegg_id = NA_character_,
             confidence = if_else(best_hit$ppm_error < 2, "high", "medium"))
    } else {
      tibble(feature = feat, compound_name = NA_character_,
             hmdb_id = NA_character_, kegg_id = NA_character_,
             confidence = NA_character_)
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

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
