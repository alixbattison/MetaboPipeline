suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

# ── Main pathway enrichment entry point ──────────────────────────────────────

run_pathway_enrichment <- function(data_obj, diff_results, config, organ_dir,
                                   mapping_tbl, mz_rt_df = NULL) {
  log_info("Running pathway enrichment for: ", data_obj$organ_name)
  pe_dir <- file.path(organ_dir, "pathway_enrichment")
  dir.create(pe_dir, recursive = TRUE, showWarnings = FALSE)

  # Gather all significantly differential features (Class A + B)
  sig_class_a <- purrr::map_dfr(diff_results$pairwise, function(df) {
    df %>% filter(significant) %>%
      select(feature, log2FC, adj_p, group1, group2, direction)
  }) %>% distinct(feature, .keep_all = TRUE)

  sig_class_b <- purrr::map_dfr(diff_results$class_b, function(df) {
    df %>% filter(significant) %>%
      select(feature, direction, group1, group2)
  }) %>% distinct(feature, .keep_all = TRUE)

  sig_features <- unique(c(sig_class_a$feature, sig_class_b$feature))
  background   <- data_obj$feature_ids

  if (length(sig_features) == 0) {
    log_warn("No significant features — skipping pathway enrichment.")
    return(invisible(NULL))
  }

  # Resolve KEGG / HMDB IDs from mapping table
  sig_mapped <- mapping_tbl %>%
    filter(feature %in% sig_features,
           confidence %in% c("high", "medium"))

  bg_mapped <- mapping_tbl %>%
    filter(feature %in% background,
           confidence %in% c("high", "medium"))

  all_results <- list()

  # ── KEGG ORA ─────────────────────────────────────────────────────────────
  kegg_dir <- file.path(pe_dir, "kegg")
  dir.create(kegg_dir, showWarnings = FALSE)
  kegg_res <- run_kegg_ora(sig_mapped, bg_mapped, config, kegg_dir, data_obj$organ_name)
  all_results$kegg <- kegg_res

  # ── Reactome ORA ─────────────────────────────────────────────────────────
  react_dir <- file.path(pe_dir, "reactome")
  dir.create(react_dir, showWarnings = FALSE)
  react_res <- run_reactome_ora(sig_mapped, bg_mapped, config, react_dir, data_obj$organ_name)
  all_results$reactome <- react_res

  # ── WikiPathways ORA ──────────────────────────────────────────────────────
  wiki_dir <- file.path(pe_dir, "wikipathways")
  dir.create(wiki_dir, showWarnings = FALSE)
  wiki_res <- run_wikipathways_ora(sig_mapped, bg_mapped, config, wiki_dir, data_obj$organ_name)
  all_results$wikipathways <- wiki_res

  # ── HMDB pathway ORA ─────────────────────────────────────────────────────
  hmdb_dir <- file.path(pe_dir, "hmdb")
  dir.create(hmdb_dir, showWarnings = FALSE)
  hmdb_res <- run_hmdb_ora(sig_mapped, bg_mapped, config, hmdb_dir, data_obj$organ_name)
  all_results$hmdb <- hmdb_res

  # ── LIPID MAPS enrichment ─────────────────────────────────────────────────
  lm_dir <- file.path(pe_dir, "lipidmaps")
  dir.create(lm_dir, showWarnings = FALSE)
  lm_res <- run_lipidmaps_enrichment(sig_mapped, bg_mapped, config, lm_dir,
                                      data_obj$organ_name)
  all_results$lipidmaps <- lm_res

  # ── mummichog-style m/z enrichment (if m/z provided) ─────────────────────
  if (!is.null(mz_rt_df)) {
    mz_dir <- file.path(pe_dir, "mummichog")
    dir.create(mz_dir, showWarnings = FALSE)
    sig_mz_df <- mz_rt_df %>% filter(feature %in% sig_features)
    bg_mz_df  <- mz_rt_df
    mz_res <- run_mummichog(sig_mz_df, bg_mz_df, config, mz_dir, data_obj$organ_name)
    all_results$mummichog <- mz_res
  }

  # ── Cross-database summary ────────────────────────────────────────────────
  cross_db <- build_cross_db_summary(all_results)
  if (!is.null(cross_db) && nrow(cross_db) > 0) {
    utils::write.csv(cross_db,
                     file.path(pe_dir, "cross_database_summary.csv"),
                     row.names = FALSE)
    p_cross <- plot_cross_db(cross_db, data_obj$organ_name)
    if (!is.null(p_cross)) {
      save_plot(p_cross, file.path(pe_dir, "cross_database_dotplot.pdf"),
                width = 12, height = 8)
    }
  }

  log_info("  Pathway enrichment complete")
  invisible(all_results)
}

# ── KEGG ORA ─────────────────────────────────────────────────────────────────

run_kegg_ora <- function(sig_mapped, bg_mapped, config, out_dir, organ_name) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    log_warn("clusterProfiler or org.Mm.eg.db not installed — skipping KEGG ORA.")
    return(NULL)
  }

  kegg_sig <- sig_mapped %>% filter(!is.na(kegg_id)) %>% pull(kegg_id)
  kegg_bg  <- bg_mapped  %>% filter(!is.na(kegg_id)) %>% pull(kegg_id)

  if (length(kegg_sig) < 1) {
    log_warn("Fewer than 1 KEGG-mapped significant features — skipping KEGG ORA.")
    return(NULL)
  }

  res <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = kegg_sig,
      universe      = kegg_bg,
      organism      = config$organism_kegg,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      minGSSize     = config$min_pathway_size,
      maxGSSize     = config$max_pathway_size
    ),
    error = function(e) {
      log_warn("KEGG ORA failed: ", e$message)
      NULL
    }
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    log_info("  KEGG: no enriched pathways found.")
    return(NULL)
  }

  res_df <- as.data.frame(res)
  utils::write.csv(res_df, file.path(out_dir, "kegg_ora_results.csv"), row.names = FALSE)

  p <- plot_enrichment_dotplot(res_df, paste("KEGG —", organ_name), "KEGG")
  save_plot(p, file.path(out_dir, "kegg_dotplot.pdf"), width = 10, height = 7)

  # KEGG pathway diagrams for top 5
  if (requireNamespace("pathview", quietly = TRUE)) {
    top_pathways <- head(res_df$ID, 5)
    sig_fc <- sig_mapped %>%
      filter(!is.na(kegg_id)) %>%
      select(kegg_id) %>%
      mutate(fc = 1) %>%
      deframe()
    for (pw in top_pathways) {
      tryCatch(
        pathview::pathview(
          gene.data  = sig_fc,
          pathway.id = sub("mmu", "", pw),
          species    = config$organism_kegg,
          out.suffix = gsub("[^A-Za-z0-9]", "_", organ_name),
          kegg.dir   = out_dir
        ),
        error = function(e) log_warn("pathview failed for ", pw, ": ", e$message)
      )
    }
  }

  res_df %>% mutate(database = "KEGG")
}

# ── Reactome ORA ─────────────────────────────────────────────────────────────

run_reactome_ora <- function(sig_mapped, bg_mapped, config, out_dir, organ_name) {
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    log_warn("ReactomePA not installed — skipping Reactome ORA.")
    return(NULL)
  }

  # ReactomePA works with Entrez IDs; attempt conversion from KEGG IDs
  kegg_sig <- sig_mapped %>% filter(!is.na(kegg_id)) %>% pull(kegg_id)
  kegg_bg  <- bg_mapped  %>% filter(!is.na(kegg_id)) %>% pull(kegg_id)

  entrez_sig <- kegg_to_entrez(kegg_sig)
  entrez_bg  <- kegg_to_entrez(kegg_bg)

  if (length(entrez_sig) < 1) {
    log_warn("Fewer than 1 Entrez-mapped features — skipping Reactome ORA.")
    return(NULL)
  }

  res <- tryCatch(
    ReactomePA::enrichPathway(
      gene          = entrez_sig,
      universe      = entrez_bg,
      organism      = config$organism_common,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      minGSSize     = config$min_pathway_size,
      maxGSSize     = config$max_pathway_size
    ),
    error = function(e) {
      log_warn("Reactome ORA failed: ", e$message)
      NULL
    }
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    log_info("  Reactome: no enriched pathways found.")
    return(NULL)
  }

  res_df <- as.data.frame(res)
  utils::write.csv(res_df, file.path(out_dir, "reactome_ora_results.csv"), row.names = FALSE)
  p <- plot_enrichment_dotplot(res_df, paste("Reactome —", organ_name), "Reactome")
  save_plot(p, file.path(out_dir, "reactome_dotplot.pdf"), width = 10, height = 7)

  res_df %>% mutate(database = "Reactome")
}

# ── WikiPathways ORA ──────────────────────────────────────────────────────────

run_wikipathways_ora <- function(sig_mapped, bg_mapped, config, out_dir, organ_name) {
  if (!requireNamespace("rWikiPathways", quietly = TRUE) ||
      !requireNamespace("clusterProfiler", quietly = TRUE)) {
    log_warn("rWikiPathways not installed — skipping WikiPathways ORA.")
    return(NULL)
  }

  kegg_sig <- sig_mapped %>% filter(!is.na(kegg_id)) %>% pull(kegg_id)
  kegg_bg  <- bg_mapped  %>% filter(!is.na(kegg_id)) %>% pull(kegg_id)
  entrez_sig <- kegg_to_entrez(kegg_sig)
  entrez_bg  <- kegg_to_entrez(kegg_bg)

  if (length(entrez_sig) < 1) {
    log_warn("Fewer than 1 Entrez-mapped features — skipping WikiPathways ORA.")
    return(NULL)
  }

  wp_gmt <- tryCatch(
    rWikiPathways::downloadPathwayArchive(
      organism  = "Mus musculus",
      format    = "gmt",
      destpath  = tempdir()
    ),
    error = function(e) { log_warn("WikiPathways download failed: ", e$message); NULL }
  )
  if (is.null(wp_gmt)) return(NULL)

  wp_df <- tryCatch(
    clusterProfiler::read.gmt(wp_gmt),
    error = function(e) NULL
  )
  if (is.null(wp_df)) return(NULL)

  res <- tryCatch(
    clusterProfiler::enricher(
      gene          = entrez_sig,
      universe      = entrez_bg,
      TERM2GENE     = wp_df,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      minGSSize     = config$min_pathway_size,
      maxGSSize     = config$max_pathway_size
    ),
    error = function(e) {
      log_warn("WikiPathways ORA failed: ", e$message)
      NULL
    }
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    log_info("  WikiPathways: no enriched pathways found.")
    return(NULL)
  }

  res_df <- as.data.frame(res)
  utils::write.csv(res_df, file.path(out_dir, "wikipathways_ora_results.csv"), row.names = FALSE)
  p <- plot_enrichment_dotplot(res_df, paste("WikiPathways —", organ_name), "WikiPathways")
  save_plot(p, file.path(out_dir, "wikipathways_dotplot.pdf"), width = 10, height = 7)

  res_df %>% mutate(database = "WikiPathways")
}

# ── HMDB pathway ORA ─────────────────────────────────────────────────────────

run_hmdb_ora <- function(sig_mapped, bg_mapped, config, out_dir, organ_name) {
  hmdb_sig <- sig_mapped %>% filter(!is.na(hmdb_id)) %>% pull(hmdb_id)
  hmdb_bg  <- bg_mapped  %>% filter(!is.na(hmdb_id)) %>% pull(hmdb_id)

  if (length(hmdb_sig) < 1) {
    log_warn("Fewer than 1 HMDB-mapped features — skipping HMDB pathway ORA.")
    return(NULL)
  }

  pathways <- fetch_hmdb_pathways(hmdb_bg)
  if (is.null(pathways) || nrow(pathways) == 0) {
    log_warn("Could not retrieve HMDB pathway data.")
    return(NULL)
  }

  # Fisher's exact test per pathway
  res <- purrr::map_dfr(unique(pathways$pathway_name), function(pw) {
    pw_members  <- pathways$hmdb_id[pathways$pathway_name == pw]
    in_sig_in   <- sum(hmdb_sig %in% pw_members)
    in_sig_out  <- length(hmdb_sig) - in_sig_in
    out_sig_in  <- sum(hmdb_bg[!hmdb_bg %in% hmdb_sig] %in% pw_members)
    out_sig_out <- length(hmdb_bg) - length(hmdb_sig) - out_sig_in
    ct  <- matrix(c(in_sig_in, in_sig_out, out_sig_in, out_sig_out), nrow = 2)
    pval <- tryCatch(fisher.test(ct, alternative = "greater")$p.value,
                     error = function(e) NA_real_)
    tibble(
      pathway      = pw,
      n_sig        = in_sig_in,
      n_pathway    = length(pw_members),
      p_value      = pval
    )
  }) %>%
    filter(!is.na(p_value), n_sig > 0) %>%
    mutate(p_adjust = p.adjust(p_value, method = "BH"),
           GeneRatio = paste0(n_sig, "/", length(hmdb_sig)),
           BgRatio   = paste0(n_pathway, "/", length(hmdb_bg))) %>%
    filter(p_adjust < 0.05) %>%
    arrange(p_adjust)

  if (nrow(res) == 0) {
    log_info("  HMDB: no enriched pathways found.")
    return(NULL)
  }

  utils::write.csv(res, file.path(out_dir, "hmdb_ora_results.csv"), row.names = FALSE)
  p <- plot_enrichment_dotplot(
    res %>% rename(Description = pathway, p.adjust = p_adjust, Count = n_sig),
    paste("HMDB —", organ_name), "HMDB"
  )
  save_plot(p, file.path(out_dir, "hmdb_dotplot.pdf"), width = 10, height = 7)

  res %>% mutate(database = "HMDB")
}

fetch_hmdb_pathways <- function(hmdb_ids) {
  # Query HMDB REST API for pathway associations — no cap, all mapped IDs used.
  # Results are cached to avoid re-querying on subsequent runs.
  cache_file <- here::here("data", ".mapping_cache", "hmdb_pathways_cache.rds")
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
  cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()

  unique_ids  <- unique(na.omit(hmdb_ids))
  uncached    <- unique_ids[!unique_ids %in% names(cache)]

  if (length(uncached) > 0) {
    log_info("    Fetching HMDB pathways for ", length(uncached), " compounds ...")
    for (i in seq_along(uncached)) {
      hid  <- uncached[i]
      if (i %% 50 == 0)
        log_info("    HMDB pathway progress: ", i, "/", length(uncached))
      url  <- paste0("https://hmdb.ca/metabolites/", hid, ".json")
      resp <- tryCatch(httr::GET(url, httr::timeout(10)), error = function(e) NULL)
      if (is.null(resp) || httr::status_code(resp) != 200) {
        cache[[hid]] <- tibble(hmdb_id = hid, pathway_name = character(0))
        next
      }
      parsed <- tryCatch(
        jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                           simplifyVector = FALSE),
        error = function(e) NULL
      )
      pathways <- if (!is.null(parsed)) parsed$pathways else NULL
      if (is.null(pathways) || length(pathways) == 0) {
        cache[[hid]] <- tibble(hmdb_id = hid, pathway_name = character(0))
      } else {
        cache[[hid]] <- purrr::map_dfr(pathways, function(pw) {
          tibble(hmdb_id = hid, pathway_name = pw$name %||% NA_character_)
        })
      }
      Sys.sleep(0.1)
    }
    saveRDS(cache, cache_file)
  }

  purrr::map_dfr(unique_ids, function(hid) cache[[hid]])
}

# ── mummichog-style m/z enrichment ───────────────────────────────────────────

run_mummichog <- function(sig_mz_df, bg_mz_df, config, out_dir, organ_name) {
  if (!all(c("feature", "mz") %in% names(sig_mz_df))) {
    log_warn("mz_rt_df missing required columns. Skipping mummichog.")
    return(NULL)
  }

  log_info("  Running mummichog-style m/z enrichment")

  # Fetch KEGG compound masses as reference
  ref_masses <- fetch_kegg_compound_masses(config$organism_kegg)
  if (is.null(ref_masses) || nrow(ref_masses) == 0) {
    log_warn("Could not fetch KEGG compound masses for mummichog.")
    return(NULL)
  }

  all_adducts <- c(config$adducts$pos, config$adducts$neg)

  # Match each observed m/z to a set of possible compound IDs
  match_mz_to_compounds <- function(mz_vec) {
    purrr::map(mz_vec, function(obs_mz) {
      compound_hits <- character(0)
      for (adduct_mass in all_adducts) {
        neutral <- obs_mz - adduct_mass
        ppm_err <- abs((ref_masses$mass - neutral) / neutral * 1e6)
        hits    <- ref_masses$compound_id[ppm_err <= config$ppm_tolerance]
        compound_hits <- c(compound_hits, hits)
      }
      unique(compound_hits)
    })
  }

  sig_hits <- match_mz_to_compounds(sig_mz_df$mz)
  bg_hits  <- match_mz_to_compounds(bg_mz_df$mz)

  sig_compounds <- unique(unlist(sig_hits))
  bg_compounds  <- unique(unlist(bg_hits))

  if (length(sig_compounds) < 1) {
    log_warn("Fewer than 1 m/z-matched compounds — skipping mummichog.")
    return(NULL)
  }

  # Fetch KEGG pathways and run Fisher ORA
  kegg_pathways <- fetch_kegg_pathways(config$organism_kegg)
  if (is.null(kegg_pathways)) return(NULL)

  res <- purrr::map_dfr(names(kegg_pathways), function(pw_id) {
    pw_members  <- kegg_pathways[[pw_id]]$compounds
    in_sig_in   <- sum(sig_compounds %in% pw_members)
    if (in_sig_in == 0) return(NULL)
    in_sig_out  <- length(sig_compounds) - in_sig_in
    out_sig_in  <- sum(bg_compounds[!bg_compounds %in% sig_compounds] %in% pw_members)
    out_sig_out <- length(bg_compounds) - length(sig_compounds) - out_sig_in
    ct   <- matrix(c(in_sig_in, in_sig_out, out_sig_in, out_sig_out), nrow = 2)
    pval <- tryCatch(fisher.test(ct, alternative = "greater")$p.value,
                     error = function(e) NA_real_)
    tibble(pathway_id   = pw_id,
           pathway_name = kegg_pathways[[pw_id]]$name,
           n_sig        = in_sig_in,
           n_pathway    = length(pw_members),
           p_value      = pval)
  }) %>%
    filter(!is.na(p_value)) %>%
    mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
    filter(p_adjust < 0.05) %>%
    arrange(p_adjust)

  if (nrow(res) == 0) {
    log_info("  mummichog: no enriched pathways found.")
    return(NULL)
  }

  utils::write.csv(res, file.path(out_dir, "mummichog_results.csv"), row.names = FALSE)
  p <- plot_enrichment_dotplot(
    res %>% rename(Description = pathway_name, p.adjust = p_adjust, Count = n_sig),
    paste("mummichog (m/z-based) —", organ_name), "mummichog"
  )
  save_plot(p, file.path(out_dir, "mummichog_dotplot.pdf"), width = 10, height = 7)

  res %>% mutate(database = "mummichog")
}

# ── Helpers ───────────────────────────────────────────────────────────────────

kegg_to_entrez <- function(kegg_ids) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(character(0))
  unique_ids <- unique(na.omit(kegg_ids))
  if (length(unique_ids) == 0) return(character(0))

  entrez <- purrr::map_chr(unique_ids, function(kid) {
    info <- tryCatch(KEGGREST::keggGet(kid)[[1]], error = function(e) NULL)
    if (is.null(info)) return(NA_character_)
    dblinks <- info$DBLINKS
    entrez_entry <- grep("^NCBI-GeneID", dblinks, value = TRUE)
    if (length(entrez_entry) == 0) return(NA_character_)
    sub("NCBI-GeneID: ", "", entrez_entry[1])
  })
  na.omit(entrez)
}

fetch_kegg_compound_masses <- function(organism) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  tryCatch({
    all_cpds <- KEGGREST::keggList("compound")
    cpd_ids  <- names(all_cpds)[seq_len(min(500, length(all_cpds)))]
    purrr::map_dfr(cpd_ids, function(cid) {
      info <- tryCatch(KEGGREST::keggGet(cid)[[1]], error = function(e) NULL)
      if (is.null(info)) return(NULL)
      mass <- as.numeric(info$EXACT_MASS)
      if (is.na(mass)) return(NULL)
      tibble(compound_id = cid, mass = mass)
    })
  }, error = function(e) {
    log_warn("Could not fetch KEGG compound masses: ", e$message)
    NULL
  })
}

fetch_kegg_pathways <- function(organism) {
  if (!requireNamespace("KEGGREST", quietly = TRUE)) return(NULL)
  tryCatch({
    pw_list <- KEGGREST::keggList("pathway", organism)
    purrr::map(names(pw_list), function(pw_id) {
      info <- tryCatch(KEGGREST::keggGet(pw_id)[[1]], error = function(e) NULL)
      if (is.null(info)) return(NULL)
      list(name = info$NAME, compounds = names(info$COMPOUND))
    }) %>% setNames(names(pw_list))
  }, error = function(e) {
    log_warn("Could not fetch KEGG pathways: ", e$message)
    NULL
  })
}

# ── Dot plot ──────────────────────────────────────────────────────────────────

plot_enrichment_dotplot <- function(res_df, title, db_name, top_n = 20) {
  if (!"Description" %in% names(res_df)) return(NULL)
  if (!"p.adjust" %in% names(res_df))   return(NULL)
  if (!"Count" %in% names(res_df) && "n_sig" %in% names(res_df)) {
    res_df$Count <- res_df$n_sig
  }
  if (!"Count" %in% names(res_df)) return(NULL)

  plot_df <- res_df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Description  = stringr::str_wrap(Description, 40),
      Description  = factor(Description, levels = rev(Description)),
      neg_log10_fdr = -log10(p.adjust)
    )

  ggplot(plot_df, aes(x = neg_log10_fdr, y = Description, size = Count,
                      colour = neg_log10_fdr)) +
    geom_point() +
    scale_colour_gradient(low = "#fee08b", high = "#d73027",
                          name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Feature count", range = c(2, 8)) +
    labs(title = title,
         x     = expression(-log[10](adjusted~p-value)),
         y     = NULL) +
    theme_metabo() +
    theme(axis.text.y = element_text(size = 8))
}

# ── Cross-database summary ────────────────────────────────────────────────────

build_cross_db_summary <- function(all_results) {
  non_null <- Filter(Negate(is.null), all_results)
  if (length(non_null) == 0) return(NULL)

  dbs_with_desc <- purrr::keep(non_null, ~ "Description" %in% names(.))
  if (length(dbs_with_desc) == 0) return(NULL)

  combined <- purrr::map_dfr(dbs_with_desc, function(df) {
    df %>% select(any_of(c("Description", "p.adjust", "Count", "database",
                            "pathway", "pathway_name"))) %>%
      rename(any_of(c(Description = "pathway", Description = "pathway_name")))
  })

  combined %>%
    group_by(Description) %>%
    summarise(
      databases       = paste(sort(unique(database)), collapse = "; "),
      n_databases     = n_distinct(database),
      min_fdr         = min(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_databases), min_fdr)
}

plot_cross_db <- function(cross_db, organ_name, top_n = 30) {
  if (nrow(cross_db) == 0) return(NULL)

  plot_df <- cross_db %>%
    slice_head(n = top_n) %>%
    mutate(
      Description = stringr::str_wrap(Description, 45),
      Description = factor(Description, levels = rev(Description))
    )

  ggplot(plot_df, aes(x = Description, y = n_databases, fill = -log10(min_fdr))) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "#fee08b", high = "#d73027",
                        name = expression(-log[10](best~FDR))) +
    labs(title = paste("Cross-database pathway overlap —", organ_name),
         x = NULL, y = "Number of databases with enrichment") +
    theme_metabo()
}

# ── LIPID MAPS enrichment ─────────────────────────────────────────────────────
# Two-level enrichment:
#   1. Lipid class level  (e.g. Glycerophospholipids vs Sphingolipids)
#   2. Lipid sub-class level (e.g. Phosphatidylcholines vs Lysophosphatidylcholines)

run_lipidmaps_enrichment <- function(sig_mapped, bg_mapped, config,
                                      out_dir, organ_name) {
  # Guard: columns may be absent if mapping failed entirely
  if (!all(c("lipid_class", "lmid") %in% names(sig_mapped))) {
    log_info("  LIPID MAPS: mapping columns absent — skipping.")
    return(NULL)
  }

  sig_lm <- sig_mapped %>% filter(!is.na(lipid_class) | !is.na(lmid))
  bg_lm  <- bg_mapped  %>% filter(!is.na(lipid_class) | !is.na(lmid))

  if (nrow(sig_lm) < 3) {
    log_info("  LIPID MAPS: fewer than 3 lipid features mapped — skipping.")
    return(NULL)
  }

  log_info("  LIPID MAPS: ", nrow(sig_lm), " significant lipid features, ",
           nrow(bg_lm), " background lipid features")

  results <- list()

  # ── Class-level ORA ───────────────────────────────────────────────────────
  class_res <- lipidmaps_ora(
    sig_classes = sig_lm$lipid_class,
    bg_classes  = bg_lm$lipid_class,
    level_name  = "Lipid class",
    n_sig_total = nrow(sig_mapped),
    n_bg_total  = nrow(bg_mapped)
  )
  if (!is.null(class_res) && nrow(class_res) > 0) {
    utils::write.csv(class_res,
                     file.path(out_dir, "lipidmaps_class_ora.csv"),
                     row.names = FALSE)
    p_class <- plot_lipidmaps(class_res,
                               paste("LIPID MAPS class —", organ_name))
    save_plot(p_class, file.path(out_dir, "lipidmaps_class_dotplot.pdf"),
              width = 10, height = 6)
    results$class <- class_res
  }

  # ── Sub-class-level ORA ───────────────────────────────────────────────────
  sub_res <- lipidmaps_ora(
    sig_classes = sig_lm$lipid_subclass,
    bg_classes  = bg_lm$lipid_subclass,
    level_name  = "Lipid sub-class",
    n_sig_total = nrow(sig_mapped),
    n_bg_total  = nrow(bg_mapped)
  )
  if (!is.null(sub_res) && nrow(sub_res) > 0) {
    utils::write.csv(sub_res,
                     file.path(out_dir, "lipidmaps_subclass_ora.csv"),
                     row.names = FALSE)
    p_sub <- plot_lipidmaps(sub_res,
                             paste("LIPID MAPS sub-class —", organ_name))
    save_plot(p_sub, file.path(out_dir, "lipidmaps_subclass_dotplot.pdf"),
              width = 10, height = 7)
    results$subclass <- sub_res
  }

  # ── Lipid class composition bar chart ─────────────────────────────────────
  if (nrow(sig_lm) > 0 && any(!is.na(sig_lm$lipid_class))) {
    p_comp <- plot_lipid_composition(sig_lm, bg_lm, organ_name)
    save_plot(p_comp, file.path(out_dir, "lipid_class_composition.pdf"),
              width = 10, height = 6)
  }

  if (length(results) == 0) return(NULL)

  combined <- purrr::map_dfr(results, function(df) {
    df %>% mutate(database = "LIPID MAPS",
                  Description = category,
                  p.adjust    = fdr,
                  Count       = n_sig)
  })
  invisible(combined)
}

lipidmaps_ora <- function(sig_classes, bg_classes, level_name,
                           n_sig_total, n_bg_total) {
  sig_classes <- na.omit(sig_classes)
  bg_classes  <- na.omit(bg_classes)
  if (length(sig_classes) < 2) return(NULL)

  all_cats <- unique(c(sig_classes, bg_classes))

  res <- purrr::map_dfr(all_cats, function(cat) {
    in_sig_in   <- sum(sig_classes == cat)
    if (in_sig_in == 0) return(NULL)
    in_sig_out  <- length(sig_classes) - in_sig_in
    out_sig_in  <- sum(bg_classes[!bg_classes %in% sig_classes] == cat)
    out_sig_out <- length(bg_classes) - length(sig_classes) - out_sig_in
    ct   <- matrix(c(in_sig_in, in_sig_out, out_sig_in, out_sig_out), nrow = 2)
    pval <- tryCatch(fisher.test(ct, alternative = "greater")$p.value,
                     error = function(e) NA_real_)
    tibble(
      category   = cat,
      level      = level_name,
      n_sig      = in_sig_in,
      n_bg       = sum(bg_classes == cat),
      pct_sig    = round(100 * in_sig_in / n_sig_total, 1),
      p_value    = pval
    )
  }) %>%
    filter(!is.na(p_value)) %>%
    mutate(fdr = p.adjust(p_value, method = "BH")) %>%
    filter(fdr < 0.05) %>%
    arrange(fdr)

  if (nrow(res) == 0) return(NULL)
  res
}

plot_lipidmaps <- function(res, title) {
  plot_df <- res %>%
    mutate(
      category     = stringr::str_wrap(category, 35),
      category     = factor(category, levels = rev(category)),
      neg_log10_fdr = -log10(fdr)
    )

  ggplot(plot_df, aes(x = neg_log10_fdr, y = category,
                      size = n_sig, colour = neg_log10_fdr)) +
    geom_point() +
    scale_colour_gradient(low = "#fee08b", high = "#d73027",
                          name = expression(-log[10](FDR))) +
    scale_size_continuous(name = "Significant\nfeatures", range = c(3, 9)) +
    labs(title = title,
         x     = expression(-log[10](adjusted~p-value)),
         y     = NULL) +
    theme_metabo() +
    theme(axis.text.y = element_text(size = 9))
}

plot_lipid_composition <- function(sig_lm, bg_lm, organ_name) {
  comp <- bind_rows(
    sig_lm %>% filter(!is.na(lipid_class)) %>%
      count(lipid_class) %>% mutate(group = "Significant"),
    bg_lm %>% filter(!is.na(lipid_class)) %>%
      count(lipid_class) %>% mutate(group = "Background")
  ) %>%
    group_by(group) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  ggplot(comp, aes(x = group, y = prop, fill = lipid_class)) +
    geom_col(position = "stack", colour = "white", linewidth = 0.3) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_brewer(palette = "Set3", name = "Lipid class") +
    labs(title = paste("Lipid class composition —", organ_name),
         x = NULL, y = "Proportion of features") +
    theme_metabo()
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
