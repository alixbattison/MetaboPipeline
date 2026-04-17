# MetaboPipeline — main entry point
#
# Usage:
#   Rscript run_pipeline.R
#
# Place your Excel files in the 'data/' directory (or change config$input_dir).
# Optionally supply an m/z + RT file per organ (see mz_rt_files map below).
# Results are written to 'results/{organ_name}/'.

# ── Locate the pipeline root reliably in both RStudio and Rscript ─────────────
if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable() &&
    nzchar(rstudioapi::getActiveDocumentContext()$path)) {
  .pipeline_root <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  .pipeline_root <- tryCatch(
    dirname(normalizePath(sys.frame(0)$ofile)),
    error = function(e) getwd()
  )
}
setwd(.pipeline_root)
message("Pipeline root: ", .pipeline_root)

suppressPackageStartupMessages(library(here))
here::i_am("run_pipeline.R")   # anchors here() to this file's directory

# ── Load configuration and all modules ───────────────────────────────────────
source(here("config.R"))
source(here("R/utils.R"))
source(here("R/01_load_data.R"))
source(here("R/02_id_mapping.R"))
source(here("R/03_qc.R"))
source(here("R/04_pca.R"))
source(here("R/05_differential.R"))
source(here("R/06_heatmap.R"))
source(here("R/07_pathway_enrichment.R"))

suppressPackageStartupMessages({
  library(purrr)
  library(dplyr)
})

set.seed(config$seed)

# ── Optional: m/z + RT data per file ─────────────────────────────────────────
# If you have exported m/z and RT from MetaboScape, create a CSV per organ with
# columns: feature, mz, rt
# Add entries here mapping the Excel base-name (without extension) to the CSV path.
# Example:
#   mz_rt_files <- list(
#     Plasma = "data/Plasma_mz_rt.csv",
#     Liver  = "data/Liver_mz_rt.csv"
#   )
mz_rt_files <- list()

# ── Discover input files ──────────────────────────────────────────────────────
excel_files <- list.files(
  path       = here(config$input_dir),
  pattern    = "\\.(xlsx|xls)$",
  full.names = TRUE
)

if (length(excel_files) == 0) {
  stop("No Excel files found in '", config$input_dir, "/'. ",
       "Please place your data files there and re-run.")
}

log_info("Found ", length(excel_files), " Excel file(s) to process:")
for (f in excel_files) log_info("  ", basename(f))

# ── Process each file ─────────────────────────────────────────────────────────
run_single_file <- function(file_path) {
  organ_name <- tools::file_path_sans_ext(basename(file_path))
  log_info("══════════════════════════════════════════")
  log_info("Processing: ", organ_name)
  log_info("══════════════════════════════════════════")

  organ_dir <- make_output_dirs(here(config$output_dir), organ_name)

  # Stage 1: Load data
  data_obj <- tryCatch(
    load_metabolomics_data(file_path, config),
    error = function(e) {
      log_error("Failed to load ", basename(file_path), ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(data_obj)) return(invisible(NULL))

  # Load m/z + RT if available for this organ
  mz_rt_df <- NULL
  if (organ_name %in% names(mz_rt_files)) {
    mz_rt_path <- mz_rt_files[[organ_name]]
    mz_rt_df   <- tryCatch(
      read.csv(mz_rt_path, stringsAsFactors = FALSE),
      error = function(e) {
        log_warn("Could not load m/z+RT file '", mz_rt_path, "': ", e$message)
        NULL
      }
    )
  }

  # Stage 2: ID mapping
  mapping_tbl <- tryCatch(
    run_id_mapping(data_obj, config, organ_dir, mz_rt_df = mz_rt_df),
    error = function(e) {
      log_warn("ID mapping failed: ", e$message, " — continuing without mapping.")
      tibble::tibble(feature = data_obj$feature_ids,
                     compound_name = NA_character_,
                     hmdb_id = NA_character_, kegg_id = NA_character_,
                     mapping_method = "none", confidence = "none")
    }
  )

  # Stage 3: QC
  tryCatch(
    run_qc(data_obj, config, organ_dir),
    error = function(e) log_warn("QC failed: ", e$message)
  )

  # Stage 4: PCA
  tryCatch(
    run_pca(data_obj, config, organ_dir),
    error = function(e) log_warn("PCA failed: ", e$message)
  )

  # Stage 5: Differential abundance
  diff_results <- tryCatch(
    run_differential(data_obj, config, organ_dir, mapping_tbl),
    error = function(e) {
      log_warn("Differential analysis failed: ", e$message)
      NULL
    }
  )

  # Stage 6: Heatmaps
  if (!is.null(diff_results)) {
    tryCatch(
      run_heatmap(data_obj, diff_results, config, organ_dir, mapping_tbl),
      error = function(e) log_warn("Heatmap generation failed: ", e$message)
    )

    # Stage 7: Pathway enrichment
    tryCatch(
      run_pathway_enrichment(data_obj, diff_results, config, organ_dir,
                             mapping_tbl, mz_rt_df),
      error = function(e) log_warn("Pathway enrichment failed: ", e$message)
    )
  }

  log_info("Finished: ", organ_name)
  log_info("Results written to: ", organ_dir)
  invisible(NULL)
}

# ── Run ───────────────────────────────────────────────────────────────────────
start_time <- proc.time()

for (f in excel_files) {
  run_single_file(f)
}

elapsed <- round((proc.time() - start_time)[["elapsed"]])
log_info("══════════════════════════════════════════")
log_info("Pipeline complete. Total time: ", elapsed, "s")
log_info("Results directory: ", here(config$output_dir))
