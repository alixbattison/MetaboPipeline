# Install all packages required by MetaboPipeline
# Run this script once before running the pipeline for the first time.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "stringr",
  "here",
  "ggplot2",
  "ggrepel",
  "RColorBrewer",
  "viridis",
  "patchwork",
  "factoextra",
  "FactoMineR",
  "dunn.test",
  "UpSetR",
  "httr",
  "jsonlite",
  "scales",
  "openxlsx",
  "pheatmap",
  "grid",
  "gridExtra"
)

bioc_packages <- c(
  "limma",
  "ComplexHeatmap",
  "circlize",
  "KEGGREST",
  "clusterProfiler",
  "ReactomePA",
  "rWikiPathways",
  "org.Mm.eg.db",
  "pathview",
  "AnnotationDbi"
)

missing_cran <- cran_packages[!cran_packages %in% installed.packages()[, "Package"]]
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, dependencies = TRUE)
}

missing_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[, "Package"]]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
}

message("All packages installed successfully.")
