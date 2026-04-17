# Pipeline configuration — edit these values to match your experiment

config <- list(

  # ── Paths ────────────────────────────────────────────────────────────────────
  input_dir  = "data",        # Directory containing Excel files
  output_dir = "results",     # Root output directory

  # ── Experimental groups ───────────────────────────────────────────────────────
  # Ordered vector of group labels in priority order (longer/more-specific first
  # so the regex does not match "C" inside "NC" or "TC").
  groups = c("NC", "NN", "TN", "TC", "C"),

  # Human-readable labels for plots (must match groups above, same order)
  group_labels = c(
    NC = "No Torpor Cancer",
    NN = "No Torpor No Cancer",
    TN = "Torpor No Cancer",
    TC = "Torpor Cancer",
    C  = "Control"
  ),

  # Colours assigned to each group for all plots
  group_colours = c(
    NC = "#E41A1C",
    NN = "#377EB8",
    TN = "#4DAF4A",
    TC = "#FF7F00",
    C  = "#984EA3"
  ),

  # ── Presence / absence thresholds ────────────────────────────────────────────
  # A feature is considered "present" in a group if the proportion of non-zero
  # replicates exceeds this value (strictly greater than).
  presence_threshold = 0.5,

  # ── Differential abundance ────────────────────────────────────────────────────
  fdr_threshold = 0.05,       # Adjusted p-value cut-off
  fc_threshold  = 1.0,        # Absolute log2 fold-change cut-off for Class A

  # ── Heatmap ──────────────────────────────────────────────────────────────────
  heatmap_top_n = 50,         # Max features to show per heatmap

  # ── ID mapping ───────────────────────────────────────────────────────────────
  ppm_tolerance = 5,          # Mass tolerance for m/z-based mapping (ppm)
  # Adducts considered in positive and negative mode
  adducts = list(
    pos = c("[M+H]+"  = 1.007276,
            "[M+Na]+" = 22.989218,
            "[M+K]+"  = 38.963158),
    neg = c("[M-H]-"  = -1.007276,
            "[M+FA-H]-" = 44.997655)
  ),

  # ── Pathway enrichment ───────────────────────────────────────────────────────
  organism_kegg   = "mmu",    # KEGG organism code for mouse
  organism_common = "mouse",
  min_pathway_size = 3,       # Minimum pathway gene/metabolite set size
  max_pathway_size = 500,

  # ── Reproducibility ──────────────────────────────────────────────────────────
  seed = 42
)
