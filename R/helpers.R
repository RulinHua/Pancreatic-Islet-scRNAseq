# Small project-local helpers used by the analysis scripts.

project_palette <- function(n) {
  # Portable replacement for the author's local palette helper.
  # It preserves a stable categorical palette without requiring a private script.
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required for project_palette().")
  }
  scales::hue_pal()(n)
}
