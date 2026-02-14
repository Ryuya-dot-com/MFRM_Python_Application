#' mfrmr: Flexible Many-Facet Rasch Modeling in R
#'
#' @description
#' `mfrmr` provides estimation, diagnostics, and reporting utilities for
#' many-facet Rasch models (MFRM) without relying on FACETS/TAM/sirt backends.
#'
#' @details
#' Recommended workflow:
#'
#' 1. Fit model with [fit_mfrm()]
#' 2. Compute diagnostics with [diagnose_mfrm()]
#' 3. Run residual PCA with [analyze_residual_pca()] if needed
#' 4. Estimate interactions with [estimate_bias()]
#' 5. Build narrative/report outputs with [build_apa_outputs()] and [build_visual_summaries()]
#'
#' Data interface:
#' - Input analysis data is long format (one row per observed rating).
#' - Packaged simulation data is available via [load_mfrmr_data()] or `data()`.
#'
#' @name mfrmr-package
"_PACKAGE"
