#' mfrmr Workflow and Method Map
#'
#' @description
#' Quick reference for end-to-end `mfrmr` analysis and for checking which
#' output objects support `summary()` and `plot()`.
#'
#' @section Typical workflow:
#' 1. Fit a model with [fit_mfrm()].
#' 2. (Optional) Use [run_mfrm_facets()] or [mfrmRFacets()] for a
#'    FACETS-style one-shot workflow wrapper.
#' 3. Build diagnostics with [diagnose_mfrm()].
#' 4. (Optional) Estimate interaction bias with [estimate_bias()].
#' 5. Generate reporting bundles:
#'    [apa_table()], [build_fixed_reports()], [build_visual_summaries()].
#' 6. (Optional) Audit FACETS parity with `facets_parity_report()`.
#' 7. Use `summary()` for compact text checks and `plot()` (or dedicated plot
#'    helpers) for base-R visual diagnostics.
#'
#' @section Interpreting output:
#' This help page is a map, not an estimator:
#' - use it to decide function order,
#' - confirm which objects have `summary()`/`plot()` defaults,
#' - and identify when dedicated helper functions are needed.
#'
#' @section Objects with default `summary()` and `plot()` routes:
#' - `mfrm_fit`: `summary(fit)` and `plot(fit, ...)`.
#' - `mfrm_diagnostics`: `summary(diag)`; plotting via dedicated helpers
#'   such as [plot_unexpected()], [plot_displacement()], [plot_qc_dashboard()].
#' - `mfrm_bias`: `summary(bias)` and [plot_bias_interaction()].
#' - `mfrm_data_description`: `summary(ds)` and `plot(ds, ...)`.
#' - `mfrm_anchor_audit`: `summary(aud)` and `plot(aud, ...)`.
#' - `mfrm_facets_run`: `summary(run)` and `plot(run, type = c("fit", "qc"), ...)`.
#' - `apa_table`: `summary(tbl)` and `plot(tbl, ...)`.
#' - `mfrm_apa_outputs`: `summary(apa)` for compact diagnostics of report text.
#' - `mfrm_threshold_profiles`: `summary(profiles)` for preset threshold grids.
#' - `mfrm_bundle` families:
#'   `summary()` and class-aware `plot(bundle, ...)`.
#'   Key bundle classes now also use class-aware `summary(bundle)`:
#'   `mfrm_unexpected`, `mfrm_fair_average`, `mfrm_displacement`,
#'   `mfrm_interrater`, `mfrm_facets_chisq`, `mfrm_bias_interaction`,
#'   `mfrm_rating_scale`, `mfrm_category_structure`, `mfrm_category_curves`,
#'   `mfrm_measurable`, `mfrm_unexpected_after_bias`, `mfrm_output_bundle`,
#'   `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`,
#'   `mfrm_iteration_report`, `mfrm_subset_connectivity`,
#'   `mfrm_facet_statistics`, `mfrm_parity_report`.
#'
#' @section `plot.mfrm_bundle()` coverage:
#' Default dispatch now covers:
#' - `mfrm_unexpected`, `mfrm_fair_average`, `mfrm_displacement`
#' - `mfrm_interrater`, `mfrm_facets_chisq`, `mfrm_bias_interaction`
#' - `mfrm_bias_count`, `mfrm_fixed_reports`, `mfrm_visual_summaries`
#' - `mfrm_category_structure`, `mfrm_category_curves`, `mfrm_rating_scale`
#' - `mfrm_measurable`, `mfrm_unexpected_after_bias`, `mfrm_output_bundle`
#' - `mfrm_residual_pca`, `mfrm_specifications`, `mfrm_data_quality`
#' - `mfrm_iteration_report`, `mfrm_subset_connectivity`, `mfrm_facet_statistics`
#' - `mfrm_parity_report`
#'
#' For unknown bundle classes, use dedicated plotting helpers or custom base-R
#' plots from component tables.
#'
#' @seealso [fit_mfrm()], [run_mfrm_facets()], [mfrmRFacets()],
#'   [diagnose_mfrm()], [estimate_bias()],
#'   [summary.mfrm_fit()], `summary(diag)`,
#'   `summary()`, [plot.mfrm_fit()], `plot()`
#'
#' @examples
#' toy <- expand.grid(
#'   Person = paste0("P", 1:4),
#'   Rater = paste0("R", 1:2),
#'   Criterion = c("Content", "Organization", "Language"),
#'   stringsAsFactors = FALSE
#' )
#' toy$Score <- (
#'   as.integer(factor(toy$Person)) +
#'   2 * as.integer(factor(toy$Rater)) +
#'   as.integer(factor(toy$Criterion))
#' ) %% 3
#'
#' fit <- fit_mfrm(
#'   toy,
#'   person = "Person",
#'   facets = c("Rater", "Criterion"),
#'   score = "Score",
#'   method = "JML",
#'   maxit = 25
#' )
#' class(summary(fit))
#'
#' diag <- diagnose_mfrm(fit, residual_pca = "none")
#' class(summary(diag))
#'
#' t4 <- unexpected_response_table(fit, diagnostics = diag, top_n = 10)
#' class(summary(t4))
#' p <- plot(t4, draw = FALSE)
#'
#' @name mfrmr_workflow_methods
NULL
