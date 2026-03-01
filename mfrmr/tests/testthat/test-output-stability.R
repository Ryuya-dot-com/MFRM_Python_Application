# --------------------------------------------------------------------------
# test-output-stability.R
# Regression guard: output structures must remain stable across versions.
# --------------------------------------------------------------------------

# === 5.1 fit_mfrm output structure ========================================

test_that("fit_mfrm returns all required components", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
  # Required top-level components
  expect_true(all(c("summary", "facets", "steps", "config", "prep", "opt") %in%
    names(fit)))
  # Person table structure
  expect_true(is.data.frame(fit$facets$person))
  expect_true(all(c("Person", "Estimate") %in% names(fit$facets$person)))
  # Others table structure
  expect_true(is.data.frame(fit$facets$others))
  expect_true(all(c("Facet", "Level", "Estimate") %in%
    names(fit$facets$others)))
  # Summary structure
  expect_true(is.data.frame(fit$summary))
  expect_true("LogLik" %in% names(fit$summary))
  # Steps structure
  expect_true(is.data.frame(fit$steps))
  expect_true(all(c("Step", "Estimate") %in% names(fit$steps)))
  # Column types
  expect_type(fit$facets$others$Estimate, "double")
  expect_type(fit$facets$others$Level, "character")
  expect_type(fit$facets$others$Facet, "character")
})

# === 5.2 MML output has SD column =========================================

test_that("MML fit includes person SD column", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "MML", maxit = 30))
  expect_true("SD" %in% names(fit$facets$person))
  expect_true(all(fit$facets$person$SD > 0))
})

# === 5.3 diagnose_mfrm output structure ===================================

test_that("diagnose_mfrm returns all required components", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  expect_s3_class(dx, "mfrm_diagnostics")
  # Required components
  required <- c("obs", "measures", "overall_fit", "reliability")
  expect_true(all(required %in% names(dx)))
  # measures columns
  expect_true(all(c("Facet", "Level", "Estimate", "SE", "Infit", "Outfit") %in%
    names(dx$measures)))
  # reliability columns
  expect_true(all(c("Facet", "Separation", "Reliability") %in%
    names(dx$reliability)))
  # obs should be data.frame
  expect_true(is.data.frame(dx$obs))
})

# === 5.4 Table builder output stability ====================================
#
# Shared fixture: fit + diagnostics computed once for table builder tests.

.stab_d   <- mfrmr:::sample_mfrm_data(seed = 42)
.stab_fit <- suppressWarnings(fit_mfrm(.stab_d, "Person",
  c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
.stab_dx  <- diagnose_mfrm(.stab_fit, residual_pca = "none")

test_that("unexpected_response_table returns table + summary", {
  res <- unexpected_response_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("table", "summary") %in% names(res)))
  expect_true(is.data.frame(res$table))
  expect_true(is.data.frame(res$summary))
})

test_that("fair_average_table returns stacked + by_facet", {
  res <- fair_average_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("stacked", "by_facet") %in% names(res)))
  expect_true(is.data.frame(res$stacked))
  expect_true(is.list(res$by_facet))
})

test_that("displacement_table returns table + summary", {
  res <- displacement_table(.stab_fit, diagnostics = .stab_dx,
    anchored_only = FALSE)
  expect_true(all(c("table", "summary") %in% names(res)))
  expect_true(is.data.frame(res$table))
  expect_true(is.data.frame(res$summary))
})

test_that("rating_scale_table returns category_table + threshold_table", {
  res <- rating_scale_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("category_table", "threshold_table") %in% names(res)))
  expect_true(is.data.frame(res$category_table))
  expect_true(is.data.frame(res$threshold_table))
})

test_that("measurable_summary_table returns summary + facet_coverage", {
  res <- measurable_summary_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("summary", "facet_coverage") %in% names(res)))
  expect_true(is.data.frame(res$summary))
  expect_true(is.data.frame(res$facet_coverage))
})

test_that("interrater_agreement_table returns pairs + summary", {
  res <- interrater_agreement_table(.stab_fit, diagnostics = .stab_dx,
    rater_facet = "Rater")
  expect_true(all(c("pairs", "summary") %in% names(res)))
  expect_true(is.data.frame(res$pairs))
  expect_true(is.data.frame(res$summary))
})

test_that("facets_chisq_table returns table + summary", {
  res <- facets_chisq_table(.stab_fit, diagnostics = .stab_dx)
  expect_true(all(c("table", "summary") %in% names(res)))
  expect_true(is.data.frame(res$table))
  expect_true(is.data.frame(res$summary))
})

# === 5.5 Extreme data - all same score =====================================

test_that("all-same-score data produces informative error", {
  d <- data.frame(
    Person = rep(paste0("P", 1:6), each = 3),
    Rater  = rep(paste0("R", 1:3), 6),
    Score  = rep(3, 18),
    stringsAsFactors = FALSE
  )
  expect_error(
    fit_mfrm(d, "Person", "Rater", "Score", method = "JML"),
    "one score"
  )
})

# === 5.6 All max score =====================================================

test_that("all-max-score data produces informative error", {
  d <- data.frame(
    Person = rep(paste0("P", 1:6), each = 3),
    Rater  = rep(paste0("R", 1:3), 6),
    Score  = rep(5, 18),
    stringsAsFactors = FALSE
  )
  expect_error(
    fit_mfrm(d, "Person", "Rater", "Score", method = "JML"),
    "one score"
  )
})

# === 5.7 Output length matches facet levels ================================

test_that("output rows match number of facet levels", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  expect_equal(nrow(fit$facets$person), length(unique(d$Person)))
  n_rater <- nrow(fit$facets$others |>
    dplyr::filter(Facet == "Rater"))
  expect_equal(n_rater, length(unique(d$Rater)))
  n_task <- nrow(fit$facets$others |>
    dplyr::filter(Facet == "Task"))
  expect_equal(n_task, length(unique(d$Task)))
})

# === 5.8 Numeric types for estimates and SEs ===============================

test_that("all Estimate and SE columns are numeric", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  # Estimate columns are numeric/double
  expect_type(fit$facets$person$Estimate, "double")
  expect_type(fit$facets$others$Estimate, "double")
  expect_type(dx$measures$Estimate, "double")
  expect_type(fit$steps$Estimate, "double")
  # SE columns are numeric and positive where present
  se_vals <- dx$measures$SE
  se_finite <- se_vals[is.finite(se_vals)]
  expect_true(all(se_finite > 0))
})

# === 5.9 Fit statistics are in expected range ==============================

test_that("fit statistics are in expected range", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  # Infit and Outfit are positive (> 0) where finite
  infit_vals <- dx$measures$Infit[is.finite(dx$measures$Infit)]
  outfit_vals <- dx$measures$Outfit[is.finite(dx$measures$Outfit)]
  expect_true(all(infit_vals > 0))
  expect_true(all(outfit_vals > 0))
  # PTMEA is between -1 and 1 where present
  if ("PTMEA" %in% names(dx$measures)) {
    ptmea_vals <- dx$measures$PTMEA[is.finite(dx$measures$PTMEA)]
    expect_true(all(ptmea_vals >= -1 & ptmea_vals <= 1))
  }
})
