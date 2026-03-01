# --------------------------------------------------------------------------
# test-bubble-chart.R
# Tests for plot_bubble() bubble chart function
# --------------------------------------------------------------------------

# ---- Helper fit object ---------------------------------------------------

local_fit <- function(envir = parent.frame()) {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
  fit
}

# ---- 1. Basic draw=FALSE structure --------------------------------------

test_that("plot_bubble returns mfrm_plot_data with draw=FALSE", {
  fit <- local_fit()
  out <- suppressWarnings(plot_bubble(fit, draw = FALSE))
  expect_s3_class(out, "mfrm_plot_data")
  expect_equal(out$name, "bubble")
  expect_true(is.data.frame(out$data$table))
  expect_true(nrow(out$data$table) > 0)
})

test_that("plot_bubble table excludes Person facet", {
  fit <- local_fit()
  out <- suppressWarnings(plot_bubble(fit, draw = FALSE))
  expect_false("Person" %in% out$data$table$Facet)
})

# ---- 2. fit_stat parameter ----------------------------------------------

test_that("plot_bubble works with Infit and Outfit", {
  fit <- local_fit()
  out_in <- suppressWarnings(plot_bubble(fit, fit_stat = "Infit", draw = FALSE))
  out_out <- suppressWarnings(plot_bubble(fit, fit_stat = "Outfit", draw = FALSE))
  expect_equal(out_in$data$fit_stat, "Infit")
  expect_equal(out_out$data$fit_stat, "Outfit")
})

# ---- 3. bubble_size parameter -------------------------------------------

test_that("plot_bubble works with all bubble_size options", {
  fit <- local_fit()
  for (bs in c("SE", "N", "equal")) {
    out <- suppressWarnings(plot_bubble(fit, bubble_size = bs, draw = FALSE))
    expect_equal(out$data$bubble_size, bs)
    expect_true(length(out$data$radius) > 0)
    expect_true(all(is.finite(out$data$radius)))
  }
})

# ---- 4. facets parameter ------------------------------------------------

test_that("plot_bubble filters by facets", {
  fit <- local_fit()
  out <- suppressWarnings(plot_bubble(fit, facets = "Rater", draw = FALSE))
  expect_true(all(out$data$table$Facet == "Rater"))
})

test_that("plot_bubble with non-existent facet errors", {
  fit <- local_fit()
  expect_error(
    suppressWarnings(plot_bubble(fit, facets = "NonExistent", draw = FALSE)),
    "No measures"
  )
})

# ---- 5. draw=TRUE produces base graphics --------------------------------

test_that("plot_bubble draws without error", {
  fit <- local_fit()
  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)
  expect_silent(suppressWarnings(plot_bubble(fit)))
})

# ---- 6. Explicit diagnostics input --------------------------------------

test_that("plot_bubble accepts diagnostics object directly", {
  fit <- local_fit()
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  out <- plot_bubble(diag, draw = FALSE)
  expect_s3_class(out, "mfrm_plot_data")
  expect_true(nrow(out$data$table) > 0)
})

test_that("plot_bubble accepts fit + diagnostics together", {
  fit <- local_fit()
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  out <- plot_bubble(fit, diagnostics = diag, draw = FALSE)
  expect_s3_class(out, "mfrm_plot_data")
})

# ---- 7. fit_range in output ----------------------------------------------

test_that("plot_bubble stores fit_range in output", {
  fit <- local_fit()
  custom_range <- c(0.7, 1.3)
  out <- suppressWarnings(plot_bubble(fit, fit_range = custom_range, draw = FALSE))
  expect_equal(out$data$fit_range, custom_range)
})

# ---- 8. CI visualization in plot.mfrm_fit -------------------------------

test_that("plot.mfrm_fit with show_ci=TRUE wright map runs without error", {
  fit <- local_fit()
  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)
  expect_silent(
    suppressWarnings(plot(fit, type = "wright", show_ci = TRUE, ci_level = 0.95))
  )
})

test_that("plot.mfrm_fit with show_ci=TRUE facet plot runs without error", {
  fit <- local_fit()
  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)
  expect_silent(
    suppressWarnings(plot(fit, type = "facet", show_ci = TRUE, ci_level = 0.95))
  )
})

test_that("plot.mfrm_fit show_ci=TRUE wright data includes SE", {
  fit <- local_fit()
  out <- suppressWarnings(plot(fit, type = "wright", show_ci = TRUE, draw = FALSE))
  loc <- out$data$locations
  expect_true("SE" %in% names(loc))
  # At least some facet-level SE should be finite
  facet_se <- loc$SE[loc$PlotType == "Facet level"]
  expect_true(any(is.finite(facet_se) & facet_se > 0))
})

test_that("plot.mfrm_fit show_ci=TRUE facet data includes SE", {
  fit <- local_fit()
  out <- suppressWarnings(plot(fit, type = "facet", show_ci = TRUE, draw = FALSE))
  facet_tbl <- out$data$facets
  expect_true("SE" %in% names(facet_tbl))
  expect_true(any(is.finite(facet_tbl$SE) & facet_tbl$SE > 0))
})

test_that("plot.mfrm_fit show_ci=FALSE does not include SE in wright data", {
  fit <- local_fit()
  out <- suppressWarnings(plot(fit, type = "wright", show_ci = FALSE, draw = FALSE))
  loc <- out$data$locations
  expect_false("SE" %in% names(loc))
})
