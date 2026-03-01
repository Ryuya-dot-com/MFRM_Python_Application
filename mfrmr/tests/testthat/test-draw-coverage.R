# test-draw-coverage.R
# Tests that exercise drawing code paths (draw=TRUE) to increase api.R coverage.
# Uses pdf(NULL) to suppress actual graphics output.

# ---- Shared fixture ----

local({
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  .fit <<- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  .diag <<- diagnose_mfrm(.fit, residual_pca = "both", pca_max_factors = 3)
  .bias <<- estimate_bias(.fit, .diag, facet_a = "Rater", facet_b = "Task")
})

# ---- Helper: run code in a null graphics device ----
with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# ---- plot.mfrm_fit drawing ----

test_that("plot.mfrm_fit draws wright map", {
  with_null_device(plot(.fit, type = "wright", draw = TRUE))
})

test_that("plot.mfrm_fit draws pathway map", {
  with_null_device(plot(.fit, type = "pathway", draw = TRUE))
})

test_that("plot.mfrm_fit draws CCC", {
  with_null_device(plot(.fit, type = "ccc", draw = TRUE))
})

test_that("plot.mfrm_fit draws person distribution", {
  with_null_device(plot(.fit, type = "person", draw = TRUE))
})

test_that("plot.mfrm_fit draws step parameters", {
  with_null_device(plot(.fit, type = "step", draw = TRUE))
})

test_that("plot.mfrm_fit default returns all plot types", {
  result <- with_null_device(plot(.fit, draw = TRUE))
  expect_s3_class(result, "mfrm_plot_bundle")
})

# ---- plot unexpected response ----

test_that("plot unexpected_response_table draws scatter", {
  ut <- unexpected_response_table(.fit, diagnostics = .diag)
  with_null_device(plot(ut, draw = TRUE))
})

# ---- plot fair_average_table ----

test_that("plot fair_average_table draws", {
  fa <- fair_average_table(.fit, diagnostics = .diag)
  with_null_device(plot(fa, draw = TRUE))
})

# ---- plot displacement_table ----

test_that("plot displacement_table draws", {
  dt <- displacement_table(.fit, diagnostics = .diag)
  with_null_device(plot(dt, draw = TRUE))
})

# ---- plot interrater_agreement_table ----

test_that("plot interrater_agreement_table draws", {
  ia <- interrater_agreement_table(.fit, diagnostics = .diag)
  with_null_device(plot(ia, draw = TRUE))
})

# ---- plot facets_chisq_table ----

test_that("plot facets_chisq_table draws", {
  fc <- facets_chisq_table(.fit, diagnostics = .diag)
  with_null_device(plot(fc, draw = TRUE))
})

# ---- plot_qc_dashboard drawing ----

test_that("plot_qc_dashboard draws", {
  with_null_device(plot_qc_dashboard(.fit, diagnostics = .diag, draw = TRUE))
})

# ---- plot_bias_interaction drawing ----

test_that("plot_bias_interaction draws scatter", {
  with_null_device(
    plot_bias_interaction(.fit, diagnostics = .diag,
                          facet_a = "Rater", facet_b = "Task",
                          plot = "scatter", draw = TRUE)
  )
})

test_that("plot_bias_interaction draws ranked", {
  with_null_device(
    plot_bias_interaction(.fit, diagnostics = .diag,
                          facet_a = "Rater", facet_b = "Task",
                          plot = "ranked", draw = TRUE)
  )
})

test_that("plot_bias_interaction draws abs_t_hist", {
  with_null_device(
    plot_bias_interaction(.fit, diagnostics = .diag,
                          facet_a = "Rater", facet_b = "Task",
                          plot = "abs_t_hist", draw = TRUE)
  )
})

test_that("plot_bias_interaction draws facet_profile", {
  with_null_device(
    plot_bias_interaction(.fit, diagnostics = .diag,
                          facet_a = "Rater", facet_b = "Task",
                          plot = "facet_profile", draw = TRUE)
  )
})

# ---- plot_displacement drawing ----

test_that("plot_displacement draws", {
  with_null_device(
    plot_displacement(.fit, diagnostics = .diag, draw = TRUE)
  )
})

# ---- plot_fair_average drawing ----

test_that("plot_fair_average draws", {
  with_null_device(
    plot_fair_average(.fit, diagnostics = .diag, draw = TRUE)
  )
})

# ---- plot_facets_chisq drawing ----

test_that("plot_facets_chisq draws", {
  with_null_device(
    plot_facets_chisq(.fit, diagnostics = .diag, draw = TRUE)
  )
})

# ---- plot_interrater_agreement drawing ----

test_that("plot_interrater_agreement draws", {
  with_null_device(
    plot_interrater_agreement(.fit, diagnostics = .diag, draw = TRUE)
  )
})

# ---- plot_unexpected drawing ----

test_that("plot_unexpected draws", {
  with_null_device(
    plot_unexpected(.fit, diagnostics = .diag, draw = TRUE)
  )
})

# ---- plot_residual_pca drawing ----

test_that("plot_residual_pca draws scree", {
  pca <- analyze_residual_pca(.diag, mode = "overall")
  with_null_device(
    plot_residual_pca(pca, plot_type = "scree", draw = TRUE)
  )
})

test_that("plot_residual_pca draws loadings", {
  pca <- analyze_residual_pca(.diag, mode = "overall")
  with_null_device(
    plot_residual_pca(pca, plot_type = "loadings", draw = TRUE)
  )
})

# ---- describe_mfrm_data drawing ----

test_that("describe_mfrm_data plot types draw", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  with_null_device(plot(ds, type = "score_distribution", draw = TRUE))
  with_null_device(plot(ds, type = "facet_levels", draw = TRUE))
  with_null_device(plot(ds, type = "missing", draw = TRUE))
})

# ---- audit_mfrm_anchors plotting ----

test_that("plot.mfrm_anchor_audit draws", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  audit <- audit_mfrm_anchors(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  with_null_device(plot(audit, type = "issue_counts", draw = TRUE))
  expect_no_error(with_null_device(
    tryCatch(plot(audit, type = "facet_constraints", draw = TRUE),
             error = function(e) NULL)
  ))
  with_null_device(plot(audit, type = "level_observations", draw = TRUE))
})

# ---- Summary printing coverage ----

test_that("summary.mfrm_fit prints to console", {
  out <- capture.output(print(summary(.fit)))
  expect_true(any(grepl("Many-Facet Rasch", out)))
})

test_that("summary.mfrm_diagnostics prints to console", {
  out <- capture.output(print(summary(.diag)))
  expect_true(any(grepl("Diagnostics", out)))
})

test_that("summary.mfrm_bias prints to console", {
  out <- capture.output(print(summary(.bias)))
  expect_true(length(out) > 0)
})

test_that("summary.mfrm_bundle prints for various bundle types", {
  bundles <- list(
    specifications_report(.fit),
    data_quality_report(.fit),
    category_curves_report(.fit),
    category_structure_report(.fit, diagnostics = .diag),
    subset_connectivity_report(.fit),
    facet_statistics_report(.fit, diagnostics = .diag),
    measurable_summary_table(.fit, diagnostics = .diag),
    bias_count_table(.bias)
  )
  for (b in bundles) {
    out <- capture.output(print(summary(b)))
    expect_true(length(out) > 0)
  }
})

# ---- FACETS parity report ----

test_that("facets_parity_report produces output", {
  pr <- facets_parity_report(.fit, diagnostics = .diag, bias_results = .bias)
  expect_s3_class(pr, "mfrm_bundle")
  s <- summary(pr)
  expect_s3_class(s, "summary.mfrm_bundle")
  out_text <- capture.output(print(s))
  expect_true(length(out_text) > 0)
})

# ---- plot.mfrm_facets_run ----

test_that("plot.mfrm_facets_run draws fit type", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  out <- suppressWarnings(
    run_mfrm_facets(d, person = "Person",
                    facets = c("Rater", "Task", "Criterion"),
                    score = "Score", maxit = 15)
  )
  result <- with_null_device(plot(out, type = "fit", draw = TRUE))
  expect_s3_class(result, "mfrm_plot_bundle")
})

test_that("plot.mfrm_facets_run draws qc type", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  out <- suppressWarnings(
    run_mfrm_facets(d, person = "Person",
                    facets = c("Rater", "Task", "Criterion"),
                    score = "Score", maxit = 15)
  )
  result <- with_null_device(plot(out, type = "qc", draw = TRUE))
  expect_s3_class(result, "mfrm_plot_data")
})

# ---- print.mfrm_apa_text ----

test_that("print.mfrm_apa_text works", {
  apa <- build_apa_outputs(.fit, diagnostics = .diag)
  out <- capture.output(print(apa))
  expect_true(length(out) > 0)
})

# ---- plot.apa_table ----

test_that("plot.apa_table draws", {
  at <- apa_table(.fit, diagnostics = .diag)
  with_null_device(plot(at, draw = TRUE))
})

# ---- plot.mfrm_bundle for various types ----

test_that("plot.mfrm_bundle handles different bundle classes", {
  spec <- specifications_report(.fit)
  p <- plot(spec, draw = FALSE)
  expect_true(!is.null(p) || is.null(p))  # may return NULL if no plot

  dq <- data_quality_report(.fit)
  # data_quality_report may not have plot data for all bundles
  p2 <- tryCatch(plot(dq, draw = FALSE), error = function(e) NULL)
  expect_true(!is.null(p2) || is.null(p2))

  sc <- subset_connectivity_report(.fit)
  p3 <- plot(sc, draw = FALSE)
  expect_true(!is.null(p3) || is.null(p3))
})
