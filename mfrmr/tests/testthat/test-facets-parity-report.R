test_that("facets_parity_report returns full-coverage audit in facets branch", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task", max_iter = 2)

  parity <- mfrmr::facets_parity_report(
    fit = fit,
    diagnostics = diag,
    bias_results = bias,
    branch = "facets"
  )

  expect_s3_class(parity, "mfrm_parity_report")
  expect_true(all(c(
    "overall", "column_summary", "column_audit", "missing_preview",
    "metric_summary", "metric_by_table", "metric_audit", "settings"
  ) %in% names(parity)))
  expect_true(is.data.frame(parity$overall))
  expect_true(is.data.frame(parity$column_audit))
  expect_true(is.data.frame(parity$metric_audit))

  expect_equal(parity$overall$ColumnMismatches[1], 0)
  expect_equal(parity$overall$ColumnMismatchRate[1], 0)
  expect_equal(parity$overall$MeanColumnCoverage[1], 1)
  expect_equal(parity$overall$MinColumnCoverage[1], 1)
  expect_equal(parity$overall$MeanColumnCoverageAvailable[1], 1)
  expect_equal(parity$overall$MinColumnCoverageAvailable[1], 1)
  expect_true(parity$overall$MetricFailed[1] <= 0)
})

test_that("facets_parity_report integrates with summary() and plot()", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task", max_iter = 2)

  parity <- mfrmr::facets_parity_report(
    fit = fit,
    diagnostics = diag,
    bias_results = bias,
    branch = "original"
  )

  s <- summary(parity)
  expect_s3_class(s, "summary.mfrm_bundle")
  expect_identical(s$summary_kind, "mfrm_parity_report")

  printed <- paste(capture.output(print(s)), collapse = "\n")
  expect_match(printed, "mfrmr FACETS Parity Summary", fixed = TRUE)

  p1 <- plot(parity, draw = FALSE)
  p2 <- plot(parity, type = "table_coverage", draw = FALSE)
  p3 <- plot(parity, type = "metric_status", draw = FALSE)
  p4 <- plot(parity, type = "metric_by_table", draw = FALSE)

  expect_s3_class(p1, "mfrm_plot_data")
  expect_s3_class(p2, "mfrm_plot_data")
  expect_s3_class(p3, "mfrm_plot_data")
  expect_s3_class(p4, "mfrm_plot_data")
})

test_that("facets_parity_report contract coverage includes unavailable rows", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20
  )

  parity <- mfrmr::facets_parity_report(
    fit = fit,
    branch = "facets",
    include_metrics = FALSE
  )

  expect_gt(parity$overall$ColumnMismatches[1], 0)
  expect_lt(parity$overall$MeanColumnCoverage[1], 1)
  expect_lt(parity$overall$MinColumnCoverage[1], 1)
  expect_equal(parity$overall$MeanColumnCoverageAvailable[1], 1)
  expect_equal(parity$overall$MinColumnCoverageAvailable[1], 1)
})
