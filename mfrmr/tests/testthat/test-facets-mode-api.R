test_that("run_mfrm_facets returns FACETS-mode workflow bundle", {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  d <- mfrmr:::sample_mfrm_data(seed = 321)

  out <- mfrmr::run_mfrm_facets(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    maxit = 20,
    top_n_interactions = 10
  )

  expect_s3_class(out, "mfrm_facets_run")
  expect_true(all(c("fit", "diagnostics", "iteration", "fair_average", "rating_scale", "run_info", "mapping") %in% names(out)))
  expect_true(is.data.frame(out$run_info))
  expect_equal(out$fit$summary$Model[[1]], "RSM")
  expect_equal(out$fit$summary$Method[[1]], "JMLE")

  out_summary <- summary(out, top_n = 5)
  expect_s3_class(out_summary, "summary.mfrm_facets_run")
  printed <- capture.output(print(out_summary))
  expect_true(any(grepl("FACETS-mode Workflow Summary", printed, fixed = TRUE)))

  p_fit <- plot(out, type = "fit", draw = FALSE)
  expect_s3_class(p_fit, "mfrm_plot_bundle")
  p_qc <- plot(out, type = "qc", draw = FALSE)
  expect_s3_class(p_qc, "mfrm_plot_data")
})

test_that("mfrmRFacets alias routes to run_mfrm_facets", {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  d <- mfrmr:::sample_mfrm_data(seed = 999)
  a <- mfrmr::run_mfrm_facets(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    maxit = 15
  )
  b <- mfrmr::mfrmRFacets(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    maxit = 15
  )

  expect_equal(a$mapping, b$mapping)
  expect_equal(a$fit$summary$Model[[1]], b$fit$summary$Model[[1]])
  expect_equal(a$fit$summary$Method[[1]], b$fit$summary$Method[[1]])
})

test_that("run_mfrm_facets accepts method/model options", {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

  d <- mfrmr:::sample_mfrm_data(seed = 2026)

  out_mml <- mfrmr::run_mfrm_facets(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    model = "RSM",
    method = "MML",
    quad_points = 7,
    maxit = 15
  )

  expect_equal(out_mml$fit$summary$Model[[1]], "RSM")
  expect_equal(out_mml$fit$summary$Method[[1]], "MML")
})
