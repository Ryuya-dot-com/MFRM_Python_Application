# test-edge-cases.R
# Tests for boundary conditions and unusual-but-valid data.
# All tests target the release code directly (no mocks).

# ---- Minimal viable data ----

test_that("fit_mfrm succeeds with minimal viable data", {
  d <- data.frame(
    Person = rep(c("P1", "P2", "P3", "P4"), each = 2),
    Rater = rep(c("R1", "R2"), 4),
    Score = c(0, 1, 1, 0, 0, 1, 1, 0)
  )
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", "Rater", "Score", method = "JML", maxit = 30)
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_true(is.data.frame(fit$summary))
  expect_true("Estimate" %in% names(fit$facets$others))
})

# ---- NA values dropped gracefully ----

test_that("fit_mfrm drops NA rows and still fits", {
  d <- data.frame(
    Person = c("P1", "P2", "P3", "P4", NA, "P5", "P6"),
    Rater = c("R1", "R2", "R1", "R2", "R1", "R2", "R1"),
    Score = c(0, 1, NA, 2, 1, 0, 1)
  )
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", "Rater", "Score", method = "JML", maxit = 30)
  )
  expect_s3_class(fit, "mfrm_fit")
})

# ---- Weight column handling ----

test_that("fit_mfrm handles observation weights correctly", {
  set.seed(42)
  d <- data.frame(
    Person = rep(paste0("P", 1:6), each = 3),
    Rater = rep(paste0("R", 1:3), 6),
    Score = sample(0:2, 18, replace = TRUE),
    W = rep(c(1, 2, 0.5), 6)
  )
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", "Rater", "Score", weight = "W", method = "JML", maxit = 30)
  )
  expect_s3_class(fit, "mfrm_fit")
})

# ---- Non-convergence detection ----

test_that("fit_mfrm warns about non-convergence with tiny maxit", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  expect_warning(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 1),
    "did not fully converge"
  )
})

# ---- PCM model path ----

test_that("fit_mfrm PCM mode works with step_facet", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             model = "PCM", method = "JML", step_facet = "Task", maxit = 20)
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_equal(fit$summary$Model[[1]], "PCM")
})

# ---- MML path with person estimates ----

test_that("fit_mfrm MML produces person EAP estimates with SD", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "MML", maxit = 30, quad_points = 7)
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_true("SD" %in% names(fit$facets$person))
  expect_true(all(fit$facets$person$SD > 0))
})

# ---- Diagnostics on fitted model ----

test_that("diagnose_mfrm produces all expected components", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  expect_s3_class(diag, "mfrm_diagnostics")
  expect_true(all(c("obs", "measures", "overall_fit", "reliability",
                     "unexpected", "fair_average", "displacement",
                     "interrater", "facets_chisq") %in% names(diag)))
})

test_that("diagnose_mfrm with PCA produces eigenvalue output", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit, residual_pca = "both", pca_max_factors = 3)
  expect_true(!is.null(diag$residual_pca_overall))
  expect_true(!is.null(diag$residual_pca_by_facet))
})

# ---- describe_mfrm_data ----

test_that("describe_mfrm_data works with minimal data", {
  d <- data.frame(
    Person = c("P1", "P2", "P3"),
    Rater = c("R1", "R2", "R1"),
    Score = c(0, 1, 2)
  )
  ds <- describe_mfrm_data(d, "Person", "Rater", "Score")
  expect_s3_class(ds, "mfrm_data_description")
  expect_equal(ds$overview$Observations, 3)
})

test_that("describe_mfrm_data summary and print work", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  s <- summary(ds)
  expect_s3_class(s, "summary.mfrm_data_description")
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
})

test_that("describe_mfrm_data plot types work", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  p1 <- plot(ds, type = "score_distribution", draw = FALSE)
  expect_s3_class(p1, "mfrm_plot_data")
  p2 <- plot(ds, type = "facet_levels", draw = FALSE)
  expect_s3_class(p2, "mfrm_plot_data")
  p3 <- plot(ds, type = "missing", draw = FALSE)
  expect_s3_class(p3, "mfrm_plot_data")
})

# ---- audit_mfrm_anchors ----

test_that("audit_mfrm_anchors works without anchors", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  result <- audit_mfrm_anchors(
    d, "Person", c("Rater", "Task", "Criterion"), "Score"
  )
  expect_s3_class(result, "mfrm_anchor_audit")
})

test_that("audit_mfrm_anchors detects issues with bad anchors", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  bad_anchors <- data.frame(
    Facet = c("Rater", "NonExistent"),
    Level = c("R1", "X"),
    Anchor = c(0.5, -0.5)
  )
  result <- audit_mfrm_anchors(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    anchors = bad_anchors
  )
  expect_true(sum(result$issue_counts$N, na.rm = TRUE) > 0)
})

# ---- mfrmRFacets alias ----

test_that("mfrmRFacets alias produces identical output to run_mfrm_facets", {
  d <- mfrmr:::sample_mfrm_data(seed = 77)
  out1 <- suppressWarnings(
    run_mfrm_facets(d, person = "Person",
                    facets = c("Rater", "Task", "Criterion"),
                    score = "Score", maxit = 10)
  )
  out2 <- suppressWarnings(
    mfrmRFacets(d, person = "Person",
                facets = c("Rater", "Task", "Criterion"),
                score = "Score", maxit = 10)
  )
  expect_equal(out1$fit$summary$LogLik, out2$fit$summary$LogLik)
  expect_s3_class(out2, "mfrm_facets_run")
})

# ---- Threshold profiles ----

test_that("mfrm_threshold_profiles returns all three profiles", {
  tp <- mfrm_threshold_profiles()
  expect_s3_class(tp, "mfrm_threshold_profiles")
  profiles <- tp$profiles
  expect_true(is.list(profiles))
  expect_true(all(c("strict", "standard", "lenient") %in% names(profiles)))
  for (p in profiles) {
    expect_true("n_obs_min" %in% names(p))
  }
})

# ---- build_visual_summaries with different profiles ----

test_that("build_visual_summaries works with all threshold profiles", {
  skip_on_cran()
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)

  for (profile in c("strict", "standard", "lenient")) {
    vs <- build_visual_summaries(fit, diagnostics = diag,
                                 threshold_profile = profile)
    expect_true(is.list(vs))
    expect_true("warning_map" %in% names(vs) || "summary_map" %in% names(vs))
  }
})
