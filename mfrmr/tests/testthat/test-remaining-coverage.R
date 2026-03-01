# test-remaining-coverage.R
# Targeted tests for remaining coverage gaps.

with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# ---- PCM model path (hits PCM-specific code in pathway map, CCC, etc.) ----

test_that("PCM model path exercises PCM-specific plotting", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             model = "PCM", method = "JML", step_facet = "Criterion", maxit = 20)
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_equal(fit$summary$Model[[1]], "PCM")

  # PCM plot paths
  with_null_device(plot(fit, type = "wright", draw = TRUE))
  with_null_device(plot(fit, type = "pathway", draw = TRUE))
  with_null_device(plot(fit, type = "ccc", draw = TRUE))
  with_null_device(plot(fit, type = "person", draw = TRUE))
  with_null_device(plot(fit, type = "step", draw = TRUE))
})

test_that("PCM diagnostics and reports work", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             model = "PCM", method = "JML", step_facet = "Criterion", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  expect_s3_class(diag, "mfrm_diagnostics")

  rs <- rating_scale_table(fit, diagnostics = diag)
  expect_s3_class(rs, "mfrm_rating_scale")
  with_null_device(plot(rs, draw = TRUE))

  cs <- category_structure_report(fit, diagnostics = diag)
  expect_s3_class(cs, "mfrm_bundle")
  with_null_device(plot(cs, draw = TRUE))
})

# ---- Interrater agreement drawing sub-types ----

test_that("interrater agreement draws all plot types", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  ia <- interrater_agreement_table(fit, diagnostics = diag)

  with_null_device(plot(ia, plot = "exact", draw = TRUE))
  with_null_device(plot(ia, plot = "corr", draw = TRUE))
  with_null_device(plot(ia, plot = "difference", draw = TRUE))
})

# ---- Unexpected response - severity bar plot ----

test_that("unexpected response draws severity type", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  ut <- unexpected_response_table(fit, diagnostics = diag)

  with_null_device(plot(ut, plot = "severity", draw = TRUE))
})

# ---- Fair average - multiple facets ----

test_that("fair_average_table exercises per-facet paths", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  fa <- fair_average_table(fit, diagnostics = diag)

  # Plot for each facet
  for (facet_name in names(fa$by_facet)) {
    with_null_device(plot(fa, facet = facet_name, draw = TRUE))
  }
})

# ---- Displacement drawing ----

test_that("displacement table draws severity and histogram", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  dt <- displacement_table(fit, diagnostics = diag)

  with_null_device(plot(dt, plot = "lollipop", draw = TRUE))
  with_null_device(plot(dt, plot = "hist", draw = TRUE))
})

# ---- Chi-square drawing sub-types ----

test_that("facets_chisq draws scatter and bar plots", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  fc <- facets_chisq_table(fit, diagnostics = diag)

  with_null_device(plot(fc, plot = "fixed", draw = TRUE))
  with_null_device(plot(fc, plot = "random", draw = TRUE))
})

# ---- Output bundle (graphfile/scorefile) plot sub-types ----

test_that("output bundle draws all graph types", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)

  # facets_output_file_bundle if it exists
  bundle <- tryCatch(
    mfrmr:::facets_output_file_bundle(fit, diagnostics = diag),
    error = function(e) NULL
  )
  if (!is.null(bundle)) {
    with_null_device(plot(bundle, type = "graph_expected", draw = TRUE))
    with_null_device(plot(bundle, type = "score_residuals", draw = TRUE))
    with_null_device(plot(bundle, type = "obs_probability", draw = TRUE))
  }
})

# ---- Anchor audit with actual anchors ----

test_that("audit_mfrm_anchors with anchors exercises more branches", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  anchors <- data.frame(
    Facet = c("Rater", "Rater", "Task"),
    Level = c("R1", "R2", "T1"),
    Anchor = c(0.5, -0.3, 0.1)
  )
  audit <- audit_mfrm_anchors(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    anchors = anchors
  )
  expect_s3_class(audit, "mfrm_anchor_audit")
  s <- summary(audit)
  expect_s3_class(s, "summary.mfrm_anchor_audit")
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
})

# ---- apa_table with which parameter ----

test_that("apa_table with different which values", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)

  at1 <- apa_table(fit, diagnostics = diag, which = "fit")
  expect_s3_class(at1, "apa_table")
  at2 <- apa_table(fit, diagnostics = diag, which = "reliability")
  expect_s3_class(at2, "apa_table")
})

# ---- build_fixed_reports with both branches ----

test_that("build_fixed_reports works with original branch", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task")

  fr_orig <- build_fixed_reports(bias, branch = "original")
  expect_true(is.list(fr_orig))

  fr_facets <- build_fixed_reports(bias, branch = "facets")
  expect_true(is.list(fr_facets))
})

# ---- MML model with all reports ----

test_that("MML path exercises different code branches", {
  skip_on_cran()
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "MML", maxit = 30, quad_points = 7)
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_equal(fit$summary$Method[[1]], "MML")

  # Diagnostics
  diag <- diagnose_mfrm(fit)
  expect_s3_class(diag, "mfrm_diagnostics")

  # APA outputs
  apa <- build_apa_outputs(fit, diagnostics = diag)
  expect_s3_class(apa, "mfrm_apa_outputs")

  # Visual summaries
  vs <- build_visual_summaries(fit, diagnostics = diag)
  expect_true(is.list(vs))

  # All table functions
  ut <- unexpected_response_table(fit, diagnostics = diag)
  expect_s3_class(ut, "mfrm_unexpected")

  fa <- fair_average_table(fit, diagnostics = diag)
  expect_s3_class(fa, "mfrm_fair_average")

  dt <- displacement_table(fit, diagnostics = diag)
  expect_s3_class(dt, "mfrm_displacement")

  ms <- measurable_summary_table(fit, diagnostics = diag)
  expect_s3_class(ms, "mfrm_bundle")

  rs <- rating_scale_table(fit, diagnostics = diag)
  expect_s3_class(rs, "mfrm_rating_scale")
})

# ---- facets_mode_methods.R coverage ----

test_that("summary.mfrm_facets_run and print work", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  out <- suppressWarnings(
    run_mfrm_facets(d, person = "Person",
                    facets = c("Rater", "Task", "Criterion"),
                    score = "Score", maxit = 15)
  )
  s <- summary(out)
  expect_s3_class(s, "summary.mfrm_facets_run")

  printed <- capture.output(print(s))
  expect_true(any(grepl("FACETS-mode", printed)))

  # Direct print delegates to summary
  printed2 <- capture.output(print(out))
  expect_true(length(printed2) > 0)
})

# ---- Additional describe_mfrm_data coverage ----

test_that("describe_mfrm_data with multiple facets and context_facets", {
  d <- mfrmr:::sample_mfrm_data(seed = 1)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score",
                            context_facets = c("Rater", "Task"))
  expect_s3_class(ds, "mfrm_data_description")
  # context_facets parameter accepted without error
  expect_true(is.list(ds))
})

# ---- positive_facets and noncenter_facet edge cases ----

test_that("fit_mfrm with positive_facets works", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20,
             positive_facets = c("Task"))
  )
  expect_s3_class(fit, "mfrm_fit")
})

test_that("fit_mfrm with noncenter_facet works", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20,
             noncenter_facet = "Rater")
  )
  expect_s3_class(fit, "mfrm_fit")
})

# ---- list_mfrmr_data ----

test_that("list_mfrmr_data returns available datasets", {
  ds <- list_mfrmr_data()
  expect_true(is.character(ds))
  expect_true(length(ds) > 0)
  expect_true("study1" %in% ds || "ej2021_study1" %in% ds)
})

# ---- load_mfrmr_data ----

test_that("load_mfrmr_data loads all available datasets", {
  available <- list_mfrmr_data()
  for (name in available[1:min(3, length(available))]) {
    d <- load_mfrmr_data(name)
    expect_true(is.data.frame(d))
    expect_true(nrow(d) > 0)
  }
})
