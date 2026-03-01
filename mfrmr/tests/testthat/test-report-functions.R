# test-report-functions.R
# Tests for report-building and table functions in isolation.
# Uses a shared fixture fit + diagnostics for efficiency.

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

# ---- specifications_report ----

test_that("specifications_report returns a bundle", {
  spec <- specifications_report(.fit)
  expect_s3_class(spec, "mfrm_bundle")
  s <- summary(spec)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- estimation_iteration_report ----

test_that("estimation_iteration_report returns a bundle", {
  iter <- suppressWarnings(estimation_iteration_report(.fit))
  expect_s3_class(iter, "mfrm_bundle")
})

# ---- data_quality_report ----

test_that("data_quality_report returns a bundle", {
  dq <- data_quality_report(.fit)
  expect_s3_class(dq, "mfrm_bundle")
  s <- summary(dq)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- category_curves_report ----

test_that("category_curves_report returns a bundle", {
  cc <- category_curves_report(.fit)
  expect_s3_class(cc, "mfrm_bundle")
})

# ---- category_structure_report ----

test_that("category_structure_report returns a bundle", {
  cs <- category_structure_report(.fit, diagnostics = .diag)
  expect_s3_class(cs, "mfrm_bundle")
})

# ---- subset_connectivity_report ----

test_that("subset_connectivity_report returns a bundle", {
  sc <- subset_connectivity_report(.fit)
  expect_s3_class(sc, "mfrm_bundle")
})

# ---- facet_statistics_report ----

test_that("facet_statistics_report returns a bundle", {
  fs <- facet_statistics_report(.fit, diagnostics = .diag)
  expect_s3_class(fs, "mfrm_bundle")
  s <- summary(fs)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- unexpected_response_table ----

test_that("unexpected_response_table produces valid output", {
  ut <- unexpected_response_table(.fit, diagnostics = .diag)
  expect_s3_class(ut, "mfrm_unexpected")
  expect_true(all(c("table", "summary", "thresholds") %in% names(ut)))
  s <- summary(ut)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(ut, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- fair_average_table ----

test_that("fair_average_table produces valid output", {
  fa <- fair_average_table(.fit, diagnostics = .diag)
  expect_s3_class(fa, "mfrm_fair_average")
  expect_true("stacked" %in% names(fa))
  s <- summary(fa)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(fa, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- displacement_table ----

test_that("displacement_table produces valid output", {
  dt <- displacement_table(.fit, diagnostics = .diag)
  expect_s3_class(dt, "mfrm_displacement")
  s <- summary(dt)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(dt, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- measurable_summary_table ----

test_that("measurable_summary_table produces valid output", {
  ms <- measurable_summary_table(.fit, diagnostics = .diag)
  expect_s3_class(ms, "mfrm_bundle")
  s <- summary(ms)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- rating_scale_table ----

test_that("rating_scale_table produces valid output", {
  rs <- rating_scale_table(.fit, diagnostics = .diag)
  expect_s3_class(rs, "mfrm_rating_scale")
  s <- summary(rs)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- bias_count_table ----

test_that("bias_count_table produces valid output", {
  bc <- bias_count_table(.bias)
  expect_s3_class(bc, "mfrm_bundle")
  s <- summary(bc)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- unexpected_after_bias_table ----

test_that("unexpected_after_bias_table produces valid output", {
  ub <- unexpected_after_bias_table(.fit, bias_results = .bias, diagnostics = .diag)
  expect_s3_class(ub, "mfrm_bundle")
})

# ---- interrater_agreement_table ----

test_that("interrater_agreement_table produces valid output", {
  ia <- interrater_agreement_table(.fit, diagnostics = .diag)
  expect_s3_class(ia, "mfrm_interrater")
  s <- summary(ia)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(ia, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- facets_chisq_table ----

test_that("facets_chisq_table produces valid output", {
  fc <- facets_chisq_table(.fit, diagnostics = .diag)
  expect_s3_class(fc, "mfrm_facets_chisq")
  s <- summary(fc)
  expect_s3_class(s, "summary.mfrm_bundle")
  p <- plot(fc, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- bias_interaction_report ----

test_that("bias_interaction_report produces valid output", {
  bi <- bias_interaction_report(.fit, diagnostics = .diag,
                                facet_a = "Rater", facet_b = "Task")
  expect_s3_class(bi, "mfrm_bundle")
  s <- summary(bi)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ---- build_apa_outputs ----

test_that("build_apa_outputs produces structured APA text", {
  apa <- build_apa_outputs(.fit, diagnostics = .diag)
  expect_s3_class(apa, "mfrm_apa_outputs")
  expect_true("report_text" %in% names(apa))
  expect_true(nchar(apa$report_text) > 50)
  s <- summary(apa)
  expect_s3_class(s, "summary.mfrm_apa_outputs")
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
})

test_that("build_apa_outputs with bias produces extended text", {
  apa <- build_apa_outputs(.fit, diagnostics = .diag, bias = .bias)
  expect_true(nchar(apa$report_text) > 100)
})

# ---- build_fixed_reports ----

test_that("build_fixed_reports produces text reports", {
  fr <- build_fixed_reports(.bias)
  expect_true(is.list(fr))
  expect_true(length(fr) > 0)
})

# ---- build_visual_summaries ----

test_that("build_visual_summaries produces warning and summary maps", {
  vs <- build_visual_summaries(.fit, diagnostics = .diag)
  expect_true(is.list(vs))
  expect_true("warning_map" %in% names(vs) || "summary_map" %in% names(vs))
})

# ---- apa_table ----

test_that("apa_table produces structured output", {
  at <- apa_table(.fit, diagnostics = .diag)
  expect_s3_class(at, "apa_table")
  s <- summary(at)
  expect_s3_class(s, "summary.apa_table")
  out <- capture.output(print(s))
  expect_true(length(out) > 0)
})

# ---- analyze_residual_pca ----

test_that("analyze_residual_pca produces eigenvalue and loading output", {
  pca <- analyze_residual_pca(.diag, mode = "both")
  expect_s3_class(pca, "mfrm_residual_pca")
  s <- summary(pca)
  expect_s3_class(s, "summary.mfrm_bundle")
})

test_that("analyze_residual_pca accepts fit object directly", {
  pca <- analyze_residual_pca(.fit, mode = "overall")
  expect_s3_class(pca, "mfrm_residual_pca")
})

# ---- plot_residual_pca ----

test_that("plot_residual_pca produces plot bundles", {
  pca <- analyze_residual_pca(.diag, mode = "overall")
  p_scree <- plot_residual_pca(pca, plot_type = "scree", draw = FALSE)
  expect_s3_class(p_scree, "mfrm_plot_data")
})

# ---- plot.mfrm_fit specific types ----

test_that("plot.mfrm_fit supports all named types", {
  p_wright <- plot(.fit, type = "wright", draw = FALSE)
  expect_s3_class(p_wright, "mfrm_plot_data")

  p_pathway <- plot(.fit, type = "pathway", draw = FALSE)
  expect_s3_class(p_pathway, "mfrm_plot_data")

  p_ccc <- plot(.fit, type = "ccc", draw = FALSE)
  expect_s3_class(p_ccc, "mfrm_plot_data")

  p_person <- plot(.fit, type = "person", draw = FALSE)
  expect_s3_class(p_person, "mfrm_plot_data")

  p_step <- plot(.fit, type = "step", draw = FALSE)
  expect_s3_class(p_step, "mfrm_plot_data")
})

# ---- plot_qc_dashboard ----

test_that("plot_qc_dashboard returns a plot bundle", {
  p <- plot_qc_dashboard(.fit, diagnostics = .diag, draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ---- make_anchor_table ----

test_that("make_anchor_table extracts anchors from fitted model", {
  at <- make_anchor_table(.fit)
  expect_true(is.data.frame(at))
  expect_true(all(c("Facet", "Level") %in% names(at)))
})

test_that("make_anchor_table includes persons when requested", {
  at <- make_anchor_table(.fit, include_person = TRUE)
  expect_true("Person" %in% at$Facet || nrow(at) > 0)
})

# ---- Formatting helpers (internal) ----

test_that("py_style_format converts Python-style format strings", {
  fmt <- mfrmr:::py_style_format
  expect_equal(fmt("{:.2f}", 3.14159), "3.14")
  expect_equal(fmt("{:.0f}", 42.7), "43")
})

test_that("fmt_num formats numbers correctly", {
  fn <- mfrmr:::fmt_num
  expect_equal(fn(3.14159, 2), "3.14")
  expect_equal(fn(NA, 2), "NA")
})

test_that("fmt_count formats integers correctly", {
  fc <- mfrmr:::fmt_count
  expect_equal(fc(42), "42")
  expect_equal(fc(NA), "NA")
})

test_that("fmt_pvalue formats p-values correctly", {
  fp <- mfrmr:::fmt_pvalue
  expect_true(grepl("< .001", fp(0.0001)))
  expect_true(grepl("= ", fp(0.05)))
})
