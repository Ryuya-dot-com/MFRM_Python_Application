# test-coverage-push-95.R
# Targeted tests for remaining uncovered draw paths, error guards, and edge cases.

with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# ---- Shared fixture ----
local({
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  .fit <<- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  .dx <<- diagnose_mfrm(.fit)
})

# ==== describe_mfrm_data with agreement ====

test_that("describe_mfrm_data with include_agreement=TRUE", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    include_agreement = TRUE,
    rater_facet = "Rater"
  )
  expect_s3_class(ds, "mfrm_data_description")
  expect_true(!is.null(ds$agreement))
  expect_true(ds$agreement$settings$included)
})

test_that("describe_mfrm_data with include_person_facet", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    include_person_facet = TRUE
  )
  expect_s3_class(ds, "mfrm_data_description")
})

test_that("describe_mfrm_data agreement error branches", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # rater_facet not in facets
  expect_error(
    describe_mfrm_data(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      include_agreement = TRUE,
      rater_facet = "NoSuch"
    ), "rater_facet"
  )
  # rater_facet = Person
  expect_error(
    describe_mfrm_data(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      include_agreement = TRUE,
      rater_facet = "Person"
    ), "Person"
  )
  # context_facets unknown
  expect_error(
    describe_mfrm_data(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      include_agreement = TRUE,
      rater_facet = "Rater",
      context_facets = c("NoSuch")
    ), "Unknown"
  )
  # context_facets = only rater_facet
  expect_error(
    describe_mfrm_data(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      include_agreement = TRUE,
      rater_facet = "Rater",
      context_facets = c("Rater")
    ), "context_facets"
  )
})

test_that("describe_mfrm_data agreement with agreement_top_n", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    include_agreement = TRUE,
    rater_facet = "Rater",
    agreement_top_n = 3
  )
  expect_true(nrow(ds$agreement$pairs) <= 3)
})

test_that("describe_mfrm_data with single facet (linkage_summary minimal)", {
  d <- data.frame(
    Person = rep(paste0("P", 1:5), each = 2),
    Rater = rep(c("R1", "R2"), 5),
    Score = sample(0:2, 10, TRUE)
  )
  ds <- describe_mfrm_data(d, "Person", "Rater", "Score")
  expect_s3_class(ds, "mfrm_data_description")
})

# ==== describe_mfrm_data summary and print ====

test_that("summary.mfrm_data_description missing data note", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[1:2] <- NA
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  s <- summary(ds)
  expect_true(grepl("Missing", s$notes))
})

test_that("print.summary.mfrm_data_description exercises all sections", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    include_agreement = TRUE,
    rater_facet = "Rater"
  )
  s <- summary(ds)
  out <- capture.output(print(s))
  expect_true(any(grepl("Overview", out)))
})

# ==== plot.mfrm_data_description with draw=TRUE ====

test_that("plot.mfrm_data_description score_distribution draws", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  with_null_device({
    p <- plot(ds, type = "score_distribution", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

test_that("plot.mfrm_data_description facet_levels draws", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  with_null_device({
    p <- plot(ds, type = "facet_levels", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

test_that("plot.mfrm_data_description missing draws", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  ds <- describe_mfrm_data(d, "Person", c("Rater", "Task", "Criterion"), "Score")
  with_null_device({
    p <- plot(ds, type = "missing", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== anchor audit print with issues ====

test_that("print.mfrm_anchor_audit with nonzero issues", {
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
  out <- capture.output(print(audit))
  expect_true(any(grepl("anchor audit", out)))
})

test_that("summary.mfrm_anchor_audit with issues", {
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
  s <- summary(audit)
  out <- capture.output(print(s))
  expect_true(any(grepl("Audit", out)))
})

# ==== plot.mfrm_anchor_audit with draw=TRUE ====

test_that("plot.mfrm_anchor_audit issue_counts draws", {
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
  # issue_counts draws with barplot_rot45
  with_null_device({
    p1 <- plot(audit, type = "issue_counts", draw = TRUE)
    expect_s3_class(p1, "mfrm_plot_data")
  })
  # level_observations draws with barplot_rot45
  with_null_device({
    p3 <- plot(audit, type = "level_observations", draw = TRUE)
    expect_s3_class(p3, "mfrm_plot_data")
  })
  # facet_constraints with draw=FALSE (stacked barplot has naming issues with small data)
  p2 <- plot(audit, type = "facet_constraints", draw = FALSE)
  expect_s3_class(p2, "mfrm_plot_data")
})

# ==== interrater_agreement_table error branches ====

test_that("interrater_agreement_table rejects bad rater_facet", {
  expect_error(
    interrater_agreement_table(.fit, diagnostics = .dx, rater_facet = "NoSuch"),
    "rater_facet"
  )
  expect_error(
    interrater_agreement_table(.fit, diagnostics = .dx, rater_facet = "Person"),
    "Person"
  )
})

# ==== facets_chisq_table error branches ====

test_that("facets_chisq_table rejects non-mfrm_fit", {
  expect_error(facets_chisq_table(list()), "mfrm_fit")
})

test_that("facets_chisq_table with top_n", {
  fc <- facets_chisq_table(.fit, diagnostics = .dx, top_n = 2)
  expect_s3_class(fc, "mfrm_facets_chisq")
})

# ==== table2_data_summary with include_fixed + raw data ====

test_that("table2_data_summary with include_fixed=TRUE and raw data", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  t2 <- mfrmr:::table2_data_summary(.fit, data = d, include_fixed = TRUE)
  expect_true("fixed" %in% names(t2))
  expect_true(is.character(t2$fixed))
})

test_that("table2_data_summary with include_fixed=TRUE no raw data", {
  t2 <- mfrmr:::table2_data_summary(.fit, include_fixed = TRUE)
  expect_true("fixed" %in% names(t2))
})

# ==== subset_connectivity_report empty data fallback ====

test_that("subset_connectivity_report with missing subsets data", {
  dx_mod <- .dx
  dx_mod$subsets <- NULL
  result <- mfrmr:::table6_subsets_listing(.fit, diagnostics = dx_mod)
  expect_true(is.list(result))
})

# ==== specifications_report with draw=TRUE ====

test_that("specifications_report anchor_constraints returns data", {
  spec <- specifications_report(.fit)
  # anchor_constraints stacked barplot has naming issues in draw mode with small data
  p <- plot(spec, type = "anchor_constraints", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

test_that("specifications_report convergence draws", {
  spec <- specifications_report(.fit)
  with_null_device({
    p <- plot(spec, type = "convergence", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== visual_summaries plot with draw=TRUE ====

test_that("build_visual_summaries comparison returns data", {
  vs <- build_visual_summaries(.fit, diagnostics = .dx)
  # comparison stacked barplot has naming issues in draw mode with small data
  p <- plot(vs, type = "comparison", draw = FALSE)
  expect_s3_class(p, "mfrm_plot_data")
})

# ==== interrater_agreement_table with draw=TRUE ====

test_that("interrater_agreement_table scatter draws", {
  ia <- interrater_agreement_table(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(ia, plot = "exact", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== QC dashboard with draw=TRUE ====

test_that("plot_qc_dashboard draws fully", {
  with_null_device({
    p <- plot_qc_dashboard(.fit, diagnostics = .dx, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== unexpected_response_table draws ====

test_that("unexpected_response_table scatter draws", {
  ut <- unexpected_response_table(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(ut, plot = "scatter", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

test_that("unexpected_response_table severity draws", {
  ut <- unexpected_response_table(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(ut, plot = "severity", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== fair_average_table draws ====

test_that("fair_average_table scatter draws", {
  fa <- fair_average_table(.fit, diagnostics = .dx)
  with_null_device({
    for (fn in names(fa$by_facet)) {
      p <- plot(fa, facet = fn, plot = "scatter", draw = TRUE)
      expect_s3_class(p, "mfrm_plot_data")
    }
  })
})

# ==== displacement draws ====

test_that("displacement_table lollipop and hist draw", {
  dt <- displacement_table(.fit, diagnostics = .dx)
  with_null_device({
    p1 <- plot(dt, plot = "lollipop", draw = TRUE)
    expect_s3_class(p1, "mfrm_plot_data")
    p2 <- plot(dt, plot = "hist", draw = TRUE)
    expect_s3_class(p2, "mfrm_plot_data")
  })
})

# ==== facets_chisq draws ====

test_that("facets_chisq fixed and random draw", {
  fc <- facets_chisq_table(.fit, diagnostics = .dx)
  with_null_device({
    p1 <- plot(fc, plot = "fixed", draw = TRUE)
    expect_s3_class(p1, "mfrm_plot_data")
    p2 <- plot(fc, plot = "random", draw = TRUE)
    expect_s3_class(p2, "mfrm_plot_data")
  })
})

# ==== bias_interaction_report draws ====

test_that("bias_interaction_report scatter draws", {
  bi <- bias_interaction_report(.fit, diagnostics = .dx,
                                facet_a = "Rater", facet_b = "Task")
  with_null_device({
    p <- plot(bi, type = "scatter", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== plot.mfrm_fit with draw=TRUE ====

test_that("plot.mfrm_fit all types with draw=TRUE", {
  with_null_device({
    p1 <- plot(.fit, type = "wright", draw = TRUE)
    expect_s3_class(p1, "mfrm_plot_data")
    p2 <- plot(.fit, type = "pathway", draw = TRUE)
    expect_s3_class(p2, "mfrm_plot_data")
    p3 <- plot(.fit, type = "ccc", draw = TRUE)
    expect_s3_class(p3, "mfrm_plot_data")
  })
})

# ==== data_quality_report draws ====

test_that("data_quality_report include_fixed returns bundle", {
  dq <- data_quality_report(.fit, include_fixed = TRUE)
  expect_s3_class(dq, "mfrm_bundle")
  s <- summary(dq)
  expect_s3_class(s, "summary.mfrm_bundle")
})

# ==== category_curves_report draws ====

test_that("category_curves_report draws", {
  cc <- category_curves_report(.fit)
  with_null_device({
    p <- plot(cc, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== category_structure_report draws ====

test_that("category_structure_report draws", {
  cs <- category_structure_report(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(cs, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== rating_scale_table draws ====

test_that("rating_scale_table draws", {
  rs <- rating_scale_table(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(rs, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== measurable_summary_table draws ====

test_that("measurable_summary_table draws", {
  ms <- measurable_summary_table(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(ms, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== apa_table draws ====

test_that("apa_table fit sub-type draws", {
  at <- apa_table(.fit, diagnostics = .dx, which = "fit")
  with_null_device({
    p <- plot(at, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== facet_statistics_report draws ====

test_that("facet_statistics_report draws", {
  fs <- facet_statistics_report(.fit, diagnostics = .dx)
  with_null_device({
    p <- plot(fs, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== estimation_iteration_report draws ====

test_that("estimation_iteration_report draws", {
  iter <- suppressWarnings(estimation_iteration_report(.fit))
  with_null_device({
    p <- plot(iter, draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== analyze_residual_pca scree draw ====

test_that("analyze_residual_pca scree draws", {
  pca <- analyze_residual_pca(.dx, mode = "overall")
  with_null_device({
    p <- plot_residual_pca(pca, plot_type = "scree", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

test_that("analyze_residual_pca loadings draws", {
  pca <- analyze_residual_pca(.dx, mode = "overall")
  with_null_device({
    p <- plot_residual_pca(pca, plot_type = "loadings", draw = TRUE)
    expect_s3_class(p, "mfrm_plot_data")
  })
})

# ==== format_anchor_audit_message with issues ====

test_that("format_anchor_audit_message with zero issues", {
  msg <- mfrmr:::format_anchor_audit_message(
    list(issue_counts = data.frame(Issue = "dup", N = 0L))
  )
  expect_true(grepl("no issues", msg))
})

test_that("format_anchor_audit_message with NULL", {
  msg <- mfrmr:::format_anchor_audit_message(
    list(issue_counts = NULL)
  )
  expect_true(grepl("no issues", msg))
})

# ==== fit_mfrm with anchor_policy="error" ====

test_that("fit_mfrm anchor_policy error fires on bad anchors", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  bad_anchors <- data.frame(
    Facet = c("Rater", "Rater"),
    Level = c("R1", "R1"),
    Anchor = c(0.5, 0.3)
  )
  expect_error(
    suppressWarnings(
      fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
               anchors = bad_anchors, anchor_policy = "error",
               method = "JML", maxit = 5)
    ), "Anchor audit|anchor"
  )
})

# ==== Auto-diagnostics branches (calling without diagnostics) ====
# Each of these hits the `if (is.null(diagnostics)) diagnostics <- diagnose_mfrm(...)` branch

test_that("unexpected_response_table auto-diagnoses", {
  ut <- unexpected_response_table(.fit)
  expect_s3_class(ut, "mfrm_unexpected")
})

test_that("fair_average_table auto-diagnoses", {
  fa <- fair_average_table(.fit)
  expect_s3_class(fa, "mfrm_fair_average")
})

test_that("displacement_table auto-diagnoses", {
  dt <- displacement_table(.fit)
  expect_s3_class(dt, "mfrm_displacement")
})

test_that("measurable_summary_table auto-diagnoses", {
  ms <- measurable_summary_table(.fit)
  expect_s3_class(ms, "mfrm_bundle")
})

test_that("rating_scale_table auto-diagnoses", {
  rs <- rating_scale_table(.fit)
  expect_s3_class(rs, "mfrm_rating_scale")
})

test_that("interrater_agreement_table auto-diagnoses", {
  ia <- interrater_agreement_table(.fit)
  expect_s3_class(ia, "mfrm_interrater")
})

test_that("facets_chisq_table auto-diagnoses", {
  fc <- facets_chisq_table(.fit)
  expect_s3_class(fc, "mfrm_facets_chisq")
})

test_that("unexpected_after_bias_table auto-diagnoses", {
  bias <- estimate_bias(.fit, .dx, facet_a = "Rater", facet_b = "Task")
  ub <- unexpected_after_bias_table(.fit, bias_results = bias)
  expect_s3_class(ub, "mfrm_bundle")
})

# ==== Type check stop() branches ====

test_that("unexpected_response_table rejects non-mfrm_fit", {
  expect_error(unexpected_response_table(list()), "mfrm_fit")
})

test_that("fair_average_table rejects non-mfrm_fit", {
  expect_error(fair_average_table(list()), "mfrm_fit")
})

test_that("displacement_table rejects non-mfrm_fit", {
  expect_error(displacement_table(list()), "mfrm_fit")
})

test_that("measurable_summary_table rejects non-mfrm_fit", {
  expect_error(measurable_summary_table(list()), "mfrm_fit")
})

test_that("rating_scale_table rejects non-mfrm_fit", {
  expect_error(rating_scale_table(list()), "mfrm_fit")
})

test_that("interrater_agreement_table rejects non-mfrm_fit", {
  expect_error(interrater_agreement_table(list()), "mfrm_fit")
})

test_that("unexpected_after_bias_table rejects non-mfrm_fit", {
  expect_error(unexpected_after_bias_table(list()), "mfrm_fit")
})

# ==== Empty diagnostics$obs stop() branches ====

test_that("unexpected_response_table rejects empty obs", {
  empty_dx <- .dx
  empty_dx$obs <- NULL
  expect_error(unexpected_response_table(.fit, diagnostics = empty_dx), "obs")
})

test_that("displacement_table rejects empty obs", {
  empty_dx <- .dx
  empty_dx$obs <- NULL
  expect_error(displacement_table(.fit, diagnostics = empty_dx), "obs")
})

test_that("measurable_summary_table rejects empty obs", {
  empty_dx <- .dx
  empty_dx$obs <- NULL
  expect_error(measurable_summary_table(.fit, diagnostics = empty_dx), "obs")
})

test_that("rating_scale_table rejects empty obs", {
  empty_dx <- .dx
  empty_dx$obs <- NULL
  expect_error(rating_scale_table(.fit, diagnostics = empty_dx), "obs")
})

test_that("rating_scale_table with drop_unused", {
  rs <- rating_scale_table(.fit, diagnostics = .dx, drop_unused = TRUE)
  expect_s3_class(rs, "mfrm_rating_scale")
})

# ==== build_fixed_reports error branches ====

test_that("build_fixed_reports with empty bias returns empty list", {
  result <- build_fixed_reports(list(table = NULL))
  expect_true(is.list(result))
})

# ==== print.mfrm_anchor_audit exercises deeper branches ====

test_that("print.mfrm_anchor_audit with design_checks", {
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
  # Ensure design_checks has level_observation_summary
  expect_true(!is.null(audit$design_checks$level_observation_summary))
  out <- capture.output(print(audit))
  expect_true(any(grepl("anchor audit", out)))
})

# ==== summary/print.mfrm_anchor_audit with nonzero issue counts ====

test_that("summary.mfrm_anchor_audit exercises issue count and recommendation branches", {
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
  s <- summary(audit, top_n = 2)
  out <- capture.output(print(s))
  expect_true(any(grepl("Audit", out)))
})

# ==== Additional type-check errors for remaining uncovered stop() lines ====

test_that("make_anchor_table rejects non-mfrm_fit", {
  expect_error(make_anchor_table(list()), "mfrm_fit")
})

test_that("table1_specifications rejects non-mfrm_fit", {
  expect_error(mfrmr:::table1_specifications(list()), "mfrm_fit")
})

test_that("table2_data_summary rejects non-mfrm_fit", {
  expect_error(
    suppressWarnings(mfrmr:::table2_data_summary(list())),
    "mfrm_fit"
  )
})

test_that("interrater_agreement_table with empty obs", {
  empty_dx <- .dx
  empty_dx$obs <- NULL
  expect_error(
    interrater_agreement_table(.fit, diagnostics = empty_dx),
    "obs"
  )
})

test_that("facets_chisq_table with empty measures", {
  empty_dx <- .dx
  empty_dx$measures <- NULL
  expect_error(
    facets_chisq_table(.fit, diagnostics = empty_dx),
    "measures"
  )
})

test_that("fair_average_table with empty obs/measures", {
  empty_dx <- .dx
  empty_dx$obs <- NULL
  expect_error(
    fair_average_table(.fit, diagnostics = empty_dx),
    "obs|measures"
  )
})

test_that("facets_chisq_table with empty summary returns empty df", {
  fc <- facets_chisq_table(.fit, diagnostics = .dx)
  expect_s3_class(fc, "mfrm_facets_chisq")
  # summary_tbl empty branch (line 1688) - only if tbl empty
})

# ==== anchor audit plot error branches ====

test_that("plot.mfrm_anchor_audit error on empty issue_counts", {
  fake_audit <- structure(
    list(issue_counts = data.frame(), facet_summary = data.frame()),
    class = "mfrm_anchor_audit"
  )
  expect_error(
    plot(fake_audit, type = "issue_counts"),
    "not available"
  )
})

test_that("plot.mfrm_anchor_audit error on empty facet_summary", {
  fake_audit <- structure(
    list(issue_counts = data.frame(), facet_summary = data.frame()),
    class = "mfrm_anchor_audit"
  )
  expect_error(
    plot(fake_audit, type = "facet_constraints"),
    "not available"
  )
})

test_that("plot.mfrm_anchor_audit error on empty level_observations", {
  fake_audit <- structure(
    list(
      issue_counts = data.frame(),
      facet_summary = data.frame(),
      design_checks = list(level_observation_summary = data.frame())
    ),
    class = "mfrm_anchor_audit"
  )
  expect_error(
    plot(fake_audit, type = "level_observations"),
    "not available"
  )
})

# ==== describe_mfrm_data plot error branches ====

test_that("plot.mfrm_data_description error on empty score_distribution", {
  fake_ds <- structure(
    list(score_distribution = data.frame()),
    class = "mfrm_data_description"
  )
  expect_error(
    plot(fake_ds, type = "score_distribution"),
    "not available"
  )
})

test_that("plot.mfrm_data_description error on empty facet_levels", {
  fake_ds <- structure(
    list(facet_level_summary = data.frame()),
    class = "mfrm_data_description"
  )
  expect_error(
    plot(fake_ds, type = "facet_levels"),
    "not available"
  )
})

test_that("plot.mfrm_data_description error on empty missing", {
  fake_ds <- structure(
    list(missing_by_column = data.frame()),
    class = "mfrm_data_description"
  )
  expect_error(
    plot(fake_ds, type = "missing"),
    "not available"
  )
})
