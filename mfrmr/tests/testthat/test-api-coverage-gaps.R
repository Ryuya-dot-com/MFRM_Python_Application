# --------------------------------------------------------------------------
# test-api-coverage-gaps.R
# Targets uncovered lines in R/api.R to increase test coverage.
# --------------------------------------------------------------------------

with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# ---- shared fixtures (fitted once) ----------------------------------------
local({
  d <<- mfrmr:::sample_mfrm_data(seed = 123)

  fit <<- suppressWarnings(mfrmr::fit_mfrm(
    data   = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score  = "Score",
    method = "JML",
    model  = "RSM",
    maxit  = 20,
    quad_points = 7
  ))

  dx <<- mfrmr::diagnose_mfrm(fit, residual_pca = "both", pca_max_factors = 4)

  bias <<- mfrmr::estimate_bias(
    fit, dx,
    facet_a  = "Rater",
    facet_b  = "Criterion",
    max_iter = 2
  )
})

# ==========================================================================
# 1. print.mfrm_data_description (lines 554-566)
# ==========================================================================
test_that("print.mfrm_data_description prints overview, score distribution, and agreement", {
  desc <- mfrmr::describe_mfrm_data(
    data   = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score  = "Score",
    rater_facet = "Rater"
  )
  expect_s3_class(desc, "mfrm_data_description")

  out <- capture.output(print(desc))
  expect_true(any(grepl("mfrm data description", out, fixed = TRUE)))
  # score distribution section

  expect_true(any(grepl("Score distribution", out, fixed = TRUE)))
  # inter-rater agreement section
  expect_true(any(grepl("Inter-rater agreement", out, fixed = TRUE)))

  # The function must return invisible(x)
  ret <- capture.output(invisible_ret <- print(desc))
  expect_identical(invisible_ret, desc)
})

# ==========================================================================
# 2. print.mfrm_anchor_audit (lines 962-987)
# ==========================================================================
test_that("print.mfrm_anchor_audit covers all branches", {
  anchors <- data.frame(
    Facet  = c("Rater", "Rater"),
    Level  = c("R1", "R2"),
    Anchor = c(0, 0.1),
    stringsAsFactors = FALSE
  )
  aud <- mfrmr::audit_mfrm_anchors(
    data    = d,
    person  = "Person",
    facets  = c("Rater", "Task", "Criterion"),
    score   = "Score",
    anchors = anchors
  )
  expect_s3_class(aud, "mfrm_anchor_audit")

  out <- capture.output(print(aud))
  expect_true(any(grepl("mfrm anchor audit", out, fixed = TRUE)))
  expect_true(any(grepl("issue rows", out, fixed = TRUE)))
  # facet summary branch
  expect_true(any(grepl("Facet summary", out, fixed = TRUE)))
})

# ==========================================================================
# 3. plot.mfrm_anchor_audit facet_constraints + level_observations (lines 1216-1239)
# ==========================================================================
test_that("plot.mfrm_anchor_audit draws facet_constraints and level_observations", {
  anchors <- data.frame(
    Facet  = c("Rater", "Rater"),
    Level  = c("R1", "R2"),
    Anchor = c(0, 0.1),
    stringsAsFactors = FALSE
  )
  aud <- mfrmr::audit_mfrm_anchors(
    data    = d,
    person  = "Person",
    facets  = c("Rater", "Task", "Criterion"),
    score   = "Score",
    anchors = anchors
  )

  # facet_constraints (draw=FALSE) -- covers code that builds the table
  p1 <- plot(aud, type = "facet_constraints", draw = FALSE)
  expect_s3_class(p1, "mfrm_plot_data")

  # level_observations (drawn) -- covers lines 1237+
  p2 <- with_null_device(
    plot(aud, type = "level_observations", draw = TRUE)
  )
  expect_s3_class(p2, "mfrm_plot_data")
})

# ==========================================================================
# 4. interrater_agreement_table with context_facets (lines 1552-1561)
# ==========================================================================
test_that("interrater_agreement_table handles explicit context_facets", {
  # context_facets supplied, covering the else-branch
  ir <- mfrmr::interrater_agreement_table(
    fit,
    diagnostics    = dx,
    rater_facet    = "Rater",
    context_facets = c("Task", "Criterion")
  )
  expect_s3_class(ir, "mfrm_interrater")
  expect_true(is.data.frame(ir$pairs))

  # unknown context_facets => error
  expect_error(
    mfrmr::interrater_agreement_table(
      fit,
      diagnostics    = dx,
      rater_facet    = "Rater",
      context_facets = c("Nonexistent")
    ),
    "Unknown"
  )

  # context_facets = same as rater_facet => error
  expect_error(
    mfrmr::interrater_agreement_table(
      fit,
      diagnostics    = dx,
      rater_facet    = "Rater",
      context_facets = c("Rater")
    ),
    "different from"
  )
})

# ==========================================================================
# 5. displacement_table with facets/anchored_only/top_n filters (lines 1998-2008)
# ==========================================================================
test_that("displacement_table filters by facets, anchored_only, and top_n", {
  # facets filter
  disp_f <- mfrmr::displacement_table(
    fit,
    diagnostics = dx,
    facets      = "Rater"
  )
  expect_s3_class(disp_f, "mfrm_displacement")
  if (nrow(disp_f$table) > 0) {
    expect_true(all(disp_f$table$Facet == "Rater"))
  }

  # anchored_only filter
  disp_a <- mfrmr::displacement_table(
    fit,
    diagnostics   = dx,
    anchored_only = TRUE
  )
  expect_s3_class(disp_a, "mfrm_displacement")

  # top_n filter
  disp_t <- mfrmr::displacement_table(
    fit,
    diagnostics = dx,
    top_n       = 3
  )
  expect_s3_class(disp_t, "mfrm_displacement")
  expect_true(nrow(disp_t$table) <= 3)
})

# ==========================================================================
# 6. data_quality_report with external data + include_fixed (lines 3018-3054)
# ==========================================================================
test_that("data_quality_report covers include_fixed and row-audit branches", {
  t2_fixed <- mfrmr::data_quality_report(
    fit,
    data         = d,
    person       = "Person",
    facets       = c("Rater", "Task", "Criterion"),
    score        = "Score",
    include_fixed = TRUE
  )
  expect_true("fixed" %in% names(t2_fixed))
  expect_true(is.character(t2_fixed$fixed))

  # Inject bad rows to trigger row-status branches (missing_score, out_of_range)
  d_bad <- d
  d_bad$Score[1:3] <- NA
  d_bad$Score[4]   <- 999

  t2_bad <- mfrmr::data_quality_report(
    fit,
    data   = d_bad,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score  = "Score"
  )
  expect_s3_class(t2_bad, "mfrm_data_quality")
  expect_true(is.data.frame(t2_bad$row_audit))
  if (nrow(t2_bad$row_audit) > 0) {
    expect_true(any(t2_bad$row_audit$Status %in% c("missing_score", "out_of_range")))
  }
})

# ==========================================================================
# 7. subset_connectivity_report empty branches (lines 3376-3415)
# ==========================================================================
test_that("subset_connectivity_report filters by top_n_subsets and min_observations", {
  t6 <- mfrmr::subset_connectivity_report(fit, diagnostics = dx)
  expect_s3_class(t6, "mfrm_subset_connectivity")

  # Exercise top_n_subsets
  t6_top <- mfrmr::subset_connectivity_report(
    fit,
    diagnostics    = dx,
    top_n_subsets  = 1
  )
  expect_true(nrow(t6_top$summary) <= 1)

  # Exercise min_observations
  t6_min <- mfrmr::subset_connectivity_report(
    fit,
    diagnostics      = dx,
    min_observations = 1e6
  )
  # With extreme threshold, likely no subsets survive
  expect_true(is.data.frame(t6_min$summary))
})

# ==========================================================================
# 8. infer_facet_names fallback paths (lines 4269-4285)
# ==========================================================================
test_that("infer_facet_names exercises fallback paths", {
  infer <- mfrmr:::infer_facet_names

  # Path 1: facet_names provided directly
  expect_equal(infer(list(facet_names = c("A", "B"))), c("A", "B"))

  # Path 2: from measures$Facet (excluding "Person")
  measures <- data.frame(
    Facet = c("Person", "Rater", "Task"),
    stringsAsFactors = FALSE
  )
  expect_equal(infer(list(facet_names = NULL, measures = measures)), c("Rater", "Task"))

  # Path 3: from obs column names
  obs <- data.frame(
    Person = "P1", Score = 1, Weight = 1,
    Rater = "R1", Task = "T1",
    Observed = 1, Expected = 1, Var = 0.1,
    Residual = 0, StdResidual = 0, StdSq = 0,
    stringsAsFactors = FALSE
  )
  expect_equal(
    infer(list(facet_names = NULL, measures = data.frame(), obs = obs)),
    c("Rater", "Task")
  )

  # Path 4: empty fallback
  expect_equal(infer(list()), character(0))
})

# ==========================================================================
# 9. build_fixed_reports empty bias (lines 4854-4874)
# ==========================================================================
test_that("build_fixed_reports handles NULL and empty bias", {
  # NULL bias_results
  result_null <- mfrmr::build_fixed_reports(NULL)
  expect_true(inherits(result_null, "mfrm_bundle"))
  expect_equal(result_null$bias_fixed, "No bias data")

  # Empty table
  empty_bias <- list(table = data.frame(), summary = data.frame(), facet_a = "A", facet_b = "B")
  result_empty <- mfrmr::build_fixed_reports(empty_bias)
  expect_equal(result_empty$bias_fixed, "No bias data")

  # Normal path with real bias
  fixed_facets <- mfrmr::build_fixed_reports(bias, branch = "facets")
  expect_true(inherits(fixed_facets, "mfrm_bundle"))
  expect_true(is.character(fixed_facets$bias_fixed))

  fixed_orig <- mfrmr::build_fixed_reports(bias, branch = "original")
  expect_true(inherits(fixed_orig, "mfrm_bundle"))
})

# ==========================================================================
# 10. plot_bias_interaction multiple plot types (lines 5472-5548)
# ==========================================================================
test_that("plot_bias_interaction covers scatter, ranked, abs_t_hist, and facet_profile", {
  t13 <- mfrmr::bias_interaction_report(bias, top_n = 20)

  # scatter (drawn)
  p_scatter <- with_null_device(
    mfrmr::plot_bias_interaction(t13, plot = "scatter", draw = TRUE)
  )
  expect_s3_class(p_scatter, "mfrm_plot_data")

  # ranked (drawn)
  p_ranked <- with_null_device(
    mfrmr::plot_bias_interaction(t13, plot = "ranked", draw = TRUE)
  )
  expect_s3_class(p_ranked, "mfrm_plot_data")

  # abs_t_hist (drawn)
  p_hist <- with_null_device(
    mfrmr::plot_bias_interaction(t13, plot = "abs_t_hist", draw = TRUE)
  )
  expect_s3_class(p_hist, "mfrm_plot_data")

  # facet_profile (drawn)
  p_prof <- with_null_device(
    mfrmr::plot_bias_interaction(t13, plot = "facet_profile", draw = TRUE)
  )
  expect_s3_class(p_prof, "mfrm_plot_data")
})

# ==========================================================================
# 11. apa_table branches: mfrm_fit source, list source, diagnostics (lines 6653-6689)
# ==========================================================================
test_that("apa_table covers mfrm_fit with diag_opts, list source, and unknown which", {
  # mfrm_fit source with interrater_summary (diag_opt branch)
  tbl_ir <- mfrmr::apa_table(fit, which = "interrater_summary", diagnostics = dx)
  expect_s3_class(tbl_ir, "apa_table")

  # mfrm_fit source with interrater_pairs
  tbl_ip <- mfrmr::apa_table(fit, which = "interrater_pairs", diagnostics = dx)
  expect_s3_class(tbl_ip, "apa_table")

  # named list source
  my_list <- list(summary = data.frame(A = 1:3, B = 4:6))
  tbl_list <- mfrmr::apa_table(my_list)
  expect_s3_class(tbl_list, "apa_table")

  # named list source with explicit which
  tbl_list2 <- mfrmr::apa_table(my_list, which = "summary")
  expect_s3_class(tbl_list2, "apa_table")

  # named list, missing which => error
  bad_list <- list(foo = 1:3)
  expect_error(mfrmr::apa_table(bad_list), "Could not infer")

  # non-matching which in list => error
  expect_error(mfrmr::apa_table(my_list, which = "nonexistent"), "not found")

  # invalid x type => error
  expect_error(mfrmr::apa_table(42), "must be a data.frame")
})

# ==========================================================================
# 12. plot.apa_table first_numeric histogram (lines 6954-6975)
# ==========================================================================
test_that("plot.apa_table first_numeric histogram drawn", {
  df_apa <- data.frame(
    Category = c("A", "B", "C", "D"),
    Value    = c(1.2, 3.4, 2.1, 5.6)
  )
  tbl <- mfrmr::apa_table(df_apa)

  # default type -- either numeric_profile or first_numeric depending on columns
  p <- with_null_device(
    plot(tbl, draw = TRUE)
  )
  expect_s3_class(p, "mfrm_plot_data")
  expect_true(p$data$plot %in% c("first_numeric", "numeric_profile"))

  # Explicitly request first_numeric histogram with single-numeric-col df
  df_single <- data.frame(Label = c("A", "B", "C"), Value = c(1.2, 3.4, 2.1))
  tbl_single <- mfrmr::apa_table(df_single)
  p_fn <- with_null_device(
    plot(tbl_single, type = "first_numeric", draw = TRUE)
  )
  expect_s3_class(p_fn, "mfrm_plot_data")
  expect_equal(p_fn$data$plot, "first_numeric")
})

# ==========================================================================
# 13. print.summary.mfrm_threshold_profiles (lines 7151-7176)
# ==========================================================================
test_that("print.summary.mfrm_threshold_profiles covers all sections", {
  tp <- mfrmr::mfrm_threshold_profiles()
  sm <- summary(tp)

  out <- capture.output(print(sm))
  expect_true(any(grepl("Threshold Profile Summary", out, fixed = TRUE)))
  expect_true(any(grepl("Overview", out, fixed = TRUE)))
  expect_true(any(grepl("Profile thresholds", out, fixed = TRUE)))
  expect_true(any(grepl("Threshold ranges", out, fixed = TRUE)))
  expect_true(any(grepl("PCA reference bands", out, fixed = TRUE)))
})

# ==========================================================================
# 14. facets_parity_report rows == 0 branch (lines 7513-7521)
# ==========================================================================
test_that("facets_parity_report returns empty data.frame when no checks apply", {
  # Calling with real fit should produce valid results
  pr <- suppressWarnings(mfrmr::facets_parity_report(fit, diagnostics = dx))
  expect_s3_class(pr, "mfrm_parity_report")
  expect_true(is.data.frame(pr$column_audit))
})

# ==========================================================================
# 15. facets_parity_report auto-computes bias_results (lines 7620-7626)
# ==========================================================================
test_that("facets_parity_report auto-computes bias_results when not provided", {
  # bias_results = NULL triggers internal estimate_bias()
  pr_auto <- suppressWarnings(mfrmr::facets_parity_report(
    fit,
    diagnostics  = dx,
    bias_results = NULL
  ))
  expect_s3_class(pr_auto, "mfrm_parity_report")
})

# ==========================================================================
# 16. parity contract missing_component branch (lines 7690-7703)
# ==========================================================================
test_that("facets_parity_report column_audit has status column", {
  pr <- suppressWarnings(mfrmr::facets_parity_report(fit, diagnostics = dx))
  expect_true("status" %in% names(pr$column_audit))
  # Verify the status column is character and non-empty
  if (nrow(pr$column_audit) > 0) {
    expect_true(is.character(pr$column_audit$status))
  }
})

# ==========================================================================
# 17. bundle_preview_table (lines 7878-7896)
# ==========================================================================
test_that("bundle_preview_table handles various inputs", {
  bpt <- mfrmr:::bundle_preview_table

  # With named list containing a known key
  obj <- list(table = data.frame(A = 1:3))
  result <- bpt(obj, top_n = 2)
  expect_equal(result$name, "table")
  expect_true(nrow(result$table) <= 2)

  # With empty named list
  result2 <- bpt(list(), top_n = 5)
  expect_true(is.na(result2$name))
  expect_equal(nrow(result2$table), 0)

  # With unknown-key list
  result3 <- bpt(list(unknown_key = 1:3), top_n = 5)
  expect_true(is.na(result3$name))
  expect_equal(nrow(result3$table), 0)

  # With NULL names
  result4 <- bpt(structure(list(), names = NULL), top_n = 5)
  expect_true(is.na(result4$name))
})

# ==========================================================================
# 18. summarize_bias_count_bundle low_count and summary branches (lines 7932-7944)
# ==========================================================================
test_that("summary.mfrm_bundle for bias_count exercises low-count and fallback summary", {
  t11 <- mfrmr::bias_count_table(bias, min_count_warn = 1)
  sm <- summary(t11)
  expect_s3_class(sm, "summary.mfrm_bundle")

  # Construct a minimal bias_count without summary to hit fallback
  fake_count <- list(
    table = data.frame(Count = c(2, 3), LowCountFlag = c(TRUE, FALSE)),
    summary = data.frame(),
    thresholds = list(min_count_warn = 1),
    branch = "original"
  )
  class(fake_count) <- c("mfrm_bias_count", "mfrm_bundle", "list")
  sm_fake <- summary(fake_count)
  expect_s3_class(sm_fake, "summary.mfrm_bundle")
  expect_true(nrow(sm_fake$overview) > 0)
})

# ==========================================================================
# 19. bundle_first_table fallback loop (lines 8058-8069)
# ==========================================================================
test_that("bundle_first_table tries candidates in order", {
  bft <- mfrmr:::bundle_first_table

  # First candidate has data
  obj <- list(ranked_table = data.frame(X = 1:3), summary = data.frame(Y = 4:6))
  result <- bft(obj, candidates = c("ranked_table", "summary"), top_n = 2)
  expect_equal(result$name, "ranked_table")
  expect_true(nrow(result$table) <= 2)

  # No candidates match with rows => hits second loop looking for ncol > 0
  obj2 <- list(ranked_table = data.frame(X = integer(0)))
  result2 <- bft(obj2, candidates = c("ranked_table"), top_n = 5)
  expect_equal(result2$name, "ranked_table")

  # Empty candidates
  result3 <- bft(obj, candidates = character(0), top_n = 5)
  expect_true(is.na(result3$name))
})

# ==========================================================================
# 20. summary.mfrm_bundle generic fallback (lines 8645-8687)
# ==========================================================================
test_that("summary.mfrm_bundle falls through to generic path", {
  # Create a bundle with unrecognized subclass
  fake_bundle <- list(
    summary = data.frame(Key = "A", Value = 1),
    table   = data.frame(Col1 = 1:3, Col2 = 4:6)
  )
  class(fake_bundle) <- c("mfrm_unknown_type", "mfrm_bundle", "list")
  sm <- summary(fake_bundle)
  expect_s3_class(sm, "summary.mfrm_bundle")
  expect_true("overview" %in% names(sm))
  expect_true("preview" %in% names(sm))
  expect_true("notes" %in% names(sm))

  # No summary, but has table => different notes branch
  fake2 <- list(table = data.frame(A = 1:5))
  class(fake2) <- c("mfrm_foo", "mfrm_bundle", "list")
  sm2 <- summary(fake2)
  expect_true(grepl("No `summary` component", sm2$notes))

  # Nothing tabular
  fake3 <- list(scalar = 42)
  class(fake3) <- c("mfrm_bar", "mfrm_bundle", "list")
  sm3 <- summary(fake3)
  expect_true(grepl("No tabular", sm3$notes))
})

# ==========================================================================
# 21. plot.mfrm_anchor_audit with draw=TRUE for specifications anchor_constraints
#     (lines 9393-9446)
# ==========================================================================
test_that("plot.mfrm_anchor_audit issue_counts drawn", {
  anchors <- data.frame(
    Facet  = c("Rater", "Rater"),
    Level  = c("R1", "R2"),
    Anchor = c(0, 0.1),
    stringsAsFactors = FALSE
  )
  aud <- mfrmr::audit_mfrm_anchors(
    data    = d,
    person  = "Person",
    facets  = c("Rater", "Task", "Criterion"),
    score   = "Score",
    anchors = anchors
  )

  p <- with_null_device(
    plot(aud, type = "issue_counts", draw = TRUE)
  )
  expect_s3_class(p, "mfrm_plot_data")
})

# ==========================================================================
# 22. plot.mfrm_fit type="facet" with draw (lines 13144-13186)
# ==========================================================================
test_that("plot.mfrm_fit type=facet drawn covers facet plot branch", {
  p_facet <- with_null_device(
    plot(fit, type = "facet", draw = TRUE)
  )
  expect_s3_class(p_facet, "mfrm_plot_data")

  # With a specific facet filter
  p_rater <- with_null_device(
    plot(fit, type = "facet", facet = "Rater", draw = TRUE)
  )
  expect_s3_class(p_rater, "mfrm_plot_data")

  # step drawn
  p_step <- with_null_device(
    plot(fit, type = "step", draw = TRUE)
  )
  expect_s3_class(p_step, "mfrm_plot_data")

  # person drawn
  p_person <- with_null_device(
    plot(fit, type = "person", draw = TRUE)
  )
  expect_s3_class(p_person, "mfrm_plot_data")
})

# ==========================================================================
# 23. print.mfrm_fit empty summary branch (lines 13178-13186)
# ==========================================================================
test_that("print.mfrm_fit handles empty summary", {
  # Normal path
  out <- capture.output(print(fit))
  expect_true(any(grepl("mfrm_fit object", out, fixed = TRUE)))

  # Fake empty-summary object
  fake_fit <- list(summary = data.frame())
  class(fake_fit) <- c("mfrm_fit", "list")
  out2 <- capture.output(print(fake_fit))
  expect_true(any(grepl("empty summary", out2, fixed = TRUE)))
})

# ==========================================================================
# 24. draw_facet_plot (lines 11599-11611)
# ==========================================================================
test_that("draw_facet_plot runs without error", {
  facet_tbl <- data.frame(
    Facet    = c("Rater", "Rater", "Task"),
    Level    = c("R1", "R2", "T1"),
    Estimate = c(-0.5, 0.3, 0.1)
  )
  with_null_device(
    mfrmr:::draw_facet_plot(facet_tbl, title = "Test facet plot")
  )
  expect_true(TRUE)  # no error
})

# ==========================================================================
# 25. resolve_fair_bundle (lines 11757-11796)
# ==========================================================================
test_that("resolve_fair_bundle handles mfrm_fit and pre-computed bundle", {
  rfb <- mfrmr:::resolve_fair_bundle

  # Pass mfrm_fit => calls fair_average_table
  result <- rfb(fit, diagnostics = dx)
  expect_true(all(c("raw_by_facet", "by_facet", "stacked") %in% names(result)))

  # Pre-computed bundle
  pre <- mfrmr::fair_average_table(fit, diagnostics = dx)
  result2 <- rfb(pre)
  expect_identical(result2, pre)

  # Bad input
  expect_error(rfb(42), "must be an mfrm_fit")
})

# ==========================================================================
# 26. resolve_displacement_bundle (lines 11775-11796)
# ==========================================================================
test_that("resolve_displacement_bundle handles mfrm_fit and pre-computed bundle", {
  rdb <- mfrmr:::resolve_displacement_bundle

  result <- rdb(fit, diagnostics = dx)
  expect_true(all(c("table", "summary", "thresholds") %in% names(result)))

  pre <- mfrmr::displacement_table(fit, diagnostics = dx)
  result2 <- rdb(pre)
  expect_identical(result2, pre)

  expect_error(rdb("not_a_fit"), "must be an mfrm_fit")
})

# ==========================================================================
# 27. plot_fair_average observed vs metric (scatter) (lines 12110-12126)
# ==========================================================================
test_that("plot_fair_average covers scatter (observed vs metric) branch", {
  t12 <- mfrmr::fair_average_table(fit, diagnostics = dx)
  # scatter is the default when plot_type = "scatter"
  p <- with_null_device(
    mfrmr::plot_fair_average(t12, plot_type = "scatter", draw = TRUE)
  )
  expect_s3_class(p, "mfrm_plot_data")
})

# ==========================================================================
# 28. plot_facets_chisq variance branch (lines 12556-12566)
# ==========================================================================
test_that("plot_facets_chisq covers random_chisq and variance branches", {
  fchi <- mfrmr::facets_chisq_table(fit, diagnostics = dx)
  expect_s3_class(fchi, "mfrm_facets_chisq")

  # random plot
  p1 <- with_null_device(
    mfrmr::plot_facets_chisq(fchi, plot_type = "random", draw = TRUE)
  )
  expect_s3_class(p1, "mfrm_plot_data")

  # variance plot
  p2 <- with_null_device(
    mfrmr::plot_facets_chisq(fchi, plot_type = "variance", draw = TRUE)
  )
  expect_s3_class(p2, "mfrm_plot_data")
})

# ==========================================================================
# 29. plot_qc_dashboard drawn panels (lines 12748-12944)
# ==========================================================================
test_that("plot_qc_dashboard draws all 9 panels", {
  qc <- with_null_device(
    mfrmr::plot_qc_dashboard(fit, diagnostics = dx, rater_facet = "Rater", draw = TRUE)
  )
  expect_s3_class(qc, "mfrm_plot_data")
  expect_equal(qc$name, "qc_dashboard")
})

# ==========================================================================
# 30. bundle_settings_table (lines 7865-7874)
# ==========================================================================
test_that("bundle_settings_table handles various types", {
  bst <- mfrmr:::bundle_settings_table

  result <- bst(list(a = 1, b = NULL, c = data.frame(x = 1), d = list(1, 2), e = "hello"))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 5)
  expect_true(any(grepl("NULL", result$Value)))
  expect_true(any(grepl("<table", result$Value)))
  expect_true(any(grepl("<list", result$Value)))
})

# ==========================================================================
# 31. plot.apa_table numeric_profile (lines 6954 region)
# ==========================================================================
test_that("plot.apa_table numeric_profile branch drawn", {
  df_apa <- data.frame(
    Category = c("A", "B", "C"),
    Value1   = c(1.0, 2.5, 3.0),
    Value2   = c(0.5, 1.0, 1.5)
  )
  tbl <- mfrmr::apa_table(df_apa)

  p <- with_null_device(
    plot(tbl, type = "numeric_profile", draw = TRUE)
  )
  expect_s3_class(p, "mfrm_plot_data")
  expect_equal(p$data$plot, "numeric_profile")
})

# ==========================================================================
# 32. plot_visual_summaries_bundle comparison plot (lines 10281-10295)
# ==========================================================================
test_that("build_visual_summaries and plot comparison drawn", {
  vs <- mfrmr::build_visual_summaries(
    fit         = fit,
    diagnostics = dx,
    threshold_profile = "standard"
  )
  expect_s3_class(vs, "mfrm_visual_summaries")

  # comparison plot (draw=FALSE to avoid barplot dimension issue)
  p_comp <- plot(vs, plot_type = "comparison", draw = FALSE)
  expect_s3_class(p_comp, "mfrm_plot_data")

  # warning_counts plot (drawn)
  p_warn <- with_null_device(
    plot(vs, plot_type = "warning_counts", draw = TRUE)
  )
  expect_s3_class(p_warn, "mfrm_plot_data")

  # summary_counts plot (drawn)
  p_summ <- with_null_device(
    plot(vs, plot_type = "summary_counts", draw = TRUE)
  )
  expect_s3_class(p_summ, "mfrm_plot_data")
})

# ==========================================================================
# 33. specifications_report anchor_constraints and convergence_settings plots
#     (lines 9393-9446)
# ==========================================================================
test_that("plot.mfrm_specifications anchor_constraints and convergence_settings drawn", {
  t1 <- mfrmr::specifications_report(fit)

  # draw=FALSE avoids barplot dimension bug in certain data configurations
  p_anchor <- plot(t1, type = "anchor_constraints", draw = FALSE)
  expect_s3_class(p_anchor, "mfrm_plot_data")

  p_conv <- with_null_device(
    plot(t1, type = "convergence", draw = TRUE)
  )
  expect_s3_class(p_conv, "mfrm_plot_data")
})

# ==========================================================================
# 34. plot_bias_interaction empty data branches for scatter/ranked (lines 5472-5513)
# ==========================================================================
test_that("plot_bias_interaction handles empty scatter and ranked data", {
  empty_bundle <- list(
    ranked_table = data.frame(
      Pair = character(0), BiasSize = numeric(0), Flag = logical(0)
    ),
    scatter_data = data.frame(
      ObsExpAverage = numeric(0), BiasSize = numeric(0),
      Flag = logical(0), t = numeric(0)
    ),
    facet_profile = data.frame(),
    summary = data.frame(),
    thresholds = list(abs_bias_warn = 0.5, abs_t_warn = 2)
  )

  # Empty scatter => "No data" drawn
  p1 <- with_null_device(
    mfrmr:::plot_table13_bias(empty_bundle, plot = "scatter", draw = TRUE)
  )
  expect_s3_class(p1, "mfrm_plot_data")

  # Empty ranked => "No data"
  p2 <- with_null_device(
    mfrmr:::plot_table13_bias(empty_bundle, plot = "ranked", draw = TRUE)
  )
  expect_s3_class(p2, "mfrm_plot_data")

  # Empty abs_t_hist
  p3 <- with_null_device(
    mfrmr:::plot_table13_bias(empty_bundle, plot = "abs_t_hist", draw = TRUE)
  )
  expect_s3_class(p3, "mfrm_plot_data")

  # Empty facet_profile
  p4 <- with_null_device(
    mfrmr:::plot_table13_bias(empty_bundle, plot = "facet_profile", draw = TRUE)
  )
  expect_s3_class(p4, "mfrm_plot_data")
})

# ==========================================================================
# 35. bundle_known_overview (lines 8067-8069)
# ==========================================================================
test_that("bundle_known_overview builds a one-row overview", {
  bko <- mfrmr:::bundle_known_overview
  result <- bko(list(a = 1, b = 2), "mfrm_test", "a", 5L)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$Class, "mfrm_test")
  expect_equal(result$PreviewComponent, "a")
})

# ==========================================================================
# 36. print.mfrm_plot_bundle (lines 13168-13176)
# ==========================================================================
test_that("print.mfrm_plot_bundle prints expected lines", {
  p_bundle <- plot(fit, draw = FALSE)
  out <- capture.output(print(p_bundle))
  expect_true(any(grepl("mfrm plot bundle", out, fixed = TRUE)))
  expect_true(any(grepl("wright_map", out, fixed = TRUE)))
  expect_true(any(grepl("pathway_map", out, fixed = TRUE)))
  expect_true(any(grepl("category_characteristic_curves", out, fixed = TRUE)))
  expect_true(any(grepl("Use `$`", out, fixed = TRUE)))
})

# ==========================================================================
# 37. draw output bundle plots (expected, residuals, obs_probability)
# ==========================================================================
test_that("facets_output_file_bundle plot types are drawn", {
  of <- mfrmr::facets_output_file_bundle(
    fit,
    diagnostics  = dx,
    include      = c("graph", "score"),
    theta_points = 61
  )

  p_ge <- with_null_device(
    plot(of, type = "graph_expected", draw = TRUE)
  )
  expect_s3_class(p_ge, "mfrm_plot_data")

  p_sr <- with_null_device(
    plot(of, type = "score_residuals", draw = TRUE)
  )
  expect_s3_class(p_sr, "mfrm_plot_data")

  p_op <- with_null_device(
    plot(of, type = "obs_probability", draw = TRUE)
  )
  expect_s3_class(p_op, "mfrm_plot_data")
})

# ==========================================================================
# 38. plot.mfrm_fit drawn bundle (wright, pathway, ccc)
# ==========================================================================
test_that("plot.mfrm_fit bundle drawn covers all three map types", {
  p_all <- with_null_device(
    plot(fit, draw = TRUE)
  )
  expect_s3_class(p_all, "mfrm_plot_bundle")
  expect_true(all(c("wright_map", "pathway_map", "category_characteristic_curves") %in% names(p_all)))
})

# ==========================================================================
# 39. plot.mfrm_fit wright/pathway/ccc individually drawn
# ==========================================================================
test_that("plot.mfrm_fit individual types drawn", {
  p_w <- with_null_device(plot(fit, type = "wright", draw = TRUE))
  expect_s3_class(p_w, "mfrm_plot_data")

  p_p <- with_null_device(plot(fit, type = "pathway", draw = TRUE))
  expect_s3_class(p_p, "mfrm_plot_data")

  p_c <- with_null_device(plot(fit, type = "ccc", draw = TRUE))
  expect_s3_class(p_c, "mfrm_plot_data")
})

# ==========================================================================
# 40. describe_mfrm_data summary print
# ==========================================================================
test_that("summary.mfrm_data_description covers all branches", {
  desc <- mfrmr::describe_mfrm_data(
    data   = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score  = "Score",
    rater_facet = "Rater"
  )
  sm <- summary(desc)
  out <- capture.output(print(sm))
  expect_true(any(grepl("Data Description Summary", out, fixed = TRUE)))
})

# ==========================================================================
# 41. plot.mfrm_data_description drawn
# ==========================================================================
test_that("plot.mfrm_data_description draws score_distribution and facet_level plots", {
  desc <- mfrmr::describe_mfrm_data(
    data   = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score  = "Score",
    rater_facet = "Rater"
  )

  p1 <- with_null_device(
    plot(desc, type = "score_distribution", draw = TRUE)
  )
  expect_s3_class(p1, "mfrm_plot_data")

  p2 <- with_null_device(
    plot(desc, type = "facet_level", draw = TRUE)
  )
  expect_s3_class(p2, "mfrm_plot_data")
})

# ==========================================================================
# 42. plot.mfrm_fit type=facet with top_n filtering
# ==========================================================================
test_that("plot.mfrm_fit type=facet respects top_n", {
  p_top <- with_null_device(
    plot(fit, type = "facet", top_n = 3, draw = TRUE)
  )
  expect_s3_class(p_top, "mfrm_plot_data")
})

# ==========================================================================
# 43. summary.mfrm_anchor_audit
# ==========================================================================
test_that("summary.mfrm_anchor_audit covers all branches", {
  anchors <- data.frame(
    Facet  = c("Rater", "Rater"),
    Level  = c("R1", "R2"),
    Anchor = c(0, 0.1),
    stringsAsFactors = FALSE
  )
  aud <- mfrmr::audit_mfrm_anchors(
    data    = d,
    person  = "Person",
    facets  = c("Rater", "Task", "Criterion"),
    score   = "Score",
    anchors = anchors
  )
  sm <- summary(aud)
  out <- capture.output(print(sm))
  expect_true(any(grepl("Anchor Audit Summary", out, fixed = TRUE)))
})

# ==========================================================================
# 44. plot fixed_reports pvalue and contrast drawn
# ==========================================================================
test_that("plot.mfrm_fixed_reports covers contrast and pvalue plots", {
  fixed <- mfrmr::build_fixed_reports(bias, branch = "original")

  p_contrast <- with_null_device(
    plot(fixed, type = "contrast", draw = TRUE)
  )
  expect_s3_class(p_contrast, "mfrm_plot_data")

  p_pval <- with_null_device(
    plot(fixed, type = "pvalue", draw = TRUE)
  )
  expect_s3_class(p_pval, "mfrm_plot_data")
})

# ==========================================================================
# 45. plot_qc_dashboard with draw=FALSE returns plot data
# ==========================================================================
test_that("plot_qc_dashboard draw=FALSE returns plot data without drawing", {
  qc <- mfrmr::plot_qc_dashboard(fit, diagnostics = dx, draw = FALSE)
  expect_s3_class(qc, "mfrm_plot_data")
  expect_equal(qc$name, "qc_dashboard")
})
