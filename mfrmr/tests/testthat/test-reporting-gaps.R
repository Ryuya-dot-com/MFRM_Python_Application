# test-reporting-gaps.R
# Targeted tests for uncovered branches in reporting.R and facets_mode files.

# ---- py_style_format fallback (line 20) ----

test_that("py_style_format returns character for non-matching format", {
  pf <- mfrmr:::py_style_format
  # No {:.Xf} match -> fallback to as.character()
  expect_equal(pf("{}", 42), "42")
  expect_equal(pf("plain", 3.14), "3.14")
})

# ---- format_fixed_width_table edge cases (lines 109-111) ----

test_that("format_fixed_width_table handles non-data.frame sections", {
  ffw <- mfrmr:::format_fixed_width_table
  # Empty data.frame returns "No data"
  expect_equal(ffw(data.frame(), columns = "X"), "No data")
  # NULL returns "No data"
  expect_equal(ffw(NULL, columns = "X"), "No data")
})

# ---- fmt_pvalue NA path (line 225) ----

test_that("fmt_pvalue returns NA for non-finite input", {
  fp <- mfrmr:::fmt_pvalue
  expect_equal(fp(NA), "NA")
  expect_equal(fp(Inf), "NA")
  expect_equal(fp(NaN), "NA")
})

# ---- safe_residual_pca with NULL (line 244) ----

test_that("safe_residual_pca returns NULL for NULL diagnostics", {
  srp <- mfrmr:::safe_residual_pca
  expect_null(srp(NULL))
})

# ---- resolve_warning_thresholds edge cases (line 317) ----

test_that("resolve_warning_thresholds falls back to standard for invalid profile", {
  rwt <- mfrmr:::resolve_warning_thresholds
  result <- rwt(threshold_profile = "nonexistent_profile")
  expect_equal(result$profile_name, "standard")
})

test_that("resolve_warning_thresholds applies custom threshold overrides", {
  rwt <- mfrmr:::resolve_warning_thresholds
  custom <- list(infit_lo = 0.3)
  result <- rwt(thresholds = custom, threshold_profile = "standard")
  expect_equal(result$thresholds$infit_lo, 0.3)
})

# ---- extract_pca functions edge cases (lines 362-397) ----

test_that("extract_overall_pca_first returns NULL for empty table", {
  ef <- mfrmr:::extract_overall_pca_first
  expect_null(ef(NULL))
  expect_null(ef(list(overall_table = NULL)))
  expect_null(ef(list(overall_table = data.frame())))
})

test_that("extract_overall_pca_second returns NULL for empty or single-row table", {
  es <- mfrmr:::extract_overall_pca_second
  expect_null(es(NULL))
  expect_null(es(list(overall_table = data.frame())))
  # Single component (no PC2)
  tbl <- data.frame(Component = 1, Eigenvalue = 2.5, Proportion = 0.5)
  expect_null(es(list(overall_table = tbl)))
})

test_that("extract_facet_pca_first returns empty data.frame for NULL input", {
  efp <- mfrmr:::extract_facet_pca_first
  expect_equal(nrow(efp(NULL)), 0)
  expect_equal(nrow(efp(list(by_facet_table = NULL))), 0)
  expect_equal(nrow(efp(list(by_facet_table = data.frame()))), 0)
  # Missing required columns
  expect_equal(nrow(efp(list(by_facet_table = data.frame(X = 1)))), 0)
})

# ---- collapse_apa_paragraph edge cases (lines 400-406) ----

test_that("collapse_apa_paragraph handles empty input and small width", {
  cap <- mfrmr:::collapse_apa_paragraph
  expect_equal(cap(character(0)), "")
  expect_equal(cap(NULL), "")
  expect_equal(cap(c("", " ")), "")
  # Width below minimum gets reset to 92
  result <- cap("A short sentence.", width = 10)
  expect_true(nchar(result) > 0)
})

# ---- summarize_anchor_constraints edge cases (lines 435, 440, 445) ----

test_that("summarize_anchor_constraints handles missing columns", {
  sac <- mfrmr:::summarize_anchor_constraints
  # Minimal config with no anchors
  config <- list(
    noncenter_facet = "none",
    facet_names = c("Rater"),
    facet_levels = list(Rater = c("R1", "R2")),
    dummy_facets = character(0)
  )
  result <- sac(config)
  expect_true(grepl("noncenter facet", result))
  expect_true(grepl("none", result))
})

# ---- summarize_step_estimates edge cases (lines 487-536) ----

test_that("summarize_step_estimates handles NULL step table", {
  sse <- mfrmr:::summarize_step_estimates
  expect_equal(sse(NULL), "Step/threshold estimates were not available.")
  expect_equal(sse(data.frame()), "Step/threshold estimates were not available.")
})

# ---- summarize_top_misfit_levels edge cases (lines 539-559) ----

test_that("summarize_top_misfit_levels handles NULL table", {
  stm <- mfrmr:::summarize_top_misfit_levels
  expect_equal(stm(NULL), "Top misfit levels were not available.")
  expect_equal(stm(data.frame()), "Top misfit levels were not available.")
})

# ---- build_apa_text with context options (lines 610-658) ----

test_that("build_apa_text respects context parameters", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)

  # With context options that trigger the assessment/setting branches
  apa <- build_apa_outputs(fit, diagnostics = diag, context = list(
    assessment = "essay writing",
    setting = "a university course",
    rater_training = "two-hour calibration",
    raters_per_response = "2",
    scale_desc = "holistic 0-3"
  ))
  expect_true(grepl("essay writing", apa$report_text))
  expect_true(grepl("university course", apa$report_text))
  expect_true(grepl("calibration", apa$report_text))
  expect_true(grepl("holistic", apa$report_text))
})

# ---- build_apa_table_figure_note_map edge cases (lines 859-943) ----

test_that("build_apa_table_figure_note_map produces note map", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)

  note_map <- mfrmr:::build_apa_table_figure_note_map(
    res = fit, diagnostics = diag, bias_results = NULL
  )
  expect_true(is.list(note_map))
  expect_true("table1" %in% names(note_map))
  expect_true("table2" %in% names(note_map))
  expect_true(grepl("Facet summary", note_map$table1))
})

# ---- build_visual_warning_map measures not available (lines 1079-1083) ----

test_that("build_visual_warning_map handles missing measures", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  # Null out measures to hit the missing measures branch
  diag$measures <- NULL

  wm <- mfrmr:::build_visual_warning_map(fit, diag)
  expect_true(is.list(wm))
  expect_true(any(grepl("not available", unlist(wm))))
})

# ---- build_visual_summary_map pathway/PCA unavailable (lines 1275, 1395, 1426) ----

test_that("build_visual_summary_map handles missing step and PCA data", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)

  sm <- mfrmr:::build_visual_summary_map(fit, diag)
  expect_true(is.list(sm))
})

# ---- facets_mode_api.R: normalize_spec_input (line 20) ----

test_that("normalize_spec_input converts valid data.frame", {
  nsi <- mfrmr:::normalize_spec_input
  result <- nsi(data.frame(Facet = "R", Level = "R1", Anchor = 0.5), "anchors")
  expect_true(is.data.frame(result))
})

# ---- facets_mode_api.R: weight not found (lines 47-48) ----

test_that("infer_facets_mode_mapping rejects missing weight column", {
  toy <- data.frame(Person = 1:3, Score = 1:3, Rater = c("A", "B", "C"))
  expect_error(
    mfrmr:::infer_facets_mode_mapping(toy, person = "Person",
                                       score = "Score", weight = "NoSuch"),
    "Weight column not found"
  )
})

# ---- facets_mode_api.R: auto-detect facet columns (lines 53-56) ----

test_that("infer_facets_mode_mapping auto-detects facet columns", {
  toy <- data.frame(Person = 1:5, Score = 1:5, Rater = paste0("R", 1:5),
                    Task = paste0("T", 1:5))
  result <- mfrmr:::infer_facets_mode_mapping(toy, person = "Person", score = "Score")
  expect_true("Rater" %in% result$facets)
  expect_true("Task" %in% result$facets)
})

# ---- facets_mode_api.R: no facets detected (line 60) ----

test_that("infer_facets_mode_mapping errors when no facet columns found", {
  # 3 columns but person and score take up 2, leaving only 1 for facets
  # Actually, we need a case where all remaining columns are blocked
  toy <- data.frame(Person = 1:3, Score = 1:3, Weight = c(1, 1, 1))
  expect_error(
    mfrmr:::infer_facets_mode_mapping(toy, person = "Person",
                                       score = "Score", weight = "Weight"),
    "No facet columns"
  )
})

# ---- facets_mode_methods.R: round_numeric_frame edge (line 4) ----

test_that("round_numeric_frame handles empty data.frame", {
  rnf <- mfrmr:::round_numeric_frame
  result <- rnf(data.frame())
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

# ---- facets_mode_methods.R: summary error path (line 54) ----

test_that("summary.mfrm_facets_run rejects non-mfrm_facets_run input", {
  expect_error(
    summary.mfrm_facets_run(list()),
    "mfrm_facets_run"
  )
})

# ---- facets_mode_methods.R: plot error path (line 157) ----

test_that("plot.mfrm_facets_run rejects non-mfrm_facets_run input", {
  expect_error(
    plot.mfrm_facets_run(list()),
    "mfrm_facets_run"
  )
})

# ---- build_apa_text with bias (line 876) ----

test_that("build_apa_table_figure_note_map includes bias note when provided", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task")

  note_map <- mfrmr:::build_apa_table_figure_note_map(
    res = fit, diagnostics = diag, bias_results = bias
  )
  expect_true(grepl("Bias", note_map$table4))
})

# ---- build_bias_fixed_text with facet_a only (line 159) ----

test_that("build_bias_fixed_text handles single facet label", {
  bfrt <- mfrmr:::build_bias_fixed_text
  result <- bfrt(
    table_df = data.frame(A = 1:2, B = 3:4),
    summary_df = NULL,
    chi_df = NULL,
    facet_a = "Rater",
    facet_b = NULL,
    interaction_label = NULL,
    columns = c("A", "B"),
    formats = list()
  )
  expect_true(grepl("Rater", result))
  expect_true(is.character(result))
})

# ---- PCM model APA text (line 653) ----

test_that("build_apa_text includes PCM step structure text", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             model = "PCM", method = "JML", step_facet = "Criterion", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  apa <- build_apa_outputs(fit, diagnostics = diag)
  # Should mention PCM and step structure
  expect_true(grepl("PCM", apa$report_text))
})

# ---- Weighted model APA text (line 658) ----

test_that("build_apa_text includes weight text", {
  set.seed(42)
  d <- data.frame(
    Person = rep(paste0("P", 1:8), each = 3),
    Rater = rep(paste0("R", 1:3), 8),
    Score = sample(0:2, 24, replace = TRUE),
    W = rep(c(1, 2, 0.5), 8)
  )
  fit <- suppressWarnings(
    fit_mfrm(d, "Person", "Rater", "Score", weight = "W",
             method = "JML", maxit = 20)
  )
  diag <- diagnose_mfrm(fit)
  apa <- build_apa_outputs(fit, diagnostics = diag)
  expect_true(grepl("[Ww]eight", apa$report_text))
})

# ---- build_pca_reference_text (line 329) ----

test_that("build_pca_reference_text produces formatted reference text", {
  refs <- mfrmr:::warning_threshold_profiles()
  bprt <- mfrmr:::build_pca_reference_text
  result <- bprt(refs$pca_reference_bands)
  expect_true(grepl("Literature bands", result))
  expect_true(grepl("EV >=", result))
})
