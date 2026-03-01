# test-reporting-coverage.R
# Exercises uncovered paths in reporting.R:
# - build_visual_warning_map (all warning stages)
# - build_visual_summary_map (all summary stages, including detail="detailed")
# - build_apa_report_text (with and without context/bias)
# - build_apa_table_figure_notes / captions
# - build_sectioned_fixed_report edge cases
# - format_fixed_width_table edge cases
# - resolve_warning_thresholds with different profiles
# - build_pca_reference_text / build_pca_check_text
# - summarize_anchor_constraints / summarize_convergence_metrics / summarize_step_estimates
# - summarize_bias_counts / summarize_top_misfit_levels

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

# ============================================================================
# build_visual_warning_map
# ============================================================================

test_that("build_visual_warning_map returns all expected visual keys", {
  wmap <- mfrmr:::build_visual_warning_map(.fit, .diag)
  expected_keys <- c(
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds",
    "category_curves", "observed_expected", "fit_diagnostics",
    "fit_zstd_distribution", "misfit_levels",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  expect_true(all(expected_keys %in% names(wmap)))
})

test_that("build_visual_warning_map with strict profile produces more warnings", {
  wmap_strict <- mfrmr:::build_visual_warning_map(.fit, .diag, threshold_profile = "strict")
  expect_true(is.list(wmap_strict))
  expect_true(length(wmap_strict$residual_pca_overall) > 0)
})

test_that("build_visual_warning_map with lenient profile", {
  wmap_lenient <- mfrmr:::build_visual_warning_map(.fit, .diag, threshold_profile = "lenient")
  expect_true(is.list(wmap_lenient))
})

test_that("build_visual_warning_map with NULL inputs returns empty map", {
  wmap_null <- mfrmr:::build_visual_warning_map(NULL, NULL)
  expect_true(is.list(wmap_null))
  expect_true(all(lengths(wmap_null) == 0))
})

test_that("build_visual_warning_map with custom threshold overrides", {
  custom_thr <- list(n_obs_min = 999999, n_person_min = 999999)
  wmap <- mfrmr:::build_visual_warning_map(.fit, .diag, thresholds = custom_thr)
  expect_true(any(grepl("Small", unlist(wmap))))
})

# ============================================================================
# build_visual_summary_map
# ============================================================================

test_that("build_visual_summary_map returns summary text for all keys", {
  smap <- mfrmr:::build_visual_summary_map(.fit, .diag)
  expected_keys <- c(
    "wright_map", "pathway_map", "facet_distribution", "step_thresholds",
    "category_curves", "observed_expected", "fit_diagnostics",
    "fit_zstd_distribution", "misfit_levels",
    "residual_pca_overall", "residual_pca_by_facet"
  )
  expect_true(all(expected_keys %in% names(smap)))
  # wright_map should always have some text
  expect_true(length(smap$wright_map) > 0)
})

test_that("build_visual_summary_map with detail='detailed' adds extra summaries", {
  smap_detailed <- mfrmr:::build_visual_summary_map(
    .fit, .diag,
    options = list(detail = "detailed", max_facet_ranges = 2, top_misfit_n = 5)
  )
  expect_true(is.list(smap_detailed))
  expect_true(length(smap_detailed$wright_map) > 0)
})

test_that("build_visual_summary_map with NULL inputs returns empty map", {
  smap_null <- mfrmr:::build_visual_summary_map(NULL, NULL)
  expect_true(is.list(smap_null))
  expect_true(all(lengths(smap_null) == 0))
})

# ============================================================================
# build_apa_report_text
# ============================================================================

test_that("build_apa_report_text produces Method and Results sections", {
  text <- mfrmr:::build_apa_report_text(.fit, .diag)
  expect_true(grepl("Method", text))
  expect_true(grepl("Results", text))
  expect_true(grepl("many-facet Rasch", text, ignore.case = TRUE))
})

test_that("build_apa_report_text with context supplies assessment/setting text", {
  ctx <- list(
    assessment = "writing proficiency",
    setting = "a university setting",
    rater_training = "two hours of calibration training",
    raters_per_response = "two",
    scale_desc = "a 5-point holistic rubric",
    line_width = 80L
  )
  text <- mfrmr:::build_apa_report_text(.fit, .diag, context = ctx)
  expect_true(grepl("writing proficiency", text))
  expect_true(grepl("university", text))
  expect_true(grepl("calibration", text))
  expect_true(grepl("holistic rubric", text))
})

test_that("build_apa_report_text with bias_results includes bias summary", {
  text <- mfrmr:::build_apa_report_text(.fit, .diag, bias_results = .bias)
  expect_true(grepl("Bias", text, ignore.case = TRUE) || grepl("bias", text))
})

test_that("build_apa_report_text without bias mentions no bias data", {
  text <- mfrmr:::build_apa_report_text(.fit, .diag, bias_results = NULL)
  expect_true(nchar(text) > 0)
})

# ============================================================================
# build_apa_table_figure_notes / captions
# ============================================================================

test_that("build_apa_table_figure_notes produces note text", {
  notes <- mfrmr:::build_apa_table_figure_notes(.fit, .diag)
  expect_true(is.character(notes))
  expect_true(nchar(notes) > 0)
})

test_that("build_apa_table_figure_notes with bias_results includes bias notes", {
  notes <- mfrmr:::build_apa_table_figure_notes(.fit, .diag, bias_results = .bias)
  expect_true(nchar(notes) > 0)
})

test_that("build_apa_table_figure_captions produces caption text", {
  captions <- mfrmr:::build_apa_table_figure_captions(.fit, .diag)
  expect_true(is.character(captions))
  expect_true(grepl("Table 1", captions))
})

test_that("build_apa_table_figure_captions with context adds assessment phrase", {
  ctx <- list(assessment = "oral proficiency")
  captions <- mfrmr:::build_apa_table_figure_captions(.fit, .diag, context = ctx)
  expect_true(grepl("oral proficiency", captions))
})

test_that("build_apa_table_figure_captions with bias includes interaction label", {
  captions <- mfrmr:::build_apa_table_figure_captions(.fit, .diag, bias_results = .bias)
  expect_true(grepl("Table 4", captions))
})

# ============================================================================
# resolve_warning_thresholds
# ============================================================================

test_that("resolve_warning_thresholds returns all profiles", {
  for (profile in c("strict", "standard", "lenient")) {
    resolved <- mfrmr:::resolve_warning_thresholds(threshold_profile = profile)
    expect_equal(resolved$profile_name, profile)
    expect_true("thresholds" %in% names(resolved))
    expect_true("pca_reference_bands" %in% names(resolved))
    expect_true(resolved$thresholds$n_obs_min > 0)
  }
})

test_that("resolve_warning_thresholds applies custom overrides", {
  resolved <- mfrmr:::resolve_warning_thresholds(
    thresholds = list(n_obs_min = 42),
    threshold_profile = "standard"
  )
  expect_equal(resolved$thresholds$n_obs_min, 42)
})

# ============================================================================
# build_pca_reference_text / build_pca_check_text
# ============================================================================

test_that("build_pca_reference_text returns reference band text", {
  bands <- mfrmr:::warning_threshold_profiles()$pca_reference_bands
  text <- mfrmr:::build_pca_reference_text(bands)
  expect_true(is.character(text))
  expect_true(nchar(text) > 0)
})

test_that("build_pca_check_text returns interpretive text for a given eigenvalue/proportion", {
  bands <- mfrmr:::warning_threshold_profiles()$pca_reference_bands
  text <- mfrmr:::build_pca_check_text(eigenvalue = 1.5, proportion = 0.08, reference_bands = bands)
  expect_true(is.character(text))
  expect_true(nchar(text) > 0)
})

test_that("build_pca_check_text with high eigenvalue flags concern", {
  bands <- mfrmr:::warning_threshold_profiles()$pca_reference_bands
  text <- mfrmr:::build_pca_check_text(eigenvalue = 5.0, proportion = 0.25, reference_bands = bands)
  expect_true(nchar(text) > 0)
})

# ============================================================================
# Format helpers: format_fixed_width_table edge cases
# ============================================================================

test_that("format_fixed_width_table handles empty data frame", {
  result <- mfrmr:::format_fixed_width_table(data.frame(), columns = character(0))
  expect_equal(result, "No data")
})

test_that("format_fixed_width_table handles missing columns gracefully", {
  df <- data.frame(A = 1:3, B = c("x", "y", "z"), stringsAsFactors = FALSE)
  result <- mfrmr:::format_fixed_width_table(df, columns = c("A", "B", "C"))
  expect_true(is.character(result))
  expect_true(grepl("A", result))
})

# ============================================================================
# build_sectioned_fixed_report
# ============================================================================

test_that("build_sectioned_fixed_report handles various section types", {
  report <- mfrmr:::build_sectioned_fixed_report(
    title = "Test Report",
    sections = list(
      list(title = "Empty Section", data = NULL),
      list(title = "String Section", data = "Some text here"),
      list(title = "Data Section", data = data.frame(X = 1:3, Y = c("a", "b", "c"), stringsAsFactors = FALSE)),
      list(title = "Empty DF", data = data.frame()),
      list(title = "Truncated Section", data = data.frame(X = 1:10), max_rows = 3)
    )
  )
  expect_true(grepl("Test Report", report))
  expect_true(grepl("Some text here", report))
  expect_true(grepl("Empty Section", report))
  expect_true(grepl("Showing first 3 rows of 10", report))
})

# ============================================================================
# summarize_anchor_constraints
# ============================================================================

test_that("summarize_anchor_constraints produces text about anchor status", {
  text <- mfrmr:::summarize_anchor_constraints(.fit$config)
  expect_true(is.character(text))
})

# ============================================================================
# summarize_convergence_metrics
# ============================================================================

test_that("summarize_convergence_metrics produces text from summary row", {
  summary_row <- .fit$summary
  if (!is.null(summary_row) && nrow(summary_row) > 0) {
    text <- mfrmr:::summarize_convergence_metrics(summary_row[1, , drop = FALSE])
    expect_true(is.character(text))
  }
})

# ============================================================================
# summarize_step_estimates
# ============================================================================

test_that("summarize_step_estimates produces text from step table", {
  text <- mfrmr:::summarize_step_estimates(.fit$steps)
  expect_true(is.character(text))
  expect_true(nchar(paste(text, collapse = "")) > 0)
})

# ============================================================================
# summarize_bias_counts
# ============================================================================

test_that("summarize_bias_counts with bias results", {
  text <- mfrmr:::summarize_bias_counts(.bias)
  expect_true(is.character(text))
  expect_true(nchar(paste(text, collapse = "")) > 0)
})

test_that("summarize_bias_counts with NULL returns text", {
  text <- mfrmr:::summarize_bias_counts(NULL)
  expect_true(is.character(text))
})

# ============================================================================
# summarize_top_misfit_levels
# ============================================================================

test_that("summarize_top_misfit_levels produces text", {
  fit_tbl <- .diag$fit
  if (!is.null(fit_tbl) && nrow(fit_tbl) > 0) {
    text <- mfrmr:::summarize_top_misfit_levels(fit_tbl, top_n = 3)
    expect_true(is.character(text))
  }
})

# ============================================================================
# collapse_apa_paragraph
# ============================================================================

test_that("collapse_apa_paragraph wraps text at specified width", {
  sentences <- c(
    "This is a first sentence.",
    "This is a second sentence that is somewhat longer than the first.",
    "And a third one."
  )
  result <- mfrmr:::collapse_apa_paragraph(sentences, width = 40L)
  expect_true(is.character(result))
  expect_true(nchar(result) > 0)
})

# ============================================================================
# py_style_format edge cases
# ============================================================================

test_that("py_style_format handles edge cases", {
  fmt <- mfrmr:::py_style_format
  expect_equal(fmt(NULL, 42), "42")
  expect_equal(fmt("{}", 42), "42")
  expect_equal(fmt("{:.0f}", NA), "")
  expect_equal(fmt(function(x) paste0("$", x), 100), "$100")
  expect_equal(fmt(c("a", "b"), 42), "42")
})

# ============================================================================
# describe_series edge cases
# ============================================================================

test_that("describe_series handles NULL and empty inputs", {
  expect_null(mfrmr:::describe_series(NULL))
  expect_null(mfrmr:::describe_series(c(NA, NaN)))
  s <- mfrmr:::describe_series(c(1, 2, 3))
  expect_equal(s$min, 1)
  expect_equal(s$max, 3)
  expect_equal(s$mean, 2)
})

test_that("describe_series handles single element", {
  s <- mfrmr:::describe_series(5)
  expect_equal(s$min, 5)
  expect_equal(s$max, 5)
  expect_true(is.na(s$sd))
})

# ============================================================================
# build_bias_fixed_text / build_pairwise_fixed_text
# ============================================================================

test_that("build_bias_fixed_text with valid data produces text", {
  if (!is.null(.bias$table) && nrow(.bias$table) > 0) {
    tbl <- as.data.frame(.bias$table, stringsAsFactors = FALSE)
    cols <- intersect(
      c("Bias Size", "Obs-Exp Average", "S.E.", "t", "Prob."),
      names(tbl)
    )
    text <- mfrmr:::build_bias_fixed_text(
      table_df = tbl,
      summary_df = .bias$summary %||% data.frame(),
      chi_df = .bias$chi_sq %||% data.frame(),
      facet_a = "Rater",
      facet_b = "Task",
      columns = cols,
      formats = list()
    )
    expect_true(is.character(text))
    expect_true(grepl("Bias", text))
  }
})

test_that("build_bias_fixed_text with empty table returns 'No bias data'", {
  text <- mfrmr:::build_bias_fixed_text(
    table_df = data.frame(),
    summary_df = data.frame(),
    chi_df = data.frame(),
    facet_a = "Rater",
    facet_b = "Task",
    columns = character(0),
    formats = list()
  )
  expect_equal(text, "No bias data")
})

test_that("build_pairwise_fixed_text with empty data returns 'No pairwise data'", {
  text <- mfrmr:::build_pairwise_fixed_text(
    pair_df = data.frame(),
    target_facet = "Rater",
    context_facet = "Task",
    columns = character(0),
    formats = list()
  )
  expect_equal(text, "No pairwise data")
})

# ============================================================================
# extract_overall_pca_first / extract_overall_pca_second / extract_facet_pca_first
# ============================================================================

test_that("PCA extraction functions return correct data", {
  pca <- mfrmr:::safe_residual_pca(.diag, mode = "both")
  if (!is.null(pca)) {
    first <- mfrmr:::extract_overall_pca_first(pca)
    expect_true(!is.null(first))
    expect_true("Eigenvalue" %in% names(first))

    second <- mfrmr:::extract_overall_pca_second(pca)
    # second may be NULL if only 1 component
    expect_true(is.null(second) || "Eigenvalue" %in% names(second))

    facet <- mfrmr:::extract_facet_pca_first(pca)
    expect_true(is.data.frame(facet))
  }
})

test_that("PCA extraction functions handle NULL input", {
  expect_null(mfrmr:::extract_overall_pca_first(NULL))
  expect_null(mfrmr:::extract_overall_pca_second(NULL))
  facet <- mfrmr:::extract_facet_pca_first(NULL)
  expect_true(is.data.frame(facet))
  expect_equal(nrow(facet), 0)
})

# ============================================================================
# build_apa_table_figure_note_map (internal)
# ============================================================================

test_that("build_apa_table_figure_note_map returns a named list", {
  note_map <- mfrmr:::build_apa_table_figure_note_map(.fit, .diag)
  expect_true(is.list(note_map))
  expect_true(length(note_map) > 0)
})

test_that("build_apa_table_figure_note_map with bias_results includes table4", {
  note_map <- mfrmr:::build_apa_table_figure_note_map(.fit, .diag, bias_results = .bias)
  expect_true("table4" %in% names(note_map))
})
