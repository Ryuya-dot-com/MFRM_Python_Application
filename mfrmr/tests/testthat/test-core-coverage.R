# test-core-coverage.R
# Exercises uncovered code paths in mfrm_core.R including:
# - Math helpers: logsumexp, weighted_mean, gauss_hermite_normal, center_sum_zero
# - Facet constraint building and expansion
# - zstd_from_mnsq (both exact and Wilson-Hilferty)
# - safe_cor / weighted_mean_safe / infer_default_rater_facet
# - calc_step_order / category_warnings_text
# - get_extreme_levels
# - make_union_find / calc_subsets
# - normalize_anchor_df / normalize_group_anchor_df
# - read_flexible_table
# - calc_reliability / calc_facets_chisq
# - extract_anchor_tables
# - resolve_pcm_step_facet / sanitize_noncenter_facet / sanitize_dummy_facets
# - build_facet_signs
# - Various compute/estimate helpers

# ---- Shared fixture ----

local({
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  .fit <<- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 20)
  )
  .diag <<- diagnose_mfrm(.fit, residual_pca = "both", pca_max_factors = 3)
})

# ============================================================================
# Math helpers
# ============================================================================

test_that("logsumexp is numerically stable", {
  lse <- mfrmr:::logsumexp
  expect_equal(lse(c(0, 0)), log(2))
  # Large values should not overflow
  result <- lse(c(1000, 1001))
  expect_true(is.finite(result))
  expect_true(result > 1000)
})

test_that("weighted_mean handles edge cases", {
  wm <- mfrmr:::weighted_mean
  expect_equal(wm(c(1, 2, 3), c(1, 1, 1)), 2)
  expect_true(is.na(wm(c(NA, NA), c(1, 1))))
  expect_true(is.na(wm(c(1, 2), c(0, 0))))
  expect_true(is.na(wm(c(1, 2), c(-1, -1))))
  # Weighted case
  expect_equal(wm(c(0, 10), c(1, 1)), 5)
})

test_that("get_weights returns weights from data frame", {
  gw <- mfrmr:::get_weights
  df_with_weight <- data.frame(X = 1:3, Weight = c(2, 0, NA))
  w <- gw(df_with_weight)
  expect_equal(w, c(2, 0, 0))
  df_no_weight <- data.frame(X = 1:3)
  w2 <- gw(df_no_weight)
  expect_equal(w2, c(1, 1, 1))
})

test_that("gauss_hermite_normal returns correct output for n=1", {
  gh <- mfrmr:::gauss_hermite_normal
  result <- gh(1)
  expect_equal(result$nodes, 0)
  expect_equal(result$weights, 1)
})

test_that("gauss_hermite_normal returns correct output for n > 1", {
  gh <- mfrmr:::gauss_hermite_normal
  result <- gh(5)
  expect_length(result$nodes, 5)
  expect_length(result$weights, 5)
  # Weights should sum to approximately 1
  expect_equal(sum(result$weights), 1, tolerance = 1e-10)
})

test_that("gauss_hermite_normal errors for n < 1", {
  expect_error(mfrmr:::gauss_hermite_normal(0), "n >= 1")
})

test_that("center_sum_zero centers values", {
  csz <- mfrmr:::center_sum_zero
  result <- csz(c(1, 2, 3))
  expect_equal(mean(result), 0)
  expect_equal(csz(numeric(0)), numeric(0))
})

test_that("expand_facet works for various n_levels", {
  ef <- mfrmr:::expand_facet
  expect_equal(ef(c(), 1), 0)
  result <- ef(c(0.5), 2)
  expect_equal(result, c(0.5, -0.5))
})

# ============================================================================
# zstd_from_mnsq
# ============================================================================

test_that("zstd_from_mnsq works with Wilson-Hilferty approximation", {
  zstd <- mfrmr:::zstd_from_mnsq
  result <- zstd(1.0, 100, whexact = FALSE)
  expect_true(is.finite(result))
  # MnSq = 1 should give Z near 0
  expect_true(abs(result) < 0.5)
})

test_that("zstd_from_mnsq works with exact approximation", {
  zstd <- mfrmr:::zstd_from_mnsq
  result <- zstd(1.0, 100, whexact = TRUE)
  expect_true(is.finite(result))
  expect_equal(result, 0)
})

test_that("zstd_from_mnsq handles vectorized input", {
  zstd <- mfrmr:::zstd_from_mnsq
  result <- zstd(c(0.5, 1.0, 1.5), 100, whexact = FALSE)
  expect_length(result, 3)
  expect_true(all(is.finite(result)))
})

test_that("zstd_from_mnsq handles edge cases", {
  zstd <- mfrmr:::zstd_from_mnsq
  expect_equal(length(zstd(numeric(0), numeric(0))), 0)
  result <- zstd(c(1, NA), c(100, 100))
  expect_true(is.finite(result[1]))
  expect_true(is.na(result[2]))
})

# ============================================================================
# safe_cor
# ============================================================================

test_that("safe_cor handles basic correlation", {
  sc <- mfrmr:::safe_cor
  r <- sc(c(1, 2, 3, 4), c(1, 2, 3, 4))
  expect_equal(r, 1, tolerance = 1e-10)
})

test_that("safe_cor handles constant vector", {
  sc <- mfrmr:::safe_cor
  r <- sc(c(1, 1, 1), c(1, 2, 3))
  expect_true(is.na(r))
})

test_that("safe_cor handles all NA", {
  sc <- mfrmr:::safe_cor
  r <- sc(c(NA, NA, NA), c(1, 2, 3))
  expect_true(is.na(r))
})

test_that("safe_cor with weights", {
  sc <- mfrmr:::safe_cor
  r <- sc(c(1, 2, 3, 4), c(1, 2, 3, 4), w = c(1, 1, 1, 1))
  expect_true(abs(r - 1) < 1e-10)
  # All zero weights
  r2 <- sc(c(1, 2, 3), c(1, 2, 3), w = c(0, 0, 0))
  expect_true(is.na(r2))
})

# ============================================================================
# infer_default_rater_facet
# ============================================================================

test_that("infer_default_rater_facet finds rater-like facets", {
  ifr <- mfrmr:::infer_default_rater_facet
  expect_equal(ifr(c("Task", "Rater", "Criterion")), "Rater")
  expect_equal(ifr(c("Task", "Judge", "Criterion")), "Judge")
  expect_equal(ifr(c("Task", "Grader", "Criterion")), "Grader")
  expect_equal(ifr(c("Task", "Reader", "Criterion")), "Reader")
  expect_equal(ifr(c("Task", "Scorer")), "Scorer")
  expect_equal(ifr(c("Task", "Assessor")), "Assessor")
  expect_equal(ifr(c("Task", "Evaluator")), "Evaluator")
  # No match: returns first
  expect_equal(ifr(c("Task", "Criterion")), "Task")
  expect_null(ifr(character(0)))
})

# ============================================================================
# calc_step_order
# ============================================================================

test_that("calc_step_order returns ordering info for step table", {
  cso <- mfrmr:::calc_step_order
  step_tbl <- .fit$steps
  result <- cso(step_tbl)
  expect_true(is.data.frame(result))
  expect_true("Ordered" %in% names(result))
  expect_true("Spacing" %in% names(result))
})

test_that("calc_step_order handles NULL input", {
  result <- mfrmr:::calc_step_order(NULL)
  expect_equal(nrow(result), 0)
})

# ============================================================================
# category_warnings_text
# ============================================================================

test_that("category_warnings_text returns warning messages", {
  cwt <- mfrmr:::category_warnings_text
  cat_tbl <- mfrmr:::calc_category_stats(.diag$obs, res = .fit)
  text <- cwt(cat_tbl)
  expect_true(is.character(text))
  expect_true(nchar(text) > 0)
})

test_that("category_warnings_text with step_tbl checks ordering", {
  cwt <- mfrmr:::category_warnings_text
  cat_tbl <- mfrmr:::calc_category_stats(.diag$obs, res = .fit)
  step_order <- mfrmr:::calc_step_order(.fit$steps)
  text <- cwt(cat_tbl, step_tbl = step_order)
  expect_true(is.character(text))
})

test_that("category_warnings_text with NULL returns no diagnostics message", {
  text <- mfrmr:::category_warnings_text(NULL)
  expect_true(grepl("No category diagnostics", text))
})

# ============================================================================
# get_extreme_levels
# ============================================================================

test_that("get_extreme_levels identifies extreme scorers", {
  gel <- mfrmr:::get_extreme_levels
  result <- gel(.diag$obs, .fit$config$facet_names,
                .fit$prep$rating_min, .fit$prep$rating_max)
  expect_true(is.list(result))
  expect_true(all(.fit$config$facet_names %in% names(result)))
})

# ============================================================================
# make_union_find and calc_subsets
# ============================================================================

test_that("make_union_find creates a working union-find structure", {
  uf <- mfrmr:::make_union_find(c("A", "B", "C", "D"))
  expect_equal(uf$find("A"), "A")
  uf$union("A", "B")
  expect_equal(uf$find("A"), uf$find("B"))
  uf$union("C", "D")
  expect_equal(uf$find("C"), uf$find("D"))
  # A and C should be in different sets
  expect_false(uf$find("A") == uf$find("C"))
  uf$union("A", "C")
  expect_equal(uf$find("A"), uf$find("D"))
})

test_that("calc_subsets returns subset info from observation data", {
  result <- mfrmr:::calc_subsets(.diag$obs, .fit$config$facet_names)
  expect_true(is.list(result))
  expect_true("summary" %in% names(result))
  expect_true("nodes" %in% names(result))
  expect_true(is.data.frame(result$summary))
  expect_true(is.data.frame(result$nodes))
})

test_that("calc_subsets handles empty input", {
  result <- mfrmr:::calc_subsets(data.frame(), character(0))
  expect_equal(nrow(result$summary), 0)
})

# ============================================================================
# normalize_anchor_df / normalize_group_anchor_df
# ============================================================================

test_that("normalize_anchor_df normalizes a valid anchor table", {
  naf <- mfrmr:::normalize_anchor_df
  df <- data.frame(Facet = "Rater", Level = "R1", Anchor = 0.5, stringsAsFactors = FALSE)
  result <- naf(df)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
  expect_equal(result$Anchor, 0.5)
})

test_that("normalize_anchor_df handles NULL and empty", {
  naf <- mfrmr:::normalize_anchor_df
  result_null <- naf(NULL)
  expect_equal(nrow(result_null), 0)
  result_empty <- naf(data.frame())
  expect_equal(nrow(result_empty), 0)
})

test_that("normalize_anchor_df handles missing columns", {
  naf <- mfrmr:::normalize_anchor_df
  df <- data.frame(X = 1:3, Y = 4:6)
  result <- naf(df)
  expect_equal(nrow(result), 0)
})

test_that("normalize_group_anchor_df normalizes a valid group anchor table", {
  ngaf <- mfrmr:::normalize_group_anchor_df
  df <- data.frame(
    Facet = "Rater", Level = "R1", Group = "G1", GroupValue = 0.0,
    stringsAsFactors = FALSE
  )
  result <- ngaf(df)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("normalize_group_anchor_df handles NULL and missing columns", {
  ngaf <- mfrmr:::normalize_group_anchor_df
  expect_equal(nrow(ngaf(NULL)), 0)
  expect_equal(nrow(ngaf(data.frame(X = 1:3))), 0)
})

# ============================================================================
# read_flexible_table
# ============================================================================

test_that("read_flexible_table reads tab-separated text", {
  rft <- mfrmr:::read_flexible_table
  text <- "A\tB\n1\t2\n3\t4"
  result <- rft(text, NULL)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  expect_true(all(c("A", "B") %in% names(result)))
})

test_that("read_flexible_table reads comma-separated text", {
  rft <- mfrmr:::read_flexible_table
  text <- "A,B\n1,2\n3,4"
  result <- rft(text, NULL)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
})

test_that("read_flexible_table reads semicolon-separated text", {
  rft <- mfrmr:::read_flexible_table
  text <- "A;B\n1;2\n3;4"
  result <- rft(text, NULL)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
})

test_that("read_flexible_table handles empty and NULL input", {
  rft <- mfrmr:::read_flexible_table
  result_null <- rft(NULL, NULL)
  expect_equal(nrow(result_null), 0)
  result_empty <- rft("", NULL)
  expect_equal(nrow(result_empty), 0)
  result_spaces <- rft("   ", NULL)
  expect_equal(nrow(result_spaces), 0)
})

# ============================================================================
# resolve_pcm_step_facet / sanitize_noncenter_facet / sanitize_dummy_facets / build_facet_signs
# ============================================================================

test_that("resolve_pcm_step_facet resolves correctly", {
  rpsf <- mfrmr:::resolve_pcm_step_facet
  expect_null(rpsf("RSM", NULL, c("Rater", "Task")))
  expect_equal(rpsf("PCM", NULL, c("Rater", "Task")), "Rater")
  expect_equal(rpsf("PCM", "Task", c("Rater", "Task")), "Task")
  expect_error(rpsf("PCM", "Missing", c("Rater", "Task")), "not among")
})

test_that("sanitize_noncenter_facet returns valid facet or Person", {
  snf <- mfrmr:::sanitize_noncenter_facet
  expect_equal(snf("Rater", c("Rater", "Task")), "Rater")
  expect_equal(snf("Person", c("Rater", "Task")), "Person")
  expect_equal(snf("Missing", c("Rater", "Task")), "Person")
  expect_equal(snf(NULL, c("Rater", "Task")), "Person")
})

test_that("sanitize_dummy_facets returns valid intersections", {
  sdf <- mfrmr:::sanitize_dummy_facets
  expect_equal(sdf(NULL, c("Rater", "Task")), character(0))
  expect_equal(sdf(c("Rater", "Missing"), c("Rater", "Task")), "Rater")
  expect_equal(sdf(c("Person"), c("Rater", "Task")), "Person")
})

test_that("build_facet_signs assigns correct signs", {
  bfs <- mfrmr:::build_facet_signs
  result <- bfs(c("Rater", "Task"), positive_facets = "Task")
  expect_equal(result$signs[["Rater"]], -1)
  expect_equal(result$signs[["Task"]], 1)
  expect_equal(result$positive_facets, "Task")
})

test_that("build_facet_signs with no positives defaults all to -1", {
  bfs <- mfrmr:::build_facet_signs
  result <- bfs(c("Rater", "Task"))
  expect_equal(result$signs[["Rater"]], -1)
  expect_equal(result$signs[["Task"]], -1)
})

# ============================================================================
# calc_reliability
# ============================================================================

test_that("calc_reliability produces reliability/separation by facet", {
  # Build a measure table mimicking diagnostics$measures
  measures <- .diag$measures
  if (!is.null(measures) && nrow(measures) > 0) {
    result <- mfrmr:::calc_reliability(measures)
    expect_true(is.data.frame(result))
    expect_true("Reliability" %in% names(result))
    expect_true("Separation" %in% names(result))
    expect_true("Strata" %in% names(result))
    expect_true(nrow(result) > 0)
  }
})

# ============================================================================
# calc_facets_chisq
# ============================================================================

test_that("calc_facets_chisq produces chi-square statistics", {
  measures <- .diag$measures
  if (!is.null(measures) && nrow(measures) > 0) {
    result <- mfrmr:::calc_facets_chisq(measures)
    expect_true(is.data.frame(result))
    expect_true("FixedChiSq" %in% names(result))
    expect_true("FixedProb" %in% names(result))
    expect_true(nrow(result) > 0)
  }
})

test_that("calc_facets_chisq handles NULL input", {
  result <- mfrmr:::calc_facets_chisq(NULL)
  expect_equal(nrow(result), 0)
})

# ============================================================================
# extract_anchor_tables
# ============================================================================

test_that("extract_anchor_tables returns anchors and groups", {
  result <- mfrmr:::extract_anchor_tables(.fit$config)
  expect_true(is.list(result))
  expect_true("anchors" %in% names(result))
  expect_true("groups" %in% names(result))
  expect_true(is.data.frame(result$anchors))
  expect_true(is.data.frame(result$groups))
})

# ============================================================================
# build_facet_constraint / count_facet_params / expand_facet_with_constraints
# ============================================================================

test_that("build_facet_constraint with anchors", {
  bfc <- mfrmr:::build_facet_constraint
  result <- bfc(
    levels = c("R1", "R2", "R3"),
    anchors = c(R1 = 0.5),
    centered = TRUE
  )
  expect_true(is.list(result))
  expect_equal(result$levels, c("R1", "R2", "R3"))
  expect_true(!is.na(result$anchors["R1"]))
  expect_true(is.na(result$anchors["R2"]))
})

test_that("build_facet_constraint with groups", {
  bfc <- mfrmr:::build_facet_constraint
  result <- bfc(
    levels = c("R1", "R2", "R3"),
    groups = c(R1 = "G1", R2 = "G1"),
    group_values = c(G1 = 0.0),
    centered = FALSE
  )
  expect_true(is.list(result))
  expect_equal(result$groups["R1"], c(R1 = "G1"))
})

test_that("count_facet_params counts free parameters", {
  cfp <- mfrmr:::count_facet_params
  spec <- mfrmr:::build_facet_constraint(
    levels = c("R1", "R2", "R3"),
    anchors = c(R1 = 0.5),
    centered = TRUE
  )
  n <- cfp(spec)
  expect_true(is.integer(n) || is.numeric(n))
  expect_true(n >= 0)
})

# ============================================================================
# guess_col / truncate_label / facet_report_id
# ============================================================================

test_that("guess_col finds matching column", {
  gc <- mfrmr:::guess_col
  cols <- c("Person_ID", "Rater_Score", "Task_Name")
  expect_equal(gc(cols, c("rater")), "Rater_Score")
  expect_equal(gc(cols, c("nonexistent")), cols[1])
  expect_equal(gc(character(0), "rater"), character(0))
})

test_that("truncate_label truncates long strings", {
  tl <- mfrmr:::truncate_label
  expect_equal(tl("short", 28), "short")
  long <- paste(rep("x", 50), collapse = "")
  result <- tl(long, 10)
  expect_true(nchar(result) <= 13)  # allows for "..."
})

test_that("facet_report_id generates report ID", {
  fri <- mfrmr:::facet_report_id
  result <- fri("Rater")
  expect_true(is.character(result))
  expect_true(nchar(result) > 0)
})

# ============================================================================
# category probability functions
# ============================================================================

test_that("category_prob_rsm returns valid probabilities", {
  cpr <- mfrmr:::category_prob_rsm
  step_cum <- c(0, -0.5, 0.5)
  probs <- cpr(eta = 0, step_cum = step_cum)
  expect_true(all(probs >= 0))
  expect_equal(sum(probs), 1, tolerance = 1e-10)
})

test_that("category_prob_pcm returns valid probabilities", {
  cpp <- mfrmr:::category_prob_pcm
  step_cum_mat <- matrix(c(0, -0.5, 0.5, 0, -0.3, 0.3), nrow = 2, byrow = TRUE)
  probs <- cpp(eta = 0, step_cum_mat = step_cum_mat, criterion_idx = 1)
  expect_true(all(probs >= 0))
  expect_equal(sum(probs), 1, tolerance = 1e-10)
})

# ============================================================================
# expected_score_table / compute_obs_table
# ============================================================================

test_that("expected_score_table produces a table from fit", {
  est <- mfrmr:::expected_score_table(.fit)
  expect_true(is.data.frame(est))
  expect_true("Expected" %in% names(est))
})

test_that("compute_obs_table produces an observation table from fit", {
  obs_tbl <- mfrmr:::compute_obs_table(.fit)
  expect_true(is.data.frame(obs_tbl))
  expect_true(all(c("Observed", "Expected") %in% names(obs_tbl)))
})

# ============================================================================
# compute_prob_matrix / compute_scorefile / compute_residual_file
# ============================================================================

test_that("compute_prob_matrix returns a matrix-like structure", {
  pm <- mfrmr:::compute_prob_matrix(.fit)
  expect_true(!is.null(pm))
})

test_that("compute_scorefile produces a data frame", {
  sf <- mfrmr:::compute_scorefile(.fit)
  expect_true(is.data.frame(sf))
  expect_true(nrow(sf) > 0)
})

test_that("compute_residual_file produces a data frame", {
  rf <- mfrmr:::compute_residual_file(.fit)
  expect_true(is.data.frame(rf))
  expect_true(nrow(rf) > 0)
})

# ============================================================================
# PCM model test (ensure different code path)
# ============================================================================

test_that("fit_mfrm works with PCM model", {
  d <- mfrmr:::sample_mfrm_data(seed = 99)
  fit_pcm <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", model = "PCM", step_facet = "Criterion",
             maxit = 15)
  )
  expect_s3_class(fit_pcm, "mfrm_fit")
  expect_equal(fit_pcm$config$model, "PCM")
  # Check PCM-specific output
  expect_true(nrow(fit_pcm$steps) > 0)
})

# ============================================================================
# MML method test (ensure MML code path)
# ============================================================================

test_that("fit_mfrm works with MML method", {
  d <- mfrmr:::sample_mfrm_data(seed = 88)
  fit_mml <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "MML", maxit = 15, quad_points = 5)
  )
  expect_s3_class(fit_mml, "mfrm_fit")
  expect_equal(fit_mml$config$method, "MML")
})

# ============================================================================
# ensure_positive_definite
# ============================================================================

test_that("ensure_positive_definite handles already positive definite matrix", {
  epd <- mfrmr:::ensure_positive_definite
  mat <- diag(3)
  result <- epd(mat)
  expect_equal(result, mat)
})

test_that("ensure_positive_definite attempts smoothing for non-PD matrix", {
  epd <- mfrmr:::ensure_positive_definite
  # Create a non-positive-definite matrix
  mat <- matrix(c(1, 2, 2, 1), nrow = 2)
  result <- epd(mat)
  expect_true(is.matrix(result))
})

# ============================================================================
# weighted_mean_safe
# ============================================================================

test_that("weighted_mean_safe delegates to weighted_mean", {
  wms <- mfrmr:::weighted_mean_safe
  expect_equal(wms(c(1, 2, 3), c(1, 1, 1)), 2)
})

# ============================================================================
# Edge case: positive_facets in fit_mfrm
# ============================================================================

test_that("fit_mfrm with positive_facets changes sign interpretation", {
  d <- mfrmr:::sample_mfrm_data(seed = 77)
  fit_pos <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 15, positive_facets = "Criterion")
  )
  expect_s3_class(fit_pos, "mfrm_fit")
  expect_true("Criterion" %in% fit_pos$config$positive_facets)
})

# ============================================================================
# Edge case: noncenter_facet
# ============================================================================

test_that("fit_mfrm with noncenter_facet adjusts centering", {
  d <- mfrmr:::sample_mfrm_data(seed = 66)
  fit_nc <- suppressWarnings(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
             method = "JML", maxit = 15, noncenter_facet = "Rater")
  )
  expect_s3_class(fit_nc, "mfrm_fit")
  expect_equal(fit_nc$config$noncenter_facet, "Rater")
})

# ============================================================================
# calc_bias_pairwise
# ============================================================================

test_that("calc_bias_pairwise produces pairwise comparisons", {
  bias_tbl <- .diag$bias_interactions
  if (!is.null(bias_tbl) && nrow(bias_tbl) > 0) {
    result <- mfrmr:::calc_bias_pairwise(bias_tbl, "Rater", "Task")
    expect_true(is.data.frame(result))
  }
})

# ============================================================================
# compute_bias_adjustment_vector / compute_prob_matrix_with_bias
# ============================================================================

test_that("compute_bias_adjustment_vector returns adjustments with bias results", {
  bias_results <- estimate_bias(.fit, .diag, facet_a = "Rater", facet_b = "Task")
  adj <- mfrmr:::compute_bias_adjustment_vector(.fit, bias_results)
  expect_true(is.numeric(adj))
  expect_equal(length(adj), nrow(.diag$obs))
})

test_that("compute_bias_adjustment_vector returns zeros without bias", {
  adj <- mfrmr:::compute_bias_adjustment_vector(.fit, NULL)
  expect_true(is.numeric(adj))
  expect_true(all(adj == 0))
})

test_that("compute_prob_matrix_with_bias produces output", {
  bias_results <- estimate_bias(.fit, .diag, facet_a = "Rater", facet_b = "Task")
  pm <- mfrmr:::compute_prob_matrix_with_bias(.fit, bias_results)
  expect_true(!is.null(pm))
})
