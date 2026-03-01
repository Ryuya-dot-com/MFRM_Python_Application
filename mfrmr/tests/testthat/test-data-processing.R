# test-data-processing.R
# Tests for the data processing layer -- the most important bug hotspot.
# Covers: facet ID conversion, factor vs character equivalence, label
# permutation invariance, missing data, score recoding, anchors, weights,
# and edge cases.

# ---------------------------------------------------------------------------
# 1.1  Facet ID conversion accuracy
# ---------------------------------------------------------------------------

test_that("character facet labels are preserved in fitted output", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  # Rater labels

  rater_levels <- sort(unique(d$Rater))
  fitted_rater <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::pull(Level) |>
    sort()
  expect_equal(fitted_rater, rater_levels)

  # Task labels
  task_levels <- sort(unique(d$Task))
  fitted_task <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::pull(Level) |>
    sort()
  expect_equal(fitted_task, task_levels)

  # Criterion labels
  crit_levels <- sort(unique(d$Criterion))
  fitted_crit <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::pull(Level) |>
    sort()
  expect_equal(fitted_crit, crit_levels)

  # Person labels preserved
  expect_setequal(fit$facets$person$Person, unique(d$Person))
})

test_that("facet level ordering follows alphabetical (factor default)", {

  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  rater_tbl <- fit$facets$others |>
    dplyr::filter(Facet == "Rater")
  expect_equal(rater_tbl$Level, sort(unique(d$Rater)))
})

# ---------------------------------------------------------------------------
# 1.2  Factor vs character equivalence
# ---------------------------------------------------------------------------

test_that("factor vs character facet columns produce identical estimates", {
  d <- mfrmr:::sample_mfrm_data(seed = 100)
  d_factor <- d |>
    dplyr::mutate(Rater = factor(Rater), Task = factor(Task))

  fit_char <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))
  fit_factor <- suppressWarnings(fit_mfrm(d_factor, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))

  expect_equal(fit_char$summary$LogLik, fit_factor$summary$LogLik, tolerance = 1e-6)

  for (facet in c("Rater", "Task", "Criterion")) {
    est_c <- fit_char$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    est_f <- fit_factor$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    expect_equal(est_c, est_f, tolerance = 1e-8,
      label = paste("factor vs char for", facet))
  }
})

test_that("factor with explicit level ordering matches character", {
  d <- mfrmr:::sample_mfrm_data(seed = 55)
  d_fac <- d |>
    dplyr::mutate(Rater = factor(Rater, levels = c("R3", "R1", "R2")))

  fit_chr <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))
  fit_fac <- suppressWarnings(fit_mfrm(d_fac, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))

  expect_equal(fit_chr$summary$LogLik, fit_fac$summary$LogLik, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 1.3  Integer ID vs string ID equivalence
# ---------------------------------------------------------------------------

test_that("integer facet IDs and string IDs produce equivalent estimates", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  rater_map <- c(R1 = "1", R2 = "2", R3 = "3")
  d_int <- d |> dplyr::mutate(Rater = rater_map[Rater])

  fit_str <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))
  fit_int <- suppressWarnings(fit_mfrm(d_int, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 40))

  expect_equal(fit_str$summary$LogLik, fit_int$summary$LogLik, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 1.4  Label permutation invariance (CRITICAL)
# ---------------------------------------------------------------------------

test_that("shuffling rater labels preserves estimates after remapping", {
  d <- mfrmr:::sample_mfrm_data(seed = 200)
  fit_orig <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  shuffle <- c(R1 = "RaterC", R2 = "RaterA", R3 = "RaterB")
  d_shuf <- d |> dplyr::mutate(Rater = shuffle[Rater])
  fit_shuf <- suppressWarnings(fit_mfrm(d_shuf, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  rev_map <- c(RaterC = "R1", RaterA = "R2", RaterB = "R3")
  est_orig <- fit_orig$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  est_remapped <- fit_shuf$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::mutate(OrigLevel = rev_map[Level]) |>
    dplyr::arrange(OrigLevel) |>
    dplyr::pull(Estimate)

  expect_equal(unname(est_orig), unname(est_remapped), tolerance = 0.01)
})

test_that("shuffling task labels preserves estimates after remapping", {
  d <- mfrmr:::sample_mfrm_data(seed = 200)
  fit_orig <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  shuffle <- c(T1 = "TaskD", T2 = "TaskB", T3 = "TaskA", T4 = "TaskC")
  d_shuf <- d |> dplyr::mutate(Task = shuffle[Task])
  fit_shuf <- suppressWarnings(fit_mfrm(d_shuf, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  rev_map <- c(TaskD = "T1", TaskB = "T2", TaskA = "T3", TaskC = "T4")
  est_orig <- fit_orig$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  est_remapped <- fit_shuf$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::mutate(OrigLevel = rev_map[Level]) |>
    dplyr::arrange(OrigLevel) |>
    dplyr::pull(Estimate)

  expect_equal(unname(est_orig), unname(est_remapped), tolerance = 0.01)
})

test_that("shuffling criterion labels preserves estimates after remapping", {
  d <- mfrmr:::sample_mfrm_data(seed = 200)
  fit_orig <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  shuffle <- c(C1 = "CritZ", C2 = "CritX", C3 = "CritY")
  d_shuf <- d |> dplyr::mutate(Criterion = shuffle[Criterion])
  fit_shuf <- suppressWarnings(fit_mfrm(d_shuf, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 100))

  rev_map <- c(CritZ = "C1", CritX = "C2", CritY = "C3")
  est_orig <- fit_orig$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  est_remapped <- fit_shuf$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::mutate(OrigLevel = rev_map[Level]) |>
    dplyr::arrange(OrigLevel) |>
    dplyr::pull(Estimate)

  expect_equal(unname(est_orig), unname(est_remapped), tolerance = 0.01)
})

# ---------------------------------------------------------------------------
# 1.5  Missing data handling
# ---------------------------------------------------------------------------

test_that("NA in person column: rows dropped, fit succeeds", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Person[1:5] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), nrow(d))
})

test_that("NA in score column: rows dropped, fit succeeds", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score[10:20] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), nrow(d))
})

test_that("5% random NA insertion: fit still succeeds, prep data rows < original", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  n_orig <- nrow(d)
  set.seed(999)
  na_idx <- sample(seq_len(n_orig), size = ceiling(0.05 * n_orig))
  d$Score[na_idx] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
    "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), n_orig)
})

test_that("all NA after filtering raises explicit error", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score <- NA
  expect_error(
    suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"),
      "Score", method = "JML", maxit = 30)),
    "No valid observations"
  )
})

# ---------------------------------------------------------------------------
# 1.6  Score category handling
# ---------------------------------------------------------------------------

test_that("unused intermediate category recoded to contiguous (keep_original=FALSE)", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Only keep scores 1, 3, 5 (drop 2 and 4)
  d_gap <- d |> dplyr::filter(Score %in% c(1, 3, 5))
  fit <- suppressWarnings(fit_mfrm(d_gap, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30, keep_original = FALSE))

  expect_s3_class(fit, "mfrm_fit")
  # After recoding, the data should have contiguous scores
  observed_scores <- sort(unique(fit$prep$data$Score))
  expected_contiguous <- seq(min(observed_scores), max(observed_scores))
  expect_equal(observed_scores, expected_contiguous)
})

test_that("keep_original=TRUE preserves original score codes", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Only keep scores 1, 3, 5
  d_gap <- d |> dplyr::filter(Score %in% c(1, 3, 5))
  fit <- suppressWarnings(fit_mfrm(d_gap, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30, keep_original = TRUE))

  expect_s3_class(fit, "mfrm_fit")
  observed_scores <- sort(unique(fit$prep$data$Score))
  expect_true(all(observed_scores %in% c(1, 3, 5)))
})

test_that("explicit rating_min/rating_max override auto-detection", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_auto <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30))
  fit_explicit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30,
    rating_min = 1, rating_max = 5))

  # With the same data and explicit range matching the data, results should match
  expect_equal(fit_auto$summary$LogLik, fit_explicit$summary$LogLik, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 1.7  Constraint (anchor) settings
# ---------------------------------------------------------------------------

test_that("anchor R2=0 forces Rater R2 estimate to zero", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchor_tbl <- data.frame(
    Facet = "Rater", Level = "R2", Anchor = 0,
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60,
    anchors = anchor_tbl))

  r2_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R2") |>
    dplyr::pull(Estimate)

  expect_equal(unname(r2_est), 0, tolerance = 1e-8)
})

test_that("invalid anchor (non-existent level) with warn policy triggers warning", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchor_tbl <- data.frame(
    Facet = "Rater", Level = "R99_NONEXISTENT", Anchor = 0,
    stringsAsFactors = FALSE
  )
  expect_warning(
    fit_mfrm(d, "Person",
      c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 30,
      anchors = anchor_tbl,
      anchor_policy = "warn"),
    regex = "."
  )
})

test_that("group anchors constrain group mean to target", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Group R1 and R3 together with target mean = 0
  group_tbl <- data.frame(
    Facet = c("Rater", "Rater"),
    Level = c("R1", "R3"),
    Group = c("G1", "G1"),
    GroupValue = c(0, 0),
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60,
    group_anchors = group_tbl))

  r1_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R1") |>
    dplyr::pull(Estimate)
  r3_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater", Level == "R3") |>
    dplyr::pull(Estimate)

  group_mean <- mean(c(r1_est, r3_est))
  expect_equal(group_mean, 0, tolerance = 1e-4)
})

# ---------------------------------------------------------------------------
# 1.8  Weight handling
# ---------------------------------------------------------------------------

test_that("zero-weight rows are excluded from estimation", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Weight <- 1
  d$Weight[1:50] <- 0

  fit_wt <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30, weight = "Weight"))

  expect_lt(nrow(fit_wt$prep$data), nrow(d))
})

test_that("explicit weight column produces different results from unit weights", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  set.seed(777)
  d$Weight <- sample(c(1, 2, 3), nrow(d), replace = TRUE)

  fit_unit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 40))
  fit_wt <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 40, weight = "Weight"))

  # LogLik should differ when weights are non-uniform
  expect_false(isTRUE(all.equal(fit_unit$summary$LogLik, fit_wt$summary$LogLik,
    tolerance = 1e-4)))
})

# ---------------------------------------------------------------------------
# 1.9  Edge cases
# ---------------------------------------------------------------------------

test_that("duplicate column names in data raise error with 'duplicate'", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d2 <- cbind(d, Rater = d$Rater)
  # This data.frame now has two columns named "Rater"
  expect_error(
    suppressWarnings(fit_mfrm(d2, "Person",
      c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 20)),
    regex = "[Dd]uplic"
  )
})

test_that("large number of facet levels (20+ persons) does not crash", {
  set.seed(1234)
  n_person <- 25
  d_large <- expand.grid(
    Person = paste0("P", sprintf("%02d", seq_len(n_person))),
    Rater = paste0("R", 1:2),
    Task = paste0("T", 1:2),
    stringsAsFactors = FALSE
  )
  d_large$Score <- sample(1:3, nrow(d_large), replace = TRUE)

  fit <- suppressWarnings(fit_mfrm(d_large, "Person",
    c("Rater", "Task"), "Score", method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_equal(nrow(fit$facets$person), n_person)
})

test_that("non-numeric score column (character) is coerced correctly", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Score <- as.character(d$Score)

  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 30))

  expect_s3_class(fit, "mfrm_fit")
  expect_true(is.integer(fit$prep$data$Score) || is.numeric(fit$prep$data$Score))
})

test_that("prepare_mfrm_data returns expected structure", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  prep <- mfrmr:::prepare_mfrm_data(d, "Person",
    c("Rater", "Task", "Criterion"), "Score")

  expect_true(is.list(prep))
  expect_true("data" %in% names(prep))
  expect_true("levels" %in% names(prep))
  expect_true("facet_names" %in% names(prep))
  expect_true("rating_min" %in% names(prep))
  expect_true("rating_max" %in% names(prep))
  expect_equal(prep$facet_names, c("Rater", "Task", "Criterion"))
  expect_true("Person" %in% names(prep$levels))
  expect_true("Rater" %in% names(prep$levels))
})

test_that("duplicate person/facets/score arguments raise error", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_error(
    mfrmr:::prepare_mfrm_data(d, "Person", c("Person", "Task"), "Score"),
    regex = "[Dd]uplic"
  )
})
