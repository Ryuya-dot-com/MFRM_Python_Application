# --------------------------------------------------------------------------
# test-exception-regression.R
# Regression guard: exception handling must remain stable across versions.
# --------------------------------------------------------------------------

# === 6.1 Unused category ===================================================

test_that("unused intermediate category does not crash estimation", {
  set.seed(42)
  d <- data.frame(
    Person = rep(paste0("P", 1:10), each = 3),
    Rater  = rep(paste0("R", 1:3), 10),
    Score  = sample(c(1, 3, 5), 30, replace = TRUE),
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(d, "Person", "Rater", "Score",
    method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
})

# === 6.2 Constant rater (all same score) ===================================

test_that("constant rater does not crash estimation", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  # Make R3 always score 3
  d$Score[d$Rater == "R3"] <- 3
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
})

# === 6.3 Single-level facet ================================================

test_that("single-level facet is handled gracefully", {
  d <- data.frame(
    Person = rep(paste0("P", 1:8), each = 2),
    Rater  = rep("OnlyRater", 16),
    Task   = rep(c("T1", "T2"), 8),
    Score  = c(1, 2, 2, 3, 1, 3, 2, 1, 3, 2, 1, 2, 2, 3, 1, 2),
    stringsAsFactors = FALSE
  )
  fit <- suppressWarnings(fit_mfrm(d, "Person", c("Rater", "Task"), "Score",
    method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |> dplyr::pull(Estimate)
  expect_equal(unname(rater_est), 0, tolerance = 1e-8)
})

# === 6.4 Scattered NA handling =============================================

test_that("scattered NAs are dropped and estimation succeeds", {
  d <- mfrmr:::sample_mfrm_data(seed = 300)
  set.seed(300)
  na_idx <- sample(nrow(d), size = round(0.05 * nrow(d)))
  d$Score[na_idx] <- NA
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
  expect_lt(nrow(fit$prep$data), nrow(d))
})

# === 6.5 Error message regression tests ====================================

test_that("fit_mfrm error messages are stable", {
  expect_error(fit_mfrm(42, "P", "R", "S"), "data.frame")
  expect_error(fit_mfrm(data.frame(), "P", "R", "S"), "zero rows")
  expect_error(fit_mfrm(data.frame(P = 1), 123, "R", "S"), "character")
  expect_error(fit_mfrm(data.frame(P = 1), "P", character(0), "S"), "facet")
  expect_error(fit_mfrm(data.frame(P = 1, R = 1, S = 1), "P", "R", "S",
    maxit = -1), "positive")
})

# === 6.6 Non-convergence warning ===========================================

test_that("non-convergence warning message is stable", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_warning(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", maxit = 1),
    "converge"
  )
})

# === 6.7 diagnose_mfrm input guards =======================================

test_that("diagnose_mfrm rejects non-mfrm_fit input", {
  expect_error(diagnose_mfrm(list()), "mfrm_fit")
  expect_error(diagnose_mfrm(42), "mfrm_fit")
})

# === 6.8 estimate_bias input guards ========================================

test_that("estimate_bias rejects invalid inputs", {
  expect_error(estimate_bias(list()), "mfrm_fit")
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score", method = "JML", maxit = 20))
  dx <- diagnose_mfrm(fit, residual_pca = "none")
  # Non-existent facet returns an empty result (not an error)
  res <- estimate_bias(fit, dx, facet_a = "NonExistent", facet_b = "Rater")
  expect_true(is.list(res))
  expect_equal(length(names(res)), 0)
})

# === 6.9 Extremely small dataset ==========================================

test_that("extremely small dataset estimation attempt", {
  d <- data.frame(
    Person = c("A", "A", "B", "B"),
    Rater  = c("X", "Y", "X", "Y"),
    Score  = c(1, 2, 2, 3),
    stringsAsFactors = FALSE
  )
  # Should either succeed or give an informative error/warning, not crash
  result <- tryCatch(
    suppressWarnings(fit_mfrm(d, "Person", "Rater", "Score",
      method = "JML", maxit = 20)),
    error = function(e) e
  )
  expect_true(inherits(result, "mfrm_fit") || inherits(result, "error"))
})

# === 6.10 Zero-weight rows excluded =======================================

test_that("zero-weight rows are excluded from estimation", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  d$Weight <- 1
  d$Weight[1:10] <- 0
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    weight = "Weight", method = "JML", maxit = 30))
  expect_s3_class(fit, "mfrm_fit")
  expect_equal(nrow(fit$prep$data), nrow(d) - 10)
})

# === 6.11 Invalid anchor warning ===========================================

test_that("anchor with non-existent level produces warning", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchors <- data.frame(
    Facet = "Rater", Level = "NonExistent", Anchor = 0,
    stringsAsFactors = FALSE
  )
  expect_warning(
    fit_mfrm(d, "Person",
      c("Rater", "Task", "Criterion"), "Score",
      anchors = anchors, method = "JML", maxit = 30),
    "anchor|Anchor"
  )
})

# === 6.12 Anchor policy "error" mode =======================================

test_that("anchor_policy error stops on invalid anchors", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  anchors <- data.frame(
    Facet = "Rater", Level = "NonExistent", Anchor = 0,
    stringsAsFactors = FALSE
  )
  expect_error(
    fit_mfrm(d, "Person", c("Rater", "Task", "Criterion"), "Score",
      anchors = anchors, anchor_policy = "error", method = "JML", maxit = 30)
  )
})
