# test-numerical-validation.R
# Direct tests of internal mathematical functions for correctness and stability.
# All tests target the release code (mfrm_core.R) through the ::: accessor.

# ---- logsumexp ----

test_that("logsumexp computes log(sum(exp(x))) correctly", {
  lse <- mfrmr:::logsumexp
  expect_equal(lse(log(c(1, 2, 3))), log(6), tolerance = 1e-10)
  expect_equal(lse(c(0, 0)), log(2), tolerance = 1e-10)
  expect_equal(lse(0), 0, tolerance = 1e-10)
})

test_that("logsumexp is numerically stable for large values", {
  lse <- mfrmr:::logsumexp
  # Large positive values (would overflow naive exp)
  expect_equal(lse(c(1000, 1001)), 1001 + log(1 + exp(-1)), tolerance = 1e-10)
  # Large negative values (would underflow naive exp)
  expect_equal(lse(c(-1000, -999)), -999 + log(1 + exp(-1)), tolerance = 1e-10)
})

# ---- weighted_mean ----

test_that("weighted_mean handles standard cases", {
  wm <- mfrmr:::weighted_mean
  expect_equal(wm(c(1, 2, 3), c(1, 1, 1)), 2)
  expect_equal(wm(c(1, 3), c(1, 3)), 2.5)
})

test_that("weighted_mean handles edge cases safely", {
  wm <- mfrmr:::weighted_mean
  expect_true(is.na(wm(c(NA, NA), c(1, 1))))
  expect_true(is.na(wm(c(1, 2), c(0, 0))))
  expect_true(is.na(wm(numeric(0), numeric(0))))
  # Inf is not finite, so it's excluded; result is just the finite value
  expect_equal(wm(c(Inf, 1), c(1, 1)), 1)
})

# ---- get_weights ----

test_that("get_weights extracts Weight column when present", {
  gw <- mfrmr:::get_weights
  df <- data.frame(X = 1:3, Weight = c(1, 2, 0.5))
  expect_equal(gw(df), c(1, 2, 0.5))
})

test_that("get_weights returns ones when Weight column absent", {
  gw <- mfrmr:::get_weights
  df <- data.frame(X = 1:3, Y = 4:6)
  expect_equal(gw(df), c(1, 1, 1))
})

test_that("get_weights zeroes invalid weights", {
  gw <- mfrmr:::get_weights
  df <- data.frame(X = 1:3, Weight = c(1, -2, NA))
  result <- gw(df)
  expect_equal(result[1], 1)
  expect_equal(result[2], 0)
  expect_equal(result[3], 0)
})

# ---- gauss_hermite_normal ----

test_that("gauss_hermite_normal returns valid quadrature for n=1", {
  gh <- mfrmr:::gauss_hermite_normal
  q1 <- gh(1)
  expect_equal(q1$nodes, 0)
  expect_equal(q1$weights, 1)
})

test_that("gauss_hermite_normal weights sum to 1 for standard normal", {
  gh <- mfrmr:::gauss_hermite_normal
  for (n in c(3, 7, 15, 21)) {
    q <- gh(n)
    expect_equal(length(q$nodes), n)
    expect_equal(length(q$weights), n)
    expect_equal(sum(q$weights), 1, tolerance = 1e-10,
                 label = paste("n =", n))
  }
})

test_that("gauss_hermite_normal nodes are symmetric around zero", {
  gh <- mfrmr:::gauss_hermite_normal
  q <- gh(15)
  sorted <- sort(q$nodes)
  expect_equal(sorted, -rev(sorted), tolerance = 1e-10)
})

# ---- center_sum_zero ----

test_that("center_sum_zero produces zero-sum vector", {
  csz <- mfrmr:::center_sum_zero
  x <- c(1, 2, 3, 10)
  result <- csz(x)
  expect_equal(sum(result), 0, tolerance = 1e-15)
  expect_equal(result, x - mean(x), tolerance = 1e-15)
})

test_that("center_sum_zero handles empty input", {
  expect_equal(mfrmr:::center_sum_zero(numeric(0)), numeric(0))
})

# ---- expand_facet ----

test_that("expand_facet produces sum-to-zero constraint", {
  ef <- mfrmr:::expand_facet
  result <- ef(c(0.5, -0.3), 3)
  expect_equal(sum(result), 0, tolerance = 1e-15)
  expect_equal(length(result), 3)
})

test_that("expand_facet returns zero for single level", {
  expect_equal(mfrmr:::expand_facet(numeric(0), 1), 0)
})

# ---- zstd_from_mnsq ----

test_that("zstd_from_mnsq gives near-zero for perfect fit (MnSq=1)", {
  zstd <- mfrmr:::zstd_from_mnsq
  result <- zstd(1.0, 100)
  expect_true(abs(result) < 0.5)
})

test_that("zstd_from_mnsq gives positive for overfit (MnSq>1)", {
  zstd <- mfrmr:::zstd_from_mnsq
  result <- zstd(2.0, 100)
  expect_true(result > 0)
})

test_that("zstd_from_mnsq returns NA for NA or zero df", {
  zstd <- mfrmr:::zstd_from_mnsq
  expect_true(is.na(zstd(NA, 100)))
  expect_true(is.na(zstd(1.0, NA)))
})

# ---- category_prob_rsm ----

test_that("category_prob_rsm probabilities sum to 1 per observation", {
  cp <- mfrmr:::category_prob_rsm
  eta <- c(-2, -1, 0, 1, 2)
  step_cum <- c(0, -0.5, 0.5)  # 3 categories (0, 1, 2)
  probs <- cp(eta, step_cum)
  expect_equal(nrow(probs), 5)
  expect_equal(ncol(probs), 3)
  row_sums <- rowSums(probs)
  expect_equal(row_sums, rep(1, 5), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

test_that("category_prob_rsm responds to eta direction", {
  cp <- mfrmr:::category_prob_rsm
  step_cum <- c(0, 0, 0)  # 3 categories
  # High eta should favor higher categories
  probs_high <- cp(5, step_cum)
  probs_low <- cp(-5, step_cum)
  expect_true(probs_high[1, 3] > probs_high[1, 1])
  expect_true(probs_low[1, 1] > probs_low[1, 3])
})

# ---- loglik_rsm ----

test_that("loglik_rsm returns finite non-positive values", {
  ll <- mfrmr:::loglik_rsm
  eta <- c(0, 0.5, -0.5)
  score_k <- c(0L, 1L, 2L)
  step_cum <- c(0, -0.3, 0.3)
  result <- ll(eta, score_k, step_cum)
  expect_true(is.finite(result))
  expect_true(result <= 0)
})

# ---- category_prob_pcm ----

test_that("category_prob_pcm probabilities sum to 1", {
  cp <- mfrmr:::category_prob_pcm
  eta <- c(-1, 0, 1)
  # step_cum_mat: rows = criterion levels, cols = categories
  step_cum_mat <- matrix(c(0, -0.5, 0.5,
                            0, -0.3, 0.3), nrow = 2, byrow = TRUE)
  step_idx <- c(1L, 2L, 1L)
  probs <- cp(eta, step_cum_mat, step_idx)
  expect_equal(nrow(probs), 3)
  row_sums <- rowSums(probs)
  expect_equal(row_sums, rep(1, 3), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

# ---- JML vs MML directional consistency ----

test_that("JML and MML produce correlated facet estimates", {
  skip_on_cran()
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit_jml <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 30
  ))
  fit_mml <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 30, quad_points = 15
  ))

  jml_rater <- fit_jml$facets$others |>
    dplyr::filter(.data$Facet == "Rater") |>
    dplyr::arrange(.data$Level)
  mml_rater <- fit_mml$facets$others |>
    dplyr::filter(.data$Facet == "Rater") |>
    dplyr::arrange(.data$Level)

  cor_val <- cor(jml_rater$Estimate, mml_rater$Estimate)
  expect_true(cor_val > 0.9, label = paste("JML-MML rater correlation:", round(cor_val, 3)))
})

# ---- build_param_sizes / split_params round-trip ----

test_that("split_params correctly partitions parameter vector", {
  sp <- mfrmr:::split_params
  sizes <- list(theta = 3, Rater = 2, steps = 2)
  par <- c(0.1, 0.2, 0.3, 0.5, -0.5, 0.7, -0.7)
  result <- sp(par, sizes)
  expect_equal(result$theta, c(0.1, 0.2, 0.3))
  expect_equal(result$Rater, c(0.5, -0.5))
  expect_equal(result$steps, c(0.7, -0.7))
})

test_that("split_params handles zero-size correctly", {
  sp <- mfrmr:::split_params
  sizes <- list(theta = 0, Rater = 2)
  par <- c(0.5, -0.5)
  result <- sp(par, sizes)
  expect_equal(result$theta, numeric(0))
  expect_equal(result$Rater, c(0.5, -0.5))
})

# ---- sample_mfrm_data reproducibility ----

test_that("sample_mfrm_data is reproducible with same seed", {
  d1 <- mfrmr:::sample_mfrm_data(seed = 999)
  d2 <- mfrmr:::sample_mfrm_data(seed = 999)
  expect_identical(d1, d2)
  expect_true(is.data.frame(d1))
  expect_true(all(c("Person", "Rater", "Task", "Criterion", "Score") %in% names(d1)))
})
