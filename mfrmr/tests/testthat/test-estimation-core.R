# --------------------------------------------------------------------------
# test-estimation-core.R
# Tests the likelihood and estimation core of the mfrmr package.
# --------------------------------------------------------------------------

# ---- 2.1  RSM category probabilities sum to 1 and are non-negative --------

test_that("category_prob_rsm rows sum to 1 and are non-negative", {
  step_cum <- c(0, -0.5, 0.3, 0.7)
  etas <- c(-3, -1, 0, 1, 3)
  probs <- mfrmr:::category_prob_rsm(etas, step_cum)
  expect_equal(rowSums(probs), rep(1, length(etas)), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

# ---- 2.2  PCM category probabilities sum to 1 and are non-negative --------

test_that("category_prob_pcm rows sum to 1 and are non-negative", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  etas <- c(-1, 0, 1, 0.5)
  crit_idx <- c(1L, 2L, 1L, 2L)
  probs <- mfrmr:::category_prob_pcm(etas, step_cum_mat, crit_idx)
  expect_equal(rowSums(probs), rep(1, length(etas)), tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

# ---- 2.3a  Monotonicity -- higher eta -> higher expected score (RSM) ------

test_that("higher eta produces monotonically higher expected score under RSM", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  etas <- seq(-3, 3, by = 0.5)
  probs <- mfrmr:::category_prob_rsm(etas, step_cum)
  k_vals <- 0:(ncol(probs) - 1)
  expected <- as.vector(probs %*% k_vals)
  diffs <- diff(expected)
  expect_true(all(diffs >= -1e-10),
              info = "Expected score must be non-decreasing in eta")
})

# ---- 2.3b  Monotonicity -- higher eta -> higher expected score (PCM) ------

test_that("higher eta produces monotonically higher expected score under PCM", {
  step_cum_mat <- matrix(c(0, 0.3, 0.8, 1.5,
                           0, 0.1, 0.6, 1.2), nrow = 2, byrow = TRUE)
  etas <- seq(-3, 3, by = 0.5)
  # All observations belong to criterion 1
  crit_idx <- rep(1L, length(etas))
  probs <- mfrmr:::category_prob_pcm(etas, step_cum_mat, crit_idx)
  k_vals <- 0:(ncol(probs) - 1)
  expected <- as.vector(probs %*% k_vals)
  diffs <- diff(expected)
  expect_true(all(diffs >= -1e-10),
              info = "Expected score must be non-decreasing in eta (PCM)")
})

# ---- 2.4  Numerical stability for extreme eta values ---------------------

test_that("extreme eta values do not produce NaN or Inf", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  extreme_etas <- c(-50, -10, 0, 10, 50)
  probs <- mfrmr:::category_prob_rsm(extreme_etas, step_cum)
  expect_true(all(is.finite(probs)))
  expect_equal(rowSums(probs), rep(1, length(extreme_etas)), tolerance = 1e-10)
})

# ---- 2.5a  RSM log-likelihood is finite and non-positive -----------------

test_that("RSM log-likelihood is finite and non-positive", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  score_k <- c(0L, 1L, 2L, 3L, 1L)
  eta <- c(-1, 0, 0.5, 1, -0.5)
  ll <- mfrmr:::loglik_rsm(eta, score_k, step_cum)
  expect_true(is.finite(ll))
  expect_lte(ll, 0)
})

# ---- 2.5b  PCM log-likelihood is finite and non-positive -----------------

test_that("PCM log-likelihood is finite and non-positive", {
  step_cum_mat <- matrix(c(0, -0.3, 0.2, 0.5,
                           0, -0.1, 0.4, 0.8), nrow = 2, byrow = TRUE)
  score_k <- c(0L, 1L, 2L, 3L)
  eta <- c(-1, 0, 0.5, 1)
  crit_idx <- c(1L, 2L, 1L, 2L)
  ll <- mfrmr:::loglik_pcm(eta, score_k, step_cum_mat, crit_idx)
  expect_true(is.finite(ll))
  expect_lte(ll, 0)
})

# ---- 2.6  JML vs MML consistency ----------------------------------------

test_that("JML and MML facet estimates are highly correlated", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)

  fit_jml <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 50, quad_points = 7
  ))
  fit_mml <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 50, quad_points = 7
  ))

  for (facet in c("Rater", "Task", "Criterion")) {
    est_jml <- fit_jml$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    est_mml <- fit_mml$facets$others |>
      dplyr::filter(Facet == facet) |>
      dplyr::arrange(Level) |>
      dplyr::pull(Estimate)
    r <- cor(est_jml, est_mml)
    expect_gt(r, 0.9,
              label = paste("JML-MML correlation for", facet))
  }
  expect_true(is.finite(fit_jml$summary$LogLik))
  expect_true(is.finite(fit_mml$summary$LogLik))
})

# ---- 2.7a  Convergence -- optim converges with adequate maxit ------------

test_that("optim converges with adequate maxit", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 100, quad_points = 7
  ))
  expect_true(fit$summary$Converged)
  expect_true(is.finite(fit$summary$LogLik))
  expect_lt(fit$summary$LogLik, 0)
})

# ---- 2.7b  Convergence -- maxit=1 produces warning ----------------------

test_that("maxit=1 produces a convergence warning", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_warning(
    fit_mfrm(
      d, "Person", c("Rater", "Task", "Criterion"), "Score",
      method = "JML", model = "RSM", maxit = 1, quad_points = 7
    ),
    "converge|converg",
    ignore.case = TRUE
  )
})

# ---- 2.7c  Convergence -- more iterations give better LogLik -------------

test_that("more iterations produce equal or better log-likelihood", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_short <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 5, quad_points = 7
  ))
  fit_long <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 100, quad_points = 7
  ))
  # More iterations should give equal or better (less negative) log-likelihood
  expect_gte(fit_long$summary$LogLik, fit_short$summary$LogLik - 0.1)
})

# ---- 2.8  Step parameters are finite ------------------------------------

test_that("RSM step estimates are finite", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 50, quad_points = 7
  ))
  expect_true(all(is.finite(fit$steps$Estimate)))
  # After centering, steps should sum to approximately 0
  expect_equal(sum(fit$steps$Estimate), 0, tolerance = 1e-6)
})

# ---- 2.9a  EAP person estimates are all finite (MML) --------------------

test_that("all EAP person estimates are finite", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 50, quad_points = 7
  ))
  expect_true(all(is.finite(fit$facets$person$Estimate)))
})

# ---- 2.9b  EAP SDs are all positive (MML) -------------------------------

test_that("all EAP person SDs are positive", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 50, quad_points = 7
  ))
  expect_true("SD" %in% names(fit$facets$person))
  expect_true(all(fit$facets$person$SD > 0))
})

# ---- 2.10a  ZSTD: mnsq = 1.0 -> ZSTD near 0 ----------------------------

test_that("ZSTD is near zero when mnsq equals 1", {
  z <- mfrmr:::zstd_from_mnsq(1.0, df = 30)
  expect_equal(z, 0, tolerance = 0.3)
})

# ---- 2.10b  ZSTD: mnsq > 1 -> ZSTD > 0 ---------------------------------

test_that("ZSTD is positive when mnsq exceeds 1", {
  z <- mfrmr:::zstd_from_mnsq(1.5, df = 30)
  expect_gt(z, 0)
})

# ---- 2.10c  ZSTD: mnsq < 1 -> ZSTD < 0 ---------------------------------

test_that("ZSTD is negative when mnsq is below 1", {
  z <- mfrmr:::zstd_from_mnsq(0.7, df = 30)
  expect_lt(z, 0)
})

# ---- 2.11a  Gauss-Hermite: weights sum to 1 -----------------------------

test_that("Gauss-Hermite weights sum to 1", {
  for (n in c(5, 11, 21)) {
    gh <- mfrmr:::gauss_hermite_normal(n)
    expect_equal(sum(gh$weights), 1, tolerance = 1e-10,
                 label = paste("GH weights sum for n =", n))
  }
})

# ---- 2.11b  Gauss-Hermite: nodes are symmetric around 0 -----------------

test_that("Gauss-Hermite nodes are symmetric around 0", {
  for (n in c(5, 11, 21)) {
    gh <- mfrmr:::gauss_hermite_normal(n)
    nodes_sorted <- sort(gh$nodes)
    nodes_neg_sorted <- sort(-gh$nodes)
    expect_equal(nodes_sorted, nodes_neg_sorted, tolerance = 1e-10,
                 label = paste("GH node symmetry for n =", n))
  }
})

# ---- 2.11c  Gauss-Hermite: E[X] approx 0, E[X^2] approx 1 for N(0,1) ---

test_that("Gauss-Hermite recovers N(0,1) moments", {
  gh <- mfrmr:::gauss_hermite_normal(21)
  ex <- sum(gh$weights * gh$nodes)
  ex2 <- sum(gh$weights * gh$nodes^2)
  expect_equal(ex, 0, tolerance = 0.01,
               label = "E[X] should be ~0 under standard normal")
  expect_equal(ex2, 1, tolerance = 0.01,
               label = "E[X^2] should be ~1 under standard normal")
})

# ---- 2.12  LogLik non-deterioration across optimization ------------------

test_that("log-likelihood does not deteriorate with more iterations (MML)", {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit_short <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 5, quad_points = 7
  ))
  fit_long <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "MML", model = "RSM", maxit = 100, quad_points = 7
  ))
  expect_gte(fit_long$summary$LogLik, fit_short$summary$LogLik - 0.1)
})

# ---- Additional: logsumexp stability ------------------------------------

test_that("logsumexp handles large values without overflow", {
  x <- c(1000, 1001, 999)
  result <- mfrmr:::logsumexp(x)
  expect_true(is.finite(result))
  # Should be close to log(exp(1000) + exp(1001) + exp(999)) ~= 1001.41
  expect_equal(result, log(exp(0) + exp(1) + exp(-1)) + 1000, tolerance = 1e-10)
})

test_that("logsumexp handles very negative values", {
  x <- c(-1000, -999, -1001)
  result <- mfrmr:::logsumexp(x)
  expect_true(is.finite(result))
})

# ---- Additional: RSM probability boundary behavior -----------------------

test_that("RSM at extreme positive eta concentrates on highest category", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  probs <- mfrmr:::category_prob_rsm(50, step_cum)
  # At very high eta, highest category should dominate

  expect_gt(probs[1, ncol(probs)], 0.99)
})

test_that("RSM at extreme negative eta concentrates on lowest category", {
  step_cum <- c(0, 0.5, 1.0, 1.5)
  probs <- mfrmr:::category_prob_rsm(-50, step_cum)
  # At very low eta, lowest category should dominate
  expect_gt(probs[1, 1], 0.99)
})
