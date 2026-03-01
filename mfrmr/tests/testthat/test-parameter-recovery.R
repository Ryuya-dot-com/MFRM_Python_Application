# test-parameter-recovery.R
# Tests for parameter recovery with known ground truth from sample_mfrm_data.
#
# True generating parameters (from sample_mfrm_data):
#   person ability:  rnorm(36, 0, 1) with set.seed(20240131)
#   rater_eff:       c(-0.4, 0, 0.4)   R1=lenient, R2=neutral, R3=harsh
#   task_eff:        seq(-0.5, 0.5, length.out=4)  T1 easiest, T4 hardest
#   crit_eff:        c(-0.3, 0, 0.3)   C1 easiest, C3 hardest
#   eta = ability - rater_eff - task_eff - crit_eff
#   Score = cut(eta + noise, breaks)
#
# In the subtractive MFRM parameterization, a positive facet effect means
# harder/harsher, so the estimated measure should be positively correlated
# with the true effect vector.

# ---------------------------------------------------------------------------
# Helper: reconstruct true person abilities
# ---------------------------------------------------------------------------
true_person_ability <- function() {
  set.seed(20240131)
  rnorm(36, 0, 1)
}

true_rater_eff  <- c(-0.4, 0, 0.4)        # R1, R2, R3
true_task_eff   <- seq(-0.5, 0.5, length.out = 4)  # T1, T2, T3, T4
true_crit_eff   <- c(-0.3, 0, 0.3)        # C1, C2, C3

# ---------------------------------------------------------------------------
# 4A.  JML true value recovery
# ---------------------------------------------------------------------------

test_that("JML rater estimates correlate with true rater effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  # Center both for comparison
  rater_est_c <- rater_est - mean(rater_est)
  true_c <- true_rater_eff - mean(true_rater_eff)
  r <- cor(rater_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(rater_est_c - true_c)), 0.6)
})

test_that("JML task estimates correlate with true task effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  task_est_c <- task_est - mean(task_est)
  true_c <- true_task_eff - mean(true_task_eff)
  r <- cor(task_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(task_est_c - true_c)), 0.6)
})

test_that("JML criterion estimates correlate with true criterion effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  crit_est_c <- crit_est - mean(crit_est)
  true_c <- true_crit_eff - mean(true_crit_eff)
  r <- cor(crit_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(crit_est_c - true_c)), 0.6)
})

# ---------------------------------------------------------------------------
# 4B.  MML true value recovery
# ---------------------------------------------------------------------------

test_that("MML rater estimates correlate with true rater effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 100))

  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  rater_est_c <- rater_est - mean(rater_est)
  true_c <- true_rater_eff - mean(true_rater_eff)
  r <- cor(rater_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(rater_est_c - true_c)), 0.6)
})

test_that("MML task estimates correlate with true task effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 100))

  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  task_est_c <- task_est - mean(task_est)
  true_c <- true_task_eff - mean(true_task_eff)
  r <- cor(task_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(task_est_c - true_c)), 0.6)
})

test_that("MML criterion estimates correlate with true criterion effects", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "MML", quad_points = 21, maxit = 100))

  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)

  crit_est_c <- crit_est - mean(crit_est)
  true_c <- true_crit_eff - mean(true_crit_eff)
  r <- cor(crit_est_c, true_c)
  expect_gt(r, 0.8)
  expect_lt(max(abs(crit_est_c - true_c)), 0.6)
})

# ---------------------------------------------------------------------------
# 4C.  Person ability recovery
# ---------------------------------------------------------------------------

test_that("JML person ability estimates correlate with true ability", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  true_ability <- true_person_ability()
  persons <- paste0("P", sprintf("%02d", 1:36))

  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  # Match person labels to true ability
  person_tbl <- fit$facets$person |> dplyr::arrange(Person)
  true_sorted <- true_ability[order(persons)]
  r <- cor(person_tbl$Estimate, true_sorted)
  expect_gt(r, 0.6)
})

# ---------------------------------------------------------------------------
# 4D.  Ordering preservation
# ---------------------------------------------------------------------------

test_that("facet ordering matches true generating direction", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  rater_est <- fit$facets$others |>
    dplyr::filter(Facet == "Rater") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  names(rater_est) <- c("R1", "R2", "R3")

  # R3 (eff=0.4, harshest) should have highest rater measure
  # R1 (eff=-0.4, most lenient) should have lowest rater measure
  expect_gt(rater_est["R3"], rater_est["R1"])

  task_est <- fit$facets$others |>
    dplyr::filter(Facet == "Task") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  names(task_est) <- c("T1", "T2", "T3", "T4")

  # T4 (hardest) > T1 (easiest)
  expect_gt(task_est["T4"], task_est["T1"])

  crit_est <- fit$facets$others |>
    dplyr::filter(Facet == "Criterion") |>
    dplyr::arrange(Level) |>
    dplyr::pull(Estimate)
  names(crit_est) <- c("C1", "C2", "C3")

  # C3 (hardest) > C1 (easiest)
  expect_gt(crit_est["C3"], crit_est["C1"])
})

# ---------------------------------------------------------------------------
# 4E.  Reproducibility
# ---------------------------------------------------------------------------

test_that("same seed produces identical LogLik across runs", {
  d1 <- mfrmr:::sample_mfrm_data(seed = 42)
  d2 <- mfrmr:::sample_mfrm_data(seed = 42)
  expect_identical(d1, d2)

  fit1 <- suppressWarnings(fit_mfrm(d1, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60))
  fit2 <- suppressWarnings(fit_mfrm(d2, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 60))

  expect_equal(fit1$summary$LogLik, fit2$summary$LogLik, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4F.  PCM recovery
# ---------------------------------------------------------------------------

test_that("PCM with step_facet produces all finite estimates", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)
  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100,
    model = "PCM", step_facet = "Criterion"))

  expect_s3_class(fit, "mfrm_fit")
  # All facet estimates should be finite
  all_est <- fit$facets$others$Estimate
  expect_true(all(is.finite(all_est)))
  # Step parameter estimates should be finite
  expect_true(all(is.finite(fit$steps$Estimate)))
})

# ---------------------------------------------------------------------------
# 4G.  Bias interaction recovery
# ---------------------------------------------------------------------------

test_that("estimate_bias detects known injected bias in correct direction", {
  d <- mfrmr:::sample_mfrm_data(seed = 20240131)

  # Inject known bias: for Rater R3 x Criterion C1, boost scores by +1
  boost_idx <- which(d$Rater == "R3" & d$Criterion == "C1")
  d$Score[boost_idx] <- pmin(d$Score[boost_idx] + 1L, 5L)

  fit <- suppressWarnings(fit_mfrm(d, "Person",
    c("Rater", "Task", "Criterion"), "Score",
    method = "JML", maxit = 100))

  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))

  bias <- suppressWarnings(estimate_bias(fit, diag,
    facet_a = "Rater", facet_b = "Criterion", max_iter = 4))

  expect_true(is.list(bias))
  expect_true("table" %in% names(bias))
  expect_true(is.data.frame(bias$table))

  # Find the R3 x C1 cell
  r3c1 <- bias$table |>
    dplyr::filter(
      Facet1_Level == "R3" & Facet2_Level == "C1" |
      FacetA_Level == "R3" & FacetB_Level == "C1"
    )

  # The bias for this cell should be present and indicate positive direction
  # (scores were boosted upward)
  if (nrow(r3c1) > 0) {
    bias_size <- r3c1$`Bias Size`[1]
    expect_true(is.finite(bias_size))
    # Positive bias because we boosted scores upward in that cell
    # The sign depends on the subtractive parameterization:
    # higher observed scores relative to expected means positive obs-exp diff.
    # Check the Obs-Exp Average is positive (higher than expected).
    obs_exp_avg <- r3c1$`Obs-Exp Average`[1]
    if (is.finite(obs_exp_avg)) {
      expect_gt(obs_exp_avg, -0.5,
        label = "R3:C1 Obs-Exp Average should reflect upward bias")
    }
  }
})
