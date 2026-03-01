tolerance <- 1e-8

safe_num <- function(x) suppressWarnings(as.numeric(x))

test_that("FACETS-style metric contracts hold for key summary tables", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Task", max_iter = 2)

  t4 <- mfrmr::unexpected_response_table(fit, diagnostics = diag, top_n = 20)
  s4 <- t4$summary[1, , drop = FALSE]
  total_obs <- safe_num(s4$TotalObservations)
  unexpected_n <- safe_num(s4$UnexpectedN)
  unexpected_pct <- safe_num(s4$UnexpectedPercent)
  if (is.finite(total_obs) && total_obs > 0) {
    expect_equal(unexpected_pct, 100 * unexpected_n / total_obs, tolerance = tolerance)
  }

  t10 <- mfrmr::unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20)
  s10 <- t10$summary[1, , drop = FALSE]
  baseline <- safe_num(s10$BaselineUnexpectedN)
  after <- safe_num(s10$AfterBiasUnexpectedN)
  reduced <- safe_num(s10$ReducedBy)
  reduced_pct <- safe_num(s10$ReducedPercent)
  if (all(is.finite(c(baseline, after, reduced)))) {
    expect_equal(reduced, baseline - after, tolerance = tolerance)
  }
  if (is.finite(baseline) && baseline > 0 && is.finite(reduced_pct)) {
    expect_equal(reduced_pct, 100 * reduced / baseline, tolerance = 1e-6)
  }

  t11 <- mfrmr::bias_count_table(bias, branch = "facets")
  s11 <- t11$summary[1, , drop = FALSE]
  cells <- safe_num(s11$Cells)
  low <- safe_num(s11$LowCountCells)
  low_pct <- safe_num(s11$LowCountPercent)
  if (is.finite(cells) && cells > 0 && is.finite(low_pct)) {
    expect_equal(low_pct, 100 * low / cells, tolerance = 1e-6)
  }

  t3 <- mfrmr::estimation_iteration_report(fit, max_iter = 5)
  s3 <- t3$summary[1, , drop = FALSE]
  expect_true(is.logical(s3$FinalConverged) || s3$FinalConverged %in% c(TRUE, FALSE))
  expect_true(safe_num(s3$FinalIterations) >= 1)
  expect_true(safe_num(s3$ReplayRows) >= 1)
})

test_that("FACETS-style range contracts hold for agreement, fit, displacement, and rating scale", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)
  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  t7agree <- mfrmr::interrater_agreement_table(fit, diagnostics = diag)
  if (nrow(t7agree$summary) > 0) {
    s <- t7agree$summary[1, , drop = FALSE]
    exact <- safe_num(s$ExactAgreement)
    expected_exact <- safe_num(s$ExpectedExactAgreement)
    adjacent <- safe_num(s$AdjacentAgreement)
    expect_true(exact >= -tolerance && exact <= 1 + tolerance)
    expect_true(expected_exact >= -tolerance && expected_exact <= 1 + tolerance)
    expect_true(adjacent >= -tolerance && adjacent <= 1 + tolerance)
  }

  t7chisq <- mfrmr::facets_chisq_table(fit, diagnostics = diag)
  tbl7 <- t7chisq$table
  if (nrow(tbl7) > 0) {
    fp <- safe_num(tbl7$FixedProb)
    rp <- safe_num(tbl7$RandomProb)
    expect_true(all(fp[is.finite(fp)] >= -tolerance & fp[is.finite(fp)] <= 1 + tolerance))
    expect_true(all(rp[is.finite(rp)] >= -tolerance & rp[is.finite(rp)] <= 1 + tolerance))
  }

  disp <- mfrmr::displacement_table(fit, diagnostics = diag)
  sdisp <- disp$summary[1, , drop = FALSE]
  levels_n <- safe_num(sdisp$Levels)
  anchored <- safe_num(sdisp$AnchoredLevels)
  flagged <- safe_num(sdisp$FlaggedLevels)
  flagged_anch <- safe_num(sdisp$FlaggedAnchoredLevels)
  expect_true(anchored <= levels_n + tolerance)
  expect_true(flagged <= levels_n + tolerance)
  expect_true(flagged_anch <= anchored + tolerance)

  t81 <- mfrmr::rating_scale_table(fit, diagnostics = diag)
  s81 <- t81$summary[1, , drop = FALSE]
  cats <- safe_num(s81$Categories)
  used <- safe_num(s81$UsedCategories)
  expect_true(used <= cats + tolerance)

  tt <- t81$threshold_table
  if (nrow(tt) > 1 && "GapFromPrev" %in% names(tt)) {
    gaps <- safe_num(tt$GapFromPrev)
    expect_equal(isTRUE(s81$ThresholdMonotonic), !any(gaps[is.finite(gaps)] < -tolerance))
  }
})
