test_that("plot APIs accept title/palette/label customization", {
  d <- mfrmr:::sample_mfrm_data(seed = 321)

  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "JML",
    model = "RSM",
    maxit = 20,
    quad_points = 7
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  expect_no_error(
    plot(
      fit,
      type = "wright",
      draw = FALSE,
      title = "Custom Wright",
      palette = c(facet_level = "#1f78b4", step_threshold = "#d95f02"),
      label_angle = 45
    )
  )

  expect_no_error(
    mfrmr::plot_unexpected(
      fit,
      diagnostics = diag,
      plot_type = "severity",
      draw = FALSE,
      main = "Custom Unexpected",
      palette = c(higher = "#d95f02", lower = "#1b9e77", bar = "#2b8cbe"),
      label_angle = 45
    )
  )

  expect_no_error(
    mfrmr::plot_interrater_agreement(
      fit,
      diagnostics = diag,
      rater_facet = "Rater",
      plot_type = "exact",
      draw = FALSE,
      main = "Custom Inter-rater",
      palette = c(ok = "#2b8cbe", flag = "#cb181d", expected = "#08519c"),
      label_angle = 45
    )
  )

  bias3 <- mfrmr::estimate_bias(
    fit,
    diag,
    interaction_facets = c("Rater", "Task", "Criterion"),
    max_iter = 2
  )
  t13_3 <- mfrmr::bias_interaction_report(bias3, top_n = 20)

  expect_no_error(
    mfrmr::plot_bias_interaction(
      t13_3,
      plot = "facet_profile",
      draw = FALSE,
      main = "Custom Higher-Order Bias Profile",
      palette = c(normal = "#2b8cbe", flag = "#cb181d", profile = "#756bb1"),
      label_angle = 45
    )
  )
})
