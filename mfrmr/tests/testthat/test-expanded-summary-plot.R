test_that("expanded bundle classes dispatch through plot.mfrm_bundle", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      2 * as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- mfrmr::fit_mfrm(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "both")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)

  targets <- list(
    list(obj = mfrmr::measurable_summary_table(fit, diag), type = "facet_coverage"),
    list(obj = mfrmr::unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20), type = "comparison"),
    list(obj = mfrmr::facets_output_file_bundle(fit, diagnostics = diag, include = c("graph", "score")), type = "graph_expected"),
    list(obj = mfrmr::analyze_residual_pca(diag, mode = "both"), type = "overall_scree"),
    list(obj = mfrmr::specifications_report(fit), type = "facet_elements"),
    list(obj = mfrmr::data_quality_report(fit, data = toy, person = "Person", facets = c("Rater", "Criterion"), score = "Score"), type = "row_audit"),
    list(obj = mfrmr::estimation_iteration_report(fit, max_iter = 5), type = "residual"),
    list(obj = mfrmr::subset_connectivity_report(fit), type = "subset_observations"),
    list(obj = mfrmr::facet_statistics_report(fit), type = "means")
  )

  for (tr in targets) {
    plt <- plot(tr$obj, type = tr$type, draw = FALSE)
    expect_s3_class(plt, "mfrm_plot_data")
  }
})

test_that("expanded bundle classes provide class-aware summary output", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      2 * as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- mfrmr::fit_mfrm(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "both")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)

  targets <- list(
    list(obj = mfrmr::measurable_summary_table(fit, diag), kind = "mfrm_measurable"),
    list(obj = mfrmr::unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20), kind = "mfrm_unexpected_after_bias"),
    list(obj = mfrmr::facets_output_file_bundle(fit, diagnostics = diag, include = c("graph", "score")), kind = "mfrm_output_bundle"),
    list(obj = mfrmr::analyze_residual_pca(diag, mode = "both"), kind = "mfrm_residual_pca"),
    list(obj = mfrmr::specifications_report(fit), kind = "mfrm_specifications"),
    list(obj = mfrmr::data_quality_report(fit, data = toy, person = "Person", facets = c("Rater", "Criterion"), score = "Score"), kind = "mfrm_data_quality"),
    list(obj = mfrmr::estimation_iteration_report(fit, max_iter = 5), kind = "mfrm_iteration_report"),
    list(obj = mfrmr::subset_connectivity_report(fit), kind = "mfrm_subset_connectivity"),
    list(obj = mfrmr::facet_statistics_report(fit), kind = "mfrm_facet_statistics")
  )

  for (tr in targets) {
    out <- summary(tr$obj)
    expect_s3_class(out, "summary.mfrm_bundle")
    expect_identical(out$summary_kind, tr$kind)
    expect_true(is.data.frame(out$overview))
    expect_true(nrow(out$overview) == 1)
  }
})

test_that("core bundle classes now return class-aware summary_kind", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      2 * as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- mfrmr::fit_mfrm(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)

  targets <- list(
    list(obj = mfrmr::unexpected_response_table(fit, diagnostics = diag, top_n = 10), kind = "mfrm_unexpected"),
    list(obj = mfrmr::fair_average_table(fit, diagnostics = diag), kind = "mfrm_fair_average"),
    list(obj = mfrmr::displacement_table(fit, diagnostics = diag), kind = "mfrm_displacement"),
    list(obj = mfrmr::interrater_agreement_table(fit, diagnostics = diag), kind = "mfrm_interrater"),
    list(obj = mfrmr::facets_chisq_table(fit, diagnostics = diag), kind = "mfrm_facets_chisq"),
    list(obj = mfrmr::bias_interaction_report(bias), kind = "mfrm_bias_interaction"),
    list(obj = mfrmr::rating_scale_table(fit, diagnostics = diag), kind = "mfrm_rating_scale"),
    list(obj = mfrmr::category_structure_report(fit, diagnostics = diag), kind = "mfrm_category_structure"),
    list(obj = mfrmr::category_curves_report(fit, theta_points = 101), kind = "mfrm_category_curves")
  )

  for (tr in targets) {
    out <- summary(tr$obj)
    expect_s3_class(out, "summary.mfrm_bundle")
    expect_identical(out$summary_kind, tr$kind)
    expect_true(is.data.frame(out$overview))
    expect_true(nrow(out$overview) == 1)
  }
})

test_that("print.summary.mfrm_bundle uses class-aware titles", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      2 * as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  fit <- mfrmr::fit_mfrm(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    method = "JML",
    maxit = 25
  )
  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")

  t4 <- mfrmr::unexpected_response_table(fit, diagnostics = diag, top_n = 10)
  t4_sum <- summary(t4)
  out_t4 <- paste(capture.output(print(t4_sum)), collapse = "\n")
  expect_match(out_t4, "mfrmr Unexpected Response Summary", fixed = TRUE)
  expect_match(out_t4, "Threshold summary", fixed = TRUE)

  rs <- mfrmr::rating_scale_table(fit, diagnostics = diag)
  rs_sum <- summary(rs)
  out_rs <- paste(capture.output(print(rs_sum)), collapse = "\n")
  expect_match(out_rs, "mfrmr Rating Scale Summary", fixed = TRUE)
  expect_match(out_rs, "Category/threshold summary", fixed = TRUE)
})

test_that("data description and anchor audit support summary and plot", {
  toy <- expand.grid(
    Person = paste0("P", 1:4),
    Rater = paste0("R", 1:2),
    Criterion = c("Content", "Organization", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      2 * as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  ds <- mfrmr::describe_mfrm_data(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score"
  )
  ds_sum <- summary(ds)
  expect_s3_class(ds_sum, "summary.mfrm_data_description")
  ds_plot <- plot(ds, type = "score_distribution", draw = FALSE)
  expect_s3_class(ds_plot, "mfrm_plot_data")

  aud <- mfrmr::audit_mfrm_anchors(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score"
  )
  aud_sum <- summary(aud)
  expect_s3_class(aud_sum, "summary.mfrm_anchor_audit")
  aud_plot <- plot(aud, type = "issue_counts", draw = FALSE)
  expect_s3_class(aud_plot, "mfrm_plot_data")
})
