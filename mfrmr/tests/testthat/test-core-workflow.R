test_that("core fit/diagnostics workflow runs", {
  d <- mfrmr:::sample_mfrm_data(seed = 123)

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

  expect_s3_class(fit, "mfrm_fit")
  expect_true("summary" %in% names(fit))

  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "both", pca_max_factors = 4)
  expect_true("overall_fit" %in% names(diag))
  expect_true("residual_pca_overall" %in% names(diag))
  expect_true("residual_pca_by_facet" %in% names(diag))

  pca <- mfrmr::analyze_residual_pca(diag, mode = "both", pca_max_factors = 4)
  expect_true(is.data.frame(pca$overall_table))
  expect_true(is.data.frame(pca$by_facet_table))
  expect_gt(nrow(pca$overall_table), 0)

  p1 <- mfrmr::plot_residual_pca(pca, mode = "overall", plot_type = "scree")
  p2 <- mfrmr::plot_residual_pca(pca, mode = "facet", facet = "Rater", plot_type = "loadings", top_n = 5)
  expect_s3_class(p1, "plotly")
  expect_s3_class(p2, "plotly")

  vis <- mfrmr::build_visual_summaries(fit, diag)
  expect_true("figure10" %in% names(vis$warning_map))
  expect_true("figure11" %in% names(vis$warning_map))
  expect_true("figure10" %in% names(vis$summary_map))
  expect_true("figure11" %in% names(vis$summary_map))
  expect_match(paste(vis$warning_map$figure10, collapse = " "), "Threshold profile: standard", fixed = TRUE)
  expect_match(paste(vis$summary_map$figure10, collapse = " "), "Literature bands", fixed = TRUE)

  vis_strict <- mfrmr::build_visual_summaries(fit, diag, threshold_profile = "strict")
  expect_match(paste(vis_strict$warning_map$figure10, collapse = " "), "Threshold profile: strict", fixed = TRUE)

  apa <- mfrmr::build_apa_outputs(fit, diag)
  expect_match(apa$table_figure_notes, "Figure 10", fixed = TRUE)
  expect_match(apa$table_figure_notes, "Figure 11", fixed = TRUE)
  expect_match(apa$table_figure_captions, "Residual PCA Scree", fixed = TRUE)
  expect_match(apa$report_text, "Literature bands", fixed = TRUE)

  profiles <- mfrmr::mfrm_threshold_profiles()
  expect_true("profiles" %in% names(profiles))
  expect_true(all(c("strict", "standard", "lenient") %in% names(profiles$profiles)))
})

test_that("packaged extdata includes baseline and iterative-calibrated CSVs", {
  ext <- system.file("extdata", package = "mfrmr")
  expect_true(nzchar(ext))

  files <- sort(list.files(ext))
  expect_true("eckes_jin_2021_study1_sim.csv" %in% files)
  expect_true("eckes_jin_2021_study2_sim.csv" %in% files)
  expect_true("eckes_jin_2021_study1_itercal_sim.csv" %in% files)
  expect_true("eckes_jin_2021_study2_itercal_sim.csv" %in% files)
})

test_that("packaged data objects are available via data() and loader", {
  expect_true("study1" %in% mfrmr::list_mfrmr_data())
  expect_true("combined_itercal" %in% mfrmr::list_mfrmr_data())

  data("ej2021_study1", package = "mfrmr", envir = environment())
  expect_true(exists("ej2021_study1"))
  expect_true(is.data.frame(ej2021_study1))
  expect_true(all(c("Study", "Person", "Rater", "Criterion", "Score") %in% names(ej2021_study1)))

  d2 <- mfrmr::load_mfrmr_data("study1")
  expect_true(is.data.frame(d2))
  expect_equal(nrow(d2), nrow(ej2021_study1))
})
