test_that("core fit/diagnostics workflow runs", {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

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
  fit_summary <- summary(fit)
  expect_s3_class(fit_summary, "summary.mfrm_fit")
  expect_true(all(c("overview", "facet_overview", "person_overview", "step_overview") %in% names(fit_summary)))
  printed_summary <- capture.output(summary(fit))
  expect_true(any(grepl("Many-Facet Rasch Model Summary", printed_summary, fixed = TRUE)))
  expect_true(any(grepl("Facet overview", printed_summary, fixed = TRUE)))
  p_fit_default <- plot(fit, draw = FALSE)
  expect_s3_class(p_fit_default, "mfrm_plot_bundle")
  expect_true(all(c("wright_map", "pathway_map", "category_characteristic_curves") %in% names(p_fit_default)))
  expect_s3_class(p_fit_default$wright_map, "mfrm_plot_data")
  expect_s3_class(p_fit_default$pathway_map, "mfrm_plot_data")
  expect_s3_class(p_fit_default$category_characteristic_curves, "mfrm_plot_data")
  printed_bundle <- capture.output(print(p_fit_default))
  expect_true(any(grepl("mfrm plot bundle", printed_bundle, fixed = TRUE)))

  p_fit_wright <- plot(fit, type = "wright", draw = FALSE)
  p_fit_pathway <- plot(fit, type = "pathway", draw = FALSE)
  p_fit_ccc <- plot(fit, type = "ccc", draw = FALSE)
  p_fit_person <- plot(fit, type = "person", draw = FALSE)
  p_fit_step <- plot(fit, type = "step", draw = FALSE)
  expect_s3_class(p_fit_wright, "mfrm_plot_data")
  expect_s3_class(p_fit_pathway, "mfrm_plot_data")
  expect_s3_class(p_fit_ccc, "mfrm_plot_data")
  expect_s3_class(p_fit_person, "mfrm_plot_data")
  expect_s3_class(p_fit_step, "mfrm_plot_data")

  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "both", pca_max_factors = 4)
  expect_s3_class(diag, "mfrm_diagnostics")
  expect_true("overall_fit" %in% names(diag))
  expect_true("residual_pca_overall" %in% names(diag))
  expect_true("residual_pca_by_facet" %in% names(diag))
  expect_true("unexpected" %in% names(diag))
  expect_true("fair_average" %in% names(diag))
  expect_true("displacement" %in% names(diag))
  expect_true("interrater" %in% names(diag))
  expect_true("facets_chisq" %in% names(diag))
  expect_true(is.data.frame(diag$unexpected$table))
  expect_true(is.data.frame(diag$fair_average$stacked))
  expect_true(is.data.frame(diag$displacement$table))
  expect_true(is.data.frame(diag$interrater$pairs))
  expect_true(is.data.frame(diag$facets_chisq))
  diag_summary <- summary(diag)
  expect_s3_class(diag_summary, "summary.mfrm_diagnostics")
  expect_true(all(c("overview", "overall_fit", "reliability", "top_fit", "flags") %in% names(diag_summary)))
  printed_diag <- capture.output(summary(diag))
  expect_true(any(grepl("Many-Facet Rasch Diagnostics Summary", printed_diag, fixed = TRUE)))

  t4 <- mfrmr::unexpected_response_table(fit, diagnostics = diag, abs_z_min = 1.5, prob_max = 0.4, top_n = 15)
  expect_s3_class(t4, "mfrm_unexpected")
  expect_true(all(c("table", "summary", "thresholds") %in% names(t4)))
  expect_true(is.data.frame(t4$table))
  expect_true(is.data.frame(t4$summary))
  t4_summary <- summary(t4)
  expect_s3_class(t4_summary, "summary.mfrm_bundle")
  t4_plot <- plot(t4, draw = FALSE)
  expect_s3_class(t4_plot, "mfrm_plot_data")

  t12 <- mfrmr::fair_average_table(fit, diagnostics = diag, udecimals = 2)
  expect_s3_class(t12, "mfrm_fair_average")
  expect_true(all(c("by_facet", "stacked", "raw_by_facet", "settings") %in% names(t12)))
  expect_true(is.data.frame(t12$stacked))
  expect_gt(nrow(t12$stacked), 0)
  t12_summary <- summary(t12)
  expect_s3_class(t12_summary, "summary.mfrm_bundle")
  t12_plot <- plot(t12, draw = FALSE)
  expect_s3_class(t12_plot, "mfrm_plot_data")

  disp <- mfrmr::displacement_table(fit, diagnostics = diag, anchored_only = FALSE)
  expect_s3_class(disp, "mfrm_displacement")
  expect_true(all(c("table", "summary", "thresholds") %in% names(disp)))
  expect_true(is.data.frame(disp$table))
  expect_true(is.data.frame(disp$summary))
  disp_summary <- summary(disp)
  expect_s3_class(disp_summary, "summary.mfrm_bundle")
  disp_plot <- plot(disp, draw = FALSE)
  expect_s3_class(disp_plot, "mfrm_plot_data")

  t5 <- mfrmr::measurable_summary_table(fit, diagnostics = diag)
  expect_s3_class(t5, "mfrm_measurable")
  expect_true(all(c("summary", "facet_coverage", "category_stats", "subsets") %in% names(t5)))
  expect_true(is.data.frame(t5$summary))
  expect_true(is.data.frame(t5$category_stats))
  t5_summary <- summary(t5)
  expect_s3_class(t5_summary, "summary.mfrm_bundle")

  t1 <- mfrmr::specifications_report(fit, title = "Toy run")
  expect_s3_class(t1, "mfrm_specifications")
  expect_true(all(c("header", "data_spec", "facet_labels", "output_spec", "convergence_control", "anchor_summary") %in% names(t1)))
  expect_true(is.data.frame(t1$header))
  expect_true(is.data.frame(t1$facet_labels))
  t1_summary <- summary(t1)
  expect_s3_class(t1_summary, "summary.mfrm_bundle")
  expect_false("fixed" %in% names(t1))
  t1_fixed <- mfrmr::specifications_report(fit, title = "Toy run", include_fixed = TRUE)
  expect_true("fixed" %in% names(t1_fixed))
  expect_true(is.character(t1_fixed$fixed))
  t1_alias <- mfrmr::specifications_report(fit, title = "Toy run")
  expect_true(is.data.frame(t1_alias$header))
  expect_equal(names(t1_alias), names(t1))

  t2 <- mfrmr::data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score"
  )
  expect_s3_class(t2, "mfrm_data_quality")
  expect_true(all(c("summary", "model_match", "row_audit", "unknown_elements", "category_counts") %in% names(t2)))
  expect_true(is.data.frame(t2$summary))
  expect_true(is.data.frame(t2$model_match))
  t2_summary <- summary(t2)
  expect_s3_class(t2_summary, "summary.mfrm_bundle")
  expect_false("fixed" %in% names(t2))
  t2_fixed <- mfrmr::data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    include_fixed = TRUE
  )
  expect_true("fixed" %in% names(t2_fixed))
  expect_true(is.character(t2_fixed$fixed))
  t2_alias <- mfrmr::data_quality_report(
    fit,
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score"
  )
  expect_true(is.data.frame(t2_alias$summary))
  expect_equal(names(t2_alias), names(t2))

  t3 <- mfrmr::estimation_iteration_report(fit, max_iter = 3, include_prox = TRUE)
  expect_s3_class(t3, "mfrm_iteration_report")
  expect_true(all(c("table", "summary", "settings") %in% names(t3)))
  expect_true(is.data.frame(t3$table))
  expect_true(is.data.frame(t3$summary))
  t3_summary <- summary(t3)
  expect_s3_class(t3_summary, "summary.mfrm_bundle")
  expect_false("fixed" %in% names(t3))
  t3_fixed <- mfrmr::estimation_iteration_report(fit, max_iter = 3, include_prox = TRUE, include_fixed = TRUE)
  expect_true("fixed" %in% names(t3_fixed))
  expect_true(is.character(t3_fixed$fixed))
  t3_alias <- mfrmr::estimation_iteration_report(fit, max_iter = 3, include_prox = TRUE)
  expect_true(is.data.frame(t3_alias$table))
  expect_equal(names(t3_alias), names(t3))

  t6 <- mfrmr::subset_connectivity_report(fit, diagnostics = diag)
  expect_s3_class(t6, "mfrm_subset_connectivity")
  expect_true(all(c("summary", "listing", "nodes", "settings") %in% names(t6)))
  expect_true(is.data.frame(t6$summary))
  expect_true(is.data.frame(t6$listing))
  expect_true(is.data.frame(t6$nodes))
  expect_true("ObservationPercent" %in% names(t6$summary))
  t6_summary <- summary(t6)
  expect_s3_class(t6_summary, "summary.mfrm_bundle")
  t6_alias <- mfrmr::subset_connectivity_report(fit, diagnostics = diag)
  expect_true(is.data.frame(t6_alias$summary))
  expect_equal(names(t6_alias), names(t6))

  t62 <- mfrmr::facet_statistics_report(fit, diagnostics = diag)
  expect_s3_class(t62, "mfrm_facet_statistics")
  expect_true(all(c("table", "ranges", "settings") %in% names(t62)))
  expect_true(is.data.frame(t62$table))
  expect_true(is.data.frame(t62$ranges))
  expect_true("Ruler" %in% names(t62$table))
  t62_summary <- summary(t62)
  expect_s3_class(t62_summary, "summary.mfrm_bundle")
  t62_alias <- mfrmr::facet_statistics_report(fit, diagnostics = diag)
  expect_true(is.data.frame(t62_alias$table))
  expect_equal(names(t62_alias), names(t62))

  t8 <- mfrmr::rating_scale_table(fit, diagnostics = diag)
  expect_s3_class(t8, "mfrm_rating_scale")
  expect_true(all(c("category_table", "threshold_table", "summary") %in% names(t8)))
  expect_true(is.data.frame(t8$category_table))
  expect_true(is.data.frame(t8$summary))
  t8_summary <- summary(t8)
  expect_s3_class(t8_summary, "summary.mfrm_bundle")
  t8_plot <- plot(t8, draw = FALSE)
  expect_s3_class(t8_plot, "mfrm_plot_data")

  t8b <- mfrmr::category_structure_report(fit, diagnostics = diag)
  expect_s3_class(t8b, "mfrm_category_structure")
  expect_true(all(c("category_table", "mode_peaks", "mode_boundaries", "median_thresholds", "mean_halfscore_points", "settings") %in% names(t8b)))
  expect_true(is.data.frame(t8b$category_table))
  expect_true(is.data.frame(t8b$mode_peaks))
  expect_true(is.data.frame(t8b$median_thresholds))
  t8b_summary <- summary(t8b)
  expect_s3_class(t8b_summary, "summary.mfrm_bundle")
  t8b_plot <- plot(t8b, draw = FALSE)
  expect_s3_class(t8b_plot, "mfrm_plot_data")
  expect_false("fixed" %in% names(t8b))
  t8b_fixed <- mfrmr::category_structure_report(fit, diagnostics = diag, include_fixed = TRUE)
  expect_true("fixed" %in% names(t8b_fixed))
  expect_true(is.character(t8b_fixed$fixed))
  t8b_alias <- mfrmr::category_structure_report(fit, diagnostics = diag)
  expect_true(is.data.frame(t8b_alias$category_table))
  expect_equal(names(t8b_alias), names(t8b))

  t8c <- mfrmr::category_curves_report(fit, theta_points = 101)
  expect_s3_class(t8c, "mfrm_category_curves")
  expect_true(all(c("graphfile", "graphfile_syntactic", "probabilities", "expected_ogive", "settings") %in% names(t8c)))
  expect_true(is.data.frame(t8c$graphfile))
  expect_true(is.data.frame(t8c$graphfile_syntactic))
  expect_true(any(grepl("^Prob:", names(t8c$graphfile))))
  expect_true(any(grepl("^Prob_", names(t8c$graphfile_syntactic))))
  t8c_summary <- summary(t8c)
  expect_s3_class(t8c_summary, "summary.mfrm_bundle")
  t8c_plot <- plot(t8c, draw = FALSE)
  expect_s3_class(t8c_plot, "mfrm_plot_data")
  expect_false("fixed" %in% names(t8c))
  t8c_fixed <- mfrmr::category_curves_report(fit, theta_points = 101, include_fixed = TRUE)
  expect_true("fixed" %in% names(t8c_fixed))
  expect_true(is.character(t8c_fixed$fixed))
  t8c_alias <- mfrmr::category_curves_report(fit, theta_points = 101)
  expect_true(is.data.frame(t8c_alias$graphfile))
  expect_equal(names(t8c_alias), names(t8c))

  of <- mfrmr::facets_output_file_bundle(
    fit,
    diagnostics = diag,
    include = c("graph", "score"),
    theta_points = 81
  )
  expect_s3_class(of, "mfrm_output_bundle")
  expect_true(all(c("graphfile", "graphfile_syntactic", "scorefile", "settings") %in% names(of)))
  expect_true(is.data.frame(of$graphfile))
  expect_true(is.data.frame(of$scorefile))
  expect_true("ObsProb" %in% names(of$scorefile))
  of_summary <- summary(of)
  expect_s3_class(of_summary, "summary.mfrm_bundle")
  expect_false("graphfile_fixed" %in% names(of))
  expect_false("scorefile_fixed" %in% names(of))

  of_fixed <- mfrmr::facets_output_file_bundle(
    fit,
    diagnostics = diag,
    include = c("graph", "score"),
    theta_points = 81,
    include_fixed = TRUE
  )
  expect_true("graphfile_fixed" %in% names(of_fixed))
  expect_true("scorefile_fixed" %in% names(of_fixed))
  expect_true(is.character(of_fixed$graphfile_fixed))
  expect_true(is.character(of_fixed$scorefile_fixed))

  write_dir <- tempfile("mfrmr-output-")
  dir.create(write_dir, recursive = TRUE)
  of_written <- mfrmr::facets_output_file_bundle(
    fit,
    diagnostics = diag,
    include = c("graph", "score"),
    theta_points = 61,
    write_files = TRUE,
    output_dir = write_dir,
    file_prefix = "toy",
    overwrite = TRUE
  )
  expect_true("written_files" %in% names(of_written))
  expect_true(is.data.frame(of_written$written_files))
  expect_gt(nrow(of_written$written_files), 0)
  expect_true(all(file.exists(of_written$written_files$Path)))

  bias <- mfrmr::estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion", max_iter = 2)
  expect_s3_class(bias, "mfrm_bias")
  expect_true(all(c("table", "summary", "chi_sq", "facet_a", "facet_b") %in% names(bias)))
  expect_true(is.data.frame(bias$table))
  bias_summary <- summary(bias)
  expect_s3_class(bias_summary, "summary.mfrm_bias")
  expect_true(all(c("overview", "chi_sq", "final_iteration", "top_rows") %in% names(bias_summary)))
  printed_bias <- capture.output(summary(bias))
  expect_true(any(grepl("Many-Facet Rasch Bias Summary", printed_bias, fixed = TRUE)))

  t13 <- mfrmr::bias_interaction_report(bias, top_n = 20)
  expect_s3_class(t13, "mfrm_bias_interaction")
  expect_true(all(c("ranked_table", "scatter_data", "summary", "thresholds", "facet_a", "facet_b") %in% names(t13)))
  expect_true(is.data.frame(t13$ranked_table))
  expect_true(is.data.frame(t13$scatter_data))
  t13_summary <- summary(t13)
  expect_s3_class(t13_summary, "summary.mfrm_bundle")
  t13_plot <- plot(t13, draw = FALSE)
  expect_s3_class(t13_plot, "mfrm_plot_data")
  t13_alias <- mfrmr::bias_interaction_report(bias, top_n = 20)
  expect_true(is.data.frame(t13_alias$ranked_table))
  expect_equal(names(t13_alias), names(t13))

  p13 <- mfrmr::plot_bias_interaction(t13, plot = "scatter", draw = FALSE)
  expect_s3_class(p13, "mfrm_plot_data")
  expect_equal(p13$name, "table13_bias")
  p13_alias <- mfrmr::plot_bias_interaction(t13, plot = "scatter", draw = FALSE)
  expect_s3_class(p13_alias, "mfrm_plot_data")
  expect_equal(p13_alias$name, "table13_bias")

  t11 <- mfrmr::bias_count_table(bias, min_count_warn = 1)
  expect_s3_class(t11, "mfrm_bias_count")
  expect_equal(t11$branch, "original")
  expect_true(all(c("table", "by_facet_a", "by_facet_b", "summary", "thresholds") %in% names(t11)))
  expect_true(is.data.frame(t11$table))
  expect_true(is.data.frame(t11$summary))
  t11_summary <- summary(t11)
  expect_s3_class(t11_summary, "summary.mfrm_bundle")
  t11_plot <- plot(t11, draw = FALSE)
  expect_s3_class(t11_plot, "mfrm_plot_data")

  t11_facets <- mfrmr::bias_count_table(bias, min_count_warn = 1, branch = "facets", fit = fit)
  expect_s3_class(t11_facets, "mfrm_bias_count")
  expect_equal(t11_facets$branch, "facets")
  expect_equal(t11_facets$style, "facets_manual")
  expect_true("Observd Count" %in% names(t11_facets$table))
  t11_facets_summary <- summary(t11_facets)
  expect_s3_class(t11_facets_summary, "summary.mfrm_bundle")

  t10 <- mfrmr::unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20)
  expect_s3_class(t10, "mfrm_unexpected_after_bias")
  expect_true(all(c("table", "summary", "thresholds", "facets") %in% names(t10)))
  expect_true(is.data.frame(t10$table))
  expect_true(is.data.frame(t10$summary))
  t10_summary <- summary(t10)
  expect_s3_class(t10_summary, "summary.mfrm_bundle")

  bias3 <- mfrmr::estimate_bias(
    fit,
    diag,
    interaction_facets = c("Rater", "Task", "Criterion"),
    max_iter = 2
  )
  expect_s3_class(bias3, "mfrm_bias")
  expect_equal(as.integer(bias3$interaction_order), 3L)
  expect_equal(as.character(bias3$interaction_mode), "higher_order")
  expect_true(all(c("Facet1_Level", "Facet2_Level", "Facet3_Level") %in% names(bias3$table)))
  bias3_summary <- summary(bias3)
  expect_s3_class(bias3_summary, "summary.mfrm_bias")
  expect_true("Pair" %in% names(bias3_summary$top_rows))

  t13_3 <- mfrmr::bias_interaction_report(bias3, top_n = 20)
  expect_s3_class(t13_3, "mfrm_bias_interaction")
  expect_true("facet_profile" %in% names(t13_3))
  expect_true(is.data.frame(t13_3$facet_profile))
  expect_true(all(c("Facet", "Level", "MeanAbsBias", "FlagRate") %in% names(t13_3$facet_profile)))

  p13_profile <- mfrmr::plot_bias_interaction(t13_3, plot = "facet_profile", draw = FALSE)
  expect_s3_class(p13_profile, "mfrm_plot_data")
  expect_equal(p13_profile$data$plot, "facet_profile")

  t10_3 <- mfrmr::unexpected_after_bias_table(fit, bias3, diagnostics = diag, top_n = 20)
  expect_s3_class(t10_3, "mfrm_unexpected_after_bias")
  expect_true(is.data.frame(t10_3$table))

  fixed3 <- mfrmr::build_fixed_reports(bias3)
  expect_s3_class(fixed3, "mfrm_fixed_reports")
  expect_true(is.character(fixed3$pairwise_fixed))
  expect_true(grepl("2-way interactions", fixed3$pairwise_fixed, fixed = TRUE))
  expect_true(is.data.frame(fixed3$pairwise_table))
  expect_equal(nrow(fixed3$pairwise_table), 0)
  fixed3_summary <- summary(fixed3)
  expect_s3_class(fixed3_summary, "summary.mfrm_bundle")

  fixed_facets <- mfrmr::build_fixed_reports(bias, branch = "facets")
  expect_s3_class(fixed_facets, "mfrm_fixed_reports")
  expect_equal(fixed_facets$branch, "facets")
  fixed_original <- mfrmr::build_fixed_reports(bias, branch = "original")
  expect_s3_class(fixed_original, "mfrm_fixed_reports")
  expect_equal(fixed_original$branch, "original")
  fixed_plot <- plot(fixed_original, draw = FALSE)
  expect_s3_class(fixed_plot, "mfrm_plot_data")

  ir <- mfrmr::interrater_agreement_table(fit, diagnostics = diag, rater_facet = "Rater")
  expect_s3_class(ir, "mfrm_interrater")
  expect_true(all(c("summary", "pairs", "settings") %in% names(ir)))
  expect_true(is.data.frame(ir$summary))
  expect_true(is.data.frame(ir$pairs))
  ir_summary <- summary(ir)
  expect_s3_class(ir_summary, "summary.mfrm_bundle")
  ir_plot <- plot(ir, draw = FALSE)
  expect_s3_class(ir_plot, "mfrm_plot_data")

  chi <- mfrmr::facets_chisq_table(fit, diagnostics = diag)
  expect_s3_class(chi, "mfrm_facets_chisq")
  expect_true(all(c("table", "summary", "thresholds") %in% names(chi)))
  expect_true(is.data.frame(chi$table))
  expect_true(is.data.frame(chi$summary))
  chi_summary <- summary(chi)
  expect_s3_class(chi_summary, "summary.mfrm_bundle")
  chi_plot <- plot(chi, draw = FALSE)
  expect_s3_class(chi_plot, "mfrm_plot_data")

  p_unexp <- mfrmr::plot_unexpected(fit, diagnostics = diag, abs_z_min = 1.5, prob_max = 0.4, top_n = 10, draw = FALSE)
  p_unexp2 <- mfrmr::plot_unexpected(t4, draw = FALSE)
  p_fair <- mfrmr::plot_fair_average(fit, diagnostics = diag, draw = FALSE)
  p_fair2 <- mfrmr::plot_fair_average(t12, draw = FALSE)
  p_disp <- mfrmr::plot_displacement(fit, diagnostics = diag, anchored_only = FALSE, draw = FALSE)
  p_disp2 <- mfrmr::plot_displacement(disp, draw = FALSE)
  p_ir <- mfrmr::plot_interrater_agreement(fit, diagnostics = diag, rater_facet = "Rater", draw = FALSE)
  p_ir2 <- mfrmr::plot_interrater_agreement(ir, draw = FALSE)
  p_chi <- mfrmr::plot_facets_chisq(fit, diagnostics = diag, draw = FALSE)
  p_chi2 <- mfrmr::plot_facets_chisq(chi, draw = FALSE)
  p_qc <- mfrmr::plot_qc_dashboard(fit, diagnostics = diag, draw = FALSE, top_n = 10)
  expect_s3_class(p_unexp, "mfrm_plot_data")
  expect_s3_class(p_unexp2, "mfrm_plot_data")
  expect_s3_class(p_fair, "mfrm_plot_data")
  expect_s3_class(p_fair2, "mfrm_plot_data")
  expect_s3_class(p_disp, "mfrm_plot_data")
  expect_s3_class(p_disp2, "mfrm_plot_data")
  expect_s3_class(p_ir, "mfrm_plot_data")
  expect_s3_class(p_ir2, "mfrm_plot_data")
  expect_s3_class(p_chi, "mfrm_plot_data")
  expect_s3_class(p_chi2, "mfrm_plot_data")
  expect_s3_class(p_qc, "mfrm_plot_data")
  expect_true(all(c("unexpected", "fair_average", "displacement", "interrater", "facets_chisq", "reliability") %in% names(p_qc$data)))

  pca <- mfrmr::analyze_residual_pca(diag, mode = "both", pca_max_factors = 4)
  pca_from_fit <- mfrmr::analyze_residual_pca(fit, mode = "both", pca_max_factors = 4)
  expect_true(is.data.frame(pca$overall_table))
  expect_true(is.data.frame(pca$by_facet_table))
  expect_gt(nrow(pca$overall_table), 0)
  expect_true(is.data.frame(pca_from_fit$overall_table))
  expect_gt(nrow(pca_from_fit$overall_table), 0)

  p1 <- mfrmr::plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
  p2 <- mfrmr::plot_residual_pca(pca, mode = "facet", facet = "Rater", plot_type = "loadings", top_n = 5, draw = FALSE)
  p3 <- mfrmr::plot_residual_pca(fit, mode = "overall", plot_type = "scree", draw = FALSE)
  expect_s3_class(p1, "mfrm_plot_data")
  expect_s3_class(p2, "mfrm_plot_data")
  expect_s3_class(p3, "mfrm_plot_data")

  vis <- mfrmr::build_visual_summaries(fit, diag)
  expect_s3_class(vis, "mfrm_visual_summaries")
  expect_true("residual_pca_overall" %in% names(vis$warning_map))
  expect_true("residual_pca_by_facet" %in% names(vis$warning_map))
  expect_true("residual_pca_overall" %in% names(vis$summary_map))
  expect_true("residual_pca_by_facet" %in% names(vis$summary_map))
  expect_match(paste(vis$warning_map$residual_pca_overall, collapse = " "), "Threshold profile: standard", fixed = TRUE)
  expect_match(paste(vis$summary_map$residual_pca_overall, collapse = " "), "Literature bands", fixed = TRUE)
  vis_summary <- summary(vis)
  expect_s3_class(vis_summary, "summary.mfrm_bundle")
  vis_plot <- plot(vis, draw = FALSE)
  expect_s3_class(vis_plot, "mfrm_plot_data")

  vis_strict <- mfrmr::build_visual_summaries(fit, diag, threshold_profile = "strict")
  expect_match(paste(vis_strict$warning_map$residual_pca_overall, collapse = " "), "Threshold profile: strict", fixed = TRUE)
  vis_facets <- mfrmr::build_visual_summaries(fit, diag, branch = "facets")
  expect_s3_class(vis_facets, "mfrm_visual_summaries")
  expect_equal(vis_facets$branch, "facets")
  expect_true(all(c("Visual", "FACETS") %in% names(vis_facets$crosswalk)))

  apa <- mfrmr::build_apa_outputs(fit, diag)
  expect_s3_class(apa, "mfrm_apa_outputs")
  expect_s3_class(apa$report_text, "mfrm_apa_text")
  expect_match(apa$table_figure_notes, "Residual PCA scree", fixed = TRUE)
  expect_match(apa$table_figure_notes, "Residual PCA by facet", fixed = TRUE)
  expect_match(apa$table_figure_captions, "Residual PCA Scree", fixed = TRUE)
  expect_match(apa$report_text, "Literature bands", fixed = TRUE)
  expect_match(apa$report_text, "Optimization", fixed = TRUE)
  expect_match(apa$report_text, "Constraint settings:", fixed = TRUE)
  expect_match(apa$report_text, "Step/threshold summary:", fixed = TRUE)
  expect_match(apa$report_text, "Largest misfit", fixed = TRUE)
  printed_apa_text <- capture.output(print(apa$report_text))
  expect_true(any(grepl("Method\\.", printed_apa_text)))
  apa_summary <- summary(apa)
  expect_s3_class(apa_summary, "summary.mfrm_apa_outputs")
  expect_true(is.data.frame(apa_summary$overview))
  expect_true("report_text" %in% apa_summary$components$Component)
  apa_wrapped <- mfrmr::build_apa_outputs(fit, diag, context = list(line_width = 60))
  expect_true(grepl("Method\\.\\n\\n", apa_wrapped$report_text))
  expect_true(grepl("\\n\\nResults\\.\\n\\n", apa_wrapped$report_text))

  at <- mfrmr::apa_table(fit, which = "summary", caption = "Model summary")
  expect_s3_class(at, "apa_table")
  expect_true(is.data.frame(at$table))
  printed_table <- capture.output(print(at))
  expect_true(any(grepl("Model summary", printed_table, fixed = TRUE)))
  at_summary <- summary(at)
  expect_s3_class(at_summary, "summary.apa_table")
  at_plot <- plot(at, draw = FALSE)
  expect_s3_class(at_plot, "mfrm_plot_data")
  at_ir <- mfrmr::apa_table(fit, which = "interrater_pairs", diagnostics = diag)
  expect_s3_class(at_ir, "apa_table")
  expect_true(is.data.frame(at_ir$table))
  at_facets <- mfrmr::apa_table(fit, which = "summary", branch = "facets")
  expect_s3_class(at_facets, "apa_table")
  expect_equal(at_facets$branch, "facets")
  expect_match(at_facets$caption, "FACETS-aligned table", fixed = TRUE)

  profiles <- mfrmr::mfrm_threshold_profiles()
  expect_s3_class(profiles, "mfrm_threshold_profiles")
  expect_true("profiles" %in% names(profiles))
  expect_true(all(c("strict", "standard", "lenient") %in% names(profiles$profiles)))
  profiles_summary <- summary(profiles)
  expect_s3_class(profiles_summary, "summary.mfrm_threshold_profiles")
  expect_true("ThresholdCount" %in% names(profiles_summary$overview))
  expect_true("standard" %in% names(profiles_summary$thresholds))
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

test_that("legacy numbered API names are internal (not exported)", {
  old_opt <- options(lifecycle_verbosity = "quiet")
  on.exit(options(old_opt), add = TRUE)

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
    maxit = 20
  )

  expect_error(
    mfrmr::table1_specifications(fit),
    "not an exported object"
  )
  old_t1 <- mfrmr:::table1_specifications(fit)
  new_t1 <- mfrmr::specifications_report(fit)
  expect_equal(names(old_t1), names(new_t1))

  expect_error(
    mfrmr::table8_curves_export(fit, theta_points = 101),
    "not an exported object"
  )
  old_t8 <- mfrmr:::table8_curves_export(fit, theta_points = 101)
  new_t8 <- mfrmr::category_curves_report(fit, theta_points = 101)
  expect_equal(names(old_t8), names(new_t8))
})

test_that("descriptive and anchor-audit helpers run", {
  toy <- expand.grid(
    Person = paste0("P", 1:6),
    Rater = paste0("R", 1:3),
    Criterion = c("Content", "Language"),
    stringsAsFactors = FALSE
  )
  toy$Score <- (
    as.integer(factor(toy$Person)) +
      as.integer(factor(toy$Rater)) +
      as.integer(factor(toy$Criterion))
  ) %% 3

  ds <- mfrmr::describe_mfrm_data(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score"
  )
  expect_s3_class(ds, "mfrm_data_description")
  expect_true(all(c("overview", "score_distribution", "linkage_summary") %in% names(ds)))
  expect_gt(nrow(ds$score_distribution), 0)
  expect_true("agreement" %in% names(ds))
  expect_true(is.list(ds$agreement))
  expect_true(is.data.frame(ds$agreement$summary))
  expect_true(is.data.frame(ds$agreement$pairs))
  expect_gt(nrow(ds$agreement$summary), 0)

  anchors <- data.frame(
    Facet = c("Rater", "Rater", "Rater", "UnknownFacet"),
    Level = c("R1", "R1", "R999", "X1"),
    Anchor = c(0.0, 0.1, 0.2, 0.3),
    stringsAsFactors = FALSE
  )
  group_anchors <- data.frame(
    Facet = c("Rater", "Rater", "Rater"),
    Level = c("R2", "R2", "R3"),
    Group = c("G1", "G2", "G1"),
    GroupValue = c(0, NA, 0.2),
    stringsAsFactors = FALSE
  )

  aud <- mfrmr::audit_mfrm_anchors(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    anchors = anchors,
    group_anchors = group_anchors,
    min_common_anchors = 4,
    min_obs_per_element = 20,
    min_obs_per_category = 8
  )
  expect_s3_class(aud, "mfrm_anchor_audit")
  expect_true(is.data.frame(aud$anchors))
  expect_true(is.data.frame(aud$group_anchors))
  expect_true(is.list(aud$design_checks))
  expect_true(is.list(aud$thresholds))
  expect_equal(aud$thresholds$min_common_anchors, 4L)
  expect_equal(aud$thresholds$min_obs_per_element, 20)
  expect_equal(aud$thresholds$min_obs_per_category, 8)
  expect_true(is.data.frame(aud$design_checks$level_observation_summary))
  expect_true(is.data.frame(aud$design_checks$category_counts))
  expect_true(any(aud$issue_counts$Issue == "duplicate_anchors" & aud$issue_counts$N > 0))
  expect_true(any(aud$issue_counts$Issue == "unknown_anchor_facets" & aud$issue_counts$N > 0))
  expect_true(any(aud$issue_counts$Issue == "unknown_anchor_levels" & aud$issue_counts$N > 0))

  expect_error(
    mfrmr::fit_mfrm(
      data = toy,
      person = "Person",
      facets = c("Rater", "Criterion"),
      score = "Score",
      anchors = anchors,
      group_anchors = group_anchors,
      method = "JML",
      maxit = 15,
      min_common_anchors = 4,
      anchor_policy = "error"
    ),
    "Anchor audit detected"
  )

  fit <- mfrmr::fit_mfrm(
    data = toy,
    person = "Person",
    facets = c("Rater", "Criterion"),
    score = "Score",
    anchors = anchors,
    group_anchors = group_anchors,
    method = "JML",
    maxit = 15,
    min_common_anchors = 4,
    anchor_policy = "silent"
  )
  expect_s3_class(fit, "mfrm_fit")
  expect_true("anchor_audit" %in% names(fit$config))

  anchor_tbl <- mfrmr::make_anchor_table(fit)
  expect_true(is.data.frame(anchor_tbl))
  expect_true(all(c("Facet", "Level", "Anchor") %in% names(anchor_tbl)))
  expect_gt(nrow(anchor_tbl), 0)
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

test_that("MML + PCM path runs and returns diagnostics", {
  d <- mfrmr:::sample_mfrm_data(seed = 101)

  fit <- mfrmr::fit_mfrm(
    data = d,
    person = "Person",
    facets = c("Rater", "Task", "Criterion"),
    score = "Score",
    method = "MML",
    model = "PCM",
    maxit = 10,
    quad_points = 7
  )

  expect_s3_class(fit, "mfrm_fit")
  expect_true(is.data.frame(fit$summary))
  expect_identical(as.character(fit$summary$Method[1]), "MML")
  expect_identical(as.character(fit$summary$Model[1]), "PCM")

  diag <- mfrmr::diagnose_mfrm(fit, residual_pca = "none")
  expect_s3_class(diag, "mfrm_diagnostics")
  expect_true(is.data.frame(diag$overall_fit))
  expect_true(is.data.frame(diag$fit))
})
