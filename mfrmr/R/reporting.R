# Reporting and narrative helpers (FACETS-style / APA-style)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

py_style_format <- function(fmt, value) {
  if (is.na(value)) return("")
  if (is.null(fmt)) return(as.character(value))
  if (is.function(fmt)) return(as.character(fmt(value)))
  if (!is.character(fmt) || length(fmt) != 1) return(as.character(value))
  if (identical(fmt, "{}")) return(as.character(value))

  m <- regexec("^\\{:\\.(\\d+)f\\}$", fmt)
  g <- regmatches(fmt, m)[[1]]
  if (length(g) == 2) {
    digits <- as.integer(g[2])
    return(sprintf(paste0("%.", digits, "f"), as.numeric(value)))
  }
  as.character(value)
}

format_fixed_width_table <- function(df,
                                     columns,
                                     formats = list(),
                                     right_align = NULL,
                                     max_col_width = 16,
                                     min_col_width = 6) {
  if (is.null(df) || nrow(df) == 0) return("No data")

  if (is.null(right_align)) {
    right_align <- columns[columns %in% names(df)][vapply(columns[columns %in% names(df)], function(col) {
      is.numeric(df[[col]])
    }, logical(1))]
  }

  str_cols <- list()
  widths <- list()

  for (col in columns) {
    if (!col %in% names(df)) {
      vals <- rep("", nrow(df))
    } else {
      vals <- vapply(df[[col]], function(v) {
        if (is.na(v)) return("")
        py_style_format(formats[[col]], v)
      }, character(1))
    }

    str_cols[[col]] <- vals
    max_len <- max(c(nchar(col), nchar(vals)), na.rm = TRUE)
    widths[[col]] <- max(min_col_width, min(max_len, max_col_width))
  }

  pad <- function(col, text) {
    text <- substr(text, 1, widths[[col]])
    if (col %in% right_align) {
      stringr::str_pad(text, widths[[col]], side = "left")
    } else {
      stringr::str_pad(text, widths[[col]], side = "right")
    }
  }

  header <- paste(vapply(columns, function(col) pad(col, col), character(1)), collapse = " ")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    paste(vapply(columns, function(col) pad(col, str_cols[[col]][i]), character(1)), collapse = " ")
  }, character(1))

  paste(c(header, rows), collapse = "\n")
}

build_bias_fixed_text <- function(table_df,
                                  summary_df,
                                  chi_df,
                                  facet_a,
                                  facet_b,
                                  columns,
                                  formats) {
  if (is.null(table_df) || nrow(table_df) == 0) return("No bias data")

  fixed_table <- format_fixed_width_table(table_df, columns, formats = formats, max_col_width = 18)
  lines <- c(paste0("Bias/Interaction: ", facet_a, " x ", facet_b), "", fixed_table)

  if (!is.null(summary_df) && nrow(summary_df) > 0) {
    lines <- c(lines, "", "Summary")
    for (i in seq_len(nrow(summary_df))) {
      row <- summary_df[i, , drop = FALSE]
      line <- paste0(
        as.character(row$Statistic),
        ": Bias Size=", fmt_num(row$`Bias Size`, 2),
        ", Obs-Exp Avg=", fmt_num(row$`Obs-Exp Average`, 2),
        ", Model S.E.=", fmt_num(row$`S.E.`, 2)
      )
      lines <- c(lines, line)
    }
  }

  if (!is.null(chi_df) && nrow(chi_df) > 0) {
    chi <- chi_df[1, , drop = FALSE]
    line <- if (!is.na(chi$FixedChiSq)) {
      paste0(
        "Fixed (all = 0) chi-squared: ", sprintf("%.2f", as.numeric(chi$FixedChiSq)),
        "  d.f.: ", ifelse(is.na(chi$FixedDF), "", as.character(as.integer(round(chi$FixedDF)))),
        "  significance (probability): ", sprintf("%.4f", as.numeric(chi$FixedProb))
      )
    } else {
      "Fixed (all = 0) chi-squared: N/A"
    }
    lines <- c(lines, "", line)
  }

  paste(lines, collapse = "\n")
}

build_pairwise_fixed_text <- function(pair_df, target_facet, context_facet, columns, formats) {
  if (is.null(pair_df) || nrow(pair_df) == 0) return("No pairwise data")
  fixed_table <- format_fixed_width_table(pair_df, columns, formats = formats, max_col_width = 18)
  paste(
    paste0("Bias/Interaction Pairwise Report: Target=", target_facet, "  Context=", context_facet),
    "",
    fixed_table,
    sep = "\n"
  )
}

to_float <- function(value) {
  out <- suppressWarnings(as.numeric(value))
  ifelse(length(out) == 0, NA_real_, out)
}

fmt_count <- function(value) {
  val <- to_float(value)
  if (!is.finite(val)) return("NA")
  if (abs(val - round(val)) < 1e-6) return(as.character(as.integer(round(val))))
  sprintf("%.0f", val)
}

fmt_num <- function(value, decimals = 2) {
  val <- to_float(value)
  if (!is.finite(val)) return("NA")
  sprintf(paste0("%.", decimals, "f"), val)
}

fmt_pvalue <- function(value) {
  val <- to_float(value)
  if (!is.finite(val)) return("NA")
  if (val < 0.001) return("< .001")
  paste0("= ", sprintf("%.3f", val))
}

describe_series <- function(series) {
  if (is.null(series)) return(NULL)
  arr <- suppressWarnings(as.numeric(series))
  arr <- arr[is.finite(arr)]
  if (length(arr) == 0) return(NULL)
  list(
    min = min(arr),
    max = max(arr),
    mean = mean(arr),
    sd = if (length(arr) > 1) stats::sd(arr) else NA_real_
  )
}

safe_residual_pca <- function(diagnostics, mode = "both", pca_max_factors = 10L) {
  if (is.null(diagnostics)) return(NULL)
  tryCatch(
    analyze_residual_pca(
      diagnostics = diagnostics,
      mode = mode,
      pca_max_factors = pca_max_factors
    ),
    error = function(e) NULL
  )
}

warning_threshold_profiles <- function() {
  list(
    profiles = list(
      strict = list(
        n_obs_min = 200,
        n_person_min = 50,
        low_cat_min = 15,
        min_facet_levels = 4,
        misfit_ratio_warn = 0.08,
        missing_fit_ratio_warn = 0.15,
        zstd2_ratio_warn = 0.08,
        zstd3_ratio_warn = 0.03,
        expected_var_min = 0.30,
        pca_first_eigen_warn = 1.5,
        pca_first_prop_warn = 0.10
      ),
      standard = list(
        n_obs_min = 100,
        n_person_min = 30,
        low_cat_min = 10,
        min_facet_levels = 3,
        misfit_ratio_warn = 0.10,
        missing_fit_ratio_warn = 0.20,
        zstd2_ratio_warn = 0.10,
        zstd3_ratio_warn = 0.05,
        expected_var_min = 0.20,
        pca_first_eigen_warn = 2.0,
        pca_first_prop_warn = 0.10
      ),
      lenient = list(
        n_obs_min = 60,
        n_person_min = 20,
        low_cat_min = 5,
        min_facet_levels = 2,
        misfit_ratio_warn = 0.15,
        missing_fit_ratio_warn = 0.30,
        zstd2_ratio_warn = 0.15,
        zstd3_ratio_warn = 0.08,
        expected_var_min = 0.10,
        pca_first_eigen_warn = 3.0,
        pca_first_prop_warn = 0.20
      )
    ),
    pca_reference_bands = list(
      eigenvalue = c(
        critical_minimum = 1.4,
        caution = 1.5,
        common = 2.0,
        strong = 3.0
      ),
      proportion = c(
        minor = 0.05,
        caution = 0.10,
        strong = 0.20
      )
    )
  )
}

resolve_warning_thresholds <- function(thresholds = NULL, threshold_profile = "standard") {
  refs <- warning_threshold_profiles()
  profile_name <- tolower(as.character(threshold_profile %||% "standard"))
  if (!profile_name %in% names(refs$profiles)) profile_name <- "standard"
  active <- refs$profiles[[profile_name]]
  if (!is.null(thresholds) && length(thresholds) > 0) {
    active[names(thresholds)] <- thresholds
  }
  list(
    profile_name = profile_name,
    thresholds = active,
    pca_reference_bands = refs$pca_reference_bands
  )
}

build_pca_reference_text <- function(reference_bands) {
  eigen <- reference_bands$eigenvalue
  prop <- reference_bands$proportion
  paste0(
    "Literature bands: EV >= ", fmt_num(eigen[["critical_minimum"]], 1), " (critical minimum), ",
    ">= ", fmt_num(eigen[["caution"]], 1), " (caution), ",
    ">= ", fmt_num(eigen[["common"]], 1), " (common), ",
    ">= ", fmt_num(eigen[["strong"]], 1), " (strong); ",
    "variance >= ", fmt_num(100 * prop[["minor"]], 0), "% (minor), ",
    ">= ", fmt_num(100 * prop[["caution"]], 0), "% (caution), ",
    ">= ", fmt_num(100 * prop[["strong"]], 0), "% (strong)."
  )
}

build_pca_check_text <- function(eigenvalue, proportion, reference_bands) {
  eigen <- to_float(eigenvalue)
  prop <- to_float(proportion)
  if (!is.finite(eigen) && !is.finite(prop)) {
    return("Current PC1 threshold checks are unavailable.")
  }

  e <- reference_bands$eigenvalue
  p <- reference_bands$proportion
  checks <- c(
    paste0("EV>=", fmt_num(e[["caution"]], 1), ":", ifelse(is.finite(eigen) && eigen >= e[["caution"]], "Y", "N")),
    paste0("EV>=", fmt_num(e[["common"]], 1), ":", ifelse(is.finite(eigen) && eigen >= e[["common"]], "Y", "N")),
    paste0("EV>=", fmt_num(e[["strong"]], 1), ":", ifelse(is.finite(eigen) && eigen >= e[["strong"]], "Y", "N")),
    paste0("Var>=", fmt_num(100 * p[["caution"]], 0), "%:", ifelse(is.finite(prop) && prop >= p[["caution"]], "Y", "N")),
    paste0("Var>=", fmt_num(100 * p[["strong"]], 0), "%:", ifelse(is.finite(prop) && prop >= p[["strong"]], "Y", "N"))
  )
  paste0("Current PC1 checks: ", paste(checks, collapse = ", "), ".")
}

extract_overall_pca_first <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$overall_table) || nrow(pca_obj$overall_table) == 0) return(NULL)
  tbl <- pca_obj$overall_table
  idx <- which.min(suppressWarnings(as.numeric(tbl$Component)))
  if (!length(idx) || !is.finite(idx)) return(NULL)
  tbl[idx[1], , drop = FALSE]
}

extract_overall_pca_second <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$overall_table) || nrow(pca_obj$overall_table) == 0) return(NULL)
  tbl <- pca_obj$overall_table
  comp <- suppressWarnings(as.numeric(tbl$Component))
  idx <- which(comp == 2)
  if (length(idx) == 0) return(NULL)
  tbl[idx[1], , drop = FALSE]
}

extract_facet_pca_first <- function(pca_obj) {
  if (is.null(pca_obj) || is.null(pca_obj$by_facet_table) || nrow(pca_obj$by_facet_table) == 0) {
    return(data.frame())
  }

  tbl <- pca_obj$by_facet_table
  if (!"Facet" %in% names(tbl) || !"Component" %in% names(tbl)) return(data.frame())

  split_tbl <- split(tbl, as.character(tbl$Facet))
  out <- lapply(split_tbl, function(df) {
    comp <- suppressWarnings(as.numeric(df$Component))
    idx <- which.min(comp)
    if (!length(idx) || !is.finite(idx)) return(NULL)
    df[idx[1], , drop = FALSE]
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.frame())
  out <- dplyr::bind_rows(out)
  out[order(suppressWarnings(as.numeric(out$Eigenvalue)), decreasing = TRUE), , drop = FALSE]
}

build_apa_report_text <- function(res, diagnostics, bias_results = NULL, context = list(), whexact = FALSE) {
  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  prep <- res$prep
  config <- res$config

  n_obs <- if (!is.null(summary)) to_float(summary$N) else NA_real_
  n_person <- if (!is.null(summary)) to_float(summary$Persons) else nrow(res$facets$person)
  n_cat <- if (!is.null(summary)) to_float(summary$Categories) else to_float(config$n_cat)
  rating_min <- to_float(prep$rating_min)
  rating_max <- to_float(prep$rating_max)

  facet_names <- config$facet_names
  facet_counts <- vapply(facet_names, function(f) length(config$facet_levels[[f]]), numeric(1))
  facets_text <- if (length(facet_counts) > 0) {
    paste(paste0(names(facet_counts), " (n = ", fmt_count(facet_counts), ")"), collapse = ", ")
  } else {
    "no additional facets"
  }

  assessment <- trimws(as.character(context$assessment %||% ""))
  setting <- trimws(as.character(context$setting %||% ""))
  rater_training <- trimws(as.character(context$rater_training %||% ""))
  raters_per_response <- trimws(as.character(context$raters_per_response %||% ""))
  scale_desc <- trimws(as.character(context$scale_desc %||% ""))

  method_sentences <- character(0)
  if (nzchar(assessment)) {
    method_sentences <- c(method_sentences, if (nzchar(setting)) {
      paste0("The analysis focused on ", assessment, " in ", setting, ".")
    } else {
      paste0("The analysis focused on ", assessment, ".")
    })
  }

  method_sentences <- c(
    method_sentences,
    paste0(
      "A many-facet Rasch model (MFRM) was fit to ", fmt_count(n_obs),
      " observations from ", fmt_count(n_person),
      " persons scored on a ", fmt_count(n_cat),
      "-category scale (", fmt_count(rating_min), "-", fmt_count(rating_max), ")."
    ),
    if (length(facet_names) > 0) {
      paste0("The design included facets for ", facets_text, ".")
    } else {
      "No additional facets beyond Person were modeled."
    }
  )

  if (nzchar(scale_desc)) method_sentences <- c(method_sentences, paste0("The rating scale was described as ", scale_desc, "."))
  if (nzchar(rater_training)) method_sentences <- c(method_sentences, paste0("Raters received ", rater_training, "."))
  if (nzchar(raters_per_response)) method_sentences <- c(method_sentences, paste0("Each response was scored by ", raters_per_response, " raters on average."))

  model <- config$model
  method <- config$method
  model_sentence <- paste0("The ", model, " specification was estimated using ", method, " in the native R MFRM package.")
  if (identical(model, "PCM") && !is.null(config$step_facet) && nzchar(config$step_facet)) {
    model_sentence <- paste0(model_sentence, " The step structure varied by ", config$step_facet, ".")
  }
  method_sentences <- c(method_sentences, model_sentence)

  if (!is.null(config$weight_col) && nzchar(config$weight_col)) {
    method_sentences <- c(method_sentences, "Observation weights were applied as frequency counts.")
  }

  method_text <- paste0("Method.\n", paste(method_sentences, collapse = " "))

  results_sentences <- character(0)
  cat_tbl <- if (!is.null(diagnostics)) calc_category_stats(diagnostics$obs, res = res, whexact = whexact) else tibble::tibble()
  step_order <- calc_step_order(res$steps)
  unused <- if (nrow(cat_tbl) > 0) sum(cat_tbl$Count == 0, na.rm = TRUE) else 0
  low_count <- if (nrow(cat_tbl) > 0) sum(cat_tbl$Count < 10, na.rm = TRUE) else 0
  disordered <- if (nrow(step_order) > 0) step_order |> dplyr::filter(Ordered == FALSE) else tibble::tibble()
  usage_label <- if (unused == 0 && low_count == 0) "adequate" else "uneven"
  threshold_text <- if (nrow(disordered) == 0) "thresholds were ordered" else paste0("thresholds were disordered for ", fmt_count(nrow(disordered)), " step(s)")

  results_sentences <- c(results_sentences,
    paste0("Category usage was ", usage_label, " (unused categories = ", fmt_count(unused),
           ", low-count categories = ", fmt_count(low_count), "), and ", threshold_text, ".")
  )

  person_stats <- describe_series(res$facets$person$Estimate)
  if (!is.null(person_stats)) {
    results_sentences <- c(results_sentences,
      paste0("Person measures ranged from ", fmt_num(person_stats$min), " to ", fmt_num(person_stats$max),
             " logits (M = ", fmt_num(person_stats$mean), ", SD = ", fmt_num(person_stats$sd), ").")
    )
  }

  if (!is.null(res$facets$others) && nrow(res$facets$others) > 0) {
    for (facet in facet_names) {
      df_f <- res$facets$others |> dplyr::filter(Facet == facet)
      stats_f <- describe_series(df_f$Estimate)
      if (!is.null(stats_f)) {
        results_sentences <- c(results_sentences,
          paste0(facet, " measures ranged from ", fmt_num(stats_f$min), " to ", fmt_num(stats_f$max),
                 " logits (M = ", fmt_num(stats_f$mean), ", SD = ", fmt_num(stats_f$sd), ").")
        )
      }
    }
  }

  overall_fit <- if (!is.null(diagnostics$overall_fit) && nrow(diagnostics$overall_fit) > 0) diagnostics$overall_fit[1, , drop = FALSE] else NULL
  if (!is.null(overall_fit)) {
    infit <- to_float(overall_fit$Infit)
    outfit <- to_float(overall_fit$Outfit)
    fit_label <- if (is.finite(infit) && is.finite(outfit) && infit >= 0.5 && infit <= 1.5 && outfit >= 0.5 && outfit <= 1.5) "acceptable" else "elevated"
    results_sentences <- c(results_sentences,
      paste0("Overall fit was ", fit_label, " (infit MnSq = ", fmt_num(infit), ", outfit MnSq = ", fmt_num(outfit), ").")
    )
  }

  fit_tbl <- diagnostics$fit
  if (!is.null(fit_tbl) && nrow(fit_tbl) > 0) {
    misfit <- with(fit_tbl, (Infit < 0.5) | (Infit > 1.5) | (Outfit < 0.5) | (Outfit > 1.5))
    results_sentences <- c(results_sentences,
      paste0(fmt_count(sum(misfit, na.rm = TRUE)), " of ", fmt_count(nrow(fit_tbl)), " elements exceeded the 0.5-1.5 fit range.")
    )
  }

  rel_tbl <- diagnostics$reliability
  if (!is.null(rel_tbl) && nrow(rel_tbl) > 0) {
    rel_lines <- vapply(seq_len(nrow(rel_tbl)), function(i) {
      row <- rel_tbl[i, , drop = FALSE]
      paste0(row$Facet, " reliability = ", fmt_num(row$Reliability), " (separation = ", fmt_num(row$Separation), ").")
    }, character(1))
    results_sentences <- c(results_sentences, paste(rel_lines, collapse = " "))
  }

  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_overall_2 <- extract_overall_pca_second(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  pca_reference_text <- build_pca_reference_text(warning_threshold_profiles()$pca_reference_bands)

  if (!is.null(pca_overall_1)) {
    ev1 <- to_float(pca_overall_1$Eigenvalue)
    pr1 <- to_float(pca_overall_1$Proportion) * 100
    if (!is.null(pca_overall_2)) {
      ev2 <- to_float(pca_overall_2$Eigenvalue)
      results_sentences <- c(
        results_sentences,
        paste0(
          "Residual PCA (overall standardized residual matrix) showed PC1 eigenvalue = ",
          fmt_num(ev1), " (", fmt_num(pr1, 1), "% variance), with PC2 eigenvalue = ", fmt_num(ev2), "."
        )
      )
    } else {
      results_sentences <- c(
        results_sentences,
        paste0(
          "Residual PCA (overall standardized residual matrix) showed PC1 eigenvalue = ",
          fmt_num(ev1), " (", fmt_num(pr1, 1), "% variance)."
        )
      )
    }
  } else {
    results_sentences <- c(results_sentences, "Residual PCA was not available for this run.")
  }

  if (nrow(pca_facet_1) > 0) {
    top <- pca_facet_1[1, , drop = FALSE]
    results_sentences <- c(
      results_sentences,
      paste0(
        "Facet-specific residual PCA showed the largest first-component signal in ",
        as.character(top$Facet), " (eigenvalue = ", fmt_num(top$Eigenvalue),
        ", ", fmt_num(100 * to_float(top$Proportion), 1), "% variance)."
      )
    )
  }
  results_sentences <- c(results_sentences, pca_reference_text)

  if (!is.null(bias_results) && !is.null(bias_results$table) && nrow(bias_results$table) > 0) {
    bias_tbl <- bias_results$table |> dplyr::filter(is.finite(t))
    if (nrow(bias_tbl) > 0) {
      idx <- which.max(abs(bias_tbl$t))
      row <- bias_tbl[idx, , drop = FALSE]
      results_sentences <- c(results_sentences,
        paste0("Bias analysis for ", bias_results$facet_a, " x ", bias_results$facet_b,
               " showed a largest contrast of ", fmt_num(row$`Bias Size`),
               " logits (t = ", fmt_num(row$t), ", p ", fmt_pvalue(row$`Prob.`), ").")
      )
    }
  } else {
    results_sentences <- c(results_sentences, "Bias analysis was not estimated in this run.")
  }

  results_text <- paste0("Results.\n", paste(results_sentences, collapse = " "))
  paste0(method_text, "\n\n", results_text)
}

build_apa_table_figure_note_map <- function(res, diagnostics, bias_results = NULL, context = list(), whexact = FALSE) {
  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  prep <- res$prep
  config <- res$config

  n_obs <- if (!is.null(summary)) to_float(summary$N) else NA_real_
  n_person <- if (!is.null(summary)) to_float(summary$Persons) else nrow(res$facets$person)
  n_cat <- if (!is.null(summary)) to_float(summary$Categories) else to_float(config$n_cat)
  rating_min <- to_float(prep$rating_min)
  rating_max <- to_float(prep$rating_max)
  model <- config$model
  method <- config$method
  rater_facet <- trimws(as.character(context$rater_facet %||% ""))

  cat_tbl <- calc_category_stats(diagnostics$obs, res = res, whexact = whexact)
  step_order <- calc_step_order(res$steps)
  unused <- if (nrow(cat_tbl) > 0) sum(cat_tbl$Count == 0, na.rm = TRUE) else 0
  low_count <- if (nrow(cat_tbl) > 0) sum(cat_tbl$Count < 10, na.rm = TRUE) else 0
  disordered <- if (nrow(step_order) > 0) step_order |> dplyr::filter(Ordered == FALSE) else tibble::tibble()
  threshold_text <- if (nrow(disordered) == 0) "ordered" else paste0("disordered in ", fmt_count(nrow(disordered)), " step(s)")

  overall_fit <- if (!is.null(diagnostics$overall_fit) && nrow(diagnostics$overall_fit) > 0) diagnostics$overall_fit[1, , drop = FALSE] else NULL
  infit <- if (!is.null(overall_fit)) to_float(overall_fit$Infit) else NA_real_
  outfit <- if (!is.null(overall_fit)) to_float(overall_fit$Outfit) else NA_real_

  rel_tbl <- diagnostics$reliability
  rater_rel <- NULL
  if (nzchar(rater_facet) && !is.null(rel_tbl) && nrow(rel_tbl) > 0) {
    match <- rel_tbl |> dplyr::filter(Facet == rater_facet)
    if (nrow(match) > 0) rater_rel <- match[1, , drop = FALSE]
  }

  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_overall_2 <- extract_overall_pca_second(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  pca_reference_text <- build_pca_reference_text(warning_threshold_profiles()$pca_reference_bands)

  note_map <- list()
  note_map$table1 <- paste0(
    "Table 1. Facet summary\n",
    "Note. Measures are reported in logits; higher values indicate more of the modeled trait for that facet. ",
    "SE = standard error; MnSq = mean-square fit. ",
    "Model = ", model, "; estimation = ", method,
    "; N = ", fmt_count(n_obs), " observations from ", fmt_count(n_person),
    " persons on a ", fmt_count(n_cat), "-category scale (", fmt_count(rating_min), "-", fmt_count(rating_max), ")."
  )

  note_map$table2 <- paste0(
    "Table 2. Rating scale diagnostics\n",
    "Note. Category counts and thresholds summarize scale functioning. Thresholds were ", threshold_text,
    "; unused categories = ", fmt_count(unused), "; low-count categories (< 10) = ", fmt_count(low_count), "."
  )

  fit_sentence <- if (is.finite(infit) && is.finite(outfit)) {
    paste0("Overall fit: infit MnSq = ", fmt_num(infit), ", outfit MnSq = ", fmt_num(outfit), ".")
  } else {
    "Overall fit indices are reported as mean infit and outfit MnSq."
  }

  note_map$table3 <- paste0(
    "Table 3. Fit and reliability summary\n",
    "Note. Separation and reliability are based on observed variance and measurement error. ",
    fit_sentence
  )

  if (!is.null(rater_rel)) {
    note_map$table3 <- paste0(
      note_map$table3,
      " Rater facet (", rater_facet, ") reliability = ", fmt_num(rater_rel$Reliability),
      ", separation = ", fmt_num(rater_rel$Separation), "."
    )
  }

  if (!is.null(bias_results) && !is.null(bias_results$table) && nrow(bias_results$table) > 0) {
    sig <- bias_results$table |> dplyr::filter(as.numeric(`Prob.`) < 0.05)
    note_map$table4 <- paste0(
      "Table 4. Bias/interaction effects\n",
      "Note. Bias contrasts are in logits and represent observed minus expected scores with main effects held fixed. ",
      "Significant interactions (p < .05) = ", fmt_count(nrow(sig)), "."
    )
  } else {
    note_map$table4 <- paste0(
      "Table 4. Bias/interaction effects\n",
      "Note. Bias contrasts are in logits and represent observed minus expected scores with main effects held fixed."
    )
  }

  note_map$figure1 <- "Figure 1. Wright map\nNote. Persons and facet elements are located on a common logit scale; higher values indicate higher ability or greater severity/difficulty depending on facet orientation."
  note_map$figure2 <- "Figure 2. Pathway map (measure vs. infit t)\nNote. Points show element measures and their standardized infit values. Extreme |ZSTD| values flag potential misfit."
  note_map$figure3 <- "Figure 3. Facet estimate distribution\nNote. Distributions summarize severity/difficulty spread within each facet."
  note_map$figure4 <- "Figure 4. Step/threshold estimates\nNote. Step ordering should generally increase; disordered thresholds suggest category structure issues."
  note_map$figure5 <- "Figure 5. Category probability curves\nNote. Curves show the most probable category across the latent continuum; well-functioning categories show distinct peaks in order."
  note_map$figure6 <- "Figure 6. Observed vs expected scores\nNote. Points summarize mean observed and expected scores by bin; deviations from the diagonal suggest local misfit."
  note_map$figure7 <- "Figure 7. Fit diagnostics (Infit vs Outfit)\nNote. Each point represents an element within a facet. Values near 1.0 indicate expected fit; values substantially above 1.0 suggest misfit."
  note_map$figure8 <- "Figure 8. Fit ZSTD distribution\nNote. Distributions of standardized fit help identify unusually large residuals across facets."
  note_map$figure9 <- "Figure 9. Misfit levels\nNote. Levels are ranked by maximum |ZSTD| to highlight potentially problematic elements."

  if (!is.null(pca_overall_1)) {
    overall_tail <- paste0(
      " PC1 eigenvalue = ", fmt_num(pca_overall_1$Eigenvalue),
      " (", fmt_num(100 * to_float(pca_overall_1$Proportion), 1), "% variance)."
    )
    if (!is.null(pca_overall_2)) {
      overall_tail <- paste0(overall_tail, " PC2 eigenvalue = ", fmt_num(pca_overall_2$Eigenvalue), ".")
    }
    note_map$figure10 <- paste0(
      "Figure 10. Residual PCA scree (overall)\n",
      "Note. Eigenvalues are from PCA of the person x facet-combination standardized residual correlation matrix.",
      overall_tail,
      " ",
      pca_reference_text
    )
  } else {
    note_map$figure10 <- paste0(
      "Figure 10. Residual PCA scree (overall)\n",
      "Note. Overall residual PCA was not available for this run. ",
      pca_reference_text
    )
  }

  if (nrow(pca_facet_1) > 0) {
    top <- head(pca_facet_1, 3)
    labels <- vapply(seq_len(nrow(top)), function(i) {
      paste0(
        top$Facet[i], ": ", fmt_num(top$Eigenvalue[i]), " (",
        fmt_num(100 * to_float(top$Proportion[i]), 1), "%)"
      )
    }, character(1))

    note_map$figure11 <- paste0(
      "Figure 11. Residual PCA by facet\n",
      "Note. Each facet is analyzed using a person x facet-level standardized residual matrix. ",
      "Largest PC1 signals: ", paste(labels, collapse = "; "), ". ",
      pca_reference_text
    )
  } else {
    note_map$figure11 <- paste0(
      "Figure 11. Residual PCA by facet\n",
      "Note. Facet-specific residual PCA was not available for this run. ",
      pca_reference_text
    )
  }

  note_map
}

build_apa_table_figure_notes <- function(res, diagnostics, bias_results = NULL, context = list(), whexact = FALSE) {
  note_map <- build_apa_table_figure_note_map(
    res = res,
    diagnostics = diagnostics,
    bias_results = bias_results,
    context = context,
    whexact = whexact
  )

  ordered_keys <- c(
    "table1", "table2", "table3", "table4",
    "figure1", "figure2", "figure3", "figure4", "figure5",
    "figure6", "figure7", "figure8", "figure9", "figure10", "figure11"
  )
  paste(vapply(ordered_keys[ordered_keys %in% names(note_map)], function(k) note_map[[k]], character(1)), collapse = "\n\n")
}

build_apa_table_figure_captions <- function(res, diagnostics, bias_results = NULL, context = list()) {
  assessment <- trimws(as.character(context$assessment %||% ""))
  facet_pair <- if (!is.null(bias_results) && !is.null(bias_results$facet_a) && !is.null(bias_results$facet_b)) {
    paste0(bias_results$facet_a, " x ", bias_results$facet_b)
  } else {
    ""
  }

  assessment_phrase <- if (nzchar(assessment)) paste0(" for ", assessment) else ""

  blocks <- c(
    paste0("Table 1\nFacet Summary (Measures, SE, Fit, Reliability)", assessment_phrase),
    "Table 2\nRating Scale Diagnostics (Category Counts and Thresholds)",
    "Table 3\nFit and Reliability Summary",
    if (nzchar(facet_pair)) paste0("Table 4\nBias/Interaction Effects for ", facet_pair) else "Table 4\nBias/Interaction Effects",
    paste0("Figure 1\nWright Map of Person and Facet Measures", assessment_phrase),
    "Figure 2\nPathway Map (Measure vs. Infit t)",
    "Figure 3\nFacet Estimate Distribution",
    "Figure 4\nStep/Threshold Estimates",
    "Figure 5\nCategory Probability Curves",
    "Figure 6\nObserved vs. Expected Scores",
    "Figure 7\nFit Diagnostics (Infit vs Outfit)",
    "Figure 8\nFit ZSTD Distribution",
    "Figure 9\nMisfit Levels (Max |ZSTD|)",
    "Figure 10\nResidual PCA Scree (Overall)",
    "Figure 11\nResidual PCA by Facet"
  )

  paste(blocks, collapse = "\n\n")
}

build_visual_warning_map <- function(res,
                                     diagnostics,
                                     whexact = FALSE,
                                     thresholds = NULL,
                                     threshold_profile = "standard") {
  # Figure-level warning text is accumulated independently for figure1..figure11.
  warnings <- stats::setNames(replicate(11, character(0), simplify = FALSE), paste0("figure", 1:11))
  if (is.null(res) || is.null(diagnostics)) return(warnings)

  # Stage 1: Resolve active threshold profile (strict/standard/lenient + overrides).
  resolved <- resolve_warning_thresholds(thresholds = thresholds, threshold_profile = threshold_profile)
  active <- resolved$thresholds
  profile_name <- resolved$profile_name
  pca_reference_text <- build_pca_reference_text(resolved$pca_reference_bands)

  n_obs_min <- active$n_obs_min %||% 100
  n_person_min <- active$n_person_min %||% 30
  low_cat_min <- active$low_cat_min %||% 10
  min_facet_levels <- active$min_facet_levels %||% 3
  misfit_ratio_warn <- active$misfit_ratio_warn %||% 0.10
  missing_fit_ratio_warn <- active$missing_fit_ratio_warn %||% 0.20
  zstd2_ratio_warn <- active$zstd2_ratio_warn %||% 0.10
  zstd3_ratio_warn <- active$zstd3_ratio_warn %||% 0.05
  expected_var_min <- active$expected_var_min %||% 0.20
  pca_first_eigen_warn <- active$pca_first_eigen_warn %||% 2.0
  pca_first_prop_warn <- active$pca_first_prop_warn %||% 0.10

  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  n_obs <- if (!is.null(summary)) to_float(summary$N) else NA_real_
  n_person <- if (!is.null(res$facets$person)) nrow(res$facets$person) else 0

  # Stage 2: Sample-size and design warnings (Wright/pathway/observed-expected plots).
  if (is.finite(n_obs) && n_obs < n_obs_min) {
    warnings$figure1 <- c(warnings$figure1, paste0("Small number of observations (N = ", fmt_count(n_obs), " < ", fmt_count(n_obs_min), ")."))
    warnings$figure2 <- c(warnings$figure2, paste0("Small number of observations (N = ", fmt_count(n_obs), " < ", fmt_count(n_obs_min), "); ZSTD values may be volatile."))
    warnings$figure6 <- c(warnings$figure6, paste0("Small number of observations (N = ", fmt_count(n_obs), " < ", fmt_count(n_obs_min), "); bin averages may be noisy."))
  }

  if (n_person < n_person_min) {
    warnings$figure1 <- c(warnings$figure1, paste0("Small person sample (n = ", fmt_count(n_person), " < ", fmt_count(n_person_min), "); interpret spread cautiously."))
  }

  facet_levels <- res$config$facet_levels
  small_facets <- names(facet_levels)[vapply(facet_levels, length, integer(1)) < min_facet_levels]
  if (length(small_facets) > 0) {
    warnings$figure1 <- c(warnings$figure1, paste0("Facets with very few levels: ", paste(small_facets, collapse = ", "), "."))
    warnings$figure3 <- c(warnings$figure3, paste0("Facet distributions are based on few levels: ", paste(small_facets, collapse = ", "), "."))
  }

  # Stage 3: Rating-scale warnings.
  cat_tbl <- calc_category_stats(diagnostics$obs, res = res, whexact = whexact)
  if (nrow(cat_tbl) > 0) {
    unused <- sum(cat_tbl$Count == 0, na.rm = TRUE)
    low_count <- sum(cat_tbl$Count < low_cat_min, na.rm = TRUE)
    if (unused > 0) warnings$figure5 <- c(warnings$figure5, paste0("Unused categories detected (n = ", fmt_count(unused), ")."))
    if (low_count > 0) warnings$figure5 <- c(warnings$figure5, paste0("Low-count categories (< ", fmt_count(low_cat_min), ") detected (n = ", fmt_count(low_count), ")."))
  }

  step_order <- calc_step_order(res$steps)
  if (nrow(step_order) > 0) {
    disordered <- step_order |> dplyr::filter(Ordered == FALSE)
    if (nrow(disordered) > 0) {
      warnings$figure4 <- c(warnings$figure4, paste0("Disordered thresholds detected (n = ", fmt_count(nrow(disordered)), ")."))
      warnings$figure5 <- c(warnings$figure5, "Disordered thresholds can distort category curves.")
    }
  }

  # Stage 4: Fit-based warnings.
  measures <- diagnostics$measures
  if (is.null(measures) || nrow(measures) == 0) {
    warnings$figure2 <- c(warnings$figure2, "Fit statistics are not available for this run.")
    warnings$figure7 <- c(warnings$figure7, "Fit statistics are not available for this run.")
    warnings$figure8 <- c(warnings$figure8, "ZSTD distributions are not available for this run.")
    warnings$figure9 <- c(warnings$figure9, "Misfit ranking requires fit statistics.")
    return(warnings)
  }

  infit <- suppressWarnings(as.numeric(measures$Infit))
  outfit <- suppressWarnings(as.numeric(measures$Outfit))
  infit_z <- suppressWarnings(as.numeric(measures$InfitZSTD))
  outfit_z <- suppressWarnings(as.numeric(measures$OutfitZSTD))

  valid_fit <- is.finite(infit) & is.finite(outfit)
  if (length(valid_fit) > 0) {
    missing_ratio <- 1 - mean(valid_fit)
    if (is.finite(missing_ratio) && missing_ratio >= missing_fit_ratio_warn) {
      warnings$figure7 <- c(warnings$figure7, paste0("Fit statistics missing for ", sprintf("%.0f", missing_ratio * 100), "% of elements."))
    }
  }

  misfit <- (infit < 0.5) | (infit > 1.5) | (outfit < 0.5) | (outfit > 1.5)
  misfit_ratio <- mean(misfit, na.rm = TRUE)
  if (is.finite(misfit_ratio) && misfit_ratio > misfit_ratio_warn) {
    warnings$figure7 <- c(warnings$figure7, paste0("High proportion of misfit elements (", sprintf("%.0f", misfit_ratio * 100), "%)."))
  }

  zstd <- pmax(abs(infit_z), abs(outfit_z), na.rm = TRUE)
  zstd <- zstd[is.finite(zstd)]
  if (length(zstd) > 0) {
    prop2 <- mean(zstd >= 2)
    prop3 <- mean(zstd >= 3)
    if (prop2 > zstd2_ratio_warn) warnings$figure8 <- c(warnings$figure8, paste0("Large share of |ZSTD| >= 2 (", sprintf("%.0f", prop2 * 100), "%)."))
    if (prop3 > zstd3_ratio_warn) warnings$figure8 <- c(warnings$figure8, paste0("Notable |ZSTD| >= 3 (", sprintf("%.0f", prop3 * 100), "%)."))
  }

  obs <- diagnostics$obs
  if (!is.null(obs) && nrow(obs) > 0 && "Expected" %in% names(obs)) {
    exp_var <- stats::var(suppressWarnings(as.numeric(obs$Expected)), na.rm = TRUE)
    if (is.finite(exp_var) && exp_var < expected_var_min) {
      warnings$figure6 <- c(warnings$figure6, "Expected scores have limited spread; trends may be muted.")
    }
  }

  # Stage 5: Residual PCA warnings (overall and by facet).
  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  warnings$figure10 <- c(
    warnings$figure10,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )
  warnings$figure11 <- c(
    warnings$figure11,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )

  if (is.null(pca_overall_1)) {
    warnings$figure10 <- c(warnings$figure10, "Overall residual PCA is not available.")
  } else {
    warnings$figure10 <- c(
      warnings$figure10,
      build_pca_check_text(
        eigenvalue = pca_overall_1$Eigenvalue,
        proportion = pca_overall_1$Proportion,
        reference_bands = resolved$pca_reference_bands
      )
    )
    if (to_float(pca_overall_1$Eigenvalue) > pca_first_eigen_warn) {
      warnings$figure10 <- c(
        warnings$figure10,
        paste0("Overall residual PCA PC1 eigenvalue is high (", fmt_num(pca_overall_1$Eigenvalue), ").")
      )
    }
    if (to_float(pca_overall_1$Proportion) > pca_first_prop_warn) {
      warnings$figure10 <- c(
        warnings$figure10,
        paste0("Overall residual PCA PC1 explains ", fmt_num(100 * to_float(pca_overall_1$Proportion), 1), "% variance.")
      )
    }
  }

  if (nrow(pca_facet_1) == 0) {
    warnings$figure11 <- c(warnings$figure11, "Facet-specific residual PCA is not available.")
  } else {
    top <- pca_facet_1[1, , drop = FALSE]
    warnings$figure11 <- c(
      warnings$figure11,
      paste0(
        "Top facet PC1 (", as.character(top$Facet), "): ",
        build_pca_check_text(
          eigenvalue = top$Eigenvalue,
          proportion = top$Proportion,
          reference_bands = resolved$pca_reference_bands
        )
      )
    )
    flagged <- pca_facet_1[
      (suppressWarnings(as.numeric(pca_facet_1$Eigenvalue)) > pca_first_eigen_warn) |
        (suppressWarnings(as.numeric(pca_facet_1$Proportion)) > pca_first_prop_warn),
      ,
      drop = FALSE
    ]
    if (nrow(flagged) > 0) {
      warnings$figure11 <- c(
        warnings$figure11,
        paste0("Facet residual PCA shows stronger PC1 signal in: ", paste(flagged$Facet, collapse = ", "), ".")
      )
    }
  }

  warnings
}

build_visual_summary_map <- function(res,
                                     diagnostics,
                                     whexact = FALSE,
                                     options = NULL,
                                     thresholds = NULL,
                                     threshold_profile = "standard") {
  # Summary map mirrors the warning map but provides descriptive, non-binary text.
  options <- options %||% list()
  detail <- tolower(as.character(options$detail %||% "standard"))
  max_facet_ranges <- as.integer(options$max_facet_ranges %||% 4)
  top_misfit_n <- as.integer(options$top_misfit_n %||% 3)
  include_top_misfit <- top_misfit_n > 0
  resolved <- resolve_warning_thresholds(thresholds = thresholds, threshold_profile = threshold_profile)
  active <- resolved$thresholds
  profile_name <- resolved$profile_name
  pca_reference_text <- build_pca_reference_text(resolved$pca_reference_bands)
  pca_first_eigen_warn <- active$pca_first_eigen_warn %||% 2.0
  pca_first_prop_warn <- active$pca_first_prop_warn %||% 0.10

  summaries <- stats::setNames(replicate(11, character(0), simplify = FALSE), paste0("figure", 1:11))
  if (is.null(res) || is.null(diagnostics)) return(summaries)

  # Stage 1: Global design summary.
  summary <- if (!is.null(res$summary) && nrow(res$summary) > 0) res$summary[1, , drop = FALSE] else NULL
  n_obs <- if (!is.null(summary)) to_float(summary$N) else NA_real_
  n_person <- if (!is.null(res$facets$person)) nrow(res$facets$person) else 0

  if (is.finite(n_obs)) summaries$figure1 <- c(summaries$figure1, paste0("Observations: N = ", fmt_count(n_obs), "."))
  summaries$figure1 <- c(summaries$figure1, paste0("Persons: n = ", fmt_count(n_person), "."))

  person_stats <- describe_series(res$facets$person$Estimate)
  if (!is.null(person_stats)) {
    summaries$figure1 <- c(summaries$figure1,
      paste0("Person range ", fmt_num(person_stats$min), " to ", fmt_num(person_stats$max),
             " (M = ", fmt_num(person_stats$mean), ", SD = ", fmt_num(person_stats$sd), ").")
    )
    if (detail == "detailed" && is.finite(person_stats$sd)) {
      summaries$figure1 <- c(summaries$figure1, paste0("Person spread (SD) = ", fmt_num(person_stats$sd), " logits."))
    }
  }

  if (!is.null(res$facets$others) && nrow(res$facets$others) > 0) {
    facet_stats <- c()
    for (facet in unique(res$facets$others$Facet)) {
      df <- res$facets$others |> dplyr::filter(Facet == facet)
      st <- describe_series(df$Estimate)
      if (!is.null(st)) {
        facet_stats <- c(facet_stats, paste0(facet, ": n = ", fmt_count(nrow(df)), ", range ", fmt_num(st$min), " to ", fmt_num(st$max)))
      }
    }
    if (length(facet_stats) > 0) {
      summaries$figure1 <- c(summaries$figure1, paste0("Facet ranges: ", paste(head(facet_stats, max_facet_ranges), collapse = "; "), "."))
      if (length(facet_stats) > max_facet_ranges) summaries$figure1 <- c(summaries$figure1, "Additional facets omitted for brevity.")
    }
  }

  # Stage 2: Fit and category summaries.
  measures <- diagnostics$measures
  if (!is.null(measures) && nrow(measures) > 0) {
    infit_z <- suppressWarnings(as.numeric(measures$InfitZSTD))
    valid <- measures[is.finite(infit_z), , drop = FALSE]
    if (nrow(valid) > 0) {
      count2 <- sum(abs(infit_z) >= 2, na.rm = TRUE)
      summaries$figure2 <- c(summaries$figure2, paste0("Elements with |Infit ZSTD| >= 2: ", fmt_count(count2), " of ", fmt_count(nrow(valid)), "."))

      if (include_top_misfit) {
        tmp <- valid
        tmp$Abs <- abs(suppressWarnings(as.numeric(tmp$InfitZSTD)))
        top <- tmp |> dplyr::arrange(dplyr::desc(Abs)) |> dplyr::slice_head(n = top_misfit_n)
        if (nrow(top) > 0) {
          labels <- vapply(seq_len(nrow(top)), function(i) {
            paste0(top$Facet[i], ": ", truncate_label(top$Level[i], 20), " (|Z|=", fmt_num(top$Abs[i]), ")")
          }, character(1))
          summaries$figure2 <- c(summaries$figure2, paste0("Largest |Z|: ", paste(labels, collapse = "; "), "."))
        }
      }
    }
  }

  if (!is.null(res$facets$others) && nrow(res$facets$others) > 0) {
    summaries$figure3 <- c(summaries$figure3, "Distributions show the spread of severity/difficulty within each facet.")
  }

  step_tbl <- res$steps
  if (!is.null(step_tbl) && nrow(step_tbl) > 0) {
    step_order <- calc_step_order(step_tbl)
    disordered <- step_order |> dplyr::filter(Ordered == FALSE)
    summaries$figure4 <- c(summaries$figure4, paste0("Steps estimated: ", fmt_count(nrow(step_tbl)), "."))
    summaries$figure4 <- c(summaries$figure4, paste0("Disordered steps: ", fmt_count(nrow(disordered)), "."))
  }

  cat_tbl <- calc_category_stats(diagnostics$obs, res = res, whexact = whexact)
  if (nrow(cat_tbl) > 0) {
    used <- sum(cat_tbl$Count > 0, na.rm = TRUE)
    total <- nrow(cat_tbl)
    max_pct <- suppressWarnings(max(cat_tbl$Percent, na.rm = TRUE))
    summaries$figure5 <- c(summaries$figure5, paste0("Categories used: ", fmt_count(used), " of ", fmt_count(total), "."))
    if (is.finite(max_pct)) summaries$figure5 <- c(summaries$figure5, paste0("Largest category share: ", fmt_num(max_pct, 1), "%."))
  }

  obs <- diagnostics$obs
  if (!is.null(obs) && nrow(obs) > 0 && all(c("Observed", "Expected") %in% names(obs))) {
    resid <- suppressWarnings(as.numeric(obs$Observed) - as.numeric(obs$Expected))
    if ("Weight" %in% names(obs)) {
      w <- suppressWarnings(as.numeric(obs$Weight))
      w <- ifelse(is.finite(w) & w > 0, w, 0)
      mean_resid <- if (sum(w) > 0) sum(resid * w, na.rm = TRUE) / sum(w) else NA_real_
      mae <- if (sum(w) > 0) sum(abs(resid) * w, na.rm = TRUE) / sum(w) else NA_real_
    } else {
      mean_resid <- mean(resid, na.rm = TRUE)
      mae <- mean(abs(resid), na.rm = TRUE)
    }
    summaries$figure6 <- c(summaries$figure6, paste0("Mean residual: ", fmt_num(mean_resid), "."))
    summaries$figure6 <- c(summaries$figure6, paste0("Mean absolute residual: ", fmt_num(mae), "."))
  }

  if (!is.null(measures) && nrow(measures) > 0) {
    infit <- suppressWarnings(as.numeric(measures$Infit))
    outfit <- suppressWarnings(as.numeric(measures$Outfit))
    ok <- is.finite(infit) & is.finite(outfit)
    if (any(ok)) {
      misfit <- (infit < 0.5) | (infit > 1.5) | (outfit < 0.5) | (outfit > 1.5)
      summaries$figure7 <- c(summaries$figure7, paste0("Misfit elements (0.5-1.5 rule): ", fmt_count(sum(misfit, na.rm = TRUE)), " of ", fmt_count(sum(ok)), "."))
      if (detail == "detailed") {
        summaries$figure7 <- c(summaries$figure7, paste0("Mean infit = ", fmt_num(mean(infit, na.rm = TRUE)), ", mean outfit = ", fmt_num(mean(outfit, na.rm = TRUE)), "."))
      }
    }

    zstd <- pmax(abs(suppressWarnings(as.numeric(measures$InfitZSTD))), abs(suppressWarnings(as.numeric(measures$OutfitZSTD))), na.rm = TRUE)
    zstd_valid <- zstd[is.finite(zstd)]
    if (length(zstd_valid) > 0) {
      summaries$figure8 <- c(summaries$figure8, paste0("|ZSTD| >= 2: ", fmt_count(sum(zstd_valid >= 2)), "."))
      summaries$figure8 <- c(summaries$figure8, paste0("|ZSTD| >= 3: ", fmt_count(sum(zstd_valid >= 3)), "."))

      if (include_top_misfit) {
        tmp <- measures
        tmp$AbsZSTD <- zstd
        top <- tmp[is.finite(tmp$AbsZSTD), , drop = FALSE] |>
          dplyr::arrange(dplyr::desc(AbsZSTD)) |>
          dplyr::slice_head(n = top_misfit_n)
        if (nrow(top) > 0) {
          labels <- vapply(seq_len(nrow(top)), function(i) {
            paste0(top$Facet[i], ": ", truncate_label(top$Level[i], 20), " (|Z|=", fmt_num(top$AbsZSTD[i]), ")")
          }, character(1))
          summaries$figure9 <- c(summaries$figure9, paste0("Top misfit: ", paste(labels, collapse = "; "), "."))
        }
      }
    }
  }

  # Stage 3: Residual PCA summaries with threshold profile context.
  pca_obj <- safe_residual_pca(diagnostics, mode = "both")
  pca_overall_1 <- extract_overall_pca_first(pca_obj)
  pca_overall_2 <- extract_overall_pca_second(pca_obj)
  pca_facet_1 <- extract_facet_pca_first(pca_obj)
  summaries$figure10 <- c(
    summaries$figure10,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )
  summaries$figure11 <- c(
    summaries$figure11,
    paste0(
      "Threshold profile: ", profile_name,
      " (PC1 EV >= ", fmt_num(pca_first_eigen_warn, 1),
      ", variance >= ", fmt_num(100 * pca_first_prop_warn, 0), "%)."
    ),
    pca_reference_text
  )

  if (!is.null(pca_overall_1)) {
    summaries$figure10 <- c(
      summaries$figure10,
      paste0(
        "Overall residual PCA PC1: eigenvalue = ", fmt_num(pca_overall_1$Eigenvalue),
        ", variance = ", fmt_num(100 * to_float(pca_overall_1$Proportion), 1), "%."
      )
    )
    summaries$figure10 <- c(
      summaries$figure10,
      build_pca_check_text(
        eigenvalue = pca_overall_1$Eigenvalue,
        proportion = pca_overall_1$Proportion,
        reference_bands = resolved$pca_reference_bands
      )
    )
    if (!is.null(pca_overall_2)) {
      summaries$figure10 <- c(
        summaries$figure10,
        paste0("Overall residual PCA PC2: eigenvalue = ", fmt_num(pca_overall_2$Eigenvalue), ".")
      )
    }
  } else {
    summaries$figure10 <- c(summaries$figure10, "Overall residual PCA unavailable.")
  }

  if (nrow(pca_facet_1) > 0) {
    show_n <- min(nrow(pca_facet_1), ifelse(detail == "detailed", 5L, 3L))
    top <- pca_facet_1[seq_len(show_n), , drop = FALSE]
    summaries$figure11 <- c(
      summaries$figure11,
      paste0(
        "Top facet PC1 (", as.character(top$Facet[1]), "): ",
        build_pca_check_text(
          eigenvalue = top$Eigenvalue[1],
          proportion = top$Proportion[1],
          reference_bands = resolved$pca_reference_bands
        )
      )
    )
    labels <- vapply(seq_len(nrow(top)), function(i) {
      paste0(
        top$Facet[i], "=", fmt_num(top$Eigenvalue[i]),
        " (", fmt_num(100 * to_float(top$Proportion[i]), 1), "%)"
      )
    }, character(1))
    summaries$figure11 <- c(
      summaries$figure11,
      paste0("Facet residual PCA PC1 signals: ", paste(labels, collapse = "; "), ".")
    )
    if (nrow(pca_facet_1) > show_n) {
      summaries$figure11 <- c(summaries$figure11, "Additional facets omitted for brevity.")
    }
  } else {
    summaries$figure11 <- c(summaries$figure11, "Facet-specific residual PCA unavailable.")
  }

  summaries
}
