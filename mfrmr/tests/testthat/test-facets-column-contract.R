contract_path <- function() {
  installed <- system.file("references", "facets_column_contract.csv", package = "mfrmr")
  if (nzchar(installed) && file.exists(installed)) {
    return(installed)
  }
  source_path <- testthat::test_path("..", "..", "inst", "references", "facets_column_contract.csv")
  if (file.exists(source_path)) {
    return(source_path)
  }
  source_path
}

split_required_columns <- function(x) {
  parts <- strsplit(as.character(x), "|", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

column_token_present <- function(token, columns) {
  token <- as.character(token)
  if (!nzchar(token)) return(TRUE)
  if (endsWith(token, "*")) {
    prefix <- substr(token, 1L, nchar(token) - 1L)
    return(any(startsWith(columns, prefix)))
  }
  token %in% columns
}

test_that("FACETS column contract file is available and valid", {
  path <- contract_path()
  expect_true(file.exists(path))

  contract <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_true(is.data.frame(contract))
  expect_true(nrow(contract) > 0)
  expect_true(all(c("table_id", "function_name", "object_id", "component", "required_columns") %in% names(contract)))
  expect_true(all(nzchar(contract$required_columns)))
})

test_that("FACETS column contract is satisfied by current outputs", {
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

  outputs <- list(
    t1 = mfrmr::specifications_report(fit),
    t2 = mfrmr::data_quality_report(
      fit,
      data = d,
      person = "Person",
      facets = c("Rater", "Task", "Criterion"),
      score = "Score"
    ),
    t3 = mfrmr::estimation_iteration_report(fit, max_iter = 5),
    t4 = mfrmr::unexpected_response_table(fit, diagnostics = diag, top_n = 20),
    t5 = mfrmr::measurable_summary_table(fit, diagnostics = diag),
    t6 = mfrmr::subset_connectivity_report(fit, diagnostics = diag),
    t62 = mfrmr::facet_statistics_report(fit, diagnostics = diag),
    t7chisq = mfrmr::facets_chisq_table(fit, diagnostics = diag),
    t7agree = mfrmr::interrater_agreement_table(fit, diagnostics = diag),
    t81 = mfrmr::rating_scale_table(fit, diagnostics = diag),
    t8bar = mfrmr::category_structure_report(fit, diagnostics = diag),
    t8curves = mfrmr::category_curves_report(fit, theta_points = 101),
    out = mfrmr::facets_output_file_bundle(fit, diagnostics = diag, include = c("graph", "score"), theta_points = 81),
    t14 = mfrmr::build_fixed_reports(bias, branch = "facets"),
    t10 = mfrmr::unexpected_after_bias_table(fit, bias, diagnostics = diag, top_n = 20),
    t11 = mfrmr::bias_count_table(bias, branch = "facets"),
    t12 = mfrmr::fair_average_table(fit, diagnostics = diag),
    t13 = mfrmr::bias_interaction_report(bias)
  )

  contract <- utils::read.csv(contract_path(), stringsAsFactors = FALSE)

  rows <- vector("list", nrow(contract))
  for (i in seq_len(nrow(contract))) {
    row <- contract[i, , drop = FALSE]
    obj <- outputs[[row$object_id]]
    expect_false(is.null(obj), info = paste("Unknown object_id:", row$object_id))

    comp <- obj[[row$component]]
    expect_true(is.data.frame(comp), info = paste("Component is not a data.frame:", row$object_id, "$", row$component, sep = ""))
    cols <- names(comp)

    tokens <- split_required_columns(row$required_columns)
    present <- vapply(tokens, column_token_present, logical(1), columns = cols)
    missing <- tokens[!present]

    rows[[i]] <- data.frame(
      table_id = row$table_id,
      function_name = row$function_name,
      object_id = row$object_id,
      component = row$component,
      required_n = length(tokens),
      present_n = sum(present),
      coverage = if (length(tokens) > 0) sum(present) / length(tokens) else 1,
      missing = paste(missing, collapse = " | "),
      stringsAsFactors = FALSE
    )
  }

  audit <- dplyr::bind_rows(rows)

  expect_true(all(audit$coverage == 1), info = paste(
    apply(audit[audit$coverage < 1, c("table_id", "component", "missing")], 1, paste, collapse = " :: "),
    collapse = "\n"
  ))

  parity_summary <- data.frame(
    Components = nrow(audit),
    FullCoverage = sum(audit$coverage == 1),
    MeanCoverage = mean(audit$coverage),
    MinCoverage = min(audit$coverage),
    stringsAsFactors = FALSE
  )
  expect_equal(parity_summary$FullCoverage, parity_summary$Components)
  expect_equal(parity_summary$MeanCoverage, 1)
  expect_equal(parity_summary$MinCoverage, 1)
})
