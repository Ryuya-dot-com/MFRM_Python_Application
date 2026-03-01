# --------------------------------------------------------------------------
# test-export-mfrm.R
# Tests for export_mfrm(), as.data.frame.mfrm_fit()
# --------------------------------------------------------------------------

# ---- Helper fit object ---------------------------------------------------

local_fit <- function(envir = parent.frame()) {
  d <- mfrmr:::sample_mfrm_data(seed = 42)
  fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
  fit
}

# ---- 1. as.data.frame.mfrm_fit returns correct structure ----------------

test_that("as.data.frame.mfrm_fit returns data.frame with expected columns", {
  fit <- local_fit()
  df <- as.data.frame(fit)
  expect_s3_class(df, "data.frame")
  expect_true(all(c("Facet", "Level", "Estimate") %in% names(df)))
  expect_true(nrow(df) > 0)
})

test_that("as.data.frame.mfrm_fit includes both Person and other facets", {
  fit <- local_fit()
  df <- as.data.frame(fit)
  facets_in <- unique(df$Facet)
  expect_true("Person" %in% facets_in)
  expect_true(any(facets_in != "Person"))
})

test_that("as.data.frame.mfrm_fit Level column is character", {
  fit <- local_fit()
  df <- as.data.frame(fit)
  expect_type(df$Level, "character")
})

# ---- 2. export_mfrm writes CSV files ------------------------------------

test_that("export_mfrm writes person and facet CSVs", {
  fit <- local_fit()
  tmpdir <- tempfile("export_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  result <- export_mfrm(fit, output_dir = tmpdir, prefix = "test",
                         tables = c("person", "facets"))
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(file.exists(result$Path)))

  person_csv <- utils::read.csv(result$Path[result$Table == "person"])
  expect_true(nrow(person_csv) > 0)
  expect_true("Estimate" %in% names(person_csv))
})

test_that("export_mfrm writes summary CSV", {
  fit <- local_fit()
  tmpdir <- tempfile("export_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  result <- export_mfrm(fit, output_dir = tmpdir, prefix = "test",
                         tables = "summary")
  expect_equal(nrow(result), 1)
  expect_true(file.exists(result$Path[1]))

  summary_csv <- utils::read.csv(result$Path[1])
  expect_true(nrow(summary_csv) >= 1)
})

test_that("export_mfrm writes step CSV", {
  fit <- local_fit()
  tmpdir <- tempfile("export_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  result <- export_mfrm(fit, output_dir = tmpdir, prefix = "test",
                         tables = "steps")
  # Steps may or may not exist depending on model
  # but the function should not error
  expect_s3_class(result, "data.frame")
})

# ---- 3. export_mfrm with diagnostics ------------------------------------

test_that("export_mfrm enriches facets when diagnostics provided", {
  fit <- local_fit()
  diag <- suppressWarnings(diagnose_mfrm(fit, residual_pca = "none"))
  tmpdir <- tempfile("export_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  result <- export_mfrm(fit, diagnostics = diag, output_dir = tmpdir,
                         prefix = "enriched", tables = c("facets", "measures"))
  expect_true(nrow(result) >= 2)

  facet_csv <- utils::read.csv(
    result$Path[result$Table == "facets"]
  )
  # Should have enrichment columns when diagnostics provided
  enrichment_cols <- intersect(c("SE", "Infit", "Outfit"), names(facet_csv))
  expect_true(length(enrichment_cols) > 0)
})

# ---- 4. overwrite protection --------------------------------------------

test_that("export_mfrm refuses to overwrite by default", {
  fit <- local_fit()
  tmpdir <- tempfile("export_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  export_mfrm(fit, output_dir = tmpdir, prefix = "ow", tables = "person")
  expect_error(
    export_mfrm(fit, output_dir = tmpdir, prefix = "ow", tables = "person"),
    "already exists"
  )
})

test_that("export_mfrm overwrites when overwrite = TRUE", {
  fit <- local_fit()
  tmpdir <- tempfile("export_test_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  export_mfrm(fit, output_dir = tmpdir, prefix = "ow", tables = "person")
  result <- export_mfrm(fit, output_dir = tmpdir, prefix = "ow",
                          tables = "person", overwrite = TRUE)
  expect_equal(nrow(result), 1)
})

# ---- 5. Input validation ------------------------------------------------

test_that("export_mfrm rejects non-mfrm_fit input", {
  expect_error(export_mfrm("not_a_fit"), "mfrm_fit")
})

test_that("export_mfrm rejects unknown table names", {
  fit <- local_fit()
  expect_error(export_mfrm(fit, tables = "invalid_table"), "Unknown table")
})
