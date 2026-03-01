# test-error-handling.R
# Tests that input validation produces clear, actionable error messages.
# Each test targets the release code directly (no mocks).

# ---- fit_mfrm input validation ----

test_that("fit_mfrm rejects non-data.frame input", {
  expect_error(
    fit_mfrm(list(a = 1), "P", "R", "S"),
    "data\\.frame.*Convert"
  )
  expect_error(
    fit_mfrm("not a frame", "P", "R", "S"),
    "data\\.frame"
  )
})

test_that("fit_mfrm rejects empty data.frame", {
  expect_error(
    fit_mfrm(data.frame(), "P", "R", "S"),
    "zero rows"
  )
})

test_that("fit_mfrm rejects non-character person argument", {
  toy <- data.frame(P = 1:3, R = 1:3, S = 1:3)
  expect_error(
    fit_mfrm(toy, 1, "R", "S"),
    "person.*character string"
  )
  expect_error(
    fit_mfrm(toy, c("P", "Q"), "R", "S"),
    "person.*single"
  )
  expect_error(
    fit_mfrm(toy, "", "R", "S"),
    "person.*non-empty"
  )
})

test_that("fit_mfrm rejects invalid facets argument", {
  toy <- data.frame(P = 1:3, R = 1:3, S = 1:3)
  expect_error(
    fit_mfrm(toy, "P", character(0), "S"),
    "facets.*character vector"
  )
  expect_error(
    fit_mfrm(toy, "P", 1, "S"),
    "facets.*character vector"
  )
})

test_that("fit_mfrm rejects invalid score argument", {
  toy <- data.frame(P = 1:3, R = 1:3, S = 1:3)
  expect_error(
    fit_mfrm(toy, "P", "R", 1),
    "score.*character string"
  )
  expect_error(
    fit_mfrm(toy, "P", "R", ""),
    "score.*non-empty"
  )
})

test_that("fit_mfrm rejects invalid weight argument", {
  toy <- data.frame(P = 1:3, R = 1:3, S = 1:3)
  expect_error(
    fit_mfrm(toy, "P", "R", "S", weight = 1),
    "weight.*character string"
  )
})

test_that("fit_mfrm rejects invalid numeric parameters", {
  toy <- data.frame(P = 1:3, R = 1:3, S = 1:3)
  expect_error(fit_mfrm(toy, "P", "R", "S", maxit = -1), "maxit.*positive")
  expect_error(fit_mfrm(toy, "P", "R", "S", reltol = -0.1), "reltol.*positive")
  expect_error(fit_mfrm(toy, "P", "R", "S", quad_points = 0), "quad_points.*positive")
})

test_that("fit_mfrm rejects missing columns in data", {
  toy <- data.frame(A = 1:5, B = 1:5, C = 1:5)
  expect_error(
    fit_mfrm(toy, "Person", "Rater", "Score"),
    "not found in data.*Person"
  )
})

test_that("fit_mfrm rejects duplicate column references", {
  toy <- data.frame(A = 1:5, B = 1:5, C = 1:5)
  expect_error(
    fit_mfrm(toy, "A", "A", "C"),
    "distinct.*duplicates"
  )
})

# ---- diagnose_mfrm input validation ----

test_that("diagnose_mfrm rejects non-mfrm_fit input", {
  expect_error(
    diagnose_mfrm(list()),
    "mfrm_fit"
  )
  expect_error(
    diagnose_mfrm("not a fit"),
    "mfrm_fit"
  )
})

# ---- estimate_bias input validation ----

test_that("estimate_bias rejects invalid fit object", {
  expect_error(
    estimate_bias(list(), list(obs = data.frame()), facet_a = "R", facet_b = "C"),
    "mfrm_fit"
  )
})

test_that("estimate_bias rejects invalid diagnostics object", {
  fake_fit <- structure(list(), class = c("mfrm_fit", "list"))
  expect_error(
    estimate_bias(fake_fit, list(), facet_a = "R", facet_b = "C"),
    "diagnose_mfrm"
  )
})

test_that("estimate_bias requires interaction facets specification", {
  fake_fit <- structure(list(), class = c("mfrm_fit", "list"))
  fake_diag <- list(obs = data.frame(x = 1))
  expect_error(
    estimate_bias(fake_fit, fake_diag),
    "interaction_facets.*facet_a.*facet_b"
  )
})

# ---- prepare_mfrm_data edge cases ----

test_that("prepare_mfrm_data rejects all-NA data after filtering", {
  toy <- data.frame(Person = c("P1", "P2"), Rater = c("R1", "R2"), Score = c(NA, NA))
  expect_error(
    fit_mfrm(toy, "Person", "Rater", "Score"),
    "No valid observations"
  )
})

test_that("prepare_mfrm_data rejects single-category scores", {
  toy <- data.frame(
    Person = c("P1", "P2", "P3", "P4"),
    Rater = c("R1", "R2", "R1", "R2"),
    Score = c(1, 1, 1, 1)
  )
  expect_error(
    fit_mfrm(toy, "Person", "Rater", "Score"),
    "Only one score category"
  )
})

test_that("prepare_mfrm_data rejects zero-weight-only data", {
  toy <- data.frame(
    Person = c("P1", "P2"),
    Rater = c("R1", "R2"),
    Score = c(0, 1),
    W = c(0, 0)
  )
  expect_error(
    fit_mfrm(toy, "Person", "Rater", "Score", weight = "W"),
    "No valid observations"
  )
})

# ---- FACETS mode API validation ----

test_that("run_mfrm_facets rejects non-data.frame", {
  expect_error(
    run_mfrm_facets(list(a = 1)),
    "data\\.frame"
  )
})

test_that("run_mfrm_facets rejects empty data", {
  expect_error(
    run_mfrm_facets(data.frame()),
    "empty"
  )
})

test_that("infer_facets_mode_mapping rejects too few columns", {
  expect_error(
    run_mfrm_facets(data.frame(A = 1, B = 2)),
    "at least 3 columns"
  )
})

test_that("infer_facets_mode_mapping rejects missing person column", {
  toy <- data.frame(A = 1:3, B = 1:3, C = 1:3)
  expect_error(
    run_mfrm_facets(toy, person = "NonExistent"),
    "Person column not found"
  )
})

test_that("infer_facets_mode_mapping rejects missing score column", {
  toy <- data.frame(A = 1:3, B = 1:3, C = 1:3)
  expect_error(
    run_mfrm_facets(toy, score = "NonExistent"),
    "Score column not found"
  )
})

test_that("infer_facets_mode_mapping rejects missing facet columns", {
  toy <- data.frame(Person = 1:3, Score = 1:3, R = 1:3)
  expect_error(
    run_mfrm_facets(toy, person = "Person", facets = c("Missing1", "Missing2"), score = "Score"),
    "Facet column.*not found"
  )
})

test_that("normalize_spec_input rejects non-data.frame anchors", {
  toy <- mfrmr:::sample_mfrm_data(seed = 1)
  expect_error(
    run_mfrm_facets(toy, person = "Person",
                    facets = c("Rater", "Task", "Criterion"),
                    score = "Score", anchors = "bad_input", maxit = 5),
    "data\\.frame"
  )
})

# ---- Gauss-Hermite validation ----

test_that("gauss_hermite_normal rejects n < 1", {
  expect_error(
    mfrmr:::gauss_hermite_normal(0),
    "quadrature.*n >= 1"
  )
})

# ---- PCM step_facet validation ----

test_that("resolve_pcm_step_facet rejects invalid step facet", {
  expect_error(
    mfrmr:::resolve_pcm_step_facet("PCM", "Invalid", c("Rater", "Task")),
    "step_facet.*Invalid.*not among.*Rater.*Task"
  )
})
