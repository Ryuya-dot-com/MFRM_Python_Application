# Extracted from test-identifiability-constraints.R:195

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "mfrmr", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- mfrmr:::sample_mfrm_data(seed = 42)
fit <- suppressWarnings(fit_mfrm(
    d, "Person", c("Rater", "Task", "Criterion"), "Score",
    dummy_facets = "Criterion",
    method = "JML", model = "RSM", maxit = 40, quad_points = 7
  ))
