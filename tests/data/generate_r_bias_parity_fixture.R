# Regenerate tests/data/r_bias_parity_output.json from the deterministic
# CSV at tests/data/r_bias_parity_input.csv via mfrmr 0.2.0.
#
# Usage from the repository root:
#   Rscript tests/data/generate_r_bias_parity_fixture.R
#
# The output JSON is the reference the Python parity test
# (test_bias_estimation_matches_r_reference_within_tolerance) compares
# against. Run this after either side of the math contract changes
# (the bias inference code, the mfrmr R package version, or the
# deterministic CSV input) to refresh the reference.

suppressMessages({
  library(mfrmr)
  library(jsonlite)
})

input_path <- "tests/data/r_bias_parity_input.csv"
output_path <- "tests/data/r_bias_parity_output.json"

stopifnot(file.exists(input_path))

df <- read.csv(input_path, stringsAsFactors = FALSE)
cat("read", nrow(df), "rows from", input_path, "\n")

fit <- fit_mfrm(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  rating_min = 0L,
  rating_max = 2L,
  model = "GPCM",
  step_facet = "Criterion",
  method = "MML",
  quad_points = 5L,
  maxit = 25L,
  reltol = 1e-4
)
cat("fit converged:", fit$summary$Converged[1], "; log_lik:",
    sprintf("%.6f", fit$summary$LogLik[1]), "\n")

dx <- diagnose_mfrm(fit, residual_pca = "none", diagnostic_mode = "legacy")
bias <- estimate_bias(fit, dx, facet_a = "Rater", facet_b = "Criterion")
tbl <- bias$table

# Project the columns the Python parity test asserts against.
keep <- c(
  "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
  "Bias Size", "S.E.", "t", "d.f.", "Prob.",
  "LR ChiSq", "LR d.f.", "LR Prob.",
  "Profile CI Lower", "Profile CI Upper",
  "Profile CI Level", "Profile CI Status",
  "Likelihood Basis", "InferenceTier",
  "Infit", "Outfit", "ObsN"
)
keep <- keep[keep %in% colnames(tbl)]
out <- tbl[, keep]

fixture <- list(
  metadata = list(
    r_version = R.version.string,
    mfrmr_version = as.character(packageVersion("mfrmr")),
    fit_summary = list(
      model = "GPCM",
      method = "MML",
      step_facet = "Criterion",
      quad_points = 5L,
      maxit = 25L,
      reltol = 1e-4,
      converged = as.logical(fit$summary$Converged[1]),
      log_lik = as.numeric(fit$summary$LogLik[1])
    ),
    n_cells = nrow(out)
  ),
  cells = out
)
write_json(fixture, output_path,
           auto_unbox = TRUE, na = "null", digits = 17, pretty = TRUE)
cat("wrote", nrow(out), "cells to", output_path, "\n")
cat("R version:", R.version.string, "\n")
cat("mfrmr version:", as.character(packageVersion("mfrmr")), "\n")
