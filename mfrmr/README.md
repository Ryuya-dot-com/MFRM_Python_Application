# mfrmr

[![GitHub](https://img.shields.io/badge/GitHub-mfrmr-181717?logo=github)](https://github.com/Ryuya-dot-com/mfrmr)
[![R-CMD-check](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/pkgdown.yaml)
[![test-coverage](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Ryuya-dot-com/mfrmr/actions/workflows/test-coverage.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Native R package for flexible many-facet Rasch model (MFRM) estimation without TAM/sirt backends.

## Current scope

- Flexible facet count (`facets = c(...)`)
- Estimation methods: `MML` (default) and `JML` (`JMLE` internally)
- Models: `RSM`, `PCM`
- FACETS-style one-shot wrapper (`run_mfrm_facets()`, alias `mfrmRFacets()`)
- Bias/interaction iterative estimation (FACETS-style)
- Optional fixed-width text reports for console/log audits (FACETS-style)
- APA-style narrative output helpers (`build_apa_outputs()`)
- Visual warning and summary maps (`build_visual_summaries()`)
- Residual PCA for unidimensionality checks (`overall` / `facet` / `both`)
- TAM-style descriptive data snapshot (`describe_mfrm_data()`)
- Anchor audit / normalization helper (`audit_mfrm_anchors()`)
- Anchor export helper for linking workflows (`make_anchor_table()`)

## Installation

```r
# GitHub (development version)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("Ryuya-dot-com/mfrmr")

# CRAN (after release)
# install.packages("mfrmr")
```

## Core workflow (recommended)

1. Fit model: `fit_mfrm()`
2. Diagnostics: `diagnose_mfrm()`
3. Optional residual PCA: `analyze_residual_pca()`
4. Optional interaction bias: `estimate_bias()`
5. Reporting: `apa_table()`, `build_apa_outputs()`, `build_visual_summaries()`
6. Optional FACETS-style parity audit: `facets_parity_report()`
7. Reproducible inspection: `summary()` and `plot(..., draw = FALSE)`

## Help-page navigation

All exported help pages are now aligned with two practical sections:

- `Interpreting output`
- `Typical workflow`

Recommended entry points:

- `?mfrmr-package` (package overview)
- `?mfrmr_workflow_methods` (method map for `summary()` / `plot()`)
- `?fit_mfrm`, `?diagnose_mfrm`, `?run_mfrm_facets`
- `?build_apa_outputs`, `?build_visual_summaries`, `?apa_table`

## Quick start (stepwise)

```r
library(mfrmr)

data("ej2021_study1", package = "mfrmr")
df <- ej2021_study1
desc <- describe_mfrm_data(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score"
)
aud <- audit_mfrm_anchors(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  min_common_anchors = 5,
  min_obs_per_element = 30,
  min_obs_per_category = 10
)
fit <- fit_mfrm(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",   # default
  model = "RSM",
  anchor_policy = "warn"
)
fit_s <- summary(fit)
fit_p <- plot(fit, draw = FALSE)

diag <- diagnose_mfrm(fit, residual_pca = "both")
diag_s <- summary(diag)
pca <- analyze_residual_pca(diag, mode = "both")
p_scree <- plot_residual_pca(pca, mode = "overall", plot_type = "scree", draw = FALSE)
p_load <- plot_residual_pca(pca, mode = "facet", facet = "Rater", plot_type = "loadings", draw = FALSE)

bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion")
bias_s <- summary(bias)
bias_p <- plot_bias_interaction(bias, draw = FALSE)

fixed <- build_fixed_reports(bias)
fixed_s <- summary(fixed)

apa <- build_apa_outputs(
  fit,
  diag,
  context = list(line_width = 92)
)
apa_s <- summary(apa)
cat(apa$report_text)

warn <- build_visual_summaries(fit, diag)
warn_s <- summary(warn)
warn_p <- plot(warn, type = "comparison", draw = FALSE)

warn_strict <- build_visual_summaries(fit, diag, threshold_profile = "strict")
profiles <- mfrm_threshold_profiles()
profiles_s <- summary(profiles)

tbl <- apa_table(fit, which = "summary", branch = "facets")
tbl_s <- summary(tbl)
tbl_p <- plot(tbl, draw = FALSE)

parity <- facets_parity_report(fit, diagnostics = diag, branch = "facets")
parity_s <- summary(parity)

anchors_next_run <- make_anchor_table(fit)
```

`build_apa_outputs()` and `build_visual_summaries()` automatically include
residual PCA narratives when diagnostics contain residual PCA results.
Residual PCA summaries now show literature-based multi-threshold bands
(eigenvalue and explained variance) with configurable profiles.

`report_text` (default) now includes convergence metrics, anchor/constraint
summary, threshold ordering details, top misfit levels, and bias-screen counts.

`build_apa_outputs()` returns an `mfrm_apa_outputs` object and supports
`summary(apa)` for compact coverage checks.

`mfrm_threshold_profiles()` returns an `mfrm_threshold_profiles` object and supports
`summary(profiles)` for side-by-side threshold comparison.

## FACETS-style one-shot wrapper

```r
run <- run_mfrm_facets(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "JML",
  model = "RSM"
)
run_s <- summary(run)
run_fit_plot <- plot(run, type = "fit", draw = FALSE)
run_qc_plot <- plot(run, type = "qc", draw = FALSE)

# Alias
run2 <- mfrmRFacets(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score"
)
```

## Public API map

Model and diagnostics:

- `fit_mfrm()`, `run_mfrm_facets()`, `mfrmRFacets()`
- `diagnose_mfrm()`, `analyze_residual_pca()`
- `estimate_bias()`, `bias_count_table()`

Table/report outputs:

- `specifications_report()`, `data_quality_report()`, `estimation_iteration_report()`
- `subset_connectivity_report()`, `facet_statistics_report()`
- `measurable_summary_table()`, `rating_scale_table()`
- `category_structure_report()`, `category_curves_report()`
- `unexpected_response_table()`, `unexpected_after_bias_table()`
- `fair_average_table()`, `displacement_table()`
- `interrater_agreement_table()`, `facets_chisq_table()`
- `facets_output_file_bundle()`, `facets_parity_report()`
- `bias_interaction_report()`, `build_fixed_reports()`
- `apa_table()`, `build_apa_outputs()`, `build_visual_summaries()`

Plots and QA dashboards:

- `plot_unexpected()`, `plot_fair_average()`, `plot_displacement()`
- `plot_interrater_agreement()`, `plot_facets_chisq()`
- `plot_bias_interaction()`, `plot_residual_pca()`, `plot_qc_dashboard()`

Data and anchor utilities:

- `describe_mfrm_data()`, `audit_mfrm_anchors()`, `make_anchor_table()`
- `mfrm_threshold_profiles()`, `list_mfrmr_data()`, `load_mfrmr_data()`

Legacy FACETS-style numbered names are internal and not exported.

## FACETS reference mapping

See:

- `inst/references/FACETS_manual_mapping.md`
- `inst/references/CODE_READING_GUIDE.md` (for developers/readers)

## Packaged synthetic datasets

Installed at `system.file("extdata", package = "mfrmr")`:

- `eckes_jin_2021_study1_sim.csv`
- `eckes_jin_2021_study2_sim.csv`
- `eckes_jin_2021_combined_sim.csv`
- `eckes_jin_2021_study1_itercal_sim.csv`
- `eckes_jin_2021_study2_itercal_sim.csv`
- `eckes_jin_2021_combined_itercal_sim.csv`

The same datasets are also packaged in `data/` and can be loaded with:

```r
data("ej2021_study1", package = "mfrmr")
# or
df <- load_mfrmr_data("study1")
```

## Citation

```r
citation("mfrmr")
```
