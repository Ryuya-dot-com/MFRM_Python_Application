# Code Reading Guide (for Human Developers)

This guide is written for intermediate R users who want to understand and modify `mfrmr`.

## 1) Start here: public API

Read `/R/api.R` first.

- `fit_mfrm()` is the entry point for estimation.
- `diagnose_mfrm()` builds diagnostics from a fitted model.
- `analyze_residual_pca()` / `plot_residual_pca()` handle residual PCA.
- `estimate_bias()` / `build_fixed_reports()` produce FACETS-style interaction reports.
- `build_apa_outputs()` / `build_visual_summaries()` produce narrative output.

If you only need to use the package, you can stop at `api.R` + the vignettes.

## 2) Core estimation flow

Core numeric logic is in `/R/mfrm_core.R`.

`mfrm_estimate()` is intentionally structured in five stages:

1. normalize inputs and model options
2. resolve facet/model configuration
3. build optimization structures
4. run optimization (JMLE or MML)
5. convert optimized parameters to readable tables

Helper functions around `mfrm_estimate()` are separated so each stage can be debugged independently.

## 3) Data assumptions

Input data is **long format** (one row = one rating event).

Required columns at runtime:

- person id (`person` argument)
- one or more facet columns (`facets` argument)
- observed score (`score` argument)
- optional frequency weight (`weight` argument)

Packaged example datasets can be loaded via `data()` or `load_mfrmr_data()`.

## 4) Diagnostics flow

Diagnostics are computed in `mfrm_diagnostics()` (`/R/mfrm_core.R`).

Key components:

- observed/expected tables
- fit statistics (Infit/Outfit and ZSTD)
- reliability/separation
- category and threshold diagnostics
- optional residual PCA (overall and by facet)

## 5) Reporting flow

Text/report builders live in `/R/reporting.R`.

- APA builders combine model summary + diagnostics into prose and figure/table notes.
- visual summary/warning builders produce figure-specific text for UIs.
- threshold profiles (`strict/standard/lenient`) are resolved in one place to avoid hidden constants.

## 6) Sign conventions

Facet orientation is controlled by `positive_facets`.

- facets listed in `positive_facets` use `+1`
- other facets use `-1`

This sign map is used consistently when converting between logits, fair averages, and report values.

## 7) Where to edit safely

If you need to change behavior:

- estimation algorithm: edit `mfrm_core.R` around log-likelihood/optimization helpers
- output wording or thresholds: edit `reporting.R`
- user-facing function signatures: edit `api.R` and regenerate docs

After changes, always run:

1. `Rscript -e 'testthat::test_dir("mfrmr/tests/testthat")'`
2. `R CMD build mfrmr`
3. `R CMD check --no-manual mfrmr_0.1.0.tar.gz`
