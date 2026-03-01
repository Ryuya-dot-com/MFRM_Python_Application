# mfrmr 0.1.0

## Initial release

- Native R implementation of many-facet Rasch model (MFRM) estimation without TAM/sirt backends.
- Supports arbitrary facet counts with `fit_mfrm()` and method selection (`MML` default, `JML`).
- Includes FACETS-style bias/interaction iterative estimation via `estimate_bias()`.
- Provides fixed-width report helpers (`build_fixed_reports()`).
- Adds APA-style narrative output (`build_apa_outputs()`).
- Adds visual warning summaries (`build_visual_summaries()`) with configurable threshold profiles.
- Implements residual PCA diagnostics and visualization (`analyze_residual_pca()`, `plot_residual_pca()`).
- Bundles Eckes & Jin (2021)-inspired synthetic Study 1/2 datasets in both `data/` and `inst/extdata/`.

## Package operations and publication readiness

- Added GitHub Actions CI for cross-platform `R CMD check`.
- Added `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, and `SECURITY.md`.
- Added citation metadata (`inst/CITATION`, `CITATION.cff`).
- Expanded README with explicit installation and citation instructions.
