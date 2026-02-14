# mfrmr

Native R package for flexible many-facet Rasch model (MFRM) estimation without TAM/sirt backends.

## Current scope

- Flexible facet count (`facets = c(...)`)
- Estimation methods: `MML` (default) and `JML` (`JMLE` internally)
- Models: `RSM`, `PCM`
- Bias/interaction iterative estimation (FACETS-style)
- Fixed-width text reports (FACETS-style)
- APA-style narrative output helpers
- Visual warning and summary maps
- Residual PCA for unidimensionality checks (`overall` / `facet` / `both`)

## Quick start

```r
library(mfrmr)

data("ej2021_study1", package = "mfrmr")
df <- ej2021_study1
fit <- fit_mfrm(
  data = df,
  person = "Person",
  facets = c("Rater", "Criterion"),
  score = "Score",
  method = "MML",   # default
  model = "RSM"
)

diag <- diagnose_mfrm(fit)
pca <- analyze_residual_pca(diag, mode = "both")
p_scree <- plot_residual_pca(pca, mode = "overall", plot_type = "scree")
p_load <- plot_residual_pca(pca, mode = "facet", facet = "Rater", plot_type = "loadings")
bias <- estimate_bias(fit, diag, facet_a = "Rater", facet_b = "Criterion")
fixed <- build_fixed_reports(bias)
apa <- build_apa_outputs(fit, diag)
warn <- build_visual_summaries(fit, diag)
warn_strict <- build_visual_summaries(fit, diag, threshold_profile = "strict")
profiles <- mfrm_threshold_profiles()
```

`build_apa_outputs()` and `build_visual_summaries()` automatically include
residual PCA narratives when diagnostics contain residual PCA results.
Residual PCA summaries now show literature-based multi-threshold bands
(eigenvalue and explained variance) with configurable profiles.

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
