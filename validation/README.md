# Validation Notes

This directory is for optional validation artifacts. The Streamlit app does not
call R packages or external MFRM engines at runtime.

## Validation Stance

The goal is implementation credibility, not unconditional numerical identity
with another package. FACETS, TAM, sirt, mirt, and this app can differ in
constraints, parameterization, quadrature, optimizer details, latent-variance
treatment, and omitted likelihood constants.

Safe claim:

- The app targets functional parity for key MFRM workflows.

Unsafe claim without an archived parity report:

- The app is an exact replacement for FACETS, TAM, sirt, mirt, or `mfrmr`.

## Generate the Fixture

Use the built-in fixture exporter to generate deterministic comparison files:

```bash
python streamlit_app.py --export-parity-fixture validation/generated/parity_fixture
```

The generated folder includes:

- Python outputs for JMLE RSM, JMLE PCM, MML RSM, MML latent regression, and MML GPCM scenarios.
- sirt-specific person-rater response fixture files: `sirt_rater_facets_response.csv` and `sirt_rater_facets_items.csv`.
- `python_parity_manifest.csv`.
- `cross_package_validation_plan.csv`.
- `cross_package_parameterization_notes.csv`.
- `cross_package_tolerance_policy.csv`.
- `r_crosscheck_scaffold.R`.
- notes explaining why exact equality is not expected.

Generated files are intentionally ignored by Git under `validation/generated/`.

## Cross-Package Validation Matrix

| Python area | External reference | What to compare | What not to claim |
| --- | --- | --- | --- |
| Data preparation and category support | all packages | row counts, score map, category support, missingness pattern | parameter equality before score support matches |
| JMLE RSM/PCM | FACETS-like fixed-effect workflows and ConQuest-style calibration | centered estimates, ordering, step behavior, fit outliers | exact equality without matching constraints and extreme-score handling |
| MML RSM | TAM `tam.mml.mfr` style faceted MML | convergence, ordering, EAP direction, broad fit behavior | direct likelihood equality unless quadrature, priors, variance, constants, and constraints match |
| Latent regression | TAM latent regression; mirt `mixedmirt` latent-regression concepts | covariate coding, coefficient direction, group-level EAP shifts | TAM/mirt identity while this app uses fixed `population_prior_sd` behavior |
| GPCM | mirt item-level GPCM checks | positive slopes, slope ordering, EAP ordering, item-level response behavior | raw slope equality without a facet-to-item parameterization map |
| Rater facets | sirt `rm.facets` | rater severity ordering, convergence, fit-log availability | exact equality when slope options and constraints differ |
| Plausible values | TAM, mirt, and sirt plausible-value or factor-score routines | distribution means, variances, covariate trends | draw-by-draw plausible-value equality |
| Anchor/linking | FACETS/TAM/ConQuest-style anchor workflows | anchor count, connectedness, hard-constraint drift | linking claims without common-scale evidence |
| Strict marginal diagnostics | TAM/mirt residual and fit diagnostics where applicable | directional flags and sparse-cell warnings | proof of model truth |

The exporter writes a machine-readable version of this plan to
`cross_package_validation_plan.csv`.

## Default Tolerance Policy

- Exact row counts, rating-category support, score recoding, and missingness flags should match.
- For comparable centered fixed-effect measures, start review at absolute differences greater than about 0.05 logits. Tighten or relax only after documenting package constraints and identification constants.
- For rank ordering of comparable Rasch-family effects, review rank correlations below 0.95.
- For item-level GPCM checks, review rank correlations below 0.90 because the comparison is not the same arbitrary-facet parameterization.
- For latent-regression coefficients, review sign reversals or standardized differences greater than about 0.10.
- For plausible-value distributions, review mean shifts greater than about 0.10 logits or variance ratios outside 0.80 to 1.25.
- Treat log-likelihood differences as interpretable only when category coding, quadrature, priors, constraints, and omitted constants are aligned.
- Do not compare slopes, latent variances, anchored constants, or package-specific nuisance terms as if they were identical parameters without a parameterization map.

## External Reference Roles

- TAM: faceted MML design, latent regression, EAP, and multifacet reference checks.
- mirt: GPCM, EAP/factor-score, plausible-value, and broader IRT diagnostic reference checks.
- sirt: rater-facet, hierarchical rater-model, and plausible-value reference checks.
- mfrmr: functional capability reference for migration completeness.

## Official Documentation Touchpoints

- TAM `tam.mml.mfr`: faceted MML design with `formulaA`, `formulaY`, facets, constraints, and variance controls. See the CRAN TAM manual: https://cran.r-project.org/web/packages/TAM/TAM.pdf
- mirt `mirt`: itemtype choices include Rasch/PCM-style and GPCM-style models; `quadpts`, `dentype`, and optimizer choices affect comparability. See the mirt reference: https://philchalmers.github.io/mirt/reference/mirt.html
- mirt `fscores`: EAP is the default factor-score method, and plausible-value support is exposed through `plausible.draws`. See: https://philchalmers.github.io/mirt/reference/fscores.html
- mirt `mixedmirt`: latent-regression inputs are modeled through `lr.fixed` / `lr.random` style arguments. See: https://philchalmers.github.io/mirt/docs/reference/mixedmirt.html
- sirt `rm.facets`: rater-facet models use person-rater rows and optional item/rater intercept and slope settings. See the CRAN sirt manual: https://cran.r-project.org/web/packages/sirt/sirt.pdf
- sirt plausible-value tools: distributional checks should be preferred over draw-by-draw equality. See the CRAN sirt manual: https://cran.r-project.org/web/packages/sirt/sirt.pdf

## Reporting Rule

When reporting a cross-package check, archive:

- the generated fixture folder;
- the external package versions;
- the exact R script used;
- the parameterization map;
- the tolerance table;
- the final comparison table; and
- a short note explaining which differences are expected.

Do not report exact cross-package parity unless those artifacts are included.
