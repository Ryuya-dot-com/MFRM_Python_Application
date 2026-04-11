# Validation Notes

This directory is for optional validation artifacts. The Streamlit app does not
call R packages or external MFRM engines at runtime.

Use the built-in fixture exporter to generate deterministic comparison files:

```bash
python streamlit_app.py --export-parity-fixture validation/generated/parity_fixture
```

The generated folder includes:

- Python outputs for JMLE RSM, JMLE PCM, MML RSM, and MML latent regression scenarios
- sirt-specific person-rater response fixture files: `sirt_rater_facets_response.csv` and `sirt_rater_facets_items.csv`
- `python_parity_manifest.csv`
- `r_crosscheck_scaffold.R`
- notes explaining why exact equality is not expected

Reference roles:

- TAM: faceted MML and latent-regression reference
- mirt: GPCM, EAP/factor-score, and plausible-value reference
- sirt: rater-facet and hierarchical rater-model reference
- mfrmr: functional capability reference

Interpretation rule:

- Passing cross-checks supports implementation credibility.
- Parameter values may differ because of constraints, parameterization, scoring support, quadrature, optimizer settings, and latent variance treatment.
- Do not claim exact package parity unless a specific fixture, tolerance, and report are included with the analysis.

Suggested tolerance policy:

- Exact row counts, rating-category support, score recoding, and missingness flags should match.
- For comparable centered fixed-effect measures, start review at absolute differences greater than about 0.05 logits, then tighten or relax only after documenting the package constraints and identification constants.
- Do not compare slopes, latent variances, anchored constants, or package-specific nuisance terms as if they were identical parameters without a parameterization map.
- Treat log-likelihood differences as interpretable only when category coding, quadrature, priors, constraints, and omitted constants are aligned.
