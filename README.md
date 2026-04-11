# MFRM Streamlit

Standalone Python beta application for Many-Facet Rasch Model estimation in Streamlit.

This app is designed to run without `mfrmr`, `rpy2`, `Rscript`, FACETS, TAM, sirt, or mirt at runtime. Those tools are used only as external methodological references for validation and interpretation.

## Status

- Release status: public beta / research preview
- Runtime engine: standalone Python
- Primary entrypoint: `streamlit_app.py`
- Intended use: exploratory analysis, teaching, reporting support, and research workflow prototyping
- Not intended as: a validated drop-in replacement for FACETS, TAM, sirt, mirt, or `mfrmr`

Before using results for high-stakes scoring, placement, certification, employment, or institutional decisions, cross-check the analysis with an established workflow and document the model assumptions.

## Data Privacy

Rating data can contain person IDs, rater IDs, school or institution identifiers, subgroup labels, and other sensitive information.

For confidential data, run the app locally on a controlled machine:

```bash
streamlit run streamlit_app.py
```

If the app is deployed to Streamlit Community Cloud or another hosted platform, uploaded data may be processed on remote infrastructure. Confirm your institutional data-handling requirements before upload. Prefer removing direct identifiers, limiting uploaded columns to analysis variables, and avoiding hosted deployments for regulated or confidential datasets.

The app is intended to keep uploaded data in memory during the session unless the user explicitly downloads or exports outputs.

## Install

Use a virtual environment:

```bash
cd /Users/ryuya/Library/CloudStorage/Dropbox/MFRM_Application/MFRM_App/MFRM_Streamlit
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

For local verification and CI-equivalent smoke tests, install the development dependencies:

```bash
python -m pip install -r requirements-dev.txt
```

## Run

```bash
streamlit run streamlit_app.py
```

## Verify

```bash
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --self-test
python -m pytest tests/test_app_smoke.py
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
python streamlit_app.py --export-parity-fixture validation/generated/parity_fixture
```

Optional Make shortcuts:

```bash
make verify
make clean
make run
```

## Data Format

Use long-format rating data:

```csv
Person,Rater,Task,Criterion,Score
P01,R1,T1,C1,4
P01,R2,T1,C1,5
P02,R1,T2,C2,2
P02,R2,T2,C2,3
```

Required:

- Person column
- Score column with ordered integer categories
- One or more facet columns such as rater, task, criterion, prompt, form, or occasion

Supported score handling includes:

- Intended rating boundaries, such as 1 to 5
- Zero-count intermediate categories, such as a 1-5 scale where category 3 is absent
- Optional recoding of non-consecutive observed categories
- Rejection of fractional categories with a clear error

## Model Scope

Implemented in the standalone Python engine:

- RSM
- PCM
- bounded GPCM
- JMLE
- MML via EM, Direct, Hybrid, and Auto engine selection
- latent regression via constrained `population_formula`
- fixed user-set population prior SD for the current MML population model path
- EAP posterior scoring
- plausible values
- strict marginal diagnostics
- residual PCA diagnostics
- bias / local interaction screening
- anchor audit and linking review
- anchor drift and equating-chain summaries
- prediction for fitted, held-out, and scenario rows
- simulation and design evaluation
- final-report readiness checklist
- downloadable tables, configuration, scripts, and method appendix
- reproducibility/config fingerprints for analysis exports

## Beginner Reporting Workflow

After fitting a model, inspect results in this order:

1. Convergence: do not interpret final measures until the optimizer converges.
2. Category functioning: check sparse categories, monotonic average measures, and threshold order.
3. Reliability / separation: confirm whether the design supports stable person and facet ordering.
4. Wright map targeting: check whether person locations and facet difficulty/severity ranges overlap.
5. Fit diagnostics: review large standardized residuals and misfitting elements.
6. Bias / local interaction: treat flags as review prompts, not automatic proof of bias.
7. PCA / dimensionality: check whether residual structure suggests a second dimension.
8. Anchor / linking review: check connectedness and anchor stability before comparing runs or groups.
9. Strict marginal diagnostics: use for final MML reports when feasible.
10. Final-report readiness: use the generated checklist before writing conclusions.

## Statistical Caveats

- GPCM is not a strict Rasch model. Its slope parameters change the interpretation of invariance and should be reported explicitly.
- The current latent regression path uses the app's documented population prior SD behavior. Do not describe it as identical to TAM unless the model, quadrature, variance treatment, and constraints have been checked.
- Bias and differential functioning outputs are screening tools unless linking, common-scale evidence, sample size, and precision support stronger claims.
- Cross-package equality is not expected by default because FACETS, TAM, sirt, mirt, and this app can use different parameterizations, constraints, latent variance handling, and optimization details.
- Treat external R cross-check outputs as validation evidence, not as a runtime dependency.
- Reproducibility fingerprints help confirm that a report came from the same data/settings, but they are not a privacy guarantee or encryption method.

## Caching

- `st.cache_resource` is used only for the read-only core function namespace.
- `st.cache_data` is used for built-in sample data, bundled anchor templates/guidelines, and short-lived export bytes keyed by reproducibility fingerprints.
- Uploaded or pasted rating data are not cached as a resource. Export-byte caches are bounded and short-lived, but confidential analyses should still be run locally.

## External Reference Roles

- TAM: faceted MML design, latent regression, EAP, and multifacet reference checks.
- mirt: GPCM, EAP/factor scores, plausible values, and broader IRT diagnostics.
- sirt: rater-facet and hierarchical rater model reference checks.
- mfrmr: functional capability reference for the Python migration target.

## Continuous Integration

The app-specific workflow is in:

```text
.github/workflows/python-streamlit.yml
```

It runs:

- dependency installation
- `python -m py_compile streamlit_app.py`
- `python streamlit_app.py --doctor`
- `python streamlit_app.py --self-test`
- `python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv`
- Streamlit AppTest smoke check
- parity fixture export smoke check

If this directory is used as a standalone GitHub repository, the workflow will be discovered normally. If it remains a subdirectory inside a larger repository, copy or mirror the workflow into the repository root `.github/workflows/` directory.

## License

MIT License. See `LICENSE`.

## Publication Note

`origin` is set to `https://github.com/Ryuya-dot-com/MFRM_Python_Application.git`. That repository already had a separate `main` history when this standalone Streamlit distribution was prepared. The intended publication path is an explicit replacement of `origin/main` with this standalone app history.

## Files

```text
Makefile
streamlit_app.py
requirements.txt
requirements-dev.txt
LICENSE
README.md
CHANGELOG.md
CONTRIBUTING.md
SECURITY.md
MFRM_STREAMLIT_RELEASE_PLAN.md
RELEASE_CHECKLIST.md
anchor_templates_and_guideline/
.streamlit/config.toml
.github/workflows/python-streamlit.yml
.github/ISSUE_TEMPLATE/
```

Generated files such as `__pycache__`, benchmark CSVs, local logs, and validation output directories are ignored.
