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

## Preview

![MFRM Streamlit sample-data result overview](docs/images/app-data-overview.png)

The screenshot uses the built-in synthetic sample data. It highlights the
default guided sidebar, the visible data-privacy warning, the input preview,
the post-estimation success status, and the beginner-oriented result tabs.

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

For a first local smoke check, keep **Sample data (built-in)** selected, leave
the guided defaults unchanged, click **Run FACETS-mode estimation**, then open
the **What should I look at first?**, **Data**, **Visuals**, and **Help** tabs.

For deployment notes, including Streamlit Community Cloud settings and privacy gates for hosted use, see `DEPLOYMENT.md`.

## Verify

```bash
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --release-check
python streamlit_app.py --self-test
python -m pytest tests
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
python streamlit_app.py --export-demo-report validation/generated/demo_report
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
- visual interpretation checklist for beginner-friendly figure reading
- latent-regression covariate type preview for numeric IDs/codes that may need categorical coding
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
11. Manuscript claim guide: check what is safe to claim, what requires a caveat, and what should not be claimed yet.

The final-report readiness checklist, first-read guide, and generated report
text use the same main thresholds: about 5% or fewer observation residuals
with `|z| >= 2`, person reliability at least 0.80 when person separation is the
goal, and residual PCA first eigenvalue below 2.0 for a clean screen.

For MML latent regression, inspect the covariate type preview before fitting.
Integer-like columns such as `GradeCode = 1, 2, 3` are flagged so you can decide
whether they are continuous predictors or category labels that should be forced
categorical.

The app and demo export also include a manuscript claim guide, public-beta
limitations, and release readiness tables. Use these to keep public claims
aligned with what the standalone Python engine currently supports.
The same public-beta bundle includes `mfrmr_015_migration_coverage.csv`, which
maps the local mfrmr 0.1.5 feature surface to Python support, boundaries, and
next validation actions.

The Visuals tab also includes a downloadable visual interpretation checklist.
It maps each figure to the first signal to read, the review trigger, and the
recommended next action for beginners.
It also includes a visual method evidence table that links each plot family to
its Rasch/MFRM diagnostic role and explains the app's readability rules.
For PCM and bounded GPCM, category probability curves can be read either as an
averaged overview or for each selected step-facet level. The downloadable table
bundle also includes long-form curve data for all available curve scopes.

To generate a synthetic beginner-facing report without uploading data:

```bash
python streamlit_app.py --export-demo-report validation/generated/demo_report
```

Open `validation/generated/demo_report/MFRM_Demo_Report.html` first, then read
`final_report_readiness.csv`, `manuscript_claim_guide.csv`,
`visual_interpretation_checklist.csv`, `visual_method_evidence.csv`,
`public_beta_limitations.csv`, `mfrmr_015_migration_coverage.csv`,
`public_release_readiness.csv`, and the interactive category curves in
`figures_html/`.

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

## Figure Exports

Interactive HTML figure exports are always attempted. Static PNG exports use
Plotly/Kaleido and require Chrome or Chromium in the runtime when using
`kaleido >= 1.0`; if that browser dependency is unavailable, the app falls back
to interactive HTML figures instead of blocking the analysis.

## External Reference Roles

- TAM: faceted MML design, latent regression, EAP, and multifacet reference checks.
- mirt: GPCM, EAP/factor scores, plausible values, and broader IRT diagnostics.
- sirt: rater-facet and hierarchical rater model reference checks.
- mfrmr: functional capability reference for the Python migration target.

For the cross-package validation matrix and tolerance policy, see `validation/README.md`.
For the latest archived local R smoke status, see
`validation/R_CROSSCHECK_STATUS.md`.
The parity fixture also includes official documentation touchpoints, an
artifact checklist, an external validation report template, and the mfrmr 0.1.5
migration coverage table so external checks can be archived without making
R packages a runtime dependency.

## Continuous Integration

The app-specific workflow is in:

```text
.github/workflows/python-streamlit.yml
```

It runs:

- dependency installation
- `python -m py_compile streamlit_app.py`
- `python streamlit_app.py --doctor`
- `python streamlit_app.py --release-check`
- `python streamlit_app.py --self-test`
- `python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv`
- Streamlit AppTest smoke check
- demo report export smoke check
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
DEPLOYMENT.md
CHANGELOG.md
CONTRIBUTING.md
SECURITY.md
MFRM_STREAMLIT_RELEASE_PLAN.md
RELEASE_CHECKLIST.md
anchor_templates_and_guideline/
.streamlit/config.toml
.github/workflows/python-streamlit.yml
.github/ISSUE_TEMPLATE/
docs/images/app-data-overview.png
```

Generated files such as `__pycache__`, benchmark CSVs, local logs, and validation output directories are ignored.
