# MFRM Streamlit Release Plan

## Goal

- Build a clean, standalone Streamlit distribution of the current Python MFRM app in this directory.
- Keep the application independent from `mfrmr`, `rpy2`, `Rscript`, FACETS, TAM, sirt, and mirt at runtime.
- Treat TAM, sirt, mirt, and `mfrmr 0.1.5` as external statistical references for validation language and parity checks, not as engines called by the app.
- Make the app suitable for a public beta release, while avoiding any claim of exact numerical parity with established statistical packages until cross-package tolerances are documented.

## Directory Scope

- Target directory: `/Users/ryuya/Library/CloudStorage/Dropbox/MFRM_Application/MFRM_App/MFRM_Streamlit`
- Primary app entrypoint: `streamlit_app.py`
- Documentation: `README.md`
- Dependency file: `requirements.txt`
- CI: `.github/workflows/python-streamlit.yml`
- Optional Streamlit config: `.streamlit/config.toml`
- Validation artifacts: `validation/`
- Tests or smoke wrappers: `tests/`
- Repository hygiene: `.gitattributes`, `Makefile`, `RELEASE_CHECKLIST.md`
- Do not copy `__pycache__`, local benchmark outputs, virtual environments, or temporary parity outputs.

## Release Positioning

- Public label: `standalone Python beta`
- Safe claim: "Functional parity target with key MFRM workflows."
- Avoid claim: "Exact replacement for FACETS/TAM/sirt/mirt/mfrmr."
- Required caveat: results should be cross-checked before high-stakes operational decisions.
- Required privacy caveat: uploaded rating data may contain examinee, rater, or institution information. Prefer local execution for confidential data. If deployed to Streamlit Community Cloud or another hosted service, users must confirm their data handling obligations before uploading.

## Phase 1: Clean Copy

- [x] Copy `MFRM_Python_Application/streamlit_app.py` to `MFRM_Streamlit/streamlit_app.py`.
- [x] Copy anchor templates/guidelines into `anchor_templates_and_guideline/`.
- [x] Create a fresh `requirements.txt` in the target directory.
- [x] Exclude `__pycache__` and any generated local artifacts.
- [x] Add a local `.gitignore` if the target directory is expected to be versioned independently.
- [x] Run `python3 -m py_compile streamlit_app.py`.
- [x] Run `python3 streamlit_app.py --self-test`.

## Phase 2: README and Public Documentation

- [x] Create a new README that matches the actual entrypoint: `streamlit run streamlit_app.py`.
- [x] Document install and run steps for local execution.
- [x] Document the current model scope:
  - [x] RSM
  - [x] PCM
  - [x] bounded GPCM
  - [x] JMLE
  - [x] MML EM / Direct / Hybrid / Auto
  - [x] latent regression / `population_formula`
  - [x] plausible values
  - [x] prediction
  - [x] simulation and design evaluation
  - [x] anchor audit / linking review
  - [x] strict marginal diagnostics
- [x] Add "What beginners should inspect first" with a short report workflow.
- [x] Add "Known limitations" with statistical caveats:
  - [x] fixed user-set population prior SD unless variance estimation is explicitly implemented
  - [x] bounded GPCM scope and identification conventions
  - [x] bias and DFF outputs as screening unless linking and precision support stronger claims
  - [x] parity fixture outputs are directional checks, not exact equality guarantees
- [x] Add a privacy and data-governance section.
- [x] Add a verification section listing `py_compile`, `--self-test`, benchmark smoke, and optional external R cross-check scaffold.

## Phase 3: App CI

- [x] Add a Python GitHub Actions workflow in `.github/workflows/python-streamlit.yml`.
- [x] CI steps:
  - [x] checkout
  - [x] set up Python 3.11 or 3.12
  - [x] install `requirements.txt`
  - [x] run `python -m py_compile streamlit_app.py`
  - [x] run `python streamlit_app.py --doctor`
  - [x] run `python streamlit_app.py --self-test`
  - [x] run a lightweight benchmark smoke if runtime is acceptable
  - [x] run a Streamlit AppTest smoke once the app can run without manual inputs
- [x] Keep R cross-checks optional because they add package installation time and are not runtime dependencies.

## Phase 4: Data Privacy Warning

- [x] Add prominent in-app warning near file upload / paste input.
- [x] State that confidential rating data should be run locally.
- [x] State that hosted deployments may process data on remote infrastructure.
- [x] Make the warning practical: remove direct identifiers where possible, avoid uploading unnecessary columns, and export/delete local results according to institutional policy.
- [x] Mirror the same warning in README.
- [x] Ensure no uploaded data are written to disk unless the user explicitly exports/downloads files.

## Phase 5: Streamlit Caching

- [x] Use `st.cache_data` for deterministic, trusted, pure data transformations when the cache key can represent the input safely.
- [x] Candidate cache targets:
  - [x] built-in sample data generation
  - [x] static anchor template loading
  - [x] method/help text assembly
  - [x] expensive deterministic figure/table exports when keyed by an explicit config fingerprint
- [x] Be cautious caching uploaded user data because Streamlit `st.cache_data` serializes returned values. If used, keep cache entries bounded with `ttl` and `max_entries`.
- [x] Use `st.cache_resource` only for the read-only core namespace; do not cache uploaded user data as a resource.
- [x] Add a small cache invalidation helper based on a reproducibility/config fingerprint.
- [x] Run self-tests after every cache-related change to prevent stale output.

## Phase 6: Statistical Refinements

- [x] Strengthen validation language against external references:
  - [x] TAM: MML/faceted design and latent regression reference
  - [x] mirt: GPCM, EAP/factor scores, plausible values reference
  - [x] sirt: rater-facet and hierarchical rater model reference
  - [x] mfrmr 0.1.5: functional capability reference
- [x] Add or improve a `validation/` fixture workflow:
  - [x] deterministic Python fixture generation
  - [x] R scaffold for TAM/mirt/sirt checks
  - [x] README explaining expected parameterization differences
  - [x] numeric tolerance policy where comparable
- [x] Add a sirt-specific rater-facet fixture instead of relying only on the current arbitrary long-facet fixture.
- [x] Add clear reporting language for cases where exact package parity is not expected.
- [x] Review GPCM documentation and UI warning so users understand it is outside strict Rasch invariance.
- [x] Review latent regression documentation so users understand the population distribution and prior SD assumptions.
- [x] Review strict marginal diagnostics thresholds and label them as diagnostic screens rather than proof of model truth.

## Phase 7: Final Public-Beta Gate

- [x] Fresh target directory contains no `__pycache__`.
- [x] `streamlit_app.py` compiles.
- [x] `--doctor` passes.
- [x] `--self-test` passes.
- [x] App starts with `streamlit run streamlit_app.py`.
- [x] README run command is correct.
- [x] Requirements install cleanly in a fresh virtual environment.
- [x] CI file exists and matches the target directory layout.
- [x] Data privacy warning appears in both app and README.
- [x] Release language says beta / functional parity, not exact replacement.
- [x] Known limitations are explicit.
- [x] Validation scaffold exists or is documented as optional.

## Phase 8: Git Repository Setup

- [x] Initialize this directory as an independent Git repository on branch `main`.
- [x] Commit the standalone Streamlit app and support files.
- [x] Add `.gitattributes` for consistent text/binary handling.
- [x] Add `Makefile` shortcuts for local verification and cleanup.
- [x] Add `RELEASE_CHECKLIST.md` for public beta release gates.
- [x] Add `CHANGELOG.md`, `CONTRIBUTING.md`, and `SECURITY.md`.
- [x] Add GitHub issue templates that warn against sharing confidential rating data.
- [ ] Add a GitHub remote after choosing the target repository URL.
- [ ] Choose and add a license before public distribution.

## Immediate Next Step

Completed setup path for this iteration:

1. [x] Copy the app and templates into this clean directory.
2. [x] Create fresh `requirements.txt`, `.gitignore`, README, and CI.
3. [x] Apply privacy warning and safe caching changes.
4. [x] Add statistical validation refinements.
5. [x] Run compile, doctor, self-test, AppTest smoke, benchmark smoke, parity fixture export smoke, optional R scaffold smoke, fresh-venv install/self-test, and Streamlit health check.
6. [x] Initialize local Git history and add repository hygiene files.
7. [x] Add public-beta support documentation and issue templates.
