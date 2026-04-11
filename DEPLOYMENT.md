# Deployment Guide

This app is safest for confidential rating data when run locally. Hosted deployment is appropriate for public demonstrations, teaching data, synthetic data, or non-confidential research prototypes.

## Local Deployment

Use local execution for confidential, regulated, or institutionally restricted data:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
streamlit run streamlit_app.py
```

Before sharing results, remove direct identifiers where possible and keep only the columns required for the analysis.

## Streamlit Community Cloud

Use these settings when deploying from GitHub:

- Repository: `Ryuya-dot-com/MFRM_Python_Application`
- Branch: `main`
- Main file path: `streamlit_app.py`
- Python dependencies: `requirements.txt`
- Development-only dependencies: `requirements-dev.txt` is for CI and local verification, not the hosted runtime.
- Secrets: none are required by the app. Do not add confidential data to `.streamlit/secrets.toml`.

The root `requirements.txt` is intentional because Streamlit Community Cloud installs app dependencies from a requirements file in the repository root or next to the entrypoint file. A `packages.txt` file is not currently needed for core analysis because the app falls back to interactive HTML figure exports if Chrome or Chromium is unavailable for Kaleido PNG export.

Set the Python version deliberately in the Streamlit Community Cloud Advanced settings. After deployment, changing the Python version requires deleting and redeploying the app, so record the chosen version with the app release notes.

## Privacy Gate Before Hosted Deployment

Do not use a hosted deployment for:

- confidential operational scoring data
- person-level records with direct identifiers
- regulated or contract-restricted assessment data
- data whose sharing policy has not been confirmed

For public demos, use synthetic data, the built-in example, or fully de-identified teaching data.

## Preflight Checks

Run these checks before deployment:

```bash
make verify
make clean
```

Then confirm:

- `git status --short` is clean.
- `validation/generated/` is absent.
- `__pycache__` and `.pytest_cache` are absent.
- `.streamlit/secrets.toml` is not tracked.
- README and in-app privacy warnings still appear.

## Official References

- Streamlit Community Cloud app dependencies: https://docs.streamlit.io/deploy/streamlit-community-cloud/deploy-your-app/app-dependencies
- Streamlit Community Cloud Python version management: https://docs.streamlit.io/deploy/streamlit-community-cloud/manage-your-app/upgrade-python
