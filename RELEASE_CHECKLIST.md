# Release Checklist

Use this checklist before tagging or pushing a public beta release.

## Scope

- [ ] Confirm this repository is being used as an independent repo, not accidentally committed as ordinary files inside the parent `MFRM_Application` repo.
- [ ] If the parent repo should reference this project, decide between a Git submodule and a parent-level ignore rule.
- [ ] Confirm the release label remains `standalone Python beta` unless the validation scope has changed.
- [ ] Confirm the README still states that the app is not an exact replacement for FACETS, TAM, sirt, mirt, or `mfrmr`.

## Privacy

- [ ] Confirm no example or generated files contain confidential rating data.
- [ ] Confirm generated outputs under `validation/generated/` are absent before committing.
- [ ] Confirm `.streamlit/secrets.toml` is not tracked.
- [ ] Confirm the in-app and README privacy warnings are still visible.

## Deployment

- [ ] Review `DEPLOYMENT.md` before hosted deployment.
- [ ] Confirm hosted demos use synthetic, built-in, or fully de-identified data.
- [ ] Record the Streamlit Community Cloud Python version selected in Advanced settings.

## Statistical Validation

- [ ] Run `make verify` or the equivalent commands below.
- [ ] Review `validation/README.md` before claiming cross-package parity.
- [ ] If external R checks are reported, archive the generated `r_crosscheck_status.csv` and note package versions.
- [ ] Do not claim exact numerical parity unless the fixture, tolerances, and parameterization map are included.

## Commands

```bash
python -m pip install -r requirements-dev.txt
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --self-test
python -m pytest tests/test_app_smoke.py
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
python streamlit_app.py --export-parity-fixture validation/generated/parity_fixture
rm -rf .pytest_cache validation/generated
find . -type d -name __pycache__ -prune -exec rm -rf {} +
```

## GitHub Release

- [x] Choose a license before public distribution.
- [x] Add a remote with `git remote add origin <repo-url>`.
- [x] Because `origin/main` already has a separate history, choose a publication path: branch/PR, new repository, or explicit replacement.
- [x] Explicit replacement of `origin/main` approved by the project owner.
- [x] Push with `git push --force-with-lease origin main`.
- [x] Confirm GitHub Actions pass on Python 3.11 and 3.12.
- [x] Add a tag only after CI passes.
- [x] Create GitHub prerelease `v0.1.0-beta`.
