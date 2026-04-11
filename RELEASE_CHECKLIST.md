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

## Statistical Validation

- [ ] Run `make verify` or the equivalent commands below.
- [ ] Review `validation/README.md` before claiming cross-package parity.
- [ ] If external R checks are reported, archive the generated `r_crosscheck_status.csv` and note package versions.
- [ ] Do not claim exact numerical parity unless the fixture, tolerances, and parameterization map are included.

## Commands

```bash
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --self-test
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
python streamlit_app.py --export-parity-fixture validation/generated/parity_fixture
python -c 'from streamlit.testing.v1 import AppTest; at = AppTest.from_file("streamlit_app.py").run(timeout=30); assert not at.exception'
rm -rf __pycache__ .pytest_cache validation/generated
```

## GitHub Release

- [ ] Choose a license before public distribution.
- [ ] Add a remote with `git remote add origin <repo-url>`.
- [ ] Push with `git push -u origin main`.
- [ ] Confirm GitHub Actions pass on Python 3.11 and 3.12.
- [ ] Add a tag only after CI passes.
