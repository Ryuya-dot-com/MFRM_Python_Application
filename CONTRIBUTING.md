# Contributing

This repository is a standalone Python beta for MFRM analysis in Streamlit. Contributions should preserve the app's current scope: a single-file application that does not call `mfrmr`, `rpy2`, `Rscript`, FACETS, TAM, sirt, or mirt at runtime.

## Before Editing

- Run `python streamlit_app.py --doctor` to confirm the local environment.
- Do not commit generated files such as `__pycache__`, `.pytest_cache`, or `validation/generated/`.
- Do not add confidential rating data, real person IDs, rater IDs, institution names, subgroup labels, or proprietary assessment data.
- Keep privacy warnings and beta/parity caveats visible in user-facing documentation.

## Verification

Use:

```bash
python -m pip install -r requirements-dev.txt
make verify
make clean
```

Equivalent commands:

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

## Statistical Changes

For changes affecting estimation, diagnostics, prediction, simulation, anchoring, or reporting:

- Add or update a self-test where feasible.
- Update `validation/README.md` if the external-comparison scope changes.
- Keep external package comparisons directional unless the fixture, tolerance, and parameterization map are documented.
- Document differences from FACETS, TAM, sirt, mirt, and `mfrmr` rather than presenting the app as an exact replacement.

## Documentation Changes

When changing the UI or supported workflow, update at least one of:

- `README.md`
- `MFRM_STREAMLIT_RELEASE_PLAN.md`
- `RELEASE_CHECKLIST.md`
- `CHANGELOG.md`

## Pull Request Checklist

- [ ] No confidential data or generated artifacts are committed.
- [ ] `python streamlit_app.py --doctor` passes.
- [ ] `python streamlit_app.py --self-test` passes.
- [ ] Streamlit AppTest smoke passes.
- [ ] Statistical caveats and privacy warnings remain accurate.
