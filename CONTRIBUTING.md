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
python -m pytest tests
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

## Localization (i18n)

User-visible strings should resolve through the ``t()`` helper so the
sidebar language picker (English / Japanese) keeps working. New strings:

- Add the English source under `locales/en.json` first; the locale code
  in `_metadata.code` must match the filename stem.
- Mirror the same dotted key in `locales/ja.json` with the Japanese
  translation. `tests/test_i18n_parity.py` enforces key, placeholder, and
  non-empty parity in CI; treat a failing parity test as a release
  blocker.
- Use `snake_case` for key segments and at most three dot-separated
  levels (e.g. `tab.diagnostics.fit.header`); deeper trees should be
  refactored into a new section.
- Keep MFRM/statistical jargon in English (``MFRM``, ``RSM``, ``PCM``,
  ``GPCM``, ``JMLE``, ``MML``, ``Infit``, ``Outfit``, ``MNSQ``, ``ZSTD``,
  ``Person``, ``Rater``, ``Criterion``, ``Item``, ``Task``, ``Region``).
  Translate generic UI vocabulary (``preview``, ``upload``, ``data``,
  ``pipeline``, ``diagnostics``, ``sub-tab``).
- Reuse Python ``str.format`` placeholders verbatim across locales
  (e.g. ``{release_label}``); never translate placeholder names.
- Do not commit translation tooling artefacts (``.po``, ``.mo``,
  ``poedit`` caches). The repo is plain JSON to keep Streamlit Community
  Cloud and GitHub Pages deployments dependency-free.
- Indent locale JSON with 2 spaces, UTF-8 without BOM, and end with a
  trailing newline.

### Plotly figure text

All Plotly figure text (chart `title`, `xaxis_title`, `yaxis_title`,
`annotation_text`, legend entries, hover templates) stays in English
regardless of locale. Streamlit re-renders figures every rerun; keeping
them in English ensures publication exports (PNG/SVG/PDF via Kaleido or
matplotlib) reproduce the wording cited in the manuscript and supports
side-by-side comparison with FACETS / Winsteps screenshots.

When the same string is used in both a Plotly figure and an adjacent
`st.markdown` / `st.warning` label, hold both forms — see
`_show_pca_panel` in `streamlit_app.py` for the canonical pattern:

```python
if mode == "overall":
    label_plot = "Overall standardized residuals"           # English, for Plotly title
    label_ui = t("dimensionality.label_overall")            # Translated, for st.markdown
else:
    label_plot = f"Residuals for facet: {facet_name}"
    label_ui = t("dimensionality.label_facet_template", facet_name=facet_name)
```

### Statistical term translations

The following terms have established Japanese translations in the
psychometric / language-testing literature and are translated rather
than kept in English in `ja.json`:

| English source | JA translation |
|---|---|
| eigenvalue | 固有値 |
| component (in PCA contexts) | 主成分 |
| residual | 残差 |
| variance explained | 説明分散 |
| dimension(ality) | 次元 (性) |

Conversely, model names, estimator acronyms, Greek letters in math
notation, and Rasch-specific facet names (Person, Rater, Task,
Criterion, Item, Facet, Score) stay in English. When in doubt, prefer
keeping the term in English and revisit during translation review.

## Pull Request Checklist

- [ ] No confidential data or generated artifacts are committed.
- [ ] `python streamlit_app.py --doctor` passes.
- [ ] `python streamlit_app.py --self-test` passes.
- [ ] Streamlit AppTest smoke passes.
- [ ] Statistical caveats and privacy warnings remain accurate.
