# Contributing

This repository contains the Streamlit Web application and a Tkinter desktop
preview for MFRM analysis.
Reusable statistical, evidence, and export computation belongs under
`mfrm_app/`; `streamlit_app.py` and `desktop_app.py` are UI entrypoints and may retain thin
compatibility wrappers. Core code must not call external estimation engines.

## Run from source

Use the branch or tag for the version you intend to inspect. A published Web
tag, the default branch and a desktop preview branch may contain different
features. These instructions assume a checkout containing `desktop_app.py`
and `requirements.txt`; GitHub's source ZIP is not a standalone executable.
Users of a supplied application ZIP do not need Python or these steps.

Use Python 3.11 or 3.12. Desktop startup also requires Tkinter/Tk. From the
source directory, create an isolated environment and install the requirements.
The initial installation needs internet access.

On Windows (PowerShell):

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe desktop_app.py
```

On macOS or Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python desktop_app.py
```

Keep the terminal open while using the source application. For subsequent
desktop launches, run the last command in the same environment. For the Web
interface, run `python -m streamlit run streamlit_app.py` instead, using the
environment's Python executable. See the [Windows guide](docs/desktop_user_guide.en.md)
for the sample workflow. Source startup does not establish support for a
packaged app on a particular OS.

## Before Editing

- When setting up or changing the runtime environment, use
  `python streamlit_app.py --doctor` to check it.
- Do not commit generated files such as `__pycache__`, `.pytest_cache`, or `validation/generated/`.
- Do not add confidential rating data, real person IDs, rater IDs, institution names, subgroup labels, or proprietary assessment data.
- Keep privacy warnings and statistical claim boundaries visible in user-facing documentation.
- Read `ROADMAP.md` before expanding product scope. EGA, external-engine
  execution, result ingestion, and cross-engine interchange are not core work.

## Verification

Choose checks for the changed behavior and dependencies. For documentation-only
changes, verify facts, links, language consistency and whether the text is
bundled; do not repeat numerical experiments or rebuild unchanged computation.
For application changes, run affected tests and the required integration checks.
The general Web verification sequence is:

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
make apptest
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
make clean
```

`make clean` removes development caches while preserving retained validation
evidence, source snapshots, backups, and both Mac and Windows virtual environments.
In a Dropbox workspace, finish synchronization before bulk cleanup.

## Statistical Changes

For changes affecting estimation, diagnostics, prediction, simulation,
anchoring, or reporting:

- Add or update a self-test where feasible.
- Put reusable computation in a focused `mfrm_app/` module without importing
  Streamlit.
- Return a stable ReasonCode when a prerequisite fails or a result is not
  assessable.
- Add deterministic native-Python fixtures and document the interpretation
  boundary of every automatic recommendation.
- Keep computation state and conclusion-stability state separate; neither is a
  pass/fail label for a person, rater, institution, or assessment.
- Prespecify sensitivity variants by exact `AnalysisID`, not a display label;
  use an executable versioned rule and restore decisions with their complete
  evidence bundle so missing or altered references fail closed.

## Legacy Compatibility

External-engine code generators, posterior-result viewers, and handoff
archives are frozen compatibility surfaces and are not the design basis for new
features. The parity-fixture CLI and Make target have been removed. The normal
`make apptest`, `make verify`, and GitHub workflow exclude tests marked
`legacy_compat`; those tests preserve dormant code only while the explicit
retain/deprecate/remove decision in `ROADMAP.md` remains open. They are not
release evidence. Do not expand or reconnect these surfaces without a separate
product-scope decision.

For compatibility-maintenance work only, the isolated tests can be inspected
with `python -m pytest -m legacy_compat tests`. This command is not part of the
standalone release gate.

## Documentation Changes

When changing the UI or supported workflow, update the current instructions
that users encounter in the README, in-app Help, or bundled guide, and record
the change in `CHANGELOG.md`. Update `ROADMAP.md` when scope or priorities
change, and `RELEASE_CHECKLIST.md` when release requirements change.

`MFRM_STREAMLIT_RELEASE_PLAN.md`, `docs/v0.2.0_beta_plan.md`, and
`docs/advanced_models_guide.md` are historical records. Updating a historical
plan does not update current user documentation. Preserve dated findings and
frozen evidence; make corrections or superseding decisions explicit. Check
whether linked evidence is tracked in Git or retained only locally before
presenting it as available to readers of a fresh checkout.

## Localization (i18n)

In the Web interface, user-visible strings should resolve through the ``t()`` helper so the
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

Desktop interface strings use `mfrm_app/desktop_i18n.py` and the desktop
translation helper. Keep its English/Japanese wording consistent with the
Web descriptions without translating user-supplied column names or data.

### Plotly figure text

Publication figure text keeps its established English labels so exported
PNG/SVG/PDF figures can be checked against manuscript descriptions. Interactive
figures used only for navigation or comparing saved analyses may follow the UI
language; keep their axis meaning and numerical values identical across locales.
Do not silently change publication figure wording as part of a UI-copy edit.

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
| regularization | 正則化 |
| ceiling effect | 天井効果 |
| floor effect | 床効果 |

Conversely, model names, estimator acronyms, Greek letters in math
notation, and Rasch-specific facet names (Person, Rater, Task,
Criterion, Item, Facet, Score) stay in English. When in doubt, prefer
keeping the term in English and revisit during translation review.

## Pull Request Checklist

- [ ] No confidential data or generated artifacts are committed.
- [ ] Checks cover the changed behavior, affected saved outputs and required integration paths.
- [ ] Relevant runtime, self-test, Streamlit AppTest or desktop checks pass; their scope and any untested conditions are stated.
- [ ] Statistical caveats and privacy warnings remain accurate.
