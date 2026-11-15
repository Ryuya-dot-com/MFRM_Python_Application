# mfrmpy: Many-facet Rasch analysis in Python

`mfrmpy` brings the supported MFRM Desktop calculations to Python scripts and
notebooks. It uses the same Python engine as the existing MFRM Web and desktop
applications. It does not require R, FACETS, Streamlit or Tkinter for analysis.

The initial API is a **pre-release**. Available analyses are RSM/PCM with JMLE
or fixed-population-SD MML, diagnostics, conditional bias, selected structural
intervals, a separate joint-interaction model, and balanced observed-score G/D
studies. Numerical eligibility checks and interpretation limits remain in
force. Availability on PyPI is not statistical validation or certification.

## Installation

For a published pre-release:

```bash
python -m pip install --pre mfrmpy
```

For Matplotlib figures and offline HTML reports:

```bash
python -m pip install --pre 'mfrmpy[figures]'
```

During release preparation, install the built wheel instead:

```bash
python -m pip install dist/mfrmpy-0.1.0a1-py3-none-any.whl
```

Python 3.11 and 3.12 are the initial validation targets. Core dependencies are
NumPy, pandas and SciPy. This distribution provides a library, not the
standalone desktop executable or the Web app launcher. Those applications
remain available from the
[application repository](https://github.com/Ryuya-dot-com/MFRM_Python_Application).

## First analysis

Ratings must be in long form: one row per observed rating, with separate
person, facet and score columns. Specify the full intended integer rating
range, even when an endpoint is absent from the observations. Missing ratings
are not zero scores. The bundled example is synthetic.

```python
import mfrmpy

data = mfrmpy.load_example()
analysis = mfrmpy.analyze(
    data,
    person_col="Person",
    facet_cols=["Rater", "Task", "Criterion"],
    score_col="Score",
    rating_min=0,
    rating_max=5,
    model="RSM",
    method="JMLE",
)

for title, frame, note in mfrmpy.tables(analysis):
    print(title)
    print(frame.head())
    print(note)

mfrmpy.save_project("writing.mfrm", analysis)
restored = mfrmpy.load_project("writing.mfrm")  # no refitting
```

Use `pandas.read_csv("ratings.csv", dtype={"Person": str, "Rater": str})`
for your own data, adjusting the ID columns as needed. Explicit string types
preserve leading zeroes. The input/design audit remains part of the returned
analysis and tables; inspect it before interpreting estimates.

## Models and supplementary requests

`analyze` is the shared application service, not a second implementation of
the estimators. It returns a retained-analysis dictionary and does not display
an interface. Its MFRM settings are:

| Setting | Meaning / default |
|---|---|
| `person_col`, `facet_cols`, `score_col` | Required column names; `facet_cols` is a list |
| `rating_min`, `rating_max` | Required integer rubric endpoints, minimum < maximum |
| `model`, `method` | `"RSM"` / `"PCM"`; `"JMLE"` / `"MML"`; defaults RSM and JMLE |
| `step_facet` | Required for PCM; one selected facet with separate step structures |
| `maxit`, `reltol` | Maximum iterations (400) and stopping tolerance (1e-6) |
| `quad_points`, `population_prior_sd` | Fixed-SD MML quadrature count (15) and normal population SD (1.0) |
| `anchors`, `group_anchors` | Fixed-level / group-mean anchor records under the existing identification rules |
| `compute_pca` | Optional residual diagnostics; default `False` |
| `umean`, `uscale` | Reporting origin (0) and positive scale multiplier (1); canonical logits are retained |
| `progress` | Optional callback receiving phase names, not percentages |

Keep `keep_original=True`, `estimate_population_sd=False` and `mml_engine="EM"`.
Other choices, bounded GPCM and latent regression are outside this initial
public interface. A numerical termination flag alone is not evidence that
inference is valid; the returned output qualifications also matter.

Example requests, supplied **when fitting**, are:

```python
analysis = mfrmpy.analyze(
    data,
    person_col="Person",
    facet_cols=["Rater", "Task", "Criterion"],
    score_col="Score",
    rating_min=0,
    rating_max=5,
    model="PCM",
    step_facet="Criterion",
    method="MML",
    population_prior_sd=1.0,
    bias_requests=[{"facet_a": "Rater", "facet_b": "Task"}],
    interval_requests=[{"facet": "Rater", "difference": ["R1", "R2"], "level": 0.95}],
    joint_settings={"facet_a": "Rater", "facet_b": "Task"},
)
```

Replace the contrast labels with levels that actually occur in your data.
Conditional bias holds the main calibration fixed and does not establish
statistical significance or fairness. Joint interaction fits a **separate**
fixed-N(0,1), unit-weight model and retains its result separately. Selected
structural intervals require fixed-SD, unit-weight RSM/PCM and suitable
identification, termination, boundaries, information and quadrature checks.
They are pointwise asymptotic intervals, with no finite-sample coverage
guarantee. Holm adjustment of joint tests does not make intervals simultaneous.

Unavailable/failed interval or joint requests remain recorded with a reason;
they do not silently erase a completed baseline fit. For details, see the
[analysis scope](https://github.com/Ryuya-dot-com/MFRM_Python_Application/blob/main/docs/analysis_scope.md).

## G study and D study without an MFRM fit

These analyze observed scores in a complete, balanced design with one or two
random facets, at least two levels of each factor and one observation per cell.
The score need not be an integer for G/D analysis. For example:

```python
import pandas as pd

ratings = pd.DataFrame([
    {"Person": str(p), "Rater": str(r), "Task": str(t),
     "Score": 0.4*p + 0.3*r + 0.2*t + 0.1*((p+r+t) % 3)}
    for p in range(8) for r in range(3) for t in range(2)
])
gd_analysis = mfrmpy.analyze(
    ratings,
    workflow="gd",
    gd_settings={
        "object_facet": "Person",
        "random_facets": ["Rater", "Task"],
        "score_col": "Score",
        "target_population": "Comparable persons, raters and tasks",
        "selected_plan": {"Rater": 4, "Task": 2},
    },
)
for title, frame, note in mfrmpy.tables(gd_analysis):
    print(title, frame, note, sep="\n")
```

Unsupported designs return an unavailable result with a reason. G is for
relative decisions and Phi for absolute decisions. These are not MFRM person
reliability. G/Phi intervals and multivariate G-theory are unavailable. A fixed
conditioning facet can be specified through `gd_settings["fixed_condition"]`
under the existing desktop rules; it does not add another random facet.

## Results, plots and saved projects

The stable public surface is the names in `mfrmpy.__all__`. The bundled
`mfrm_app` modules are implementation details, including experimental methods
outside the supported scope. Do not rely on their private functions or treat
their presence as a supported public API.

| Function | Result |
|---|---|
| `analyze(data, **settings)` | Retained analysis dictionary, including inputs, settings, identities and qualification |
| `tables(analysis)` | List of `(title, pandas.DataFrame, note)` tuples, with current output restrictions applied |
| `load_example()` | A new DataFrame containing the bundled synthetic writing example |
| `save_project(path, analysis)` | Atomically write a validated `.mfrm` archive |
| `load_project(path)` | Validate and restore an archive without refitting; malformed archives raise `ProjectError` |
| `plot(analysis, kind="wright", **options)` | Matplotlib Figure; requires `figures` extra |
| `report_html(analysis, **options)` | Standalone HTML string; requires `figures` extra |

`analysis["result"]` and `analysis["diagnostics"]` retain raw MFRM output;
`analysis["gd_analysis"]` retains G/D results. Treat the nested implementation
schema as provisional in this pre-release. For reporting, use `tables`, keeping
the notes and `Use`/`Reason` information. A blank interval or SE is unavailable,
not zero. General MML tables, person-score uncertainty and automatic equivalence
claims retain their existing inferential restrictions.

With the optional figures dependency installed:

```python
from pathlib import Path

step_facet = analysis["fit_settings"].get("step_facet")
level = str(data[step_facet].iloc[0]) if step_facet else None
figure = mfrmpy.plot(analysis, "wright", level=level)
figure.savefig("wright.png", dpi=300)
Path("report.html").write_text(mfrmpy.report_html(analysis), encoding="utf-8")
```

Other figure kinds include `"pathway"`, `"category"`, `"bias"`, `"interval"`,
`"joint"`, `"gd"` and `"pca"` when the required results are present. PCM Wright
maps and category curves require `level` (a step-facet level); supplementary requests use their
zero-based `index`. Retained results are reused; plotting and reopening must
not recompute estimation. Unsupported requests can produce a reason-bearing
figure rather than a numerical result.

Saved `.mfrm` projects use the existing desktop/Web archive format. New source
or dependency versions can cause current qualification checks to withhold old
interval output without altering the archived estimates. Record
`mfrmpy.__version__` and `analysis["source"]` when reporting an analysis.
The package API version and the shared engine version are distinct.

## Privacy, license and citation

Calculations run locally. Saved projects and exported reports can retain
ratings and person/rater IDs; store and share them according to your research
requirements. Use synthetic data for bug reports.

MIT licensed; see
[LICENSE](https://github.com/Ryuya-dot-com/MFRM_Python_Application/blob/main/LICENSE)
and the
[responsible-use notice](https://github.com/Ryuya-dot-com/MFRM_Python_Application/blob/main/LICENSE_NOTICE.md).
Cite Ryuya Komuro, *mfrmpy*, the installed package version, and
https://github.com/Ryuya-dot-com/MFRM_Python_Application, together with the
statistical methods actually used. Source and issue reporting are maintained
in that repository.
