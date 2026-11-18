# mfrmpy: Many-facet Rasch analysis in Python

`mfrmpy` brings the supported MFRM Desktop calculations to Python scripts and
notebooks. It uses the same Python engine as the existing MFRM Web and desktop
applications. It does not require R, FACETS, Streamlit or Tkinter for analysis.

Version **0.2.0** refines output qualification, four-series precision reporting,
low-d.f. and boundary fit safeguards, and CSV/archive provenance while retaining
full read compatibility with 0.1.0 `.mfrm` projects. Available analyses are RSM/PCM with JMLE
or fixed-population-SD MML, diagnostics, conditional bias, selected structural
intervals, a separate joint-interaction model, and balanced observed-score G/D
studies. Numerical eligibility checks and interpretation limits remain in
force. Availability on PyPI is not statistical validation or certification.

## Installation

Install from PyPI:

```bash
python -m pip install mfrmpy
```

For Matplotlib figures and offline HTML reports:

```bash
python -m pip install 'mfrmpy[figures]'
```

To upgrade an earlier installation, including an alpha release:

```bash
python -m pip install --upgrade mfrmpy
```

Python 3.11 and 3.12 are the initial validation targets. Core dependencies are
NumPy >=1.26, pandas >=3.0, SciPy >=1.10 and PyArrow >=16.0. pandas 3 and
PyArrow preserve the string storage and missing-value types recorded in Web
and desktop `.mfrm` projects, including their content fingerprints. This
distribution provides a library, not the
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
preflight = mfrmpy.check_design(
    data,
    person="Person",
    facets=["Rater", "Task", "Criterion"],
    score="Score",
    rating_min=0,
    rating_max=5,
)
print(preflight["summary"])

analysis = mfrmpy.analyze(
    data,
    person="Person",
    facets=["Rater", "Task", "Criterion"],
    score="Score",
    rating_min=0,
    rating_max=5,
    model="RSM",
    method="JMLE",
)

print(mfrmpy.summary(analysis))
for title, frame, note in mfrmpy.tables(analysis, kind="measures"):
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
| `person_col` (`person`), `facet_cols` (`facets`), `score_col` (`score`) | Required column names; `facet_cols` / `facets` is a list |
| `rating_min`, `rating_max` | Required integer rubric endpoints (`int` or `numpy.integer`), minimum < maximum |
| `model`, `method` | `"RSM"` / `"PCM"`; `"JMLE"` / `"MML"`; defaults RSM and JMLE |
| `step_facet` | Required for PCM; one selected facet with separate step structures |
| `maxit`, `reltol` | Maximum iterations (400) and stopping tolerance (1e-6) |
| `quad_points`, `population_prior_sd` | Fixed-SD MML quadrature count (default `31` for MML) and normal population SD (1.0) |
| `noncenter_facet`, `positive_facets` | Non-centered reference facet (default `"Person"`) and optional list of facets given additive `+1` leniency/easiness orientation |
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
| `check_design(data, **settings)` | Optimization-free preflight dict checking MFRM connectivity / category support or G/D balance |
| `analyze(data, **settings)` | Retained analysis dictionary, including inputs, settings, identities and qualification |
| `summary(analysis)` | One-screen DataFrame of design counts, qualification state, hold reasons and available table/plot kinds |
| `tables(analysis, kind=None)` | List of `(title, pandas.DataFrame, note)` tuples (optionally filtered by `kind`, e.g. `"measures"`, `"steps"`/`"thresholds"`, `"reliability"`, `"agreement"`, `"bias"`, `"intervals"`, `"joint"`, `"gd"`) |
| `export_csv(analysis, directory)` | Atomically write qualified CSV tables and `00_manifest.csv` (with current and saved reasons) |
| `load_example()` | A new DataFrame containing the bundled synthetic writing example |
| `save_project(path, analysis)` | Atomically write a validated `.mfrm` archive |
| `load_project(path)` | Validate and restore an archive without refitting; malformed archives raise `ProjectError` |
| `plot(analysis, kind="wright", **options)` | Matplotlib Figure; requires `figures` extra |
| `report_html(analysis, **options)` | Standalone HTML string (`facet`/`statistic` or `pathway_facet`/`pathway_statistic`); requires `figures` extra |

`analysis["result"]` and `analysis["diagnostics"]` retain raw MFRM output;
`analysis["gd_analysis"]` retains G/D results. The nested implementation
schema is internal; use `.mfrm` files to retain analyses across sessions rather
than serializing or modifying this dictionary yourself. For reporting, use
`tables` or `export_csv`, keeping the notes and `SEUse`/`CIUse`/`Use`/`StatusReason`
fields:

- **Distinct run qualification states:** `stat.jmle.identifiability_audit_size_limit`
  (exact dense rank audit skipped due to matrix size limit), `stat.jmle.eta_rank_deficient`
  (proven structural rank deficiency / disconnected subsets), and
  `stat.estimation.optimizer_termination_unverified` (non-convergence) are reported
  separately rather than merged under one hold code.
- **Exploratory SEs and intervals (`screening_only` / `SupportsFormalInference = False`):**
  Converged JMLE (including size-limited dense audits without proven disconnection)
  and converged fixed-SD MML display finite `ModelSE`, `RealSE = ModelSE * sqrt(max(Infit, 1))`
  and exploratory normal bands with explicit `Basis`, `SEUse = "screening_only"`,
  `CIUse = "screening_only"` and `SupportsFormalInference = False`. Non-converged,
  structurally rank-deficient, free-SD MML, fixed-constraint and extreme-endpoint
  rows withhold SEs and intervals (`""`, not zero).
- **Four-series reliability and counts:** `Separation and Reliability by Facet`
  reports all four combinations of `SEMode` (`model` vs `fit_adjusted`) and
  `DistributionBasis` (`sample` ddof=1 vs `population` ddof=0), alongside `Levels`,
  `EstimateAvailable`, `SEAvailable`, `ExcludedLevels`, `AggregationTarget`,
  `Use` and `SupportsFormalInference = False`.
- **Fit and boundary safeguards (`FitStatus` / `FitReason` / `ObservedCategories`):**
  Following `mfrmr` 0.3.0, ZSTD is withheld (`""`) when degrees of freedom are
  below 1 (preventing Wilson–Hilferty sign inversion) or when a scope includes
  JMLE extreme-score boundary proxies; single-observation free JMLE levels
  withhold zero-d.f. mean squares (`unavailable_zero_df`); unbounded endpoint
  facet levels, zero-count internal/boundary score categories (`stat.category.*`),
  and disconnected MML cohorts (`stat.mml.population_assumption_linked`) withhold
  unqualified estimates, SEs, specified intervals, joint interactions and bias;
  and `ObservedCategories`, agreement notes and G/D metadata flag constant-score
  levels explicitly.
- **Shared re-display and saved provenance:** `tables`, `plot`, `report_html` and
  `export_csv` apply the same qualification rules to fresh and reopened `.mfrm`
  projects while retaining original archived values and save-time reasons (`Saved*`
  columns, `00_manifest.csv` and HTML archive notes).
- **Cross-package numerical alignment (`mfrmr` 0.3.0 / `TAM` / `ConQuest`):** Fixed-SD
  MML defaults to `quad_points=31` (matching `mfrmr` 0.3.0 so standard designs satisfy
  the `q` vs `2q-1` sensitivity gate for `Specified intervals` and `Joint interaction`,
  while earlier `.mfrm` archives retain their stored `quad_points`), and `analyze` /
  `check_design` accept `noncenter_facet` (e.g. `noncenter_facet="Task"` to match
  `TAM` `constraint="cases"` and `ConQuest` `lconstraints=cases`) and
  `positive_facets`. `Step Thresholds` (`tables(analysis, "steps")`) reports centered
  adjacent-category thresholds, observed transition counts, disorder flags, and
  exploratory `screening_only` SEs/bands (full observed-information delta-method SEs
  including the sum-to-zero final step for MML; analytical conditional step-information
  SEs for JMLE). `Specified intervals` includes explicit `SupportsFormalInference = False`
  and `Reason` columns, and conditional `Bias` (`tables(analysis, "bias")`) reports
  `ObsScore`, `ExpScore`, `ObsExpAvg`, `BiasSize`, `SE`, `Z`, `LR_ChiSq`, and 1-d.f.
  conditional likelihood-ratio `LR_Prob` under `Use = "screening_only"`
  (`SupportsFormalInference = False`), withholding uncertainty under calibration or
  cell holds. Fixed-SD MML facet SEs propagate full joint observed-information
  covariance `sqrt(diag(J I^{-1} J^T))` across all levels including constrained
  sum-to-zero levels (matching `mfrmr` 0.3.0, whereas `TAM::tam.se()` and `ConQuest`
  `stderr=quick` compute diagonal curvature `1 / sqrt(diag(I))`, and `ConQuest`
  5.47.5 `stderr=empirical` matches on free coordinates but repeats the preceding SE
  on `*`-constrained levels); and MML Person reliability uses empirical EAP variance
  `Var(EAP) / (Var(EAP) + mean(SD^2))` (matching `TAM::EAP.rel` under
  `DistributionBasis = "population"`).

With the optional figures dependency installed:

```python
from pathlib import Path

step_facet = analysis["fit_settings"].get("step_facet")
level = str(data[step_facet].iloc[0]) if step_facet else None
figure = mfrmpy.plot(analysis, "wright", level=level)
figure.savefig("wright.png", dpi=300)
Path("report.html").write_text(mfrmpy.report_html(analysis), encoding="utf-8")
```

Other figure kinds include:
- Core maps & curves: `"pathway"` (alias `"fit_pathway"`; `statistic="Infit"|"Outfit"|"InfitZSTD"|"OutfitZSTD"`), `"category"` (alias `"ccc"`), `"expected"` (alias `"tcc"`), and `"information"` (category Fisher information and conditional SE).
- Facet & rubric profiles: `"facet"`, `"rater_severity"`, `"step"` (alias `"thresholds"`), `"bubble"`, and `"fair_average"` (`statistic="FairM"|"FairZ"`).
- Residual & supplementary diagnostics: `"unexpected"`, `"residual_qq"`, `"pca"` (alias `"residual_pca"`), `"pca_loadings"`, `"bias"`, `"bias_ranked"`, `"interval"`, `"joint"`, and `"gd"`.

PCM Wright maps and category curves require `level` (a step-facet level), while `"step"`, `"expected"`, and `"information"` accept `level` for a single step-facet level or `level=None` to overlay all step-facet levels; supplementary requests use their
zero-based `index`. Retained results are reused; plotting and reopening must
not recompute estimation. Unsupported requests can produce a reason-bearing
figure rather than a numerical result.

Saved `.mfrm` projects use the existing desktop/Web archive schema with
deterministic outer gzip compression (`\x1f\x8b`) by default; `load_project`
automatically reads both gzip-compressed and earlier uncompressed `.mfrm`
archives without changing `content_fingerprint` (up to 128 MiB on disk and
512 MiB uncompressed). Older `0.1.0` releases read uncompressed archives only.
New source or dependency versions can cause current qualification checks to
withhold old interval output without altering the archived estimates. Record
`mfrmpy.__version__` and `analysis["source"]` when reporting an analysis.
The package API version and the shared engine version are distinct.

## API and compatibility policy

The supported API is the exported functions documented above, `ProjectError`
and `__version__`. `analyze` accepts the settings and workflows documented in
this guide. `ProjectError` is a `ValueError` subclass for invalid/unsupported
project contents; filesystem errors such as a missing file remain normal
Python I/O errors.

Within **0.1.x**, patch releases preserve these public names, documented
arguments and return types. Additive options and bug fixes are allowed. Before
1.0, an incompatible public API change requires a new minor release, such as
0.2.0, with migration instructions in the changelog. Private `mfrm_app` helpers,
nested raw-result fields, table display order and exact wording are outside
this compatibility promise. Keep interpretation notes with exported results.

Version 0.1.0 reads the supported schema-1 RSM/JMLE and schema-2 `.mfrm` archives,
including those written by 0.1.0a2 and the corresponding Web/desktop workflows.
The 0.1.x series will retain read compatibility with these supported archives;
malformed or unsupported records are rejected without refitting. Future archive
changes will be documented, and older applications may not read newer formats.

Numerical bug fixes or updated qualification rules can change reported results
or withhold an output; such changes belong in the changelog. For reproducible
research, retain the original `.mfrm` file and record the package, source and
dependency versions. The first stable release uses the same calculation and
archive code as 0.1.0a2, with the same analysis scope and interpretation limits.

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
