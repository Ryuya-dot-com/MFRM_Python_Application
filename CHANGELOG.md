# Changelog

Changes to the Python library, Web and desktop applications are summarized by version below.
Historical entries describe changes at that time, not current availability or
validation of every output. Some advanced or uncertainty features were later
restricted or removed from the core workflow. Consult
[analysis scope](docs/analysis_scope.md) and the current output's notes before
using a result. Complete source changes remain in Git history; scientific
methods and validation findings belong in their dedicated technical records.

## mfrmpy 0.2.0 — 2026-11-18

- Separate JMLE size-limited dense rank audits (`stat.jmle.identifiability_audit_size_limit`),
  proven structural rank deficiency (`stat.jmle.eta_rank_deficient`) and
  unverified optimizer termination (`stat.estimation.optimizer_termination_unverified`)
  across qualification snapshots, tables, figure cautions, HTML reports and CSV
  manifests.
- Surface finite conditional and fixed-SD MML `ModelSE`, fit-adjusted `RealSE`
  and exploratory normal bands with explicit `Basis`, `SEUse = "screening_only"`,
  `CIUse = "screening_only"` and `SupportsFormalInference = False`, while keeping
  non-converged, rank-deficient, free-SD, fixed-constraint and extreme-endpoint
  cells withheld.
- Expand `Separation and Reliability by Facet` to the four `SEMode` (`model` /
  `fit_adjusted`) × `DistributionBasis` (`sample` / `population`) series with
  `Levels`, `EstimateAvailable`, `SEAvailable`, `ExcludedLevels`,
  `AggregationTarget` and `SupportsFormalInference = False`.
- Adopt `mfrmr` 0.3.0 fit safeguards (`df >= 1` ZSTD guard on both engine and
  FACETS conventions, JMLE extreme-person boundary-proxy ZSTD guard, and complete
  row variance/residual checks) with explicit `FitStatus` and `FitReason` columns.
- Add `mfrmpy.export_csv(analysis, directory)` and apply one shared qualification
  contract across `tables`, `plot`, `report_html`, `export_csv` and `.mfrm`
  reload, preserving raw archived values and save-time reasons (`Saved*` columns
  and `qualification_at_save`) without refitting or altering archive fingerprints.
- Compress newly written `.mfrm` archives with deterministic outer gzip (`\x1f\x8b`,
  single-pass serialization, ~6x smaller output) while reading both gzip-compressed
  and earlier uncompressed `.mfrm` archives with bounded streaming decompression
  guards, and generate Web `.mfrm` downloads on demand.
- Expand `mfrmpy.plot` and `.mfrm` report/viewer figures with zero-refit `mfrmr`
  diagnostic plots (`"facet"`, `"rater_severity"`, `"step"`/`"thresholds"`,
  `"bubble"`, `"expected"`/`"tcc"`, `"information"`, `"fair_average"`,
  `"unexpected"`, `"residual_qq"`, `"pca_loadings"`, `"bias_ranked"`, plus
  `"Outfit"`/`"OutfitZSTD"` on `"pathway"`), preserving qualification holds and
  `presentation="explained"|"manuscript"` caption separation.
- Harden boundary and degenerate-input safeguards in alignment with `mfrmr` 0.3.0:
  gate zero-count internal/boundary score categories and single-category PCM
  step-facet levels (`stat.category.*`), unpenalized MML endpoint recession
  levels (`stat.mml.extreme_free_facet_location`), disconnected MML rating
  cohorts across general SEs, specified intervals, joint interactions and bias
  (`stat.mml.population_assumption_linked`), single-observation free JMLE fit
  mean squares (`unavailable_zero_df`), `ObservedCategories` / G/D and agreement
  zero-variance notes, and `SavedReliability` provenance hygiene.
- Consolidate the public `mfrmpy` API and shared internal renderers: add
  optimization-free `mfrmpy.check_design`, one-screen `mfrmpy.summary`, and
  `mfrmpy.tables(analysis, kind=...)` selectors; accept `numpy.integer` /
  `numpy.bool_` scalars, `person`/`facets`/`score` column aliases on `analyze`
  and `check_design`, and `facet`/`statistic` aliases on `report_html`; and
  unify table/category-curve helpers across `basic_analysis`, `desktop_analysis`,
  `visual_maps` and `desktop_outputs`.
- Align defaults, step/bias tables, and scale-orientation options with 4-engine
  verification (`mfrmpy` / `mfrmr` 0.3.0 / `TAM` 4.3.25 / `ConQuest` 5.47.5):
  raise the default fixed-SD MML quadrature order in `analyze` to `quad_points=31`
  (retaining archived `quad_points` in reopened `.mfrm` projects) and surface
  `SupportsFormalInference = False` and `Reason` in `Specified intervals`; add the
  zero-refit `Step Thresholds` table (`kind="steps"` / `"thresholds"`) with MML
  full observed-information delta-method SEs, JMLE conditional step-information SEs,
  transition counts and disorder flags; expose `noncenter_facet` and `positive_facets`
  across `analyze`, `check_design`, `summary` and `.mfrm` validation/portability; and
  enrich the `mfrmpy` conditional `Bias` screening table (`ObsScore`, `ExpScore`,
  `ObsExpAvg`, `BiasSize`, `SE`, `Z`, `LR_ChiSq`, `LR_Prob`, `Use="screening_only"`,
  `SupportsFormalInference=False`) while preserving the Web APA publication contract.

## mfrmpy 0.1.0 — 2026-11-15

- Release the first stable Python library; normal installation uses
  `python -m pip install mfrmpy` without `--pre`. Existing installations can
  update with `python -m pip install --upgrade mfrmpy`.
- Define the supported public API and patch-release compatibility within
  0.1.x, including saved-project read compatibility. Document how incompatible
  API changes and numerical corrections will be communicated.
- Retain the 0.1.0a2 calculations, dependencies, public call signatures and
  archive format. Earlier alpha and hosted Web projects keep their analysis
  identities, qualified tables and notes when restored without refitting.

## mfrmpy 0.1.0a2 — 2026-11-15

- Include PyArrow >=16.0 and require pandas >=3.0 / NumPy >=1.26 so a normal
  installation can restore the exact string storage and missing-value types
  saved by the hosted Web app. The first package's numeric-only dependencies
  could not restore Arrow-backed projects; pandas 2 also lacks the required
  string missing-value constructor argument.
- Add installed-package checks for Python/Arrow string storage with both
  `pd.NA` and NaN, preserving leading-zero identifiers, literal `"nan"` text,
  saved fingerprints and zero-refit restoration. Calculations and the archive
  format are unchanged.

## mfrmpy 0.1.0a1 — 2026-11-15

- Publish [`mfrmpy` 0.1.0a1](https://pypi.org/project/mfrmpy/0.1.0a1/) as an initial
  pre-release Python library sharing the existing
  desktop/Web calculation engine. Add a documented public interface, bundled
  synthetic data, qualified tables, .mfrm save/reload, and optional figures/HTML.
- Package shared reporting text without depending on a repository checkout.
  Numerical algorithms and statistical interpretation boundaries are unchanged.
- Add wheel/source-distribution validation and installed-package CI, with
  publication available only through an explicit manual request on main.
- Update the Web workflow control budget for the existing saved-project entry.
- Check Web archive interoperability through completion instead of treating
  a 45-second AppTest limit as a functional failure. Statistical assertions
  and the application runtime remain unchanged.

## Unified macOS & Windows Desktop preview — 2026-10-31

The [macOS (Apple Silicon) and Windows 11 (x64) preview release](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-20261031)
hosts both standalone ZIP archives directly on GitHub Releases (retaining displayed app version `0.2.19-beta`):

- Embed custom MFRM application icons (`.icns` / `.ico`), Windows `VSVersionInfo` metadata, uncompressed binaries (`upx=False`), and GitHub Actions Sigstore SLSA build provenance on the Windows archive.
- Vectorize core `mfrm_app/engine.py` sufficient-statistic and Hessian calculations with exact numerical equivalence, and unify univariate `G/D-Study` computation under `mfrm_app/generalizability.py` (closed-form arbitrary integer D-Study projections, centered-tensor orthogonal ANOVA residual sums of squares, fixed-facet conditioning for balanced ≥3-facet data, preflight balance checks, and streamlined bilingual result notes).

## Windows Desktop preview — 2026-10-30

The [Windows 11 x64 preview](https://github.com/Ryuya-dot-com/MFRM_Python_Application/releases/tag/desktop-windows-20261030)
includes Python and retains the displayed app version 0.2.19-beta. The dated
desktop tag identifies this distribution separately from the earlier Web release.

- Integrate the complete initial desktop workflow: RSM/PCM with JMLE or fixed-SD
  MML, Wright/Pathway/category figures, conditional bias, specified intervals,
  joint interactions and supported G/D studies in one application.
- Add English/Japanese controls, explanatory/manuscript figures with separate
  captions, and HTML/CSV/PNG/SVG exports. Save extended projects while retaining
  read support for the older basic format. Older applications cannot read the
  extended format; preserve original projects during updates.
- Correct Windows worker completion and failure handling, including pipe EOF
  (`WinError 109`), cancellation and retry. Preserve results until worker exit
  status is known.
- Preserve observation indices in score-probability calculations after score
  rows are excluded, including gapped, permuted, duplicate and text indices.
- Correct FACETS score-file parsing: keep `PtMea`, `PtBis` and `PtBisI` distinct.
  The application's existing PTMEA definition remains different from FACETS
  PtMea; see the [bounded comparison](validation/FACETS_WINDOWS_DESKTOP_REVIEW_20261029.md).
- Correct pandas reporting compatibility without changing numeric values or
  unavailable-output handling. Isolate Streamlit AppTest state so later spawned
  desktop workers use the correct process entrypoint.
- Provide bilingual desktop guides and a self-contained ZIP startup guide;
  distinguish source archives from application downloads, and explain saved
  project compatibility and output-specific interpretation limits.
- Reorganize the public roadmap, release checklist and deployment guidance by
  their reader and purpose. Keep release changes separate from development
  diaries and statistical study records.

The Windows application is an unsigned preview. Its release lists download
details, tested environments and limitations. The macOS package is not yet
publicly released.
Existing inferential restrictions and supported-design conditions remain in
force. A published Web/source version does not establish desktop acceptance.

## 0.2.19-beta — 2026-10-28

This release updates the Web application and Python source; it does not include
standalone desktop binaries.

- Add specified intervals for selected fixed-SD RSM/PCM facet levels and
  prespecified differences, using full structural covariance and target-specific
  termination, identification, boundary, information and quadrature checks.
- Save the same interval record in tables, figures, CSV/JSON/HTML and
  Word/PDF/HTML reports. Reload matching interval records without refitting;
  reject incompatible, stale or inconsistent records.
- Retain intervals computed later in opted-in history for the exact analysis.
  Preserve metadata-only history and access to saved results after inputs change.
- Improve bilingual input, analysis, report and restoration guidance. Preserve
  non-default optional-analysis settings, uploaded inputs and cache identity.
- Use basic reports for the main results and detailed reports for methods and
  writing supplements. Improve table titles, units, notes, pagination and CSV
  selection, and include result-aware reporting guidance.

<a id="unreleased--directory-review-corrections"></a>

### Output and interpretation corrections

- Carry failed-fit, person-boundary and calibration restrictions through live,
  cached and reopened tables, figures, reliability, scoring and report prose.
  Preserve valid structural results without promoting unsupported person values.
- Restrict general MML uncertainty, free-population-SD uncertainty, information
  criteria/model recommendations, conditional bias significance and automatic
  equivalence conclusions to their documented scope. Close report and stale-
  result paths that could bypass those restrictions.
- Preserve the estimator's actual stopping settings and distinguish optimization
  error from quadrature sensitivity. Offer explicit finer-quadrature refitting
  without silently changing the saved original analysis.
- Require comparable data, likelihood bases, constraints and population settings
  when retaining model comparisons. Numerical agreement alone does not establish
  empirical model validity or estimator superiority.
- Restrict G/D ANOVA to supported complete, balanced, unit-weight one-/two-random-
  facet designs; explain unsupported designs and distinguish signed variance
  components from the values used in coefficients.
- Explain Strata, fit statistics, scale direction, fixed anchoring versus joint
  calibration, and interaction diagnostics without equating them with actual
  groups, exclusion decisions or assessment fairness.
- Correct category-curve fallback rendering and preserve interpretation notes
  when saving plots and reports. Hosted static PNG generation depends on the
  renderer; interactive HTML remains the offered fallback when unavailable.

Specified intervals are pointwise model-based approximations, not a finite-
sample coverage guarantee. Their availability does not remove restrictions on
Person-score uncertainty, free-SD inference, equivalence or other held outputs.

## 0.2.18-beta — 2026-09-15

- Improve Japanese guidance across setup, results, diagnostics and exports.
- Add a Home action that preserves data, settings and results.
- Add current-fit comparison exports and optional offline PCM reference scripts.
- Add a read-only fixed-SD PCM returned-point check that distinguishes optimizer
  stationarity from quadrature sensitivity and does not release inference holds.

## 0.2.17-beta — 2026-09-13

- Consolidate the Web application with simpler onboarding, English/Japanese
  guidance, preserved workspace context and clearer report/export routes.
- Keep display-density controls accessible and repair citation preparation.
- Preserve output protections, optimizer-coordinate identification and relevant
  numerical checks. Development-only studies do not establish broad inference.
- Make normal verification work from a clean checkout using synthetic fixtures;
  separately retained study evidence remains outside that test selection.

## 0.2.15-beta - 2026-06-17

- Expand rating-design, G/D planning, model-comparison, parameter-recovery,
  residual-PCA, bias and uncertainty-review outputs. Later interpretation
  restrictions supersede any earlier presentation of those outputs as ready.
- Add claim-to-evidence and manuscript/figure bundles, reference exports,
  publication table styling and result-aware reporting guidance.
- Expand English/Japanese controls and spreadsheet guidance. Add Excel,
  Parquet and JSON inputs, column-role detection and wide-table conversion.
- Improve category information and GPCM fair-average calculations, including
  handling of nonfinite predictors. Retain their method-specific limitations.
- Clarify licensing, external-reference roles and the distinction between native
  computation and optional external code-generation/compatibility surfaces.

## 0.2.14-beta - 2026-04-18

- Add sparse music peer-rating and large-scale writing scenarios and improve
  scenario selection/readiness messages.
- Correct categorical aggregation that expanded unobserved combinations and
  crashed sparse-design agreement calculations.
- Permit singleton-facet input to proceed with an appropriate readiness message
  and avoid selecting singleton columns as default facets.

## 0.2.13-beta - 2026-04-18

- Add a binary reading/testlet-format example and guidance distinguishing fixed
  task effects from random testlet models; core binary score handling already
  existed.
- Prefer exact column-name matches so `Score` is not confused with `Scorer`;
  reject duplicate columns after role mapping.
- Guard singleton variance calculations and improve multi-author citation parsing.

## 0.2.12-beta - 2026-04-18

- Add complete sample-to-results UI checks and small-input boundary cases so
  valid generated data also exercises dependent widgets and output rendering.

## 0.2.11-beta - 2026-04-18

- Fix forest-plot, misfit-ranking and bias-heatmap controls whose bounds could
  crash when facets contained only a few levels.

## 0.2.10-beta - 2026-04-18

- Add Matplotlib fallbacks for publication figures when the Plotly static-image
  renderer is unavailable.
- Simplify one-click Word/PDF/HTML downloads and make figure failure handling
  explicit rather than silently omitting expected report figures.

## 0.2.9-beta - 2026-04-17

- Add ingestion-time outlier review and broaden checks for existing outputs.
- Reduce repetitive sidebar instructions while preserving readiness messages.

## 0.2.8-beta - 2026-04-17

- Simplify data-source selection, add scenario switching and display the active
  dataset's context without requiring a separate expander.

## 0.2.7-beta - 2026-04-17

- Add selectable synthetic assessment scenarios with descriptions, references
  and scenario-specific downloadable CSV names.

## 0.2.6-beta - 2026-04-17

- Expand contextual Help and apply essential/full display selection to Help
  and Report views without removing the detailed routes.
- Make replayable configuration fields explicit and preserve stable tab keys.

## 0.2.5-beta - 2026-04-17

- Add contextual Help popovers and an essential/full display-density control.

## 0.2.4-beta - 2026-04-17

- Add observation-coverage and category-usage plots, forest plots, residual Q–Q
  plots and measure ECDFs. Add posterior diagnostic plots to the historical
  optional viewer; their existence does not qualify the underlying inference.

## 0.2.3-beta - 2026-04-17

- Align readiness and reporting checklists with the outputs actually available.
- Repair outdated navigation references and explain why a diagnostic was skipped.
- Expand Help, troubleshooting and documentation of optional upload/download
  surfaces.

## 0.2.2-beta - 2026-04-17

- Add analysis context, glossary access and consistent measure-table ordering.
- Improve diagnostic reference-line visibility and preserve current analyses
  when visiting the historical posterior viewer.
- Correct column ordering for empty tables and improve generated-code downloads.

## 0.2.1-beta - 2026-04-17

- Add downloadable sample CSVs, upload-size checks, history-deletion confirmation
  and a read-only configuration preview.
- Connect chart-reading guides to plots and show readiness with text as well as
  color. Improve keyboard focus, reduced-motion behavior and narrow layouts.
- Improve output-format guidance, rounding and history memory use.

## 0.2.0-beta - 2026-04-17

- Add historical download-only advanced-model code generators and an optional
  posterior viewer; these were not native inference engines in the application.
- Add publication-document exports, chart guidance, readiness/history controls
  and more specific analysis-failure explanations.
- Extend PCA, bias and reliability diagnostics with explicit unavailable cases.

## 0.1.2-beta - 2026-04-12

- Add broader category-curve exports, synthetic demonstration reports and
  manuscript/interpretation guidance with explicit limitations.
- Add publication figure bundles and HTML fallback when PNG rendering is
  unavailable; align report readiness and displayed interpretation thresholds.
- Preserve uploaded streams during previews and invalidate predictions when
  uploaded content changes rather than relying only on filename and size.

## 0.1.1-beta - 2026-04-12

- Add local/hosted deployment guidance and a contract for optional comparisons
  with external software, distinguishing their targets and methods.

## 0.1.0-beta - 2026-04-12

- Introduce the standalone Streamlit application with RSM/PCM, bounded GPCM,
  JMLE/MML and advanced analysis workflows, subject to their documented limits.
- Include privacy guidance, reproducibility metadata, environment/self-test
  commands, synthetic examples, anchor guidance and MIT licensing.
- Do not claim exact numerical parity with FACETS, TAM, sirt, mirt or mfrmr.
