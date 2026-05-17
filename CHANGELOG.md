# Changelog

All notable changes to this standalone Streamlit distribution should be recorded here.

## Unreleased

### Added

- **Zotero-informed method reference audit**. A safe read-only audit of
  the local Zotero/BibTeX library was used to expand the in-app APA
  reference library with Rasch, fit, local-dependence, external-package,
  and recent MFRM extension anchors. New `build_method_reference_audit()`
  exports `method_reference_audit.csv`, mapping model, estimation,
  fit/person-fit, residual PCA, bias/local interaction, simulation,
  and external-validation surfaces to citation tokens, BibTeX/RIS export
  status, manuscript use, and claim boundaries. Current-run reference
  rows now feed the APA report, manuscript template, and
  `claim_to_evidence_matrix.csv`, so citations and statistical claim
  boundaries stay aligned with the analyses actually computed. Quick,
  download, demo bundles, and the manuscript binder include these tables.
- **APA sentence evidence audit**. Generated APA-style prose now has an
  exportable `apa_report_sentence_audit.csv` that maps each draft
  Method/Results sentence to the manuscript claim area, output evidence,
  evidence files, suggested citations, citation boundaries, and copy
  decision. The audit is visible in the APA Report and Claim Guide areas
  and is included in quick/download/demo/manuscript-binder exports.
- **Final-result manuscript handoff export**. Downloads now include
  `mfrm_manuscript_handoff.md` and
  `manuscript_handoff_checklist.csv`, a concise file-opening,
  download, submission-preflight, privacy, and archiving guide that
  links the Publication Document, table bundle, figure bundle, method
  appendix, manuscript template, config JSON, app-engine runner, and
  OSF package. The same assets are included in the OSF ZIP and
  synthetic demo report export.
- **Claim-to-evidence matrix and manuscript binder**. The Downloads
  tab and demo export now include `claim_to_evidence_matrix.csv`,
  which maps each manuscript area to claim status, gate status,
  exported tables, figures, readiness checks, reviewer questions,
  caveats, and archive files. A curated
  `MFRM_Manuscript_Binder.zip` / `MFRM_Demo_Manuscript_Binder.zip`
  collects the handoff, evidence matrix, action/gate tables, method
  appendix, manuscript template, config, and figure manifest for
  coauthor review and submission preparation.
- **Visual evidence binder**. Figure downloads now include
  `MFRM_Visual_Evidence_Binder.zip`, with the actual PNG/HTML figure
  files, `visual_evidence_map.csv`, `visual_caption_drafts.md`,
  `figure_manifest.csv`, visual interpretation and method-evidence
  tables, and the claim-to-evidence matrix. The demo report exports
  `MFRM_Demo_Visual_Evidence_Binder.zip` and
  `visual_evidence_map.csv` as well.
- **SE/CI coverage diagnostic for ADEMP recovery**.
  `evaluate_parameter_recovery(...)` now preserves `SE_Method`,
  `SE_Status`, `CI_Method`, `CI_Status`, and `UncertaintyCaution` in the
  recovery long table, and attaches a separate `coverage_report`.
  `build_se_ci_coverage_report(...)` and `evaluate_se_ci_coverage(...)`
  expose nominal coverage, empirical Coverage95, coverage error,
  SE/CI status summaries, interpretation, and recommended action for
  manuscript caveats. The Report tab displays this diagnostic under the
  parameter-recovery panel and offers `mfrm_se_ci_coverage_report.csv`.
- **Executable MML fixed-prior-SD sensitivity**.
  `evaluate_mml_prior_sd_sensitivity(...)` now refits the current MML
  result across selected `population_prior_sd` values, reusing the
  current result for the base SD and returning `summary`,
  `measure_deltas`, optional `population_deltas`, and an interpretation
  `report`. The Fit Details tab exposes this as an on-demand screen with
  a ZIP download so final reports can document whether person EAPs,
  non-person facet measures, ranks, and latent-regression coefficients
  are stable to the fixed population variance setting.
- **Conditional bias-inference audit**.
  `build_bias_inference_audit(...)` now summarizes every exported
  bias/local-interaction facet pair with cells screened, DFF flags,
  strong-review cells, sparse cells, Holm/BH counts, practical-logit
  counts, profile-CI status summaries, inference-tier summaries, and
  connected-scale status. DFF tables now preserve the conditional
  likelihood/SE/multiplicity basis from `estimate_bias_interaction(...)`
  and add cell-level claim status plus reporting cautions. The
  first-read guide, final-readiness gate, manuscript claim guide,
  APA-style narrative, quick result bundle, Downloads ZIP, and demo
  exports now use this audit instead of raw `|t| >= 2` counts.
- **Residual PCA stability sensitivity audit**.
  `compute_pca_bundle(...)` now augments the basic residual-matrix
  checks with leave-one-column-out EV1/share/loading stability and a
  deterministic row-bootstrap audit of EV1, EV1 CV, 5th/95th percentile
  EV1, loading correlations, and EV=2/3 threshold crossings.
  `collect_pca_stability_tables(...)` exports overall and per-facet
  stability rows together as `pca_stability_all_scopes.csv`; the
  Dimensionality panel now displays and downloads the stability audit
  alongside the scree plot and loadings.
- **MML covariance and SE/CI claim audit**.
  `compute_mml_parameter_covariance(...)` now attaches an `audit` table
  with Hessian rank, rank deficiency, eigenvalue floor tolerance,
  regularized eigenvalue count, condition number, covariance diagonal
  range, claim status, and recommended action. Normal diagnostics export
  the same table as `mml_covariance_audit.csv`, and measure rows now
  carry covariance status, claim status, condition number, rank, rank
  deficiency, and regularized-eigenvalue metadata. ADEMP SE/CI coverage
  reports now add `SEBasisRisk` and `CoverageClaimStatus` so coverage is
  not cited without the matching SE basis and covariance caveat.
- **Visualization preferences and visual QA preflight**. The sidebar now
  lets users choose figure theme, label-density policy, base font size,
  label truncation length, static export width, minimum figure height,
  and caption-detail level. These settings are saved in
  `visualization_settings.json` and `visualization_settings.csv`,
  included in config JSON inspection via `visualization_preferences`,
  propagated into `figure_manifest.csv`, `visual_evidence_map.csv`,
  and caption drafts, and checked in
  `visual_qa_preflight.csv` before manuscript use.
- **BibTeX and RIS export of the cited reference library**. The Report
  tab's APA-style report now includes a "Bibliography downloads"
  expander with two buttons that emit a `.bib` (LaTeX / Zotero /
  BibDesk) or `.ris` (EndNote / Mendeley / Zotero) file containing
  exactly the references cited in the generated narrative. The
  conversion is heuristic — `_parse_apa_entry` reads each APA prose
  string in `_APA_REFERENCE_LIBRARY` and classifies it as
  `@article` / `@book` / `@incollection`, extracting authors / year /
  title / journal-or-publisher / volume / number / pages / DOI from
  the canonical APA 7 layout. RIS uses the matching `TY  - JOUR / BOOK
  / CHAP` records with one `AU  -` line per author. New helpers:
  `apa_entry_to_bibtex(key, apa)`, `apa_entry_to_ris(key, apa)`,
  `build_bibtex_from_cited(text, always_include=...)`,
  `build_ris_from_cited(text, always_include=...)`,
  `_cited_reference_keys(text, always_include=...)`. The author
  parser handles single authors, 2–4 authors with `&` joiners,
  hyphenated initials (`K.-Y.`), two-initial given names (`A. A.`),
  and `Jr.` suffixes. Titles ending in `?` / `!` are preserved (the
  Linacre 2002 RMT short note is one such case). All 66 entries in
  the live library parse without falling back to `@misc` — a
  regression test pins this so a future APA addition whose format
  breaks the parser fails the suite immediately. 24 unit tests in
  `tests/test_bibtex_ris_export.py` cover author conversion, every
  citation type, the question-mark-title edge case, the
  Winsteps.com-publisher case, proceedings entries with no
  publisher, and parity between the BibTeX bundle and the APA
  reference list. New i18n namespace `apa_report` in both
  `locales/{en,ja}.json`.
- **Configurable visualization CI level**. The classical 95 % Wald CI
  band on forest plots and EB-shrinkage error bars is now a
  user-selectable level: 50 % / 66 % / 80 % / 89 % / 90 % / 95 % /
  99 %. Sidebar control labelled "Visualization CI level" lives in
  `run_facets_mode` next to Workflow mode; selection persists in
  `st.session_state["viz_ci_level"]` and propagates to
  `_draw_measures_forest_plotly` (outer CI band + caption + title)
  and `_make_eb_shrinkage_figure` (raw and shrunk SE error bars + legend +
  title). New helpers `get_viz_ci_level()`, `_ci_z(level)`, and
  `_ci_level_pct_label(level)` centralise the level → critical value
  conversion via `Phi^{-1}((1 + c) / 2)`. The numeric `CI_Lower` /
  `CI_Upper` columns in exported measure tables stay at 95 % for
  backwards compatibility with downstream parsers — only the rendered
  visualizations change. The three pre-existing per-feature CI
  selectboxes (DIMTEST p-value, G-study bootstrap, GPCM fair-average
  SE) also broaden from `[0.90, 0.95, 0.99]` to the full
  `VIZ_CI_LEVEL_OPTIONS` list. 22 unit tests in
  `tests/test_viz_ci_level.py` pin the `Phi^{-1}` values for all seven
  canonical levels, the soft-fallback to 95 % on NaN / out-of-range
  / non-numeric input, the percent-label trailing-zero formatting,
  and the sidebar selector defaults / clamping.
- **Smart column-role detection** for the upload UI. New
  `auto_detect_column_roles(df)` helper combines a multilingual keyword
  dictionary (English + Japanese: 学生 / 受験者 / 評価者 / 採点者 / 観点 /
  項目 / 得点 / 評定 / ...) with dtype and cardinality cues to pre-fill
  the Person / Score / Rater / Criterion / Weight pickers in
  `run_facets_mode`. Score detection prefers numeric small-cardinality
  columns (≤ 12 distinct values), Person prefers high-cardinality
  non-numeric columns, and a confidence caption surfaces under each
  selectbox when the auto-detector picked the active value. Falls back
  to the legacy keyword-only `guess_col` and `pick_default_facets`
  when no Japanese / English keywords match, so existing English-only
  CSV workflows are unchanged. 19 unit tests in
  `tests/test_column_role_detection.py` cover canonical English / JA
  headers, mixed-language headers, dtype/cardinality bias,
  Weight-only-when-named, no-double-assignment, and confidence-in-unit
  interval. New i18n keys in `locales/{en,ja}.json` keep both locales
  in parity for the auto-detection captions. Companion fix:
  `_column_is_numeric` now uses a direct ratio (`coerced / non_blank
  >= threshold`) instead of `int(threshold * len(...))`, which was off
  by one on small columns and let a 2-of-4-numeric column pass the
  70 % threshold.
- Bilingual i18n scaffold (English / Japanese) for the sidebar header, app
  title, language picker, and starter messages; translation dictionaries live
  in `locales/{en,ja}.json` and resolve via the new `t()` helper with English
  fallback. Persists the choice in `st.session_state["lang"]`.
- mfrmr 0.1.6 migration coverage table and exports, extending the prior
  0.1.5 map with EB shrinkage advisory, facet sample-size adequacy,
  nesting/crossing audit, intraclass-cluster-ICC design effect,
  information curves, and misfit/weighting audit support.
- Full publication mode now computes non-person facet EB shrinkage as a
  post-hoc reporting advisory; Standard mode leaves it optional/off by
  default.
- `LICENSE_NOTICE.md` clarifies MIT licensing, commercial-use permission,
  no-warranty/as-is status, and responsible-use boundaries.
- Teacher-facing paste-data guidance now appears in the sidebar, onboarding
  text, README, and Help quick start, with a downloadable classroom template.
- Upload parsing now accepts Excel `.xlsx/.xlsm`, Apache Parquet, and JSON /
  JSON-lines in addition to CSV/TSV/TXT for rating data, person data, anchors,
  and held-out prediction rows.
- Bias/interaction output now includes a DFF-style screening table with
  sparse-cell flags, Holm adjustment, BH/FDR adjustment, practical logit
  thresholds, interpretation text, and CSV/export-bundle coverage.
- Advanced Stan-model downloads now include model-specific data templates and
  clearer scope notes for local dependence and Mixture Rasch.
- Anchor/linking output now includes an anchor/equating workflow checklist for
  current-run linking evidence and export review.
- **GPCM fair-average kernel test suite** (`tests/test_gpcm_fair_average.py`,
  13 tests). Pins reduction to the PCM/RSM softmax reference at slope = 1,
  the worked example at slope = 1.2 (1.8381511), the analytic derivative
  identity `dE[X]/d eta = a * Var[X|eta]`, invariance under the GPCM
  identification rescaling, the NaN-on-degenerate-slope contract, binary
  and singleton boundary fixtures, and the slope-1 / per-level-slope wiring
  in `calc_facets_report_tbls()`. Includes a 27-case parity fixture
  generated from the mfrmr R kernel (mfrmr 0.2.0, R 4.5.2) at 1e-10
  agreement.
- **Naming-hygiene gate** (`tests/test_naming_hygiene.py`). Parameterised
  test that scans code, tests, locales, validation docs, project docs,
  CI workflows, and project metadata for internal milestone codenames
  (numeric, single-letter, and a finite vocabulary of multi-letter forms
  such as NATO phonetics and spelled-out numbers). Public commits must
  not ship internal sprint or phase labels.
- **MML observed-information covariance** (`compute_mml_parameter_covariance`).
  For an MML fit, the structural parameter covariance is computed as the
  inverse of the observed Fisher information, evaluated as the numerical
  Jacobian of the analytical gradient returned by
  `mfrm_loglik_mml_value_grad` (O(P) gradient calls, not O(P^2) function
  calls). Inversion uses an eigendecomposition with regularization for
  near-singular spectra; a `regularized` flag and a status field
  (`"ok"`, `"regularized"`, `"not_applicable"`, `"fallback"`) propagate
  to downstream consumers so they can label the resulting SE as
  regularization-aware.
- **Structural delta-method SE / CI for the GPCM Linacre fair-average**
  (`add_gpcm_fair_average_delta_se`, `fair_average_table(fair_se=True)`).
  For each non-Person row, the routine computes
  `SE = sqrt(grad^T Cov grad)` where `Cov` is the MML observed-information
  covariance and `grad` is the finite-difference gradient of the
  fair-average scalar function with respect to the optimisation parameter
  vector. Confidence intervals are reported at the requested `ci_level`
  (default 0.95) and clipped to the rating bounds. Person rows are marked
  `"not available"` because MML EAP person estimates are not part of the
  structural Hessian; non-GPCM fits are marked `"not_applicable"`.
  Reference: standard large-sample delta method (e.g. Cramer 1946);
  observed-information identity for MML estimators (Louis 1982).
- **FACETS-style tables tab now exposes the fair-average SE / CI inline.**
  A new toggle at the top of the tab lets the user request fair-average
  standard errors and confidence intervals; a confidence-level selector
  (90 % / 95 % / 99 %) sits next to it. The toggle is disabled for non
  GPCM-MML fits with a one-line explanation. When enabled, a spinner runs
  the Hessian + delta-method evaluation once per fit, the result is cached
  in `st.session_state` keyed on a content hash of the optimisation
  parameter vector so subsequent tab interactions are instantaneous, and
  the table columns are reordered so each `Fair(M)` / `Fair(Z)` value
  sits adjacent to its S.E., CI bounds, and status flag. Method and
  citation captions are rendered below the controls when the toggle is
  on. Strings are localised in both English and Japanese.

### Changed

- **Onboarding banner is now compressed and foldable.** The 3-step
  quickstart guide moved into a Streamlit expander (default expanded for
  first-time visitors, collapsible thereafter) so the Run / Dismiss
  buttons stay visible without paying a large vertical footprint for the
  explanatory text. The redundant tip-caption row (which advertised
  past v0.2.x features and pointed at the Tutorial expander that already
  sits directly below) was dropped; the Tutorial pointer is now folded
  into the bottom of the banner body.
- **Decorative emojis removed from sample-data scenario labels and from
  the onboarding banner.** Functional status indicators (red / amber /
  green readiness colours, Wright & Linacre fit-band swatches) are
  retained because they pair text labels with a redundant colour signal
  for accessibility (WCAG 1.4.1).
- **Default sample-data scenario (`writing_essay`) description expanded.**
  Added a "why it is the default" paragraph, a "what to look at first"
  walkthrough of the result tabs in their natural reading order, and a
  "learning points specific to this data set" section that notes the
  modest rater-severity spread and the PCA sample-size caveat at n = 30.
- **`mfrmr_020_migration_coverage_table()`** documents the migration
  surface for the mfrmr 0.2.0 R package release. It inherits every
  0.1.5 / 0.1.6 row (renaming the area column), overrides the bounded-
  GPCM row to reflect this release's slope-aware fair-average and
  structural delta-method SE delivery, and appends eleven
  0.2.0-specific rows: the slope-aware fair-average and SE pair
  (Ready), the MML observed-information covariance (Ready), and nine
  Planned rows tracking the slope-aware bias inference, the FACETS
  df / ZSTD reporting layer, the network analysis triple, the G/D-
  study planner, ADEMP parameter recovery, Snijders polytomous lz*
  correction, category-specific information curves, the model-choice
  review helper, and the kable / flextable / monochrome APA presets.
  The new table is surfaced in the Help tab's coverage panel, in the
  per-run download bundle, in the disk export, and in the
  external-validation report template; the public release-readiness
  check now reports against the 0.2.0 table.
- **Category-specific information curve decomposition.** A new helper
  `build_category_information_curve_data()` returns the per-category
  contribution to the total Fisher information at every theta on the
  curve grid: `I_k(theta) = a^2 * P_k(theta) * (k - E[X | theta])^2`
  (Muraki, 1993, Eq. 16). The contributions sum to the existing
  total information curve by construction; the identity is pinned at
  machine precision across a synthetic three-level GPCM fixture. A
  new optional expander section in **Visuals -> Information
  Curves** renders the per-category breakdown as one line per
  category using a Wong (2011) colour-blind-safe palette, with the
  total information curve overlaid as a dashed line so the
  decomposition-equals-total identity is visible at a glance. A
  per-category summary table reports the peak information, the theta
  at the peak, the integrated information under the curve, and the
  category's share of the integrated total. The data is also
  downloadable as CSV. Reference: Muraki (1993) Applied Psychological
  Measurement 17(4), Eq. 16.
- **Slope-aware GPCM bias inference: information identity, likelihood-
  ratio test, and profile-likelihood confidence interval.** Three
  related fixes / additions to `estimate_bias_interaction()`:
  * **Information identity.** Under GPCM the conditional Fisher
    information for the additive bias shift `b` is
    `I(b) = sum_i a_i^2 * Var[X_i | eta_i + b] * w_i`, where `a_i` is
    the per-observation discrimination. The previous implementation
    used `I(b) = sum_i Var[X_i | eta_i + b] * w_i`, which is the
    RSM / PCM (a = 1) form; the resulting GPCM bias `S.E.` column was
    systematically too large because the `a^2` factor was missing. The
    fix is gated on `model == "GPCM"`, so RSM and PCM fits are
    byte-identical to the previous release. Reference: Muraki (1993)
    Applied Psychological Measurement 17(4), Eqs. 7 and 16.
  * **Likelihood-ratio test.** Each GPCM bias cell now reports `LR
    ChiSq`, `LR d.f.` (always 1), and `LR Prob.` from the chi-square
    pivotal `Lambda = 2 * (loglik(bias_hat) - loglik(0)) ~ chi2_1`
    (Wilks, 1938). RSM and PCM cells leave these columns as NaN with
    an explanatory `Likelihood Basis` string instead of silently
    reusing the t-based screening tier.
  * **Profile-likelihood confidence interval.** The new
    `_profile_bias_ci` helper inverts the same chi-square pivotal to
    produce a `(1 - alpha)` confidence interval for the additive bias
    shift, by solving `2 * (NLL(b) - NLL(bias_hat)) = chi2_{1, 1-alpha}`
    on both sides of `bias_hat` via Brent root-finding inside
    `[-max_abs, max_abs]`. When the likelihood never falls far enough
    inside the bracket the endpoint is returned with `Profile CI
    Status = "limited by search range"` so the caller can render the
    truncation honestly. Status values are `"ok"`, `"limited by
    search range"`, or `"not available"`. Reference: Cox (1975)
    Biometrika 62(2).
  * **UI.** The bias / interaction tab inserts `LR ChiSq`, `LR Prob.`,
    `Profile CI Lower`, `Profile CI Upper`, and `Profile CI Status`
    columns immediately after the existing t-based block when at
    least one cell carries a finite `LR ChiSq` (i.e. for GPCM fits).
    A short caption below the table cites Wilks (1938) and Cox (1975)
    and reminds readers that theta, step thresholds, slopes, and
    other facet estimates are held fixed inside the conditional
    profile. RSM / PCM fits show the unavailable caption.
  * **mfrmr 0.2.0 coverage.** `mfrmr_020_migration_coverage_table()`
    has the `GPCM bias inference - slope-aware` row and the
    `Category-specific information curves` row both updated from
    `Planned` to `Ready`. The new `PythonEvidence` cites the
    slope-aware Fisher information identity (Muraki, 1993, Eqs. 7,
    16), the LR chi-square pivotal (Wilks, 1938), the
    profile-likelihood CI inversion (Cox, 1975), and the publication-
    document integration; the `Boundary` paragraph documents the
    still-uniformly-`screening` inference tier and the conditional
    profile semantics (theta, step thresholds, slopes, and other
    facet measures held fixed). The Bounded GPCM row's `Boundary` is
    also updated to remove the obsolete "bias inference still uses
    the non-slope-aware identity pending a follow-up" caveat; the
    paragraph now points directly at the dedicated 0.2.0 bias row
    for the slope-aware machinery. The
    `external_validation_report_template` claim row for
    `mfrmr 0.2.0 migration` is updated in lock-step so the
    public-wording paragraph names slope-aware GPCM bias inference
    among the shipped pieces rather than among the roadmap items.
  * **Publication-document integration.** The Word, PDF, and HTML
    manuscript exports now embed a bias / interaction table after the
    element-measures and reliability tables. For GPCM fits the
    table includes the new ``LR ChiSq``, ``LR Prob.``, ``Profile CI
    Lower``, ``Profile CI Upper``, and ``Profile CI Status`` columns;
    RSM and PCM fits show the t-based screening columns only. Two
    new helpers ``_bias_table_for_publication`` and
    ``_bias_table_caption`` are shared across the three builders so
    the column selection and the caption stay in lock-step.
    Captions cite Muraki (1993) for the slope-aware information
    identity, Wilks (1938) for the LR chi-square pivotal, and
    Cox (1975) for the profile-likelihood confidence interval; the
    APA reference library and the citation map gain entries for
    Wilks (1938), Cox (1975), Louis (1982), Cramer (1946), and
    Muraki (1993) so a Methods narrative that mentions these
    citations in ``(Author, Year)`` form auto-resolves to a
    properly alphabetised references list.
  * **R parity.** A new test
    (`test_bias_estimation_matches_r_reference_within_tolerance`)
    runs the same GPCM MML fit and bias estimation in Python and in
    the mfrmr 0.2.0 R reference (R 4.5.2) against a shared
    deterministic CSV (`tests/data/r_bias_parity_input.csv`,
    180 rows). The R-side output (`tests/data/r_bias_parity_output.json`,
    six per-cell reference values) is generated once via `Rscript` +
    `estimate_bias()` and committed to the repository. The parity test
    verifies that Python's per-cell `Bias Size`, `S.E.`, `LR ChiSq`,
    `LR Prob.`, and both profile-likelihood CI endpoints agree with
    the R reference within tolerances calibrated to the observed
    implementation-vs-implementation convergence variance (the two
    independent MML fits converge to log-likelihoods within `~0.02`
    of each other, and the downstream cells inherit that variance);
    all observed disagreements stay well below the `0.05`-logit
    "negligible" threshold cited in Linacre's FACETS Manual. The
    categorical fields (`Profile CI Status`, `InferenceTier`) must
    agree exactly.
- **Wide-format Excel / spreadsheet upload with automatic melt to
  long format.** Excel users typically have one row per
  ``(Person, Rater)`` and one column per item being scored (e.g. C1,
  C2, C3 for criteria); the MFRM likelihood pipeline expects long
  format (one row per observation), so the upload sidebar now offers
  an inline melt step. A new ``detect_wide_format_columns(df)``
  heuristic classifies the parsed table as long or wide by counting
  numeric vs non-numeric columns: 3+ numeric columns plus 1+ id
  column flips the auto-detect to "wide" and the layout expander
  defaults to open. A new ``apply_wide_to_long_pivot(df, id_cols,
  score_cols, new_facet_name, score_col_name)`` helper performs the
  melt, dropping blank / non-numeric score cells so the likelihood
  never sees zeros for missing observations. The sidebar UI surfaces
  ID-columns / score-columns multiselects, a new-facet name input
  (default ``Criterion``), and a long-format preview after pivot. Math
  contract pinned in ``tests/test_wide_to_long_pivot.py`` (11 tests:
  detection on canonical wide / long fixtures, pivot row-count
  identity, id-column preservation, blank cell handling, refusal on
  overlapping / missing / name-collision inputs, end-to-end fit
  through ``mfrm_estimate`` on pivoted data).
- **Self-tests for the seven 0.2.0 helpers** added in this session
  (FACETS d.f. / ZSTD alignment, model-choice AICc + Akaike weights,
  G/D-study Brennan 3-way identity, design network topology, rater
  severity higher-prop convention, APA preset Markdown / HTML
  round-trip, DIMTEST nonparametric smoke). The module-level
  ``_self_test_*`` functions run at import time via
  ``run_self_tests()`` (now 55 tests, was 48) and complement the
  pytest suite by exercising the helpers on a tiny shared fixture
  ``_make_self_test_3way_rating_data()``. Each test pins a small
  closed-form contract: the Welch-Satterthwaite denominator
  identity, the Hurvich-Tsai AICc closed form and weight-sum-to-one
  identity, Brennan Eq. 3.18 G computation, node-count parity in
  the design graph, the Rater1HigherProp + Rater2HigherProp +
  TieRate = 1 identity, HTML escaping in the APA preset, and the
  Z / p-value finiteness on a small DIMTEST run.
- **R parity for the Brennan 3-way G-study decomposition** against
  lme4 REML with explicit two-way interaction terms. The R fixture
  generator now fits

      Score ~ 1 + (1 | Person) + (1 | Rater) + (1 | Criterion) +
              (1 | Person:Rater) + (1 | Person:Criterion) +
              (1 | Rater:Criterion)

  via lme4 (the three-way interaction is intentionally omitted
  because it is confounded with residual under one observation per
  cell) and dumps the variance components plus the closed-form
  G / Phi to ``tests/data/r_helper_parity_output.json``. Four new
  parity tests verify the agreement: the variance-component source
  list matches across implementations; Person is a dominant
  variance source in both; the observed (n_p, n_r, n_c) sample
  sizes round-trip; and the headline G / Phi coefficients agree to
  ~7e-2 absolute. Individual variance components can disagree at the
  zero boundary because Henderson III MoM (Python) clamps a
  negative-MoM estimate at zero and leaves the residual unchanged,
  while lme4 REML re-solves the remaining variances under the
  non-negativity constraint (Brennan 2001 p. 81; Searle, Casella &
  McCulloch 1992 Ch. 4.6) — both are valid but allocate boundary
  variance differently, which is documented in the test docstrings.
- **D-study CI bands at every planned design point**. The cluster
  bootstrap for G / Phi now propagates through each
  ``(n_facet_a, n_facet_b)`` row of the D-study forecast grid by
  reusing the same B variance-component replicates: each replicate
  re-evaluates Brennan (2001) Eq. 3.18 at every grid point so the
  marginal cost per band is one closed-form arithmetic step rather
  than another bootstrap loop. With ``bootstrap_ci=True`` the
  ``d_study`` DataFrame gains ``G_lower``, ``G_upper``,
  ``Phi_lower``, ``Phi_upper`` columns; manuscripts can now quote
  "G = 0.85, 95% CI [0.79, 0.91] at n_rater = 4, n_criterion = 6"
  directly from a single table. The observed-design row of the
  grid matches the observed-design bootstrap CI to machine
  precision (same replicate set). The Report-tab UI shows the
  bands inline; a new caption explains the cost-sharing structure.
  Math contract pinned in ``tests/test_g_d_study.py`` (CI columns
  present only when bootstrap is enabled, ordered bounds inside
  [0, 1], observed-design row matches the observed-design CI, band
  width shrinks with design size on >= 70 % of comparable pairs,
  seed determinism).
- **Nonparametric DIMTEST for essential unidimensionality** (Stout,
  1987, Psychometrika 52; polytomous adaptation per Nandakumar & Yu,
  1996). A new ``compute_dimtest_nonparametric(res, diagnostics)``
  helper tests whether the criteria measure a single latent
  dimension after the fitted MFRM main effects are accounted for.
  The test averages the rater scores per (Person, Item), partitions
  the items into an Assessment subtest (AT) and a Partitioning
  subtest (PT), stratifies persons by their PT total score, and
  computes the pooled within-stratum covariance among AT items as
  the test statistic T_L. Under essential unidimensionality T_L is
  asymptotically zero; departures signal a residual secondary
  dimension. The standard error and p-value come from a cluster
  bootstrap on persons (Efron & Tibshirani, 1993, Ch. 19) rather
  than the asymptotic null distribution; this keeps the pipeline
  purely nonparametric and side-steps the parametric bias-
  correction term T_B from Stout (1987) that would require fitting
  an auxiliary unidimensional IRT model on the AT subtest.

  AT / PT subtest selection follows Roussos & Stout (1996): the
  default ``pc2_loading_sign`` auto-split uses the PC2 loadings on
  the item facet from the existing residual-PCA pipeline. The
  auto-split inherits a known selection-bias inflation of Type I
  error (PC2 already maximises an approximately-orthogonal
  covariance direction); for a rigorous confirmatory test, pass
  ``at_items=`` and ``pt_items=`` derived from substantive
  considerations or from a separate calibration sample. The
  ``split_method`` field on the returned bundle names the active
  path.

  The Dimensionality tab gains a DIMTEST panel below the residual-
  PCA tabs with controls for the item facet, the AT / PT split
  method (auto vs user-specified multi-select), bootstrap
  replicates, and confidence level. The result reports T_L, Z,
  p-value, person count, the bootstrap CI, the AT / PT split, and a
  Reject (p < 0.05) / Retain status banner. Math contract pinned in
  ``tests/test_dimtest_nonparametric.py`` (10 tests: refusal on
  invalid input, overlap detection, unidim-fails-to-reject on a-
  priori split, 2-dim-rejects on correct split, T_L magnitude
  ordering, bootstrap determinism, bundle field completeness, the
  selection-bias caveat appears on the auto-split path, and the
  closed-form two-sided p-value identity). Three new APA references
  (Stout, 1987; Nandakumar & Yu, 1996; Roussos & Stout, 1996) are
  added to the library and citation map.
- **Bootstrap CI for G / Phi coefficients** (Efron & Tibshirani, 1993,
  Ch. 19). ``compute_generalizability_study()`` now accepts
  ``bootstrap_ci=True`` and runs a cluster-bootstrap resampling on
  persons (the object of measurement): each replicate resamples
  persons with replacement, re-runs the variance-component
  decomposition on the bootstrapped sample, and recomputes G / Phi.
  The bundle's new ``bootstrap_ci`` field reports
  ``G_lower / G_upper``, ``Phi_lower / Phi_upper``,
  ``G_replicates / Phi_replicates`` (the full replicate arrays for
  downstream diagnostics), and the active settings
  (``confidence``, ``n_bootstrap``, ``n_success``, ``method =
  "cluster_bootstrap_on_persons"``). Person-level cluster
  resampling is the standard nonparametric bootstrap for G-theory
  because it respects the within-person correlation structure that
  drives both relative and absolute error. The Report tab gains a
  Bootstrap CI expander with a checkbox, replicate count, and
  confidence-level dropdown; when enabled the panel reports the
  CI in a single caption line below the observed G / Phi metrics.
  Math contract pinned in ``tests/test_g_d_study.py`` (default-off
  behaviour, finite-ordered bands when enabled, CI widening with
  confidence level, seed determinism, [0, 1] bounds on the
  replicate arrays). One new APA reference (Efron & Tibshirani,
  1993).
- **G/D-study: full Brennan (2001) 3-way decomposition with explicit
  two-way interaction terms** for the canonical balanced p x i x j
  design with one observation per cell.
  ``compute_generalizability_study()`` now detects whether the
  design admits the full random-effects ANOVA (object_facet plus
  exactly two random facets, balanced, one observation per cell)
  and dispatches to a new
  ``_crossed_anova_two_way_variance_components`` helper that returns
  seven variance components: the three main effects, the three
  two-way interactions (``object:facetA``, ``object:facetB``,
  ``facetA:facetB``), and a Residual term that is the three-way
  interaction confounded with error (the standard one-observation-
  per-cell confounding). The method-of-moments estimators follow
  Brennan (2001, Table 3.5):

      sigma2_p   = (MS_p   - MS_pi - MS_pj + MS_pij) / (n_i n_j)
      sigma2_i   = (MS_i   - MS_pi - MS_ij + MS_pij) / (n_p n_j)
      sigma2_j   = (MS_j   - MS_pj - MS_ij + MS_pij) / (n_p n_i)
      sigma2_pi  = (MS_pi  - MS_pij) / n_j
      sigma2_pj  = (MS_pj  - MS_pij) / n_i
      sigma2_ij  = (MS_ij  - MS_pij) / n_p
      sigma2_pij = MS_pij

  Negative MoM estimates are clamped at zero per the Henderson III
  convention (Brennan 2001 p. 81). The G / Phi formulas use Brennan
  (2001) Eq. 3.18 / 3.19 with explicit interaction variances:

      sigma2(delta) = sigma2_pi / n_i + sigma2_pj / n_j +
                      sigma2_pij / (n_i n_j)
      sigma2(Delta) = sigma2(delta) + sigma2_i / n_i +
                      sigma2_j / n_j + sigma2_ij / (n_i n_j)
      G   = sigma2_p / (sigma2_p + sigma2(delta))
      Phi = sigma2_p / (sigma2_p + sigma2(Delta))

  The bundle's ``observed_coefficients.details.decomposition``
  field names the active mode (``"full_3way_brennan_eq_3_18"``
  vs ``"main_effects_only_approximation"``) so downstream
  consumers and manuscripts can name the formula they used. For
  designs with one random facet, three or more random facets, or
  unbalanced cells, the helper falls back to the existing main-
  effects-only ANOVA approximation; the existing caveat is
  preserved for that path. Math contract pinned in
  ``tests/test_g_d_study.py``: the seven-source variance-component
  layout, Brennan Eq. 3.18 closed-form G / Phi identity, and the
  decomposition label. ``mfrmr_020_migration_coverage_table()``
  boundary text updated so manuscripts no longer see "main-effects-
  only" as the default characterisation; the new boundary names
  the dispatch logic explicitly.
- **Akaike / Schwarz / AICc weights in model-choice guidance** (Burnham
  & Anderson, 2002, Eq. 2.8). Each model-choice comparison row now
  also carries an ``AICWeight``, ``AICcWeight``, and ``BICWeight``
  column. The weights are computed as ``w_i = exp(-DeltaIC_i / 2) /
  sum_j exp(-DeltaIC_j / 2)`` and admit a probability interpretation:
  ``w_i`` is the relative likelihood that model ``i`` is the
  Kullback-Leibler best of the candidate set (Burnham & Anderson,
  2002, Section 2.9). The weights sum to 1 over the finite-DeltaIC
  rows by construction; manuscripts that quote a "best model
  probability" can read it directly off the table. Math contract
  pinned in ``tests/test_model_choice_guidance.py`` (sum-to-one,
  closed-form match, best-row carries the largest weight).
- **AICc finite-sample correction in model-choice guidance** (Hurvich &
  Tsai, 1989, Eq. 1). Each model-choice comparison row now carries
  ``AICc = AIC + 2 K (K + 1) / (N - K - 1)`` alongside the existing AIC
  and BIC, plus the corresponding ``DeltaAICc`` and
  ``AICcEvidenceRatio = exp(DeltaAICc / 2)`` columns. A new
  ``AICcRecommended`` flag fires when ``N / K < 40`` (Burnham &
  Anderson, 2002, p. 66), and the tiered recommendation logic switches
  the secondary criterion from AIC to AICc whenever any candidate
  satisfies the recommendation window — small-design users no longer
  get an over-confident AIC-based ranking. The N / K ratio is reported
  per row so the applicability is recoverable from the table itself.
  Math contract pinned in
  ``tests/test_model_choice_guidance.py`` (closed-form AICc identity,
  applicability flag, ``AICc >= AIC`` bound, DeltaAICc minimum-is-zero,
  caveat-mentions-AICc). One new APA reference (Hurvich & Tsai, 1989).
- **R parity fixtures for four 0.2.0 helpers** (FACETS d.f. / ZSTD
  alignment, rater severity directed network, rater halo
  correlation network, design network connectivity). The R-side
  reference is regenerated by ``tests/data/generate_r_helper_parity_fixture.R``
  against the same deterministic CSV the existing bias-inference
  parity test consumes; the JSON output
  (``tests/data/r_helper_parity_output.json``) carries the per-facet
  d.f. / ZSTD columns from ``dx$fit`` under ``fit_df_method = "both"``,
  the directed rater pair metrics from
  ``rater_network_analysis(mode = "severity_direction")``, the
  Spearman correlations from ``rater_halo_network_analysis(method =
  "spearman")``, and the graph-level summary + per-node metrics from
  ``mfrm_network_analysis()``. Both ecosystems run JML with
  ``maxit = 200`` and ``reltol = 1e-6`` so the
  implementation-vs-implementation drift stays inside the calibrated
  0.02-log-likelihood band the existing bias-parity test
  characterises. Five new parity tests pin the agreement:
  exact counts on N / higher-counts / articulation flags, ~5e-2
  absolute tolerance on d.f., ~1.5e-1 on ZSTD, machine precision on
  the rater severity proportions, ~5e-6 on Spearman rho, and
  byte-identical topology on the design graph.

  Two **convention alignments to mfrmr** went with the R parity
  fixtures so the numerical comparison is apples-to-apples:

  * ``zstd_from_mnsq`` now guards ``df >= 1`` (matching the mfrmr
    engine convention) instead of ``df > 0``. Cells with d.f. < 1
    now report NaN ZSTD rather than a sign-flipped value driven by
    the ``1 - 2/(9 df)`` centring term. The FACETS-convention
    helper ``zstd_from_mnsq_facets`` retains ``df > 0`` with the
    explicit cap (the Wright-Masters Welch-Satterthwaite d.f. is
    fractional by design).
  * ``_rater_severity_pair_metrics`` now reports
    ``Rater1HigherProp = Rater1HigherCount / N`` (the mfrmr /
    FACETS / Eckes (2011) "share of all shared contexts" convention)
    as the primary proportion. The previous "share of disagreements"
    conditional probability remains available as
    ``Rater1HigherCondProp`` for callers that want the
    sum-to-one normalisation on the disagreement subset, and the
    new ``TieCount`` column carries the agreement count so
    ``HigherCount1 + HigherCount2 + TieCount = N`` holds by
    construction. The R parity fixture's
    ``Rater1HigherProp`` line up with the Python column at
    machine precision after the alignment.

- **Rater severity / leniency directed network and rater halo
  cross-criterion screen.** Two further network screens
  complete the mfrmr 0.2.0 network-analysis triple alongside the
  earlier design network. ``compute_rater_severity_network(res)``
  builds the directed rater-leniency graph: for every pair of raters
  with shared scoring contexts, the edge ``A -> B`` carries weight
  ``P(A > B | A != B)`` from the directional higher-count diagnostic,
  and the per-rater ``LeniencyIndex = out_mass - in_mass`` plus
  ``SeverityRank`` give a Rasch-style severity ranking that recovers
  the data-generating order exactly on the synthetic fixture (R1
  generated with the smallest ``rater_eff`` ends up at rank 1).
  ``compute_rater_halo_network(res)`` builds the
  ``(Rater, Criterion)`` node graph with Spearman edge weights and
  Bonferroni-adjusted retention; per-rater ``ReviewStatus`` of
  ``warning`` / ``review`` / ``ok`` is keyed on the configurable
  ``halo_weight_review`` and ``halo_contrast_review`` thresholds. The
  warning tier upgrades when the halo signal dominates the retained
  edges (high same-rater cross-criterion mean correlation AND either
  a clear halo / non-halo contrast or no non-halo edges retained at
  all). The Report tab gains two new panels — Rater severity and
  Rater halo — each with summary, per-rater table, expander for the
  full pair / edge metrics, and per-table CSV download. Math contract
  pinned in ``tests/test_rater_networks.py`` (16 tests covering
  refusal on invalid input, the directional-count identity
  ``Rater1HigherProp + Rater2HigherProp = 1``, severity-rank recovery
  on synthetic data, ``sum(LeniencyIndex) = 0`` under a fully crossed
  design, halo refusal on coinciding facets, warning detection on
  halo-patterned data, ``ok`` retention on clean RSM data, the
  contrast identity ``HaloNonHaloContrast = HaloMeanWeight -
  NonHaloMeanWeight``, settings round-trip, and node-count parity).
  Two new APA references (Lai, Wolfe, & Vickers, 2015; Lamprianou,
  2025) are added. The
  ``mfrmr_020_migration_coverage_table()`` row for "Network analysis
  (mfrm / rater / rater-halo)" flips from "Ready (design network);
  rater / halo screens remain follow-up work" to plain "Ready". The
  Welch test contrasting halo vs non-halo edge weights is reported
  but described as descriptive only (edge weights are clustered by
  rater and node), so ``ReviewStatus = 'warning'`` should be read as
  a follow-up trigger rather than a causal halo diagnosis.
- **Design network analysis (connectivity, articulation points, bridges,
  betweenness centrality).** A new ``compute_design_network_analysis(res)``
  helper treats the rating design as an undirected weighted graph
  (nodes = ``(Facet, Level)`` pairs, edges = co-observed in at least
  one rating row, weight = co-observation count) and computes the
  canonical connectivity diagnostics from networkx (Hagberg, Schult,
  & Swart, 2008): graph-level summary (Nodes, Edges, Components,
  LargestComponentNodes, LargestComponentShare, Density, MeanDegree,
  MeanStrength, ArticulationPoints, Bridges, Connected, Diameter,
  MeanDistance), per-node metrics (Degree, Strength, Betweenness,
  Closeness, EigenvectorCentrality) with an ``IsArticulationPoint``
  flag, per-edge bridge flag, and per-facet aggregates. The Report
  tab gains a Design Network section with the summary table, the
  per-facet aggregate, articulation-point and bridge-edge callouts
  (green-success messages when none are found), a per-node metrics
  expander, and a CSV download for the combined summary + node +
  edge + facet tables. ``min_observations`` lets the user drop weak
  edges; ``include_graph=True`` round-trips the underlying
  ``networkx.Graph`` for downstream visualisation. Math contract
  pinned in ``tests/test_design_network.py`` (14 tests covering
  refusal on invalid input, node-count identity, fully-crossed
  designs having 1 component and no cut-points, disjoint designs
  reporting Components >= 2, facet-summary aggregation parity, cut-
  node / bridge-edge subset invariants, ``min_observations`` filter,
  ``include_graph`` round-trip, and centrality [0, 1] bounds). Two
  new APA references (Csardi & Nepusz, 2006; Hagberg, Schult, &
  Swart, 2008) are added. The
  ``mfrmr_020_migration_coverage_table()`` row for "Network analysis
  (mfrm / rater / rater-halo)" flips from Planned to
  "Ready (design network); rater / halo screens remain follow-up
  work" — the directed severity / halo screens require per-rater
  pairwise comparison machinery not yet in scope.
- **APA 7 table presets and monochrome figure palette.** Two new
  helpers, ``apa_table_to_markdown(df, ...)`` and
  ``apa_table_to_html(df, ...)``, re-emit any report DataFrame as a
  manuscript-ready APA 7 snippet (bold table number, italic caption,
  header / body grid with top + header + bottom rules, ``Note.`` block
  below). Float digits are user-configurable; boolean columns render
  as ``Yes`` / ``No``; missing values render as the empty string; HTML
  output is HTML-escaped end-to-end (cells, captions, notes) so
  ``<script>`` payloads cannot break out of the cell. A companion
  monochrome palette is exposed via ``get_monochrome_palette()`` (eight
  grayscale steps from near-black to very light grey) and
  ``get_monochrome_dash_patterns()`` (six Plotly dash patterns) so
  figure exports under a print-friendly preset stay distinguishable
  even at low DPI. The Report tab gains an APA-presets panel with a
  table picker, format radio (Markdown / HTML), digits / table-number
  / caption / note controls, an inline preview, and per-format
  download buttons; the monochrome palette is previewed with colour
  swatches and dash-pattern names. Math contract pinned in
  ``tests/test_apa_presets.py`` (18 tests). The
  ``mfrmr_020_migration_coverage_table()`` row for "APA table presets
  (kable / flextable / monochrome)" flips from Planned to Ready. The
  helpers complement the existing single-file Publication Document
  rather than replace it: use the APA preset for per-table snippets
  authors paste into a longer manuscript, and the Publication
  Document for the bundled Word / PDF / HTML deliverable.
- **Generalizability theory (G-study) and D-study forecast.** A new
  ``compute_generalizability_study(res)`` helper performs a method-of-
  moments crossed-ANOVA decomposition on the rating data and reports
  the canonical generalizability-theory outputs (Cronbach, Gleser,
  Nanda, & Rajaratnam, 1972; Brennan, 2001): per-source variance
  components with proportion-of-variance shares, G (relative decision)
  and Phi (absolute decision) coefficients at the observed design,
  and a D-study forecast grid that scales G / Phi over planned
  numbers of raters and criteria. The dispatch follows the standard
  closed form: ``G = sigma2_p / (sigma2_p + sigma2_e / n_total)`` and
  ``Phi = sigma2_p / (sigma2_p + sum(sigma2_facet / n_facet) +
  sigma2_e / n_total)``. The Report tab gains a G-study / D-study
  section with the variance-component table, an observed-design
  G / Phi metric strip, the D-study forecast table, and a CSV
  download; the panel runs every page load because the math is fast
  (a few ANOVA reductions on the rating data). Math contract pinned
  in ``tests/test_g_d_study.py`` (15 tests: refusal on invalid input,
  variance-component non-negativity, required source list, proportion-
  of-variance sum = 1, closed-form G / Phi identity, ``Phi <= G``,
  D-study monotonicity in each facet sample size, D-study grid
  contains the observed design, row count = product of facet grids,
  bounds in [0, 1], and Brennan / Cronbach citation hygiene). Two
  new APA references (Brennan, 2001; Cronbach, Gleser, Nanda, &
  Rajaratnam, 1972) are added. ``mfrmr_020_migration_coverage_table()``
  flips the "G/D-study planning (mfrm_d_study)" row from Planned to
  Ready. The implementation deliberately uses the one-observation-
  per-cell main-effects-only approximation (no explicit person-by-
  facet interaction terms), which is the standard simplification when
  ``lme4`` random-effects fitting is not available; the caveat
  surfaces in the bundle's ``caveat`` field for manuscript citation.
- **Parameter recovery (ADEMP) simulation.** A new
  ``evaluate_parameter_recovery(...)`` helper runs a Monte Carlo
  parameter-recovery study under one explicit data-generating setup
  and reports the performance measures of Morris, White, and Crowther
  (2019): bias, MCSE(bias), RMSE, MCSE(RMSE), MAE, raw / aligned
  errors, Pearson correlation, mean SE, SE-availability rate, and
  95 % coverage. The data-generating mechanism samples person
  measures, rater severities, and criterion difficulties from
  mean-zero independent normals; step thresholds are equally spaced
  over a chosen logit span (RSM uses shared thresholds; PCM and GPCM
  use per-criterion perturbed thresholds). Each replicate refits the
  requested model via ``mfrm_estimate``, runs a lightweight
  ``mfrm_diagnostics`` pass to populate SE on the measures table,
  mean-aligns the location-indeterminate blocks (Person, Rater,
  Criterion) against the truth, and computes per-row error +
  coverage95. The bundle carries six fields: ``recovery``
  (long table), ``recovery_summary`` (per-parameter-type aggregates),
  ``coverage_report`` (SE/CI coverage diagnostic),
  ``rep_overview`` (per-rep convergence and timing), ``ademp``
  (Aims / Data-generating mechanism / Estimands / Methods /
  Performance measures narrative), and ``settings`` (round-trip of
  the input arguments). The Report tab gains a Parameter-recovery
  (ADEMP) section with eight controls (model, fit_method,
  replicates, seed, n_person, n_rater, n_criterion, n_cat) and
  caches the bundle in session_state keyed on the input settings.
  Math contract pinned in ``tests/test_parameter_recovery.py``
  (18 tests covering refusal on bad inputs, three location-block
  ParameterType values, long-table row counts, the mean-alignment
  identity at machine precision, Bias = mean(ErrorAligned), RMSE /
  MAE closed forms, Coverage95 closed form, SE/CI status propagation,
  explicit coverage-report construction, fixed-seed determinism,
  ADEMP narrative completeness, and correlation positivity on a
  clean RSM fit). One new APA reference (Morris, White, & Crowther,
  2019) is added to the library and citation map. The
  ``mfrmr_020_migration_coverage_table()`` row for "Parameter
  recovery (ADEMP)" flips from Planned to Ready.
- **Model-choice guidance: RSM / PCM / GPCM side-by-side comparison.**
  A new ``compute_model_choice_comparison(res)`` helper refits the two
  non-current rating-scale models on the same data and returns a
  publication-style comparison bundle. Each row carries ``Model``,
  ``FitStatus``, ``N``, ``KParams``, ``LogLik``, ``AIC`` (Akaike,
  1974), ``DeltaAIC``, ``AICEvidenceRatio``, ``BIC`` (Schwarz, 1978),
  ``DeltaBIC``, and ``BICEvidenceRatio``; the evidence ratios use
  the Akaike / Schwarz closed form ``exp(DeltaIC / 2)`` (Burnham &
  Anderson, 2002, Eq. 2.10). A nested likelihood-ratio table reports
  ``Lambda = 2 * (LL_alt - LL_null) ~ chi2_{df_alt - df_null}`` for
  every pair where both fits succeeded (RSM in PCM, PCM in GPCM, RSM
  in GPCM; Wilks, 1938). A tiered recommendation (``strong`` /
  ``moderate`` / ``weak`` / ``tie``) is keyed on the BIC gap to the
  second-best candidate; the rationale string names the per-criterion
  margin so manuscripts can quote it directly. The Report tab gains
  a "Model-choice guidance" section that exposes the comparison behind
  a "Run model-choice comparison" button; the result is cached in
  session_state keyed on a content-addressed hash of the result so
  re-loading the same fit skips the refit. Runs with anchors,
  population / latent-regression terms, or facet regularization are
  refused with an explanatory reason because their likelihoods are
  not directly comparable on a common scale. Math contract pinned in
  ``tests/test_model_choice_guidance.py`` (13 tests: invalid-input
  refusal, regularization / population-formula refusal, three-row
  output, DeltaIC = 0 at the minimum, evidence-ratio closed form,
  LR chi-square = 2 (LL_alt - LL_null), scipy parity on p-values,
  preference recovery on RSM-generated and GPCM-heterogeneous-slope
  data, caveat citation, column-order stability, fit_times zero for
  current model). Three new APA references (Akaike, 1974; Schwarz,
  1978; Burnham & Anderson, 2002) are added to the library and
  citation map. The ``mfrmr_020_migration_coverage_table()`` row
  for "Model-choice guidance (RSM / PCM / GPCM)" flips from Planned
  to Ready.
- **FACETS d.f. / ZSTD reporting alignment.** Two d.f. conventions for
  the standardised Infit / Outfit fit statistic now ship side-by-side.
  The engine convention (default, backwards compatible) uses
  ``DF_Infit = sum(Var * w)`` and ``DF_Outfit = sum(w)`` and reports
  the Wilson-Hilferty (1931) standardised value through
  ``zstd_from_mnsq``. The new FACETS / Wright-Masters convention
  (``Wright & Masters, 1982``, Eqs. 4.20 / 4.27) applies a Welch-
  Satterthwaite d.f. that captures variance heterogeneity through the
  per-observation fourth central moment
  ``M4_i = sum_k P_k * (k - E[X | theta])^4``:
  ``DF_Infit_FACETS = 2 * (sum Var * w)^2 / sum w * (M4 - Var^2)`` and
  ``DF_Outfit_FACETS = 2 * (sum w)^2 / sum w * (M4 / Var^2 - 1)``.
  ``compute_obs_table`` now adds a ``FourthCentralMoment`` column;
  ``calc_overall_fit`` and ``calc_facet_fit`` accept a ``fit_df_method``
  argument with values ``"engine"``, ``"facets"``, or ``"both"``; and
  ``zstd_from_mnsq_facets`` applies the FACETS ``+/- 9`` cap (with the
  cap exposed as a parameter so users can disable it). The Fit Details
  tab has a three-way radio that switches the on-screen fit table
  between conventions without re-fitting, plus a metadata caption that
  names the active method, ZSTD transform, and cap on every render.
  An Import-FACETS-fit-table expander accepts a CSV/TSV with ``Facet /
  Level / InfitZSTD / OutfitZSTD`` columns, joins on ``(Facet, Level)``,
  and renders a per-row delta table (Python minus imported) so users
  can verify the Python output against an external FACETS run. Math
  contract pinned in ``tests/test_facets_df_zstd_alignment.py`` (15
  tests covering the closed-form fourth moment, ``M4 >= Var^2`` bound,
  Wilson-Hilferty cap behaviour, Welch-Satterthwaite formula, and the
  three-way dispatch column shape including the engine-only backward-
  compatibility guarantee). The ``mfrmr_020_migration_coverage_table()``
  row for "FACETS df / ZSTD reporting alignment" flips from Planned to
  Ready. References: Wright & Masters (1982); Wilson-Hilferty cube-root
  transformation as cited in Linacre (2002).
- **Polytomous person-fit indices: lz (Drasgow et al., 1985) and lz*
  (Snijders, 2001).** A new ``compute_person_fit_indices(res)`` reports
  the standardised polytomous log-likelihood for every fitted person.
  Under JMLE the report carries the Snijders (2001, Eq. 16) correction:
  ``c_n = Cov[l, S] / I(theta)``, corrected variance
  ``Var[l] - Cov[l, S]^2 / I(theta)``, and ``lz* = (l - E[l] - c_n * S)
  / sqrt(corrected_var)``. Under MML / EAP the function reports the
  unadjusted ``lz`` with the status field
  ``"not_applicable_eap"`` so the unadjusted form is never silently
  treated as Snijders-corrected. The six per-observation diagnostic
  columns ``PrObserved``, ``ItemEntropy``, ``ItemVarLogP``,
  ``ItemLogPScoreCov``, ``ScoreInformation``, and
  ``ObservedScoreDerivative`` are added to ``compute_obs_table(res)``
  in a single pass over the per-observation probability matrix; the lz
  / lz* aggregation reads them by name. Reporting columns
  (``ReportIndex``, ``ReportValue``, ``ReportFlagLevel``,
  ``ReviewStatus``, ``ReviewReason``, ``ReportCaveat``) choose lz*
  when computed and fall back to lz otherwise, using the standard
  normal thresholds ``|z| > qnorm(0.975)`` (5 %) and
  ``|z| > qnorm(0.995)`` (1 %). The Measures tab surfaces a compact
  panel: a four-metric strip (persons measured, flagged at 5 %, flagged
  at 1 %, chosen report index), a flagged-rows table sorted by
  ``|ReportValue|`` when any person trips the practical threshold, the
  full per-person table behind an expander, a CSV download, and a
  method-aware caveat (``Snijders correction conditional on fitted
  non-person calibration`` vs ``uncorrected lz; interpret tail z-scores
  conservatively``). The lz / lz* indices also ship in the per-run
  downloads bundle as ``person_fit_indices``. Two new APA references
  (Drasgow, Levine & Williams, 1985; Snijders, 2001) are added to the
  library and citation map so a Methods narrative that mentions either
  in ``(Author, Year)`` form auto-resolves to a properly alphabetised
  references list. The math contract is pinned in
  ``tests/test_person_fit_indices.py`` (17 tests covering the closed-
  form lz identity, the closed-form lz* / c_n / corrected-variance
  identity, the threshold rule, status-field fallback for MML and for
  missing-column / degenerate inputs, the small JMLE end-to-end fit,
  and column-order stability). The ``mfrmr_020_migration_coverage_table()``
  row for "Polytomous person-fit correction (Snijders lz*)" flips from
  Planned to Ready with the new evidence narrative.
- **Remaining six sample-data scenarios** (large-scale writing, L2
  speaking, clinical OSCE, writing with missing, music peer-rating,
  reading testlet binary) now follow the same expanded template:
  a one-paragraph design summary, a "what this scenario teaches"
  paragraph, a tab-by-tab "what to look at first" walkthrough, and
  scenario-specific learning points (PCA at n = 120 for the severity-
  outlier scenario, PCM/GPCM step thresholds for the L2 speaking
  criterion contrast, low person separation as a feature for the
  clinical OSCE cohort, the missing-rate and connectivity readiness
  flags for the MAR scenario, the sparse Person × Rater graph for the
  round-robin peer rating, and the local-item-dependence / Stan
  TESTLET_RI hand-off for the binary reading testlet).
- **Decorative emojis removed across the app and i18n.** Help text,
  status messages, button labels, dict-keyed tab labels, and section
  comments no longer carry decorative emoji ornaments (document,
  chart, inbox, paste, microscope, abacus, test-tube, clock, shuffle,
  book, notepad, ruler, floppy, party-popper, alarm-bell, and other
  category badges). Status / severity markers are explicitly retained
  because they pair text labels with a redundant colour signal for
  WCAG 1.4.1 accessibility: red / amber / green readiness lights, the
  Wright & Linacre fit-band swatches, the success / warning / error
  severity icons, the no-severity default circle, and the structural
  yes / no markers used inside capability-matrix tables. Toast
  notifications now use consistent severity icons (`icon="✅"` for
  success, `icon="❌"` for failure) instead of celebratory or alarm-
  bell glyphs.
- Design evaluation and report/download bundles now include the new
  0.1.6-oriented audit tables so Essential-view users see warnings in the UI and
  expert users can export full CSV/Excel evidence.
- Config JSON import contract is pinned at 25 replayable settings so restored
  runs preserve the current analysis-depth and diagnostic controls.
- Public-facing README, Help, release-check, and validation wording now avoid
  machine-specific paths and internal workspace terminology.
- Reproducible script downloads now include a local batch workflow note so
  repeated refits stay outside hosted Streamlit.

### Fixed

- **GPCM Linacre fair-average on non slope-facet rows** (Person, and any
  non-step facet such as Rater or Task under a `step_facet = Criterion`
  fit) is now evaluated at slope = 1, the geometric-mean discrimination
  under the sum-to-zero log-slope identification. The previous code
  substituted the arithmetic mean of the estimated slopes, which by the
  arithmetic-mean / geometric-mean inequality is always >= 1 and so
  reported a systematically over-discriminated fair average there. The
  step-facet's own rows continue to use that level's own discrimination
  and that level's own step thresholds. Reported FairM / FairZ values
  in the FACETS-style report table will change for GPCM fits; RSM and
  PCM fits, and non slope-facet rows under PCM, are unaffected.
  Reference: Muraki (1992), Eqs. 2-3 and 10; Muraki (1993), Eqs. 7
  and 16; Linacre, FACETS Manual section "Fair Average". Cross-checked
  against the mfrmr R reference (0.2.0) at 1e-10 across 27 fixture
  cases.
- `expected_score_from_eta()` now returns `np.nan` for a non-finite or
  non-positive slope instead of silently falling back to slope = 1, so
  a degenerate GPCM fit can no longer masquerade as a healthy one
  through the fair-average path.
- Data source radio now shows the 7 built-in samples plus Paste and Upload
  exactly once; Paste/Upload are no longer repeated after every sample option.

## 0.2.14-beta - 2026-04-18

User feedback after v0.2.13: documentation said "four scenarios" in
several places after a fifth was added; the readiness panel blocked
the binary testlet scenario with a hard [ISSUE] for its single-scorer
facet (working as designed but UX-hostile); the Run history area
needed a FACETS-style one-click bundle download; and two more
scenarios (missing data + peer-rating) were requested. This release
addresses all five points and — once again — the v0.2.12 end-to-end
gate caught two real bugs in CI before users would have hit them.

### Added

- **Two new sample scenarios** (registry now has 7 total):
  - **Writing with missing** (80×4×2×3, ~1,632 obs after 15 %
    MAR drop) — exercises MFRM's tolerance of missing-at-random
    observations, the readiness-panel missing-rate flag, and the
    connectivity check. (Little & Rubin, 2002)
  - **Music peer-rating** (120 musicians × 2 cyclic raters × 2
    pieces × 4 criteria, 1,920 obs) — sparse round-robin design
    where every musician is both an examinee and a rater. Person_i
    rated by Person_(i+1) and Person_(i+2). 240 unique pairs out
    of 14,400 possible. (Myford & Wolfe, 2004; Linacre, 2007)
- **`render_quick_results_download()`** — FACETS-style one-click
  results bundle surfaced right below the Run history panel.
  Offers ZIP (CSVs) and Excel (multi-sheet) variants. Users
  can download the full results bundle without discovering the
  Downloads sub-tab. Built via new helper
  `build_result_bundle_frames()`.
- **`_generate_mfrm_peer_rating_data()`** — sparse cyclic
  generator for round-robin peer-assessment data. The first
  generator in the registry that does NOT assume a fully-crossed
  design.
- **`missing_rate` parameter support** in the scenario registry —
  any scenario can declare a stochastic MAR drop applied after
  generation (separate RNG stream so the response seed is
  unaffected).
- **`tests/test_v0214_scenarios.py`** (12 tests) covering the
  missing-rate drop, the peer-rating generator's structural
  invariants (no self-rating, exactly 2 raters per examinee,
  shared ID pool), the result-bundle helper, and the singleton-
  facet readiness downgrade.

### Changed

- **Readiness check for singleton facets downgraded** from
  [ISSUE] to [CAUTION]. MFRM centers a singleton facet to 0
  without error; a hard block was overkill for legitimate single-
  scorer / single-form designs (e.g., binary reading testlets).
  The detail message now explains the model behaviour.
- **`pick_default_facets` post-filter**: columns with < 2 unique
  levels are excluded from the auto-selected facet set so users
  do not have to deselect them manually for every scenario.
- **Onboarding banner is now dynamic** — scenario list builds from
  `SAMPLE_DATA_SCENARIOS`, no more hardcoded "four built-in
  scenarios" text. Adding a scenario will not desynchronise docs.
- **Sidebar Data source caption** uses the same dynamic count.
- **Inline `_self_test_sample_data_scenarios`** allows ±5 %
  tolerance on row count for scenarios with `missing_rate > 0`.
- **README.md** scenario list expanded to mention all five built-
  ins (was four).

### Fixed (caught by the v0.2.12 end-to-end gate before release)

- **`calc_interrater_agreement` `observed=False` groupby
  expansion** crashed `pd.insert` with
  `ValueError: Length of values (960) does not match length of
  index (57600)` on the music peer-rating scenario. The
  Categorical-grouped agg was materialising the full 120 ×
  N_context Cartesian even though only 240 cells had data. Fix:
  `observed=True` on both group-by sites so empty cells are
  skipped — correct behaviour for sparse designs and identical for
  fully-crossed ones.
- **`guess_col` greedy substring match** (carryover from v0.2.13):
  reading-testlet column auto-detect picked `Scorer` for the
  `score` pattern, which collided after rename. Fix landed in
  v0.2.13.

### Verification

- 41 inline self-tests pass.
- **226 pytest tests pass** (was 200) — +26: 12 v0.2.14 tests, 7
  binary/testlet (carried), and 7 across the existing files for
  the new scenarios.
- End-to-end AppTest covers all 7 scenarios × Run without any
  exception (≈ 3 min CI runtime).

## 0.2.13-beta - 2026-04-18

Binary / testlet data support, plus two bugs the new end-to-end
AppTest gate (introduced v0.2.12) caught before release. Demonstrates
the gate working as intended: issues surface in CI, not production.

### Added

- **`reading_testlet_binary` sample scenario** — 100 examinees × 1
  scorer × 6 passages × 4 items = 2,400 binary (0/1) observations.
  Grounded in the testlet literature (Wainer & Kiely, 1987;
  Bradlow, Wainer & Wang, 1999). Columns are
  `Person / Scorer / Text / Item / Score`, the tidy-long shape that
  reading / listening / SJT tests produce.
- **Help → Model Capability → "Testlet-format Data"** section
  (~120 lines). Compares MFRM to classical testlet models
  (Bradlow+ 1999) and the bi-factor equivalent (Rijmen, 2010;
  DeMars, 2006). Includes a recommendation flow for when to use
  MFRM's fixed-effect Text facet vs. the **TESTLET_RI** /
  **TESTLET_BIFACTOR** Stan downloads for random-effects
  modelling of local item dependence.
- **Five new APA 7 references**:
  Bradlow, Wainer & Wang (1999), DeMars (2006), Rijmen (2010),
  Wainer & Kiely (1987), Wang, Bradlow & Wainer (2002). All
  wired into `_APA_REFERENCE_LIBRARY` + `_CITATION_TO_KEY`.
- **`tests/test_binary_testlet.py`** (7 tests) — schema, JMLE
  convergence on 2,400-row binary data, citation resolution,
  Help-doc content pin, generator determinism.

### Fixed (bugs the new end-to-end gate caught before release)

- **`guess_col` greedy substring match.** Column auto-detection on
  a testlet dataset (`Person, Scorer, Text, Item, Score`) picked
  `Scorer` as the Score column because `"score" in "scorer"` is
  True. After rename-collision, `pd.to_numeric(df["Score"])` hit
  a DataFrame with two "Score" columns and raised `arg must be a
  list, tuple, 1-d array, or Series`. Fix: `guess_col` now makes
  two passes — exact case-insensitive match first, substring match
  only if nothing exact was found.
- **`prepare_mfrm_data` silent rename-collision** — the same class
  of bug could still surface if a user-uploaded dataset had a
  facet literally named "Score" or "Person". Added a
  `df.columns.duplicated()` check with a clear error message.
- **`np.nanvar(..., ddof=1)` runtime warning** on singleton-facet
  groups (the new binary testlet has a single scorer). Pre-check
  `len(finite_est) >= 2` to avoid the numpy warning on
  zero-variance slices.
- **APA citation regex** rejected 3-author forms like
  `(Bradlow, Wainer & Wang, 1999)`. Loosened to accept any
  `(Surname…, YYYY)` shape.

### Unchanged on purpose

- Binary support was *already* present in the core
  RSM/PCM/GPCM estimator — a Score column with 0/1 values and
  `n_cat = 2` produces a single Rasch–Andrich threshold, which
  is mathematically equivalent to the 1PL Rasch model. This
  release just documents and tests that pathway and adds the
  testlet-shaped scenario that exercises it.
- For discrimination-aware binary analysis (2PL), the existing
  sidebar **Advanced models → IRT_2PL_BINARY** Stan generator
  remains the recommended download-only path.

### Verification

- 41 inline self-tests pass.
- **200 pytest tests** pass (was 186) — +14 including 7 new
  binary/testlet tests and the updated regression guards.
- End-to-end AppTest covers all 5 scenarios × Run without
  exception (≈ 2 min CI runtime).

## 0.2.12-beta - 2026-04-18

Process hotfix responding to a user question after the v0.2.11
widget-crash ship-and-patch cycle: *"Why did this happen and aren't
there other bugs like it?"*

Honest root-cause: v0.2.7 added four sample scenarios including
Clinical OSCE (3 competencies) and the Writing Essay default (4
raters). Tests verified the generator functions but never ran the
generated data through the full UI pipeline. This release installs
the CI gate that would have caught the v0.2.10 bug and completes the
systematic audit of related hazards.

### Added

- **`tests/test_end_to_end_scenarios.py`** (5 tests). For every key
  in `SAMPLE_DATA_SCENARIOS`, boot the app via
  `streamlit.testing.v1.AppTest`, switch the data-source radio to
  that scenario, click **Run FACETS-mode estimation**, and assert
  `at.exception` stays empty across the full post-run tab render.
  This would have caught the v0.2.10 `StreamlitValueBelowMinError`
  before release. Total end-to-end runtime ≈ 2 min on a CI runner.
- **Three additional boundary-case tests** in
  `tests/test_small_data_widgets.py`: 1-row forest plot, 1-row
  misfit ranking, 2×2 bias heatmap. These cover user-provided data
  shapes that do not match any built-in scenario but are plausible
  in real research (e.g., a pilot design with 2 raters).

### Audited (no additional fixes required)

A 5-pattern systematic sweep of the whole monolith found no further
reachable crashes:

1. **Data-dependent widget bounds** — 18 call sites scanned. The 3
   unsafe ones were already fixed in v0.2.11; the remaining 15 use
   literal bounds only.
2. **Statistical calls on small samples** (`std(ddof=1)`,
   `quantile`, `percentile`) — 7 call sites. All either guarded by
   `if n > 1` / `if levels > 1` or downstream of a `len(elements) <
   2` return (Facet Equivalence, Rater Dashboard, Q-Q plot).
3. **Division by data length** — 8 call sites. All guarded by
   `if n_ok == 0 / len(data) > 0 / max(total, 1)` or similar.
4. **Facet Equivalence chi-square** — guarded by
   `len(elements) < 2: return` at L16948 so `df_chi = n_elem - 1 ≥
   1`.
5. **Q-Q residual plot `.std(ddof=1)`** — guarded by `residuals.size
   < 20: return` at L18763.

### Process improvement

- Scenario additions now require an end-to-end AppTest gate before
  release. Unit tests for the generator alone are not sufficient:
  the bug class that shipped was "valid DataFrame → widget
  constraint violation", which only the full pipeline exercises.
- Total pytest count: 178 → **186** (+5 end-to-end, +3 boundary).
- 41 inline self-tests continue to pass.

## 0.2.11-beta - 2026-04-18

Widget small-data hotfix. A user reported a hard crash with
`StreamlitValueBelowMinError` on the Visuals tab when running against
a facet with fewer than 5 elements (e.g. the Clinical OSCE scenario
has 3 competencies). The forest-plot widget's `min_value=5` was
incompatible with `value=min(40, len(sub))` when `len(sub) < 5`. An
audit found two sibling widgets with the same data-dependent
min/value hazard.

### Fixed

- **`_draw_measures_forest_plotly`** — when the selected facet has
  fewer than 6 rows, skip the "Show top N" number input entirely
  and plot every row. When it has ≥ 6, clamp the default to the
  widget's own range so `value >= min_value` always holds.
- **`_draw_misfit_ranking`** — the `st.slider("Number of elements
  to show")` had `min_value=5, max_value=len(df)`; on fit tables
  with `len(df) < 5` this flipped to max<min and triggered
  `StreamlitAPIException`. Now auto-skips the widget for tiny
  tables (< 6 rows).
- **Bias heatmap axis-label slider** — `max_axis_default` could
  exceed `max_value` when the pivot was narrow. Now all three of
  (min, max, default) are clamped to the pivot's actual extent,
  and the slider is skipped entirely for ≤ 8-label pivots.

### Added

- **`tests/test_small_data_widgets.py`** (3 tests): directly invokes
  `_draw_measures_forest_plotly` and `_draw_misfit_ranking` on
  3-row fixtures, and runs the full AppTest with the Clinical OSCE
  scenario. Pass condition = no `StreamlitValueBelowMinError` /
  `ValueAboveMaxError` / `APIException` raised. These are the
  cheapest guards against the class of bug that just shipped.

### Unchanged on purpose

- All v0.2.10 one-click PDF / Word / HTML downloads and matplotlib
  figure fallbacks continue to work unchanged.
- Non-affected widgets (maxit, reltol, rating_min, rating_max,
  number-of-latent-classes, minimum-common-anchors, etc.) all have
  literal min/max/value bounds and are not at risk.

## 0.2.10-beta - 2026-04-18

Publication-Document hotfix. On Streamlit Community Cloud, the
Publication PDF shipped without any plots because kaleido 1.x
requires a Chrome binary that Cloud does not install — every figure
call silently returned `None` and the PDF ended up as a
tables-and-text document. This release adds a matplotlib fallback
that never needs Chrome, and streamlines the download UI so users
can grab the PDF in one click.

### Added

- **Matplotlib fallback for publication figures.** Four new helpers
  render the core figures directly from the result / diagnostics
  dicts, bypassing the Plotly → kaleido → Chrome chain entirely:
  - `_mpl_wright_map_png()`
  - `_mpl_fit_scatter_png()`
  - `_mpl_category_probability_png()`
  - `_mpl_facet_distribution_png()`
- **Hybrid selector** `_plotly_or_matplotlib_png()` tries the Plotly
  export first (for parity with the on-screen look when Chrome is
  available) and falls back to matplotlib when it returns `None`. So
  the Publication Document always ships with its core figures,
  regardless of host environment.
- **`matplotlib>=3.8`** added to `requirements.txt`.
- **New pytest file** `tests/test_publication_figures.py` (10 tests)
  pins: every matplotlib helper produces valid PNG bytes, the
  payload integration function emits all 4 figures, the full PDF
  builder's output contains compressed image streams, and empty
  inputs never raise.
- **Strengthened `_self_test_publication_document_pdf`** — the test
  now enforces `len(pdf) > 100 KB` and `b"FlateDecode" in pdf` so
  the regression where figures silently vanished will fire the next
  time it happens.

### Changed

- **One-click Publication Document downloads.** The previous
  two-step "Generate → wait → Download" flow is replaced by a
  single prominent download button per format. The bytes are built
  eagerly on tab entry with a visible spinner, and cached in
  `session_state` keyed by a run-level fingerprint so subsequent
  sidebar reruns reuse the cache rather than rebuilding.
- **PDF button promoted to primary action** (`type="primary"`) and
  relabelled "Download PDF (with plots)" to make its content
  explicit after the figure-regression fix.
- Consistent button labels across formats: "Download Word (.docx)",
  "Download PDF (with plots)", "Download HTML".

### Unchanged on purpose

- The Plotly export path is still tried first when kaleido + Chrome
  are available; it produces PNGs that match the on-screen look
  exactly. The matplotlib fallback only fires when that path
  returns `None`.
- Figure content, captions, and ordering remain identical to v0.2.9.
- Word / HTML builders unchanged apart from sharing the same
  one-click UI pattern; their existing embedding paths already
  benefited from the shared `_publication_figure_payloads()` fix.

## 0.2.9-beta - 2026-04-17

Quality hotfix responding to a UX re-audit. Three gaps identified in
the v0.2.8 scoring review are addressed in this release:
verbose captions that over-explained features, missing ingestion-time
outlier detection, and thin pytest coverage.

### Added

- **Ingestion-time outlier detector** (`detect_data_outliers`).
  Screens uploaded / pasted / sample data for six anomaly types
  before the user clicks Run, and appends findings to the
  readiness panel with appropriate severity:
  1. Out-of-range scores (when rating_min / rating_max are set)
  2. Negative score values (rating scales are conventionally ≥ 0)
  3. Zero-variance persons (all-identical responses)
  4. Zero-variance facet elements (rater / task / item with
     constant scores — catches classic central-tendency raters)
  5. Extreme-frequency persons via Tukey fence (k=3)
  6. Ceiling / floor saturation (> 90 % of scores at scale extremes)
- **Six new pytest files** covering previously untested surfaces.
  Total pytest count: 10 → 165.
  - `tests/test_sample_scenarios.py` (29 tests) — every scenario's
    shape, determinism, citation resolution, PCA-readiness
  - `tests/test_data_outliers.py` (12 tests) — every anomaly type,
    clean-data-zero-findings invariant, defensive empty-data handling
  - `tests/test_readiness_report.py` (9 tests) — severity logic,
    outlier integration, overall-max invariant
  - `tests/test_apa_references.py` (12 tests) — library integrity,
    citation-map round-trip, key-naming convention, sorted output
  - `tests/test_help_popover_library.py` (58 tests) — 18 topics ×
    required-field check, v0.2.6 additions, unknown-key no-op
  - `tests/test_response_model_generators.py` (29 tests) — seven
    Stan generators + dispatcher, brace balance, required blocks
  - `tests/test_config_whitelist.py` (6 tests) — 23-key contract,
    filter drop-through, critical-keys always present
- **`_self_test_data_outlier_detection`** inline self-test (#41)
  pins the six detection types against a synthetic dirty fixture.

### Changed

- **Verbose sidebar captions trimmed.** Twelve captions and `help=`
  strings shortened to 1-sentence form. Full explanations moved to
  Help → Interpretation Guide rather than repeated at every
  widget. Affected widgets: Column mapping caption, Facet selector
  help, Workflow mode help, Model radio help, Advanced-models
  caption, population_formula help, standardize-numeric help,
  categorical-terms help, noncenter_facet help, dummy_facets help,
  positive_facets help, Analysis depth help, residual-PCA help,
  strict-marginal help, anchor-templates caption.

### Unchanged on purpose

- The readiness panel's existing status UI absorbs outlier
  findings without any rendering changes — it was already designed
  for extensible check lists.
- Existing self-tests, sample scenarios, popover library, and
  Essential-mode tab filters all unchanged.

## 0.2.8-beta - 2026-04-17

Sample-scenario UX hotfix. The v0.2.7 scenario registry hid the
switcher behind a two-step flow: users had to first select
"Sample data (built-in)" on a radio, then spot a selectbox below
labelled "Sample scenario". Early testers missed that the selectbox
existed and assumed the app shipped a single sample. This release
pulls the scenario choice out of hiding.

### Changed

- **Flat single-radio data-source picker.** The sidebar now shows
  one radio with six options: the four built-in scenarios (with
  label + observation count) followed by Paste CSV/TSV and Upload
  file. No more two-step "Sample data → Sample scenario" flow.
- **Inline info card** replaces the collapsed expander as the
  primary surface for scenario metadata. The selected scenario's
  name, design dimensions, observation count, and one-line
  diagnostic hint are always visible in the sidebar directly under
  the radio.
- **"Try another scenario" quick-switch buttons.** When any scenario
  is loaded, the sidebar now shows one-click buttons for each of the
  other three scenarios so users can hop between designs without
  scrolling back to the radio. Each button is tooltipped with the
  scenario's short description.

### Added

- **Main-area data banner** (`render_loaded_data_banner`). When a
  built-in scenario is loaded, a one-line callout appears above the
  readiness panel naming the scenario, its dimensions, and pointing
  at the sidebar switcher. Clears automatically when the user
  switches to Paste or Upload so the banner never implies a
  user-provided file came from the literature.
- **Onboarding banner updated** to list all four built-in
  scenarios (Writing essay / Large-scale writing / L2
  speaking / Clinical OSCE) with observation counts, so the
  quickstart flow advertises the scenarios before the user even
  reaches the sidebar.

### Unchanged on purpose

- The "Run with sample data" onboarding quickstart button still
  fires the default writing-essay scenario — first-timers still get
  a small dataset before they meet the switcher.
- All four scenarios, their parameter values, and their APA 7
  references are identical to v0.2.7. Only the UI wrapping changed.

## 0.2.7-beta - 2026-04-17

Sample-data enrichment. Replaces the single hardcoded 960-observation
demo with a registry of four literature-grounded synthetic scenarios,
each sized and parameterised for a different diagnostic emphasis.
Users who want to see residual PCA, bias heatmaps, or the rater-
severity Wright-map callout at publication-realistic scale can now
pick a scenario that actually exercises those paths.

### Added

- **Four built-in sample-data scenarios** selectable from a new
  `Sample scenario` sidebar selectbox:
  - **Writing essay** (30×4×2×4, 960 obs) — the v0.1+ default,
    retained unchanged. (Eckes, 2011; McNamara, 1996; Linacre, 1989)
  - **Large-scale writing** (120×4×2×3, 2,880 obs) — sized for
    stable residual-PCA; one rater deliberately injected at +1.6
    logits of severity so the bias heatmap and misfit ranking have a
    legible signal. (Myford & Wolfe, 2003, 2004; Engelhard, 1994,
    2013; Smith, 2002)
  - **L2 speaking** (80×3×3×5, 3,600 obs) — analytic rubric with
    five criteria; Pronunciation is hardest by design. Tighter rater
    severity spread typical of trained panels.
    (McNamara, 1996; Bachman & Palmer, 1996; Luoma, 2004)
  - **Clinical OSCE** (60×4×5×3, 3,600 obs) — five stations ×
    three competencies on a compact 4-point scale; station difficulty
    is the dominant source of variation. (Downing & Yudkowsky, 2009;
    Tavakol & Dennick, 2011; Wolfe & Song, 2015)
- **"About this dataset" sidebar expander** renders the chosen
  scenario's description, design dimensions, observation count, and
  APA 7 reference list built from the same `_APA_REFERENCE_LIBRARY`
  the Publication Document uses.
- **Five new APA 7 references** added to support the new scenarios:
  Bachman & Palmer (1996), Downing & Yudkowsky (2009), Engelhard
  (2013), Luoma (2004), Tavakol & Dennick (2011). All integrated
  into `_CITATION_TO_KEY` so they are pickupable by the narrative
  citation scanner.
- **`_self_test_sample_data_scenarios`** (test #40) pins every
  scenario's row-count / score-bounds / citation resolution /
  determinism contract and enforces backward compatibility of the
  legacy `sample_mfrm_data()` signature.

### Changed

- `sample_mfrm_data(seed)` is now a thin wrapper over the new
  `sample_mfrm_data_by_key(DEFAULT_SAMPLE_SCENARIO_KEY, seed)`. The
  default scenario still produces the same 960-row writing-essay
  demo as v0.1–v0.2.6, so existing screenshots, saved configs, and
  onboarding tours continue to match.
- New module-level helper `_generate_mfrm_rsm_from_params(params, seed)`
  is the single RSM-family generator every scenario shares; each
  scenario only owns its parameter dict, not the sampling loop.
- Sample-data CSV download filename now encodes the chosen scenario
  (`mfrm_sample_<key>.csv`) so users can re-download each scenario
  without overwriting the previous one.

### Unchanged on purpose

- `cached_sample_mfrm_data(seed)` is retained with the same
  signature for existing call sites; `cached_sample_mfrm_data_by_key`
  is the new per-scenario cache.
- The "Run with sample data" onboarding banner still fires the
  default writing-essay scenario — first-time users get the familiar
  small dataset before they discover the selectbox.

## 0.2.6-beta - 2026-04-17

Progressive-disclosure polish + self-test expansion. Continues the
cognitive-load reduction started in v0.2.5: every major diagnostic
chart now has a just-in-time help popover, the Essential view-density
mode is extended to Help and Report (not just Visuals), and four new
contract-level self-tests pin previously manually-tested code paths.

### Added

- **Help popover library expanded from 13 → 18 topics.** New entries:
  `misfit_ranking`, `facet_distribution`, `rater_agreement`,
  `threshold_map`, `zstd_distribution`. Each follows the existing
  What it shows / How to read / Watch for schema.
- **Five more charts wired to popovers:** ZSTD distribution histogram,
  Top-misfit ranking, Facet element distribution, Step / Threshold
  Ordering, and Inter-Element Agreement. Combined with v0.2.5's
  Wright Map / Scree / Fit-scatter / Category-probability /
  Bias-heatmap wiring, every chart in the core results tabs now
  has a one-click "How to read" affordance next to its subheader.
- **Essential mode extended to Help and Report sub-tabs.** Help drops
  from 10 → 7 tabs in Essential (hides Rating Scale Guide, Model
  Capability, Public Beta). Report drops from 10 → 6 sub-tabs in
  Essential (hides Manuscript Template, Method Appendix, Facet
  Equivalence, Stan Code). Full mode (sidebar toggle) restores all
  sub-tabs for publication-depth work.
- **Four new inline self-tests** raise the total from 35 → 39:
  `_self_test_posterior_load_netcdf` (synthetic 2-chain × 100-draw
  round-trip via `az.from_dict` → `.to_netcdf()` → uploader wrapper,
  portable across arviz 0.x / 1.x),
  `_self_test_posterior_load_cmdstan_csvs` (two minimal CmdStan CSVs
  → `_posterior_load_cmdstan_csvs` → shape assertions),
  `_self_test_config_json_import_whitelist` (23-key frozenset +
  critical-key check + filter drop-through),
  `_self_test_run_history_clear_confirmation` (push N+2 runs, verify
  cap, verify clear + confirm flag independence).
- **`_self_test_essential_mode_tab_filters`** pins the exact
  Essential-visible label sets for the Help and Report tabs so a
  future tab rename or reorder surfaces as a failing test rather
  than a silently-hidden sub-tab.

### Changed

- `_CONFIG_JSON_IMPORT_WHITELIST` extracted from an inline `set` inside
  the Report tab's import handler to a module-level `frozenset` near
  `_RUN_HISTORY_KEY`, so the new whitelist test pins the same object
  the handler uses. No behavioural change — same 23 keys, same filter.
- `show_help_section` and `show_report_section` now build their tab
  lists dynamically based on the View density setting, with a
  `{label: tab}` dict replacing positional `help_tabs[4]` indexing so
  hiding a tab in the middle no longer shifts downstream indices.

### Unchanged on purpose

- All 10 Help tabs and 10 Report sub-tabs remain reachable in Full
  mode — Essential only hides the *default* view.
- Popover library topics from v0.2.5 are untouched; v0.2.6 only
  appends new entries.

## 0.2.5-beta - 2026-04-17

Cognitive-load hotfix. The v0.2.0 → v0.2.4 feature stack added 23
diagnostic visualisations and roughly 40 new widgets, tabs, and
expanders, crossing the threshold where the default view started to
overwhelm first-time users. This release layers three fixes on top
without removing any existing feature.

### Added

- **Just-in-time help popovers** backed by `_HELP_POPOVER_LIBRARY`
  with 13 three-section entries (What it shows / How to read /
  Watch for). The new helper `render_help_popover(topic_key)`
  mounts an `st.popover` button next to chart titles so users get
  focused guidance without leaving the current tab. Wired into the
  Pathway Map, Forest plot, Q-Q plot, ECDF, Observation-coverage
  heatmap, Category-usage bar, Posterior trace, and Posterior
  Rhat/ESS bars. On Streamlit < 1.32 the helper falls back to an
  expander so help is still reachable.
- **View density toggle** in the sidebar (radio: Essential / Full,
  default Essential). Essential hides the three advanced Visuals
  sub-tabs (Forest / Q-Q / ECDF) and shows a compact caption
  pointing users at the Full mode when they need publication-depth
  diagnostics. State survives reruns via
  `st.session_state["app_view_density"]`.
- Inline self-test `_self_test_help_popover_library` (#34) pins the
  13 required topic keys and enforces minimum prose length so the
  library cannot silently empty out.

### Changed

- **First-read guide is collapsed by default** (was unfolded on
  every rerun). Returning users see the result tabs immediately;
  first-time users still see the affordance at the top.
- **Visual interpretation roadmap expander** in the Visuals tab
  defaults to collapsed (was `expanded=True`), so the four core
  diagnostic sub-tabs are the first thing users see.

### Unchanged on purpose

- Readiness panel in Warning / Issue state still renders expanded
  (alerts must stay visible).
- Anchor-issues, submission action plan "first actions", and
  copy-edit wording repair expanders stay expanded (all three are
  actionable blockers before final reporting).

## 0.2.4-beta - 2026-04-17

Visualization-coverage hotfix. Six new diagnostic plots land across
the Data tab, Visuals tab, and Posterior Viewer, closing the gaps
identified by the post-v0.2.3 chart-coverage audit. All new plots
are Plotly-based, respect the v0.2.2 WCAG palette (primary teal
`#0d7a5a`, warn `#c24e00`, neutral `#666666`), and degrade to an
st.info notice when their required data is missing rather than
rendering broken canvases.

### Added — Data tab

- **Observation coverage heatmap** (`_draw_data_coverage_heatmap`).
  Person × first-facet grid, blue cells = ≥ 1 observation, white =
  none. Capped at 80 persons for readability; fires a warning when
  > 40 % of displayed cells are empty (connectivity risk).
- **Category usage bar chart** (`_draw_category_usage_bar`). Raw
  observed frequency per score category, with an optional "split
  by facet" stack. Categories below `max(5, 1%)` of the total are
  called out in a warning so under-used rating categories are
  caught before they cascade into disordered thresholds downstream
  (Linacre 2004).

### Added — Visuals tab (three new sub-tabs)

- **Forest plot** (`_draw_measures_forest_plotly`). Per-facet
  forest of element measures with dot = Estimate, thick bar =
  50 % CI (± 0.67 · SE), thin bar = 95 % CI (± 1.96 · SE). The
  frequentist analogue of the Bayesian forest already in
  Posterior Viewer; covers JMLE / MML users who do not cross into
  Stan.
- **Q-Q plot** (`_draw_residual_qq_plotly`). Standardized residuals
  vs N(0, 1) theoretical quantiles with 45° reference. Requires
  ≥ 20 residuals. Heavy tails flag extreme-score misfit; S-shape
  flags possible residual structure.
- **ECDF of measures** (`_draw_measure_ecdf_plotly`). Empirical
  cumulative distribution of person and facet-element measures on
  a shared logit axis. Flat stretches = measurement gaps; steep
  jumps = clusters. Complements the Wright Map binned histogram.

### Added — Posterior Viewer

- **Rhat / ESS bar chart tab** (`_posterior_rhat_ess_figure`).
  Two-panel horizontal bar chart with dashed reference at
  Rhat = 1.01 and ESS = 400. Bars turn red (#c24e00) when either
  threshold is crossed so convergence issues jump out without
  reading the full summary table row by row. Reuses
  `_posterior_compute_summary` so numbers stay in sync.

## 0.2.3-beta - 2026-04-17

Third-pass hotfix landing the remaining findings from a second,
post-v0.2.2 comprehensive gap audit. Focus: reconciling in-app text
and metadata tables with the feature surface that actually ships.

### Fixed

- `public_release_readiness_table` row "Simulation validation
  templates" had Status predicate reading from a 4-key scripts
  dict while Evidence sentence quoted a 3-row inventory. Aligned
  both on the inventory so "X rows are documented" and the
  Ready / Review flag cannot drift apart.
- `_render_reporting_checklist` two cells were hardcoded to
  unchecked (PCA of residuals, Average measures per category)
  regardless of whether the diagnostic was actually available.
  Wire them to the real data.
- 11 stale tab-reference strings updated to match the v0.2.0 Report
  meta-grouping: `Report → Stan Code` → `Report → Exports → Stan
  Code`, `Report → Tables` → `Report → Tables & checks → Tables`,
  `Report → APA Report` → `Report → Reports → APA Report`, etc.

### Added

- Help tab Quick Start + Troubleshooting documenting every v0.2.x
  feature and diagnostic (shipped in a separate commit on the same
  release branch).
- `README.md` "What's new in v0.2.0 → v0.2.2" covers the four
  feature tracks plus both post-release hotfix passes.
- `SECURITY.md` Posterior Viewer Uploads and Advanced-Model Stan
  Code Downloads sections document the two v0.2.x file-upload /
  file-download surfaces distinct from rating-data handling.
- `public_beta_limitations_table` rows for Advanced models (Stan,
  download-only), Posterior Viewer (upload), and Publication
  Document (Word / PDF / HTML).
- `mfrmr_015_migration_coverage_table` rows for Advanced response
  models (Stan generators only) and Posterior Viewer (diagnostic).
- `visual_method_evidence_table` row for Posterior trace / ridge /
  pair / forest diagnostics with primary references.
- `external_reference_documentation_table` rows for arviz and
  cmdstan / cmdstanpy / rstan so reviewers have explicit hooks.
- Dimensionality readiness row now concatenates
  `diagnostics['pca_reason']` into its Evidence sentence, so the
  readiness table, manuscript claim guide, publication gate, and
  submission action plan all explain WHY PCA was skipped (small
  data / single item-combination / disabled by Analysis depth /
  >95% NaN).
- Onboarding-banner tip line carries a one-line v0.2.x feature
  catalogue (Publication Document / Posterior Viewer / Advanced
  models / Run history) so new users discover them before clicking
  "Got it".
- Self-tests `_self_test_publication_document_html` (#32) and
  `_self_test_diagnose_estimation_error_patterns` (#33).

## 0.2.2-beta - 2026-04-17

Second-pass hotfix completing the remaining MAJOR and MINOR items
from the pre-redeploy 8-agent UX audit.

### Added

- Run breadcrumb banner (`_render_run_breadcrumb`) above the critical-
  warnings row: one compact line showing model / method / convergence
  status / iterations / obs × person × facet counts / 8-char run
  fingerprint. Users flipping between history snapshots finally have
  a persistent "where am I?" signal.
- MFRM / Bayesian quick-reference glossary (`_MFRM_GLOSSARY` +
  `render_glossary_expander`) with 20 curated entries (logit, Infit,
  Outfit, MnSq, ZSTD, step facet, anchor, separation, strata, RSM,
  PCM, GPCM, SE, CI, Rhat, ESS, divergent, E-BFMI, …). Wired into
  the existing Help → Glossary tab so both the quick reference and
  the deep 50+ term glossary are one click away.
- `reorder_measure_columns(df)` + module-level priority list so every
  measure table surfaces Facet / Level / Estimate / SE / CI /
  Infit / Outfit first; Anchor, N, ReliabilityNote, etc. move to the
  tail. Applied to the Combined measures + fit table in the Measures
  tab; callers can opt in one table at a time.
- `style_fit_columns(df)` returns a pandas Styler that colour-codes
  Infit / Outfit mean-square cells by the Wright & Linacre (1994)
  bands (0.5–1.5 acceptable · 1.5–2.0 noisy · ≥ 2.0
  distorting · < 0.5 over-fit). ZSTD columns are excluded (they
  use a different interpretation scale). Applied to the Combined
  measures + fit table.
- Scree plot now draws three reference lines — EV = 1 (expected
  null), EV = 2 (caution), EV = 3 (strong secondary residual signal) —
  matching the thresholds in the chart guide and narrative.
  Previously the EV = 3 threshold was documented but invisible on
  the plot.
- Inline self-tests: `_self_test_style_fit_columns` (#29),
  `_self_test_mfrm_glossary` (#30), `_self_test_reorder_measure_columns`
  (#31).

### Changed

- Stan downloads now advertise `application/x-stan` instead of
  `text/plain`; browsers trigger a save dialog instead of inline
  preview, and Stan-aware editors pick up syntax highlighting.
- Scree plot neutral-threshold line shifted from `#999999` to the
  stronger `#666666` grey to satisfy WCAG AA contrast on white
  (was 2.85 : 1, now 5.7 : 1).
- Posterior Viewer sidebar carries a reassurance caption when a
  FACETS-mode run is in session_state, so users know switching to
  the viewer does not discard their live run.

### Fixed

- `reorder_measure_columns` short-circuited on `df.empty` (True for
  any rowless DataFrame), leaving placeholder / pre-estimation
  tables with the wrong column order. Now guarded on
  `len(df.columns) == 0` so only genuinely schema-less frames
  pass through untouched.

## 0.2.1-beta - 2026-04-17

Pre-redeploy hotfix landing the 8 highest-impact findings from a
parallel UI / UX audit (downloads, tables, charts, information
architecture, microcopy, onboarding, accessibility, performance).

### Added

- Downloadable built-in sample dataset (Download sample CSV) in the
  sidebar when "Sample data (built-in)" is selected, so users can
  inspect the exact column structure the estimator expects before
  preparing their own upload.
- File-size preflight on CSV uploads: ≥50 MB shows a st.warning,
  ≥100 MB shows a st.error, before pandas parsing consumes memory.
- Clear-history confirmation: the Clear history button now uses a
  two-step flow (first click surfaces a warning with explicit
  "Yes, delete all N snapshots" / "Cancel" buttons), preventing
  accidental one-click loss of the run stack.
- Publication Document signpost at the top of the Downloads tab, so
  users expecting a Word/PDF manuscript file are pointed at
  Report → Exports → Publication Document instead of bouncing.
- Read-only companion import for `mfrm_config.json`: upload a
  previously downloaded config file to inspect its settings as a
  table. Does not overwrite sidebar state (most fields are
  run-specific metadata rather than replayable inputs).
- Canonical rounding helper `format_measure_table()` with `3` decimals
  for measure-like columns (Estimate / SE / CI / Infit / Outfit /
  RMSE / Separation / Eigenvalue …) and `4` decimals for probability-
  like columns (Probability / Proportion / Percent / Share / Rate).
- Inline self-tests `_self_test_format_measure_table` (case #28).

### Changed

- `render_chart_guide()` (v0.2.0-beta) was dead code — the
  11-entry library was registered but never rendered. Now called
  after the scree plot, yardstick Wright map, category probability
  curves, pathway map, facet distribution, posterior trace, and
  posterior ridge, so every diagnostic plot carries a consistent
  "How to read this" expander backed by the library.
- Readiness-panel states render as textual [OK] / [CAUTION] / [ISSUE]
  labels at every render site, so users can distinguish severity
  without hue alone.
- Run-history snapshot cap lowered from 10 → 5 to cut
  session-state memory pressure on Streamlit Community Cloud
  (each deep-copy of facets_mode_output can reach 50+ MB on
  realistic datasets).
- Four most prominent empty-state messages in the Measures tab
  ("No person measures available" → "No person measures yet. Click
  Run FACETS-mode estimation …") rewritten to be action-oriented
  and explain when the table becomes available.

### Accessibility

- Global CSS additions via `_inject_desktop_readability_css`:
  - `:focus-visible` outline 3 px #0066cc at 2 px offset for
    keyboard users (WCAG 2.4.7 Focus Visible).
  - `prefers-reduced-motion` suppresses all animations /
    transitions / smooth scrolling for vestibular-sensitive users.
  - `@media (max-width: 899px)` wraps tabs and tightens the
    container padding so the wide layout does not crush tablet /
    iPad content.

## 0.2.0-beta - 2026-04-17

### Added — advanced model Stan code generators (download-only)

- Seven advanced response models registered in `_ADVANCED_RESPONSE_MODELS`:
  `DINA` (CDM; de la Torre 2009), `HRM` (Patz et al. 2002), `TESTLET_RI`
  and `TESTLET_BIFACTOR` (Bradlow 1999 / DeMars 2006), `MIXTURE_RASCH`
  (Rost 1990), `IRT_2PL_BINARY` (Birnbaum 1968), `PAIRWISE_BTL`
  (Bradley & Terry 1952). Each entry carries `family`, `binary`,
  `needs_q_matrix`, `needs_testlet_column`, and `needs_class_count`
  metadata so the sidebar can surface the right extra widgets.
- Per-model Stan code generators:
  - `generate_dina_stan_code(n_items, n_attributes)` — Bernoulli
    likelihood with slip / guess Beta(2, 10) priors, enumerated
    profile class, `log_lik`, and `profile_class` in
    `generated quantities`.
  - `generate_hrm_stan_code(n_categories)` — signal-detection
    parameters phi (accuracy, lognormal) and eta (bias), ordered
    kappa thresholds, `y_rep` posterior predictive.
  - `generate_testlet_stan_code(n_categories, bifactor)` — random
    intercept by default; bifactor adds a `theta_testlet_general`
    latent factor per testlet.
  - `generate_mixture_rasch_stan_code(n_classes)` — latent classes
    with class-specific item difficulties, Dirichlet prior on class
    probabilities.
  - `generate_2pl_binary_stan_code()` — free discrimination
    (lognormal) + difficulty, `y_rep` PPC.
  - `generate_pairwise_btl_stan_code()` — ability vector + pairwise
    Bernoulli logit.
  - `generate_advanced_model_stan_code(name, ...)` — dispatcher.
- `validate_q_matrix(df)` validates DINA Q-matrices (0/1 values,
  attribute coverage, item coverage, per-item / per-attribute row
  and column sums) and returns a structured verdict with messages.
- New sidebar section "Advanced models (Stan, download only)"
  (collapsed by default) with:
  - enable checkbox
  - model family picker (labels from the registry)
  - Q-matrix uploader (DINA), class-count input (Mixture Rasch)
  - "Generate Stan code" button that serialises the Stan program
    into session_state and exposes a Download .stan button.
- Inline self-test `_self_test_advanced_model_generators` walks every
  registered advanced model, verifies that the generated Stan code
  contains each of the four canonical blocks, has balanced braces,
  carries model-specific keyword markers (bernoulli_lpmf for DINA,
  phi for HRM, gamma_testlet for testlets, simplex + dirichlet for
  mixture, alpha_item for 2PL, ability for BTL), and exercises
  `validate_q_matrix` on well-formed and broken inputs.

### Added — UX tweaks

- Toast notification (`st.toast`) fires in addition to the persistent
  `st.status` accordion whenever estimation completes; message variants
  communicate convergence, non-convergence, and failure.
  Users who scroll away from the sidebar still get a completion signal.
- Report tab regrouped into three meta-categories:
  Reports (APA Report, Manuscript Template, Method Appendix, Claim Guide),
  Tables & checks (Tables, Reporting Checklist, Facet Equivalence, Readiness),
  Exports (Stan Code, Publication Document).
  Individual sub-section renderers are unchanged; only the surface shape is.
- Unified chart interpretation helper `render_chart_guide(chart_name)`
  backed by `_CHART_GUIDE_LIBRARY`. Every diagnostic plot can now drop
  a consistent "How to read this" expander with headline + body text
  sourced from a single library (Wright map, pathway map, category
  probability curves, threshold map, ICC, scree, facet distribution,
  reliability, rater agreement, posterior trace, posterior ridge).
- Keyboard shortcut cheat sheet (`render_keyboard_shortcuts_help`) lives
  in a collapsed sidebar expander so the shortcut surface is documented
  once and kept up to date (R rerun, C clear cache, Esc close, ? cheat
  sheet, Tab focus, Enter activate, Ctrl/Cmd+F search).
- Cache-stale detection banner now exposes a one-click **Rerun now**
  button that sets `_facets_mode_force_rerun` and triggers
  `st.rerun()`, so users don't have to scroll back to the sidebar.
- Inline self-test `_self_test_chart_guide_library` pins the chart-guide
  library key set and enforces a minimum text budget per entry so the
  explanatory text doesn't silently empty out over time.

### Added — Posterior Viewer

- New top-level app mode "Posterior Viewer (upload)" exposed via the
  sidebar radio, for inspecting externally-produced posterior draws
  without leaving the browser. Estimation itself still happens locally
  (Streamlit Cloud is used only as a visualisation frontend).
- Upload pipeline supports three formats with auto-detection by file
  extension:
  - CmdStan CSV (one or more per-chain outputs, loaded via
    `arviz.from_cmdstan`).
  - Apache Parquet (any DataFrame with `chain`, `draw`, and one
    column per parameter).
  - ArviZ NetCDF `.nc` (round-tripped via `arviz.from_netcdf`).
- Summary table with mean / sd / median / 5% / 95% / n_eff / Rhat
  (the ESS and Rhat columns come from `arviz.ess` / `arviz.rhat`) plus
  a one-click CSV download.
- HMC transition diagnostics banner with coloured metrics: divergent
  transitions count + %, max-treedepth hits, acceptance rate mean,
  step-size mean, and E-BFMI min across chains. Inline guidance triggers
  when thresholds are crossed (e.g. `adapt_delta 0.95 → 0.99` for
  divergences, `max_treedepth 10 → 12` for treedepth hits, or
  re-parameterisation for E-BFMI < 0.3).
- Plot suite (Plotly, matches the rest of the app):
  - Trace — chain-coloured trace per parameter
  - Ridge — offset KDE ridge across selected parameters
  - Pair  — scatter-matrix with divergence highlights
  - Forest — posterior mean + 50% and 95% CIs
- Parameter multi-select (defaults to the first 6 parameters) drives
  every plot and the summary table.
- Inline self-test `_self_test_posterior_viewer_loaders` round-trips
  synthetic parquet draws through the loader, summary, HMC diagnostics,
  and plot builders. Gracefully skipped when `arviz` / `pyarrow` are
  not yet installed.

### Added — publication export

- Word (.docx) publication document builder with auto-generated abstract,
  exhaustive Methods (estimator / priors / anchor policy / convergence
  criteria / unused-category handling), manuscript-style Results,
  embedded figures (Wright map, fit scatter, category probability curves,
  facet distribution), core tables (element measures, reliability), and
  an APA 7 reference list. Runs lazily via `python-docx` (new dependency).
- PDF publication document builder with the same structure rendered as a
  letter-size reportlab platypus flowable stack (hanging-indent
  references, 6-inch figure embeds, 8 pt table styling). The PDF
  excludes the reproducibility code appendix per the monolith's public
  beta convention; use the Word export for code-bearing variants.
- HTML publication document builder (`build_publication_html_bytes`):
  self-contained single-file HTML with inline CSS, base64 figure embeds,
  and a minimal markdown → HTML converter. Works offline with no external
  assets.
- APA 7 reference-list generator (`collect_cited_references`,
  `build_apa_reference_list`) keyed by the narrative's `(Author, Year)`
  citations, with an always-include anchor set so every publication
  document carries the foundational Rasch references (Andrich 1978,
  Masters 1982, Linacre 1989 / 2024, Myford & Wolfe 2003 / 2004,
  Smith 2002, Wright & Masters 1982).
- New Report tab sub-tab "Publication Document" exposing three
  download paths (Word / PDF / HTML) from identical narrative sources.
- Inline self-tests: `_self_test_apa_reference_list`,
  `_self_test_publication_document_word`,
  `_self_test_publication_document_pdf`.

### Added — PCA + bias/reliability silent-fail diagnostics

- `diagnose_pca_skip_reason` surfaces exactly why residual PCA was
  skipped (small person count, too few item-combinations, all-NaN
  columns, NaN share above 95%, exception, or an Analysis depth that
  disabled PCA) instead of the previous silent "PCA is not available"
  fallback. The Dimensionality tab reads the stored reason and shows
  it with concrete fix suggestions.
- `estimate_bias_interaction` now returns a `_skip_reason`-carrying
  dict on every early-return path (missing result / diagnostics, empty
  observation table, same-facet pair, unknown facet names, no valid
  cells after filtering). `show_bias_section` surfaces the reason
  instead of crashing on the missing `table` key.
- `calc_reliability` records a `ReliabilityNote` column explaining
  why Separation / Reliability is NaN for a given facet ("only N
  level", "all levels have identical measures", "SE zero or NaN",
  "measurement variance ≤ error variance").

### Added — UX layer (st.status, error remedy, readiness, history, comparison)

- Replace the plain `st.spinner` with an expandable `st.status`
  accordion that writes step-by-step progress (Estimate → Diagnostics
  → Report tables → Bias interactions) and the final elapsed seconds.
- Replace the generic estimation "common causes" checklist with a
  pattern-matched diagnosis → action guide
  (`diagnose_estimation_error`): singular matrix → non-centered facet,
  maxit reached → `maxit=1000`, rating scale error → Score column
  check, memory error → reduce facet cardinality, anchor error →
  fix anchor tables, etc. Falls back to the generic list for
  unrecognised failures.
- Pre-estimation data-quality readiness panel above the Run button
  (`build_readiness_report`, `render_readiness_panel`) with a
  [OK] / [CAUTION] / [ISSUE] readiness summary across eight checks (observation
  count, person count, score column dtype, facet count, per-facet
  level count, per-person coverage, column-role overlap).
- Dismissible onboarding banner with a "Run with sample data"
  one-click quickstart that fires estimation with the built-in
  sample dataset and default column mapping.
- Session-scoped run history (`record_run_in_history`,
  `render_run_history_panel`) caps the last 10 runs and offers a
  Restore button per entry plus a "Clear history" control.
- Two-run comparison panel (`render_comparison_panel`,
  `render_comparison_selector`) on top of the history stack:
  convergence / iterations / elapsed delta, element-level Pearson r
  + RMSE scatter with a 45° reference line, facet-level reliability
  diff table, and an automatic interpretation ribbon
  (near-perfect / strong / moderate / weak agreement).

### Dependencies

- Added `python-docx >= 1.0` (Word publication export).
- Added `reportlab >= 4.0` (PDF publication export).
- Added `arviz >= 0.17` (Posterior Viewer — Rhat / ESS / InferenceData I/O).
- Added `netcdf4 >= 1.6` (Posterior Viewer — NetCDF file round-trip).
- Added `pyarrow >= 15.0` (Posterior Viewer — Apache Parquet draws).

### Changed

- App version label: `0.1.2-beta` → `0.2.0-beta`.
- Report tab restructured into three meta-categories (Reports /
  Tables & checks / Exports); individual sub-section renderers
  unchanged.
- Residual-PCA skip-path surfaces a structured reason instead of the
  generic "not available" fallback.
- Estimation failure UI uses pattern-matched specific remedies instead
  of a generic four-bullet checklist.

## Unreleased

No unreleased changes yet.

## 0.1.2-beta - 2026-04-12

### Added

- Guided visual interpretation checklist in the app and downloadable report bundle.
- PCM/GPCM category probability curves can now be inspected by selected step-facet level, not only as averaged curves.
- Category probability curve exports now use the same curve builder as the Visuals tab and include long-form data for all curve scopes.
- Synthetic guided demo report export through `--export-demo-report`, including report tables, method appendix, and interactive category curves.
- Visual method evidence table plus readability safeguards for dense Wright maps, yardsticks, marginal heatmaps, and bias heatmaps.
- Latent-regression covariate type preview and export table for numeric, categorical, and integer-code review decisions.
- Public-beta limitation and release-readiness tables in the app, demo report export, and CLI release check.
- mfrmr 0.1.5 migration coverage table plus external validation artifact checklist, official-reference table, and reviewer template in the parity fixture.
- README preview screenshot generated from the real Streamlit app using built-in synthetic sample data.
- Archived local R cross-check smoke status for generated TAM, mirt, and sirt handoff fixtures.
- Result-aware manuscript claim guide in the Report tab, table downloads, and demo report export.
- Result-aware Markdown manuscript template for Methods, Results, limitations, reviewer preflight checks, and OSF/demo report exports.
- Publication gate summary that aligns APA Report conclusions with readiness checks and manuscript claim guardrails.
- Case-specific interpretation guidance for sparse categories, dimensionality, bias screens, MML marginal checks, rater reliability, and linking claims.
- Result-specific avoid/safer wording repairs in APA Report guidance, manuscript templates, and case-interpretation guidance exports.
- Submission action plan that combines publication gates, readiness checks, claim guardrails, and wording repairs into a prioritized first-read table.
- Desktop readability refinements for result-tab wrapping, wrapped guide tables, and dense or tightly spaced Wright map / yardstick labels.
- Publication-styled figure export bundle with 300 DPI PNG target, matching HTML figures, and a figure-use manifest.
- Archived simulation-validation inventory for external numerical validation planning without bundling private data or machine-specific absolute paths.
- Sanitized Python/R/Julia external Simulation validation templates, downloadable from the app, OSF/demo report exports, and parity fixture exports.

### Changed

- App version label now marks this public beta release as `0.1.2-beta`.
- Plotly/Kaleido floors now target Kaleido 1.x PNG export behavior, with HTML figure fallback when Chrome/Chromium is unavailable.
- Final-report readiness, first-read guidance, and report text now use the same residual, reliability, and PCA thresholds.
- APA Report summary wording now respects claim-guide caveats and avoids broad no-bias claims outside computed screens.
- New/held-out prediction session-state invalidation now fingerprints uploaded file content, not only file name and size.
- Uploaded CSV/TSV reads are now non-destructive, so sidebar previews do not consume the file stream before estimation.
- `make verify` now includes `python streamlit_app.py --release-check`.
- The optional R cross-check scaffold now records package versions and output file manifests for archived validation evidence.

### Validation

- `make verify` passed with 21 built-in self-tests, 10 pytest checks, release check, benchmark smoke, demo report export, and parity fixture export.
- GitHub Actions passed on Python 3.11 and 3.12 for the final pre-release template commit before version finalization.

## 0.1.1-beta - 2026-04-12

### Added

- Deployment guide for local and Streamlit Community Cloud use.
- Cross-package validation contract for TAM, sirt, mirt, FACETS-like, and ConQuest-style comparisons.

### Validation

- `make verify` passed with 18 built-in self-tests, 2 pytest checks, benchmark smoke, and parity fixture export.
- GitHub Actions passed on Python 3.11 and 3.12 for the validation-protocol commit.

## 0.1.0-beta - 2026-04-12

Initial public-beta repository setup.

### Added

- Standalone single-file Streamlit app for Many-Facet Rasch Model workflows.
- RSM, PCM, bounded GPCM, JMLE, MML EM/Direct/Hybrid/Auto, latent regression, plausible values, prediction, simulation, strict marginal diagnostics, residual PCA, anchor audit, and linking review workflows.
- Data privacy warning in the app and README.
- `st.cache_data` / `st.cache_resource` usage for built-in sample data, static bundled assets, read-only core namespace, and bounded export-byte caches.
- Reproducibility/config fingerprints in exported analysis metadata.
- `--doctor`, `--self-test`, `--benchmark-quick`, and `--export-parity-fixture` CLI paths.
- CI smoke workflow for Python 3.11 and 3.12.
- Validation scaffold for optional TAM, mirt, and sirt external checks.
- Anchor templates and a user guideline.
- MIT license file.

### Validation

- Local compile, doctor, self-test, AppTest, benchmark smoke, parity fixture export, optional R scaffold, fresh virtual environment install/self-test, and Streamlit health checks were run during repository preparation.

### Caveats

- Release status is public beta / research preview.
- Exact numerical parity with FACETS, TAM, sirt, mirt, or `mfrmr` is not claimed.
- Licensed under MIT.
