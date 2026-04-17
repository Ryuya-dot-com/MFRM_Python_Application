# Changelog

All notable changes to this standalone Streamlit distribution should be recorded here.

## Unreleased

### Added — Phase D (UX tweaks)

- Toast notification (`st.toast`) fires in addition to the persistent
  `st.status` accordion whenever estimation completes — ✅ / ⚠️ / ❌
  variants communicate convergence, non-convergence, and failure.
  Users who scroll away from the sidebar still get a completion signal.
- Report tab regrouped into three meta-categories:
  📝 Reports (APA Report, Manuscript Template, Method Appendix, Claim Guide),
  📊 Tables & checks (Tables, Reporting Checklist, Facet Equivalence, Readiness),
  💾 Exports (Stan Code, Publication Document).
  Individual sub-section renderers are unchanged; only the surface shape is.
- Unified chart interpretation helper `render_chart_guide(chart_name)`
  backed by `_CHART_GUIDE_LIBRARY`. Every diagnostic plot can now drop
  a consistent ❓ "How to read this" expander with headline + body text
  sourced from a single library (Wright map, pathway map, category
  probability curves, threshold map, ICC, scree, facet distribution,
  reliability, rater agreement, posterior trace, posterior ridge).
- Keyboard shortcut cheat sheet (`render_keyboard_shortcuts_help`) lives
  in a collapsed sidebar expander so the shortcut surface is documented
  once and kept up to date (R rerun, C clear cache, Esc close, ? cheat
  sheet, Tab focus, Enter activate, Ctrl/Cmd+F search).
- Cache-stale detection banner now exposes a one-click 🔁 **Rerun now**
  button that sets `_facets_mode_force_rerun` and triggers
  `st.rerun()`, so users don't have to scroll back to the sidebar.
- Inline self-test `_self_test_chart_guide_library` pins the chart-guide
  library key set and enforces a minimum text budget per entry so the
  explanatory text doesn't silently empty out over time.

### Added — Phase B (Posterior Viewer)

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
  - 📈 Trace — chain-coloured trace per parameter
  - 🏔 Ridge — offset KDE ridge across selected parameters
  - 🔗 Pair  — scatter-matrix with divergence highlights
  - 🌲 Forest — posterior mean + 50% and 95% CIs
- Parameter multi-select (defaults to the first 6 parameters) drives
  every plot and the summary table.
- Inline self-test `_self_test_posterior_viewer_loaders` round-trips
  synthetic parquet draws through the loader, summary, HMC diagnostics,
  and plot builders. Gracefully skipped when `arviz` / `pyarrow` are
  not yet installed.

### Added — Phase A (publication export)

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
- New Report tab sub-tab "📄 Publication Document" exposing three
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
  🟢 / 🟡 / 🔴 traffic-light summary across eight checks (observation
  count, person count, score column dtype, facet count, per-facet
  level count, per-person coverage, column-role overlap).
- Dismissible onboarding banner with a "🎯 Run with sample data"
  one-click quickstart that fires estimation with the built-in
  sample dataset and default column mapping.
- Session-scoped run history (`record_run_in_history`,
  `render_run_history_panel`) caps the last 10 runs and offers a
  Restore button per entry plus a "🗑 Clear history" control.
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

## 0.1.2-beta - 2026-04-12

### Added

- Beginner-facing visual interpretation checklist in the app and downloadable report bundle.
- PCM/GPCM category probability curves can now be inspected by selected step-facet level, not only as averaged curves.
- Category probability curve exports now use the same curve builder as the Visuals tab and include long-form data for all curve scopes.
- Synthetic beginner-facing demo report export through `--export-demo-report`, including report tables, method appendix, and interactive category curves.
- Visual method evidence table plus readability safeguards for dense Wright maps, yardsticks, marginal heatmaps, and bias heatmaps.
- Latent-regression covariate type preview and export table for numeric, categorical, and integer-code review decisions.
- Public-beta limitation and release-readiness tables in the app, demo report export, and CLI release check.
- mfrmr 0.1.5 migration coverage table plus external validation artifact checklist, official-reference table, and reviewer template in the parity fixture.
- README preview screenshot generated from the real Streamlit app using built-in synthetic sample data.
- Archived local R cross-check smoke status for generated TAM, mirt, and sirt handoff fixtures.
- Result-aware manuscript claim guide in the Report tab, table downloads, and demo report export.
- Result-aware Markdown manuscript template for Methods, Results, limitations, reviewer preflight checks, and OSF/demo report exports.
- Publication gate summary that aligns APA Report conclusions with readiness checks and manuscript claim guardrails.
- Case-specific beginner guidance for sparse categories, dimensionality, bias screens, MML marginal checks, rater reliability, and linking claims.
- Result-specific avoid/safer wording repairs in APA Report guidance, manuscript templates, and beginner case-guide exports.
- Submission action plan that combines publication gates, readiness checks, claim guardrails, and wording repairs into a prioritized first-read table.
- Desktop readability refinements for result-tab wrapping, wrapped guide tables, and dense or tightly spaced Wright map / yardstick labels.
- Publication-styled figure export bundle with 300 DPI PNG target, matching HTML figures, and a figure-use manifest.
- Local Simulation-directory validation inventory for external numerical validation planning without bundling private data or local absolute paths.
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
