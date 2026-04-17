# Changelog

All notable changes to this standalone Streamlit distribution should be recorded here.

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
  emoji + observation count) followed by Paste CSV/TSV and Upload
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
  scenarios (✏️ Writing essay / 📚 Large-scale writing / 🎙️ L2
  speaking / 🏥 Clinical OSCE) with observation counts, so the
  quickstart flow advertises the scenarios before the user even
  reaches the sidebar.

### Unchanged on purpose

- The "🎯 Run with sample data" onboarding quickstart button still
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
  - ✏️ **Writing essay** (30×4×2×4, 960 obs) — the v0.1+ default,
    retained unchanged. (Eckes, 2011; McNamara, 1996; Linacre, 1989)
  - 📚 **Large-scale writing** (120×4×2×3, 2,880 obs) — sized for
    stable residual-PCA; one rater deliberately injected at +1.6
    logits of severity so the bias heatmap and misfit ranking have a
    legible signal. (Myford & Wolfe, 2003, 2004; Engelhard, 1994,
    2013; Smith, 2002)
  - 🎙️ **L2 speaking** (80×3×3×5, 3,600 obs) — analytic rubric with
    five criteria; Pronunciation is hardest by design. Tighter rater
    severity spread typical of trained panels.
    (McNamara, 1996; Bachman & Palmer, 1996; Luoma, 2004)
  - 🏥 **Clinical OSCE** (60×4×5×3, 3,600 obs) — five stations ×
    three competencies on a compact 4-point scale; station difficulty
    is the dominant source of variation. (Downing & Yudkowsky, 2009;
    Tavakol & Dennick, 2011; Wolfe & Song, 2015)
- **"📚 About this dataset" sidebar expander** renders the chosen
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
  has a one-click "❓ How to read" affordance next to its subheader.
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
  flags potential multidimensionality.
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
  meta-grouping: `Report → Stan Code` → `Report → 💾 Exports → Stan
  Code`, `Report → Tables` → `Report → 📊 Tables & checks → Tables`,
  `Report → APA Report` → `Report → 📝 Reports → APA Report`, etc.

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
  bands (🟩 0.5–1.5 acceptable · 🟧 1.5–2.0 noisy · 🟥 ≥ 2.0
  distorting · 🟨 < 0.5 over-fit). ZSTD columns are excluded (they
  use a different interpretation scale). Applied to the Combined
  measures + fit table.
- Scree plot now draws three reference lines — EV = 1 (expected
  null), EV = 2 (caution), EV = 3 (strong secondary dimension) —
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

- Downloadable built-in sample dataset (⬇ Download sample CSV) in the
  sidebar when "Sample data (built-in)" is selected, so users can
  inspect the exact column structure the estimator expects before
  preparing their own upload.
- File-size preflight on CSV uploads: ≥50 MB shows a st.warning,
  ≥100 MB shows a st.error, before pandas parsing consumes memory.
- Clear-history confirmation: the 🗑 Clear history button now uses a
  two-step flow (first click surfaces a warning with explicit
  "Yes, delete all N snapshots" / "Cancel" buttons), preventing
  accidental one-click loss of the run stack.
- Publication Document signpost at the top of the Downloads tab, so
  users expecting a Word/PDF manuscript file are pointed at
  Report → 💾 Exports → 📄 Publication Document instead of bouncing.
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

- `render_chart_guide()` (v0.2.0-beta Phase D) was dead code — the
  11-entry library was registered but never rendered. Now called
  after the scree plot, yardstick Wright map, category probability
  curves, pathway map, facet distribution, posterior trace, and
  posterior ridge, so every diagnostic plot carries a consistent
  ❓ "How to read this" expander backed by the library.
- Readiness-panel traffic lights pair their 🟢 / 🟡 / 🔴 icons with
  textual labels [OK] / [CAUTION] / [ISSUE] at every render site,
  so users with deuteranopia / protanopia / tritanopia can
  distinguish severity without hue alone.
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

### Added — Phase C (advanced model Stan code generators, download-only)

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
- New sidebar section "🧪 Advanced models (Stan, download only)"
  (collapsed by default) with:
  - enable checkbox
  - model family picker (labels from the registry)
  - Q-matrix uploader (DINA), class-count input (Mixture Rasch)
  - 📥 "Generate Stan code" button that serialises the Stan program
    into session_state and exposes a ⬇ Download .stan button.
- Inline self-test `_self_test_advanced_model_generators` walks every
  registered advanced model, verifies that the generated Stan code
  contains each of the four canonical blocks, has balanced braces,
  carries model-specific keyword markers (bernoulli_lpmf for DINA,
  phi for HRM, gamma_testlet for testlets, simplex + dirichlet for
  mixture, alpha_item for 2PL, ability for BTL), and exercises
  `validate_q_matrix` on well-formed and broken inputs.

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

### Changed

- App version label: `0.1.2-beta` → `0.2.0-beta`.
- Report tab restructured into three meta-categories (📝 Reports /
  📊 Tables & checks / 💾 Exports); individual sub-section renderers
  unchanged.
- Residual-PCA skip-path surfaces a structured reason instead of the
  generic "not available" fallback.
- Estimation failure UI uses pattern-matched specific remedies instead
  of a generic four-bullet checklist.

## Unreleased

No unreleased changes yet.

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
