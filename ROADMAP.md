# MFRM Python Application Roadmap

- Status: active
- Last updated: 2026-08-04
- Product boundary: standalone Python

## Product thesis

The application will become a design-to-decision workbench for many-facet
rating studies. Its primary value is not the number of estimators it exposes.
Its value is the ability to show:

1. where a rating design is fragile before data collection;
2. which model and data conditions support a reported result;
3. whether a conclusion survives reasonable sensitivity checks; and
4. what evidence-based action a researcher can take next; and
5. which intended interpretations and uses remain outside the evidence that
   the application can assess.

The intended workflow is:

```text
Plan -> Check -> Estimate -> Diagnose -> Stress -> Decide -> Archive
```

The active statistical-engine remediation lane is governed separately by the
self-contained
[`docs/statistical_engine_remediation_roadmap.html`](docs/statistical_engine_remediation_roadmap.html).
As of 2026-08-04, G0 output protection and G1 exact identified
parameterization are complete; S2 inference reconstruction is active.
Until its release-qualification gate is complete, identification, inference,
bias-screening, and validation work follows that document's dependency order
and temporary output-safety boundaries. New estimator scope must not bypass
those gates.

Every step must remain executable in Python. The Streamlit interface is the
guided research surface; the statistical and evidence contracts belong in
testable modules under `mfrm_app/`.

`StudyContext` is a first-class input to interpretation, not decorative report
metadata. The same numerical analysis may support different bounded
interpretations for different populations, constructs, purposes, and stakes;
those differences change interpretation and decision identities without
silently changing the fitted AnalysisID.

## Execution roadmap: August 2026 to July 2027

This section is the operational view of the detailed delivery tracks below.
Work is selected by milestone exit gate, not by the number of UI features that
can be added. The default priority order is: prevent unsafe interpretation,
shorten the primary journey, modularize the owning contract, then expand
statistical scope only when validation can travel with it.

### Baseline — consolidated on 2026-08-04

The current baseline includes the standalone-Python runtime and release
boundary, AnalysisIdentity and EvidenceRecord foundations, bilingual Help and
terminology contracts, privacy-safe problem routes, contextual Help return,
phase-aware UX presentation policy, native release checks, and a clean native
test profile. The adversarial UX baseline is recorded in
[`docs/ux_adversarial_audit.md`](docs/ux_adversarial_audit.md).

Baseline verification:

- `make apptest`: 1,068 passed, 62 compatibility tests deselected;
- `--self-test`: 56 passed;
- `--doctor` and `--release-check`: passed;
- default sample and real-data privacy paths: AppTest verified; and
- worktree and roadmap: consolidated into a reviewable commit.

### M1 — one authoritative workflow shell (0-30 days)

Outcome: the default journey exposes one primary next action per phase while
retaining every scientific caveat and expert detail.

Deliverables:

- define one widget-independent workflow-shell record for source, setup, run,
  first-read, evidence detail, and archive phases;
- remove remaining duplicate orientation surfaces from the default result
  route while keeping stable Help targets and compatibility projections;
- keep the authoritative result navigation visible during long-page scrolling,
  with keyboard focus, Help return, reduced motion, and narrow-viewport browser
  acceptance;
- measure initial and post-run render topology so new visible controls require
  an explicit information-architecture decision; and
- add a clean-checkout CI job that runs compile, doctor, release-check,
  self-test, native AppTests, and privacy/localization gates.

Exit gate:

- one primary action is visually exposed in every default phase;
- the authoritative result selector remains operable while a long selected
  section scrolls, without covering required content on a narrow viewport;
- the sample journey and real-data setup are keyboard-completable without a
  hidden required control;
- view/language/help transitions trigger no estimation and preserve
  AnalysisID; and
- no raw response row is visible by default.

Progress on 2026-08-04:

- introduced a widget-independent six-phase workflow-shell contract with one
  named primary-action owner per phase;
- removed the duplicate top-level first-read overview from the default result
  route while retaining its full detail in the First read section;
- captured initial and post-run AppTest control budgets so additions to the
  rendered topology require an explicit information-architecture change; and
- converted the Python workflow into an explicit clean-checkout matrix with
  named workflow/privacy/localization contracts and a final worktree check.

Navigation refinement on 2026-08-04:

- removed the sidebar keyboard-shortcut cheat sheet and its maintenance/cache
  commands; the app defines no custom hotkeys;
- placed Essential and All-panels result selectors in one sticky navigation
  dock that remains available during long-page scrolling; and
- kept the narrow layout compact with a touch-scrollable one-line section
  switcher instead of a multi-row fixed overlay.

First-run route refinement on 2026-08-04:

- added the pure `mfrm_app.guidance` five-node catalog and exhaustive reducer;
  it imports neither Streamlit nor pandas and has no claim-readiness field or
  estimator callback;
- replaced the initial configuration-heavy page with one optional three-route
  decision: learn with a sample, start with user data, or continue without the
  guide;
- connected data-role review, production-path estimation, a bounded formative
  interpretation, and a non-scientific learning checkpoint into one
  skippable/resumable/restartable sample route; and
- kept guide completion independent from AnalysisID and preserved the fitted
  identity across formative answers and completion.

The browser accessibility matrix remains the next M1 acceptance item; native
AppTest topology is the required fallback gate until that browser connection
is available.

Increment verification:

- `make apptest`: 1,117 passed, 62 compatibility tests deselected;
- complete five-step sample AppTest: production fit, overclaim retry,
  learning-only completion, and stable AnalysisID passed;
- guide state/isolation AppTests: focused landing, own-data route, Japanese
  rerender, exit/resume, and existing-workspace restoration passed;
- representative initial/post-run AppTests: two scenarios passed;
- All-panels Help return: panel state and AnalysisID preserved without refit;
- `--self-test`: 56 passed;
- `--doctor` and `--release-check`: passed; and
- benchmark and demo-report generation smokes: passed.

### M2 — scalable input and modular presentation (31-90 days)

Outcome: data entry and result navigation can grow without extending the
70,000-line entrypoint or lengthening one global control list.

Deliverables:

- replace the flat source list with source class -> scenario/source detail,
  preserving stable IDs and migrating existing session/help state;
- extract source ingestion, setup/readiness, workflow shell, and result-router
  adapters into focused `mfrm_app/` modules;
- keep computation, evidence, and presentation policies independently tested;
- add rendered-route locale coverage for Japanese and English, including
  errors and advanced settings; and
- establish privacy-safe UX telemetry schemas for phase reached, stable reason
  code, rerun count, and time to first evidence, without input values or IDs.

Exit gate:

- `streamlit_app.py` is smaller after extraction and owns no new pure business
  rule;
- adding a sample scenario changes a registry, not the global page topology;
- old stable source IDs and Help returns continue to work; and
- no telemetry event can contain uploaded values, labels, paths, free text, or
  fingerprints.

### M3 — coherent evidence and design resilience (3-6 months)

Outcome: the app explains tensions between diagnostics and identifies fragile
rating designs without turning either into an automatic verdict.

Deliverables:

- normalize fit, dimensionality, local dependence, category, bias,
  reliability, and design evidence into one coherence view;
- implement deterministic leave-one-level and concentrated-missingness design
  perturbations with stable reason codes;
- provide constrained repair candidates with burden and non-claim boundaries;
- generate a human-readable Analysis Brief plus deterministic JSON sidecar;
  and
- introduce a conclusion-level sensitivity ledger that retains failed and
  unsupported specifications in its denominator.

Exit gate:

- contradictory fixtures produce an explicit tension state;
- disconnected, bridge-dependent, sparse, and imbalanced fixtures reproduce
  known behavior;
- every recommended repair states what it addresses and what it cannot prove;
  and
- no report claim cites stale, unavailable, or differently identified
  evidence.

### M4 — calibrated native validation and beta decision (6-12 months)

Outcome: public claims are backed by frozen Python-native validation evidence,
and the project can make an explicit beta-exit decision.

Deliverables:

- freeze ADEMP-style RSM, PCM, and bounded-GPCM scenario contracts;
- report bias, RMSE, interval coverage, recovery, flag rates, numerical
  failures, and Monte Carlo uncertainty for applicable conditions;
- establish checked runtime and peak-memory baselines for supported Python and
  dependency versions;
- complete the retain/deprecate/remove decision for dormant compatibility
  generators, viewers, fixtures, and locale copy; and
- publish a machine-readable validation and release-scope summary used by the
  UI, documentation, and release gate.

Exit gate:

- every advertised diagnostic or estimator claim links to a named validation
  condition or is explicitly labelled uncalibrated;
- clean-environment release gates pass on the supported version matrix;
- compatibility-only code cannot re-enter the public journey or release gate;
  and
- the public-beta label is retained or removed by a documented evidence-based
  decision, never by schedule alone.

### Operating scorecard

Review these measures at each monthly roadmap checkpoint:

| Dimension | Measure | Direction |
|---|---|---|
| Journey clarity | Visible primary actions per phase | Exactly 1 |
| Safety | Raw response rows visible by default | 0 |
| State integrity | Presentation-only transitions that refit or change AnalysisID | 0 |
| Evidence integrity | User-facing claims without current EvidenceIDs and boundaries | 0 |
| Accessibility | Sample and real-data setup tasks passing the browser matrix | Increase to 100% |
| Architecture | New pure rules added to `streamlit_app.py` | 0 |
| Reliability | Native test, doctor, self-test, and release-check failures | 0 |
| Performance | Material regressions without reviewed benchmark evidence | 0 |

Items outside the current milestone may proceed only when they do not delay its
exit gate and do not create another workflow, evidence, or state authority.

## Product boundaries

### In scope

- Native Python RSM, PCM, and bounded GPCM workflows within their documented
  identification and estimation limits.
- Native JMLE and MML estimation, uncertainty audits, prediction, simulation,
  and parameter-recovery evaluation.
- Prospective rating-design diagnostics and post-collection design audits.
- Facet-aware design, severity/leniency, halo, local-dependence, and residual
  network diagnostics.
- Residual PCA, PCA stability, DIMTEST, model-fit, category, bias/interaction,
  covariance, and reliability evidence.
- Sensitivity analysis that records when a substantive conclusion changes.
- Guided and research-oriented interfaces backed by the same statistical
  contracts.
- A structured study-context, rating-design, and scoring-process record that
  distinguishes what the application computed from what a study intends to
  infer or do.
- Reproducible Python-native reports, tables, figures, settings, fingerprints,
  and claim boundaries.

### Explicitly out of scope

- EGA or automatic discovery of a confirmatory Q matrix.
- Runtime calls to TAM, Shiny, FACETS, ConQuest, `mfrmr`, `Rscript`, `rpy2`, or
  any other external estimation engine.
- External-engine job submission, result ingestion, or cross-engine ZIP/schema
  interoperability.
- External-engine control/data-script generation or posterior-result ingestion
  as part of the core product journey.
- A remote R worker or a hosted external-solver queue.
- Confirmatory multidimensional MFRM unless it is later implemented and
  validated as a native Python model through a separate scope decision.
- Claims that a single diagnostic proves model truth, unidimensionality,
  fairness, rater quality, or suitability for a high-stakes decision.

Existing legacy export or comparison surfaces are not part of this roadmap and
will not be expanded. They will be inventoried, isolated from the default user
journey and release gates, and given an explicit retain/deprecate/remove
decision before the standalone core leaves beta. They must not become
dependencies of the core workflow.
External publications and software may still be cited as methodological
references, but they are not product integrations or release gates.

Python-native report archives may remain as application-specific download
packaging. They are not cross-engine interchange formats, and the application
will not promise that another engine can execute or return them.

## Starting point

This roadmap extends existing Python work rather than proposing a replacement
application.

| Area | Current foundation | Principal gap addressed here |
|---|---|---|
| Estimation | Native RSM/PCM/bounded-GPCM JMLE and MML paths | Conclusion-level sensitivity and native calibration |
| Dimensionality screens | Residual PCA, leave-one-column/bootstrap stability, DIMTEST | Coherent interpretation when the screens disagree |
| Rating networks | Design, severity/leniency, and halo networks | Perturbation resilience, effective-N evidence, and repair candidates |
| Uncertainty | Covariance/SE audits and prior-SD sensitivity | One versioned sensitivity ledger across supported decisions |
| Simulation | Synthetic data and parameter-recovery utilities | Immutable scenario contracts, formal runner, failure accounting, and MCSE |
| Evidence and reports | Readiness, claim-to-evidence, publication, and visual outputs | One normalized evidence source for UI, prose, and exports |
| Architecture | A growing `mfrm_app/` helper package plus a large entrypoint | Tested computation modules with compatibility wrappers |
| UX | Guided essential view, research detail, bilingual scaffolding | A shorter decision journey and evidence-linked actions |

The adversarial baseline, delivered first corrections, and task-level UX
acceptance measures are recorded in
[`docs/ux_adversarial_audit.md`](docs/ux_adversarial_audit.md). Future UX work
should reduce competing orientation surfaces against that baseline rather than
adding another local guide.

The current network functions are therefore modularized and extended, not
reimplemented. Likewise, existing PCA, DIMTEST, sensitivity, and reporting
outputs become inputs to a common evidence contract rather than parallel new
versions.

## North-star outcomes

### Defensible design

Before collecting data, a user can compare candidate assignments and identify
disconnection, single-link dependence, overloaded raters, weak overlap, and
likely failure under missing ratings or rater withdrawal.

### Explainable evidence

After estimation, every surfaced conclusion can be traced to a named evidence
record with its scope, applicability, warning state, sensitivity state, and
recommended next action.

### Robust conclusions

The application reports whether important conclusions persist across
prespecified estimator, scale, missingness, category, and leave-one-level-out
checks. It does not silently convert a fragile point result into a firm claim.

### Reproducible decisions

The application preserves the data fingerprint, resolved settings, seeds,
numerical diagnostics, evidence state, claim limits, and user-visible action
record needed to reproduce the analysis in the same Python application.

## Evidence vocabulary

Scientific status uses two different state systems, which will be kept
separate from workflow and learning progress.

### Computation state

- `AVAILABLE`: the evidence was computed and passed its applicability checks.
- `CAUTION`: the evidence was computed but has a documented limitation.
- `HOLD`: a failed prerequisite prevents the intended interpretation.
- `NOT_ASSESSABLE`: the data or model do not support the calculation.

### Conclusion stability

- `STABLE`: the conclusion persists across all required sensitivity checks.
- `CONDITIONALLY_STABLE`: it persists within an explicitly bounded subset of
  supported conditions.
- `SENSITIVE`: at least one required, credible specification changes it.
- `NOT_ASSESSED`: the required sensitivity analysis was not run.

These states are not pass/fail labels for a person, rater, institution, or
assessment. A conclusion may be computationally available and still be
sensitive.

### Workflow, learning, and claim readiness

- `workflow_complete`: the required interface actions for a bounded route were
  performed, including a route that ends with a documented failed analysis.
- `learning_complete`: the learner supplied the defined evidence of
  understanding for a guide node; merely opening a card is insufficient where
  a formative check is specified.
- `claim_ready`: every required scientific prerequisite for a particular
  ClaimID is satisfied for the current StudyContextID, AnalysisID, and
  InterpretationID.

These flags are independent. A learner may complete a node while the analysis
remains `HOLD`; a failed or unassessable run may be archived as a reproducible
no-conclusion record; and neither guide completion nor an acknowledgement may
make a substantive claim ready.

## Implementation record

### Help and audience-language shadow foundation — 2026-07-24

Delivered locally, without changing the visible Streamlit route:

- Added a pure three-layer terminology registry for task-facing, method, and
  technical vocabulary. Twenty initial review-tracked concepts retain exact
  reproducibility labels without exposing them through the Standard
  projection.
- Added immutable Help topic, section, link, target, and verified-context
  contracts plus a pure navigation reducer. Opening, changing, invalidating,
  localizing, and returning from Help carries presentation identity only and
  cannot mutate an analysis, evidence conclusion, or learning state. Route
  construction is reducer-only, including protection against forged
  `dataclasses.replace()` transitions; closing Help clears route history.
- Registered 22 task-centered Help topics, a reserved localized safe fallback,
  42 source/return targets, and 74 Help links, including exact mappings for 18
  existing popover keys. No invented topic alias is retained. All terminology
  and problem references resolve, every return focus is owned by its target,
  and the graph has no orphan nodes.
- Separated the fallback from Help-home entries and placed all 42 intended
  surface, focus, and panel identities in an independent provisional catalog.
  Consumer typos can no longer create their own valid targets; the remaining
  activation slices must bind and test each catalog entry against the real
  Streamlit surface.
- Added 11 privacy-safe user-problem records for input, mapping, estimation,
  resource, and residual-structure failures. A rendered notice can contain
  only a registered problem code, a bounded occurrence phase, and a freshly
  generated random support reference; exception text, paths, uploaded values,
  and analysis identifiers are not retained. Generic application failures use
  their own troubleshooting topic rather than an estimation-failure route.
- Added matched English and Japanese seed copy for 100 terminology keys, 34
  problem keys, 437 Help-topic keys, and 24 Help-navigation keys, including
  explicit boundaries for fit, residual structure, differential interaction,
  privacy, estimator choice, and reporting claims. Every reporting example is
  structurally paired with a warning to verify its conditions first.
- Added contract, graph, reducer, locale, exposure, privacy, and standalone-core
  tests. The seed Help copy contains no active external-runtime handoff and no
  default implementation vocabulary.
- Marked terminology and Help content `in_review`; SLA-method and UX review
  remain explicit gates before this copy replaces the visible application.
  The current ClaimBoundary/reference IDs and surface-target catalog are
  provisional contract identifiers, not completed bibliographic records or
  evidence that a live UI target exists. Those records and the remaining
  adapter bindings are separate activation gates across R2 and R3.

This was deliberately delivered first as a shadow slice. The bounded R2
activation below now consumes it for the global Help entry and first-failure
paths. The first R3 vertical slice now also activates the fit-scatter link;
legacy Help sections and every other chart popover remain compatibility
surfaces until their exact registered adapters pass the later gates.

### Persistent Help and safe-problem activation — 2026-07-24

Delivered locally:

- Added a persistent bilingual Help launcher to the application shell and a
  full main-content Help surface available before data are present. It uses
  stable topic/link IDs, renders only the Standard audience layer by default,
  checks lifecycle applicability, and excludes the system fallback from the
  22-topic picker.
- Kept the source data and settings controls alive in the sidebar while Help
  owns the main content. Paste drafts, delimiters, uploads, selected source,
  and keyed analysis settings therefore remain in Streamlit's normal widget
  lifecycle; the Help route never consumes a run or refit trigger.
- Added an in-Help topic selector, persistent English/Japanese switching, an
  explicit return control, localized return status, and a safe generic route
  for unknown, removed, or inapplicable destinations. Locale changes retain
  route identity rather than translating IDs.
- Activated privacy-safe problem notices on paste/upload, wide-to-long,
  simulation-threshold, and covariate-preview parsing, plus primary estimation
  failures. The visible notice contains registered bilingual copy, a random
  support reference, and exact registered Help actions; raw exception text,
  paths, input values, and identifiers remain absent from the public UI even
  if a legacy technical-error environment variable is present.
- Replaced visible “technical details” failure gates on these paths with
  expanded, task-facing next-action panels connected to Help.
- Added AppTests for no-data open → locale change → return, source-state
  preservation, zero estimator/refit calls, unknown-link recovery, denied and
  allowed lifecycle states, pending-trigger preservation, and sentinel
  exception redaction through both direct and simulation-threshold
  problem-to-Help paths. The zero-call guard instruments both the cached core
  namespace and direct module-level estimation/refit functions.

This is a bounded activation, not a claim that all 42 provisional targets are
live. The global start target and first problem routes are connected; legacy
popover, guided, result, report, and download adapters still require exact
surface/focus verification. Topic copy, ClaimBoundary IDs, and reference IDs
remain `in_review` until SLA-method and UX review are recorded.

### First contextual Help return slice — 2026-07-24

Delivered locally as a deliberately narrow R3 slice:

- Activated only `fit_scatter` among the 18 registered core popover aliases.
  Its button resolves the exact registered
  `link.popover.fit_scatter` -> `help.results.fit#fit-scatter` edge and the
  exact `target.results.figure.fit_scatter` /
  `focus.results.figure.fit_scatter` return pair. The other 17 aliases and the
  two legacy posterior popovers have no detailed-Help action in this slice.
- Required a current fitted context built from the retained result's
  AnalysisID and its stored sample/real provenance. A legacy output with no
  stored provenance is not inferred from the current sidebar selection, and
  missing, stale, mismatched, or malformed context fails closed.
- Added separate allowlists for presentation projection and target-owned
  scientific view state. Return restores the Essential or Full fit-details
  panel plus only `fit_df_method_method` and `fit_df_method_cap` through the
  return intent. Native MML prior-SD sensitivity draft controls and misfit
  ranking display controls use separately validated, AnalysisID-scoped durable
  session state, so Streamlit widget cleanup cannot reset them while Help owns
  the main surface or carry them into a different fitted analysis. Analysis
  settings, data, evidence, and guide state remain untouched.
- Removed the legacy imported-FACETS fit-table comparator from the default Fit
  Details journey. Its compatibility function remains isolated and is not an
  active Help dependency or standalone-core workflow; no external comparison
  upload needs to be retained by this return adapter.
- Bound the return intent to the exact link, target, focus, projection,
  AnalysisID, provenance, schema version, and validated source-view values.
  Missing or altered intent opens the reserved localized fallback instead of
  guessing a destination or reporting a successful return.
- Rechecks the current fitted phase, AnalysisID, and sample/real provenance
  before restoring any panel or source value. Browsing a static Help topic can
  still return through an intact origin intent, while an A-to-B analysis or
  provenance change fails before any A-owned value is written into B.
- Suppressed analysis triggers for the one return rerun without consuming
  pending user-initiated run/quick-start flags. The real-app AppTest verifies
  zero estimator/refit calls, the same AnalysisID, bilingual route retention,
  recovery of a `both` / `12.5` fit-display convention after Streamlit removes
  the hidden widgets, retention of native MML-sensitivity and misfit-ranking
  drafts, and consumption of the intent at the real figure marker.
- Made a new Help journey an explicit cleanup boundary for an abandoned
  `awaiting_render` intent, its one-shot trigger suppression, and an old return
  notice. A target that was not rendered therefore cannot poison a later
  global Help journey.
- Anchored the real localized Fit scatter `H3` with the registry-owned marker
  ID and, only after a valid intent reaches that rendered heading, emits one
  static focus/scroll request. The progressive-enhancement script contains no
  AnalysisID, provenance, locale copy, or input data; validates the exact
  connected `H3`; waits two animation frames; avoids stealing focus after a
  user move; and guards duplicate execution. The visible localized return
  confirmation is a non-live caption, so an empty ARIA marker and duplicate
  success-alert announcement are not introduced.
- Encoded adapter maturity as `APPTEST_ONLY` or `BROWSER_ACCEPTED`. Exactly one
  AppTest-only pilot is allowed; every other active adapter requires a stable
  retained browser-evidence reference. `fit_scatter` remains the sole pilot,
  has no browser evidence ID, and is not browser accepted. The executable
  matrix and promotion procedure are in
  [`docs/help_browser_acceptance.md`](docs/help_browser_acceptance.md).

This is not R3 completion. AppTest verifies the native heading anchor, exact
static controller body, JavaScript opt-in flag, one-shot emission, localized
caption, and fail-closed absence for invalid or stale intents; it does not run
the browser script. The exact `document.activeElement`, viewport placement,
visible focus outline, subsequent Tab order, CSP behavior, zoom, reduced-motion
behavior, and single screen-reader announcement still require the browser
acceptance pass. A remaining popover can become the next sole pilot only after
the current pilot passes those gates and its reviewed evidence reference is
recorded; defining the checks is not evidence that they passed.

### Evidence-contract foundation — 2026-07-24

Delivered:

- Added `mfrm_app/evidence.py` with versioned contracts for deterministic
  analysis identity and lineage, evidence records, sensitivity comparisons,
  sensitivity decisions, workflow dispositions, and stable reason codes.
- Fixed the distinction between computation state and conclusion stability in
  executable validation rules rather than relying on display strings.
- Added order-independent sensitivity synthesis: a conclusion is `STABLE` only
  when every required credible variant succeeds, is comparable, and leaves the
  conclusion unchanged.
- Bound each prespecified sensitivity label to an embedded, self-validating
  baseline/variant AnalysisIdentity ledger, evaluate a closed versioned rule
  against archived metrics with fail-closed three-valued logic, and
  content-fingerprint the complete aggregate outcome.
- Added strict bundle restoration: sensitivity and workflow decisions are
  accepted from JSON only when their linked plan, evidence, and source records
  reproduce the saved decision exactly; unknown fields, type-coerced arrays,
  registry-version mismatches, and unresolved EvidenceIDs are rejected.
- Added a static core-boundary test that prevents `mfrm_app/` from importing
  Streamlit, invoking subprocess/dynamic external engines, or embedding new
  TAM/sirt/mirt/mfrmr/CmdStan handoffs.
- Removed the legacy parity-fixture exporter from the CLI, `make verify`, and
  GitHub CI, and removed cross-package, Stan, and Posterior compatibility tests
  from the standalone `--self-test` registry.
- Updated contribution and release guidance so new statistical computation is
  modular, Python-native, reason-coded, and validated with native fixtures.
- Added `mfrm_app/**` and this roadmap to CI path filters so changes to the
  modular core cannot bypass the Python test workflow.
- Removed the repository-external `Simulation/engines/mfrm_python.py` loader,
  its import-time execution path, and the unreachable shared-engine dispatch.
  An entrypoint AST guard now prevents that runtime loader from returning.
- Added `mfrm_app/readiness.py` and routed final-report readiness through it.
  The original five display columns remain compatible, while every row now
  carries a deterministic AnalysisID/EvidenceID, explicit computation and
  stability states, stable ReasonCode, scope, source-artifact links, and a
  validated full EvidenceRecord export.
- Separated the different meanings previously hidden behind final-readiness
  `Review`: non-convergence is `HOLD`, an executed diagnostic requiring review
  is `CAUTION`, and skipped PCA, bias, or strict-marginal evidence is
  `NOT_ASSESSABLE`.
- Added `mfrm_app/mml_prior_sensitivity.py` and connected fixed-prior-SD MML
  refits to exact baseline/variant identities, a prespecified executable plan,
  sensitivity records, and a fail-closed SensitivityDecision. The executed
  plan and JSON contract are now included in the Python-native sensitivity
  archive, and final readiness consumes the same decision.
- Closed prior-SD edge cases: free-population-SD fits and baseline-only grids
  are rejected as inapplicable/unassessable, `max_values` fixes the plan before
  refitting, and run failure, non-convergence, missing rule metrics, and
  non-comparable measures cannot produce a false stable result.
- Hardened the sensitivity comparison against false stability: outcome flags
  require exact booleans; measure and population-term keys must be complete,
  finite, and duplicate-free; rank evidence must be evaluable for every
  eligible facet; and each returned refit must attest the planned fixed-SD,
  engine, quadrature, tolerance, iteration, no-plausible-value controls, and
  unchanged population/anchor/constraint semantics.
- Revalidate saved MML bundles by reconstructing their identities, evidence,
  plan, records, and synthesized decision before readiness or export consumes
  them; a mutated same-AnalysisID payload is ignored rather than displayed.
- Unified the MML sensitivity baseline with the final-readiness AnalysisID and
  reject stored sensitivity evidence whose baseline identity differs from the
  current fitted result.
- Added reconstructable AnalysisIdentity, EvidenceRecord, sensitivity-contract,
  settings, and result-table sidecars to the standard Python-native ZIPs. Their
  source-artifact links are tested against the actual archive contents rather
  than surviving only as in-memory DataFrame attributes.
- Reduced `--release-check` to the standalone Python contract: its JSON has
  only release metadata, native readiness, and native limitation rows, and it
  no longer generates or reports external inventories, Stan handoffs, or
  `mfrmr` migration coverage.
- Removed external simulation inventories/templates, R/Julia reproducibility
  bundles, TAM/sirt/mirt handoffs, Stan/Posterior packages, and `mfrmr`
  migration tables from the default Downloads and deterministic demo archives.
  Python engine/local-batch/self-contained-JMLE assets and the Evidence
  sidecars remain first-class native outputs.
- Removed Posterior Viewer and advanced Stan generation from the sidebar,
  Report, Help, and tutorial routes. Public Yardstick and rating-scale recode
  downloads now enter Python-only generator branches directly instead of
  constructing legacy R/Julia assets and filtering them afterward.
- Split dormant compatibility tests behind the `legacy_compat` marker. The
  standard `make apptest`, `make verify`, and GitHub workflow now test only the
  standalone Python release surface; frozen compatibility fixtures are not
  release evidence.
- Added demo-specific handoff filenames so the generated checklist and
  manuscript handoff point to the files that the deterministic demo actually
  writes.
- Applied the demo artifact profile to every generated table, Markdown/HTML
  report, visual evidence map, and binder README, and added the previously
  referenced `export_privacy_manifest.csv` as an actual synthetic-demo
  artifact. Recursive archive verification now finds no stale app-download
  names in the demo packages.

Remaining before the product-boundary completion gate:

- Record and execute the explicit retain/deprecate/remove decision for dormant
  generator definitions, unreachable locale copy, and frozen compatibility
  fixtures. Public Help, Posterior, Report, Downloads, demo, CLI, self-test, and
  CI routes are already isolated; none of the dormant surfaces may re-enter
  them.

## Delivery sequence

The order below is dependency-driven. Later outcomes should not be built on
unversioned tables or UI-only logic.

### Product boundary and modular contracts

Objective: make the standalone-Python boundary enforceable and give new work a
stable home outside the large Streamlit entrypoint.

Deliverables:

- Record this roadmap as the active product plan and retain completed release
  plans as historical documents.
- Define versioned Python contracts for analysis identity, evidence records,
  sensitivity results, decision states, and skip/hold reasons.
- Move computational network and decision helpers into focused `mfrm_app/`
  modules while retaining compatibility wrappers where tests or public imports
  require them.
- Separate UI rendering, computation, and export-frame construction.
- Update contributor guidance so new statistical work is modular by default.
- Add a boundary test that rejects runtime imports, subprocess calls, and UI
  actions that invoke external estimation engines in the core workflow.
- Classify existing legacy external-comparison surfaces as compatibility-only
  so they cannot be selected as prerequisites by new features, appear in the
  default workflow, or gate a release.
- Produce a separate compatibility decision for legacy code generators,
  posterior-result viewers, cross-engine fixtures, and handoff downloads. Do
  not remove them incidentally during statistical refactoring.

Completion gate:

- Existing public fixtures and result columns remain compatible or have an
  explicit migration note.
- All tests pass after each extraction, with no statistical change hidden in a
  modularization commit.
- Every unavailable computation returns a stable reason code rather than an
  empty result or silent exception.
- Core estimation, diagnostics, sensitivity, and reporting execute without an
  external engine installed.

### Resilient rating design

Objective: progress from a static connectivity summary to an analysis of how a
design fails and how it can be strengthened.

Deliverables:

- Maintain distinct graph views for the global multipartite design, direct
  rater-person links, rater-task/criterion coverage, and projected rater
  overlap. Do not compare density across unlike graph definitions.
- Add weighted edge strength, minimum shared cases, bridges, articulation
  points, component structure, workload imbalance, and coverage summaries.
- Add deterministic perturbation audits for:
  - leave-one-rater-out;
  - leave-one-task/item/criterion-level-out where applicable;
  - removal of singleton or weakest links;
  - prespecified random rating loss; and
  - concentrated missingness in one facet region.
- Report the proportion and identity of perturbations that disconnect the
  design or materially reduce overlap.
- Generate constrained repair candidates, such as additional shared ratings,
  while displaying burden and the exact diagnostic each candidate repairs.
- Keep connectivity and precision distinct: a repaired connection is not
  automatically an adequately precise design.

Completion gate:

- Complete, chain, star, bridge-dependent, disconnected, and imbalanced
  synthetic fixtures return their known graph properties.
- Leave-one-level-out results are deterministic and use the same validated
  observation set as the fitted analysis.
- Every proposed repair reconnects the targeted synthetic design when a repair
  is feasible, and no repair is labelled as sufficient evidence of precision.
- Exported graph definitions, thresholds, removed observations, and
  perturbation seeds reproduce the displayed result.

### Coherent diagnostic evidence

Objective: combine diagnostics without treating them as a vote or collapsing
them into an unsupported universal score.

Deliverables:

- Normalize fit, covariance, residual PCA, PCA stability, DIMTEST,
  local-dependence, category, bias/interaction, reliability, design-network,
  severity-network, and halo-network outputs into evidence records.
- Record for each item of evidence:
  - analysis identity and scope;
  - data and configuration fingerprints;
  - applicability prerequisites;
  - computation state and reason code;
  - observed statistic and uncertainty when available;
  - interpretation boundary;
  - recommended next inspection; and
  - source table or figure.
- Build an evidence-coherence view that highlights agreement, tension, and
  unresolved prerequisites. It must not use majority voting.
- Encode explicit interpretation rules for common tensions, including:
  - a residual PCA signal with unstable bootstrap loadings;
  - a DIMTEST result that conflicts with residual PCA;
  - acceptable global fit with concentrated local dependence;
  - high reliability with fragile network connectivity; and
  - a rater severity signal that disappears after residual adjustment.
- Use the same evidence records in guided UI, research UI, readiness checks,
  narrative drafts, and exports.

Completion gate:

- Contradictory synthetic evidence produces a documented tension state, not a
  false all-clear or automatic dimensionality claim.
- No final-readiness statement cites an unavailable, stale, or differently
  fingerprinted diagnostic.
- Every user-facing recommendation links to its evidence and caveat.
- Locale and export tests confirm that status meanings are consistent in
  Japanese and English.

### Conclusion robustness

Objective: evaluate the stability of the decisions users actually make, not
only the stability of individual coefficients.

Deliverables:

- Define a versioned sensitivity specification before running a grid.
- Support applicable checks across:
  - JMLE/MML and supported MML engine choices;
  - fixed population prior SD settings;
  - RSM/PCM/GPCM only when the data and comparison question permit it;
  - category handling and declared score support;
  - influential facet levels;
  - reasonable missingness perturbations; and
  - documented fit-statistic conventions.
- Track conclusion-level outcomes such as:
  - facet severity direction and rank bands;
  - misfit and bias-review flags;
  - reliability interpretation;
  - category-functioning warnings;
  - dimensionality-screen wording; and
  - final claim availability.
- Distinguish a numerical change from a substantive conclusion change.
- Preserve failed specifications in the denominator and explain their failure;
  do not summarize only successful refits.
- Produce a compact stability map plus a full specification ledger.

Completion gate:

- Fixed fixtures demonstrate each stability state and prevent stale-base-run
  reuse.
- The specification ledger reproduces every displayed comparison from frozen
  inputs and seeds.
- A failed or unsupported refit cannot strengthen the final claim.
- Rank, flag, and narrative changes use documented thresholds rather than
  post-hoc visual judgment.

### Calibrated simulation and stress testing

Objective: show under which data-generating and design conditions the app's
diagnostics and decisions behave as intended.

Deliverables:

- Complete the modular condition/scenario system for immutable, serializable,
  fingerprinted simulation specifications.
- Expand native DGP coverage for supported RSM, PCM, and bounded GPCM use cases
  without implying support for an unfitted model family.
- Add prespecified stressors for:
  - disconnected and bridge-dependent assignments;
  - unequal rater workload and restricted coverage;
  - severity spread and extreme raters;
  - halo and supported local interactions;
  - rater drift across occasions;
  - missing completely at random and concentrated design loss;
  - sparse and unused score categories; and
  - non-normal person distributions where the estimator permits evaluation.
- Evaluate bias, RMSE, interval coverage, rank recovery, flag rates,
  disconnection rates, conclusion reversals, numerical failure, and Monte
  Carlo uncertainty.
- Maintain an ADEMP-style scenario ledger that separates development checks,
  calibration runs, and frozen validation simulation runs.
- Add checkpoint/resume only as a Python-native computation feature, not as an
  external-engine exchange contract.

Completion gate:

- Repeating a scenario with the same specification reproduces task identities,
  seeds, truth, and summaries.
- Scenario additions do not change the random streams of existing frozen
  scenarios.
- Rates retain all requested replications in their denominator and report
  numerical failures separately.
- Proportions include Monte Carlo uncertainty; small development runs are not
  presented as performance evidence.
- Each public diagnostic claim is backed by a named set of frozen simulation
  conditions or is explicitly labelled uncalibrated.

### Guided decisions and research transparency

Objective: serve first-time and expert users without maintaining two different
statistical implementations.

This work is an information-architecture consolidation, not another tutorial
surface layered onto the existing onboarding, Tutorial, Guided, and Help
systems. The application already contains substantial explanatory material.
The long-term problem is to make one route and one evidence state authoritative
while rendering different levels of detail for different users.

#### P0 research, interpretation, and learning contracts

The interface will not be redesigned until the following contracts are frozen
in testable, widget-independent form. The detailed implementation contract is
maintained in
[`docs/study_context_guidance_contract.md`](docs/study_context_guidance_contract.md).

| Contract | Minimum responsibility | Boundary |
|---|---|---|
| `StudyContext` | Purpose, stakes, target population/domain, construct and subconstructs, performance mode, intended score interpretation and use, intended decision, and unsupported inferences | Does not assert validity, fairness, or fitness for use |
| `RatingDesignRecord` | Object of measurement, nuisance facets, crossed/nested/anchored/fixed/random roles, planned and observed assignment, common performances, occasions/order, double scoring, adjudication, rubric version, and missingness plan | Graph connectedness alone does not establish comparability or generalizability |
| `ScoringProcessRecord` | Rater selection, training, qualification, monitoring, rescoring, adjudication, and benchmark process | Statistical flags do not authorize punitive personnel action |
| `AnalysisSpec` | Canonical data mappings, model, estimator, identification, diagnostic, sensitivity, and seed settings independent of Streamlit widget visibility | View, guide, language, and theme are excluded |
| `GuidanceState` | Stable route/node IDs, progress, sample context, bindings, invalidation, and lifecycle events | Cannot mutate computation or conclusion state |
| `HelpTopic` / `HelpLink` | Stable concepts, audience vocabulary, applicability, claim boundaries, references, exact contextual targets, and return/focus behavior | Help navigation cannot mutate science/learning or present stale analysis context as current |
| `EvidenceIssue` | One prioritized question with all relevant EvidenceIDs, tensions, applicability, boundaries, and next actions | Never hides other `HOLD`, `CAUTION`, or required-unassessed evidence |
| Interpretation and decision records | Stable claim, limitation, action, authorship, and review-state identities linked to study and analysis context | An analysis record may document that no interpretable conclusion is available |

Planned absence and unexpected missingness are represented separately. At a
minimum the application distinguishes planned incomplete assignment,
test-taker nonresponse, not-reached responses, technical failure, rater skip
or unratable response, criterion-level absence, administrative exclusion, and
adjudication-only ratings. An unexplained `NA` is not silently treated as a
planned absent design cell.

For analytic rubrics, the study context must state whether criteria are
intended as one construct, distinct subconstructs, or a reported composite.
Treating criterion as a facet does not by itself justify an overall person
measure. If the intended claim requires a multidimensional model outside the
native scope, the claim is `HOLD` or `NOT_ASSESSABLE`, not rescued by a clean
residual screen.

All currently visible explanatory text is part of the evidence surface. Help,
tutorial, chart guides, popovers, warnings, captions, and exported prose are
therefore audited in G0 against
[`docs/interpretation_copy_audit.md`](docs/interpretation_copy_audit.md).
Later Evidence-linked copy cannot coexist with stronger legacy statements
elsewhere in the application.

The file-level remediation sequence, registry design, concrete bilingual copy,
Help lifecycle behavior, migration slices, and test changes are specified in
[`docs/help_and_terminology_remediation_plan.md`](docs/help_and_terminology_remediation_plan.md).

#### Non-negotiable UX principles

- The guide is optional, skippable from every step, restartable, and resumable
  within the session.
- Skipping the guide never deletes uploaded data, fitted results, the current
  AnalysisID, EvidenceRecords, or a draft decision note.
- Acknowledging a guide card never changes `HOLD`, `NOT_ASSESSABLE`,
  `CAUTION`, or conclusion-stability states. Only new data, settings,
  computation, or evidence may change an analysis state.
- Guided and research surfaces use the same estimator results,
  AnalysisIdentity, EvidenceRecords, sensitivity records, and reason codes.
- Tutorial progress is not analysis progress, and guide completion is not a
  certificate of MFRM competence or evidence that a real-data analysis is
  valid.
- The guide teaches the standalone Python workflow. TAM, R, Julia, external
  posterior ingestion, and cross-engine ZIP compatibility do not appear in the
  guide or become prerequisites for it.
- The sample route is allowed to teach operation and interpretation, but a
  successful sample run must explicitly state that it does not validate the
  user's data or design.
- The guide remains a thin projection over the production computation path; it
  does not maintain canned results or a separate tutorial estimator.

#### Entry routes

The landing surface offers three clear choices without asking users to label
themselves by an assumed ability level.

| Entry route | Behaviour | Guidance level |
|---|---|---|
| Learn with a sample | Continue from the landing surface into the deterministic guide | Full first-run guide |
| Start with my data | Open the real-data workflow | Contextual readiness and evidence guidance |
| Continue without the guide | Open the normal Standard view | No tutorial sequence |

The landing surface itself is the `welcome` node. It records the selected
route and moves forward in one action; users are not asked to make the same
sample-versus-data choice on a second welcome screen. A future pre-collection
route may begin from `Plan`, but it must reuse the same StudyContext and
RatingDesignRecord rather than become a fourth onboarding system.

The fixed tutorial uses a built-in sample because an arbitrary uploaded data
set may be disconnected, sparse, malformed, or non-convergent. Real data
therefore receives evidence-dependent coaching rather than being forced
through a linear tutorial. After the basic sample, an optional single-warning
practice scenario may teach one `CAUTION` or `NOT_ASSESSABLE` branch at a time.
It must not combine many artificial problems into one confusing exercise.

#### Five-step guided first-run route

Stable internal node IDs, rather than translated labels, drive routing.

| Node ID | User-facing purpose | Learning/workflow completion | Scientific effect | Primary output |
|---|---|---|---|---|
| `welcome` | Choose the sample or leave the guide | A landing route is selected | None | Bound route |
| `data_check` | Confirm roles, mappings, design, and readiness | The required concepts and current preflight outcome were reviewed | `HOLD` may remain and directs repair | Data fingerprint, design record, and readiness state |
| `estimate` | Run the bound specification and inspect its numerical outcome | The current run outcome, including failure, was reviewed | Only a successful current fit creates usable fitted evidence | Analysis identity or a versioned failure record |
| `evidence_review` | Examine the priority evidence issue and its claim boundary | The learner distinguishes one supported from one unsupported statement in the sample exercise | Required evidence may remain `HOLD`, `CAUTION`, or `NOT_ASSESSABLE` | Evidence-linked boundary and next action |
| `archive` | Save a bounded interpretation/decision draft | A supported claim or explicit no-conclusion statement, a limitation, and a next action are recorded | Archiving does not certify or improve claim readiness | Human-readable Analysis Brief plus JSON sidecar |

The visual sequence is:

```text
Landing/welcome -> Data check -> Estimate -> Read evidence -> Save record
                                      |
                                      +-- HOLD: repair, or archive why no
                                          interpretable conclusion is available
```

The guide is a five-node learning projection over the full research workflow,
not a replacement vocabulary. The mapping is shown throughout the guide so
the final transition does not create a second onboarding cliff.

| Guide node | Long-term research stages |
|---|---|
| `welcome` | Plan |
| `data_check` | Plan -> Check |
| `estimate` | Estimate |
| `evidence_review` | Diagnose -> Stress -> interpretation substep |
| `archive` | Decide -> Archive |

The canonical seven-stage workflow remains stable. `Interpret` is an explicit
substep between `Stress` and the use decision inside the `Decide` stage; its
InterpretationRecord and the subsequent DecisionRecord remain distinct
identities even when the navigation groups them on one screen.

Each guide screen has at most one primary card, three headline metrics, and one
state-changing primary call to action. `Exit guide` is persistent escape
chrome, not the screen's secondary action budget. Non-mutating Help and detail
disclosures may remain available without competing visually with the primary
action. Long tables, model theory, all-panel tours, and publication options
remain in contextual detail or Help.

#### Skip, exit, resume, and restart contract

The guide state has four lifecycle values: `NOT_STARTED`, `ACTIVE`, `SKIPPED`,
and `COMPLETED`. It also records the active node, last valid node, completed
nodes, entry route, bound data fingerprint, bound AnalysisID, reviewed
EvidenceID, and the node at which the guide was skipped.

- Before starting, the low-emphasis action is `Continue without the guide`.
- During a guide, the action is `Exit guide`, because users are skipping an
  explanatory sequence rather than statistical prerequisites.
- Exiting before estimation opens the current data/check surface.
- Exiting after estimation opens Current focus without refitting.
- Exiting while reading evidence preserves the same EvidenceRecord and detail
  target in the normal workspace.
- Resuming uses the first valid incomplete node, not blindly the previously
  displayed node.
- Restarting resets guide progress only; it does not clear input data, results,
  AnalysisIdentity, EvidenceRecords, or decision drafts.
- A skipped guide does not repeatedly reopen or display a blocking prompt in
  the same session. A quiet resume action remains in Help or the sidebar.
- Starting or resuming the sample guide while a real-data analysis exists uses
  an isolated tutorial context. It snapshots and restores the real-data
  `AnalysisSpec`, fitted identity, evidence, and draft notes; the sample may
  never overwrite them or silently become the current research analysis.
- Guide state is session-only in the initial implementation. Persistent browser
  preference is deferred until observed user friction justifies the privacy,
  accessibility, and maintenance cost.

No confirmation dialog is required for exiting because the operation is
non-destructive. The UI states that data and current results are retained.

#### Standard and Detailed views

Completing or skipping the tutorial leads to the Standard view, not directly
to every expert panel. This prevents an onboarding cliff in which a simple
guide ends at the existing high-density interface. The user-facing names avoid
using `Guided` for both the guide and a workspace, and avoid implying that only
a `Research` workspace is suitable for research.

| Stable view ID | User-facing projection | Purpose |
|---|---|---|
| `standard` | Standard view | One current evidence issue, its boundary, and next action |
| `detailed` | Detailed view | Resolved settings, full diagnostics, sensitivity ledgers, simulations, and exports |

The view selector is a presentation projection only. It never applies
estimation defaults, writes hidden widget values, changes a draft or fitted
`AnalysisSpec`, changes AnalysisID, invalidates a fit, triggers estimation, or
starts a hidden diagnostic. `AnalysisSpec` is canonical state outside
Streamlit widgets so temporarily undisplayed advanced controls retain their
values. The interface distinguishes fitted settings from an edited draft and
offers an explicit return to the fitted specification. Analysis depth remains
an Estimate setting because it changes which diagnostics and artifacts are
computed.

The canonical seven-stage navigation remains:

```text
Plan -> Check -> Estimate -> Diagnose -> Stress -> Decide -> Archive
```

Existing surfaces migrate without duplicating computation:

| Existing surface | Long-term home |
|---|---|
| Start | Plan / Check |
| First Read | Diagnose / Decide |
| Results | Estimate / Diagnose |
| Diagnostics | Diagnose / Stress |
| Figures | The relevant result or diagnostic panel |
| Report & Export | Decide / Archive |
| Learn | Persistent Help, not an analysis stage |

The Standard view starts with one Current focus card backed by an
`EvidenceIssue`, not a single convenient result. It prioritizes stale results,
non-convergence, disconnected designs, required-but-unassessable evidence,
material diagnostic cautions, evidence tensions, unassessed sensitivity,
reporting limits, and then the standard next result. It also shows the total
number of `HOLD`, `CAUTION`, and required-unassessed issues and why the current
issue was selected. Additional items remain available without competing with
the highest-priority action.

#### Information density and terminology

During the sample guide, the tutorial binds a documented sample
`AnalysisSpec`; it does not simulate a choice among display, control, and
computation presets. Standard-view sidebar content is limited to data, column
mapping, a resolved estimation summary, and the run action. Optimizer,
tolerance, anchors, visual preferences, custom diagnostic options, simulation,
and publication controls stay collapsed or appear in Detailed view. Hiding a
control never resets its canonical value.

Guide labels use plain language followed by the canonical technical term
where useful, for example:

- Model discrepancy (`Fit`);
- Residual-structure check (`Residual PCA`);
- Element differentiation (`Reliability / separation`); and
- Differential interaction screen (`facet interaction`; legacy `bias` labels
  are retained only where required for a clearly marked method reference).

Rater flags are framed as review prompts for the design, observation count,
rating examples, rubric, and scoring process. They do not label a rater as
accurate, fair, competent, or unsuitable, and they do not automatically
recommend exclusion, retraining, or personnel action.

The default run action uses standalone product language such as `Run MFRM
estimation`, not a label that implies runtime dependence on FACETS or another
external product.

#### Audience language and Help continuity

Implementation vocabulary is not automatically user vocabulary. Internal
enums, record types, IDs, cache/schema terms, and support codes remain in the
core contract and machine-readable sidecar; they are not the default wording
of the guide, Standard view, or Help landing page.

| Presentation layer | Intended language | Examples |
|---|---|---|
| Guide and Standard view | Plain task, status, boundary, and action language | “Study purpose,” “analysis settings,” “why this result is limited,” “stop and repair,” “analysis reference” |
| Detailed view and method Help | Canonical MFRM/SLA terminology with a short definition and applicability boundary | Residual PCA, Infit/Outfit, separation, differential interaction, JMLE/MML |
| Technical details, JSON, and privacy-safe support | Exact implementation names and stable identifiers | `StudyContextID`, `AnalysisID`, `EvidenceID`, `ReasonCode`, schema version |

For example, the default interface says “current evidence question,” not
`EvidenceIssue`; “interpretation draft,” not `InterpretationRecord`; and “why
this is limited,” not `ReasonCode`. A stable analysis reference may be shown
and copied in Detailed/technical details, but the user is not instructed to
reason in terms of internal class or enum names. Status codes such as `HOLD`
are paired with a localized action label such as “Stop and repair”; the code is
secondary technical detail.

Help is one connected content system, not a second tutorial or a collection of
unrelated tables. A versioned `HelpTopic` registry owns stable topic and
concept IDs, locale keys, audience layer, prerequisites, applicability,
ClaimBoundaryIDs, references, related workflow targets, and review metadata.
Every guide node, EvidenceIssue, chart popover, warning, and advanced control
links through a `HelpLink` to that registry.

Persistent static Help is reachable before data entry, during mapping and
preflight, after a failed first run, in every `HOLD`/`NOT_ASSESSABLE` state, and
after a successful fit. It never depends on a fitted result merely to render.
Only the optional “For this analysis” block requires a verified current
analysis binding.

The Help path has three depths:

```text
Inline boundary -> focused Help topic -> method/reference detail
       ^                    |
       +------ return ------+
```

- Inline Help answers what is shown, how to read it, and the most important
  limit without leaving the current task.
- `Open detailed help` is an actionable control, not a caption telling the
  user to find another page manually. It opens the exact topic/section,
  preserves the originating stage/card and current analysis context, announces
  the new heading after rerun, and provides a return action.
- Opening Help, following a HelpLink, searching, or returning changes no
  AnalysisID, setting, evidence, claim readiness, or learning-completion state.
- A missing/unknown HelpTopic fails a contract test and renders a safe generic
  fallback with a support reference; it never silently removes the Help
  control.
- Topic/target routing uses stable IDs only, never translated labels, screen
  names, or Japanese/English substring matching. The glossary has one
  canonical source for inline definitions, Help search, and exports.
- A claim-critical limitation remains inline. Help may explain it but cannot
  be the only place where it appears.
- Standard and Detailed views use the same topic and ClaimBoundaryID. Detailed
  view may reveal method/reference depth, but it does not receive a different
  substantive conclusion.

The Help home supports search and browsing by users' questions rather than
internal object names: research question; design and missingness; model and
estimation; diagnostics; interpretation and fairness; reporting and review;
reproducibility, privacy, and ethics; troubleshooting; and glossary and
references. Each focused topic states the question, prerequisites, what is
computed, what it can and cannot show, the next action, safe/unsafe reporting
language, references, related screens, and return target. Dense developer or
reviewer inventories and machine filenames move to optional technical details
or exported support/reviewer artifacts, not the guided Help path.

#### Evidence-linked teaching and claim boundaries

The evidence-review card is rendered directly from an EvidenceRecord. It shows
the computation state, what the evidence supports, what it does not establish,
the recommended action, the next inspection, and the source artifact. No
parallel tutorial copy is allowed to silently contradict the evidence ledger.

All guide-facing statements are bounded by the current data, model, and
computed scope. In particular:

- convergence means that the numerical stopping rule was met; it does not
  establish fit, validity, fairness, or conclusion stability;
- fit statistics are screening evidence, not universal pass/fail rules or
  labels for a person or rater;
- residual PCA can identify residual-structure signals but cannot prove
  unidimensionality when no large component is detected;
- orderly categories do not prove scale validity, and a disordered threshold
  does not mechanically require category collapse;
- rater reliability can indicate distinguishable severity differences; it is
  not evidence that raters agree or are interchangeable;
- dropping missing rows does not establish an ignorable missingness mechanism;
  and
- no flagged interaction within the computed pairs and detection conditions is
  not proof that bias or unfairness is absent.

Normal-state language therefore uses formulations such as `no high-priority
flag was detected under the configured screen`, not `good`, `normal`, or `no
problem`. Each card retains the applicable method and limitation note.

Before real-data upload, the application presents deployment-specific data
governance guidance covering processing location, retention/cache limits,
de-identification, consent or secondary-use authority, access control,
small-cell and intersectional re-identification risk, and data that must not be
uploaded to a public deployment. An acknowledgement records that guidance was
shown; it does not authorize processing or unlock a high-stakes claim.

#### Guidance architecture

Add `mfrm_app/guidance.py` as a pure-Python module that imports neither
Streamlit nor pandas. It owns stable node, goal, stage, completion-rule, and
target IDs plus a deterministic state reducer. It stores i18n keys, not display
strings. Streamlit adapters translate only at render time.

The authoritative model separates:

- learning state: `LOCKED`, `AVAILABLE`, `ACTIVE`, `COMPLETE`;
- workflow completion: whether the bounded route reached a terminal record;
- computation state: `AVAILABLE`, `CAUTION`, `HOLD`, `NOT_ASSESSABLE`;
- conclusion stability: `STABLE`, `CONDITIONALLY_STABLE`, `SENSITIVE`,
  `NOT_ASSESSED`; and
- claim readiness: a ClaimID-specific derivation from required evidence,
  context, applicability, review, and stability rather than a user toggle.

User acknowledgement may complete a learning node but cannot alter the latter
three scientific systems. Existing onboarding, first-run route, Guided, and
Help functions remain compatibility wrappers until their AppTests are
migrated. Routing by translated label or English substring is replaced by
explicit target IDs.

The first implementation consolidates route identity, prerequisites,
completion rules, targets, and locale keys. Long conceptual Help text remains
separate content projected from the same concept IDs; forcing every paragraph
into one catalog would create another hard-to-maintain monolith.

#### Invalidation and identity rules

| Change | Guidance effect | Analysis effect |
|---|---|---|
| Data or required mapping changes | Re-evaluate from `data_check` | New or stale AnalysisID |
| Model or estimator setting changes | Re-evaluate from `estimate` | New or stale AnalysisID |
| Analysis depth or computed diagnostic changes | Re-evaluate affected evidence | New evidence/analysis identity under the current contract |
| New fitted AnalysisID | Clear analysis-dependent review and archive completion | Bind to new result |
| Intended purpose, population, interpretation, use, or stakes change | Re-evaluate interpretation and decision nodes | Preserve AnalysisID unless the change also alters `AnalysisSpec`; create a new StudyContextID/InterpretationID/DecisionID |
| Planned or observed rating-design/scoring-process facts change | Re-evaluate `data_check` and affected claims | Preserve or invalidate AnalysisID according to whether data mappings/specification changed; always version the design/process record |
| View, guide, or language changes | Preserve guide validity | Preserve AnalysisID and the fitted/draft `AnalysisSpec` |
| Evidence card is opened | Record learning progress only | Preserve evidence state |
| Decision note changes | Update the decision artifact | Preserve AnalysisID |
| Pure display-theme change | Preserve guidance | Preserve statistical AnalysisID; update only artifact identity if needed |

A result whose data or settings no longer match the current context cannot
complete `estimate`, satisfy evidence review for the current context, or
produce a claim-ready archive. The guide returns to the first invalidated node
and explains what changed. It may still preserve an explicitly stale or
no-conclusion history record that cannot be mistaken for the current result.

#### Native analysis record

Add a versioned Python-native analysis-record contract after the study,
evidence, identity, and guidance contracts are stable. The primary user
artifact is a human-readable Markdown or HTML `MFRM_Analysis_Brief`; a
machine-readable `MFRM_Analysis_Record.json` is its deterministic sidecar.
Both derive from the same structured record and reference the same stable
StudyContextID, AnalysisID, EvidenceIDs, InterpretationID, DecisionID, and
ClaimIDs.

The structured record contains resolved settings, study/design/scoring-process
context, EvidenceRecords and unresolved evidence issues, supported and
unsupported claim entries, an explicit no-conclusion state when applicable,
the next action, authorship/review status, application/schema versions, and a
privacy declaration. Raw response rows and direct person/rater identifiers are
absent by default. Display language is artifact metadata; canonical ClaimIDs
and user-authored text are stored separately so translation does not silently
change a scientific identity.

The scientific record does not include tutorial progress, skipped-node
history, view preference, or external-engine compatibility metadata. Guide
state is not scientific evidence. Generating either artifact is distinct from
claiming that a browser download completed successfully, and the default
status is `draft_for_review`, not `validated` or `approved`.

#### Delivery slices

The UI migration is deliberately incremental. Accessibility and bilingual
meaning are acceptance criteria for every visible slice, not a final polish.

1. **G0 — Freeze contracts, behaviour, and interpretive risk.** Approve the
   StudyContext/Guidance/Identity contract, inventory every interpretation
   surface, classify unsafe, stale, or internal-only copy, map the current Help
   link graph and orphans, and add baseline tests for the built-in sample,
   AnalysisID, EvidenceRecords, stale detection, modes, locales, and exports
   without redesigning the UI.
2. **G1 — Canonical state and guidance core.** Add widget-independent
   `AnalysisSpec`, study/design/scoring-process records, the pure guidance
   catalog and reducer, the public-vocabulary map, `HelpTopic`/`HelpLink`
   registries, lifecycle and invalidation rules, locale-key parity tests, and
   compatibility wrappers. Replace G0 `STOP` copy before exposing a new
   evidence-linked route.
3. **G2 — Optional sample guide.** Replace the overlapping landing
   onboarding/Tutorial presentation with the five-step sample route, including
   skip, exit, resume, restart, formative interpretation, `HOLD`, stale, and
   non-destructive sample-context branches. Every contextual Help action opens
   the exact topic and can return to the originating node/card.
4. **G3 — View projection.** Map the existing Essential/Guided and
   Full/Advanced surfaces to Standard and Detailed views. Switching view must
   preserve fitted and draft specifications and invoke no estimator or hidden
   diagnostic.
5. **G4 — Evidence-to-decision artifacts.** Render the priority EvidenceIssue
   with its tensions and aggregate counts, collect bounded
   claim/limit/action or no-conclusion notes, and generate the human-readable
   Analysis Brief plus its privacy-safe JSON sidecar.
6. **G5 — Research and education companion.** Add an unknown-case transfer
   exercise, data-governance preflight, Methods/Results/Limitations reporting
   matrix, reviewer checklist, and versioned sample manifest without turning
   learning progress into scientific evidence or the application into an LMS.
7. **G6 — Consolidate and harden.** Move long theory to Help, remove only
   proven-unreachable duplicate guidance, complete the copy-audit remediation,
   bilingual AppTests, browser accessibility checks, and then retire
   compatibility wrappers.

Statistical changes, state-model extraction, visible UI changes, and legacy
code removal remain separate review units. Old routes are not removed until
new end-to-end gates pass.

The first G1/G2 vertical slice shipped on 2026-08-04: the pure guide catalog
and reducer, focused bilingual five-step sample journey, direct skip/exit,
resume/restart, production-estimator fit binding, formative overclaim guard,
and learning-only completion are implemented and AppTest-covered. G1 is not
complete: canonical `AnalysisSpec` and the remaining study/design/scoring
records still precede a full claim of canonical state. G2 is also not complete
until every node has exact contextual Help return, `HOLD`/stale/failure paths,
a complete existing-real-fit sample roundtrip, and browser keyboard/focus
acceptance.

#### Validation and evaluation

Unit tests cover unique IDs, valid routes and targets, completion/invalidation,
locale-key parity, one-time quick-start consumption, and the rule that user
acknowledgement cannot clear a machine `HOLD`. They also verify that view,
language, Help, learning answers, and guide lifecycle events cannot mutate the
fitted or draft `AnalysisSpec`, AnalysisID, EvidenceRecords, or claim readiness.

AppTest covers at least:

- Japanese and English sample completion;
- immediate skip by an experienced user;
- exit after estimation without refitting or losing AnalysisID;
- stale settings returning to the correct step;
- a `HOLD` path with a named repair and an archiveable no-conclusion record but
  no claim-ready continuation;
- a `CAUTION` path that carries its boundary into the decision note;
- language change at evidence review with the same AnalysisID and EvidenceID;
- resume and restart without deleting data or results;
- sample start/resume without overwriting an existing real-data context;
- Standard/Detailed projections displaying the same evidence state while
  retaining hidden advanced values and making zero estimator calls;
- independent `workflow_complete`, `learning_complete`, and `claim_ready`
  states;
- planned-absence and unexpected-missingness reason paths; and
- an evidence-tension fixture that cannot collapse into an all-clear message.
- persistent static Help access from landing/no-data, mapping/preflight,
  first-run failure, `HOLD`/`NOT_ASSESSABLE`, and post-fit states;
- one-action contextual Help routing to the exact topic and back to the same
  guide/EvidenceIssue/control without changing scientific or learning state;
- unknown, unavailable, and stale Help targets producing a visible safe
  fallback rather than a silent no-op; and
- absence of internal enum/class/file-pipeline language from the guide and
  Standard Help path, with optional technical identifiers still available in
  the declared technical layer.

Browser-level checks begin with the first visible slice and target WCAG 2.2 AA
in the supported Streamlit/browser scope. They cover keyboard-only completion
and exit, rerun focus/status announcement, logical headings and reading order,
320 CSS-pixel layouts and 200--400% zoom, long Japanese text, reduced motion,
touch-target size, non-colour state recognition, text/table alternatives for
figures, and keyboard-safe popovers/dialogs. AppTest alone is not treated as an
accessibility audit.

Initial evaluation uses deterministic tests and moderated observation rather
than production telemetry. Success is measured by whether a user can run the
sample, recover from a blocker, find Current focus, identify one supported and
one unsupported claim, and create the analysis record. A second, unfamiliar
single-problem SLA scenario measures transfer: the user must identify what is
supported, what is not supported, what to inspect next, and when supervisory
review is needed. Satisfaction or guide completion alone is not learning
evidence. If telemetry is ever
considered, it is opt-in, coarse, and excludes raw ratings, file/column names,
identifiers, fingerprints, AnalysisID, free text, and resolved settings.

#### Explicit non-goals for the initial guide

- A forced wizard or a separate top-level Tutorial tab.
- DOM-selector coach marks, automatic tutorial progression or selector-driven
  tour scrolling, hover-only instruction, or a video-first experience.
- A tour of every diagnostic, figure, export, and advanced setting.
- A second estimation or result pipeline for tutorial data.
- Universal threshold teaching, automatic validity claims, completion badges,
  or qualification certificates.
- Automatic high-stakes authorization, rater personnel decisions, or claims
  that a completed guide establishes research competence.
- A learning-management system, graded student accounts, or hidden instructor
  surveillance. Education worksheets and rubrics remain separable companion
  artifacts.
- Persistent browser tracking before a demonstrated need.
- TAM, R, Julia, external posterior, or cross-engine ZIP workflows.

Exact-source focus/scroll restoration after an explicit Help return is not
tutorial progression. It is permitted only for a registry-owned real heading,
must yield to user-moved focus, uses no smooth motion, and remains governed by
the contextual Help browser-acceptance runbook.

#### Completion gate

- Standard and Detailed views consume identical current EvidenceRecords and
  never display conflicting states for the same AnalysisID.
- The complete sample journey, immediate skip, exit at every step, resume,
  restart, stale invalidation, `HOLD`, and `CAUTION` paths pass in both locales.
- Skipping or switching view preserves data, fitted results,
  AnalysisIdentity, EvidenceRecords, and draft decision notes and does not
  trigger re-estimation.
- Every visible status has text, an interpretation boundary, a next action,
  and a method or limitation note; no critical meaning depends on colour.
- A user can identify one supported and one unsupported claim and can return
  from stale evidence to the setting that invalidated it without restarting
  the application.
- The Analysis Brief and JSON sidecar are schema-consistent, linked to the
  current study/analysis/interpretation/decision identities, privacy-safe by
  default, and contain no tutorial-progress or external-engine compatibility
  state.
- A failed, `HOLD`, or required-`NOT_ASSESSABLE` analysis can preserve a
  reviewable no-conclusion record without being presented as claim-ready.
- The default journey contains one primary action per guide screen, remains
  keyboard and narrow-screen usable, and does not eagerly render every heavy
  diagnostic panel.
- Every visible Help action reaches a registered bilingual topic in one action,
  retains a usable return target, and exposes claim-critical limits inline.
- Guide and Standard surfaces do not require users to understand internal
  class, enum, schema, manifest, runner, or raw ID names.

### Public-beta hardening

Objective: make release claims depend on native validation evidence and
operational reliability.

Deliverables:

- Replace external-parity language as a release gate with native numerical
  consistency, simulation calibration, and explicitly bounded method claims.
- Establish checked benchmark baselines for built-in scenarios and review
  material regressions in runtime or peak memory.
- Test supported Python versions in clean environments with pinned dependency
  bounds and deterministic fixtures.
- Exercise upload limits, malformed inputs, cancellation, stale state, cache
  isolation, privacy redaction, and report generation.
- Maintain a machine-readable release-scope summary covering supported models,
  estimators, diagnostics, sensitivities, and uncalibrated areas.
- Keep the public-beta label until the documented native validation matrix and
  release checklist are complete.

Completion gate:

- Compile, doctor, release-check, self-test, unit/integration, AppTest, privacy,
  localization, and benchmark gates all pass from a clean checkout.
- The release contains no confidential, machine-specific, or generated test
  data.
- README, in-app help, claim-to-evidence output, and release-scope summary agree
  on supported and unsupported functionality.
- No public claim depends on an external-engine run or an external result
  archive.

## Cross-cutting engineering rules

- `main` is the canonical public integration target. Experimental checkouts are
  donor workspaces only; pure-Python changes are selected by contract and test,
  not merged wholesale with rejected external-handoff work.
- New statistical computation belongs in `mfrm_app/`; Streamlit rendering
  calls it through a narrow contract.
- Canonical `AnalysisSpec`, StudyContext, rating-design, scoring-process,
  evidence, interpretation, decision, and guidance state live outside
  Streamlit widgets. Views are projections over those records, not alternate
  sources of truth.
- Refactoring and statistical changes are separate review units.
- Public tables are versioned contracts. Column removal or semantic change
  requires migration documentation and tests.
- Fixed seeds are part of the analysis identity, not hidden implementation
  details.
- Cache keys include all data, model, diagnostic, sensitivity, and code inputs
  that can change a result.
- Raw uploaded data and direct identifiers are not written to disk or included
  in reports unless the user explicitly requests an allowed export.
- Real-data entry states the deployment-specific processing and retention
  boundary before upload. Education or support artifacts exclude raw rows,
  column names, identifiers, fingerprints, free text, and person/rater labels
  by default.
- Every warning has a stable reason code, user-facing explanation, and
  recommended action.
- Every automatic recommendation states what it repairs and what it does not
  establish.
- No threshold is promoted from a descriptive convention to a universal
  acceptance rule without calibration and documentation.
- No help text, chart guide, popover, caption, success message, export prose,
  or manuscript template may express a stronger claim than the governing
  EvidenceRecord and StudyContext permit.
- Statistical rater or subgroup flags are investigation prompts. Automated
  exclusion, retraining, punishment, fairness verdicts, or high-stakes use are
  outside the core decision contract.

## Prioritized work inventory

One serial queue would postpone user, interpretation, and accessibility risk
until after most statistical work. Delivery therefore proceeds through three
parallel tracks with a shared P0 gate. Sequence is strict within a track, not
between tracks. A visible feature may ship only when its applicable scientific,
experience, and research-practice gates all pass.

### Shared P0 gate

| Outcome | Prerequisite | Principal evidence of completion |
|---|---|---|
| Canonical scope, evidence, and identity contracts | None | Boundary and contract tests |
| StudyContext, rating-design, scoring-process, guidance, and `AnalysisSpec` contracts | Scope/identity decisions | Approved schemas, transition/invalidation matrix, and no widget ownership |
| Full interpretation-copy inventory and risk classification | Evidence vocabulary | Every locale/help/chart/export surface has an owner and gate |
| Audience terminology and Help connection graph | Guidance/claim-boundary decisions | Surface-exposure matrix plus zero-orphan topic/link audit |
| Priority user tasks and SLA research scenarios | Study-context contract | Scenario matrix covering planning, first analysis, reproduction, teaching, and review |

### Track A — Scientific evidence

| Sequence | Outcome | Prerequisite | Principal evidence of completion |
|---:|---|---|---|
| A1 | Modular network and decision computation | Evidence contracts | Compatibility fixtures |
| A2 | Design perturbation and fragility audit | Modular network layer | Known-graph tests |
| A3 | Evidence-coherence computation | Normalized evidence records | Conflict fixtures |
| A4 | Conclusion-level sensitivity ledger | Analysis identity contract | Stability-state fixtures |
| A5 | Design repair candidates | Fragility audit | Feasibility and burden tests |
| A6 | Modular simulation condition system | Analysis identity contract | Serialization/fingerprint tests |
| A7 | Facet-aware stress scenarios | Simulation conditions | Recovery and flag-rate summaries |
| A8 | Native validation matrix | Frozen scenarios | Calibration reports with MCSE |

### Track B — Guided experience

| Sequence | Outcome | Prerequisite | Principal evidence of completion |
|---:|---|---|---|
| B1 | Canonical guidance state, audience vocabulary, Help registry, and safe-copy remediation | Shared P0 | Reducer, terminology/help-link/copy contracts, and locale tests |
| B2 | Optional five-node sample route with contextual Help return | Guidance and isolated tutorial context | Bilingual skip/exit/resume/stale/formative/Help AppTest journeys |
| B3 | Standard/Detailed view projections | Widget-independent `AnalysisSpec` | State-consistency and zero-computation view-switch tests |
| B4 | Current EvidenceIssue and recovery actions | Evidence-coherence computation | Conflict, aggregate-count, and repair-routing fixtures |
| B5 | Human-readable Analysis Brief and JSON sidecar | Interpretation/decision identities | Schema, privacy, language, and deterministic-generation tests |
| B6 | Accessibility and legacy-route consolidation | Proven replacement journeys | Browser matrix and unreachable-route evidence |

### Track C — SLA research and learning practice

| Sequence | Outcome | Prerequisite | Principal evidence of completion |
|---:|---|---|---|
| C1 | Study/rating/scoring records and RQ-to-design preflight | Shared P0 | Holistic/analytic and crossed/nested scenario tests |
| C2 | Missingness provenance and data-governance preflight | C1 | Reason-code, redaction, and hosted/local boundary tests |
| C3 | Claim-to-evidence and Methods/Results/Limitations matrix | Evidence and interpretation contracts | Reviewer fixtures and bounded prose tests |
| C4 | Formative sample, unfamiliar transfer case, and education companion | B2 and C3 | Learner evidence, versioned worksheet/rubric, and no-state-mutation tests |
| C5 | Privacy-safe support/reviewer bundle | B5 and C3 | Blind handoff and disclosure tests |

Public-beta hardening is a cross-track release gate. It does not require every
future item in all tracks, but every advertised claim and journey must have its
scientific validation, safe interpretation, accessibility, privacy, and clean-
environment evidence complete.

## Research validation matrix

The validation program will cross the conditions below selectively rather than
claiming that one giant factorial experiment is always feasible.

| Dimension | Representative conditions |
|---|---|
| Response model | RSM, PCM, bounded GPCM within supported scope |
| Estimator | JMLE; supported MML engine choices; fixed/free population SD where implemented |
| Performance mode | Synthetic writing, speaking, and human-scored short response within declared applicability |
| Rubric structure | Holistic; analytic criteria intended as one construct; analytic subconstructs/composites that trigger bounded or unsupported claims |
| Assignment | Fully crossed, balanced incomplete, chain, bridge-dependent, disconnected |
| Facet relation | Crossed, nested, anchored, fixed/random role declarations, and partial confounding |
| Scoring operation | Single, double, adjudicated, rescored, common-performance, and monitoring designs |
| Rater structure | Few/many levels, balanced/imbalanced workload, narrow/wide severity, session/order/drift, training cohort |
| Score support | Well populated, sparse middle category, unused declared category, extreme concentration |
| Missingness | None, planned incomplete assignment, nonresponse/not-reached, rater skip/unratable, technical/administrative loss, rater/task/subgroup-concentrated loss |
| Dependence | None, supported halo/local interaction, criterion/performance clustering, occasion drift |
| Differential interaction | Rater x task/criterion/proficiency/subgroup with overlap, multiplicity, sparse cells, and small-cell suppression |
| Person distribution | Reference normal and documented non-normal sensitivity conditions |
| Sample information | Person count, observations per level, and pairwise overlap varied separately |
| Resampling unit | Person, performance, task, rater, and justified clusters rather than rating-row bootstrap by default |

For each applicable cell, the report must identify estimand, requested and
successful replication counts, numerical failures, Monte Carlo uncertainty,
and the exact claim the result can or cannot support.

## Release decision rules

A feature is complete only when all applicable items below are satisfied:

- The statistical question and non-goals are documented.
- Applicability checks fail closed with stable reason codes.
- Deterministic unit fixtures and at least one end-to-end scenario exist.
- Sensitivity to missingness or sparse overlap is tested when relevant.
- UI, export, narrative, and evidence-ledger outputs use the same result.
- Japanese and English interfaces preserve the same statistical meaning.
- Privacy and cache behavior are reviewed.
- Performance is measured against a checked baseline.
- README/help/claim boundaries are updated.
- Current Help, tutorial, chart guides, popovers, captions, success messages,
  exported prose, and manuscript templates do not contradict or exceed the
  governing EvidenceRecord and StudyContext claim boundary.
- Every contextual Help control reaches a registered bilingual topic/section
  and returns to its source; stale or mismatched analysis context cannot appear
  as a current dynamic explanation.
- Guide and Standard surfaces use the approved audience vocabulary. Raw
  enum/class/schema/backend/environment/runner terms and IDs appear only in
  declared technical details or machine-readable artifacts.
- Production errors reveal no raw exception text, backend symbol, environment
  override, uploaded value, or identifier; they provide a stable user reason,
  reversible action, and privacy-safe support reference.
- Optional guidance can be skipped, exited, resumed, and restarted without
  changing or deleting the current analysis state.
- Workflow completion, learning completion, and ClaimID-specific readiness are
  tested as independent states; acknowledgement cannot unlock scientific or
  high-stakes use.
- Pure presentation changes, including language and view, do not alter
  AnalysisID or trigger re-estimation.
- Standard/Detailed switching preserves all hidden draft and fitted settings
  and invokes no estimator or newly requested diagnostic.
- A `HOLD` or required `NOT_ASSESSABLE` state can be archived honestly as a
  no-conclusion record but cannot produce an unbounded supported claim.
- Real-data upload and exported/support artifacts pass the applicable
  data-governance, privacy, redaction, and small-cell checks.
- Statistical flags about a person, rater, or subgroup never become an
  automated punitive or fairness decision.
- The full relevant test gate passes without suppressing unrelated failures.

## Deferred questions

The following require explicit product decisions and are not silently assumed
by this roadmap:

- Whether legacy external-comparison and advanced-model handoff surfaces should
  be deprecated or retained as clearly separated compatibility tools.
- Whether the large Streamlit entrypoint should ultimately become a thin shell
  or retain some stable compatibility wrappers.
- What native evidence would be required before considering a confirmatory
  multidimensional Python model.
- Whether guide-completion preference should persist beyond one session; this
  requires observed user need and a separate privacy/accessibility decision.
- Whether an explicit, privacy-safe learning checkpoint should be exported and
  imported for multi-session teaching before any browser persistence is added.
- Whether `Interpret` should later become an eighth top-level navigation stage;
  the InterpretationRecord remains separate even while it is grouped inside
  `Decide`.
- Which domain profiles beyond the initial synthetic SLA cases deserve
  versioned applicability rules rather than generic prose.
- When the legacy seven-section navigation and compatibility wrappers can be
  removed after Standard/Detailed view migration.
- What calibrated criteria, if any, justify removing the public-beta label.

Until those decisions are made, they cannot block the standalone core or be
used to expand its claims.
