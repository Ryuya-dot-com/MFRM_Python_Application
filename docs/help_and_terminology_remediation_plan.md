# Help and audience-language remediation plan

- Status: R1 contracts complete; bounded R2 active; fit-scatter is the sole local R3 `APPTEST_ONLY` pilot
- Last updated: 2026-07-24
- Product boundary: standalone Python
- Governing roadmap: [`../ROADMAP.md`](../ROADMAP.md)
- Semantic contract:
  [`study_context_guidance_contract.md`](study_context_guidance_contract.md)
- Copy-risk inventory:
  [`interpretation_copy_audit.md`](interpretation_copy_audit.md)

## Purpose

This plan turns the audience-language and Help contracts into file-level,
testable implementation work. It defines what changes, how it changes, the
order in which it changes, and the evidence required before legacy Help or
copy is removed.

The change is not a search-and-replace exercise. It separates three concerns:

1. scientific and SLA/measurement language that users need to learn;
2. application state and implementation language needed for reproducibility;
   and
3. navigation between the current task, a bounded explanation, and a safe
   return target.

The numerical estimators, EvidenceRecord schemas, AnalysisID computation, and
machine-readable exports do not change merely because their user-facing
projection changes.

## Implementation status — 2026-07-24

`R1` is implemented as a reversible contract layer:

- `mfrm_app/terminology.py` owns the three audience projections and their
  review-tracked definitions, aliases, boundaries, and Help references;
- `mfrm_app/help_contract.py`, `help_topics.py`, and `help_navigation.py` own
  the immutable Help graph, exact source/return identities, current-context
  verification, reserved safe fallback, strict review metadata, and
  presentation-only reducer;
- `mfrm_app/user_problems.py` owns 11 coarse problem codes, reversible actions,
  Help destinations, bounded occurrence phases, and freshly generated random
  privacy-safe support references; and
- both locale files contain the complete registered term, problem, topic, and
  navigation key sets under parity and exposure tests.

The contract registry currently contains 22 Help-home topics, one system-only
fallback topic, 42 independently declared target/focus/panel locations, 74
links, 437 bilingual topic keys, and 24 bilingual navigation keys. Every
report-wording example is paired structurally with a verify-before-use guard.
The generic unexpected-failure route is separate from estimation failure.

All new term and Help content remains marked `in_review`. Contract
completeness, bilingual parity, and automated boundary checks are not a
substitute for the planned SLA-method and UX review before copy is marked
reviewed.

The current ClaimBoundary and reference catalogs are provisional ID
allowlists, not reviewed ClaimBoundary records or bibliographic citations.
The surface catalog likewise declares intended adapter locations; it does not
prove that the current Streamlit application exposes them. The remaining R2
and R3 activation slices must resolve those records, bind each claimed
surface/focus/panel ID, and pass rendered lifecycle tests before that surface
or method copy is treated as active or reviewed.

Topic `lifecycle_states` are adapter applicability metadata; the pure
navigation reducer does not infer an application lifecycle. All built-in
topics currently remain explicitly universal (`all`). The R2 renderer calls
the contract applicability check, and its AppTest replaces one topic with a
restricted fixture to verify both denied and allowed lifecycle outcomes.

The bounded `R2` adapter now consumes these contracts for the global start
route and the first parse/estimation problem paths, including table reshaping,
custom thresholds, and covariate preview. Help owns the full main
content while the source controls remain rendered in the sidebar; this keeps
data-source state and settings alive without consuming an estimator or refit
trigger. AppTests cover no-data open, locale change, return, unknown fallback,
lifecycle denial, source-state preservation, pending triggers, and
exception-sentinel redaction. Their zero-call guard instruments both cached
core functions and direct module-level estimation/refit paths.

The first bounded `R3` adapter now additionally consumes the exact
`link.popover.fit_scatter` edge. It opens
`help.results.fit#fit-scatter` from the real fit-scatter figure and returns to
the exact registered target/focus pair. The adapter requires the retained
AnalysisID and stored sample/real provenance, restores the Essential or Full
fit-details projection, and separately restores only the figure-owning
`fit_df_method_method` and `fit_df_method_cap` values through an independent
key/value allowlist. Native MML prior-SD sensitivity drafts and misfit-ranking
view controls use separately validated, AnalysisID-scoped durable session
state so main-surface widget cleanup cannot silently reset them or transfer
them to a new fit. The imported-FACETS table comparator has been removed from
the default Fit Details journey and remains only as an isolated legacy
compatibility function. Missing, altered, stale, or unverifiable intent fails
closed; fitted phase, AnalysisID, and sample/real provenance are rechecked
before any return value is restored. A new Help journey clears an abandoned
return intent and its transient status/suppression state.

The adapter is explicitly `APPTEST_ONLY`, with no browser evidence ID. A
one-pilot activation gate permits this first slice while requiring every other
active adapter to be `BROWSER_ACCEPTED` with a stable retained evidence
reference. The required cases and promotion boundary are defined in
[`help_browser_acceptance.md`](help_browser_acceptance.md); the existence of
that runbook is not a passing browser result.

This is not completion of R2 or R3. Guided routes and the remaining
provisional surface/focus pairs are not yet authoritative; of the 18
registered core popover aliases, only `fit_scatter` has an active full-guide
edge. The two legacy posterior popovers remain outside the core migration.
The cached Python core namespace and source controls may be loaded while Help
is open so Streamlit can retain their widget state, but no estimation or refit
is initiated by Help navigation.

## Decisions fixed by this plan

1. The application has one authoritative Help content system. Tutorial,
   Quick Start, chart popovers, warnings, glossary, and result Help are
   projections of that system rather than independent prose stores.
2. Static Help is available before data input and remains available through
   mapping, preflight, execution, first failure, `HOLD`, stale results, and a
   completed fit.
3. Help opens as a full main-content surface from a persistent app-level
   control. A modal is not the primary route. This avoids a second navigation
   hierarchy and gives long method content a stable heading and reading order.
4. A chart popover remains a short inline explanation. Its final control opens
   the exact detailed topic and records a stable return target.
5. Help routing uses stable IDs. English titles, Japanese labels, screen-name
   substrings, and free-text search results never act as routing identity.
6. Guide/Standard, Detailed/method, and Technical are content layers. Technical
   is not a competing analysis mode and never implies a stronger result.
7. A Help action changes no data, analysis setting, AnalysisID,
   EvidenceRecord, diagnostic request, claim state, or learning-completion
   state.
8. Active core Help contains no workflow to TAM, R, Julia, Stan, Posterior
   Viewer, or a cross-engine package. Publications and software may appear as
   bounded method references, not as execution or validation actions.
9. Unknown or unavailable Help topics remain visible as a localized safe
   fallback and fail a registry test. They never silently disappear.
10. No public runtime toggle is introduced to maintain two competing Help
    systems. Migration is reversible at the commit/slice level until the
    legacy path is deleted.

## Current-to-target change map

The `Current behavior` column records the pre-remediation baseline used to
define this plan; it is not a live inventory of every local change. As of the
status date, `HELP-01` and `HELP-03` are active, and `ERR-01` is active only for
paste/upload parsing and primary estimation failures. Other rows remain
partial or planned unless the implementation-status section says otherwise.

| Work ID | Current behavior | Target behavior | Main implementation mechanism |
|---|---|---|---|
| `LANG-01` | Guide and Help expose class names, enums, filenames, backend terms, and environment variables | Guide/Standard use task, meaning, limit, and action language; exact terms remain in declared technical surfaces | Versioned terminology registry plus surface allowlist |
| `LANG-02` | The persistent product badge and run breadcrumb foreground runtime, source commit, and input fingerprint | The default badge states the supported Python product boundary; build and input-match details move to an explicit technical/reproducibility disclosure | Audience projection for app identity and current-run identity |
| `HELP-01` | Full Help is reached only through fitted-result routes | Static Help is opened from the app shell before every data/result-dependent return | Persistent Help launcher and top-level Help route |
| `HELP-02` | Popover footer says to find another Help panel | One action opens the exact topic/section and returns to the source | `HelpLink` registry and session-scoped route reducer |
| `HELP-03` | Unknown popover topics silently return | Localized fallback, Help-home action, and failing contract test | Fail-closed topic resolution |
| `HELP-04` | English titles and translated substrings drive routing | Stable topic, section, source, and target IDs drive routing | Topic and target registries |
| `HELP-05` | Help is a long panel containing workflow, developer, audit, and reviewer tables | Help home is organized by research task; technical/reviewer material is an explicit deeper layer or companion artifact | Topic template and audience projection |
| `HELP-06` | Static and current-analysis explanations are mixed | Static method content always renders; dynamic content renders only after identity verification | `HelpContextBinding` validation |
| `ERR-01` | Raw exception type/message and a developer environment variable are shown | Localized problem, preserved-state statement, reversible action, exact Help link, and privacy-safe support reference | Stable `UserProblem` records and server-side logging |
| `COPY-01` | Quick Start teaches EvidenceRecord, AnalysisID, runner, and file names | Quick Start teaches the shortest safe task sequence | New topic copy; technical reproducibility disclosure remains optional |
| `COPY-02` | Analysis Workflow teaches internal state names and a large audit inventory | Analysis path explains evidence availability, limits, and next actions in user language | Topic projections over evidence states |
| `COPY-02A` | The estimation tutorial tells users to archive record classes, raw IDs, and an engine contract | The tutorial explains what must be saved and why; exact artifact names appear only after a technical-detail request | Reuse the repeat-analysis Help topic instead of hard-coded tutorial prose |
| `COPY-03` | Residual-PCA errors expose function names, NumPy details, and raw exceptions | Reason-specific, bounded messages state whether the screen was not run, unsupported by the data, or failed numerically | Problem/reason mapping with no raw exception interpolation |
| `COPY-04` | Resource preflight advertises a safety-override environment variable | Hosted limit, retained state, supported scope reductions, and local-run option are explained | Localized resource problem plus target actions |
| `COPY-05` | Downloads describe hidden panels, reruns, package generation, and raw manifest names | Downloads describe what will be prepared, expected cost, and what is included/excluded | User-facing artifact labels plus technical filename disclosure |
| `COPY-06` | Active Help sends users to external engines or Posterior Viewer | Native Python settings, diagnostics, sensitivity, and archives are the only core actions | Remove active calls to action; retain bounded citations only |
| `COPY-07` | “Full table” mixes substantive results with contract columns | Substantive detailed table and technical audit table/export are separate | Column exposure registry |
| `COPY-08` | Glossary definitions live in parallel Python and locale structures | Inline definition, Help glossary, search aliases, and export derive from one term registry | Canonical concept/term source |
| `COPY-09` | Keyboard Help teaches Streamlit rerun/cache mechanics and a FACETS-mode button name | No dedicated shortcut surface; native control labels, visible focus, and persistent navigation carry the interaction contract | Remove the sidebar projection and retain browser acceptance for ordinary Tab/Enter operation |

## Target architecture

```text
source card / warning / setting / chart
                 |
                 v
          HelpLink registry
       (stable source + topic + return)
                 |
                 v
          Help route reducer --------------------+
                 |                               |
                 v                               |
          HelpTopic registry                     |
     (static content + claim boundaries)         |
                 |                               |
                 +--> verified context? -- yes --+--> “For this analysis”
                 |                    \\-- no -----> static explanation only
                 v
        Standard / Detailed / Technical projection
                 |
                 v
     return target registry -> restored panel/card/focus
```

The registries and reducer are pure Python. Streamlit translates locale keys,
stores the session-scoped navigation state, and renders controls. Pandas,
Streamlit, estimator imports, network calls, and external runtimes are not
allowed in the contract modules.

## File-level implementation

### New pure modules

| File | Responsibility | Must not own |
|---|---|---|
| `mfrm_app/terminology.py` | `AudienceLayer`, `TermSpec`, surface exposure, aliases, definitions, prohibited uses, and registry validation | Streamlit rendering, locale prose, analysis state |
| `mfrm_app/help_contract.py` | Immutable `HelpTopic`, `HelpLink`, and `HelpTarget` schemas plus graph and field validation | Widget keys, translated display strings, route mutation, diagnostics |
| `mfrm_app/help_topics.py` | Reviewed topic/link/target registry instances, stable sections, concept/claim-boundary references, and legacy exact aliases during migration | Display-string routing or runtime analysis values |
| `mfrm_app/help_navigation.py` | `HelpRouteState`, context verification, and the pure open/change/return/invalidate reducer | Streamlit session access or scientific/learning state |
| `mfrm_app/user_problems.py` | Stable problem codes, severity, locale keys, safe actions, Help topic, and privacy-safe support-reference construction | Rendering raw exception text or deciding scientific claim readiness |

These modules follow the immutable, versioned style already used by
`mfrm_app/evidence.py`. They expose deterministic validation and serialization
where an exported contract is needed, but Help navigation state itself is not
inserted into AnalysisSpec or EvidenceRecords.

### Existing files

| File | Required change |
|---|---|
| `mfrm_app/help_popovers.py` | Become a temporary compatibility projection from registered Help topics. Remove the independent Japanese-only topic store after every caller uses the registry. |
| `streamlit_app.py` application shell | Render the persistent Help control and active Help surface after locale/title/privacy setup. While Help owns the main content, keep source controls rendered in the sidebar and block every run/refit trigger. |
| `streamlit_app.py` product/run identity | Replace runtime/commit/fingerprint-first captions with a short product boundary and current-run summary; expose build and match values only under registered technical details. |
| `streamlit_app.py` Help section | Replace `show_help_section()` title-string branching with topic-ID selection and projections. Keep a narrow compatibility wrapper only while callers migrate. |
| `streamlit_app.py` contextual Help | Replace `_HELP_POPOVER_LIBRARY` lookup and caption-only footer with a registered inline projection plus an actionable detailed-Help button. |
| `streamlit_app.py` guided routing | Replace `guided_section_id_for_target()` translated-label and substring matching with exact `HelpTarget` mappings. |
| `streamlit_app.py` errors | Replace `render_exception_details()` on public paths with `render_user_problem()`. Raw exceptions remain in confidential logs only; no environment variable enables a public exception panel. |
| `streamlit_app.py` PCA reasons | Store stable reason/problem codes instead of exception names/messages in locale arguments. Completed results remain available when a diagnostic fails. |
| `streamlit_app.py` tables/reports | Project user labels and substantive columns for Standard/Detailed; retain exact contract columns in technical audit CSV/JSON. |
| `streamlit_app.py` tutorial/navigation | Keep repeat-analysis guidance task-focused; remove shortcut instructions and keep result navigation visible during scrolling. |
| `locales/en.json`, `locales/ja.json` | Add matched `terms`, `help_nav`, `help_topics`, and `problems` keys. Remove an old key only after no production/test reference remains. |
| `tests/` | Add contract, route, lifecycle, copy exposure, privacy, and AppTest coverage described below. |

### Current code touchpoints

| Current location | Exact modification | Acceptance evidence |
|---|---|---|
| [`streamlit_app.py`](../streamlit_app.py#L68369), `main()` | Add the global Help launcher and active-route render immediately after the localized app shell; keep the source sidebar alive but suppress source main content and analysis triggers | Help AppTest preserves no-data input/settings and records zero estimator/refit calls |
| [`streamlit_app.py`](../streamlit_app.py#L22965), `render_app_scope_badges()` | Default caption becomes a localized product/version/boundary summary; source commit and runtime wording move to technical details | Standard-surface copy scan has no `source commit`/runtime implementation label |
| [`streamlit_app.py`](../streamlit_app.py#L2651), `_render_run_breadcrumb()` | Keep model, estimator, convergence, and data-size orientation; move fingerprint to reproducibility details | Restored runs remain distinguishable without a raw match value in Standard |
| [`streamlit_app.py`](../streamlit_app.py#L22280), `_standalone_estimation_tutorial_markdown()` | Replace AnalysisID/EvidenceRecord/engine-contract instructions with the user purpose of saving settings, evidence, and limits; link to repeat-analysis Help | Tutorial and Quick Start use the same topic/ClaimBoundaryIDs |
| [`streamlit_app.py`](../streamlit_app.py#L23521), `_render_compact_dataframe()` | Require both compact and detailed column allowlists; never infer “detailed” as every remaining DataFrame column | Contract columns cannot appear merely by expanding a substantive table |
| [`streamlit_app.py`](../streamlit_app.py#L42399), `_render_final_readiness_section()` | Standard table shows localized check, state, evidence summary, next action, and required status; enum/code columns remain in technical CSV/JSON | AppTest shows no raw computation/stability/reason fields; export schema test remains stable |
| [`streamlit_app.py`](../streamlit_app.py), `render_user_problem()` and `render_exception_details()` | Keep confidential logging capability but replace public rendering with `UserProblem`; a legacy environment variable cannot re-enable raw exception output | Sentinel exception content is absent from the DOM and support output even when the legacy variable is set |
| [`streamlit_app.py`](../streamlit_app.py#L23424), paste/upload parse failures | Render `input.parse_failed`, format/delimiter action, exact Help link, and support reference; do not open a raw technical expander | Paste and upload failure AppTests preserve input and show the same bounded problem contract |
| [`streamlit_app.py`](../streamlit_app.py#L32958), `_ESTIMATION_ERROR_PATTERNS` | Convert `(diagnosis, action, keyword)` strings into match rule -> stable problem code; locales and targets own visible remediation | Pattern tests assert codes and targets rather than English prose fragments |
| [`streamlit_app.py`](../streamlit_app.py#L34389), estimation exception path | Render the classified problem, preserved-state statement, return-to-settings action, and troubleshooting link; keep unclassified errors generic | First-failure Help works without refit and no raw exception is displayed |
| [`streamlit_app.py`](../streamlit_app.py#L21915), residual-PCA exception path | Store stable diagnostic reason/problem codes without `{error_type}` or `{error}`; retain other completed diagnostics | PCA failure fixture renders a bounded reason in both locales |
| [`streamlit_app.py`](../streamlit_app.py#L29199), `render_estimation_resource_preflight()` | Replace environment-variable bypass copy and internal metric labels with user consequence, safe reductions, retained state, and Help action | Hosted-limit test finds no bypass variable and can return to the exact setting |
| [`streamlit_app.py`](../streamlit_app.py#L26209), `guided_section_id_for_target()` | Stop using it on new paths; replace substring inference with registered target IDs. A temporary adapter accepts only enumerated exact legacy aliases | Locale/title-renaming tests do not change the resolved section |
| [`streamlit_app.py`](../streamlit_app.py#L53343), `show_help_section()` | Turn into a compatibility wrapper over registered topic/category projections; remove English-title `if` identity and developer inventory from the default path | Topic graph is zero-orphan and Help selection works after display-title changes |
| [`streamlit_app.py`](../streamlit_app.py#L67147), `render_help_popover()` | Resolve an explicit HelpLink; render shows/reads/limits inline; add detailed-Help button; unknown IDs render safe fallback | Every active caller has one valid link and one return/focus target |
| [`streamlit_app.py`](../streamlit_app.py), removed `render_keyboard_shortcuts_help()` | Remove the shortcut surface and its cache/rerun guidance; do not replace it with custom hotkeys | Initial AppTest has no shortcut expander and both result projections retain native focusable navigation |
| [`streamlit_app.py`](../streamlit_app.py#L58713), `_render_downloads()` | Replace implementation-mechanics captions, give privacy manifest a user-purpose label, and keep filenames/renderers in technical details | Sharing-mode AppTest names included/excluded content without raw manifest vocabulary |
| [`mfrm_app/help_popovers.py`](../mfrm_app/help_popovers.py#L1) | During migration, generate legacy inline fields from registered topic sections; then remove the independent translation mapping | One canonical content source supplies both locales and every projection |
| [`locales/en.json`](../locales/en.json#L798) and [`locales/ja.json`](../locales/ja.json#L798) | Remove function/NumPy/exception placeholders from visible dimensionality reasons | Placeholder parity passes and no raw exception placeholder remains in a public key |
| [`locales/en.json`](../locales/en.json#L986) and [`locales/ja.json`](../locales/ja.json#L986) | Replace the hosted-limit bypass instruction | Resource-limit copy has consequence, choices, and exact Help target |
| [`locales/en.json`](../locales/en.json#L2357) and [`locales/ja.json`](../locales/ja.json#L2357) | Replace hidden-panel/rerun/package-generation explanation and raw manifest reference | Download copy states what is prepared and what is safe to share |

### Application-shell order

The top-level order becomes:

```text
configure page
initialize language and presentation-only state
render title, product boundary, privacy notice, and persistent Help control
if Help is active:
    render the selected Help topic and return controls
    stop this render without reading data or invoking analysis
render optional guide/onboarding
read or select data
if no data:
    show the no-data action and stop
render preflight and analysis workspace
```

This order makes Help independent of fitted output. Opening Help from a result
does not discard the result; returning rerenders the preserved result state and
does not press the estimation control.

## Contract details

### `AudienceLayer` and `TermSpec`

The initial layers are:

| Stable layer | Purpose | Examples |
|---|---|---|
| `standard` | Task, status, meaning, boundary, and next action | “Analysis settings,” “why this is limited,” “pause and fix” |
| `method` | Canonical SLA/measurement term with definition, applicability, and provenance | residual PCA, Infit/Outfit, separation, JMLE/MML |
| `technical` | Support and machine-contract identity | `AnalysisID`, `ReasonCode`, schema version, exact filename |

`TermSpec` contains at least:

| Field | Meaning |
|---|---|
| `concept_id` | Stable concept identity |
| `standard_label_key` | Plain localized label |
| `method_label_key` | Canonical localized method label |
| `technical_label` | Exact machine term when one exists |
| `definition_key` | Bounded definition |
| `aliases_keys` | English, Japanese, acronym, and corrected legacy search aliases |
| `allowed_layers` | Surfaces on which the term may appear |
| `prohibited_claim_keys` | Interpretations the term must not imply |
| `related_help_topic_ids` | Topics that explain the concept |
| `version_and_review` | Content version, owner role, last review, and review state |

Scientific terms are not treated as implementation leakage. MFRM, RSM, PCM,
GPCM, JMLE, MML, logit, Infit, Outfit, residual PCA, and separation remain
available with definitions and boundaries. Backend function names, cache keys,
payload/schema vocabulary, environment variables, and raw IDs are not part of
the Standard conceptual model.

### `HelpTopic`

Each topic uses stable locale keys rather than embedding prose in Python and
contains:

| Field | Requirement |
|---|---|
| `help_topic_id` | Stable language-independent ID |
| `title_key`, `summary_key` | Localized heading and search summary |
| `concept_ids` | Concepts explained by the topic |
| `audience_layers` | Available Standard/method/technical projections |
| `lifecycle_states` | States in which the topic is available; static topics normally include all states |
| `applicability` | Models, estimators, designs, rubric structures, and context conditions |
| `prerequisite_keys` | What must be known before applying the explanation |
| `computed_key` | What the application computes |
| `can_show_key`, `cannot_show_key` | Claim boundary shown together |
| `next_check_key`, `next_action_key` | Evidence inspection and reversible action |
| `safe_report_key`, `avoid_report_key` | Bounded reporting example and wording to avoid |
| `claim_boundary_ids` | Required scientific boundaries |
| `reference_ids` | Versioned method references |
| `related_target_ids` | Valid application destinations |
| `search_alias_keys` | Bilingual and corrected legacy terms |
| `version_and_review` | Owner role, reviewers, version, date, and state |

### `HelpLink`, `HelpTarget`, and route state

`HelpLink` contains only stable relationships:

| Field | Requirement |
|---|---|
| `help_link_id` | Unique link ID |
| `source_target_id` | Exact calling card, chart, warning, setting, or report section |
| `help_topic_id`, `section_id` | Exact destination |
| `return_target_id`, `return_focus_id` | Exact return location and focus marker |
| `context_policy` | `static_only`, `optional_current`, or `required_current` |

`HelpTarget` maps a stable target to presentation state. The Streamlit adapter,
not the registry, translates that state into existing widget keys such as
`guided_essential_section`, `guided_diagnostics_panel`,
`guided_figures_panel`, `main_results_panel`, and `downloads_panel`.

`HelpRouteState` is session-scoped presentation state and contains the selected
topic/section, origin, return target/focus, locale-independent history, and an
optional context binding. It is explicitly excluded from AnalysisSpec,
AnalysisID, EvidenceRecords, interpretation/decision readiness, and learning
completion.

### Context binding

Static “About this method” content always renders. “For this analysis” renders
only when all supplied current references match the current application state:

- sample versus real-data context;
- StudyContextID when one exists;
- AnalysisID;
- EvidenceID or EvidenceIssueID;
- ClaimBoundaryIDs; and
- fitted-versus-draft status.

A stale, absent, or mismatched binding suppresses the dynamic paragraph and
states that general method Help is being shown. It does not try to reconstruct
old analysis state, compute a missing diagnostic, or silently bind to the most
recent run. Runtime context references remain in session state and are not put
in URLs, search text, telemetry, or support references.

### `UserProblem`

The problem contract contains:

| Field | Requirement |
|---|---|
| `problem_code` | Stable coarse code unrelated to exception class names |
| `severity` | Information, caution, or blocked action |
| `title_key`, `body_key` | Localized problem and preserved-state explanation |
| `action_keys` | One or more reversible actions |
| `action_target_ids` | Exact controls/screens for those actions |
| `help_topic_id` | Exact troubleshooting topic |
| `support_reference` | Random incident reference; never a hash of data, filename, identifier, or exception text |

Exception matching may remain behind the adapter, but it returns a stable
problem code rather than user prose. The full exception is written to the
configured server/local log with the same support reference. Hosted UI output
contains no exception message, stack, path, uploaded value, column value,
environment-variable instruction, or identifier.

Initial problem mapping:

| Problem code | User meaning | Primary target |
|---|---|---|
| `input.parse_failed` | The pasted/uploaded table could not be read | input format and delimiter |
| `input.mapping_invalid` | Required analysis roles are not mapped | column mapping |
| `input.no_usable_rows` | No usable rating rows remain | missingness/data audit |
| `estimation.identification_failed` | The current design/constraints do not identify the requested fit | design and constraints |
| `estimation.nonconvergence` | Estimation stopped before the stopping criterion was met | estimation settings and data audit |
| `estimation.rating_scale_invalid` | Scores do not define the requested category structure | score/category mapping |
| `estimation.anchor_invalid` | An anchor does not match the current data/specification | anchor settings |
| `resource.hosted_limit` | The requested work exceeds the hosted safety envelope | scope reduction or controlled local use |
| `diagnostic.residual_pca_unavailable` | The diagnostic was not requested or its data prerequisites are absent | residual-structure Help |
| `diagnostic.residual_pca_failed` | The residual matrix could not be analyzed reliably | residual-structure troubleshooting |
| `problem.unexpected` | The step failed for an unclassified reason | safe generic troubleshooting |

## Help information architecture

The Help home uses six task-centered categories:

1. Get started / はじめる
2. Choose an analysis / 分析を選ぶ
3. Understand results / 結果を理解する
4. Resolve a problem / 問題を解決する
5. Prepare a report / 報告を準備する
6. Methods and glossary / 方法と用語集

Initial stable topics are:

| Topic ID | Main question |
|---|---|
| `help.get_started.overview` | What is the shortest safe path through the app? |
| `help.data.long_format` | What must one row represent? |
| `help.data.mapping` | Which columns play which roles? |
| `help.data.missingness` | Which ratings were planned absent versus unexpectedly missing? |
| `help.design.coverage` | Is the observed rating assignment connected and informative for the intended comparison? |
| `help.run.large_analysis` | Why is this run blocked or expensive, and what can be changed safely? |
| `help.run.estimation_failed` | What remains available after an estimation failure, and what should be checked next? |
| `help.run.nonconvergence` | What does non-convergence mean and not mean? |
| `help.results.first_read` | Which evidence should be read first? |
| `help.results.measures_targeting` | How are measures, uncertainty, and targeting read? |
| `help.results.fit` | What do fit screens show and not show? |
| `help.results.categories` | What do category use and thresholds show? |
| `help.results.residual_structure` | What can residual PCA and related screens detect? |
| `help.results.rater_evidence` | How are severity, agreement, and differentiation kept distinct? |
| `help.results.differential_interaction` | What does an interaction flag mean, and why is it not a fairness verdict? |
| `help.report.claim_limits` | Which statements are supported, limited, or unavailable? |
| `help.downloads.privacy` | What is included or excluded from a sharing package? |
| `help.downloads.repeat_analysis` | What must be saved to repeat this analysis in Python? |
| `help.methods.rsm_pcm_gpcm` | How do supported response models differ? |
| `help.methods.jmle_mml` | How do supported estimators and person-distribution choices differ? |
| `help.glossary` | What do the canonical SLA/measurement terms mean? |

Search uses localized aliases and concept IDs but returns localized titles and
summaries. Searching an unsafe legacy phrase such as “perfect fit” or “bias”
routes to the corrected bounded topic; it does not preserve the unsafe claim.

### Existing popover migration

| Existing key(s) | Registered topic/section |
|---|---|
| `scree` | `help.results.residual_structure#residual-pca` |
| `fit_scatter`, `pathway_map`, `misfit_ranking`, `zstd_distribution`, `qq_residuals` | `help.results.fit` with figure-specific sections |
| `category_probability`, `category_usage`, `threshold_map` | `help.results.categories` with figure-specific sections |
| `coverage_heatmap` | `help.design.coverage#observed-coverage` |
| `wright_map`, `ecdf_measures`, `forest_measures`, `facet_distribution` | `help.results.measures_targeting` with figure-specific sections |
| `bias_heatmap`, `classical_dif` | `help.results.differential_interaction` with method-specific sections |
| `rater_agreement` | `help.results.rater_evidence#agreement` |
| `mml_person_sd` | `help.methods.jmle_mml#person-distribution-sd` |
| `posterior_trace`, `posterior_rhat_ess` | No core migration; isolate with the legacy posterior route and remove from the active core registry |

The compatibility function may accept an old key during migration, but it
resolves through an explicit alias map. A substring or translated-title match
is never permitted.

## Lifecycle behavior

| Application state | Help behavior | Return target |
|---|---|---|
| Landing/no data | Static Help, data format, privacy, and sample route are available; no “For this analysis” block | input start |
| Data loaded/pre-run | Static Help plus safe summaries of current mapping/preflight when explicitly linked | originating mapping or preflight card |
| Preflight blocked | Explain the limit, retained state, safe scope reductions, and what the block does not imply | exact blocked warning/control |
| Running | Static Help only; opening Help does not cancel, restart, or change the requested work | progress surface |
| First failure | Preserve the localized problem and support reference; Help opens without rerunning | error heading and next-action controls |
| `HOLD` | Explain the failed prerequisite, affected claim, still-available evidence, and repair | originating evidence row/card |
| `CAUTION` | Explain the limitation and next check without converting it to a pass/fail result | originating figure/table/card |
| Fitted/current | Add current-analysis explanation only after identity verification | originating result target |
| Stale/mismatched | Show static content and an explicit stale-context note; do not show a current-analysis assertion | stale warning or source target |

The return adapter restores the Standard/Detailed projection and selected
section/panel. A target adapter may additionally restore only independently
allowlisted, type/range-validated state owned by that exact source. In the
first fit-scatter slice this means the fit d.f./ZSTD display method and cap.
Native MML-sensitivity drafts and misfit-ranking view controls persist through
their own validated backing state rather than being placed in the return
intent. The external-comparison uploader is not part of the active standalone
Fit Details journey. The adapter anchors the real localized Fit scatter `H3`,
renders a non-live localized return caption, and emits a registry-owned static
focus/scroll request only after the exact target consumes a valid intent. The
script contains no AnalysisID, provenance, locale copy, or input data and
fails silently without changing scientific or presentation state. AppTest can
verify the heading anchor, controller body, JavaScript opt-in, one-shot
emission, and fail-closed absence, but not execution. Exact keyboard focus,
viewport placement, focus outline, Tab order, CSP behavior, zoom,
reduced-motion behavior, and screen-reader announcement remain browser
acceptance requirements.

## Concrete audience copy

### Term projections

| Technical form | Standard English | Standard Japanese | Method/technical disclosure |
|---|---|---|---|
| `AnalysisID` | Analysis reference | 分析参照番号 | Full ID is copyable only in technical/reproducibility details |
| `EvidenceRecord` | Evidence summary | 根拠の要約 | Exact record and ID remain in JSON/technical audit |
| `ReasonCode` | Why this is limited | 制限される理由 | Exact code remains in support/JSON |
| `ComputationState` | Evidence availability | 根拠の確認状況 | Enum remains technical |
| `StabilityState` | Conclusion across planned checks | 計画した確認での結論 | Enum remains technical |
| app-engine runner | Re-run this analysis in Python | この分析をPythonで再実行 | Script filename is secondary technical detail |
| fingerprint | Input match value | 入力照合値 | Algorithm/full value stays technical |
| manifest | Included and excluded files | 含まれるファイル・除外されたファイル | Exact manifest filename stays technical |
| cache/session state | Kept for this browser session | このブラウザセッション中に保持 | Implementation mechanism stays technical |

### State projections

| Internal state | English label | Japanese label | Required accompanying text |
|---|---|---|---|
| `AVAILABLE` | Available | 利用可能 | Computed/applicable is not proof of a broader claim |
| `CAUTION` | Available with caution | 注意付きで利用可能 | Name the limitation and next check |
| `HOLD` | Pause and fix before interpreting | 修正するまで解釈を保留 | Name the prerequisite and repair |
| `NOT_ASSESSABLE` | Cannot be assessed with the current design | 現在の設計では評価不能 | Distinguish from not yet assessed |
| `STABLE` | Stable across the planned checks | 計画した確認で安定 | Name the prespecified checks |
| `CONDITIONALLY_STABLE` | Stable under the stated conditions | 条件付きで安定 | Name the conditions |
| `SENSITIVE` | Conclusion changes across settings | 設定により結論が変化 | Name the influential setting/check |
| `NOT_ASSESSED` | Stability not checked | 安定性は未確認 | Provide the action that requests the check |

### Quick Start replacement

English:

> Before reporting, review every item marked “Needs attention.” Then save the
> tables you need, the analysis settings, and the reporting checklist from
> Downloads.

Japanese:

> 報告前に「要確認」の項目をすべて見直します。その後、「ダウンロード」
> から必要な表、分析設定、報告前チェックリストを保存します。

Quick Start does not teach record class names, raw IDs, runner terminology, or
filenames. An optional reproducibility disclosure explains their purpose after
the user asks for technical detail.

### Analysis path replacement

English:

> The app distinguishes four situations: evidence is available, available
> with caution, interpretation should pause until a problem is fixed, or the
> question cannot be assessed with the current design. Before reporting a
> conclusion, confirm that it belongs to the analysis currently shown.

Japanese:

> 本アプリは、結果を「利用可能」「注意付きで利用可能」「問題を修正する
> まで解釈を保留」「現在の設計では評価できない」に分けて示します。結論
> を報告する前に、現在表示している分析の結果であることを確認してください。

### Estimation failure replacement

| Element | English | Japanese |
|---|---|---|
| Heading | Estimation did not finish | 推定を完了できませんでした |
| Body | No result was created for this attempt. Your data and settings are still available. Try the suggested action below. | 今回の実行結果は作成されていません。データと設定は保持されています。以下の対処を確認してください。 |
| Help action | Open Help: estimation did not finish | ヘルプを開く：推定を完了できない場合 |
| Return action | Return to settings | 設定に戻る |
| Support action | Show support reference | 問い合わせ用情報を表示 |

### Hosted resource-limit replacement

English:

> This analysis is larger than this hosted session can run safely. Reduce the
> data or analysis scope, or run the app in a controlled local environment.

Japanese:

> この分析は、現在のホスト環境で安全に実行できる規模を超えています。
> データまたは分析範囲を縮小するか、管理されたローカル環境で実行して
> ください。

No environment-variable bypass appears in the user interface.

### Residual-PCA reason replacements

| Situation | English | Japanese |
|---|---|---|
| Not requested | Residual PCA was not included in this run. Request the diagnostic and run the analysis again if it is required for the study question. | この実行では残差PCAを実施していません。研究課題に必要な場合は診断を選択して再実行してください。 |
| Insufficient comparable profiles | Residual PCA is not available because too few persons have comparable residual profiles. | 比較可能な残差パターンを持つPersonが少ないため、残差PCAを利用できません。 |
| Sparse comparisons | Too many residual comparisons are missing. Review observation coverage and connectedness before rerunning. | 残差間の比較に必要な観測が不足しています。再実行前に観測カバレッジと連結性を確認してください。 |
| Numerical failure | The residual correlation matrix could not be analyzed reliably. Other completed results remain available. | 残差相関行列を安定して分析できませんでした。すでに完了した他の結果は引き続き利用できます。 |

The bounded no-signal statement is “This screen did not identify a strong
secondary residual contrast” / 「このスクリーニングでは、強い副次的な
残差コントラストは確認されませんでした」. It never says that residual
PCA proved or confirmed unidimensionality.

### Downloads replacement

| Surface | English | Japanese |
|---|---|---|
| Panel selector | Choose a download category. Files are prepared when you open that category. | ダウンロードする種類を選んでください。ファイルは、その項目を開いたときに準備されます。 |
| Sharing mode | Sharing mode is on. Individual-response and person-level files are left out. Review “What is included in this download” before sharing. | 共有用モードが有効です。個別反応およびPerson単位のファイルは除外されます。共有前に「このダウンロードに含まれるもの」を確認してください。 |
| Manifest label | Download contents and privacy check | ダウンロード内容とプライバシー確認 |
| Figure limitation | PNG files are unavailable in this environment. Interactive HTML figures remain available. | この環境ではPNGを作成できません。インタラクティブHTMLは利用できます。 |

The manifest filename, renderer, browser dependency, and file-generation
mechanics remain optional technical detail.

## Implementation slices

### `R0` — Freeze reachable behavior and risk

Changes:

- complete the copy manifest with surface, lifecycle, audience layer, owner,
  Help topic, return target, and risk;
- record current no-data, pre-run, failure, `HOLD`, and post-fit Help behavior;
- add characterization tests before moving renderers; and
- classify legacy-only locale keys and functions without deleting them.

Gate:

- every reachable Help/copy surface is present in the manifest;
- current failures are captured by tests rather than accepted as target
  behavior; and
- no estimator or output schema changes are included.

### `R1` — Add pure contracts and safe seed content

Changes:

- add terminology, Help, target, route, and user-problem contracts;
- register the initial topics required for landing, data mapping, resource
  limits, estimation failure, result first read, residual structure, and
  download privacy;
- add matched locale keys; and
- validate IDs, fields, references, aliases, layers, and claim boundaries.

Gate:

- registries are deterministic, bilingual, zero-orphan, and Streamlit-free;
- all links and targets resolve; and
- unsafe legacy copy is not reused as seed content.

### `R2` — Make Help persistent and errors safe

Changes:

- insert the Help launcher/active route before data-dependent returns;
- render safe static topics without AnalysisID;
- replace the three public `render_exception_details()` call sites with
  `UserProblem` rendering;
- remove environment-variable instructions and raw PCA error interpolation;
  and
- connect parse, preflight, and estimation failures to exact Help topics.

Gate:

- Help opens from no-data, pre-run, preflight block, and first-failure states;
- synthetic exception text containing names, paths, columns, and values is
  absent from rendered UI;
- data/settings remain intact; and
- Help open/return causes zero estimator calls.

### `R3` — Connect result surfaces and return paths

Changes:

- migrate each current popover through the explicit alias table;
- add “Open full guide” and localized return controls;
- replace translated/substring route resolution with target IDs;
- verify optional current-analysis bindings; and
- add stale/sample-real mismatch handling.

Gate:

- every contextual control reaches one exact registered section in one action;
- unknown/unavailable topics show a safe fallback;
- Standard and Detailed reach the same concept/ClaimBoundaryIDs; and
- return restores the original panel/card and passes keyboard/browser checks.

Current local activation record (not an R3 completion claim):

- active popover set: exactly `{fit_scatter}`;
- acceptance state: `fit_scatter = APPTEST_ONLY`, browser evidence ID = none,
  and AppTest-only pilot set = exactly `{fit_scatter}`;
- exact route: `link.popover.fit_scatter` ->
  `help.results.fit#fit-scatter` ->
  `target.results.figure.fit_scatter` /
  `focus.results.figure.fit_scatter`;
- verified by AppTest: English/Japanese retention, current AnalysisID and
  sample/real provenance binding, Essential/Full return projections,
  fit-display method/cap restoration after widget cleanup, native
  MML-sensitivity and misfit-ranking draft retention, zero estimator or refit
  calls, pending-trigger preservation, positive return after browsing another
  static topic, and pre-restore fail-closed behavior for missing, malformed,
  stale, cross-analysis, or abandoned intents; plus the same stable native
  `H3` anchor in both locales, exact static JavaScript controller, one-shot
  emission, non-live returned caption, and no controller for rejected intents;
- intentionally inactive: the other 17 registered core popover aliases and
  both legacy posterior popovers; and
- still pending: real-browser `activeElement`, viewport/toolbar placement,
  focus outline and following Tab order, CSP, heading/reading order, responsive
  zoom, reduced-motion, and single-announcement screen-reader acceptance.

No second AppTest-only adapter may enter the active set. After every case in
[`help_browser_acceptance.md`](help_browser_acceptance.md) passes, the reviewed
record receives a stable `browser.*` reference, `fit_scatter` may be promoted
to `BROWSER_ACCEPTED`, and one next exact popover may become the sole pilot in
a separate change.

### `R4` — Consolidate Help content and glossary

Changes:

- split the large Analysis Workflow into the six task-centered categories;
- migrate Quick Start, method, troubleshooting, reporting, and privacy copy;
- generate popover, Help glossary, inline definitions, and search aliases from
  canonical topic/term sources;
- move developer/reviewer tables and filenames to explicit technical detail or
  companion artifacts; and
- remove active external-engine/posterior calls to action.

Gate:

- no default Help path contains unapproved internal terms or external workflow
  actions;
- every visible heading and body is localized;
- Help content never exceeds the governing EvidenceRecord/ClaimBoundary; and
- only one glossary source remains.

### `R5` — Separate research detail from technical audit

Changes:

- split substantive detailed tables from technical contract columns;
- revise human-readable reports to describe method, settings, evidence,
  uncertainty, and limitations without class/function names;
- retain exact IDs and schemas in machine-readable sidecars and optional
  reproducibility appendices; and
- label privacy manifests and repeat-analysis scripts by user purpose first.

Gate:

- machine contracts remain byte/schema compatible unless separately versioned;
- public reports contain no unexplained raw IDs or implementation function
  names; and
- privacy allowlists cover every user-facing and technical artifact separately.

### `R6` — Remove compatibility sources

Changes:

- remove `_HELP_POPOVER_LIBRARY`, parallel glossary data, English-title Help
  routing, translated substring routing, and unused unsafe locale bodies only
  after their production reference count reaches zero;
- replace source-inspection tests with contract/behavior tests where possible;
  and
- document retained legacy routes explicitly or remove them under the product
  boundary.

Gate:

- full bilingual, accessibility, privacy, standalone-boundary, and regression
  suites pass;
- the topic/link/target graph has zero orphan or unreachable core nodes; and
- a clean installation follows every advertised Help path.

## Test changes

### New tests

| Test file | Required coverage |
|---|---|
| `tests/test_terminology_contract.py` | Unique concept IDs, layer exposure, locale keys, aliases, prohibited claims, and technical-only terms |
| `tests/test_help_contract.py` | Topic/link/target uniqueness, required fields, valid sections, zero orphans, references, and ClaimBoundary coverage |
| `tests/test_help_navigation.py` | Open, topic change, locale change, return, stale invalidation, sample/real isolation, and no scientific-state mutation |
| `tests/test_user_problems.py` | Stable problem mapping, no raw exception interpolation, safe support references, and exact action targets |
| `tests/test_help_surface_copy.py` | Standard-key allowlist/denylist, active external-action absence, bounded thresholds, and matched bilingual meaning fixtures |
| `tests/test_help_app_lifecycle.py` | Streamlit AppTest journeys from no data, pre-run, blocked preflight, first failure, `HOLD`, `CAUTION`, stale, and post-fit |

### Existing tests to change

| Existing file | Migration |
|---|---|
| `tests/test_help_popover_library.py` | Test registered projections and detailed links. Replace the current unknown-key no-op expectation with visible fallback plus contract failure. |
| `tests/test_guided_essential_ui.py` | Replace Help title/source-string assertions with stable route and rendered behavior assertions. Retain zero-eager-computation checks. |
| `tests/test_i18n_parity.py` | Require every registered key and placeholder in both locales; add boundary-pair fixtures rather than key-count parity alone. |
| `tests/test_privacy_safety_guards.py` | Inject sentinel identifiers/paths/exception text and verify that errors, Help URLs/state, and support output do not disclose them. |
| `tests/test_standalone_core_boundary.py` | Verify new core modules remain Streamlit-free and active Help actions contain no external execution route. |
| `tests/test_respectful_language.py` | Continue task-centered wording enforcement across new Help content and design documents. |
| `tests/test_app_smoke.py` and end-to-end scenarios | Verify persistent Help availability and no regression to analysis execution. |

The following current expectations are explicitly reversed rather than
carried forward as compatibility requirements:

- unknown popover keys being allowed to no-op silently;
- English Help titles acting as route identity;
- translated prose or keyword substrings selecting a Guide section;
- raw ZIP/CSV/JSON filenames being required in Standard Help; and
- a missing critical Help control causing an end-to-end test to skip rather
  than fail.

Source-introspection assertions may remain temporarily to characterize legacy
code, but target behavior is tested through pure registries/reducers and
rendered AppTest journeys. Streamlit AppTest does not establish returned DOM
focus, scroll placement, CSP execution, zoom behavior, or screen-reader
announcement; those claims require the separate browser acceptance pass.
The normative matrix, evidence fields, and pass/fail rules are maintained in
[`help_browser_acceptance.md`](help_browser_acceptance.md).

### AppTest scenario matrix

Each scenario runs in English and Japanese where rendering is language
dependent:

1. no data -> open Help -> select data-format topic -> return to input;
2. data loaded -> open column-mapping Help -> return with mappings unchanged;
3. blocked preflight -> open resource Help -> return to the blocked warning;
4. injected parse/estimation failure -> open troubleshooting -> return with
   data/settings preserved and raw exception absent;
5. `HOLD` EvidenceIssue -> open exact boundary topic -> return to the same
   issue;
6. post-fit contextual popover -> open its exact section -> return to the same
   diagnostic panel and target-owned view state without refit (the first
   active fixture is fit scatter; residual PCA remains a later slice);
7. language switch inside Help -> same topic/origin IDs and translated content;
8. unknown topic -> visible fallback and Help-home action;
9. stale AnalysisID or sample/real mismatch -> static Help only; and
10. Standard -> Help -> Detailed -> return -> same evidence/claim state and no
    newly requested diagnostic.

Browser acceptance additionally covers keyboard-only operation, focus on the
exact returned heading, viewport placement and visible outline, subsequent Tab
order, heading/reading order, one heading announcement without a duplicate
live alert, CSP behavior, 320-pixel width, 200--400% zoom, reduced motion,
non-colour meaning, and touch-target size.

### Copy and privacy gates

The copy scanner operates on the manifest of reachable Standard/Help strings
and hard-coded Streamlit calls, not on machine schemas. This prevents false
failures when a valid JSON export contains `AnalysisID` while still catching a
raw ID in ordinary Help.

The Standard denylist initially covers:

- class and enum names such as `EvidenceRecord`, `ComputationState`,
  `StabilityState`, and `ReasonCode`;
- backend functions such as `compute_pca_bundle` and
  `compute_pca_overall`;
- `payload`, `schema`, `sidecar`, `fixture`, `pipeline`, cache/session-state
  implementation wording, and raw manifest/runner names;
- developer environment variables and exception/stack wording; and
- active instructions to external estimation, posterior, or cross-engine
  routes.

Allowlisted method terms remain visible only with their required definition
and boundary. Exact technical names remain available in registered technical
details and machine-readable artifacts.

## Review and merge boundaries

Implementation is split into reviewable changes:

1. pure contracts and registry tests;
2. app-shell Help route and safe user-problem rendering;
3. contextual Help/return migration;
4. topic content, terminology, and glossary consolidation;
5. report/download technical-layer separation; and
6. legacy removal and browser accessibility completion.

Statistical computation changes, EvidenceRecord/schema changes, visible copy,
navigation changes, and legacy deletion do not share one review unless one is
strictly required for the other. Each slice records the old and new behavior,
tests the no-state-mutation rule, and can be reverted without deleting user
data or changing an existing analysis identity.

## Rollback and safe-fallback rules

- `R0` and `R1` can be reverted without a data or schema migration because
  they add characterization and pure shadow contracts only.
- From `R2` onward, a failed topic renderer falls back to registered generic
  Help plus a return action. It never falls back to a raw exception panel,
  environment-variable bypass, silent no-op, or external-engine action.
- Error redaction and privacy boundaries are monotonic safety changes. They
  are not disabled by a Help-presentation rollback.
- If a temporary build-internal migration mode is required, its only allowed
  states are `shadow`, `active`, and `safe_fallback`. It is not shown to users,
  cannot change analysis state, and is removed with the compatibility wrappers.
- Machine-readable analysis/evidence schemas remain unchanged during these
  slices. Any later schema change requires its own version and migration plan.
- Each merged slice records the last safe tag/commit and its registry/content
  versions. Recovery restores a tested presentation slice, not an unreviewed
  duplicate content store.

## Definition of done

- Static Help is available from every lifecycle state and requires no fitted
  result.
- Every visible Help control has a registered bilingual topic, exact section,
  return target, and focus target.
- No Help route depends on a translated label, English title, or substring.
- Guide/Standard surfaces contain no unapproved implementation term, raw
  exception, developer environment variable, or unexplained ID.
- Method terms retain definitions, applicability, uncertainty, and claim
  boundaries rather than being removed for simplicity.
- Active core Help contains no external-engine, posterior-ingestion, or
  cross-engine-package action.
- Help open/search/locale/return operations make zero estimator calls, request
  zero new diagnostics, and change no scientific or learning state.
- Stale and sample/real context cannot appear as current-analysis Help.
- Human-readable reports and tables separate substantive research detail from
  technical contract columns; machine contracts remain reproducible.
- The glossary has one authoritative source.
- Unknown topics are visible, recoverable, and test-failing rather than silent.
- Bilingual AppTest, accessibility/browser, privacy, standalone-boundary, and
  full regression gates pass before legacy Help sources are removed.
