# Interpretation-copy safety audit

- Status: G0 initial findings; exhaustive manifest and remediation pending
- Last updated: 2026-07-24
- Governing roadmap: [`../ROADMAP.md`](../ROADMAP.md)
- Semantic contract:
  [`study_context_guidance_contract.md`](study_context_guidance_contract.md)
- Remediation specification:
  [`help_and_terminology_remediation_plan.md`](help_and_terminology_remediation_plan.md)

## Purpose

Every explanation that helps a user read a result is part of the scientific
evidence surface. A cautious EvidenceRecord cannot protect users if another
popover, success message, chart guide, tutorial paragraph, or export makes a
stronger claim.

This audit establishes the G0 inventory and release gate for:

- onboarding, Tutorial, Quick Start, and contextual Help;
- status cards, warnings, success messages, metrics, and captions;
- chart guides, popovers, table notes, and accessibility alternatives;
- Japanese and English locale content;
- hard-coded Streamlit strings;
- Analysis Brief, report, manuscript, and other export narratives; and
- legacy external-product language that conflicts with the standalone Python
  product boundary.

This document records findings only. It does not silently rewrite locale or
application code. Stable locale key/function names are the authoritative
locations; line numbers describe the 2026-07-24 working tree and may move.

## Risk classes

| Class | Meaning | Gate |
|---|---|---|
| `STOP` | Directly contradicts product scope/evidence, overstates a scientific or high-stakes conclusion, or can drive harmful action | Must be removed, replaced, or made unreachable before a new guided route exposes the surface |
| `REWRITE` | Uses a universal threshold, all-clear label, deterministic remedy, or unbounded generalization | Must be rewritten and scenario-tested before the affected surface is considered complete |
| `CONTEXT` | Core description may be retained but needs applicability, uncertainty, provenance, or decision boundary | Must receive a stable ClaimBoundaryID and bilingual review |
| `RETAIN` | Already bounded and consistent with the evidence contract | Keep under regression test; do not weaken during consolidation |

An exemption requires an owner, a documented non-core/legacy route, an explicit
boundary shown beside the text, and a test that the core journey cannot present
it as current guidance.

## Required wording contract

All replacement copy follows these rules.

1. Lead with scope: “In this dataset, fitted model, and configured screen …”
   where the boundary materially affects meaning.
2. Say what was detected, not what is universally true. Prefer “no
   high-priority signal was detected under this screen” to “good,” “normal,”
   “acceptable,” “clean,” or “no problem.”
3. Distinguish `AVAILABLE`, `CAUTION`, `HOLD`, and `NOT_ASSESSABLE` from
   `STABLE`, `CONDITIONALLY_STABLE`, `SENSITIVE`, and `NOT_ASSESSED`.
4. Treat numeric cutoffs as named screening conventions or study-specific
   decision rules. Show provenance and uncertainty; never imply a universal
   validity threshold.
5. Distinguish “not detected” from “absent,” “failed to reject” from
   “equivalent,” numerical convergence from model support, and interval overlap
   from a formal comparison.
6. Residual PCA and DIMTEST may identify signals; neither proves
   unidimensionality. An analytic rubric may require a multidimensional claim
   boundary regardless of a quiet residual screen.
7. Reliability/separation may describe distinguishability under the model; it
   does not establish agreement, interchangeability, accuracy, fairness, or
   rater quality.
8. Differential interaction is a screen. Avoid using “bias” as a verdict and
   never convert a statistical flag into automatic exclusion, retraining,
   punishment, category collapse, or fairness conclusions.
9. Dropping a missing row from the likelihood does not establish planned
   absence or an ignorable mechanism. State denominator, reason provenance,
   and what remains unassessed.
10. Logit estimates are model-based measures conditional on data, design,
    identification, fit, and linkage. Do not promise unconditional
    interval-level status or cross-facet/study comparability.
11. High-stakes use requires a declared StudyContext and evidence beyond a
    generic coefficient band. A UI acknowledgement cannot authorize it.
12. The core journey uses standalone Python language. TAM, R, Julia, FACETS,
    external posterior ingestion, Stan handoff, and cross-engine ZIP routes are
    neither prerequisites nor current product promises.
13. Guide and Standard Help use audience-facing terms. Raw enum/class names,
    schema/payload/cache/fixture language, backend function names, environment
    variables, exception types, manifests, and runner names belong only in an
    explicitly requested technical/support layer.
14. “See Help” is not a connection contract. Contextual Help provides an
    actionable link to one registered topic/section and a return path; a
    missing topic cannot silently remove the control.

## Open G0 findings

The initial manual pass contains 38 open findings: 17 `STOP`, 16 `REWRITE`,
and 5 `CONTEXT`. This is not the exhaustive manifest; automated discovery may
add locations or consolidate repeated instances under the same stable CopyID.

### Product boundary and information architecture

| ID | Risk | Surface and stable location | Finding | Required action |
|---|---|---|---|---|
| COPY-001 | `STOP` | [`locales/en.json`](../locales/en.json#L11), keys `sidebar.app_mode_facets`, `app.title`, `sidebar_advanced.run_button`; [`streamlit_app.py`](../streamlit_app.py#L68370) | The default product is named “FACETS-mode,” implying another product defines the runtime or identity | Rename the core route and run action to standalone MFRM language; isolate any historical method comparison as non-core documentation |
| COPY-002 | `STOP` | [`locales/en.json`](../locales/en.json#L1668) and [`locales/ja.json`](../locales/ja.json#L1668), key `help.quick_start_body` | Quick Start advertises Posterior Viewer, generated Stan models, external posterior ingestion, and FACETS-mode execution, all outside the roadmap's core journey | Replace with the canonical seven-stage Python workflow and five-node optional guide; do not link to external handoffs from the guide |
| COPY-003 | `REWRITE` | Same `help.quick_start_body` keys and `help.analysis_workflow_body` at line 1669 | Dense feature/version catalog, universal sample-size minima, “2+ facets,” color-led readiness, and default-model advice are presented before research purpose/design | Split operation from interpretation; route through StudyContext and RQ-to-design preflight; state conditional screening assumptions and give non-colour status text |
| COPY-004 | `STOP` | [`streamlit_app.py`](../streamlit_app.py#L66942), hard-coded `forest_measures.watch`; related Posterior help at lines 66991 onward | Contextual guidance sends users to Posterior Viewer, an excluded external-result route, and the whole library bypasses locale/evidence contracts | Remove external handoff from core guidance; migrate live topics to stable concept/boundary IDs and locale keys |

### Dimensionality, fit, categories, and measurement claims

| ID | Risk | Surface and stable location | Finding | Required action |
|---|---|---|---|---|
| COPY-005 | `STOP` | [`locales/ja.json`](../locales/ja.json#L776), key `dimensionality.interpretation_main`; compare the already-bounded English key at [`locales/en.json`](../locales/en.json#L776) | Japanese copy says fixed residual-PCA cutoffs and an elbow “support unidimensionality,” creating both overclaim and bilingual semantic mismatch | Bring Japanese meaning to parity with the bounded English residual-structure screen; add analytic-rubric and local-dependence boundary |
| COPY-006 | `REWRITE` | [`locales/en.json`](../locales/en.json#L2113) and Japanese counterpart, key `tutorial.section_scores` | Scale scores are declared interval-level and comparable across facets without model, linkage, fit, identification, uncertainty, or target-domain conditions | Describe them as model-based logit estimates under the fitted specification; bind any comparison claim to design/linkage evidence |
| COPY-007 | `STOP` | [`locales/en.json`](../locales/en.json#L2203) and Japanese counterpart, key `fit_details.guide_body`; summary keys at lines 2208--2211 and 2431 | “Perfect fit,” “acceptable,” “generally acceptable,” and “all functioning within acceptable limits” convert screens into universal judgments and can label raters | Replace with deviation descriptions, uncertainty/applicability, configured screen provenance, and non-punitive review actions |
| COPY-008 | `REWRITE` | [`streamlit_app.py`](../streamlit_app.py#L45493), hard-coded all-acceptable success; calibration plot reference line named “Perfect fit” at line 52833; chart guide at line 66908 | Hard-coded strings bypass locale parity and repeat an all-clear decision rule | Use neutral terms such as “identity/reference line” and EvidenceIssue-derived status; add both locales and a text alternative |
| COPY-009 | `REWRITE` | [`locales/en.json`](../locales/en.json#L2177), keys under `categories_steps`; [`streamlit_app.py`](../streamlit_app.py#L67038) | Counts, 0.5--1.5, 1% use, ordered thresholds, and collapse are framed as necessary quality rules or mechanical remedies | Label thresholds as configurable screens; require rubric/content review and sensitivity before any category-change suggestion; say what re-estimation does not establish |
| COPY-010 | `CONTEXT` | [`streamlit_app.py`](../streamlit_app.py#L66958), `qq_residuals`; lines 66942--66955, `forest_measures`; lines 67023--67035, `coverage_heatmap` | “Points near the line -> model fits,” non-overlapping CIs -> significant difference, “mostly blue -> well-covered,” and a 40% rule collapse limited visual cues into global conclusions | State the exact visual feature and next diagnostic; avoid formal-test claims from CI overlap and design sufficiency from colour/density alone |

### Raters, interaction/fairness, weights, and high-stakes use

| ID | Risk | Surface and stable location | Finding | Required action |
|---|---|---|---|---|
| COPY-011 | `STOP` | [`streamlit_app.py`](../streamlit_app.py#L45487), hard-coded reliability success | Low rater reliability is called “the ideal outcome” and raters “interchangeable,” although reliability does not establish agreement or accuracy | Report limited severity differentiation under the current model; require agreement/benchmark/process evidence for broader claims |
| COPY-012 | `STOP` | [`locales/en.json`](../locales/en.json#L2121), `bias_interaction.intro_caption`, no-significance/significance/chi-square messages at lines 2126--2168; hard-coded chart guide at [`streamlit_app.py`](../streamlit_app.py#L66925) | Statistical significance is labelled “significant bias,” while a nonsignificant result becomes “no significant systematic bias,” inviting fairness and rater-quality conclusions | Default to “differential interaction screen”; include multiplicity, uncertainty, sparse overlap, reference group, content review, and “not a fairness verdict” boundary |
| COPY-013 | `REWRITE` | [`locales/en.json`](../locales/en.json#L2114) and Japanese counterpart, key `tutorial.section_weight` | Weight is explained as arbitrary “importance,” including emphasizing selected rater-task combinations, without defining the likelihood, estimand, design weight, or sensitivity consequence | Document the supported weight semantics exactly; block unsupported use cases and require a sensitivity/estimand note |
| COPY-014 | `STOP` | [`locales/en.json`](../locales/en.json#L2022) and Japanese counterpart, key `simulation_validation.g_d_study_d_study_caption` | Generic G/Phi >= .8 is labelled suitable for high-stakes use and >= .7 for routine reporting | Present coefficients descriptively under the declared design; require StudyContext, uncertainty, decision loss, classification, and external validity evidence before use claims |
| COPY-015 | `REWRITE` | [`locales/en.json`](../locales/en.json#L1752) and Japanese counterpart, key `help.reporting_body` | Fixed thresholds define “good differentiation,” rater interchangeability, fit, bias, and required reporting; generated manuscript prose appears authoritative | Generate a draft reporting matrix from actual EvidenceIDs, applicability, sensitivity, and unresolved limitations; label it `draft_for_review` |

### Missingness, validation, and decision provenance

| ID | Risk | Surface and stable location | Finding | Required action |
|---|---|---|---|---|
| COPY-016 | `STOP` | [`locales/en.json`](../locales/en.json#L2117), key `tutorial.section_missing_values`; related workflow/troubleshooting copy | Missing values are treated as structurally absent/FACETS-conventional, collapsing planned non-assignment and unexpected missingness | Use explicit reason provenance and denominators; say only that specified rows do not enter the likelihood and preserve what remains unassessed |
| COPY-017 | `REWRITE` | [`locales/en.json`](../locales/en.json#L2004), key `simulation_validation.parameter_recovery_convergence_note` | Failed replications are said to be excluded so results are “not distorted,” which can hide selection from the requested simulation denominator | Retain requested/completed/successful counts, report numerical failures separately, and bound conditional recovery metrics to successful fits |
| COPY-018 | `CONTEXT` | English and Japanese analysis-workflow/sample-size/rating-scale guides, especially keys at lines 1669, 1672, and 1754 | Thresholds and remedies are often stated as universal minima, optimal ranges, or automatic category/design changes | Convert them to named literature conventions or configured screens, add applicability and uncertainty, and link actions to what they repair and do not prove |
| COPY-019 | `REWRITE` | Multiple success labels, including `fit_details.summary_all_acceptable_template`, `facet_dashboard.all_clean_success_template`, and visually green states | Positive colour and “clean/all acceptable” language can imply an overall validity/fairness all-clear even when other evidence is missing or in tension | Render the current EvidenceIssue plus aggregate `HOLD`/`CAUTION`/unassessed/tension counts; require text status independent of colour |

### Internal-language leakage and Help connectivity

| ID | Risk | Surface and stable location | Finding | Required action |
|---|---|---|---|---|
| COPY-020 | `REWRITE` | [`streamlit_app.py`](../streamlit_app.py#L53276), visible `_standalone_quick_start_markdown` and `_standalone_analysis_workflow_markdown` | Standard Help instructs users with `EvidenceRecords`, `app-engine runner`, `ComputationState`, `StabilityState`, `ReasonCode`, and `AnalysisID`; Japanese additionally mixes many untranslated implementation terms | Use the terminology registry and plain task/status/action labels; place exact IDs/enums and reproducibility internals in optional technical details |
| COPY-021 | `STOP` | [`streamlit_app.py`](../streamlit_app.py#L67147), `render_help_popover`; locale captions at [`locales/en.json`](../locales/en.json#L1665) | The popover merely says to see Help rather than opening the exact topic, and an unknown topic silently returns without a control or fallback | Introduce `HelpTopic`/`HelpLink`; provide one-action exact routing, return/focus restoration, registry tests, and a visible safe fallback |
| COPY-022 | `REWRITE` | [`streamlit_app.py`](../streamlit_app.py#L53401), Analysis Workflow Help | A first-run Help panel expands into many developer/reviewer tables, raw filenames and a reference-coverage CSV; “Visual claim guardrails” and captions are hard-coded English | Keep the first-run path task-oriented; move audit inventories/filenames to Detailed technical or reviewer artifacts; localize every visible heading |
| COPY-023 | `STOP` | [`streamlit_app.py`](../streamlit_app.py#L259), `render_exception_details` | Hosted UI still exposes raw exception type/message and tells users to set `MFRM_SHOW_TECHNICAL_ERRORS=1`, leaking implementation detail and potentially sensitive error content | Show a localized problem statement, reversible next step, and privacy-safe support reference; keep exception text/environment instructions in local developer logging only |
| COPY-024 | `STOP` | [`locales/en.json`](../locales/en.json#L798), dimensionality reason/exception keys; Japanese counterparts | Users can see `result['config']`, `StdResidual`, `compute_pca_bundle`, `compute_pca_overall`, NumPy eigendecomposition, error types, and raw exception messages | Map backend failures to stable user reason codes and plain repairs; make technical details opt-in and redacted |
| COPY-025 | `STOP` | [`locales/en.json`](../locales/en.json#L986), `resource_preflight.blocked_error`; download caption at line 2357 | UI advises an environment-variable safety override and explains hidden rerun/ZIP/Stan package generation mechanics instead of the user consequence | Do not advertise bypass flags in hosted UI; explain size/safety limits and supported next actions. Describe download cost/output, not internal render mechanics |
| COPY-026 | `CONTEXT` | [`streamlit_app.py`](../streamlit_app.py#L53343), `show_help_section` | Help routing uses English presentation phrases as IDs, has no versioned topic/claim-boundary/return registry, and cannot validate orphaned contextual links | Replace with language-independent `help_topic_id` and `section_id`, a current-analysis binding, locale aliases, review metadata, and reachability tests |
| COPY-027 | `REWRITE` | [`locales/en.json`](../locales/en.json#L691), reproducibility/claim-trace copy; first-run tutorial at [`streamlit_app.py`](../streamlit_app.py#L22280) | Guided/research guidance foregrounds fingerprints, JSON filenames, `AnalysisID`, `EvidenceRecord`, engine and reproducibility-contract wording | Explain the user purpose first (“save what is needed to reproduce this analysis”); move exact filenames/IDs to an optional reproducibility detail |
| COPY-028 | `CONTEXT` | Locale Help around [`locales/en.json`](../locales/en.json#L829), lines 1978 and 2062--2084; run-history UI | Cache, session state, bundle, hidden columns, raw complexity notation, dtype/cardinality/OOM and similar implementation terms are used where users need only cost, persistence, or repair consequences | Translate to browser-session duration, expected cost, memory limit, readable-data requirement, and next action; keep exact mechanics in administrator/developer Help |
| COPY-029 | `REWRITE` | [`streamlit_app.py`](../streamlit_app.py#L23521), compact/full table renderer; sensitivity table path around lines 35646 and 48437 | “Full table” can reveal `SensitivityPlanID`, `DecisionRuleJSON`, `ContractSchemaVersion`, and other machine-contract columns as ordinary research detail | Separate substantive full results from an explicit technical-audit table; preserve machine columns in CSV/JSON contracts |
| COPY-030 | `REWRITE` | [`locales/en.json`](../locales/en.json#L1754) and Japanese counterpart, Troubleshooting | User repair guidance exposes pattern matcher, `_skip_reason`, config, dtype, cardinality, OOM, and fallback vocabulary | Render the actual user condition and repair in plain language; retain `maxit`/`reltol` only in advanced method troubleshooting |
| COPY-031 | `REWRITE` | Report/action prose around [`streamlit_app.py`](../streamlit_app.py#L38796), lines 39691 and 42144 | Manuscript-facing prose asks users to report/archive `AnalysisIdentity`, `EvidenceRecords`, raw IDs, runner and function names as if they were substantive results | Keep app version, model, estimator, settings and limitations in Methods; move IDs/file inventory to a reproducibility appendix or technical sidecar |
| COPY-032 | `CONTEXT` | [`locales/en.json`](../locales/en.json#L2360), public-export notice | `export_privacy_manifest` is named without explaining what it contains or where to inspect it | Lead with “included/excluded files and disclosure review,” then show the filename as secondary technical detail |
| COPY-033 | `STOP` | Main flow around [`streamlit_app.py`](../streamlit_app.py#L34425) and lines 68424--68433; Help renderers at lines 27499 and 35206 | Full Help is reachable only inside post-estimation Guided/Full result routes, so users cannot reach Troubleshooting when there is no data, during mapping, or after the first failed run | Make static Help persistent from landing through all lifecycle/error states; current-analysis Help remains conditional |
| COPY-034 | `REWRITE` | [`streamlit_app.py`](../streamlit_app.py#L53343), `force_full` Help; description key at [`locales/en.json`](../locales/en.json#L771) | “Show all Help” promises rating-scale, model-capability, and public-beta content, but the implementation adds only Rating Scale Guide while other keys/bodies are orphaned | Replace display-density coupling with topic-level method/reference depth; remove or explicitly isolate orphaned content and make the control description exact |
| COPY-035 | `STOP` | Guided routing around [`streamlit_app.py`](../streamlit_app.py#L26209); Help panel IDs around line 53355 | Some routes depend on translated display labels and English/Japanese substring matching; Help panel IDs reuse English titles | Route only by stable language-independent topic/section/target IDs; add rename, alias, and locale regression tests |
| COPY-036 | `REWRITE` | Glossary sources around [`streamlit_app.py`](../streamlit_app.py#L2696) and line 53518 | Inline glossary and the large Help glossary are parallel sources despite comments implying one source, allowing definition and boundary drift | Generate inline, Help, search aliases, and export glossary from one reviewed terminology registry |
| COPY-037 | `STOP` | Locale Help bodies such as `help.interpretation_guide_body`, `reporting_body`, and `troubleshooting_body` | Instructions hard-code Full-specific tab paths that do not exist under the guide/Standard projection, producing practical dead ends rather than executable links | Store related target IDs and resolve the correct localized screen/action for the active view, with a return path |
| COPY-038 | `STOP` | Active contextual topics around [`streamlit_app.py`](../streamlit_app.py#L66820) and line 66953 | MML Help still recommends TAM/ConQuest/mirt comparison and measure Help sends users to Posterior Viewer, contradicting the standalone core boundary | Remove active external-engine/posterior calls to action; retain any methodological comparison only as clearly bounded reference text |

## Safe patterns to retain

The inventory also contains copy that already demonstrates the target style
and should become regression fixtures.

| Location | Why retain |
|---|---|
| [`locales/en.json`](../locales/en.json#L776), `dimensionality.interpretation_main` | Explicitly calls residual PCA an exploratory screen, rejects standalone proof, and gives bounded wording |
| [`locales/en.json`](../locales/en.json#L2128), `bias_interaction.dff_caption` | Calls the result a screening aid and rejects standalone proof of differential functioning |
| [`streamlit_app.py`](../streamlit_app.py#L66807), `classical_dif.watch` | Says a flag warrants review rather than automatic removal or a fairness verdict |
| Roadmap Evidence vocabulary and claim-boundary section | Separates computation from conclusion stability and lists prohibited inferences |

Retention does not exempt these surfaces from bilingual, applicability, and
accessibility review.

## Inventory manifest

G0 produces a machine-readable manifest with one row per interpretive surface.
The proposed fields are:

| Field | Purpose |
|---|---|
| `copy_id` | Stable audit identity such as `COPY-012` |
| `surface_type` | Tutorial, Help, popover, warning, success, caption, chart guide, table note, export, or manuscript |
| `source_location` | Locale key or stable function/component ID; line number is secondary |
| `locale` | `en`, `ja`, shared/hard-coded, or user-authored |
| `concept_id` | Stable statistical/research concept |
| `claim_boundary_id` | Governing bounded-claim rule |
| `audience_layer` | Guide/Standard, Detailed/method, or technical/support exposure |
| `help_topic_id` | Exact registered Help destination when the surface links onward |
| `return_target_id` | Calling card/control restored after Help |
| `evidence_requirements` | Required EvidenceIDs/types and applicability states |
| `risk_class` | `STOP`, `REWRITE`, `CONTEXT`, or `RETAIN` |
| `finding_status` | `OPEN`, `REWRITE_PLANNED`, `IMPLEMENTED`, `VERIFIED`, or `EXEMPT` |
| `owner` | Responsible module/reviewer, not an unowned TODO |
| `replacement_key` | New locale/template key where applicable |
| `test_ids` | Unit, AppTest, browser, and human-review evidence |

An inventory based only on locale keys is incomplete because current chart
guides and status messages include hard-coded English. Static discovery MUST
scan both locales, `streamlit_app.py`, `mfrm_app/` report builders, and document
templates.

## Automated gates

Create a focused contract suite, provisionally
`tests/test_interpretation_copy_contract.py`, with the following checks.

### Static checks

- Every rendered interpretive locale key and hard-coded surface is present in
  the manifest or an approved non-interpretive allowlist.
- New instances of high-risk phrases such as “perfect fit,” “all acceptable,”
  “no bias,” “proves unidimensionality,” or high-stakes coefficient bands fail
  in governed surfaces unless a reviewed exemption identifies the quotation or
  negation context.
- Core guide/navigation strings do not advertise external engines, external
  posterior ingestion, Stan handoff, or cross-engine interchange.
- Stable concept/ClaimBoundaryIDs exist in both locales; translated display
  text is never used for routing.
- Hard-coded English interpretation does not enter a bilingual journey.
- Guide/Standard Help contains no unapproved implementation enum, class,
  schema, backend function, environment-variable, raw exception, fixture,
  manifest, or runner wording.
- Every visible contextual Help control resolves to a registered topic and
  section; no unknown topic silently no-ops.
- Help topic IDs are language-independent and every rendered topic has a valid
  source/return target and ClaimBoundaryID.

Simple banned-word tests are insufficient on their own. “Does not prove
unidimensionality” must not fail because it contains “prove,” while “supports
unidimensionality” may still be too strong in context. The manifest and
scenario assertions remain authoritative.

### Scenario checks

- Quiet residual PCA with missing local-dependence evidence never produces an
  unidimensionality all-clear.
- Low rater reliability never produces agreement, interchangeability,
  accuracy, or quality language.
- No flagged differential interaction never produces a no-bias/fairness claim.
- Ordered categories never certify scale validity; disordered thresholds never
  trigger an automatic collapse instruction.
- A converged run with other `HOLD` evidence never produces “analysis good” or
  “all acceptable.”
- Planned missingness and unknown missingness produce different wording and
  claim consequences.
- Failed simulation replications remain visible in requested/completed/success
  counts.
- A high-stakes StudyContext does not unlock use from a generic reliability,
  G/Phi, fit, or acknowledgement threshold.
- Conflicting evidence renders the tension and aggregate issue counts in both
  Standard and Detailed views.
- Analysis Brief and JSON sidecar express the same ClaimIDs, statuses,
  boundaries, and unresolved limitations in both locales.
- Opening contextual Help reaches the requested section in one action,
  preserves the originating analysis/card, and returns focus without changing
  scientific or learning state.
- Stale, mismatched, or sample-versus-real Help context displays only the
  general method explanation and never presents a previous result as current.
- A production error reveals no raw exception, backend symbol, environment
  override, uploaded value, or identifier; the privacy-safe support reference
  still permits diagnosis.

## Human review gate

Automated parity cannot establish equivalent meaning. Before G2 is complete,
the affected scenarios receive independent review from:

1. an MFRM/measurement reviewer;
2. an SLA or language-assessment researcher;
3. a bilingual Japanese/English reviewer;
4. a plain-language/accessibility reviewer; and
5. for data-governance/high-stakes content, an appropriate ethics or data
   governance reviewer.

Each reviewer marks:

- what the statement supports;
- what it could be misread to support;
- the required context/evidence;
- whether it labels a person, rater, group, or institution unfairly;
- whether the recommended action is reversible and proportionate; and
- whether Japanese and English preserve the same boundary.

Disagreement is recorded as an unresolved evidence/copy issue, not erased by
majority vote.

## Remediation sequence

1. **G0 inventory:** discover all surfaces, assign stable IDs, classify every
   finding, and freeze characterization screenshots/tests. Do not redesign the
   UI in this step.
2. **G1 STOP remediation:** remove stale external-product promises and replace
   claims that can certify validity, fairness, rater quality, category changes,
   or high-stakes use. Remove raw exception/environment/backend leakage, add
   evidence-bound renderers, and add locale parity tests.
3. **G1 REWRITE/CONTEXT remediation:** convert fixed rules to named screens,
   add applicability and uncertainty, route actions through EvidenceIssue,
   implement the audience terminology registry, and register all
   `HelpTopic`/`HelpLink` targets.
4. **G2 guide gate:** expose the optional sample route only after every string
   used by that journey is `VERIFIED` in both locales and every Help action can
   open the exact topic and return to its source.
5. **G3/G4 projection and artifact gate:** verify that Standard, Detailed,
   Analysis Brief, JSON, report, and manuscript surfaces cannot disagree for
   the same IDs.
6. **G6 consolidation:** remove only proven-unreachable duplicates, keep safe
   copy regression fixtures, and close or explicitly exempt every manifest row.

## Completion criteria

The audit is complete only when:

- all rendered interpretation surfaces are in the manifest;
- no `STOP` finding is reachable in the core journey;
- every `REWRITE` and `CONTEXT` item used by a shipped surface is `VERIFIED`;
- Japanese and English pass semantic scenario review, not just key parity;
- Current focus and exported prose use the same EvidenceIssue/ClaimBoundary
  source;
- Guide and Standard Help contain no unexplained implementation language, and
  technical/support terms appear only in their declared audience layer;
- every contextual Help control has a registered bilingual target, safe
  current-context binding, accessible focus transition, and return path;
- person/rater/subgroup and high-stakes scenarios pass the non-punitive and
  use-boundary review;
- automated, AppTest, browser, and human-review evidence is linked from each
  applicable row; and
- accepted exemptions are visibly separated from the standalone core and have
  an owner and removal/review date.

Until then, existing copy is a known migration risk and must not be cited as
evidence that the planned guidance architecture has been implemented.
