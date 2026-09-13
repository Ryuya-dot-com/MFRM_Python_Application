# Study context, guidance, and identity contract

- Status: P0 contract; initial G1/G2 guidance slice implemented
- Last updated: 2026-08-04
- Product boundary: standalone Python
- Governing roadmap: [`../ROADMAP.md`](../ROADMAP.md)
- Copy-risk gate: [`interpretation_copy_audit.md`](interpretation_copy_audit.md)
- Remediation specification:
  [`help_and_terminology_remediation_plan.md`](help_and_terminology_remediation_plan.md)

## Purpose

This document freezes the semantic boundary required before restructuring the
Streamlit interface. It defines which records own research context, analysis
settings, evidence, learning progress, interpretation, decisions, and
artifacts. It is intentionally independent of Streamlit widgets and page
layout.

The contract has four goals:

1. prevent a display or tutorial action from changing a scientific result;
2. prevent a completed workflow from being mistaken for a valid claim;
3. preserve why an analysis was undertaken and what it can and cannot support;
4. give the Standard and Detailed views one canonical state to project.

`MUST`, `MUST NOT`, `SHOULD`, and `MAY` are normative terms in this document.
Examples illustrate the contract but do not expand the supported statistical
scope.

Implementation note (2026-08-04): `mfrm_app/guidance.py` now owns the pure
five-node catalog, session lifecycle, skip/exit/resume/restart transitions,
fit binding, formative learning evidence, invalidation, and learning-only
completion. The bilingual Streamlit adapter exposes the complete basic sample
journey. This does not mark G1/G2 complete: canonical `AnalysisSpec`, the
remaining research-context records, exact Help return at every guide node,
`HOLD`/stale/failure journeys, a complete existing-real-fit sample roundtrip,
and browser accessibility acceptance remain gated work.

## Non-negotiable invariants

1. `StudyContext`, `AnalysisSpec`, `AnalysisIdentity`, evidence,
   interpretation, decision, and guidance are separate records.
2. A view, language, Help, or guide event MUST NOT mutate `AnalysisSpec`, a
   fitted AnalysisID, EvidenceRecords, conclusion stability, or claim
   readiness.
3. Streamlit widget state MUST NOT be the authoritative copy of an analysis
   setting.
4. Guide completion MUST NOT certify validity, fairness, research competence,
   or suitability for a high-stakes use.
5. A `HOLD` or required `NOT_ASSESSABLE` state MAY be archived as an honest
   no-conclusion record; it MUST NOT be converted into a supported claim.
6. The same numerical analysis MAY be interpreted differently under a
   different purpose, population, construct, use, or stakes. This changes
   context, interpretation, and decision identities without changing the
   AnalysisID unless the analytic specification also changes.
7. Planned absence and unexpected missingness MUST NOT share an
   indistinguishable default code.
8. Statistical flags for a person, rater, or subgroup are investigation
   prompts, not automated personnel, bias, fairness, or quality verdicts.
9. A human-readable Analysis Brief is the primary user artifact. Its JSON
   sidecar is the machine-readable reproducibility contract.
10. Canonical records MUST contain stable IDs and values, not translated UI
    labels or English-substring routing rules.
11. Every visible contextual Help control MUST resolve to a registered topic
    and return target; critical limits remain inline and unknown topics never
    fail silently.
12. Guide and Standard surfaces MUST use the approved audience vocabulary and
    MUST NOT expose raw exceptions, backend symbols, environment overrides, or
    unexplained implementation IDs.

## Record graph

```text
StudyContextID
├── RatingDesignRecordID
├── ScoringProcessRecordID
└── planned ClaimIDs
          │
DataFingerprint + fitted AnalysisSpecID
          │
          └── AnalysisID
                ├── EvidenceIDs
                ├── EvidenceIssueIDs
                └── InterpretationID
                      └── DecisionID
                            └── ArtifactID

GuidanceState
├── may reference StudyContextID, DataFingerprint, AnalysisID, and EvidenceID
└── is never an input to any scientific identity above

HelpTopic registry <- HelpLink(source target, optional verified context,
                               return target/focus)
└── explains scientific records but is never an input to their identities
```

The record graph is directional. A guidance record may point to scientific
state so it can resume at the correct place; scientific state MUST NOT depend
on whether guidance was opened, skipped, answered, or completed.

## Common record envelope

Every versioned record SHOULD use a small common envelope.

| Field | Requirement |
|---|---|
| `schema_name` | Stable machine name, not a translated label |
| `schema_version` | Explicit semantic schema version |
| `record_id` | Content-derived or lineage-derived stable ID under the record's identity rule |
| `created_at` | ISO 8601 timestamp; excluded from content identity unless chronology is the subject |
| `app_version` | Application version that produced the record |
| `parent_ids` | Typed lineage references where applicable |
| `status` | Record-specific closed vocabulary |
| `extensions` | Versioned namespace; unknown core fields fail closed |

Canonical serialization MUST define key ordering, Unicode normalization,
number representation, missing-value representation, and list ordering before
content hashing. Display strings, timestamps that merely record export time,
and view preferences MUST NOT produce a new scientific identity.

## `StudyContext`

`StudyContext` records the interpretation-and-use argument that exists before
a result is read. It does not claim that the argument is valid.

### Required fields

| Field | Meaning |
|---|---|
| `study_context_id` | Stable ID for the normalized context |
| `purpose` | Closed value such as `METHOD_RESEARCH`, `FORMATIVE`, `SUMMATIVE`, `SELECTION`, `CERTIFICATION`, or `OTHER_DECLARED` |
| `stakes` | `LOW`, `MODERATE`, `HIGH`, or `UNDECLARED` |
| `target_population` | Population for which interpretation is intended, including exclusions or an explicit undeclared value |
| `target_domain` | Performances, tasks, and conditions to which the score is intended to generalize |
| `construct` | Named focal construct in the user's study vocabulary |
| `subconstructs` | Declared subconstructs and their relationship to the focal construct |
| `performance_mode` | For example writing, speaking, or human-scored short response |
| `rubric_structure` | `HOLISTIC`, `ANALYTIC_SINGLE_CONSTRUCT`, `ANALYTIC_SUBCONSTRUCTS`, `COMPOSITE`, or `UNDECLARED` |
| `intended_score` | Reported score/measure, direction, unit, aggregation, and target object |
| `intended_interpretation` | Bounded proposition the score is intended to support |
| `intended_use` | Action or communication for which the interpretation may be used |
| `intended_decision` | Decision rule or explicit `NONE_RESEARCH_ONLY` |
| `planned_claim_ids` | Stable claims whose prerequisites will be evaluated |
| `unsupported_inferences` | Important inferences known to be outside current evidence or product scope |
| `analysis_status` | `PLANNED_CONFIRMATORY`, `PLANNED_EXPLORATORY`, `POST_HOC_EXPLORATORY`, or `UNDECLARED` |

An undeclared required value is represented explicitly. It is not silently
filled from a domain stereotype. Claim-specific applicability decides whether
an undeclared value causes `CAUTION`, `HOLD`, or `NOT_ASSESSABLE`.

Free-text notes MAY exist locally but MUST be separated from stable structured
fields, excluded from telemetry and support bundles, and reviewed before
export because they may contain personal or confidential information.

### Scope boundary

Native MFRM evidence can contribute to scoring and limited generalization
arguments within a declared design. It does not by itself establish content
representation, response processes, external relationships, accessibility,
consequences, fairness, or the validity of a proposed use. The context record
therefore stores both planned claims and unassessed warrants.

## `RatingDesignRecord`

`RatingDesignRecord` describes the planned and observed rating design rather
than reducing design quality to graph connectedness.

### Required structure

| Field | Meaning |
|---|---|
| `rating_design_record_id` | Versioned record identity |
| `object_of_measurement` | Person, performance, product, or another declared target |
| `unit_of_observation` | What one input row represents |
| `performance_unit` | Script, response, turn, task performance, or declared equivalent |
| `facets` | Stable facet IDs with role, relation, population intent, anchor state, and level provenance |
| `relations` | Crossed, nested, partially crossed, confounded, or unknown relationships |
| `planned_assignment` | Intended assignment rule, overlap, common performances, and workload |
| `observed_assignment` | Realized counts, overlap, components, bridges, and deviations from plan |
| `scoring_multiplicity` | Single, double, variable, third rating, or adjudication-only rules |
| `occasion_and_order` | Session, occasion, order, fatigue, and timing variables when relevant |
| `anchor_plan` | Anchor persons/performances/raters/tasks, source, and maintenance rule |
| `rubric_version` | Stable rubric/scoring-guide version |
| `training_cohort` | Rater training/qualification cohort or explicit unavailable value |
| `missingness_plan` | Expected absent cells and permitted reason-code vocabulary |

Each facet declares at least:

- whether it is the measurement object, a nuisance facet, a contextual
  covariate, or record-only metadata;
- whether levels are treated as fixed, random, anchored, or unspecified for
  the intended claim;
- whether the facet is crossed, nested, partially crossed, or confounded with
  other facets; and
- whether inference concerns only observed levels or a broader population.

Few facet levels do not automatically imply stable network inference. The
record and downstream evidence retain observations per level, common ratings,
score range, extremes, leverage, overlap, and the boundary on generalizing
beyond observed levels.

## `ScoringProcessRecord`

`ScoringProcessRecord` preserves operational scoring evidence that model-fit
statistics cannot reconstruct.

| Field | Meaning |
|---|---|
| `scoring_process_record_id` | Versioned record identity |
| `selection` | How raters were recruited or assigned |
| `training` | Materials, duration, calibration examples, and rubric version |
| `qualification` | Criteria, attempt rules, and benchmark source |
| `monitoring` | Ongoing checks, frequency, and intervention policy |
| `double_scoring` | Selection rule and whether ratings are blind/independent |
| `third_rating` | Trigger and selection rule |
| `adjudication` | Procedure, authority, and whether adjudicated scores enter estimation |
| `rescoring` | Trigger, scope, overwrite/retention rule, and provenance |
| `benchmark_process` | Expert/consensus/backrating source needed for any accuracy claim |
| `non_punitive_review` | Declared rule that statistical flags begin contextual review rather than automated personnel action |

Consistency and accuracy are distinct. Without a declared benchmark process,
the application MUST NOT infer rater accuracy from severity, fit, agreement,
or consistency statistics.

## Missingness provenance

Every absent expected rating SHOULD have a reason code where the source data
permits it. The minimum vocabulary is:

| Code | Class | Example handling |
|---|---|---|
| `PLANNED_NOT_ASSIGNED` | Planned design | Excluded from the expected-rating denominator defined for completed assignments |
| `TEST_TAKER_NONRESPONSE` | Unexpected | Retained in missingness summaries; mechanism and score dependence assessed |
| `NOT_REACHED` | Unexpected | Distinguished from intentional omission |
| `TECHNICAL_FAILURE` | Unexpected | Reported by occasion/task and included in sensitivity planning |
| `RATER_SKIP` | Unexpected | Reported by rater/task and investigated |
| `UNRATABLE_RESPONSE` | Unexpected/operational | Reason and rubric rule retained |
| `CRITERION_NOT_APPLICABLE` | Planned or rule-based | Requires an explicit applicability rule |
| `ADMINISTRATIVE_EXCLUSION` | Administrative | Counts before and after exclusion retained |
| `ADJUDICATION_ONLY` | Process | Distinguished from an independent routine rating |
| `UNKNOWN_MISSING` | Unknown | Never recoded as planned absence; may block a claim |

Missingness summaries retain denominators, exclusions, and distributions by
relevant facet, score range, subgroup, and reason. A row dropped from a
likelihood is not thereby ignorable for interpretation.

## `AnalysisSpec` and `AnalysisIdentity`

`AnalysisSpec` is a canonical, serializable object owned by the application
core. It has separate `draft` and `fitted` instances.

### Included in `AnalysisSpec`

- data-role mappings and score/category handling;
- response model and facet structure used by the estimator;
- estimator, identification, constraints, anchors, weights, initialization,
  optimizer, tolerance, and iteration rules;
- uncertainty method and settings;
- requested diagnostics and analysis depth;
- sensitivity plan and comparison rules; and
- all seeds and deterministic options that can change a result.

### Excluded from `AnalysisSpec`

- Standard/Detailed view;
- guide lifecycle and learning answers;
- language, theme, chart display density, expanded panels, and selected Help;
- translated strings; and
- browser/download completion state.

An AnalysisID is derived from the canonical fitted specification, data
fingerprint, code/schema identity, and other computational inputs defined by
the existing evidence contract. Editing a draft does not rewrite the fitted
specification. It marks the current fit stale for claims that require the
draft until the user explicitly fits or restores the fitted settings.

## `EvidenceIssue`

`EvidenceIssue` groups records needed to answer one prioritized question. It
prevents one reassuring EvidenceRecord from hiding contradictory or missing
evidence.

| Field | Meaning |
|---|---|
| `evidence_issue_id` | Stable identity for issue definition plus bound evidence |
| `issue_type` | For example connectivity, convergence, dimensionality boundary, category use, interaction screen, or sensitivity |
| `priority` | Deterministic priority plus reason code |
| `required_evidence_ids` | Evidence required for the current ClaimID/context |
| `supporting_evidence_ids` | Available supporting records |
| `limiting_evidence_ids` | Cautions, holds, unassessable records, and failed variants |
| `tension_state` | `NONE`, `PRESENT`, `UNRESOLVED`, or `NOT_ASSESSABLE` |
| `claim_boundary_ids` | Stable boundaries that control rendered wording |
| `repair_targets` | Named stage/control targets that may change evidence |
| `next_actions` | Ordered, bounded actions; no automatic validity or personnel verdict |

Current focus renders the highest-priority issue and also reports totals for
all current `HOLD`, `CAUTION`, required `NOT_ASSESSABLE`, required
`NOT_ASSESSED`, and tension states.

## Interpretation, claim, and decision records

### Claim record

A ClaimID identifies a semantic claim template, not one translated sentence.
A claim instance stores:

- StudyContextID, AnalysisID, and required EvidenceIDs;
- scope and population boundaries;
- `SUPPORTED_WITHIN_SCOPE`, `UNSUPPORTED`, `WITHHELD`, or `NOT_EVALUATED`;
- conclusion-stability state and unresolved limitations;
- stable reason codes;
- optional user-authored wording stored separately from localized system copy;
  and
- author/reviewer status.

`claim_ready` is derived per ClaimID. It is true only when the context is
sufficient, every required applicability check and computation satisfies the
claim rule, required sensitivity is complete, and no blocking tension remains.
It is never a persisted checkbox.

### `InterpretationRecord`

The interpretation record answers: “What does the current evidence support
about the declared score or construct, within which boundary?” It references
ClaimIDs and explicitly lists unsupported and unassessed warrants.

### `DecisionRecord`

The decision record answers: “What action, communication, or deferral is made
for the declared use, by whom, under what review state?” A decision may be
`DEFER`, `COLLECT_MORE_EVIDENCE`, `REVISE_DESIGN`, `REPORT_WITH_BOUNDARY`, or
another versioned action. It MUST NOT imply that the application made or
authorized a high-stakes decision.

InterpretationID changes when the bounded interpretation changes. DecisionID
changes when the intended use, action, decision rule, review, or interpretation
binding changes. Neither change alone changes AnalysisID.

## `GuidanceState`

`GuidanceState` is session-scoped in the initial implementation.

| Field | Meaning |
|---|---|
| `lifecycle` | `NOT_STARTED`, `ACTIVE`, `SKIPPED`, or `COMPLETED` |
| `route_id` | `sample`, `own_data`, or `without_guide` |
| `active_node_id` | Stable current node |
| `last_valid_node_id` | Last node whose bindings remain valid |
| `node_states` | `LOCKED`, `AVAILABLE`, `ACTIVE`, or `COMPLETE` plus learning evidence |
| `bound_study_context_id` | Optional context binding |
| `bound_data_fingerprint` | Optional data binding |
| `bound_analysis_id` | Optional fitted-result binding |
| `reviewed_evidence_ids` | Evidence viewed; viewing alone does not establish understanding |
| `sample_context_id` | Isolated tutorial context when using the sample |
| `skip_origin_node_id` | Node at which the route was exited |

The landing surface is the `welcome` node. It MUST NOT route to a second screen
that asks the same entry question.

### Independent completion fields

- `workflow_complete`: a route reached its required terminal action or a
  preserved failure/no-conclusion terminal record.
- `learning_complete`: the defined learning evidence for a node or route is
  present. Where a formative check exists, opening the card is insufficient.
- `claim_ready`: derived from scientific records and never written by the
  guidance reducer.

Valid combinations include:

| Workflow | Learning | Claim | Meaning |
|---|---|---|---|
| Complete | Complete | False | Guide finished, but evidence does not support the claim |
| Complete | Incomplete | True | Analysis may be scientifically ready although tutorial learning is unfinished/skipped |
| Complete | Complete | True | Both independent requirements happen to be satisfied |
| Complete | Complete | False with `HOLD` | A no-conclusion analysis record was archived correctly |

## Guide-to-research mapping

| Guide node | Research stage(s) | Learning evidence | Scientific rule |
|---|---|---|---|
| `welcome` | Plan | Route selected | No scientific mutation |
| `data_check` | Plan, Check | Roles/design/readiness reviewed | `HOLD` remains until data/design changes |
| `estimate` | Estimate | Run outcome, including failure, interpreted | Only a current successful fit supplies fitted evidence |
| `evidence_review` | Diagnose, Stress, interpretation substep | Supported and unsupported statements distinguished in the sample | Evidence state is unchanged by the answer |
| `archive` | Decide, Archive | Bounded claim or no-conclusion, limitation, and next action recorded | Artifact creation does not make a claim ready |

The canonical product navigation remains `Plan -> Check -> Estimate ->
Diagnose -> Stress -> Decide -> Archive`. `InterpretationRecord` is a separate
identity inside the transition from Stress to the use decision even while the
UI groups it within Decide.

## Lifecycle events and effects

| Event | Guidance effect | Scientific effect |
|---|---|---|
| `START_SAMPLE` | Activate isolated sample context | Preserve current real-data state; bind documented sample specification separately |
| `START_OWN_DATA` | Route to real-data preflight | None until data/context/spec changes |
| `SKIP` / `EXIT` | Preserve progress and record origin | None |
| `RESUME` | Open first valid incomplete node | None |
| `RESTART` | Reset learning progress only | None |
| `VIEW_CHANGED` | Re-render `standard` or `detailed` projection | No spec change, estimator call, or new diagnostic |
| `LANGUAGE_CHANGED` | Re-render locale keys | No ID or evidence change |
| `HELP_OPENED` / `HELP_RETURNED` | Update topic, source, return, and focus navigation only; does not complete learning | No ID, evidence, claim, setting, or computation change |
| `FORMATIVE_ANSWERED` | Store learning evidence/feedback | No computation, stability, or claim change |
| `STUDY_CONTEXT_CHANGED` | Re-evaluate interpretation/decision nodes | New StudyContextID; preserve AnalysisID unless spec also changes |
| `DATA_OR_MAPPING_CHANGED` | Return to `data_check` | New/stale data binding and AnalysisID according to evidence contract |
| `DRAFT_SPEC_CHANGED` | Return to `estimate`; show draft/fitted difference | Preserve fitted AnalysisID; mark relevant use stale until fit/restore |
| `FIT_SUCCEEDED` | Bind current AnalysisID | Create new fitted identity/evidence lineage |
| `FIT_FAILED` | Expose repair and no-conclusion archive paths | Preserve versioned failure; no usable fitted claim created |
| `EVIDENCE_UPDATED` | Invalidate affected review completion | Recompute issue/readiness under existing evidence rules |
| `ARCHIVE_CREATED` | May complete workflow node | Create artifact identity only; no scientific state promotion |

## Invalidation matrix

| Changed input | Preserve | Invalidate or recompute |
|---|---|---|
| View, locale, theme, Help state | All scientific IDs and both specs | Rendered artifact only when its display content changes |
| Guide progress or answer | All scientific IDs and specs | Learning record only |
| Intended purpose/population/construct/use/stakes | AnalysisID and numerical evidence | StudyContextID, applicable claims, InterpretationID, DecisionID, artifacts |
| Rating-design/scoring-process narrative only | AnalysisID when analytic inputs are unchanged | Record IDs, applicability, interpretations/decisions that depend on them |
| Data or required mapping | Prior archive lineage | Current fit/evidence binding and downstream interpretation/decision |
| Draft estimator/model/identification setting | Existing fitted result | Draft-spec ID and current-use freshness |
| Successful refit | Prior immutable records | Current AnalysisID binding and downstream evidence/review/archive completion |
| Requested diagnostic/sensitivity | Unaffected evidence | Analysis/evidence identity according to existing computation contract |
| User-authored decision note | AnalysisID and evidence | DecisionID and rendered artifact |

Invalidation MUST explain what changed, which identity changed, why a node was
reopened, whether recomputation is required, and how to return to the fitted
specification.

## Standard and Detailed view contract

The stable internal IDs are `standard` and `detailed`.

- Standard view prioritizes one EvidenceIssue, its boundary, counts of other
  blocking/limiting issues, and the next action.
- Detailed view exposes resolved settings, full diagnostics, sensitivity
  ledgers, simulations, and exports.
- Both views consume the same canonical records and produce the same
  computation state for the same IDs.
- Switching view MUST preserve hidden draft values and MUST make zero estimator
  calls and zero previously unrequested diagnostic calls.
- Detailed view does not imply a stronger, more valid, or more research-worthy
  result.
- Primary execution controls MUST remain discoverable in the main content on a
  narrow screen rather than existing only in a collapsed sidebar.

## Audience vocabulary and surface exposure

Internal identifiers are required for correctness but are not the conceptual
model users must learn. One versioned terminology registry maps each stable
concept to its Standard label, Detailed/Help term, technical name, Japanese and
English aliases, definition, prohibited interpretations, and surfaces on which
it may appear.

| Internal contract | Guide/Standard wording | Detailed/Help wording | Technical-only form |
|---|---|---|---|
| `StudyContext` | Study purpose, people/domain, and intended use | Research context and interpretation/use argument | `StudyContextID` and schema |
| `RatingDesignRecord` | Rating assignment | Planned and observed rating design | Record ID/payload |
| `ScoringProcessRecord` | Rating process | Rater selection, training, qualification, monitoring, and rescoring process | Record ID/payload |
| `AnalysisSpec` | Analysis settings | Draft and fitted analysis specification | Canonical JSON/spec hash |
| `AnalysisID` | Analysis reference | Reproducible analysis identity | Full stable ID |
| `EvidenceIssue` | Current evidence question / what to check now | Evidence issue and tension | `EvidenceIssueID` |
| `ReasonCode` | Why this is limited / why action is needed | Stable limitation or applicability reason | Enum/code |
| `InterpretationRecord` | Interpretation draft | Bounded interpretation record | Record ID/payload |
| `DecisionRecord` | Use decision and next action | Decision/use record | Record ID/payload |

Scientific and SLA/measurement terms such as MFRM, logit, Infit/Outfit,
separation, residual PCA, and JMLE/MML are not implementation leakage. They MAY
appear in Detailed/Help when defined in plain language and bounded by
applicability. Standard MAY show the plain label followed by the canonical term
when it helps transfer learning.

State names follow the same rule:

| Internal state | Default user label | Required clarification |
|---|---|---|
| `AVAILABLE` | Available to review | Computed/applicable is not necessarily supported or stable |
| `CAUTION` | Review with a limitation | Name the limitation and next check |
| `HOLD` | Stop interpretation and repair | Name the prerequisite; allow a no-conclusion record |
| `NOT_ASSESSABLE` | Cannot be assessed with the current data/model | Distinguish from not yet assessed |
| `NOT_ASSESSED` | Not assessed yet | Do not render as no problem |
| `STABLE` | Maintained across the specified checks | Name the prespecified set; do not imply universal generalization |

Guide and Standard copy MUST NOT require users to understand class names,
enum names, schema versions, payloads, cache keys, manifests, fixtures,
pipelines, backend names, or runner names. Exact IDs and support codes MAY be
revealed in an optional technical-details disclosure and copied into a
privacy-safe support bundle. They MUST NOT be embedded in public URLs,
telemetry, or Help search queries.

`standard` and `detailed` remain stable internal view IDs. Their localized
labels are tested for the implication that one is statistically recommended or
more valid; Japanese may use an equivalent such as “要点表示 / 詳細表示” rather
than transliterating the internal IDs.

## `HelpTopic`, `HelpLink`, and contextual Help

Help is a versioned knowledge layer shared by the guide, Standard view,
Detailed view, popovers, warnings, EvidenceIssues, reporting, and
troubleshooting. It is not a parallel tutorial or a collection of page-local
strings.

Static Help MUST be reachable from the landing/no-data state, mapping and
preflight, first-run failure, every `HOLD`/`NOT_ASSESSABLE` path, and fitted
results. Rendering static Help MUST NOT require a current AnalysisID. The
current-analysis block is additive and appears only under a verified binding.

### `HelpTopic`

| Field | Meaning |
|---|---|
| `help_topic_id` | Stable, language-independent topic identity |
| `concept_ids` | Canonical concepts explained by the topic |
| `claim_boundary_ids` | Boundaries every rendering of the topic must preserve |
| `audience_layers` | Standard summary, Detailed/method explanation, technical/reviewer detail |
| `applicability` | Supported models, estimators, rubric structures, design states, and StudyContext conditions |
| `prerequisites` | Concepts/data/evidence needed before the topic can be applied |
| `static_content_keys` | Locale keys for method explanation independent of a current run |
| `dynamic_template_keys` | Locale keys for a current-analysis explanation bound to verified IDs |
| `related_stage_ids` | Plan/Check/Estimate/Diagnose/Stress/Decide/Archive locations |
| `related_target_ids` | Controls, evidence issues, charts, and reports that may link here |
| `search_alias_keys` | Plain-language, Japanese, English, acronym, and corrected legacy aliases |
| `method_reference_ids` | Versioned references and screening-convention provenance |
| `owner_and_review` | Content owner, reviewers, version, last review, and review state |
| `accessibility_requirements` | Text alternatives and reading-order requirements for figures, equations, and tables |

Each topic includes a stable template:

1. the research/user question;
2. prerequisites and applicable designs;
3. what the application computes;
4. what the result can show;
5. what it cannot show;
6. important alternative explanations or limiting conditions;
7. the next EvidenceIssue or check;
8. a reversible next action;
9. bounded reporting language and language to avoid;
10. method references and threshold provenance; and
11. related screens plus a return target.

Standard renders a short projection of this template. Detailed/Help may add
method and reference depth from the same topic. It MUST NOT fork a second
substantive interpretation.

### Static and current-analysis Help

“About this method” is static, citation-ready content that can be read without
an analysis. “For this analysis” is dynamic content rendered only when its
StudyContextID, AnalysisID, EvidenceIDs/EvidenceIssueID, and ClaimBoundaryIDs
match the current context. A stale, missing, sample-versus-real, or mismatched
binding MUST suppress the dynamic assertion and state that only general method
Help is being shown.

Help MUST NOT compute a previously unrequested diagnostic merely to populate a
dynamic paragraph. It reports `NOT_ASSESSED` and points to the explicit action
that would request the computation.

### `HelpLink`

| Field | Meaning |
|---|---|
| `help_link_id` | Stable link identity |
| `source_target_id` | Calling guide node, EvidenceIssue, chart, warning, setting, or artifact section |
| `help_topic_id` | Exact registered topic |
| `section_id` | Optional stable subsection within the topic |
| `context_binding` | Optional current StudyContext/Analysis/Evidence references, never raw data |
| `return_target_id` | Exact screen/card/control to restore |
| `return_focus_id` | Focus target after closing/returning |

The user-facing control is an actionable “Open detailed help” or equivalent,
not text instructing the user to find a page manually. Opening it selects the
exact topic and section, announces the new heading after Streamlit rerun, and
offers a return action. View and locale changes retain the same topic and
origin when valid.

An unknown or unavailable topic MUST fail the registry contract test. At
runtime it renders a localized generic explanation, support reference, and
safe return action; it MUST NOT silently no-op. Claim-critical limitations
remain visible on the source card even when Help is unavailable.

Help opening, search, navigation, and return events MAY be stored as
privacy-safe navigation state. They MUST NOT count as `learning_complete`,
alter scientific state, or include AnalysisID, file/column/person/rater names,
fingerprints, or free text in a URL or telemetry event.

### Help information architecture

The Help home supports search and browsing by:

- research question and intended use;
- study context, rating design, scoring process, and missingness;
- model and estimation;
- diagnostics;
- interpretation, differential interaction, and fairness boundaries;
- reporting, thesis/manuscript writing, and reviewer questions;
- reproducibility, privacy, and ethics;
- troubleshooting; and
- glossary, bilingual aliases, and references.

Developer inventories, raw filenames, schema explanations, and support/reviewer
tables belong in optional technical details or separate companion artifacts,
not the default Analysis Workflow Help panel.

Related-screen instructions are resolved from stable target IDs into the
current view's localized label; Help copy does not hard-code Full/Standard tab
names. Glossary definitions, inline term explanations, search aliases, and
exports derive from one canonical terminology source rather than parallel
Python and locale lists.

## Sample isolation and learning transfer

The deterministic sample uses the production estimator and evidence path but
has its own tutorial context. Starting, resuming, or restarting it MUST snapshot
and restore any current real-data draft/fitted specifications, identities,
evidence, and decision notes without overwrite.

The initial sample SHOULD use a synthetic holistic SLA rating study with a
clear Person/Rater/Task structure. A second unfamiliar scenario SHOULD contain
one interpretable `CAUTION` or `NOT_ASSESSABLE` issue and ask the learner to
identify:

1. what the evidence supports;
2. what it does not support;
3. what should be inspected next; and
4. when supervisory or methodological review is needed.

Learning answers MUST NOT change analysis state. The expected output, sample
data fingerprint, specification, evidence IDs, and feedback rubric are
versioned. An analytic-rubric scenario belongs in a later exercise so it can
teach criterion dependence and the boundary on an overall measure.

## SLA interpretation and ethics rules

- A declared analytic rubric is not assumed unidimensional because criteria
  were entered as a facet.
- Residual PCA and DIMTEST are bounded screens; absence of a large signal does
  not prove unidimensionality.
- Reliability/separation does not establish rater agreement,
  interchangeability, or accuracy.
- Differential-interaction evidence is a screening result, not a finding of
  bias or fairness. Subgroup purpose, reference group, overlap, multiplicity,
  uncertainty, small cells, privacy, accessibility, and consequences remain
  part of the interpretation boundary.
- Rater review proceeds from statistical flag to observation/design check,
  qualitative rating/rubric review, and a predeclared operational policy. It
  does not jump to exclusion, punishment, or retraining.
- Real-data entry presents deployment-specific processing, cache/retention,
  de-identification, consent/secondary-use, access, and prohibited-upload
  guidance before upload.
- Acknowledging governance guidance records only that it was shown. It does
  not provide ethical approval, consent, or authorization.

## Artifact separation

### Scientific artifacts

1. `MFRM_Analysis_Brief.md` or HTML: primary, human-readable, explicitly
   `draft_for_review` unless independently reviewed.
2. `MFRM_Analysis_Record.json`: deterministic sidecar with structured IDs,
   evidence, bounded claims, limitations, actions, versions, and privacy
   declaration.
3. Optional reviewer evidence pack: Methods/Results/Limitations mapping,
   planned-versus-observed design, failed sensitivity specifications, and
   unresolved warrants.

### Non-scientific companion artifacts

- education worksheet, instructor rubric, sample manifest, and learning
  checkpoint;
- privacy-safe support bundle containing app/schema version, stage/node ID,
  stable reason codes, sample-versus-real-data flag, and identity match status.

Learning/support artifacts MUST NOT be inserted into EvidenceRecords and MUST
exclude raw data, identifiers, file/column names, fingerprints, user free text,
and person/rater labels by default.

## Acceptance tests before visible migration

### Contract and identity

- Canonical serialization is deterministic and rejects unknown core fields.
- Changing locale/view/guide state leaves every scientific identity unchanged.
- Changing intended use changes StudyContextID/InterpretationID/DecisionID but
  not AnalysisID when the analytic specification is unchanged.
- Changing a computational setting preserves the old fitted record, changes
  the draft identity, and produces a new AnalysisID only after refit.
- Stable ClaimIDs have locale parity; translated wording is not used as an ID.

### State transitions

- Every lifecycle event has an exhaustive reducer test.
- `workflow_complete`, `learning_complete`, and `claim_ready` cover all valid
  independent combinations.
- Acknowledgement and formative answers cannot clear `HOLD`, change stability,
  or make a claim ready.
- A failed fit and required `NOT_ASSESSABLE` evidence can generate a valid
  no-conclusion record.
- Data/spec/context changes reopen the first affected node with a reason code.

### View and sample safety

- Standard -> Detailed -> Standard round trips preserve every hidden setting.
- View/language/Help switching invokes zero estimator calls.
- A sample journey never overwrites an existing real-data analysis or draft.
- Current focus never hides aggregate blocking/limiting counts or a registered
  evidence tension.

### Help, terminology, and error disclosure

- Every guide node, EvidenceIssue type, contextual popover, warning, setting,
  and report explanation has a valid `HelpLink` to a registered `HelpTopic` or
  an explicit no-link rationale.
- Help is reachable with no data, before estimation, after first-run failure,
  under every blocking/unassessable state, and after estimation.
- Opening Help reaches the exact section in one action and returning restores
  the source card/control and keyboard focus.
- Help open/search/return does not complete learning, mutate scientific state,
  or request a diagnostic.
- Standard and Detailed projections resolve to the same concept and
  ClaimBoundaryIDs; Detailed adds depth only.
- A stale, mismatched, or sample-versus-real binding never renders dynamic
  current-analysis Help as current.
- Unknown/unavailable topics produce a localized safe fallback and fail the
  registry test; no Help control silently disappears.
- The terminology registry covers every exposed contract/status term, and the
  guide/Standard surface contains no unapproved raw enum, class, schema,
  payload, cache, fixture, pipeline, backend, manifest, environment-variable,
  or runner language.
- Japanese and English search aliases reach the same topic; corrected legacy
  terms such as “bias” or “perfect fit” lead to the bounded explanation rather
  than preserving the unsafe claim.
- Routing uses stable topic/section/target IDs rather than translated labels or
  substring matching, and all glossary surfaces derive from one source.
- Claim-critical limits remain on the source card and are not hidden only in
  Help.
- Production errors show a localized problem, reversible next action, and
  privacy-safe support reference without raw exception text, backend symbols,
  environment overrides, uploaded values, or identifiers.

### Research and privacy

- Holistic, analytic-single-construct, analytic-subconstruct, and unsupported
  multidimensional scenarios produce the intended claim boundaries.
- Planned missingness and every unexpected reason code preserve denominators
  and cannot collapse into one anonymous `NA` state.
- Rater and subgroup flags contain non-punitive/fairness boundaries.
- Brief, JSON, education, and support artifacts pass separate disclosure
  allowlists.

### Bilingual and accessibility

- Japanese and English have the same stable concept, node, reason, boundary,
  and ClaimIDs.
- Meaning parity is reviewed, not inferred from key-count parity alone.
- Each visible slice receives browser checks for keyboard operation, focus and
  status after rerun, headings/reading order, narrow layouts/zoom, non-colour
  status, touch targets, and text/table alternatives for figures.

## Implementation sequence

1. Approve this contract and resolve any open semantic decision before adding
   Streamlit controls.
2. Freeze current behavior with characterization tests and complete the copy
   inventory.
3. Implement pure records, canonical serialization, identity rules, state
   reducer, audience terminology registry, and `HelpTopic`/`HelpLink` registry
   under `mfrm_app/` without importing Streamlit or pandas in the guidance or
   Help-contract layer.
4. Add a narrow adapter from existing session/widget state to canonical draft
   and fitted records; retain compatibility wrappers.
5. Implement exact contextual Help routing/return plus isolated sample context
   and deterministic guide journeys.
6. Project Standard and Detailed views over the same state.
7. Add EvidenceIssue, Analysis Brief, JSON sidecar, and companion artifacts.
8. Remove legacy routes only after bilingual, accessibility, privacy, and
   end-to-end replacement evidence passes.

Statistical changes, state extraction, visible UI changes, copy remediation,
and legacy removal remain separate review units.

## Deferred decisions

- Whether a privacy-safe learning checkpoint should be exportable/importable
  before any browser persistence is introduced.
- Whether `Interpret` eventually deserves its own top-level navigation stage.
- Which SLA domain profiles beyond the initial synthetic writing/speaking
  cases warrant versioned applicability rules.
- Which reviewer or education companion formats have enough observed demand
  to become maintained public contracts.

These decisions do not block the P0 record and state boundary.
