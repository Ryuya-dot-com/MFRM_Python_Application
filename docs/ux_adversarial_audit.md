# Adversarial UX Audit

- Audit date: 2026-08-04
- Scope: default Streamlit journey from landing through the first fitted result
- Product boundary: standalone Python beta

## Executive finding

The application is unusually strong at preserving statistical evidence,
warnings, reproducibility details, and advanced controls. Its main UX risk is
now the inverse of feature scarcity: too many individually reasonable aids
compete to be the user's starting point.

The default journey should optimize for one decision at a time:

```text
Choose data -> Confirm mapping and readiness -> Run -> Read the first blocker
-> Inspect one evidence surface -> Export only after claim review
```

Scientific evidence must remain available, but availability does not require
simultaneous prominence. Workflow progress, statistical readiness, learning
progress, and claim readiness must remain separate state systems.

## Current-state evidence

The automated initial render completed without an exception. The default
sample was correctly detected as 960 rows, 30 persons, three facets, and a
ready-to-run design. The same render tree also contained 7 radio groups, 12
select boxes, 4 multi-selects, 19 expanders, two competing run actions, and
multiple summaries of the same loaded sample.

After the sample fit, the default result journey contained a goal selector,
three action buttons, a first-read overview, a second section switcher, run
history, a quick download, and repeated response-data audit content. Each
surface is defensible in isolation; together they create navigation inside
navigation and make it harder to identify the single next decision.

## What is already working

- The app fails closed on stale results and keeps restored-run identity visible.
- Compact and full views share the fitted statistical result instead of fitting
  different models silently.
- Readiness, first-read, claim boundaries, and privacy-safe export contracts are
  materially stronger than a typical exploratory dashboard.
- Keyboard focus styling, reduced-motion handling, bilingual locale parity, and
  narrow-screen fallbacks already have regression coverage.
- Expensive result panels are increasingly rendered lazily.

These strengths should be retained while simplifying the route around them.

## Adversarial findings

### P0 — Warning habituation

The landing page previously showed the same strong privacy warning for a
synthetic built-in sample and a real uploaded rating file. Repeated false-alarm
severity teaches users to ignore the warning precisely when it matters.

Required rule: warning severity follows data origin. Synthetic data receives a
quiet provenance notice; paste, upload, and unknown future sources receive the
strong handling warning.

### P0 — Raw-row exposure and setup dominance

The first 20 raw rows were visible by default. For real data this increases
shoulder-surfing risk; for every data source it pushes the actual readiness and
run decision down the page. After fitting, the same setup content remained
above the result journey.

Required rule: raw rows are opt-in, and setup becomes subordinate after a
result exists. A user can always reopen the input, row audit, and mapping.

### P0 — Duplicate guidance after fitting

Pre-run readiness and row-audit content was followed by another row audit in
the guided Start section. The result surface then added a goal router, action
hub, first-read plan, and section navigation. Repetition did not add evidence;
it weakened hierarchy.

Required rule: each phase owns one primary orientation surface. Setup owns
readiness. Results own first-read priority. Detail sections own evidence.

### P0 — First-run explanation without learning state

The earlier landing combined two actions, an optional three-step explanation,
a separate terminology tutorial, and a static five-row route shown only after
fitting. None of those surfaces knew which task had been reviewed, whether the
user had exited, or whether a sample answer confused a completed fit with a
valid claim.

Current remediation: a new session stops at one three-route landing. The
optional sample path uses five stable nodes and one primary action per node;
the ordinary workspace remains directly available. A pure reducer stores
learning progress separately from AnalysisID and evidence, and the formative
question cannot promote scientific readiness. The old tutorial renderer is no
longer part of `main()`; it remains dormant only as a compatibility surface
until replacement evidence is complete.

### P1 — Data-source selection will not scale

All sample scenarios, generation, paste, and upload share one flat radio group.
This prevents a hidden scenario selector today, but each new scenario lengthens
the global setup surface.

Long-term direction: first choose the source class (example, generate, paste,
upload), then show an always-visible scenario selector for the example class.
Preserve stable internal IDs and session migration so saved help routes and
tests do not depend on display labels.

Current remediation: the four source classes now own the first decision, and
the example scenario selector is always visible immediately below its class.
The old `scenario:*`, `simulate`, `paste`, and `upload` projection remains the
analysis-facing contract. Restored old state migrates forward, privacy severity
uses the new class on the same rerun, and the last selected sample survives a
temporary move to paste, upload, or simulation.

### P1 — Navigation inside navigation

The post-fit page asks for a current goal, offers three route buttons, presents
a first-read overview, and then exposes a seven-section switcher. This can make
the user choose how to navigate before they understand what the evidence says.

Long-term direction: one persistent journey header with exactly one primary
next action. Goal selection may change the recommendation, but should not add a
second navigation system. Deep evidence remains addressable by stable section
and focus IDs.

Current remediation: the duplicate first-read overview was removed, and the
single Essential/All-panels result selector now becomes a sticky navigation
dock during long-page scrolling. The separate keyboard-shortcut cheat sheet
was removed rather than introducing another interaction system. Narrow layouts
keep the Essential selector on one touch-scrollable line instead of a tall
fixed overlay.

### P1 — Monolithic UI ownership

`streamlit_app.py` is approximately 70,000 lines. Pure contracts are already
moving under `mfrm_app/`, but substantial rendering, copy, state transitions,
and computation remain interleaved. This increases the chance that a local UX
change causes rerun, localization, privacy, or stale-state regressions.

Long-term direction: keep scientific computation and evidence records pure;
add small presentation-policy modules; move source, setup, result-shell, and
export adapters behind tested boundaries without rewriting the estimator.

### P2 — Accessibility and persistent navigation need task-level acceptance tests

CSS-level focus, motion, and sticky-position safeguards are useful but
insufficient. The next gate should cover keyboard completion of the sample
run, focus after rerun, navigation-dock position after long scrolling, narrow
viewport hierarchy, non-color status text, and accessible names for every
primary action. Browser zoom remains a reflow compatibility check; it is not a
separate application feature.

### P2 — Mixed-language technical surfaces

Locale parity ensures matching key topology, but a number of older technical
strings are still embedded directly in the entrypoint. Japanese users can
therefore encounter English during errors, advanced controls, and result
details. Translation coverage should be measured at the rendered-route level,
not only by JSON key equality.

## Changes delivered in this pass

- Added a pure `mfrm_app.ux` contract for data-origin classification and
  phase-aware presentation policy.
- Made privacy severity contextual and fail-closed for unknown future sources.
- Replaced the overlapping onboarding, three-step explanation, and visible
  terminology tutorial with one optional five-step sample route and a direct
  ordinary-workspace escape.
- Added a pure `mfrm_app.guidance` catalog and reducer for skip, exit, resume,
  restart, invalidation, fit binding, formative review, and completion without
  a claim-readiness field or estimator dependency.
- Moved raw input rows behind a collapsed disclosure.
- Consolidated input preview, sample identity, readiness, and response-row
  audit into one setup workspace.
- Moved the Guided-defaults primary Run action from the technical sidebar into
  that setup workspace after mapping and readiness. The main CTA uses
  goal-oriented bilingual wording and a compact model/method/depth summary;
  Advanced controls retain the expert sidebar route without duplicating both.
- Reduced Guided setup to the decisions needed for a first defensible run:
  Person/Score/facet roles, model, estimation method, and analysis coverage.
  Weighting, output styling, regularization, population modeling, score-scale
  overrides, identification/optimizer controls, anchors, and report scaling
  now require Advanced controls. Switching back resets hidden technical values
  to documented Guided defaults; a detected weight column is excluded from
  facet suggestions and explained rather than silently used.
- Added a versioned, fail-closed browser accessibility contract and runbook.
  The catalog expands eight tasks across English/Japanese and seven targeted
  profiles into 52 cases and 482 required evidence rows. Templates and AppTest
  cannot create a passing decision; complete version-matched browser evidence
  is required. Static preflight broadened focus styling, forced-colors and
  coarse-pointer handling, and compact-table semantics, while leaving actual
  keyboard, reflow, contrast, accessibility-tree, touch, and screen-reader
  acceptance explicitly open.
- Cleared setup content immediately after a successful fit; on later reruns it
  remains available in one collapsed panel.
- Added regression tests for the policy, privacy severity, onboarding hierarchy,
  and setup/result progressive disclosure.

## Report & Export revision — 2026-09-12

The default export entry previously placed report summaries, several audience
boards, editable checklists, claim traces, and nested report tabs before the
document and data downloads. Users had to understand the internal reporting
structure before saving a file.

The entry now starts with three tasks: **Create a report**, **Save tables &
figures**, and **Review checks**. The report task shows one format choice and
one download button: PDF for reading/printing, Word for editing, or HTML for
browser viewing. Only the selected format is built. Results text is optional;
methods, claim guidance, work notes, and the legacy reporting tools are explicit
choices under Review checks. Checks open as readable action/evidence sections,
with the highest-priority issue open first; long actions are not clipped in a
wide grid. Individual CSV export uses a searchable table selector instead of a
separate button for every table.

The inference hold and publication checks remain visible before download.
Unavailable check results are identified explicitly. The existing privacy
filter for tables remains enabled by default; the separate report document
discloses that its English text and individual results need review before
sharing. Full result content now determines document cache reuse, so a change
inside a same-sized table invalidates the saved document. Document and single
CSV downloads do not rerun the app.

AppTest covers both locales, each task, opt-in advanced material, free-SD MML
holds, a blocked or unavailable publication check, format-specific generation,
same-shape content changes, recovery after export failure, and removing a
selected person-level CSV when public export is enabled again. The table
selector starts with the analysis summary when available.

Validation: **150 tests passed**, covering the new interactions and the existing
export, readiness, publication-figure, and locale suites. In a local Chrome
152.0.7977.83 session, the built-in 960-observation RSM/JMLE example retained
fingerprint `067a6d8d` across report tasks and document formats. PDF, DOCX, HTML,
and the selected summary CSV were downloaded and their file content checked.
The default document task was visually checked at desktop width and at
400/320 CSS pixels in English and at 320 pixels in Japanese. The Japanese
check details were also inspected at 320 pixels. These views had no document
or main-panel horizontal overflow. These focused checks do not constitute acceptance of the separate
browser accessibility protocol or a first-time-user study. The diagnostic
evidence wording itself remains English; the new task and output controls
support both English and Japanese.

Follow-up observed during browser review: the existing lightweight language
switch preserves the fitted result but can return the section selector to
Start after redraw. Retaining the user's location belongs in the shared
navigation work below; this revision does not claim to resolve that behavior.

## Long-term implementation sequence

1. Execute the versioned browser matrix for keyboard/focus acceptance and
   exact contextual Help return, retaining evidence for every required row.
2. Complete browser and first-time-user acceptance for the implemented
   two-level source chooser and compact Guided setup while preserving stable
   IDs and old session state.
3. Collapse the goal router, action hub, and section navigator into one route
   model backed by the existing Help target registry.
4. Extract source/setup/result-shell renderers from the monolith, one tested
   vertical slice at a time.
5. Add task-level accessibility and rendered-locale acceptance tests.
6. Instrument privacy-safe UX events: source class, phase reached, blocked
   reason code, rerun count, and time to first interpretable evidence. Never
   record uploaded values, person identifiers, free text, or raw file names.

## Acceptance measures

- A first-time user can finish the five-step sample route without encountering
  unrelated source, model, or export controls before they are needed.
- A user with their own data can identify the next required setup action without
  opening a tutorial.
- No raw response row is visible by default.
- Synthetic data does not produce the same alert severity as user data.
- After fitting, exactly one surface names the highest-priority next check.
- Changing display density never changes the fitted analysis identity.
- Every blocked route provides one reversible action and one stable Help target.
- The sample-run journey is completable by keyboard on a narrow viewport, and
  long-page scrolling never hides the current result selector or primary
  action. Browser zoom is checked only as reflow compatibility.
