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

### P1 — Data-source selection will not scale

All sample scenarios, generation, paste, and upload share one flat radio group.
This prevents a hidden scenario selector today, but each new scenario lengthens
the global setup surface.

Long-term direction: first choose the source class (example, generate, paste,
upload), then show an always-visible scenario selector for the example class.
Preserve stable internal IDs and session migration so saved help routes and
tests do not depend on display labels.

### P1 — Navigation inside navigation

The post-fit page asks for a current goal, offers three route buttons, presents
a first-read overview, and then exposes a seven-section switcher. This can make
the user choose how to navigate before they understand what the evidence says.

Long-term direction: one persistent journey header with exactly one primary
next action. Goal selection may change the recommendation, but should not add a
second navigation system. Deep evidence remains addressable by stable section
and focus IDs.

### P1 — Monolithic UI ownership

`streamlit_app.py` is approximately 70,000 lines. Pure contracts are already
moving under `mfrm_app/`, but substantial rendering, copy, state transitions,
and computation remain interleaved. This increases the chance that a local UX
change causes rerun, localization, privacy, or stale-state regressions.

Long-term direction: keep scientific computation and evidence records pure;
add small presentation-policy modules; move source, setup, result-shell, and
export adapters behind tested boundaries without rewriting the estimator.

### P2 — Accessibility needs task-level acceptance tests

CSS-level focus and motion safeguards are useful but insufficient. The next
gate should cover keyboard completion of the sample run, focus after rerun,
200% zoom, narrow viewport hierarchy, non-color status text, and accessible
names for every primary action.

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
- Replaced explanation-first onboarding with two explicit starting choices;
  the three-step detail is now optional.
- Moved raw input rows behind a collapsed disclosure.
- Consolidated input preview, sample identity, readiness, and response-row
  audit into one setup workspace.
- Cleared setup content immediately after a successful fit; on later reruns it
  remains available in one collapsed panel.
- Added regression tests for the policy, privacy severity, onboarding hierarchy,
  and setup/result progressive disclosure.

## Long-term implementation sequence

1. Stabilize one workflow-shell contract with a single primary next action.
2. Replace the flat data-source list with a scalable two-level chooser while
   preserving stable IDs and old session state.
3. Collapse the goal router, action hub, and section navigator into one route
   model backed by the existing Help target registry.
4. Extract source/setup/result-shell renderers from the monolith, one tested
   vertical slice at a time.
5. Add task-level accessibility and rendered-locale acceptance tests.
6. Instrument privacy-safe UX events: source class, phase reached, blocked
   reason code, rerun count, and time to first interpretable evidence. Never
   record uploaded values, person identifiers, free text, or raw file names.

## Acceptance measures

- A first-time user can run the sample from the landing page with one action.
- A user with their own data can identify the next required setup action without
  opening a tutorial.
- No raw response row is visible by default.
- Synthetic data does not produce the same alert severity as user data.
- After fitting, exactly one surface names the highest-priority next check.
- Changing display density never changes the fitted analysis identity.
- Every blocked route provides one reversible action and one stable Help target.
- The sample-run journey is completable by keyboard at 200% zoom and a narrow
  viewport without losing the current phase or primary action.
