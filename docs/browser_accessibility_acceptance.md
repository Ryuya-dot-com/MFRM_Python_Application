# Browser Accessibility Acceptance Protocol

- Protocol version: `mfrm_browser_accessibility_acceptance_v1`
- Target: WCAG 2.2 AA within the supported Streamlit/browser scope
- Data boundary: synthetic fixtures only
- Current decision: **NOT_ACCEPTED — browser evidence has not been run**

## Purpose

This protocol converts browser review from an informal visual inspection into
a versioned, task-level release gate. Static source checks and Streamlit
AppTest remain useful preflight evidence, but they cannot establish keyboard,
focus, reflow, touch, contrast, accessibility-tree, or screen-reader behavior.

The executable catalog lives in
`mfrm_app/accessibility_acceptance.py`. It currently expands eight user tasks,
two locales, and seven targeted browser profiles into 52 cases. Each case has
an explicit set of required checks and evidence types.

## Fail-closed decision rule

Every required case/check pair starts as `NOT_RUN`.

- `PASS`: all evidence types declared for the check are present as portable
  relative paths, and observation, UTC time, browser build, and reviewer are
  recorded.
- `FAIL`: retained failure evidence and audit metadata are required.
- `BLOCKED`: a retained moderator note and audit metadata are required.
- `NOT_RUN`: carries no evidence and never contributes to acceptance.
- `NOT_APPLICABLE`: is not permitted for required pairs. Inapplicable checks
  are omitted when the case is generated instead.

The overall result is `ACCEPTED` only when every version-matched required row
is `PASS`. Missing, stale-fingerprint, duplicate, blocked, failed, and not-run
records all yield `NOT_ACCEPTED`.

## Browser profiles

| Profile | CSS viewport | Zoom | Input / condition |
|---|---:|---:|---|
| `desktop_keyboard` | 1440 × 900 | 100% | Keyboard only |
| `narrow_keyboard` | 320 × 568 | 100% | Keyboard only |
| `narrow_touch` | 320 × 568 | 100% | Touch / coarse pointer |
| `zoom_200_keyboard` | 640 × 450 | 200% | Keyboard only |
| `zoom_400_keyboard` | 320 × 225 | 400% | Keyboard only |
| `reduced_motion_keyboard` | 1440 × 900 | 100% | Keyboard + reduced motion |
| `screen_reader_keyboard` | 1440 × 900 | 100% | Keyboard + platform screen reader |

CSS viewport dimensions—not physical display pixels—are the governing reflow
measure. The exact browser/version is stored per evidence row. The connected
in-app Chromium surface can cover ordinary DOM, keyboard, viewport, zoom, and
screenshots. Actual screen-reader announcement cases remain blocked until a
platform screen reader can be observed; a DOM accessibility snapshot is not a
substitute for that observation.

## User tasks

1. Choose the sample-learning route from a new session.
2. Complete the five-step sample route and save its learning-only checkpoint.
3. Paste synthetic fixture data, confirm mapping/readiness, and run Guided
   defaults.
4. Load a synthetic Weight column and understand the Guided equal-weight
   explanation.
5. Switch to Advanced controls and configure a synthetic anchor table.
6. Run the sparse sample and follow its highest-priority reversible action.
7. Navigate a fitted result and obtain its ZIP archive.
8. Open contextual Help and return to its exact registered source heading.

Every task is generated in English and Japanese. Narrow/zoom profiles require
page-level horizontal-scroll and sticky-obscuration checks. Japanese cases add
long-label reflow. Keyboard, touch, reduced-motion, and screen-reader profiles
add their corresponding checks automatically.

## Evidence types

- `dom_snapshot`: current rendered structure, dimensions, focus target, and
  overflow facts needed for the criterion.
- `screenshot`: visual state at the specified viewport/zoom; crop only when a
  full-page image is also retained.
- `keyboard_trace`: ordered keys, active element/accessibility name after each
  material step, and final task state.
- `accessibility_tree`: roles, names, states, headings, and reading order.
- `moderator_note`: concise observed behavior; required for blocked conditions
  and real screen-reader/reduced-motion observations.
- `download_artifact`: synthetic-data result archive produced by the task.
- `app_identity_record`: AnalysisID and estimator-call evidence before/after a
  presentation-only action.

Evidence must remain under the acceptance directory and use relative paths.
Absolute machine paths, URLs, parent-directory traversal, raw ratings, real
identifiers, file names supplied by a participant, and free-text user data are
rejected by the contract or prohibited by this protocol.

## Execution procedure

Export a version-matched blank bundle:

```bash
python3 -m mfrm_app.accessibility_acceptance export \
  --output validation/browser_accessibility_acceptance_YYYYMMDD
```

This creates:

- `browser_acceptance_matrix.csv`
- `browser_acceptance_checks.csv`
- `browser_evidence_template.csv`
- `browser_acceptance_manifest.json`

The manifest correctly remains `NOT_ACCEPTED`. For each matrix row:

1. Start from a fresh session or the task's documented prerequisite state.
2. Set locale, CSS viewport, zoom, input mode, reduced-motion preference, and
   screen reader exactly as specified.
3. Complete the task without using an unlisted input mode.
4. Collect only the evidence types required for each check.
5. Store evidence beneath `evidence/<CaseID>/` and fill the corresponding CSV
   row. Do not change case IDs, check IDs, schema version, or catalog
   fingerprint.
6. Repeat a failed case after a fix, retaining the failed run separately rather
   than overwriting it. Only the reviewed final evidence belongs in the gate
   CSV; the failure history remains in an adjacent remediation directory.

Evaluate completed evidence:

```bash
python3 -m mfrm_app.accessibility_acceptance evaluate \
  --evidence validation/browser_accessibility_acceptance_YYYYMMDD/browser_evidence_completed.csv \
  --output validation/browser_accessibility_acceptance_YYYYMMDD/browser_acceptance_decision.json
```

The evaluator exits nonzero unless every required row is `PASS`.

## What the static preflight establishes

The repository can establish the following before browser connection:

- required Streamlit widget labels are not deliberately collapsed;
- the app defines no custom keyboard shortcuts;
- focus styling covers links, buttons, disclosures, native inputs, common ARIA
  widget roles, and positive `tabindex` values;
- forced-colors and reduced-motion preferences have explicit CSS handling;
- coarse-pointer buttons/disclosures receive a 44 CSS-pixel minimum height;
- compact HTML tables have an accessible name and explicit column/row header
  scope; and
- readiness states use text prefixes rather than color alone.

These are implementation prerequisites, not browser acceptance evidence.

## Browser-only gates that remain open

- actual focus visibility and contrast across Streamlit themes;
- logical Tab/Shift+Tab order and focus behavior after reruns;
- accessibility-tree roles, names, heading hierarchy, and live announcements;
- reflow at 320 CSS pixels and 200–400% zoom without hidden required controls;
- sticky navigation that does not obscure content;
- measured target size and operability on a coarse pointer;
- long Japanese label wrapping; and
- one real screen-reader announcement per material state transition.

No release note or readiness surface may describe these gates as passed until
the completed evidence evaluates to `ACCEPTED`.
