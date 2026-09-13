# Contextual Help browser acceptance

- Contract: `mfrm_help_browser_acceptance_v1`
- Current pilot: `fit_scatter`
- Current stage: `APPTEST_ONLY`
- Browser evidence ID: none
- Browser result: **NOT RUN**

This runbook governs the real-browser promotion of one contextual Help adapter.
It does not turn AppTest output, a checklist, or a reviewer assertion into
browser evidence. `APPTEST_ONLY` and `BROWSER_ACCEPTED` are internal release
states and must never be rendered as user-facing copy.

## Promotion boundary

Only one adapter may be the AppTest-only pilot. All other active adapters must
be browser accepted and carry a stable evidence reference matching:

```text
browser.help.<topic_key>.<retained_record>
```

The reference must be bound to the same adapter topic. The retained-record
segment must begin and end with a lowercase letter or digit; lowercase letters,
digits, `.`, `_`, `:`, and `-` may appear inside it. A syntactically valid ID is
not evidence by itself: it must resolve to the reviewed, retained record for
the exact commit and deployment under test.

The current `fit_scatter` pilot remains active so its real journey can be
tested. A second AppTest-only popover must remain inactive. `fit_scatter` must
not be changed to `BROWSER_ACCEPTED` until every required case below passes and
the evidence is reviewed and retained.

## Fixed test context

Use one immutable commit SHA and one deployment URL for a complete run. Record
the exact Streamlit, OS, browser, and assistive-technology versions. Use only
the built-in sample or a dedicated synthetic fixture; never place real ratings,
names, file paths, column names, free text, or a real-data AnalysisID in the
evidence bundle.

Start with a fitted sample analysis at a desktop viewport of `1440 x 900` and
100% zoom. Give every core case a different valid fit convention/cap so state
restoration is observable.

| Case | Locale | View | Expected source after return |
|---|---|---|---|
| `CORE-01` | English | Essential | diagnostics -> fit details |
| `CORE-02` | Japanese | Essential | diagnostics -> fit details |
| `CORE-03` | English | Full | main results -> fit details |
| `CORE-04` | Japanese | Full | main results -> fit details |

For each case, use the keyboard to follow this journey:

1. Move to the Fit scatter reading-guide popover.
2. Open the exact full guide.
3. Scroll the Help surface away from the original result position.
4. Activate the localized return control.
5. Do not click or press another key until the return probe is captured.

Opening Help or returning must not estimate, refit, consume a pending
user-initiated run, or change the fitted AnalysisID/data context.

## Exact return probe

Run this read-only probe immediately after return. Save its JSON output as part
of the case evidence.

```javascript
(() => {
  const id = "mfrm-focus-results-figure-fit-scatter";
  const nodes = [...document.querySelectorAll(`#${CSS.escape(id)}`)];
  const heading = nodes[0] || null;
  const rect = heading?.getBoundingClientRect() || null;
  const style = heading ? getComputedStyle(heading) : null;
  return {
    count: nodes.length,
    tag: heading?.tagName || null,
    active: document.activeElement === heading,
    tabindex: heading?.getAttribute("tabindex") || null,
    focusId: heading?.dataset.mfrmFocusId || null,
    focusRequest: heading?.dataset.mfrmFocusRequest || null,
    focusStatus: heading?.dataset.mfrmFocusStatus || null,
    focusHandled: heading?.dataset.mfrmFocusHandled || null,
    rect: rect && {
      top: rect.top,
      bottom: rect.bottom,
      left: rect.left,
      right: rect.right,
      width: rect.width,
      height: rect.height,
    },
    viewport: {
      width: visualViewport?.width || innerWidth,
      height: visualViewport?.height || innerHeight,
    },
    outline: style && {
      style: style.outlineStyle,
      width: style.outlineWidth,
      offset: style.outlineOffset,
      color: style.outlineColor,
    },
  };
})()
```

Every core case passes only when all of the following are true:

- exactly one target exists and it is an `H3`;
- `document.activeElement` is that target;
- `tabindex` is `-1`;
- `focusId` is `focus.results.figure.fit_scatter`;
- `focusRequest` is `return-from-help`;
- `focusStatus` is `focused` and `focusHandled` is `true`;
- the accessible role/level is heading level 3;
- the accessible name is `Fit scatter: Infit vs Outfit` in English or
  `Fit 散布図: Infit と Outfit` in Japanese;
- the localized return caption occurs exactly once and there is no returned
  success alert, empty focus group, duplicate ID, or JavaScript exception;
- the expected panel, fit convention/cap, AnalysisID match, data-context match,
  pending triggers, and unrelated analysis settings are unchanged; and
- the one-shot controller and caption are absent on the next ordinary rerun.

## Viewport, focus appearance, and Tab order

Let `V` be the visual viewport height, `O` the lowest visible bottom edge of a
fixed/sticky Streamlit toolbar touching the viewport top, and `R` the returned
heading rectangle. Require:

```text
R.top >= O + 8
R.bottom <= V - 8
abs(R.centerY - ((O + V) / 2)) <= 0.20 * (V - O)
```

The complete four-sided focus ring must be visible in a screenshot. Computed
style must have a solid outline of at least 3 CSS pixels, an offset of at least
4 CSS pixels, and a nontransparent colour.

Press `Tab` once from the focused heading. Focus must move to the next native
sequential target in the rendered DOM (the heading permalink if it is
tabbable, otherwise the reading-guide control). The heading must not enter the
ordinary Tab sequence, return to the old Help button, or create a focus trap.

## Persistent result-navigation dock

This gate has priority over adding any application-level zoom or keyboard
shortcut feature. Run it in both Essential and All-panels views after fitting
the built-in sample.

### `NAV-01` — long result section

Open a result section that is at least three visual viewports tall, scroll
until the original selector position is no longer visible, and capture the
authoritative result-navigation dock rectangle. Pass only when:

- exactly one dock is rendered for the active view;
- computed `position` is `sticky`, its top edge clears the Streamlit toolbar,
  and the entire selector remains inside the visual viewport;
- the current selection and localized caption remain visible;
- selecting another section works with pointer and ordinary keyboard
  interaction, preserves AnalysisID, and does not estimate/refit; and
- the dock does not cover an alert, focused element, required field, or the
  first heading of the selected section.

### `NAV-02` — narrow result navigation

Repeat at `390 x 844` and 320 CSS-pixel widths. In Essential view, the section
choices must remain on one horizontally touch-scrollable line rather than
wrapping into a tall overlay. The selected choice, horizontal overflow, and
focus ring must be discoverable without page-level horizontal panning. In
All-panels view, the compact select control must remain fully visible.

Record dock, toolbar, selector, and first-section-heading rectangles before
and after the section change, plus the before/after AnalysisID and estimator
call count.

## Required adversarial and compatibility cases

### `REPEAT-01` — two legitimate returns

Open and return from the same guide twice. The second journey must create a new
heading node, focus it successfully, and emit one caption/controller. The
first node's handled marker must not suppress the second legitimate journey.

### `ABORT-01` — user focus wins

With a controlled browser fixture, move focus to a visible enabled app control
after the returned heading mounts but before the controller's second animation
frame. Log calls to the original `scrollIntoView` without replacing its
behaviour.

Pass only when the user's control remains active, the heading reports
`aborted-user-focus`, the heading is not focused or scrolled, and there is no
JavaScript exception. Only this case expects `aborted-user-focus`; normal
returns require `focused`.

### `STALE-01` — identity or provenance changes

Open Help from analysis A, then replace the retained fit with analysis B or a
different sample/real context before return.

Pass only when A-owned panel/control values are not written to B, no return
controller/caption/focus is produced, the localized safe fallback is shown,
and B remains unchanged.

### `CSP-01` — deployed response policy

Capture the actual preview/public response headers, meta policy, and console.
There must be no inline-script or inline-style CSP violation, and the complete
core focus probe must pass. An absent CSP may establish compatibility with the
current environment, but it is not evidence that a CSP defence is adequate.

### `MOTION-01` — reduced motion

Enable `prefers-reduced-motion: reduce` and capture the matching media-query
result. Core conditions must still pass, computed scroll behaviour must remain
`auto`, no multi-frame smooth animation may occur, and the focus outline must
remain visible.

### `REFLOW-01/02` — browser reflow compatibility

After `NAV-01/02` pass, run the same journey at desktop 200% browser zoom. This
is a reflow compatibility check, not an application zoom feature. The guide
and return controls must work without hover, required content and the focus
ring must be reachable without page-level horizontal panning, the toolbar and
navigation dock must not obscure the heading, and the same focus/state/caption
conditions must pass.

The exact-source return scroll in this runbook is an accessibility restoration
after an explicit Help journey. It is not a tutorial coach mark or automatic
tutorial progression. Smooth scrolling and selector-driven tours remain out of
scope.

## Screen-reader gate

Run at minimum:

- VoiceOver with Safari on macOS, English and Japanese; and
- NVDA with Firefox or Chrome on Windows, English and Japanese.

Across those runs include both Essential and Full views. Pass only when the
localized real heading receives focus, its accessible name and level are
announced once, no empty group or returned-status alert receives focus, the
caption remains available to browse/virtual-cursor navigation, and the next
Tab follows the logical sequence. Do not require a screen reader's surrounding
stock phrase to match verbatim; judge the accessible name, role, level, and
focus-event count.

## Evidence record

Create one reviewed record per case with:

- case ID, UTC timestamp, operator, reviewer, result (`PASS`, `FAIL`, or
  `BLOCKED`), defect reference, commit SHA, and deployment URL;
- Streamlit, OS, browser, and assistive-technology versions;
- locale, view, viewport, visual viewport, DPR, zoom, and reduced-motion state;
- stable Help link/target/focus IDs and the return probe JSON;
- before/after presentation-state comparison using sample-safe values;
- heading/toolbar/navigation-dock rectangles, scroll-call log, console log,
  and CSP headers;
- screenshots before return, immediately after return, and with the focus ring;
- short recordings for scroll, abort, reduced-motion, and screen-reader cases;
  and
- an accessibility-tree snapshot and a concise VoiceOver/NVDA transcript.

`NOT RUN`, `BLOCKED`, and inability to reproduce a case are not passes. Do not
create a `browser_evidence_id` until all four core cases plus repeat, abort,
stale, CSP, motion, reflow, VoiceOver, and NVDA have passed.

## Promotion procedure

1. Retain the reviewed evidence outside the app's user-facing surfaces.
2. Assign a stable reference such as
   `browser.help.fit_scatter.r3.20260724`.
3. Change the adapter to `BROWSER_ACCEPTED` and set that exact reference.
4. Run the acceptance-contract, contextual Help, full test, self-test, and
   deployment-health gates.
5. In a separate change, move the single AppTest-only pilot key to the next
   exact popover. Never have two AppTest-only pilots simultaneously.
