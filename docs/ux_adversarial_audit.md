# Adversarial UX Audit

Current checkpoint, 2026-10-16: the [Help readability review](#help-readability-and-internal-copy-review--2026-10-16)
removes duplicated tutorial explanations and development wording from the
reviewed Help and model-comparison surfaces. Automated rendering covers all
eight full Help panels in both languages; browser and human acceptance remain
open. This is an integration-branch update, not a deployment or release.

Earlier checkpoint, 2026-10-16: the [comparison and fresh-process follow-up](#local-comparison-and-fresh-process-follow-up--2026-10-16)
corrects unsupported agreement claims and preserves reliability holds in the
comparison view. The local fresh-process/changed-input timing targets passed
for one small fixture. Browser, human and broader-scale acceptance remain open.

Checkpoint, 2026-10-15: bounded report reuse and bilingual history
controls have passed the [local follow-up](#local-download-and-history-refinement--2026-10-15).
The first performance candidate missed its target and is retained as a failure;
the second passed both input-specific targets. Browser and human acceptance
remain open. The following dated records retain their original scope.

Checkpoint, 2026-10-14: local Japanese-copy, documentation-reference
and report-qualification corrections have passed their recorded automated
checks. Four fitted cases now have AppTest input-to-export coverage, with
additional input-recovery and restored-result checks. Changed Word/PDF table
pages were inspected in the earlier table review. AppTest is not browser
acceptance or evidence of user understanding. See the
[change record](../CHANGELOG.md#unreleased--directory-review-corrections).
UX1 follows the [execution roadmap](../ROADMAP.md#現在の実行計画--2026-09-15)
and the protocol below. Actual-user reviews and full accessibility acceptance
remain outstanding. This document records publication of 0.2.18-beta (`b483d97`)
on September 15; this update does not verify the current deployment.

Links into `validation/generated/` and `validation/source_snapshots/` identify
local retained evidence, not files distributed in a fresh checkout. Dated
screenshots and test results cover their recorded versions and tasks only;
they do not establish acceptance of the current working tree or hosted app.

## Next acceptance slice — 2026-10-14 (planned)

Status: **PARTIAL — local AppTest journeys executed; browser and human review
NOT RUN**. The heading retains the original planning anchor. The execution
record below limits what has been checked. This protocol reuses the existing
interface and tests. It does not add a new navigation system or
promote a contextual Help adapter. That promotion still requires
[the browser acceptance protocol](help_browser_acceptance.md).

Use synthetic inputs only. Start from the four real-fit cases in
[`test_measure_report_consistency.py`](../tests/test_measure_report_consistency.py),
the diagnostic-navigation tests, and the existing Help return tests. These
fixtures exercise behavior; they do not establish statistical operating
characteristics. Record the source fingerprint and environment for local
checks; formal browser acceptance must use the immutable candidate commit and
deployment required by the browser protocol.

| Case | User task and expected outcome | Existing evidence and remaining gap |
|---|---|---|
| Ordinary JMLE | Import the synthetic table, check columns, run, find one diagnostic, explain conditional uncertainty and export a bounded descriptive report | Real fit and exports tested; complete journey and unaided explanation remain to be observed |
| Extreme-score free facet JMLE | Recognize that optimizer termination does not remove the interpretation hold; find its reason and preserve it in the saved report | Hold propagation tested; whether the user notices it and chooses an appropriate next action is unknown |
| Fixed-SD MML | Identify person uncertainty as conditional posterior SD, distinguish it from sampling SE, and find the relevant limitation | Report basis tested; comprehension and return navigation remain to be observed. Fixed-SD PCM quadrature checks have a separate restricted scope |
| Free-SD MML | Locate the inference hold, explain why reported precision is restricted, and export an honest record without a final inference | UI/table/export checks pass; whole-journey recognition and explanation remain untested |
| Saved or stale result | Open an older result with positive uncertainty flags, then change the active input/settings; distinguish historical values from the current analysis | Row-level hold override tested; full restore/invalidation/export journey must be checked without overwriting the original |
| Recoverable input error | Paste a synthetic table with an invalid score/column choice, find the explanation, correct it, and complete the intended analysis | Error/Help tests exist; instruction clarity and recovery without losing valid settings remain to be observed |

Run the applicable task sequence in both Japanese and English:

1. Enter or restore the input and confirm its provenance, column mapping and
   settings. Record the source data/configuration identity.
2. Perform one explicit analysis. Record the fitted AnalysisID and its actual
   convergence, uncertainty and interpretation states.
3. Follow the first useful diagnostic, open its available explanation and
   return. Exercise Home, language and display switches where applicable;
   verify that these presentation actions preserve the fit and do not estimate.
4. Obtain the CSV and Word/PDF/HTML reports. Trace their values, basis and hold
   reasons to the same analysis. Allow documented rounding, display-scale
   transformations and intentional blank precision cells; do not require
   technical CSV values and held publication cells to be textually identical.
5. Ask the reviewer to explain what the result permits, what it does not permit,
   and the next useful action before giving corrective guidance.

Browser work must cover the ordinary and held journeys at desktop and narrow
widths, keyboard navigation and visible focus, and state preservation after
return. Reuse the existing Help/browser matrix for its exact cases rather than
inventing a second acceptance standard. AppTest output is supporting evidence,
not a substitute for rendering, focus or assistive-technology checks.

For an initial formative review, aim for three novice and two experienced
reviewers, subject to availability; this is a working recruitment target, not
a statistically representative sample or a pass threshold. Use the same
task wording and record unaided completion, assistance, confusing wording,
interpretation errors and the reviewer's explanation. Avoid collecting real
ratings or identifying research data. A belief that convergence proves
validity, posterior SD equals sampling SE, or a held result supports a final
conclusion is a critical misunderstanding requiring a correction and another
review. Record unobserved tasks as untested, not as zero errors. With no human
participants, the understanding gate remains open.

Measure J0 during the same scenarios: cold and repeated end-to-end time,
preparation/estimation/diagnostics/rendering/export time, peak memory, data
size, settings and failures. Compare the same requested work on the same
environment. Use the results to select one bottleneck before choosing a new
backend; the performance target must be fixed before evaluating its replacement.
The [current J0 protocol](../ROADMAP.md#performance-priority-20261014) specifies
input selection, paired repetitions, cache conditions and adoption criteria.
UX1/J0 are the next work pair; further TAM comparison scope is chosen after
that cycle unless a numerical or output-qualification defect needs priority.

Retain a bounded result record with case, language/view/device, exact source,
observed behavior, assistance, misunderstanding, timing, artifact checks,
unresolved issue and follow-up outcome. Keep implementation, browser, human
understanding, statistical qualification and deployment status separate.
An observed correctness/state/claim defect blocks acceptance of the affected
journey. Fix and retest it before extending the interface. A small formative
review can support the tested journey; it cannot establish general usability
or qualify the estimator's uncertainty.

### Local execution record — 2026-10-14

[`test_analysis_journeys.py`](../tests/test_analysis_journeys.py) passed seven
AppTest cases against the local working tree. The earlier table consistency
checks remain relevant to rounding and uncertainty meaning. The matrix above
is the original acceptance protocol; the following table records execution,
not completion of its human tasks.

| Executed case | Observed result | Remaining limit |
|---|---|---|
| Ordinary JMLE; extreme-score free-facet JMLE; fixed-SD EM MML; free-SD EM MML | Paste and map synthetic data, fit in Japanese, visit fit diagnostics and Help, switch to English, return through Home, build PDF/Word/HTML, switch back to Japanese, export CSV. AnalysisID and values stay fixed; presentation actions make zero estimation/refit calls. Publication holds and CSV row statuses persist. Public CSV omits Person rows; explicit private export retains them | Four synthetic inputs, not every model/design. Automated document tests omit figures; a separate free-SD timing run includes actual figure generation. No browser or unaided task observation |
| Changed settings | Changing the iteration limit hides current-result tables and their downloads; restoring the setting reveals the same fit without estimating. Warning and rerun button are checked in both languages | This exercises one setting change, not every input/control combination |
| Incorrect score-column mapping, Japanese and English | Selecting Person as Score produces an error and no fitted output. Correcting Score allows a real fit and preserves the iteration limit and language | One recoverable error, not all parsing, file-size or missing-data cases |
| Restored free-SD result with older positive uncertainty labels | An opt-in session-history snapshot is given older `ok` labels and restored through its button. Changed active settings receive the translated historical-result notice; the actual Word report still withholds precision. The saved snapshot is unchanged and no refit occurs | Synthetic legacy labels in current session history; not an import of every historical file format or version |

The run exposed a real reuse failure: EM returned a function-local class that
could not be serialized, so both report-local calculations and document-byte
reuse fell back to rebuilding. EM now uses SciPy's existing `OptimizeResult`
representation. Actual fixed/free-SD EM reuse, changed-input invalidation and
copy isolation are tested. Six before/after fits (RSM/PCM/GPCM × fixed/free SD)
retain exactly equal parameter values, objective/SD traces, termination states,
measure and reliability tables; elapsed time is intentionally excluded.
No inference hold or numerical estimation step was changed.

The stale/restored-result notices and rerun controls now use both locale files,
with a direct explanation of which analysis the displayed output belongs to.
They no longer insert raw setting identifiers into the main warning. The
broader history panel still contains English copy and an internal G2 label;
that remains a concrete language-review task, not a completed localization.

Validation: **127 distinct related tests passed** (49 report/reuse/locale
checks, 71 engine/quadrature/adapter checks and seven journeys). The environment
doctor, release check, syntax compilation and whitespace checks passed.
The complete native suite was not rerun. These are implementation checks, not
a resumption of formal G3/MML studies or evidence of coverage/qualification.

Environment: macOS 27.0.1 arm64; Python 3.12.6; NumPy 2.4.4; pandas 3.0.3;
SciPy 1.17.1; Streamlit 1.57.0; Matplotlib 3.10.9; ReportLab 4.5.1;
python-docx 1.2.0. Tested source SHA-256 values:

```text
streamlit_app.py  cc685cd3b3adf162402b7e0537a3441565cada407c248fdabf5640d416aed7d8
locales/en.json   0b879936d55bff574923dc931cf14be7d24d616d9ed708f56c958656447d18eb
locales/ja.json   b01886f8b25c26e92a939688efaf4c8d9c1753c5de40e984837535ee42397b1c
```

No in-app browser control connection was available in this session. Desktop
and narrow layouts, keyboard/focus behavior, assistive technology and reviewer
understanding remain untested for this candidate. There was no deployment.

### Local timing observation — 2026-10-14

One before/after AppTest process used the 12-row gradient fixture (three
persons, two raters, two tasks, scores 0–2), RSM, EM, free population SD,
15 quadrature points, 100 maximum iterations, tolerance `1e-6`, starting SD 1,
and the Standard analysis scope. Both followed the same Japanese fit → Help
and English return → Home return → PDF → repeat PDF → Word → PDF → files →
measures-CSV sequence. Actual document builders and figures ran without mocks;
public export was on. The processes ran without concurrent test jobs.

| Operation (seconds, including its AppTest rerun) | Before | After |
|---|---:|---:|
| Explicit fit and first result view | 2.682 | 2.529 |
| First PDF view | 21.517 | 3.651 |
| Same PDF view again | 2.596 | 0.712 |
| First Word view | 2.225 | 2.333 |
| Return to the already-built PDF | 2.338 | 0.723 |
| Open data files | 24.059 | 7.642 |
| Select the measures CSV | 23.247 | 6.859 |
| Sum of all 18 measured UI steps | 81.689 | 27.384 |

Process peak RSS was 605.5 MiB before and 513.8 MiB after (`ru_maxrss`, converted
from macOS bytes). This includes Python, AppTest and rendering libraries;
it is not the hosted server's memory budget or browser memory. Both processes
used `/tmp/mfrm_journey_mpl` and `/tmp/mfrm_journey_cache` for writable rendering
caches. Renderer/font initialization and persistent cache warming confound the
first-PDF and total-time comparison. Repeated views, together with builder-call
tests, establish that reuse now works; this single small-data observation does
not estimate a general speedup or prove memory improvement at production scale.
The first Word build did not become faster. Neither run raised an AppTest
exception. These local timings are not an immutable release benchmark.

J0 remains partial: phase timings, representative larger inputs and repeated
measurements are still needed. CSV selection still takes about seven seconds
on this tiny input. The next performance task is to separate frame collection
from ZIP/Excel/HTML preparation and rendering costs before changing that route.
Source inspection confirms that `_render_downloads` collects the complete frame
bundle and calls ZIP/Excel/HTML cache functions in the data-table panel, including
CSV-selection reruns. It does not establish cache misses or attribute the seven
seconds to any one builder. Time fingerprinting and cache behavior separately.
The next copy task is the existing history panel; recovery notices have already
been translated, but their meaning still needs to be tested with users. Browser
operation and human understanding remain independent acceptance requirements;
these findings do not justify a new backend or relaxing statistical holds.

### Local download and history refinement — 2026-10-15

The first profile identified repeated report-frame construction as the main
cost: on the 12-rating input, `collect_download_frames` took 14.03 of the
14.35 profiled seconds inside the downloads renderer. ZIP/Excel/HTML cache
calls on CSV reselection took less than 0.001 seconds combined; they were
already hits. Profiling overhead is excluded from the comparisons below.

The first candidate reused the manuscript claim guide and claim/evidence
matrix within the existing report-local scope. Its five-pair CSV median
reduction on 960 ratings was **27.5%, below the prespecified 30% target**.
All content checks passed, but this candidate was not accepted. Its record
remains at `validation/generated/download_j0_20261014/`, including the plan,
source, profile, timings and explicit failed acceptance decision.

The remaining profile computed interrater agreement 15 times for one export,
using 20.6 of 28.8 profiled seconds. The accepted candidate also reuses the
agreement summary within the same bounded scope. It uses complete argument
fingerprints, isolated copies and the existing eight-entry limit; the scope
expires when the export ends. Estimation and agreement formulas, thresholds,
publication holds, output formats and cross-session retention are unchanged.

The second comparison retained the original targets: at least 30% reduction
in the CSV median for each input, no more than 10% regression in the first
files-view median, process RSS below 1 GiB, and a 120-second step timeout.
Each input had five alternating baseline/candidate pairs in one process.
The baseline disables only the three new reuse wrappers on identical fitted
inputs. Each first files view clears the ZIP/Excel/HTML byte caches; CSV
selection and its repeat retain them. Values below are medians (minimum–maximum)
in seconds, including the AppTest rerun, with no profiler or parallel test job.

| Input | Operation | Baseline | Candidate |
|---|---|---:|---:|
| 12 ratings, free-SD EM RSM | Open files; export-byte caches cleared | 9.336 (8.442–9.535) | 4.297 (4.207–4.741) |
| 12 ratings, free-SD EM RSM | Select measures CSV | 8.412 (7.936–8.600) | 3.758 (3.660–3.983) |
| 12 ratings, free-SD EM RSM | Repeat the same CSV view | 7.967 (7.756–8.515) | 3.703 (3.564–3.968) |
| 960 ratings, JMLE RSM | Open files; export-byte caches cleared | 22.615 (22.228–23.874) | 6.849 (6.643–7.211) |
| 960 ratings, JMLE RSM | Select measures CSV | 21.180 (21.011–21.882) | 5.863 (5.624–5.929) |
| 960 ratings, JMLE RSM | Repeat the same CSV view | 21.229 (21.048–22.666) | 5.791 (5.712–5.943) |

All 20 content comparisons passed: every collected public table and its
attributes, ZIP member bytes, HTML bytes, and Excel sheet names/cell values
matched. Workbook creation timestamps are not compared. Public measures CSVs
excluded Person rows. After each three-step sequence, the fit and AnalysisID were unchanged;
the 60 measured UI steps made zero estimation/refit calls. Both inputs passed the
unchanged performance targets. Peak process RSS was 468.1 MiB,
a shared-process resource bound, not evidence of a before/after memory saving.

The history panel now translates headings, controls, column labels and restore
availability in Japanese and English. It explains optimizer-reported convergence
and the model-selection hold without an internal G2 label. Technical CSVs
retain qualification and reason codes; the display omits the descriptive
table's preferred-direction cue because it could suggest ranking held results.
Restore, unavailable snapshots, delete/cancel/confirm, and retained current
results are covered by bilingual AppTest checks.

Validation: **224 targeted tests passed**, including real input-to-export
journeys, report/uncertainty consistency, reuse invalidation and copy isolation,
output qualification, locales, history controls and the standalone boundary.
The new agreement reuse is checked after changing ratings, weights, fitted
parameters and the selected facet. Syntax and whitespace checks passed.
Final source fingerprints:

```text
streamlit_app.py  e91b1fc66c9ef59c7e360787c9219b570471d1bb4cd2df325f9f10e6100abb06
locales/en.json  1f84ec461eba6d36457f76b9ccf463c5d8b621eba1798314954c27d774f94064
locales/ja.json  ae7c9b52d05bc1cbea80d51cdb1ecc1d5ad543b41ce38f315c4b586296f60af2
```

Local evidence is retained at `validation/generated/download_j0_20261015/`;
the preceding failed candidate is preserved separately. The manifests record
sources, inputs, environments, plans, all timings and acceptance results.
These generated directories are not distributed in a fresh checkout.

Limits: two synthetic inputs and Japanese performance runs, not broad scale
acceptance. The initial process/import/font-cache cost and changed-input timing
were not compared. Browser latency, layout, keyboard/focus and actual-user
understanding remain untested; no browser-control connection was available.
The full native suite, formal G3/MML studies and deployment were not run.
Next, examine remaining waits under fresh-process and changed-input conditions,
then selected larger/sparse inputs; review the adjacent comparison-panel copy
and complete browser/human checks separately. Numerical accuracy and coverage
qualification remain independent of this performance result.

### Local comparison and fresh-process follow-up — 2026-10-16

The adjacent comparison view had two interpretation defects, so the roadmap's
correctness-first rule took precedence over copy-only changes. A high Pearson
correlation triggered claims of equivalent measures or defensible methods;
the reliability table also bypassed the hold applied by publication reports.
For example, `(0, 1, 2)` versus `(5, 6, 7)` has correlation 1 and RMSE 5. A
correlation threshold cannot establish equality or method validity.

The view now describes differences without success/error bands based on
correlation. It summarizes matched finite estimates separately for each facet,
shows each side's element count and the comparable count, and leaves undefined
correlations blank. Duplicate or missing element identifiers prevent matching.
The scatterplot and difference table use the same finite pairs. Fixed values,
scales, weights and estimation targets still need review; this view does not
align coordinates or qualify method agreement.

Comparison reliability uses the existing publication hold logic. Held precision
cells are blank, with translated interpretation scope and a warning for each
held analysis; saved technical values remain intact. Model-selection holds
remain visible without internal stage labels. Missing data/source records and
one-sided additions of fixed-value inputs no longer count as matching inputs.
A matching record is not a statistical qualification.

Controls, column labels, plot labels, input warnings and interpretation notes
are translated into Japanese and English. Numbered history options distinguish
entries with identical timestamps and summary labels. The initial selector
checks exposed this ambiguity (2 failed, 25 passed); after correction,
**93 targeted tests passed**. These include four real fitted cases,
bilingual comparison and history controls, input-to-export journeys, uncertainty
and information-criterion holds, exact preservation of saved values, and zero
estimation/refit calls during comparison. The offset example, separate facets,
constant values, too few pairs, non-finite values, unmatched and duplicate IDs,
and missing or changed source records have direct regression checks.

Performance follow-up used the existing 12-rating free-SD EM RSM fixture and a
predeclared one-score change, five alternating pairs, and a fresh Python process
for each variant. Each process ran the real app, fitted the first input, exported,
changed the input, confirmed stale exports were hidden, explicitly refitted and
exported the new result. The baseline disables only the three report-reuse
wrappers added on October 15; comparison-view corrections are common to both.
No additional performance implementation was made in this follow-up.

The plan preceded measurement: CSV medians must improve by at least 30% for both
inputs, files-view and total measured-work medians must not regress by over 10%,
RSS must stay below 1 GiB, each step below 120 seconds and each worker below
300 seconds. All targets passed. Values are seconds, median (minimum–maximum):

| Stage | Reuse disabled | Reuse enabled |
|---|---:|---:|
| Imports in fresh Python process | 1.901 (1.664–3.588) | 1.762 (1.606–1.960) |
| First page render | 0.255 (0.222–0.283) | 0.260 (0.208–0.286) |
| Input/setup reruns | 1.364 (1.326–1.393) | 1.355 (1.319–1.389) |
| Initial fit and result render | 2.562 (2.457–2.619) | 2.491 (2.428–2.606) |
| Initial files view | 7.378 (6.914–7.501) | 3.635 (3.629–3.776) |
| Initial CSV selection | 6.712 (6.536–7.074) | 3.096 (3.016–3.275) |
| Changed input; stale exports hidden | 0.272 (0.253–0.275) | 0.266 (0.256–0.278) |
| Changed-input fit and result render | 3.282 (3.073–3.339) | 3.247 (3.154–3.476) |
| Changed-input files view | 7.493 (6.762–8.002) | 3.630 (3.459–3.995) |
| Changed-input CSV selection | 6.622 (6.218–7.016) | 3.211 (3.067–3.219) |
| Sum of measured work | 37.483 (36.371–39.857) | 23.018 (22.682–23.558) |

Peak worker RSS was 320.6 MiB, including correctness checks. It is a resource
bound, not an isolated memory-saving estimate. All 20 fitted phases preserved
public tables/attributes, ZIP members, HTML, and Excel sheet names/cells against
an uncached reference on the same fit. Workbook creation timestamps are excluded.
CSV excludes Person rows and preserves the interpretation hold. Input changes
produce a new AnalysisID; presentation preserves the current and saved values
and performs no estimation/refit. Two explicit fits per process are development
checks, not a resumption of formal research.

The total is a sum of measured computation, including imports, setup reruns,
fitting/result display and both exports. It excludes correctness checking,
human think time and process launch overhead. Fitting and result rendering are
timed together, so these results do not identify estimator-kernel speedup.
Python/Streamlit caches start fresh, but OS and persistent font caches are not
flushed. After the first export, the uncached correctness check repopulates
export caches before the changed-input phase; this is common to both variants.
The retained evidence is local to `validation/generated/comparison_ui_j0_20261016/`:
plan, runner, worker logs, fitted inputs/outputs, all timings, initial and final
test results, source snapshots and SHA-256 manifest. It is not distributed in a
fresh checkout. Source scope checking found changes only in the five comparison
helpers/renderers; estimator and interrater-agreement kernels are unchanged.

```text
streamlit_app.py  2c3524365391a91a4ae6642e0810195192c6b1e42811346b352e3235318222db
locales/en.json  ca487e0539beb1470c6a62242dbd766ce3051bc57bcdfb690f2f24ed53d9ed8c
locales/ja.json  d66eb3bc0de4c9ac3d1de608bedb89de7b39b1453c6805af0197856a1f48f8e5
tests/test_comparison_ui.py  b00de30ae9ff68bbc255633fa0fe755054a7e9c608a749cd22e00f3fc6faa1ef
```

Limits and next decision: one small input and one modification in Japanese do
not cover 960-rating fresh runs, larger/sparse inputs, or other model/weight/
anchor settings. Browser control remains unavailable in this session, so
layout, keyboard/focus and actual-user comprehension remain untested. The full
native suite, formal G3/MML studies and deployment were not run. Next select a
bounded larger/sparse input under the existing J0 criteria, and complete browser/
human checks when available. Avoid widening the model or backend scope solely
because these implementation checks passed.

### Help readability and internal-copy review — 2026-10-16

Scope: the tutorial, all eight full Help panels, the glossary, five chart
popovers, and model-comparison/bias hold explanations. The review found a
developer instruction in the glossary, raw state and validation-stage names
in guidance, and two conflicting estimation explanations shown together.
The workflow panel rendered 25 tables, including repeated advice.

Changes:

- Use the current localized estimation explanation once. Explain JMLE/MML,
  fixed or estimated population SD, EAP uncertainty, excluded rows and input
  errors against the current implementation.
- Put practical steps first. Keep the workflow's nine remaining reference
  tables in seven collapsed disclosures and remove repeated specialist
  tables from the interpretation and rating-scale panels.
- Rewrite the interpretation, rater, category, reporting and troubleshooting
  guides in both languages. Correct the reviewed glossary and chart claims:
  fit bands are review prompts, correlation is not agreement, category
  findings do not automatically justify merging, and an uncertainty label
  does not remove the corresponding inference hold.
- Explain comparison and bias restrictions without development-stage names.
  Keep technical reason codes and withheld numerical quantities in the
  model-comparison CSV. A cached recommendation is still overridden by the
  current restriction. Export filenames remain visible where users need them
  to locate evidence.

The copy revision changed 60 Japanese and 68 English locale entries. Its
application function changes are confined to eight tutorial/Help/rendering
functions; five chart-popover entries were also updated in each language.
Numerical kernels and interpretation qualification rules are unchanged by
this copy revision. Earlier uncommitted numerical/report changes are included
in the integration checkpoint and have their own records above.

Automated rendering checks cover all eight panels in Japanese and English,
one estimation tutorial, collapsed workflow details, absence of the reviewed
internal vocabulary and zero estimation/refit calls. A separate bilingual
test verifies that cached model guidance stays held, public tables hide
reason codes and p-values, and the CSV retains its technical record. The
first focused run found one remaining English implementation term and two
assertions tied to obsolete guide wording; those were corrected, and the
six affected/new cases passed on rerun.

Pre-push verification completed with Python 3.12.6:

- `make verify PYTHON=.venv/bin/python` passed compilation, environment and
  release checks, all 56 self-tests, the native test selection, the benchmark
  smoke check and demo-report generation.
- The native suite reported **3,110 passed, 2 skipped, 84 deselected** in
  1,062.36 seconds. The two optional CMLE comparisons require the unavailable
  R package `immer`; the normal selection excludes retained-study and legacy
  tests and the documented external bundles. These are implementation checks,
  not a restart of formal validation studies.
- Local documentation link targets and `git diff --check` passed. Generated
  reports, source snapshots and logs remain outside the commit.

This local run does not verify Python 3.11. The current GitHub push filter does
not include `integration/unified-app`, so pushing this branch alone does not
run the two-version CI matrix.

Source scope and hashes are retained locally in
`validation/generated/help_copy_20261016/final_change_scope.json`, with the
pre-review source under `before/`. These ignored artifacts are not included
in a fresh checkout. The earlier timing measurements apply to their recorded
source hashes; this copy review does not repeat or extend those performance
claims. Browser layout, keyboard/focus behavior and actual-user comprehension
remain untested. Formal G3 and later MML studies remain paused. Pushing the
integration branch does not publish a new beta or verify the hosted app.

## Initial audit — 2026-08-04

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

## Residual PCA and bias navigation revision — 2026-09-13

The diagnostic selector previously placed Wright Map and Visuals before bias,
and labelled residual PCA only as Dimensionality. More seriously, recommended
actions changed the top-level section without selecting the named diagnostic.
A PCA recommendation could therefore open fit details or the last-used view.
The PCA page also preferred the first rater facet, although First Read reports
the overall residual matrix.

The existing selector now starts with model fit, residual PCA, bias, and
categories. Recommended actions select their diagnostic directly; a PCA action
opens the overall scope used by First Read. Explicit facet selection remains
available. Stable panel IDs and fitted-result reuse are preserved. Programmatic
navigation no longer competes with widget defaults in Session State.

PCA starts with a short scope explanation and the plot. The reading guide and
optional DIMTEST controls are closed disclosures below the existing evidence.
The bilingual guide distinguishes residual variance from total variance and
avoids treating fixed eigenvalue cutoffs as proof of unidimensionality; see the
[Winsteps residual-dimensionality discussion](https://www.winsteps.com/winman/dimensionality.htm).
Overall and facet PCA use different residual aggregation and are explicitly
distinguished. Bias starts with the pair choice and screening results; its
technical settings explanation is inside the existing settings disclosure.
Heatmap stars are described as an unadjusted |t| >= 2 screening flag, including
the English contextual help, rather than multiplicity-adjusted significance.
Estimation, screening thresholds, stability evidence, and inference holds are
unchanged.

Validation: **153 tests passed**, covering direct PCA/bias/category routing in
both locales, overall PCA selection after a stale rater selection, fit
preservation, lazy diagnostic rendering, locale parity, help navigation, and the existing PCA stability and
bias inference/publication suites. A local Chrome check of the built-in
960-observation RSM/JMLE example retained fingerprint `067a6d8d`; its recommended
PCA destination showed the same first eigenvalue (3.74) as First Read and kept
the stability warning. The Person x Rater bias view retained 120 cells, 24
flagged cells, zero strong flags, and the pairwise-inference hold. Focused
desktop and 400-CSS-pixel checks, including Japanese PCA and bias views, supplement
AppTest; they are not acceptance of the full accessibility protocol or evidence
from first-time users.

## APA manuscript entry — 2026-09-13

Added **Start a paper** as a distinct Report & Export task. The previous
template entry was hidden in review resources and displayed a long technical
worksheet. The new entry offers an English APA Word or Markdown scaffold,
with recorded analysis facts filled in and research-context prompts left for
the author. Preview and method-literature candidates are closed disclosures.
Results prose is opt-in and reuses the existing guarded draft, including
simulation context. The top-level report holds remain visible on this route.

The template follows the user's question-to-design-to-answer editorial
priorities. Sources, author workflow, and boundaries are documented in
[the manuscript guide](apa_manuscript_template.md). The Zotero review informs
the prompts; it does not install live citation fields or add uncited references
to the article. The runtime needs no new dependency or Zotero connection.

Validation: **91 tests passed**, covering the new route and document structure,
both UI locales, opt-in draft generation, simulation-context forwarding,
existing report holds, reference integrity, and the validation contract.
The reusable Word template was rendered with the bundled LibreOffice runtime;
all seven pages were inspected, including headings, page numbers, indentation,
and the absence of inherited title decoration.
The built-in 960-rating RSM/JMLE run retained fingerprint `067a6d8d` while
opening the writing task. The Word download was inspected and contained the
recorded 960 ratings, 30 person identifiers, RSM/JMLE, and version 0.2.15-beta;
it contained no person-level rows. Desktop and English/Japanese 400-CSS-pixel
browser views were checked without horizontal page overflow. These are focused layout and interaction checks, not full browser
accessibility acceptance or a first-time-user study.

## Focused results instead of stacked guidance — 2026-09-13

Reference: [langtest.jp MFRM](https://langtest.jp/shiny/mfrm/), inspected in a
browser on 2026-09-13. Its useful pattern is direct selection of a named analysis
view, with a focused input area and one Run action. This is a UI reference;
its PCA cutoffs and bias-significance wording are not adopted as statistical
policy or treated as external validation of this application's estimators.

The previous Compact result shell placed interpretation status, a goal chooser,
an action hub, run history, global ZIP/Excel exports, and a second reading-order
brief ahead of the selected result. The shell now uses one section selector,
one run-specific interpretation notice, and one named next-check button. The
first result view opens First Read. The claim boundary remains visible; complete
check evidence and safe-output guidance are available in a closed disclosure.
The same existing priority and diagnostic-target functions determine the route.
Missing guide evidence displays a warning rather than interpretation clearance.

Removed the duplicate generic success banner and section-reading brief. Reading
guides and data-design details start closed. Run history and comparison live in
Start; Compact file exports live in Report & Export. Full display retains its
quick bundle entry. No estimation algorithms, diagnostic thresholds, inference
holds, or export privacy filters were changed. No package was added.

Validation: **214 tests passed**, including both locales, integrated section and
next-check navigation, unchanged fitted-state checks, unavailable guide evidence,
PCA stability, bias/report holds, data readiness, and the validation contract.
The browser sample retained RSM/JMLE, 960 ratings, 30 persons, 30 iterations,
fingerprint `067a6d8d`, first PCA eigenvalue 3.74, and the stability warning.
At 1365 CSS pixels the section selector and next-check button fit in the initial
viewport; previously the selector followed repeated guidance and global download controls. The
Japanese 400-CSS-pixel view had no horizontal page overflow, and Enter on the
recommended action opened overall PCA. Report & Export remained reachable and
its public-mode ZIP contained 95 tables plus the two evidence-contract JSON
assets, with a clean ZIP integrity check. Help's screen-order table was updated.
An additional Help/workspace check run passed all 35 tests (overlapping the
workspace checks above).

These are focused layout and interaction checks, not first-time-user acceptance
or completion of the full accessibility protocol. The setup sidebar still has
many controls, and some run-specific diagnostic text is English in Japanese UI.
Those are the next concrete simplification/localization targets; adding another
onboarding layer would repeat the problem addressed here.

## Sidebar: visible selections, optional editors — 2026-09-13

The sidebar now uses one source selector plus the existing scenario selector.
The sample name and a compact size/category summary stay visible; the long
scenario explanation, references and CSV download share a closed disclosure.
The sample CSV download does not rerun the app. Display density is also a closed
disclosure, while language and Help remain directly available.

Column mapping is shown as an always-visible summary of the actual Person,
Score and facet selections. Its editor starts closed for registered built-in
samples and open for other sources, including pasted and uploaded data. All
mapping widgets are still instantiated with their existing keys; closing the
editor does not remove them from Streamlit session state. An ignored weight
column and the insufficient-facet warning remain outside the editor. Model,
estimator, analysis depth and selected bias pair remain directly adjustable.

Merged the explanatory compute-plan captions and detailed settings into the
existing setup disclosure, and removed the repeated sidebar row-count banner
and Run-location caption. The main sample banner now names the example in the
selected UI language without repeating the full design. README's smoke-run
instructions were updated to the current source, Run and result controls.

Browser inspection found that Streamlit could retain an old selected-option
label after a language change even when its dropdown options had translated.
The source and depth selectors now resend their canonical selection, refreshing
the label. Standard remains the initial analysis depth, and returning from an
Advanced-only Custom plan to Guided retains the existing Standard fallback.
The pre-run depth description is localized independently of its analysis ID.

Validation: **59 tests passed**, covering source-state migration, source/scenario
round trips, both mapping-disclosure states in both languages, manual mapping,
visible facet warnings, weight handling, Advanced controls, locale label updates,
Help/sample-guide state, small-data widgets and the config whitelist. The real
Guided-versus-Advanced default run test compared person, facet and step estimates
with zero relative tolerance and absolute tolerance 1e-12. The browser sample
retained 960 ratings, RSM/JMLE, fingerprint `067a6d8d`, and PCA eigenvalue 3.74.
A separate pasted synthetic dataset opened the mapping editor automatically.
English and Japanese 400-CSS-pixel views had no horizontal page overflow; the
mapping editor remained operable in the sidebar. After the localized sample
banner and README update, the overlapping app-smoke/locale subset also passed
all 22 tests.

These changes concern presentation and state synchronization, not estimation
algorithms or statistical thresholds. The sidebar still scrolls: detailed
mapping for user data and necessary model choices remain available. This is
focused verification, not completion of the full accessibility protocol or a
first-time-user study; the previously recorded language fast-path/exact-return
acceptance work remains separate.

## Next UI refinement — 2026-09-13 (planned)

Reference: KAWAI's [September 12 post on 20 UI psychology perspectives](https://x.com/kawai_design/status/2098698136600690965).
The ordinary web reader returned 403; the author's excerpt was checked through
X's public oEmbed API. The full post/media list was then obtained through a
public mirror and all four original `pbs.twimg.com` images were inspected.
This is a design reference, not primary evidence that a psychological effect
will improve this application. Its four groups concern perceptual grouping,
memory support, attention, and interaction. The images themselves caution
against treating chunking as a universal seven-item rule, applying choice-time
laws to every complex judgment, or equating attractive design with successful
use. No reference images are copied into the product or repository.

The application already has named result routes, one next-check action,
task-based exports, visible column selections, and an optional sample guide.
The next pass refines those surfaces; it does not add twenty new components,
another onboarding system, or a separate front end. Code inspection of
`run_facets_mode()`, `_render_guided_report_export_section()`, and the retained
Stan/Viewer functions supplies the current-state basis below. These are planned
changes and checks, not a claim that the new UX has been implemented or tested.

| Priority / existing surface | Concrete refinement | Acceptance evidence |
|---|---|---|
| U0: Run and current-result state | Keep the selected data/model/method beside Run. Distinguish pending settings from the settings of the displayed fit. Put completion, error or refit-needed feedback near the triggering action, with a meaningful focus destination. | Change settings after fitting, open Help, change language and return: the old result is never presented as a newly fitted result. Display/navigation/download changes do not silently reestimate or alter AnalysisID. |
| U1: Input checks | Group the source, column summary and required corrections by task. Consolidate repeated successful counts/readiness tables into a short summary with detail on request; retain exclusions, weighting problems and blocked states beside the relevant control. | A person using their own data can identify the field to fix and the rows excluded from the likelihood. Simplification does not remove a refusal, change retained rows or hide a necessary data-handling decision. |
| U2: PCA and bias | Retain stable named routes and the existing next-check priority. Put the result, interpretation limit and specific action together; use labels such as “Inspect residual structure (PCA)” or “Check rater–task interaction”, localized consistently. | Users can locate PCA and the selected bias pair without recalling an internal tab name. Neither a PCA value nor a screening flag is presented as proof of unidimensionality, unfairness or causation. No cutoff changes. |
| U3: Report & Export | Preserve the existing document/manuscript/files/review tasks. Use content-specific download labels and show the artifact's purpose, format and data scope near its action. Place future offline Bayesian code under Files as an optional task, followed by engine selection. | The first export screen does not become an engine/language catalog. A user distinguishes a manuscript draft, a reproducibility package, a model template and actual posterior results. A template download is never labelled analysis completion. |
| U4: Remembering the current analysis | Keep a compact, consistent analysis-context summary and stable navigation order. Use the existing history only where comparison is requested. Preserve current selections across language, density and contextual Help changes. | The current fit, facet pair and result section are recoverable without reconstructing choices from memory. No new duplicate identity/state store is introduced. |
| U5: Visual hierarchy and language | Standardize spacing, alignment, heading levels and primary/secondary button roles through native Streamlit components and the existing theme. Keep one primary action per task. Finish the already recorded Japanese gaps in dynamic diagnostic text. | Meaning remains visible without color. Literal translations do not change the statistical target. Avoid warning walls, nested cards, decorative pictures and pervasive confirmation dialogs. |
| U6: Reach and focus | Keep primary actions near their objects, with adequate target size and spacing; preserve visible labels, keyboard operation, focus order and readable reflow. Reuse the existing accessibility matrix. | Check English/Japanese at narrow and desktop widths, keyboard-only use and zoom; retained warnings and selected routes remain reachable. Visual polish is not accessibility certification. |

The connection between stages should show the actual workflow: prepare → run →
inspect → report. A progress indicator must not mark interpretation or
publication complete merely because estimation ended. Optional help stays
skippable and recoverable. Significant result changes are announced at the
action/result location; animation is not required to communicate them.

For the future Bayesian route, the sequence is “prepare code → run locally →
review returned results”, with a clear indication when only code preparation is
available. A completed local run, acceptable MCMC diagnostics and a defensible
substantive claim are separate states. The sampler is not added to the main
JMLE/MML selector. See the [B0–B4 model and output gates](jax_numpyro_conquest_tam_plan.md)
before exposing either downloads or posterior uploads.

### First implementation and evaluation slice

Start with U0/U1 and the existing English/Japanese sample and pasted-data paths.
Measure the current version first: time to find Run, time to find the requested
diagnostic, wrong turns, assistance, avoidable reruns and export generation
time. Use the same synthetic inputs and retained fitting conditions. Separate
calculation time from navigation and rendering time; JAX cannot fix all three.
Then apply U2/U3 and U5/U6 to those same journeys, one surface at a time.

Use a small formative review with both first-time and experienced MFRM users;
its purpose is to find breakdowns, not estimate population-wide usability.
Tasks: prepare a sample/own-data analysis, find PCA, inspect a named bias pair,
obtain an English manuscript draft, and distinguish a future code-only export
from a fitted Bayesian result. Record unaided completion and interpretation
errors, not just satisfaction or click counts. Preserve critical controls for
missing evidence, unsupported inference, pending settings, incomplete external
runs and mismatched result files. Numerical fixtures and state guards verify
that easier navigation leaves the fitted results and their restrictions intact.
Do not launch background telemetry or a large new study framework for this review.

## License, software citation and Home — 2026-09-13 (local implementation)

The follow-up request adds a small part of U3/U4. The existing MIT license is
unchanged. A collapsed bilingual sidebar entry now contains the license link,
a copyable APA software reference and BibTeX/CFF downloads. The manuscript
template's references guidance reuses the same renderer. The user supplied
the author name **Ryuya Komuro**; the reference uses **Komuro, R.**
`CITATION.cff` is the canonical metadata, written in JSON syntax (valid YAML 1.2)
so the runtime can read it with the standard library. Its official CFF 1.2.0
schema validation passed. Citation remains a scholarly request, not an added
MIT restriction, and method citations and recorded fit versions remain distinct.
See [GitHub's citation support](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files).

Home and Return to my analysis change the displayed workspace without clearing
the input or fit. They reuse the Help path that keeps input/upload/settings
widgets alive and returns before estimation. Existing result-route selections
are retained; queued one-shot Run requests are cancelled by Home navigation.
No second fit store or background job is added. This preserves a browser
session; it is not a backup across restarts. AppTest covers navigation and
English/Japanese state changes; manual browser upload, focus and narrow-width
acceptance remain part of the existing browser matrix. This local change has
not been committed, pushed or deployed.

Local verification: the app-smoke, Help-adapter, locale-parity, workspace and
export suite passed **56 tests**. The contextual-Help, standalone-boundary and
repeated citation/locale smoke checks passed **57 tests**. The latter reuses a
real sample fit and confirms that Home → language change → return preserves
the exact output object, AnalysisID, selected result panel and fit-view values,
with zero additional estimator/refit calls. Compile and whitespace checks also
passed. Existing estimator function bodies were unchanged by AST comparison.

### Result-bound citation exports — 2026-09-13 (local follow-up)

The next slice carries the same metadata into the existing manuscript
Markdown/Word template, publication Word/PDF/HTML references, and standard
result ZIPs and manuscript binder. Shared export sidecars add
`software_citation.md`, `mfrm_software.bib` and `CITATION.cff` only when the saved
`config.app_version` matches the verified citation version. Missing or older
versions receive an explanatory Markdown note, without current-release
BibTeX/CFF files. Report version fields no longer fill missing values from the
running app. The manuscript UI applies the same check; the sidebar still offers
the current app's reference. This adds no top-level navigation or runtime dependency.

The related manuscript/export/evidence/privacy/MML/i18n suite passed **78 tests**.
A follow-up citation/publication-figure/standalone suite passed **49 tests**,
including repeated citation cases and new English/Japanese download-state checks.
These verify version matching, omission of raw rows, retention of inference
holds, reference formatting and disappearance of stale downloads after a result
change. LibreOffice rendering was inspected on all 7 manuscript-template pages
and all 10 publication-document pages; the direct PDF's reference page was also
inspected, and PDF extraction confirmed the software title's italic formatting.
The software reference fits within the page with hanging indentation.
These used a synthetic configuration fixture, with publication figures omitted;
they establish export behavior, not numerical validity or submission readiness.

The legacy full publication report still repeats drafting guidance and exposes
some Markdown emphasis markers as literal text. Its general reference styling
also needs a separate APA review; adding the software reference does not qualify
the entire report. Keep the concise manuscript template as the drafting entry
point and address this legacy-report cleanup under U3. This follow-up remains
uncommitted and unpublished.

## 日本語の見直し — 2026-09-15（ローカル実装）

利用者から日本語の不自然さを指摘されたため、主要画面の翻訳672項目と、画面に直接書かれていた件数表示を修正した。
ホーム、サンプルガイド、列の指定、分析設定、入力確認、推定結果、残差PCA、バイアス、報告・保存の案内を対象とした。
「返却点を再確認」は「推定結果を確認」、「解析深度」は「分析する内容」、「Response データ監査」は
「分析に使う回答の確認」へ変更した。指示文はです・ます調にそろえ、必要な専門用語には意味を補った。

単に語を置き換えるのではなく、何を確認し、次に何をすればよいかが分かる文章にした。
「勾配が小さい」「計算が終了した」「積分が十分に正確」を区別し、勾配の対象は総NLLと明記した。
「31点以上を報告用に使う」という古い求積の案内は、点数を変えた推定値・尤度・得点の安定性を確認する説明へ改めた。
バイアス確認の未補正の目安や、帰無仮説を棄却しなかった結果も、問題がないと断定する表現にしない。
推定式・判定条件・既定値・翻訳キー・変数の差込欄・データ列名は変更していない。
既存の数値関数が未変更であることはAST比較で確認し、変更されたPythonの関数は2つの画面表示関数に限られる。

主要UIの139テストが通過した。ヘルプ・母集団SD・比較出力の39テストも確認し、
古い日本語の単語だけを要求していた1件を修正して比較出力3件を再実行した。重複を除く確認対象は178件である。
最後の説明文の調整後にも、翻訳キーと差込欄を含む言語テスト9件が通過した。
ローカルブラウザで日本語の開始・準備画面、390px幅のホーム、固定SD PCMの計算確認を表示し、
横方向のはみ出しがないことを確認した。翻訳を更新する途中のプレビューにはキャッシュが残ったため、
再起動・再読込み後の表示を最終記録として保存した。

数値の記録やソフトウェアの識別子など、保存した原文を読む箇所には英語が残る。
全ヘルプ・全診断原文の翻訳完了や、アプリ全体の新たな受入完了を意味するものではない。
過去の研究のソースと結果はそのまま保存し、公開版への反映は行っていない。
[修正・検証記録](../validation/generated/japanese_ui_polish_20260915/qa.json)と
[変更前のソース](../validation/source_snapshots/pre_japanese_ui_polish_20260915/sources.zip)を参照。

## 求積の読み方とオフライン確認 — 2026-09-15（ローカル実装）

固定SDの数値詳細に「読み方と次の確認」を追加した。説明は操作したときだけ開き、
通常表示には、停止状態と閉じた数値詳細を維持する。計算操作と説明を開く操作を分け、
説明を読んでも再推定・求積評価・外部エンジンの起動を行わない。

案内は「表の読み方」「さらに確かめるには」「手元で確認する準備」の三つに絞った。
EAP差だけでなく事後SD差も見ること、基準行が0でも精度確認の証明にならないこと、
列ごとの最大差が同じ人とは限らないことを説明する。
同一点の積分差、再推定したことによる変化、連続積分との比較を区別し、
画面のJSONと再実行用入力の保存先を分けて示す。

比較用ZIPには、表示言語に対応する `README_quadrature_check.md` を同梱する。
内容は[日本語手順](fixed_sd_quadrature_check.ja.md)／[英語手順](fixed_sd_quadrature_check.en.md)に置く。
Pythonの開発用確認は対応するアプリのソース一式が別途必要で、ZIPだけで実行できるとは案内しない。
対応外の入力を都合よく変更しないこと、総NLL・個人得点の誤差評価や再推定候補の解釈範囲を記載した。
独立した出力メニューや新しい分析エンジンは追加していない。

確認した範囲：

- UI・言語の13テスト、比較用出力・限定PCMランナーの5テストが通過（計18）。
  表示だけでは求積評価ゼロ、操作後は3次数を評価、言語切替で追加評価なし、
  元の分析の保持、未対応・不整合入力の扱い、日英の手順書が正しいZIPへ入ることを確認した。
- 最初のUIテストでは、テスト側がpopoverのprotobufラベルの階層を誤って参照して2件失敗した。
  実際の構造に合わせて参照を修正し、上記13テストを再実行した。画面実装の失敗とは分けて記録する。
- 8人・64回答の合成PCMを実際に推定して両言語のZIPを出力し、展開した入力から記載のCLIを実行した。
  予定した6経路すべてを評価でき、元の入力・分析は不変だった。これは手順・入出力の確認で、新しい統計的資格ではない。
- 変更前とのAST照合で、変更した既存関数は数値詳細の表示と比較用ZIP作成の2つに限定された。
  尤度・勾配・推定・個人得点の関数本体は未変更。
- ローカルChromeで日英の説明、比較後の言語切替、Tab・Enterによる説明の展開、Escapeによる閉鎖と
  案内ボタンへのフォーカス復帰を点検した。狭いドック表示でも説明が折り返され、内部をスクロールして最後まで読めた。
  幅を確定した390pxの受入や、アプリ全体のモバイル・スクリーンリーダー受入は実施していない。
  確認用タブを閉じ、ローカルサーバーを停止した。

[検証記録](../validation/generated/quadrature_guidance_20260915/qa.json)、
[実際のZIPからの実行確認](../validation/generated/quadrature_guidance_20260915/offline_guide_check.json)、
[変更前ソース](../validation/source_snapshots/pre_quadrature_guidance_20260915/sources.zip)を保存した。
公開・既定の推定方法・推論の資格は変更していない。次は通常の結果画面から保存までの通し確認を行い、
利用者が三つの比較を使い分けられるかは別途確認する。

## 通常画面から求積確認と比較資料への往復 — 2026-09-15（ローカル実装）

通常の `main()` から、64行・8人の合成PCMデータを貼り付け、固定SDのMML推定を実行した。
受験者データの分析結果ではなく、操作・保存の検証例である。Direct、GH31、SD=1.1、Criterionを
ステップと非中心化の対象、Fast previewとし、AppTestではmaxit=40、ブラウザーでは400を使用した。
ブラウザーの返却点は8反復、総NLL=98.67433024092199、最大絶対勾配=2.8936931435860644e-7だった。
この値から最尤解や積分精度を認定していない。

実際に直した問題：

- 比較資料を選んだ後に「結果」を開くと、元の出力先が初期化された。既存のナビゲーション保持処理を
  ホームだけでなく通常のセクション切替にも適用した。
- 言語切替専用の中間画面が入力ウィジェットを省略し、次の描画で入力・設定が失われた。
  中間画面と再表示ボタンを削除し、現在の画面を選択した言語で描画する。
- 上の修正だけでは、実ブラウザーの次の操作で「詳細設定」が「標準設定」に戻った。
  ラジオボタンは翻訳済みラベルを送受信するが、言語切替後に現在の選択を送り直していなかった。
  言語変更のコールバックで既存の文字列値を再設定し、新しい表示名をブラウザーへ送る。
  ボタンやアップロードファイルの状態を一括で書き換える処理は加えていない。
- 外部ソフトの得点を表す `native scores` が「標準得点」と誤訳されていた。
  共通の質問文を「各ソフトの個人得点は安定しているか」に直し、画面・HTML・表の表示を揃えた。

検証の範囲：

- `test_normal_workspace_preserves_fit_and_comparison_route_across_navigation` は実際のアプリ入口・推定処理・
  ダウンロードに渡されるデータを使用する。日英の比較用ZIPの入力とソース上の手順を照合し、
  画面往復・言語切替・ホーム復帰で元の分析オブジェクト、入力の保存内容、得点表、求積確認JSONが
  不変であることを確認した。ナビゲーション中の再推定と外部エンジン起動は失敗させる。
  実際に入力を編集した場合は、古い結果・比較用ダウンロードが非表示となることも確認した。
- 最初の自動テストはナビゲーションのリセットを検出した。ブラウザーで発見したラジオボタンの問題は、
  サーバー上の選択値だけを見たテストでは見逃していた。送信メッセージの `set_value` と翻訳済みの
  `raw_value` を検査するようにし、修正前の失敗と修正後の通過を記録した。
- ホーム・ヘルプ・初期画面・日英表示を含む142テストと、比較レポート・出力の13テストが通過した。
  先行する通過記録も残しているが、最終コードに対する142件は `regression_radio_final.log` を参照する。
- ローカルChromeで通常の推定、求積確認、JSON保存、比較資料への移動を実行した。
  修正前の英語ZIPと修正後の日本語ZIPは、ともに元のJSONに記録された入力ハッシュと一致する。
  データ・設定・推定値は言語と表示修正の前後で同一だった。日本語ZIPには日本語手順と修正した質問文を含む。
  数値比較の4項目は未評価のままで、推論を許可するフラグはfalseだった。
- ネイティブ操作では古い要素ID・画面取得エラーが発生し、途中から同じタブをブラウザー拡張経由で操作した。
  これらをアプリの計算障害として扱っていない。最終のブラウザー確認範囲は下記QA記録に保存する。

[QAと保存物のハッシュ](../validation/generated/quadrature_navigation_20260915/qa.json)、
[最終の回帰テスト](../validation/generated/quadrature_navigation_20260915/regression_radio_final.log)、
[変更前ソース](../validation/source_snapshots/pre_quadrature_navigation_20260915/sources.zip)を保存した。
初見の利用者による理解、スクリーンリーダー、幅を固定したモバイル表示の受入は未実施。
優先度・信頼性の案内に英語の動的メッセージが残っており、次のUX1ではその説明とPCA・バイアス・原稿の導線を点検する。
この変更は未公開であり、数学的な検証や統計的推論の資格を追加するものではない。

## Long-term implementation sequence

1. Apply the bounded U0/U1 slice above and execute the versioned browser matrix for keyboard/focus acceptance and
   exact contextual Help return, retaining evidence for every required row.
2. Complete browser and first-time-user acceptance for the implemented
   two-level source chooser and compact Guided setup while preserving stable
   IDs and old session state.
3. Verify the simplified sidebar with first-time users and confirm that
   contextual Help and language switching return to the exact selected view.
4. Extract source/setup/result-shell renderers from the monolith, one tested
   vertical slice at a time.
5. Add task-level accessibility and rendered-locale acceptance tests.
6. Start with a local task-review record of source class, phase, blocked reason,
   reruns and time to interpretable evidence. Add product instrumentation only
   if those observations leave a concrete measurement gap. Never record uploaded
   values, person identifiers, free text or raw file names in UX records.

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
