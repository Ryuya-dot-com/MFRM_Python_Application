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
