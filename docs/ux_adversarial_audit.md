# Adversarial UX Audit

Current checkpoint, 2026-10-23: the [layout and copy rewrite](#layout-copy-rewrite-20261023)
groups actions and removes repeated instructions from the main workflow.
Privacy, inference holds and meaningful control labels remain visible where needed.

<a id="layout-copy-rewrite-20261023"></a>

## Layout and copy rewrite — 2026-10-23

The requested change is to let placement and control names explain the main
workflow. This is a UI revision; estimation and statistical qualification are
unchanged.

| Surface | Change |
|---|---|
| Landing | Two adjacent starting actions, with direct sample access as a quieter third action. Remove the introductory paragraph and repeated guide instructions. |
| Application information | Collect version, limitations, license and citation in one sidebar disclosure; remove the permanent subtitle and separate main-page disclosure. |
| Data entry | Place the table/file input beside parsing options. Keep the paste example available on demand. Combine counts and the preview in one card. |
| Preflight | Keep one successful readiness notice. Move the redundant next-action row and clean row audit into disclosures. Missing/invalid rows still expose their warnings, reasons and audit downloads. |
| Results | Open the combined measures table first, followed by the unchanged convergence diagnostics. Individual tables and reading guidance remain available in closed sections. |
| Diagnostics and figures | Start with the content selector; remove headings and captions that repeat the selected section or explain the selector. |
| Report | Group file format, detail option, sharing notice and download button. Put reference material after the action and show the three follow-up actions together. |

Column mappings, data-source selection, input labels, numerical holds and
SE/CI interpretation notes remain in place. The column auto-detection
percentages are no longer repeated beneath controls; the selected columns
remain visible and editable. Navigation labels may be visually collapsed where
their options name the destinations; meaningful data-entry labels remain visible.
The accessibility regression checks now distinguish those uses rather than
prohibiting all visually collapsed navigation labels.

The retained record is `validation/generated/ui_layout_rewrite_20261023/`.
It contains the starting source, the UI diff, test logs and Chrome captures.
AppTest covers the guide, sample/paste/upload routes, language changes,
history/state behavior, report choices, and inference holds. The input-audit
fixtures specifically check complete, missing and fractional-score data.
The first combined run retained one expected failure from old button wording;
the assertion was updated and rechecked rather than weakening its route checks.

Local Chrome checks use a separate temporary profile and synthetic sample data.
The desktop viewport is 1440×1000; the narrow viewport is 390×844. The latter
keeps the report settings and save action in one column, without document-wide
horizontal overflow. The final input layout is shown in
[the desktop capture](../validation/generated/ui_layout_rewrite_20261023/final_setup_ja_desktop.png);
the report settings in [the narrow capture](../validation/generated/ui_layout_rewrite_20261023/report_ja_mobile.png).
These local artifacts are excluded from a fresh Git checkout.

The final UI suite passed 195 tests; a separate state/qualification suite passed
60, and all 56 built-in self-tests passed. These suites overlap and their counts
are not a count of independent usability observations.
DOM/accessibility-tree checks confirm that navigation options retain their
names; AppTest checks selection retention across routes and languages.
This is a focused layout and interaction review,
not completion of the repository's full browser/assistive-technology matrix
or evidence of novice-user comprehension. The hosted app was not deployed.

## Earlier checkpoints

Earlier checkpoint, 2026-10-16: the [external-group output specification](#external-group-output-contract-20261016)
separates the model's identification from the requested comparisons, and keeps
numerical release distinct from both. Three design-only cases exercise the
proposed states. Bilingual guidance is available through the existing global
Help route; no new contextual adapter or group estimator is activated.

Earlier checkpoint, 2026-10-16: the [external-group design review](#external-group-identification-20261016)
checks three designs under two interaction specifications. Shared references
identify the additive lesson but leave group-DFF ambiguities; extra symmetry
constraints can conceal missing-cell limitations. The school lesson now names
the next design action. No new estimator or inferential clearance is released.

Earlier checkpoint, 2026-10-16: the [group assignment review](#group-assignment-review-20261016)
checks one Person/group mapping before the classical DIF screen, displays
distinct-person assignment counts, and adds the school lesson's group CSV.
Observed coverage is separate from identification or validation of a group model.

Earlier checkpoint, 2026-10-16: the [DFF learning review](#dff-learning-review-20261016)
puts the plot and a compact reading table before technical details, connects
the two-school lesson to its intended pair, and adds five bilingual reading
questions. Automated journeys do not establish novice comprehension.

Earlier checkpoint, 2026-10-16: the [MML anchor and DFF review](#mml-anchor-dff-review-20261016)
rejects unsupported Person anchors that ordinary MML previously ignored,
checks supported non-Person constraints, and clarifies the existing conditional
DFF heatmap. External-group curves and inferential acceptance remain open.

Earlier checkpoint, 2026-10-16: the [two-school linking lesson](#two-school-linking-20261016)
adds a joint-calibration example and bilingual guidance. Anchor counts and
constraint agreement no longer imply validated linking. Browser, reference
stability and human comprehension acceptance remain open.

Earlier checkpoint, 2026-10-16: the [report Q&A review](#report-faq-review-20261016)
adds four English answers and their references to all publication formats, with
Japanese guidance beside downloads. Questions cover overfit, PCM versus RSM,
Wright-map directions and bias interactions. Browser and human acceptance remain open.

Earlier checkpoint, 2026-10-16: the [metric definition review](#metric-definitions-review-20261016)
puts Strata and fit definitions next to results and in report Notes, removes
contradictory EAP aggregation, and expands the planned comprehension tasks to six.
Human comprehension and browser acceptance remain untested.

Earlier checkpoint, 2026-10-16: the [appendix review](#report-appendix-review-20261016)
preserves emphasis and list hierarchy in document exports, corrects stale
population-scale wording, and propagates calibration holds into dependent
writing guidance. Browser and human comprehension acceptance remain open.

Earlier checkpoint, 2026-10-16: the [report follow-up review](#report-followup-review-20261016)
connects missing reporting information to existing download, review and writing
screens. The formative task sheet now matches the active recruitment target;
actual human and browser acceptance remain untested.

Earlier checkpoint, 2026-10-16: the [basic-report and table review](#basic-report-table-review-20261016)
starts with a bounded report, adds bilingual reading guidance, and separates
table numbers, titles and Notes in Word/PDF/HTML. Detailed output remains
available. Browser and novice-user comprehension acceptance remain open.

Earlier checkpoint, 2026-10-16: the [report-output review](#report-output-review-20261016)
connects the first writing-assessment tables and figures to reading guidance,
fixes the category-curve image fallback, and preserves interpretation holds in
publication outputs. Full reporting, browser and human acceptance remain open.

Earlier checkpoint, 2026-10-16: the [Help readability review](#help-readability-and-internal-copy-review--2026-10-16)
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

The October 14 proposal aimed for three novice and two experienced reviewers.
It is superseded by the [October 16 task sheet](#reporting-comprehension-protocol-20261016)
and the active roadmap: five novice graduate students and two experienced
reviewers. This remains a formative sample, not a representative survey. Use the same
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

<a id="reporting-comprehension-protocol-20261016"></a>

### 報告導線の理解を確認する手順 — 2026-10-16（実施前）

本節は上のUX1手順の課題・採点方法を具体化する。参加者は**初学者の大学院生5人＋経験者2人**を
目標とし、10月14日の3人＋2人という募集案を更新する。参加者・実施担当・日程は未確定であり、
**実施0人、成功率・誤解率は未算出**。アプリの自動テスト結果を人の成績に含めない。

**実施前に固定するもの。** 候補版・入力・設定・保存文書の識別情報、ブラウザーと画面幅、言語、
課題の文面、制限時間、採点者を記録する。[作文例](apa_manuscript_template.md#writing-assessment-first-report)の
960評定の通常例と、R1の全評定を最高点にした保留例を使う。通常例も推論全般の合格例とは扱わない。
SEと事後SDの課題には、既存の固定SD MML結果から個人の行とNoteを提示し、推定法・尺度を明示する。
比較カードは[既存の出力点検](#basic-report-table-review-20261016)を再利用できるが、配布前に
対象版の値・Note・保留を照合する。新しい数値研究を始める手順ではない。

参加者が慣れている言語を記録する。5人＋2人の確認で日英両方を受入済みとはしない。
各言語の人数・課題数を分け、ブラウザーの日英4経路、既存のアクセシビリティ受入も別に残す。

#### 操作課題（初学者、20分以内）

参加者には次の課題文だけを先に示す。採点基準と回答例は、初回回答が終わるまで説明しない。
アプリ内のHelpと通常の利用案内は参照してよい。担当者が操作箇所や解釈を教えた場合は援助ありと記録する。
有効なサンプルを開いた時点から計時し、読み込み待ちも所要時間に含める。担当者による中断は別記する。

> 作文の評価者の厳しさと評定の適合について、研究会で説明する準備をしてください。
> 列と分析条件を確認して分析し、分からない用語をHelpで調べて同じ結果に戻ってください。
> 必要な結果表・図を選んで報告書を保存し、カテゴリーの使用状況もCSVで保存してください。
> 原稿を書き始めるためのファイルを開き、研究者自身が補う情報を挙げてください。

- [ ] 人数・評定数、列、モデル・推定法を確認し、明示的な分析操作で結果を得る。
- [ ] Helpから同じ分析に戻る。分析IDと設定を採点者が確認する。
- [ ] 基本レポートを取得し、評価者の推定値・適合度とNote、関連する図を指し示す。
- [ ] カテゴリーの使用度数のCSVを取得し、実測度数と確率曲線の役割を区別する。
- [ ] 原稿テンプレートを開き、自動入力される分析情報と自分で書く研究内容を挙げる。

5項目すべてを20分以内・担当者の援助なしで完了した人を操作成功と数える。
未完了・援助あり・中断を分けて残し、途中の参加者を分母から外して成功率を上げない。
機器障害で評価不能なら理由を記録し、未実施扱いの課題を補ってから目標を判定する。

#### 読解課題（操作後、時間は別記）

Strataとunderfitの判断課題を追加し、当初の4課題から6課題へ更新した。まだ実施前であり、
旧課題と新課題の成績を同じ成功率にまとめない。

| 課題文 | 正答に必要な内容（採点者用） | 重大な誤解の例 |
|---|---|---|
| R1とR4について、厳しさの推定値とInfit・Outfitから説明してください | この例の向きではR4が厳しい。厳しさとモデルへの適合を区別し、適合だけで正確さ・公平さを保証しない | 厳しい評価者だから不適合／不正、または適合しているから公平だと断定する |
| 評価者の分離信頼性が高い結果を、どのように説明しますか | 厳しさの違いを区別できる程度であり、評価者間の得点の一致度ではない。設計・不確実性に依存する | 高い分離信頼性を評価者間の高い一致として報告する |
| JMLEのSEと、このMMLの個人行にあるSDは何を表しますか | Noteから、前者は他の推定値を固定した条件付き近似、後者は個人の潜在特性の条件付き事後の広がりと説明する。後者をEAPの標本SEと呼ばない | 両方を同じ標本SEと説明し、未検証の区間精度を保証する |
| 保留例の空欄をどう説明し、報告前に何をしますか | 空欄はゼロや差の不在ではない。理由を読み、入力・符号化・割付・採点状況を点検し、確定的報告を控える | 技術CSVの値を確定的に報告する、空欄を0にする、警告を消すため評定を変更する |
| H = 3と表示されています。学生・演奏者を3群に分けてよいですか | Gから求める概念的な水準数の目安で、実在する群や境界ではない。分布・不確実性と分類目的を別に確認する | Hを丸めて群数・受賞区分を確定する |
| ピアノ審査員のInfit = 1.1、Outfit = 2.1なら何を確認しますか。MnSq = 1.2との違いは何ですか | underfitの意味と注意表示の閾値を区別する。元の評定、観測数、割付、評価基準を確認し、値だけで除外しない | underfitを能力不足・不正とみなす、1超をすべて除外する、1を完全適合の証明とする |

各課題を「正答／一部正答／誤答／未実施」で記録する。意味が同じ説明を認め、専門用語の暗唱を求めない。
6課題すべてが正答の人を読解成功とする。初学者は操作・読解それぞれ**4/5人以上**を目標とし、
成功人数／実施人数と各課題の内訳を併記する。経験者2人の成績を初学者の分母に混ぜない。
経験者は同じ資料で、論文に必要な設計・診断・留保が欠けていないかも点検する。

7人のレビューで重大な誤解が見つかった場合、該当文言・導線を修正して再確認する。
初回結果を残し、再確認の人数・版・課題・結果を別記する。再確認で0件でも、初回から0件だったとは記さない。
公式順位や不正の断定が会話に現れた場合も重大な誤解として残す。言語評価だけを試して、
ピアノ・スケートでの理解を確認したとは扱わない。

記録欄：匿名参加者ID、経験区分、言語、候補版、入力／分析ID、所要時間、完了5項目、
援助した箇所、6課題の回答原文と採点、重大な誤解、修正内容、再確認結果。
現在は全課題が未実施であり、誤解0件という実績はない。

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

<a id="report-faq-review-20261016"></a>

### 出力末尾のQ&Aと参考文献 — 2026-10-16

**実装した範囲。** 基本版・詳細版のWord・PDF・HTMLに、結果と補足資料の後、Referencesの前に
英語の4項目のQ&Aを追加した。画面では「結果の読み方：4つのQ&A」を開くと、日英で読める。
説明・具体例・次の確認をまとめ、質問と回答が改ページで分断されにくいようにした。
2項目目は利用者の訂正に従い、PCM（partial credit model）とRSM（rating scale model）の比較とした。
共通の段階構造と水準別の構造、データ量・カテゴリー使用の確認を説明し、同じ得点ラベルだけで
モデルを自動選択しない。ほかの項目では、厳しさの類似とoverfit、ファセットの向き、主効果と
交互作用を区別する。バイアスの有無の結論と、各分析の解釈保留は維持した。

**文献の扱い。** Q&Aの引用を既存の参考文献生成処理に渡し、重複を除いて収録する。
FACETS 64-bitの不適合・分離信頼性・向き・バイアス、WinstepsのRSM/PCM構造の一次資料を確認した。
モデルの基礎文献と現在の実装の検証状況は区別する。発行日不明のマニュアルはn.d.と参照日を分け、
同一著者の年付き文献より前に置く。既存のBibTeX/RISテストで、日付なし資料が未分類になる問題を
検出したため、オンラインマニュアルとして著者・題名・URL・参照日を保持するようにした。
参照年を発行年へ変換していない。

**確認したこと。** 初回95件、参考文献関連を含めた修正後150件が通過した。改ページ・並び順の
変更後78件、PCMへの訂正後は日英整合性も含む**87件**が通過した。重複を含む実行記録であり、
件数を合算しない。訂正前のPCA説明を含むキャッシュも無効化し、最後に出力関連27件が通過した。
基本／詳細×3形式で英語Q&Aと引用の一意な収録を確認し、アプリの日英表示、
初期状態の折りたたみ、旧キャッシュの無効化、再表示での結果保持を確認した。
通常／保留の保存済み合成例から3形式の基本レポートを生成し、Word・PDFを画像化して点検した。
図は今回変更していないため、確認済みの4画像から基本版の2画像を再利用した。
コンパイル、doctor、release-checkも通過した。

[実行ログ・文書・画像・出力ハッシュ](../validation/generated/report_faq_20261016/)はローカルに保持する。
初学者の読解は未実施であり、既存の6課題とは別に、Q&Aの4項目を説明できるか確認する。
公開先、実ブラウザー、狭い画面の受入、全テストの一括実行は今回の確認に含めていない。
正式G3以降の研究は再開しておらず、推論の妥当性を追加で認定する変更ではない。

<a id="metric-definitions-review-20261016"></a>

### Strata・適合度の定義と判断 — 2026-10-16

**確認した問題。** ヘルプの一部はStrataを「実際に区別できる群数」、MnSq = 1を「完全な適合」、
低い評価者の分離信頼性を「評価者が相互交換可能」と説明していた。分離度を分散比とする語彙集も
残っていた。定義、参考範囲、個別の判断を混同しないように改めた。

- 日英の結果表・適合度画面・分析比較に短い定義を置き、式、仮定例、確認事項を折りたたむ。
  報告書の保存画面、Help、語彙集も同じ説明につなぐ。定義を読む操作で再推定しない。
- Word・PDF・HTMLの表のNoteに、Strataの式・概念的な意味、underfit・overfitの定義、
  単独で群分け・除外を決めないことを残す。旧内容の保存済み文書をセッションが再利用しないよう、
  文書のキャッシュ識別も更新した。
- 分離度Gは、SEを使う場合には調整済みSD/RMSEであり、分散比はG²である。
  標本SD、負の調整済み分散の0への切り詰め、MML個人EAPの別方式を説明する。
  H = 3でも3群や境界を確定しない。G = 0のときH = 1/3になる式の挙動も説明する。
- UnderfitはMnSq > 1の方向、注意表示は0.50・1.50・2.00の別の閾値である。
  1.2は期待値より大きいが目安の範囲内となる。丸め前の境界処理や欠測の扱いは変更していない。
  ZSTDの基準値や小さい絶対値を適合の証明とする説明も修正した。

**重複計算の修正。** 報告画面の二つ目の集計は、MML個人の事後SDにもSE方式の分散差引きを使い、
既存の`calc_reliability`と矛盾していた。この重複表示を除き、保存用文書と同じ、保留を反映した
信頼性表へ統一した。中心となる推定・信頼性の計算式は変更せず、通常／保留の表示と元データの保持を確認する。
古い自動解釈にあった、低い信頼性から評価者の一致を断定する案内も除いた。

**出典の扱い。** 検索には旧32-bit版も現れたため、現行の64-bit版
[FACETS 4：信頼性・分離度・Strata](https://www.winsteps.com/facetman64/reliability.htm)、
[要素別の適合度](https://www.winsteps.com/facetman64/table7.htm)、
[不適合の点検](https://www.winsteps.com/facetman64/diagnosingmisfit.htm)を確認した（2026-10-16）。
マニュアルの群数の例だけを切り出さず、仮想分布に基づく概念的指標という説明も読む。
参考範囲を普遍的な合否基準にはしない。FACETSの母集団／標本SD、Model／Real SEと
本アプリのSE／事後SDの規約が異なりうるため、式が似ているだけで数値一致を認定しない。
今回、TAM・FACETSとの新たな推定比較や停止中の正式研究は行っていない。

**検証。** [関連81件](../validation/generated/metric_definitions_20261016/tests.log)が通過した。
新規7件は、SE・EAP・G = 0の計算例、日英の定義と折りたたみ、丸め前の閾値、
通常／保留の報告画面が同じ集計値を使うことを確認する。既存の4種類の実推定を使う出力テストも含む。
[Help・境界・操作経路の47件](../validation/generated/metric_definitions_20261016/journeys_and_help.log)も通過した。
最後のキャッシュ更新とZSTD文言の確認でも[60件が通過](../validation/generated/metric_definitions_20261016/tests_final.log)した。
これらは重複するテストを含むため、81＋47＋60を異なるテストの総数とは数えない。

[出力記録](../validation/generated/metric_definitions_20261016/inventory.json)には、保存済みの960評定の
通常例／境界保留例から作った基本レポートを記録する。変更していない図は前回の画像を再利用した。
Word・PDFの全ページの配置と、変更した表・Noteを原寸の画像で確認した。HTMLは構造・内容の自動確認である。
レポート本文は既存どおり英語、画面の説明は日英対応である。ローカルの生成物はGit管理外である。
[改ページ・注釈の確認](../validation/generated/metric_definitions_20261016/pagination_checks.json)では、
通常例はWord 5／PDF 4ページ、保留例はWord 6／PDF 5ページで、白紙は0ページだった。
4文書とも定義とNote全文が残っている。[compile・doctor・release-check](../validation/generated/metric_definitions_20261016/checks_final.log)も通過した。
今回、全件の`make verify`は実行していない。

**人による確認は未実施。** [説明用の例](apa_manuscript_template.md#strata-and-fit-decisions)と
[読解6課題](#reporting-comprehension-protocol-20261016)を用意した。合成データでの計算確認と、
初学者が自分の言葉で説明できることを区別する。学生の群分け、ピアノの除外判断を課題に加えたが、
実施0人であり、誤解率は未算出。スケートの具体的な規則や実データの妥当性を確認したとは扱わない。

<a id="report-appendix-review-20261016"></a>

### 詳細レポートの補足資料と保留の伝達 — 2026-10-16

今回の入口は、詳細版の方法・執筆の補足に残るMarkdown記号だった。保存済みの通常／保留の作文例では、
太字の記号とバッククォートがWord・PDF・HTMLに残り、入れ子の箇条書きが同じ階層で出ていた。
PDFでは見出しがページ末尾に孤立する例も実物で確認した。標準の生成例にはMarkdown表はなく、
表の処理は短い合成例と55行の例で別に確認した。

- 共通の限定したMarkdown解釈を使い、強調・斜体・コード表記、箇条書きの階層、表を各形式へ反映する。
  生のタグや不等号をエスケープし、数値・文字列・空欄・記号だけのデータ行を落とさない。
  補足の表は全行を保持し、CSVに存在するとは限らない表へ誤った省略案内を付けない。
- Wordは編集可能な見出し・リスト・表、PDFは折り返しを揃えた箇条書きと実際の表、
  HTMLは入れ子リストと見出し行を持つ表を出す。見出し・親項目を次の説明と一緒に置く。
- 「自由SDを推定していないと書く」と読める古い一律の指示を、実際の固定／推定条件と保留に
  対応する説明へ置き換えた。条件付きSE・事後SD・保留の説明は結果表と同じ根拠から作る。
  生の不確実性状態コードや、文章案の開発用表現を読みやすい文へ置き換えた。

**途中で優先した不整合。** 保留例の結果表は精度指標を空欄にする一方、補足の信頼性・マップ等には
「Claim status: Ready」と報告を促す文が残っていた。全体の保留だけで局所の案内を制約していなかった。
共有のclaim guideで、境界JMLE／自由SD MMLの校正保留を、カテゴリー、信頼性、適合・次元性、マップ、
アンカー・リンク、予測・シミュレーションの6領域へ伝えるように修正した。元の数値は保持し、
記載する場合は確認用であることを近くに明示する。保留中の局所的な解釈文例は、保留を説明する文例に絞る。
通常の条件付き近似や固定SDの留保を、校正全体の保留へ一律に変更してはいない。

入力は前回の960評定の保存済み通常／保留例を再利用した。図は今回変更しておらず、最初の描画後は
同じ画像を再利用した。生成条件とソースハッシュは[出力記録](../validation/generated/report_appendix_20261016/inventory.json)、
変更前のMarkdownは[初回点検](../validation/generated/report_appendix_20261016/baseline.json)と同じ場所に残す。
[148件の関連確認](../validation/generated/report_appendix_20261016/tests_qualified.log)は、4種類の実際の推定結果、
文書の値・保留、記号・階層・表、出力キャッシュ、報告用の集計再利用、境界・準備状態を対象に通過した。
途中の61件・118件は、最終の保留伝達を追加する前の確認として区別する。
一度、存在しないテストファイル名を指定して収集前に終了したログも残し、修正後に実行した。
[操作経路・案内の11件](../validation/generated/report_appendix_20261016/journeys.log)も通過し、
4種類の実際の分析で数値・分析ID・保留の保持と、表示操作による再推定0を確認した。
最後に状態別件数の`np.int64(...)`という表示と参考文献案内の重複した語句を修正し、
[最終文言の回帰確認38件](../validation/generated/report_appendix_20261016/tests_copy_final.log)も通過した。
[compile・doctor・release-check](../validation/generated/report_appendix_20261016/checks_final.log)を確認し、
今回も全件の`make verify`や停止中の正式研究は実行していない。

Word・PDFの補足ページを画像化して配置を確認し、強調・入れ子・保留・改ページの箇所を原寸で確認した。
[最終出力の改ページ確認](../validation/generated/report_appendix_20261016/pagination_checks.json)では、
通常例はWord・PDFとも17ページ、保留例はWordが17ページ、PDFが18ページだった。
4文書とも白紙、補足見出しだけのページ末尾、未変換の強調記号・バッククォートは0件だった。
HTMLは構造と内容の自動確認であり、実ブラウザーでの見た目・キーボード受入ではない。
文書全体のAPA準拠、全モデル・全設定の妥当性、人による理解を今回の結果から認定しない。
詳細版には長い説明や確認用の情報が残るため、初学者の入口は引き続き基本版とし、
必要な引用・説明に絞る作業は実際の読解確認と併せて進める。正式G3以降の停止は維持する。

<a id="report-followup-review-20261016"></a>

### 基本レポートの後に必要な操作 — 2026-10-16

共有の案内文が「上の『論文を書き始める』」を指していたが、詳細画面では同じ位置にその操作がない。
さらに、実際には著者が研究内容を記入するテンプレートを、完成した方法・結果の下書きのように
案内していた。共有文を場所に依存しない説明へ直し、初心者向けの「報告書を作る」には
折りたたんだ「論文にする前に確認すること」を追加した。

- 設計・除外・カテゴリー・残差を確認する表名は、現在の保存一覧と同じ翻訳を使う。
  「表・図を保存」「確認事項」「論文を書き始める」のボタンは既存の選択先へ移動する。
- 空欄をゼロや差の不在と解釈しないこと、保存・詳細化で保留が解除されないこと、
  原稿に著者が補う情報を日英で説明する。詳細な資料を最初から展開する構成には戻さない。
- [作文の利用案内](apa_manuscript_template.md#writing-assessment-first-report)に4つの自己確認課題と
  折りたたんだ回答の着眼点を追加した。開発中の受入状況の説明はこの監査記録で管理する。
- [形成的レビューの課題・採点表](#reporting-comprehension-protocol-20261016)を既存UX1手順に接続し、
  旧3人＋2人の案が現行5人＋2人の目標を上書きしないようにした。実施0人は変わらない。

ソースSHA-256: `92aec55dc002de3633e5e6ddb970f957b5daf642bb1cfda903482cf2bb56b19d`。
[関連テスト125件](../validation/generated/report_guidance_20261016/tests.log)が通過。
案内ボタンは日英×保留有無×確認結果3状態の12条件で既存の画面へ移動できることを確認した。
実際の入力を推定する4ケース（通常／境界JMLE、固定／自由SD MML）では、英語の原稿・確認事項への
往復、日本語の表保存への移動を追加確認し、値・分析ID・保留、キャッシュしたPDFを保持した。
表示操作の再推定は0。これは新しい不確実性研究や、言語教材全体の受入ではない。

文書生成の数値・書式処理は今回変更していない。[基本レポートの出力点検](#basic-report-table-review-20261016)を
維持し、今回は画面の案内と操作の接続を確認した。
[compile・doctor・release-check](../validation/generated/report_guidance_20261016/checks.log)も通過し、
[変更文書のローカル参照先](../validation/generated/report_guidance_20261016/checks.json)に欠落はなかった。
実ブラウザーの見た目・狭い画面・キーボード・人による理解は未確認。
正式G3以降の停止は維持する。全件の`make verify`は今回も再実行していない。

<a id="basic-report-table-review-20261016"></a>

### Basic report, table titles and Notes — 2026-10-16

前回の点検では出力の網羅性を優先したため、結果より先に長い説明が続き、表題にも注意書きを
詰め込んでいた。Wordの青い格子、PDFの全セル罫線、HTMLの別書式も混在していた。
[APA第7版の公式説明](https://apastyle.apa.org/style-grammar-guidelines/tables-figures/tables)を確認し、
論文全体の準拠を自動的に主張する文を削除した。今回の対象は、報告の入口と結果表の構造である。

- 画面は基本レポートを既定とし、方法を説明する→結果を読む→論文を書く、という日英の3段階を案内する。
  言語・ピアノ・スケートの問い別の着眼点は展開して読める。教材そのものの完成を意味しない。
- 基本版はファセットの推定値・適合度、信頼性の2表と、マップ・カテゴリー曲線を、利用可能な範囲で出す。
  個人別の行、相互作用表、追加の2図、方法・執筆の補足は詳細版に置く。通常の数値・元のCSVは保持する。
  基本版では不要な2図と長い補足本文を生成せず、基本／詳細の変更は保存ファイルのキャッシュを更新する。
- 各結果表は太字の番号、短い斜体の表題、水平線だけの本体、下のNoteで構成する。
  表と本文の文字サイズを揃え、区間は角括弧で1列にまとめる。全行に共通する水準・資格をNoteへ移し、
  混在する不確実性の種類は行ごとに残す。基本版のNoteは実際に表示するファセットの種類から生成する。
- 解釈保留の注意は概要・表・図に残す。保留のSE・区間と信頼性の精度指標は空欄を維持し、
  古い保存値が肯定的な状態を持っていても現在の保留を上書きしない。分離信頼性と一致度も区別する。
- 全形式で表の表示上限を50行に揃え、超過時はNoteに表示範囲と全行CSVへの案内を記す。
  列の切り捨ては行わない。Word/PDFは複数ページに渡る表の見出しを繰り返す。

入力は前節の作文の通常／保留の保存済み推定結果を再利用した。追加で固定SD MMLの既存小規模
テスト入力（3人・12評定）を用い、事後SDと構造パラメータの不確実性が混在する8列の表示を確認した。
最初のPDFでは`0.792`が`0.79`と`2`に折り返され、数値が読みにくかった。
この失敗を[回帰テスト](../tests/test_publication_table_layout.py)で再現し、列幅・余白を修正した。
Noteも意味を保って短縮し、Wordで注記だけが分断されないようにした。

描画・出力確認時のソースSHA-256: `31d8fa558a2d25b588cb3007013b1599ee0b3e86796dee2b6ce11b2fde3c79cb`。
[保存物と出力範囲](../validation/generated/apa_reporting_20261016/inventory.json)、
[210件の回帰ログ](../validation/generated/apa_reporting_20261016/regression_final_pagination.log)、
[数値の折り返しを再現した修正前の失敗](../validation/generated/apa_reporting_20261016/numeric_wrap_before.log)を保存する。
210件の関連テストは通過した。その後、PDFの最後の図の余白が白紙を作る問題を
[再現](../validation/generated/apa_reporting_20261016/blank_page_before.log)して修正し、
[出力・表・保留の40件](../validation/generated/apa_reporting_20261016/export_final.log)と
[文書生成の3自己テスト](../validation/generated/apa_reporting_20261016/document_selftests_final.log)が通過した。
先行する209件と56自己テストは、最終の列幅・Noteの短縮より前の確認として区別する。
補助的な混在例の生成スクリプトで図のスタブを戻さず自己テストを呼んだ失敗は、アプリの生成失敗とせず、
新しいプロセスでスタブなしの3自己テストを実行し直した。今回は全件の`make verify`は再実行していない。
最後にベータ機能一覧の古い「自動Abstract」「Wordは60行×12列」「PDFだけコード付録を除く」という説明も
現在の出力に合わせて修正した。最終ソースSHA-256は
`6efb56f03b04f91591ceba057ab3e4d00a1d671d09233cbf2a9bf1d6a421a6d8`。
[compile・doctor・release-check](../validation/generated/apa_reporting_20261016/checks_final.log)と
[変更文書のローカル参照168件](../validation/generated/apa_reporting_20261016/final_checks.json)を確認した。

基本版のWordは通常／保留とも4ページ、PDFは通常4・保留5ページとなった。
[改ページ確認](../validation/generated/apa_reporting_20261016/pagination_checks.json)では白紙0ページである。
基本版はWord/PDFの全ページ、詳細版は変更した結果表のページ、固定SD MMLは混在する不確実性の表を目視点検した。
表・Noteの分断、数字の途中改行、罫線、余白を確認した。HTMLは生成内容・数値・Notesの自動確認であり、
実ブラウザーでの見た目の受入を代替しない。

残る受入は実ブラウザーの狭い画面・キーボード操作、初見の大学院生による表の選択と説明、
詳細な補足資料に残るMarkdownの全体整形である。生成できたこと、表の書式、推論の適切性、
人の理解は別々に評価する。正式G3以降の停止は維持し、新しい推定法や区間の妥当性を認定していない。

<a id="report-output-review-20261016"></a>

### Reporting outputs and first writing-assessment example — 2026-10-16

出発点は`a5ce52da7`。この節は、その後の作業中の出力修正を記録する。
[報告・学習要件](mfrm_reporting_requirements.md)の初回点検として、作文の合成データで
通常例と保留例を作成し、CSV・HTML・PDF・Wordを照合した。
研究の停止範囲を変更したり、新しい区間の被覆率研究を実行したりする作業ではない。

**入力と条件。** 組込み作文例、シード20240101、30人×4評価者×2課題×4観点、960評定、
0–5の6カテゴリー。RSM・JMLE、最大100反復、アンカーなし、単位重み。
保留例はR1の全評定を最高点にした。どちらも既存の推定・診断処理で計算し、
保存した同じ結果から修正後の出力を生成した。通常例は40要素が条件付き近似、
保留例は40要素が校正の解釈保留。通常例を、全推論の資格を満たす合格例とは扱わない。

| 発見した問題 | 修正と確認 |
|---|---|
| CSVの選択肢が多数の技術的な名前だけだった | 主要14表を日英の内容名で選び、読む列と限界を添える。詳細表は切替で選べる。共有設定で選択肢が消えたときは有効な表へ戻し、元の保存名・内容・一括出力は維持 |
| 説明から参照するカテゴリー集計・段階順序が通常の出力一式に含まれなかった | `category_diagnostics.csv`と`step_order.csv`を追加。既存診断との一致、度数の合計・割合、入力の不変性を照合。内部用の全表数は通常103→105、保留99→101。公開設定では表数が異なる |
| 画像の代替処理で確率曲線が誤っていた | 長形式の`Category`・`Probability`・`CategoryValue`を独立した曲線にしていた処理を修正。カテゴリーごとにThetaとProbabilityを描く。RSM・PCM・GPCMで画面用データとの一致、ラベル、範囲、各点の確率和1を検証 |
| 公表用の推定値表から適合度が落ちていた | Infit／Outfitを追加し、期待値1と機械的な合否にしない説明を添える。SE・区間の意味と保留を維持。Word/PDFの列幅を調整 |
| 保留例の図だけを見ると制約が分からなかった | JMLEの校正保留・自由SD MMLの保留を、共通の4図の説明へ反映。図と説明を同じページに保つ |
| 相互作用の報告表に開発用状態コードが出ていた | 表の列名と計算状態を人が読める語に変更し、`Use = Withheld`を維持。古い結果・未知の状態を推論可能へ昇格しない。元の技術表は変更しない |

**読み方。** [作文評価の3表・2図](apa_manuscript_template.md#writing-assessment-first-report)に、
問い、1行の意味、条件、操作、読む箇所、短い報告文、保留例の対応を日本語で追加した。
READMEの文書案内から開ける。アプリ内には14表の個別の読み方を置いたが、
教材全文の画面内接続・英語版・5出力だけをまとめる専用ボタンは未実装。

**実物の確認。** 通常／保留の2条件で各4図と3形式の報告書を生成した。
PDFとWordを画像にし、変更した推定値・適合度・相互作用の表、確率曲線、保留の説明と改ページを確認した。
HTMLは生成と内容の照合までで、ブラウザーでの描画確認は未実施。
Wordの変換はサンドボックス内でLibreOfficeが異常終了したため、許可されたローカル変換で確認した。
アプリの依存関係には変換用の依存を追加していない。

**操作と回帰確認。** 実際の作文の保存結果を使い、通常／保留×日英の4条件で主要14表、
計56回の選択をAppTestで確認した。表示名、CSV保存ボタン、ファイル名に矛盾はなく、
例外・表示操作による再推定は0だった。入力からの教材全体の受入ではなく、表選択画面の確認である。
対象回帰テストは82件、単独の既存画面テストは88件が通過した。
途中でソースを変更した全体検証は、読込済みソースとの不一致を検知したため中断し、
その失敗ログを保存した。アプリのコードを固定してやり直した`make verify`では、
自己テスト56件が成功し、通常の回帰テストは**3,126成功・1失敗・2スキップ・84対象外**だった。
失敗した1件は相互作用の公表用表を旧列名で読むテストで、計算・技術表の保留に失敗したものではない。
期待値を新しい表示名へ合わせ、計算状態が未検証であることと解釈保留を確認する条件を維持した。
その後、`test_bias_cell_math.py`と`test_bias_publication_integration.py`の**51件がすべて成功**した。
全体検証後の変更はこのテストの期待値と文書だけで、アプリのコードは変更していない。
全件をもう一度実行した成功ログとは区別し、`make_verify_final.log`と`bias_contract_followup.log`に残す。
追加した文書内参照26件と`git diff --check`も確認し、ロードマップの履歴部分は変更していない。

保存先はローカル確認用の`validation/generated/reporting_review_20261016/`。
入力、保存結果、CSV、3形式の報告書、図、画像、実行ログ、再生成スクリプトを残す。
`final_inventory.json`に出力数・状態・ソースSHA-256を記録した。
このディレクトリはGit管理外で、リポジトリの読者が取得できる根拠とは区別する。

**残る点検。** 全モデル・全条件のR1–R8受入、教材全体の日英4経路、実ブラウザー、狭い画面、
キーボード、初学者による理解、変更後の待ち時間は未確認。
自動報告本文には長い確認リストやMarkdown記号の残存があり、論文用の簡潔な本文とは分けて整理する必要がある。
冒頭も「Validation stance」等の実装側の説明が多く、通常例と保留例で同じ定型要旨になる。
次の本文改善では、何を分析したか・結果をどこまで使えるかを冒頭に置き、技術記録を付録へ分ける。
次は本文・必要な図表・技術資料の役割を整理し、教材から必要な出力を取得できるか確認する。
今回の局所修正を、報告8領域の受入完了やCIの妥当性の証明には数えない。

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

<a id="two-school-linking-20261016"></a>

## Two schools and common recordings — 2026-10-16

**Problem.** Users need to distinguish joint calibration through shared ratings
from fixing previously estimated measures. The existing anchor checklist called
an unanchored run "Not linked" and promoted an anchor-count threshold to "Linked"
or "Ready". Neither description establishes whether two school panels can be
compared. A dummy facet also appeared as a baseline link solely because its
levels were constrained to zero.

**Change.** Add the `two_school_linking` synthetic scenario: 20 students at each
school and ten distinct common reference performances, one shared task, four
raters, three criteria and 360 scores. Local students have local raters; the
references have all four. A collapsed lesson remains available before and after
fitting and explains the question, row meaning, assignment, joint RSM/JMLE fit,
output directions and reporting limits in Japanese or English. It links to
current FACETS 64-bit equating, anchor and nesting guidance; generation values
are explicitly illustrative. The reference IDs are excluded from student school
summaries, and the app makes no validated school-contrast or causal-effect claim.

Switching scenarios now scopes facet selection to the selected sample. This
fixes a real journey failure: switching from the writing example retained Task
although the new example has just one task and its default mapping excludes it.
Language changes retain the selection and fitted result.

Anchor tables now describe count thresholds rather than validated links. The
unanchored checklist marks fixed-anchor checks not applicable and explains joint
calibration. Dummy constraints do not count as evidence of reference stability.
Anchor-chain reporting remains under Review even when all count thresholds are
met, including when older saved Strength values use "Linked". Existing technical
CSV column names remain compatible, with an explicit boundary on their meaning.
Hard-anchor agreement is described as constraint compliance, not empirical drift.
Publication-document cache versioning includes the revised reporting boundary.

**Evidence.** Retained locally under
[`two_school_linking_20261016`](../validation/generated/two_school_linking_20261016/):

- `final_tests.log`: 82 passing checks for the lesson, all sample contracts,
  anchor thresholds, locale parity and the updated registry assertion.
- `regression_tests.log`: 191 passing checks for identification, hard anchors,
  existing journeys, exports, report Q&A and uncertainty holds. The sole failure
  was the expected registry-name set, corrected and included in the final 82.
  Earlier targeted failures (wrong test button, locale key location, retained
  facet selection) and their corrections remain in the local logs.
- The lesson AppTest uses the real default fit and CSV download. Both languages
  retain identical measures and analysis identity, with zero estimation/refit
  calls after fitting; switching to pasted data removes the scenario lesson.
- `example_review.json`: seven pre-existing samples match the before snapshot
  exactly at two seeds (14/14 comparisons); the changed inline self-tests pass.
  The new 360-row example converges. The exact Person/facet design has nullity 0
  with shared references, versus nullity 1 when references are removed or given
  school-specific IDs. Person–Rater component counts are respectively 1, 2 and 2.
  Input CSV, estimates, checklist and rank summaries are retained alongside it.

**Limits and next decision.** This is a teaching-design and implementation
check, not empirical reference-stability validation, a new FACETS/TAM comparison,
or an interval-coverage study. Browser acceptance and human participants remain
zero; formal G3 and subsequent MML studies remain paused. Next, ask readers to
explain removal of common recordings, equality constraints on school means and
zero hard-anchor deviation. Then specify a fixed-anchor reuse example with
independent stability review before adding school-contrast uncertainty.

<a id="mml-anchor-dff-review-20261016"></a>
## MML anchor scope and DFF visualization — 2026-10-16

**Finding and correction.** Ordinary MML accepted a fixed Person anchor in its
audit without constraining the fitted model: changing P1 from -2 to +2 returned
identical optimization vectors, log likelihood and EAP scores. Fixed Person
values, Person group means and a dummy Person facet now stop before optimization,
with a localized explanation. Shared Person IDs can still combine common
recordings in one MML calibration. Fixed and group-mean location constraints for
non-Person facets remain supported by ordinary EM/direct/hybrid/auto engines;
the explicit joint free-SD engine continues to exclude all anchors. These are
implementation boundaries, not evidence of valid linking or TAM equivalence.

**DFF scope.** The existing Bias–Interaction screen provides cell tables,
interactive heatmaps and CSVs. Its estimates hold calibration fixed, including
MML person EAP scores; they are not a full marginal interaction fit. Updated
bilingual guidance defines DFF, explains the positive/negative score direction
beside the heatmap, and removes a contradictory significance claim from the
facet dashboard. The separate classical DIF screen uses observed-score matching
and does not supply external-group response curves. Common reference recordings
link the two-school example, but do not create cross-school ratings of local
students; school-specific rater behavior needs an appropriate crossed design.

**Evidence.** Local logs and artifacts are retained under
[`mml_anchor_dff_review_20261016`](../validation/generated/mml_anchor_dff_review_20261016/).
The before-fix probe is preserved. Regression checks cover rejected Person
constraints, returned non-Person constraints across four ordinary engines and
RSM/PCM/GPCM, Japanese/English error recovery, and existing identification,
joint-engine guards and bias-output boundaries. Two real fits of the synthetic
360-rating school example (JMLE and fixed-SD direct MML) converged; each actual
heatmap contains 12 cells, agrees with its numeric table, and has no significance
markers. Interactive HTML examples and CSVs are retained alongside the check.
Source fingerprints and test outcomes are recorded in `final_review.json`.

**Limits.** Constraint agreement does not establish convergence of every test
fit, reference stability or equating validity. Figure-data checks are not browser
or human acceptance. DFF error rates, power and interval coverage remain
unqualified; every screening claim remains withheld. No formal G3 or subsequent
MML study was started. The next scope is a learning example that distinguishes
general severity, conditional interaction and external-group DIF before any
new inferential model or school-comparison claim.

<a id="dff-learning-review-20261016"></a>
## DFF learning path — 2026-10-16

**Problem and change.** The Bias–Interaction screen exposed multiple inference
and audit tables before its plot, while the school lesson's default pair was
Person × Rater. The first view now shows the interpretation hold, an optional
reading guide, the heatmap and four columns (two facet levels, rating count and
bias). Technical calculations and CSVs remain available in a collapsed section.
The school example starts with Rater × Criterion and has a button to open the
existing diagnostic view. Explicit pair choices survive rerenders and language
changes; navigation does not launch estimation. The lesson says that it adds
no intentional DFF, so sampling variation must not be treated as known bias.

**Missingness and language.** Heatmap pivots retain rows and columns whose
estimates are all missing when other cells are estimable. Missing cells have no
hover value; a true zero stays zero and has a nondegenerate centered color range.
Counts, metric titles and hover labels follow the interface language. The
reading guide distinguishes the standard rater-severity sign from the score
direction of conditional bias, and warns against comparing color intensity
across plots with different ranges. No estimator or inference release changed.

**Reading questions for human review (not yet administered).** Record correct
answers out of five, correctness per question and any consequential
misinterpretation. These are a separate DFF module, not additional claims about
the time target of the earlier general reporting task. Participants: **0**.

| Question | Required explanation |
|---|---|
| Does an overall severe rater have DFF? | Overall severity is distinct from a combination-specific departure beyond the additive prediction |
| Does negative plotted bias mean greater leniency? | It means a shift toward lower scores than predicted; a logit is not a raw score point |
| Does a blank or unflagged cell establish no bias? | Neither establishes absence; distinguish missing evidence from an estimated zero |
| Can a dark cell be reported as significant DFF? | The plot holds calibration fixed (including MML EAPs); current tests and intervals remain unqualified |
| Do common recordings alone reveal school-specific differential treatment? | The lesson's local students have only local raters; cross-school student ratings need to be designed |

**Evidence and limits.** Local AppTest and numerical display checks are retained
in [`dff_learning_review_20261016`](../validation/generated/dff_learning_review_20261016/).
The initial failing run exposed the wrong default pair and a test-harness state
assignment error; both are retained with the subsequent corrections. Real
JMLE/MML fixtures exercise both languages, unchanged bias values and CSVs,
plot/table correspondence and continued inference holds. The full school
journey checks navigation and preservation of user settings without refitting.
`final_review.json` records source fingerprints and final test outcomes.
The final selected regression run passed **225/225** tests, including the
JMLE/MML × Japanese/English display checks and the full school lesson journey.
Browser control was not available in this session; actual browser layout and
human comprehension acceptance remain open. Formal G3 and later MML studies
remain paused. The [FACETS 64-bit pairwise-bias manual](https://www.winsteps.com/facetman64/table14.htm)
supports the distinction between overall measures and local departures, but
does not validate this application's inferential procedures.

<a id="group-assignment-review-20261016"></a>
## External groups and observed assignments — 2026-10-16

**Findings.** The classical DIF interface converted Person/group rows into a
dictionary, silently retaining the last duplicate Person. Blank strings could
remain group labels, and automatic CSV type inference could strip leading zeros
from IDs. The core and interface now share validation: Person IDs must be
nonblank and unique after string conversion; missing/whitespace-only groups
remain missing. Text IDs are preserved on CSV/Excel import and matched exactly
to response IDs. Unmatched response persons and unused group-file IDs are shown.

The old explanation that averaging removes rater severity has been corrected
in the screen's interpretation text, audit output and Help. The classical DIF
screen uses mean ratings and observed-score matching, not an MFRM severity
adjustment. A missing flag is displayed as information, not a green clearance.

**Reading path.** The two-school lesson supplies a Person/School CSV with 20 A,
20 B and ten ungrouped reference recordings. Its Rater × School count table
has 20/0 for the A raters and 0/20 for the B raters: **0/4 raters** have observed
students from both schools. The common recordings still connect calibration;
they are not assigned to either student group. In classical DIF, users can
inspect the same kind of coverage before explicitly running any screen.
Every person/facet-level pair counts once; outcomes and weight magnitudes do
not enter the count. This is not a sufficient-design or identification check.

**Computation and exports.** Coverage is aggregated over observed assignments
first. Dense matrices above 50,000 cells are not allocated. The CSV contains
only observed combinations as FacetLevel, Group and Persons, with the selected
Facet and GroupVariable recorded, whether or not
the matrix is displayed; zero cells are omitted. No new dependency, calibration
engine, group-specific response curve or inferential qualification was added.

**Evidence.** Retained under
[`group_assignment_review_20261016`](../validation/generated/group_assignment_review_20261016/)
are the before-source archive, mapping probe, test logs and final source hashes.
The first run found a duplicate-DataFrame-index failure in the new count matrix;
aggregating observed pairs directly fixes it and also bounds matrix allocation.
Tests cover synthetic school assignments, added cross-school ratings, duplicate
and missing IDs, blank groups, leading-zero IDs in CSV/Excel, outcome-blind counts,
large matrices, Japanese/English preflight and explicit-run behavior. Existing
classical-DIF point-result regressions are checked separately from formal studies.
The selected regression run passed **160 tests**. After clarifying displayed
column names, retaining column context in CSVs and handling stale facet names,
the final UI/input subset passed **31 tests** (29 repeated checks plus two new
cases: **162 unique passing checks** overall). Source hashes and the scope of
each run are recorded in `final_review.json`.
The 12-replicate operating-characteristic test is deliberately excluded; no
formal G3 or subsequent MML study was started.

**Remaining work.** The [external-group DFF specification draft](../anchor_templates_and_guideline/anchor_user_guidelines.md#external-group-dff-plan)
separates uniform interactions, group ability distributions, identification,
contrast covariance, multiplicity and the conditions for future curves. Its
three model-identification cases were **0/3 complete at this checkpoint**;
the follow-up below completes the bounded structural checks. Count-table tests
alone do not complete them. Human comprehension and actual browser acceptance remain
unperformed. The [FACETS nested-design guidance](https://www.winsteps.com/facetman64/nested.htm)
supports explaining how constraints affect interpretation; the app's proposed
model and acceptance criteria still require independent mathematical review.

<a id="external-group-identification-20261016"></a>
## External-group DFF identification — 2026-10-16

**Finding and specification correction.** With common recordings outside the
student groups and their interaction term absent, the reference recordings
define the rater baseline. Separate group-specific rater departures can each
sum to zero over raters while person locations remain free. Additionally making
each rater's A/B departures sum to zero imposes symmetry around that reference
baseline. It is a substantive model restriction, not just a change of origin.
The earlier draft's double-centering suggestion is revised accordingly.

**Three designs checked.** The small acceptance cases retain two students per
school, two references, four raters, three criteria and five score categories.
RSM step parameters are included. Under the reference-relative specification,
crossed, nested and one-missing-cell designs have rank **20/20, 16/20 and 19/20**.
Adding A/B symmetry changes those results to **17/17, 16/17 and 17/17**. Independent
rational elimination of full parameters plus explicit constraints agrees with
SVD using the app's parameter expansion in all six comparisons.

Ten declared contrasts per specification are checked for estimability, including
two within-school comparisons that must not be labeled between-school DFF.
For the nested reference-relative design, only those last two are identified.
With one absent cell, some differences between observed raters remain identified,
but individual A/B departures do not. Symmetry makes the missing-cell contrasts
identifiable through an added assumption. The
[model guide](../anchor_templates_and_guideline/anchor_user_guidelines.md#external-group-dff-identification)
documents the constraints, contrasts and a constructive probability-preserving
parameter change. Four independent null directions leave every reference
location and rater severity unchanged in the nested reference-relative model.

The actual 50-person teaching design also has four such directions; its additive
model has full rank. A single cross-rated student removes only one direction
under the reference-relative specification, although it removes the sole
remaining direction under symmetry. Neither result is a precision calculation
or a recommendation that one extra rating is sufficient.

**User guidance.** A Japanese/English note beside the lesson's assignment table
explains that linking the scale does not separate each rater's A/B difference.
It directs users toward overlapping ratings of students from both schools and
checking the intended contrasts and precision. Matrix ranks and development
identifiers are confined to technical documentation and evidence.

**Evidence and limits.** The reproducible bounded script and tests are
`validation/external_group_dff_identification_check.py` and
`tests/test_external_group_dff_identification.py`; retained design matrices,
contrast results, source hashes and logs are under
[`external_group_dff_identification_20261016`](../validation/generated/external_group_dff_identification_20261016/).
The final selected regression run passed **60 tests**, including the 14 new
structural checks, existing group-input checks, bilingual lesson navigation,
locale parity and JMLE/MML DFF reading views. The two rank methods agree for
**6/6 design/specification combinations**; both estimability methods agree for
**60/60 declared contrast checks**. These are separate denominators from the
pytest count. An AST comparison with the retained before-source confirms that
the only changed application function is the lesson/setup renderer.
The mathematical cases are **3/3 complete**, not three successful estimations or
an external expert review. RSM structural rank is distinct from finite MLEs,
precision, fit, error rates and marginal-model identification. PCM/GPCM, MML
population assumptions, a new estimator and group curves remain outside this
change. Formal G3 and subsequent MML studies remain paused; actual browser and
human comprehension acceptance remain open.

<a id="external-group-output-contract-20261016"></a>
## External-group DFF output specification — 2026-10-16

**Scope and decisions.** The [output specification](../anchor_templates_and_guideline/anchor_user_guidelines.md#external-group-dff-output-contract)
defines three design states: model identified, partially identified with an
estimable requested contrast, and the requested targets unidentified. Missing,
stale or incompatible evidence is unassessed rather than a fourth mathematical
conclusion. The fixture helper rejects invalid evidence; it is not a production
auditor. No uploaded-data classification or new interaction estimator is added.

The requested family contains eight between-school contrasts, not the two
within-school comparisons retained in the earlier mathematical audit. The
crossed, missing-cell and nested reference-relative designs therefore give
**8/8, 2/8 and 0/8** structurally estimable requested contrasts. Selecting only
the two supported comparisons from the missing-cell case gives 2/2 without
changing the partially identified model status. Selecting the four individual
rater A/B effects instead gives 0/4. The cases and these target-dependent
boundaries are executable, using the independent mathematical checks already
implemented. All three states continue to withhold point estimates, intervals,
p-values and group curves because no corresponding estimator is released.

**Display and export contract.** Name every requested comparison and its
observed support; retain unidentified targets with a reason rather than a
zero value. Separate a short table title from its model/assumption Note. Keep
rank details and machine codes out of the initial view. Contrast uncertainty
requires joint covariance; curve predictions require their own identification
check. The design-only CSV has no effect estimates or SE/CI/p-value columns.
These are future output rules exercised by local fixtures, not newly released
report tables. Existing bias and classical DIF exports retain their own scope.

**Beginner path.** The existing differential-interaction Help topic now includes
“Planning school-group comparisons” in Japanese and English. It distinguishes
each rater's school difference from the difference between raters' school
differences, explains the three possible design states and names the next
checks. The content version is increased to 1.1.0; review remains `in_review`,
and the retained earlier review date is not presented as human acceptance of
the new content. The heatmap popover directs readers to the global Help button,
the differential-interaction topic and its collapsed details.

The first navigation test found that the registered bias-heatmap contextual
edge is not activated: the current acceptance gate permits only the fit-scatter
AppTest pilot without browser evidence. The draft instruction that implied an
active direct link was corrected. No acceptance gate was weakened and no adapter
was added. The test now exercises the existing global Help route in both
languages, verifies the collapsed content and preserves the selected bias view
without estimation/refitting. A fixed locale-key count also increased by two
for the new section title and body.

The full two-school lesson journey subsequently replaced the simplified Help
harness as the acceptance check. It changes the heatmap metric, opens global
Help, reads the new section and returns in Japanese and English. This found
that the metric needed inclusion in the existing navigation-state retention
list. It also exposed live translation lookup during AppTest's widget-value
serialization: diagnostic labels are now resolved at render time, as the
outer section selector already does. The final change keeps the existing
segmented-control defaults and adds no new state-storage abstraction. The
journey checks the retained output object, analysis ID and measure table, with
estimation and refit calls blocked after the initial fit. This does not
substitute for an actual browser test.

**Evidence.** Before-source snapshots, original failed test logs, corrected
test logs, the three JSON output cases, a 24-row comparison CSV and source hashes
are retained in
[`external_group_dff_output_contract_20261016`](../validation/generated/external_group_dff_output_contract_20261016/);
the `final/` subdirectory contains the artifacts regenerated against the final
application source. Earlier logs and preliminary artifacts are retained.
The final selected regression suite passed **207 tests in 24.08 seconds**,
including all 39 structural/output-contract checks, the full bilingual school
journey, diagnostic navigation, existing heatmap mathematics and Help/locale
contracts. The earlier 168-test run preceded the stronger widget-return check;
its passes are not added to the final total. AST comparison with the retained
source confirms that the only changed application functions are the diagnostic
menu renderer and `main`; estimation and diagnostic calculation functions are
unchanged. The Help-library text changes are outside those function bodies.
The fixture checks also reject missing/nonboolean status values, inconsistent
rank/contrast evidence, mixed models or cases, empty/unknown/duplicate targets,
within-school targets presented as school DFF, and MML/PCM/GPCM substitutions.
The [estimability package author's reference](https://search.r-project.org/CRAN/refmans/estimability/html/estimability-package.html)
supports the row-space/null-space definition only; no R execution, cross-engine
parity or Rasch interval validity is inferred from it.

**Remaining work.** The three-state specification has automated acceptance,
not user acceptance. Actual browser checks, first-time-user comprehension,
production data auditing, joint estimation and uncertainty qualification remain
open. Formal G3 and subsequent MML studies remain paused.

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


<a id="gd-browser-review-20261026"></a>
## G/D browser and distribution-file review — 2026-10-26

Status: **PARTIAL — local Chrome flow and one downloaded CSV inspected;
not hosted acceptance or a complete release gate.** The inspected runtime was
`integration/unified-app` plus the retained working changes. The starting
entrypoint SHA-256 was `a36bd6dd130ebb26421e796e5871080034588f58d107c7b8c32a04e789b4158f`.
Use the local record `validation/generated/gd_browser_20261026/verification.json`
for the final source identity and the small follow-up change.

- In actual Chrome, entered the existing synthetic 4-person × 2-rater × 2-item
  control (16 observations), ran JMLE, and opened Report & Export → Review
  checks → Supporting material → All reporting tools → Tables & checks.
  The app warned that this is a small MFRM dataset. This was a deterministic
  display check, not evidence that the input supports substantive MFRM claims.
- English and Japanese G/D displays showed G = 1.000 and Phi = 0.930, retained
  the observed 2 × 2 row, and explained the interval hold. The successfully
  downloaded Japanese-view CSV retained G = 1 and Phi = 40/43 without display
  rounding, 16 analyzed rows, calculation version, input fingerprint, analysis
  ID, variance policy, and `interval_inference_allowed=False`. It had no interval
  columns and recorded zero bootstrap attempts.
- At 768 × 900, Japanese captions wrapped within the page. Browser document
  and body widths were both 768 px; wide tables used their own scrolling.
  The viewport was restored. This is not a keyboard/screen-reader or mobile
  accessibility acceptance result.
- Removing one cell hid the old analysis until explicit re-analysis. The
  15-observation result refused G/D with the Japanese missing-combination
  explanation, instead of reusing the complete-design coefficients.
- The G/D download action triggered a page rerun and returned the advanced
  report to its first tab. Added `on_click="ignore"`, matching existing document
  downloads. Six affected CSV and bilingual panel AppTests passed afterward.
  Other calculation functions were unchanged, so prior numerical tests and
  studies were not repeated.
- Chrome subsequently showed the missing-case CSV transfer as **blocked by the
  organization**. No bypass was attempted. Actual browser acquisition of that
  file, and browser confirmation of the final no-rerun download behavior,
  remain open. The prior CSV tests cover its metadata/reason and absence of
  unqualified intervals; those tests are not browser-download evidence.
  A download-event tool call also returned after an unusually long wait; that
  control-tool duration is not an app performance measurement.

A temporary copy of the tracked files plus the reviewed untracked additions
passed `--doctor` and `--release-check --json` using local Python 3.14.3.
The copy included the required `reliability.py`, `report_markup.py` and
`reporting_facts.py`; it excluded generated evidence, retained source snapshots,
local AGENTS instructions and secrets. This was a file-identified integration
snapshot, not a committed clean checkout. It preceded the download-only change;
its import/dependency evidence remains applicable, while the final renderer was
checked by the six affected tests. `release-check` marks the presence of tests
and a CI workflow as Ready: it does not run the tests or attest a successful CI
run. Python 3.11/3.12 CI, the full journey matrix and hosted verification remain
separate release requirements.

The operation confirmed a product limitation already in the roadmap: basic
G/D is buried in advanced reporting and currently requires an MFRM fit. Preserve
its explicit target, design limits and self-describing output when adding the
planned independent entry. Do not use this usability finding to expand the
supported statistical designs or duplicate the engine for Tkinter.


<a id="report-selector-followup-20261026"></a>

## Advanced reporting selection and independent D-study follow-up — 2026-10-26

Status at implementation: **PARTIAL — local behavior checks passed; browser and human acceptance were untested.** The later [bounded Chrome follow-up](#report-selector-browser-20261026) adds browser evidence; human and hosted acceptance remain open.

`show_report_section()` previously constructed three groups of nested tabs and
called all six Essential / nine full-report renderers on every page rerun,
including unselected tables and document generation. The earlier Chrome review
also observed loss of the advanced tab selection on language change. The function
now uses one localized content selector and invokes one renderer. A separate
session value preserves the choice when the widget is temporarily absent; a
choice unavailable in Essential mode falls back to APA Report. Existing report
renderers and numerical definitions are unchanged. The ordinary guided export
entry remains unchanged; this is a bounded repair of its optional advanced view,
not a completed redesign or an independent G/D entry.

The README's first reporting journey now has three steps. Its 23 technical and
writing resources remain in a collapsed, optional reference list. A bilingual
Help banner was updated to remove the obsolete nested-tab directions. Current
interpretation holds remain applicable; fewer visible controls do not qualify
estimates or intervals.

AppTest checks exercised every available content choice, Essential/full/forced
full modes, switching English to Japanese while reading Tables, leaving and
returning, and an unavailable advanced selection after changing to Essential.
They asserted one renderer call, unchanged result identity, and no MFRM refit.
The actual equivalence renderer still displays its hold when explicitly selected.
The initial affected run passed 56 tests and failed two older tests that expected
the equivalence panel to render without selection. After changing those tests to
select the panel, both passed. A final two-test locale check passed after updating
the Help directions. This is combined evidence, not a second complete run.
Logs and the before/after source are in
`validation/generated/gd_report_followup_20261026/`. No browser timing, hosted
acceptance, human comprehension, CI or publishing was performed in this follow-up.

The statistical follow-up used an unmodified archived October 21 mfrmr 0.2.4
source package. Its G-study function body and formals match the candidate used in
the preceding comparison. Its D-study directly uses the raw component matrices;
the current candidate instead calls the missing `.mfrm_mvds_components` helper.
The archived function is a historical independent reference, not a repair or
qualification of the current R candidate.

For three existing inputs with nonnegative components, all 36 common future
plans (72 G/Phi values) agree with saved Python output within the prespecified
absolute tolerance 1e-12; the maximum difference is 9.992e-16. These are correlated
projections of three small fixed inputs, not 36 independent datasets. The four
plans of a negative-component input reproduce the expected treatment difference:
with two raters, archived raw Phi is 0.9677419355 and Python zero-clipped Phi is
0.9375. No component matrices were patched to force agreement. The current R
D-study, intervals, coverage and MFRM/ConQuest equivalence remain unverified here.
The retained Python values were reused; no new Python G/D calculation or MFRM fit
was required. See the [reference scope and results](../validation/generated/gd_report_followup_20261026/reference-verification.json).


<a id="report-selector-browser-20261026"></a>

## Advanced reporting selector in Chrome — 2026-10-26

Status: **PARTIAL — the selected local browser operations passed; beginner,
hosted and complete release acceptance remain open.** This follow-up changed
documentation only. It tested the uncommitted integration source identified in
`validation/generated/report_browser_20261026/baseline.json`, using local
Python 3.14.3 and Streamlit 1.66.0 at loopback port 8518. It was not a clean,
committed release or the publicly deployed app.

The existing deterministic 16-row `gd_probe_complete_control.csv` has four
persons, two raters and two items. It is a G/D display control, not a realistic
writing-assessment dataset or a satisfactory MFRM calibration example. The
operator pasted that input and explicitly clicked Run once with Guided RSM/JMLE
defaults to reach the current report route. The displayed stopping status and
reliability caution were retained. No independent numerical comparison or study
was rerun, and this browser exercise does not qualify the fitted estimates.

| Operation | Observed result | Evidence boundary |
| --- | --- | --- |
| Open Report & Export, Review checks, then All reporting tools | One content selector exposed the nine full-report choices; selecting Tables displayed its contents | Confirms the real guided route, not just the isolated AppTest harness |
| Read the G/D panel in Tables | Observed-design G = 1.000 and Phi = 0.930; the 2 × 2 row was marked observed; scope, raw/clipped variance distinction and interval hold were visible | Display check of this control only; not MFRM reliability or interval validation |
| Switch English to Japanese | Tables remained selected as 表; the same rounded G/Phi values and the Japanese interval hold remained visible | No exact-value or AnalysisID comparison was instrumented in this browser session; prior AppTests cover result identity and refit calls |
| Go to 報告書を作る and back to 確認事項 | Both すべての報告ツール（詳細） and 表 were retained | One concrete leave/return route; not every navigation path or a new session |
| Use two Down keys and Enter from 表 | ファセットの同等性 was selected; the all-method hold and its explanation appeared | A keyboard operation, not a complete accessibility audit |
| Inspect the final viewport | Selector, Japanese heading and hold were readable without clipping in the captured viewport | Screenshot inspected in the tool transcript, not saved as a screenshot artifact; no new narrow-width check |

The native-app capture initially failed; a subsequently available Chrome
extension connection completed the pending checks. A radio-input `check()` did
not change language; the visible language label did. These control-tool events
are not app timing measurements. The final browser console query returned no
captured error entries. The server did log a Streamlit warning about setting both
the default and session value of `guided_export_task`; the tested return route
worked, but this is not a claim of a warning-free run. The extension did not
support content-file export. Browser observations are recorded in
`validation/generated/report_browser_20261026/observations.json`; raw UI traces
and the screenshot remain in the conversation. No report/CSV download was
attempted in this follow-up, so the earlier download limitation remains open.

### Findings that constrain the next UX1 changes

1. **The readiness message and allowed action disagree.** With this input the
   banner said “Not ready to run” and two rows said “Fix before running”, while
   Run was enabled and completed. In `build_readiness_report`, fewer than 100
   rows and fewer than 10 persons become `issue`; the observation explanation
   instead mentions about 30 rows. In `run_facets_mode`, the guided button's
   disabled state uses the resource preflight, independently of those advisory
   checks. Do not repair this by banning every small dataset or changing a
   threshold without a basis. Separate malformed/unsupported input from
   sparse-data and precision cautions, and make the action and explanation
   agree. Counts alone must not promise stable estimates. Test both a usable
   small teaching input and genuinely invalid input; preserve post-fit holds.
2. **The review list repeats topics and is partly untranslated.** Chrome showed
   Reliability / separation and Reliability and separation, plus two Bias /
   local interaction headings. The same English headings, actions and evidence
   remained under Japanese navigation. `build_publication_gate_summary`
   combines readiness and claim-guide rows and removes only identical whole
   rows; the guided renderer prints them directly. Similar headings therefore
   need not mean identical evidence. Group by the user's question, retaining
   every distinct reason, action and applicable hold; do not drop rows just by
   their names. Localize the display while preserving technical identifiers in
   the evidence export. No grouping implementation was attempted here.
3. **Tables still mixes reporting and new computation.** Besides report tables,
   the selected panel exposes model comparison, G/D, network analyses and a
   parameter-recovery simulator. One selector fixes nested tabs, but does not
   finish the information architecture. Keep the short basic report path;
   place optional computations behind explicit purpose-specific entry points.
   The already planned independent G/D entry must not require a preceding
   MFRM fit. This finding does not authorize a new simulation campaign or
   new statistical support.

The Python app's useful distinction is the connection between results, readable
explanations, figures and bounded report wording. Repeating the same evidence
across internal checklists can obscure that advantage. The next simplification
should reduce repeated presentation while preserving the evidence and limits.
Web and Tkinter should share that interpretation contract; desktop batch/export
capabilities are not a reason to copy the entire advanced panel into the first
screen or relax inference holds.


<a id="input-readiness-actions-20261026"></a>

## Input readiness and allowed actions — 2026-10-26

Status: **LOCAL FIX VERIFIED; hosted and human acceptance remain open.** The
first finding in the preceding browser review is addressed. This is a bounded
input-UX change, not a change to the estimator, inference qualification or
supported statistical models.

`build_readiness_report` now reserves `issue` for known input problems. Small
row/person counts, a high missing percentage, limited per-person combinations
and constant-score pattern screens request review without preventing a run.
The existing 300-row / 30-person review thresholds remain UI reminders only;
the former 100-row / 10-person prohibitive classification and claims of
sufficient size or inevitably weak reliability were removed. No sample-size
adequacy or precision threshold has been statistically validated by this fix.
The counts are input rows and distinct input person identifiers; the separate
row audit documents which responses actually enter estimation.

Both Guided and Advanced run buttons use the same report, and a run-time check
also covers internal rerun/quickstart triggers. It stops known column-role
errors (including weight), duplicate selected column names, invalid score or
weight values, no retained rows, and fewer than two retained score categories.
The last check uses the same accepted integer conversion as preprocessing.
Existing estimator validation still applies to other settings and conditions;
the readiness panel is not a complete mathematical-identification test.
Missing rows can be excluded without declaring the entire input invalid.
Post-fit interpretation and uncertainty holds remain unchanged.

The person-count action now asks users to review coverage and needed precision
instead of remapping a valid Person column. New count messages and the revised
banner/actions are localized in English and Japanese; the existing detail
table is not yet wholly localized. In Chrome the next-action table clipped
the instruction horizontally, so warning/error summaries now use wrapping
prose. The detailed checks and row-audit downloads remain available.

Validation records are in `validation/generated/input_readiness_20261026/`:

- Initial targeted suite: **59 passed**, covering readiness, rendering, the
  actual main-app controls, input-row audits and outlier screens.
- After adding the near-integer constant-score boundary and making the main-app
  estimator spy independent of the cached namespace: **11 passed**, including
  the singleton-facet compatibility check. Input blockers were also checked
  against the unchanged `prepare_mfrm_data` rejection; valid small, mostly
  missing, single-combination and constant-person-pattern examples were accepted
  by that preprocessor. These checks did not fit a new statistical study.
- The selected integration run initially had **2 passed, 1 failed**. Both
  bilingual invalid-mapping → correction → actual-fit paths passed. The normal
  input-to-export journey reached HTML creation but its test then assumed that
  the first table had an Estimate column. The current first table is observed
  design; the preceding column-schema comparison already agreed with it. The
  assertion was corrected to compare the actual first-table columns and values,
  as the Word check already did; no report generator was changed.
- After the wrapping-summary change and test correction, **10 passed**:
  rendering, input disclosure and one complete ordinary-JMLE journey through
  Help/Home, language switches, PDF/Word/HTML/CSV, result identity, stale-input
  protection and zero refits during those navigation/export operations. These
  counts overlap; they are separate logs, not a new all-suite pass.

Real Chrome checks used the same earlier 16-row G/D control at loopback port
8520. English and Japanese showed two cautions with Run enabled. Replacing one
score with 0.5 showed an input-error explanation and disabled Run in both
Guided and Advanced modes; correcting it re-enabled Run. The next action was
visually readable as a paragraph. A fresh server start, followed by a fresh
page session, confirmed the final English/Japanese copy without the locale
cache retaining an earlier sentence. No browser fit or download was requested;
actual fits and document creation above are AppTest evidence. Screenshots were
inspected in the conversation, not saved as image artifacts. The temporary tab
and server were closed. No release, CI, cross-engine rerun or human review was
performed; report-list duplication, remaining translation and advanced Tables
organization remain the existing UX1 follow-ups.

<a id="report-review-questions-20261026"></a>
## Guided report checks by question — 2026-10-26

The report-review list combined diagnostic-readiness rows with manuscript-claim
rows. Similar titles therefore appeared more than once, but their findings and
actions were not interchangeable. Deleting rows by title would lose evidence;
translating a generic question as though it were an observed finding would also
mislead the reader.

The Guided Review checks view now groups known aliases by the question they
address. Questions, all distinct current statuses and short review guidance are
localized in English and Japanese. The caption distinguishes run-specific
status from general review guidance. Each group retains every original source
row, including its area, status, evidence and action, in one initially collapsed
expander. Details have distinct numbered labels. Technical source text remains
in English where it was originally English; this is not full translation of
all findings or the advanced report panels.

Blocking groups appear first. A caveat and a claim boundary within the same
group both remain visible, rather than collapsing them into a single favorable
label. Unrecognized areas stay as separate additional checks; unrecognized
statuses request source review. The renderer also uses existing qualification
helpers to explain unestablished optimizer termination, current MML holds, or
JMLE identification/facet-boundary holds in the selected language. It does not
infer successful termination from the presence of a numerical result. Bias and
no-bias inference remains withheld, and optional prediction/simulation checks
are explicitly not prerequisites for completing the basic report.

Only the Guided report entry function changed among existing top-level Python
functions; two presentation helpers and a topic map were added. AST comparison
against the turn's starting snapshot confirms that the gate/readiness/claim
builders, estimator and uncertainty functions are unchanged. CSV schemas,
source evidence and document-generation rules are also unchanged. Existing
cross-engine evidence remains applicable within its previously stated scope;
no ConQuest, TAM or mfrmr calculation was rerun for this display change.

Evidence is retained in `validation/generated/report_review_20261026/`:

- Initial targeted run: **38 passed**, covering report-task selection, lazy
  output behavior, external-comparison routing, first-read localization,
  question grouping, mixed statuses, unknown checks and source preservation.
- After making source-detail labels unique: **32 passed**, covering the final
  renderer and current evidence/hold contracts. Actual separated JMLE fits
  retain their report hold in both languages and all original findings in the
  details. The display test rejects estimator calls while changing languages.
  These runs overlap and are not a single 70-test suite.
- Local Chrome used the built-in Writing essay example (960 ratings, 30 people,
  Rater/Task/Criterion) with default RSM/JMLE. Seven source rows were presented
  as four questions: categories, fit/residuals, bias and optional prediction.
  English/Japanese switching retained the selected review task and caveats.
  Opening the bias details showed both the descriptive screen (24 of 120 review
  flags) and the separate prohibition on bias/no-bias inference. This count is
  a technical screen, not an inferential conclusion.
- A final fresh server/session verified the numbered detail labels in both
  languages and the return from Create a report to Review checks. The two
  browser sessions each explicitly ran the example once; no additional
  statistical validation study was conducted. Screenshots were inspected in
  the conversation, not saved as image artifacts.

The first loopback server start was blocked by the filesystem/network sandbox;
a permitted loopback-only start supplied the browser checks. A direct Enter
attempt on the detail text did not open it; a mouse click did. Keyboard
acceptance is therefore not claimed. This work did not verify narrow-window
layout, actual downloaded files, real novice comprehension or a public
deployment. The browser tab and temporary server were closed afterward.
Existing Tables organization and untranslated technical details remain UX1
follow-ups; inference holds and the shared-engine desktop plan are unchanged.

<a id="optional-report-tools-20261026"></a>
## Separate current tables from optional calculations — 2026-10-26

The advanced Tables view displayed six unrelated tool panels before much of the
current analysis: model comparison, G/D, design networks, paired rater-score
differences, within-rater criterion associations and recovery simulation.
Some had execution buttons and G/D computed point results while rendering.
Fixed-SD MML convergence details also exposed a new numerical-check execution
button. These controls made reading a report look like a sequence of additional
analyses to complete.

The seven tools now share **Additional analyses (optional)**. It is available
directly through Guided **Review checks → Supporting material**, and through the
advanced report content selector. The initial choice is empty; only the selected
tool renders. English/Japanese question labels explain what each tool examines,
what it computes and the limits of the result. G/D point calculations run on
selection; all other new calculations still require their execution button.
G/D remains based on retained scores and its supported crossed-ANOVA design,
not on MFRM parameter estimates. This layout change does not implement the
separately planned standalone G/D input workflow.

Tables retains current estimation/fit/reliability/measure results, qualification,
APA formatting and existing numerical-check evidence. For fixed-SD MML it hides
only the execution button and points to the optional entry; the cached gradient,
records and hold notices remain. The Results/diagnostics route keeps its existing
explicit numerical-check button. Statistical calculation, qualification and
export definitions are unchanged. The display signatures have an optional
`allow_check` argument whose default preserves existing diagnostic callers.

Conditional display removes hidden Streamlit widget keys. Recovery settings
are therefore retained separately and restored before showing its controls,
so returning to a panel can retrieve the cache corresponding to the chosen
settings. Tool selection also survives language changes and leaving the view.
This does not fit new data or reclassify a cached result as qualified.

Evidence is under `validation/generated/report_optional_tools_20261026/`:

- Initial focused run: **61 passed, 3 failed**. The failures arose while AppTest
  serialized a translated selectbox: its formatting callback consulted session
  language outside the rendering context and could produce a label absent from
  the displayed options. Labels are now captured at render time, as in the
  existing report selector. This was a UI/verification failure, not a numerical
  test failure.
- After that correction: **8 passed**, covering every tool route, no automatic
  refit, current tables without optional launchers, bilingual fixed-SD evidence
  with/without an execution button, retained simulation settings/cache, and the
  direct Guided route with a global interpretation hold. The simulation test
  substitutes a small result for its explicitly clicked calculation; it does
  not perform a new Monte Carlo study.
- Existing returned-point UI checks: **3 passed, 2 deselected**, confirming
  explicit execution, no optimizer/native-engine refit, retained language/state,
  invalidation after input changes, and distinct widget keys across contexts.
  These small checks use their existing fixed-SD MML fixture; they do not expand
  statistical validation or repeat a cross-engine comparison. Test counts from
  these runs overlap and must not be added as an all-suite result.

Local Chrome used the built-in Writing essay example (960 ratings, 30 persons,
Rater/Task/Criterion), fitted once with default RSM/JMLE. The direct additional
entry initially showed no selected tool or run button. Selecting simulation
opened its settings only; changing replicates to 4 and switching English →
Japanese retained the choice and value. Tables retained result/formatting and
convergence content without the optional analysis launchers. Reopening the same
simulation from the advanced content selector restored 4 replicates. G/D
selection correctly reported that this three-facet example is outside its
one/two-random-facet scope and retained the interval hold.

No model-comparison or recovery-simulation execution button was clicked in the
browser. Screenshots were inspected in the conversation, not saved as files;
the temporary tab and loopback server were closed. Downloaded files, narrow
screens, keyboard acceptance, novice comprehension and public deployment were
not checked in this turn. Some existing technical labels (including recovery
setting names) remain English in Japanese mode. These are remaining UX work,
not a reason to relax inference restrictions. No commit, push, deployment or
new statistical study was performed.

<a id="optional-tool-copy-20261027"></a>
## Explain optional-tool setup in Japanese — 2026-10-27

The recovery screen still used English setup labels in Japanese mode. More
importantly, the replicate count did not explain the distinction between new
response datasets and iterations within one estimation. The default ten
replicates could also be mistaken for a precision recommendation.

All eight setup labels are now Japanese. Four tooltips explain the shared
generator/fit model, the simulation's separately chosen estimator, replicate
counts, and the seed needed for reproduction. A bilingual preview shows the
planned persons × raters × criteria rating count, the number of generate/fit
replicates, the score range starting at zero, and the complete rating design.
It explicitly refers to synthetic observations separate from the current fit.
The default remains ten, with no precision guarantee. The setup controls align
at the bottom when translated labels wrap.

Japanese G/D headings and network entry headings, descriptions, execution
buttons and help were clarified. Paired observed score differences are
distinguished from MFRM severity estimates and bias tests; criterion
associations are not treated as a causal halo diagnosis. Existing halo source
references are retained. This is entry/control copy, not full translation of
all technical result columns or dynamic error messages.

Evidence is in `validation/generated/optional_tool_copy_20261027/`:

- **14 passed** initially: optional-tool navigation, settings/cache retention,
  qualification notices, existing help readability and the new input guidance.
- **10 passed** after the final column-alignment change (overlap with the first
  run). The preview test compares an 80-rating, 10-person/2-rater/4-criterion
  setting and scores 0–4 with the existing response generator. It substitutes
  generation only for the simulation evaluator and rejects estimator calls;
  this is not a new recovery study or validation of inferential precision.
- Local Chrome: one explicit RSM/JMLE fit of the built-in Writing essay example
  (960 ratings) to reach the report route. English/Japanese setup labels and
  the 270-rating default preview were checked. The Japanese replicate tooltip
  opened and showed the distinction from optimizer iterations and the absence
  of a precision guarantee. Changing only persons from 30 to 40 updated the
  preview to 360 while the execution button remained available; the simulation
  was not started. A display rerun applied the alignment change while retaining
  the setup and current analysis. Screenshots were inspected in conversation.

Only `render_parameter_recovery_simulation` changed among Python functions;
existing generator, evaluator, estimator, exports and qualification functions
remain identical to the starting snapshot. After browser/test checks, the
preexisting author/year references were retained in the translated halo
caption; this final citation-only copy was checked as text, not rendered again.
No cross-engine comparison or statistical study was rerun. The browser tab and
loopback server were closed. Real novice comprehension, narrow-screen and full
keyboard acceptance, downloaded files, public deployment and the planned
standalone G/D entry remain outside this turn's completed work.

<a id="gd-independent-entry-20261027"></a>
### Observed-score G/D entry without a preceding MFRM fit — 2026-10-27

The welcome page and the sidebar of an existing workspace now expose
**Explore a rating plan / 評定計画を検討する**. This local implementation completes
the bounded input-to-result route for basic G/D; it is not a public deployment
or a desktop implementation. It preserves the application's distinction between
a calculation, its assumptions and permission to make an inferential claim.

The user chooses the small artificial example, pastes a long-format table, or
loads CSV/TSV; selects scores, the object and one or two random facets; describes
the population; and acknowledges equal weights and the absence of fixed facets.
Calculation is explicit. Continuous scores and identifier strings such as `001`
are preserved. No MFRM rounding, category recoding or row exclusion is applied.
A population declaration records the intended use; it does not verify sampling
or validate generalization. Unused columns are disclosed and are not adjusted for.

Both the existing fit-based adapter and `compute_generalizability_from_data`
call `_compute_generalizability_study`. They share the complete balanced crossed
design checks, signed and clipped variance components, G/Phi formulas and interval
hold. Their score bases and analysis identifiers differ. Fixed-facet models,
incomplete designs, duplicate cells, nonfinite scores, fewer than two levels,
more than two random facets and non-unit `Weight` values are not newly supported.
No external R or MFRM estimator runs on entry or on a G/D calculation.

The first results show observed G/Phi and compare a selected complete crossed
plan with the observed design, including relative and absolute **error variances
in squared score units**. All computed plans and variance-decomposition details
are expandable. The existing bounded grid is retained; it is not an unrestricted
sample-size search or a recommendation for a required number of persons.
The downloadable CSV adds those observed error variances and the selected plan
to the existing metadata/components/grid sections. It records score basis,
population declaration, fixed facets (none), data fingerprint, analysis ID,
calculation version and the interval hold. No unqualified CI is exported.

Changing input, column roles, population or scope acknowledgment hides the old
result and its download until the matching calculation exists. Language and Help
preserve the G/D setup and plan. Returning to MFRM preserves its input, settings
and fitted output. While G/D is open, existing MFRM input widgets are kept alive
in a collapsed sidebar container, and the existing auxiliary-view guard returns
before estimation. Uploaded G/D tables are retained within the session; an explicit
clear resets both the retained input and upload widget so the same file can be
loaded again. This is not persistence across browser sessions or project-file save.

Validation evidence is under `validation/generated/gd_entry_20261027/`:

- Initial G/D numerical/boundary suite: **61 passed, 2 failed**, 14.13 s. Both
  failures were test lookups of the old table position after a comparison table
  was inserted. The tests now identify the grid by `IsObservedDesign`.
- Initial new-entry/display batch: **11 passed, 3 failed**, 9.63 s. Test assertions
  incorrectly treated citation downloads as G/D downloads and a zero-valued
  call-count dictionary as an empty list. They were corrected to check the G/D
  download key and actual estimator/refit counts.
- Entry, optional-tool and guide batch: **25 passed, 1 failed**, 40.96 s.
  The guide failure is discussed below; it was not hidden by changing convergence.
- Final targeted entry/display batch: **15 passed**, 8.49 s. This covers equal
  raw/fit calculations for one/two facets, independent hand-derived error terms,
  raw-score identity and preservation, invalid inputs, lazy execution, changed-input
  invalidation, EN/JA/Help/return, selected-plan CSV, a real 12-rating MFRM fit
  retained without re-estimation, paste rejection and uploaded-table retention,
  clear and reload. The upload-widget value is supplied by a test double while
  the production CSV parser and application route run; this is not an OS file-picker
  acceptance check. The selected plan and error terms are checked in the actual
  bytes passed to the download widget, not a browser-saved file.
- Existing R-helper fixture checks: **4 passed**, 2.08 s. These are the bounded
  source-set/dominant-component/coefficient-tolerance/observed-row checks of the
  committed fixture, not a new mfrmr candidate execution or estimator equivalence.
  The earlier direct matched-ANOVA and archived D-study comparisons remain the
  evidence for their original, limited conditions. ConQuest and mfrmr were not rerun.

These are overlapping stages, not a sum of independent successful tests. AST
comparison confirms unchanged MFRM estimation, input qualification, ANOVA
components, D-grid and bootstrap functions. The one-facet coefficient helper only
adds already-computed relative/absolute error values to its details. Existing
MFRM input/render functions differ only in the sidebar container argument.
The calculation version remains `balanced_anova_v2`; the CSV schema is extended.

Chrome checked the English sample route without an MFRM fit (8 ratings,
G = 0.797, Phi = 0.654), Japanese display, selection of a 1 × 1 plan, and a
pasted four-row continuous-score control (G = 0.938, Phi = 0.857; error variances
2 and 5). Removing one of its cells produced the explicit incomplete-design
refusal and no coefficient output. Screenshots were inspected in conversation.
The final Japanese observed-design caption was checked by AppTest/text; the
long-running browser's locale cache still displayed its earlier wording.
Real novice comprehension, full keyboard/narrow-screen acceptance, native upload
and downloaded-file behavior, and public deployment remain unverified.

**Existing guide limitation discovered during regression checking — diagnosis corrected below.**
The initial description treated the summary's `Converged=False` as actual
nonconvergence. The [follow-up](#guide-termination-recovery-20261027) established
that the optimizer reported Boolean success and the summary had been qualified
without that record. The earlier nonconvergence diagnosis is withdrawn; the
observed test failure and log remain valid evidence of the summary defect.

`test_complete_sample_guide_uses_real_fit_and_keeps_analysis_identity` also fails
in isolation. Running the pre-change source snapshot reproduces it: the Writing
sample has `Converged=False`, `InferenceReady=False`, 30 iterations and gradient
norm about 0.04535, so the guide remains at Estimate with
`learning.fit_nonconvergence_reviewed`. This was a pre-existing summary/guide defect, not evidence that the G/D entry
changed numerical convergence. The follow-up corrects the missing-record call
and separately tests genuinely unsuccessful attempts without relaxing any
statistical criterion.


<a id="guide-termination-recovery-20261027"></a>
### Correct termination evidence and recoverable guide failures — 2026-10-27

The Writing sample was blocked by an incorrectly constructed summary, not by
unsuccessful optimizer termination. Before this change, a direct production fit
returned Boolean `opt.success=True` and the optimizer message
`CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH`. Its convergence table was
True, and the current termination predicate was True, but its summary was
`Converged=False`, `SavedConverged=True`, `InferenceReady=False` with
`stat.estimation.optimizer_termination_unverified`. The prior diagnosis based on
the summary alone is therefore corrected. Optimizer success still establishes
only that the configured stopping rule was met, not scientific validity or
sampling-uncertainty calibration.

`mfrm_estimate` had called `qualified_estimation_summary` with just summary and
config, before assembling the result containing the optimizer and identification
audit. The fail-closed predicate correctly rejected the incomplete object. The
call now occurs after the complete result is assembled. This fixes fresh JMLE
and MML summaries through their common entry; MML inference holds still apply
even when optimizer termination succeeds. The optimizer, iteration limits,
tolerances, parameterization, likelihoods and uncertainty formulas are unchanged.
An AST comparison after removing the moved summary assignment finds the estimator
function otherwise identical. The shared qualification predicate and summary
qualifier themselves are unchanged.

Guide advancement now consults explicit Boolean optimizer evidence instead of
trusting a summary cell or interpreting a truthy string. This controls the
learning stage, not inference permission. Successful calculations proceed to the
existing evidence/overclaim exercise and bounded learning record. On genuine
unsuccessful termination the guide remains at Estimate, explains the hold and
can open existing results and settings for review. The user can then make a
justified change and explicitly run again. Display reruns, language changes and
Help do not retry estimation. An exception remains a failed attempt; its existing
bounded problem notice, without raw exception text, is retained for subsequent
language changes. The generic fallback list that suggested tolerance relaxation
without a diagnosis was removed in favor of the existing classified problem notice.
Existing pattern-specific error remedies were not comprehensively rewritten.

Evidence under `validation/generated/guide_termination_20261027/`:

- `before_probe.log` records the optimizer/summary contradiction before editing.
  The earlier baseline failure was reused rather than rerunning the full prior
  studies. No ConQuest or mfrmr execution was needed for this record-passing fix.
- Initial guide, termination and state-machine batch: **57 passed, 3 failed**,
  49.58 s. The regular real-fit guide now completed. All three failures were old
  expectations in MML tests that never call the changed estimator: one expected
  saved caution prose to be copied; two accepted a summary-only convergence flag
  without an optimizer record. Those assertions were aligned with the existing
  fail-closed contract, retaining original inputs, current method holds and
  idempotence checks. No qualification code was changed to make them pass.
- New result/guide cases, those corrected tests and analysis journeys:
  **21 passed, 3 failed**, 85.50 s. New numerical cases cover successful and
  iteration-limited JMLE plus fixed/free-SD MML; Boolean and malformed stopping
  records are distinguished. Four actual JMLE/MML input-to-report journeys passed.
  Two new recovery tests looked for a run button that is absent in the existing
  post-fit settings view; they now use its actual sidebar Run control. The third
  failure was an older Word-table assertion assuming wholly withheld columns
  always exist; it now accepts their omission or blank cells, while still checking
  the hold note, original snapshot and zero refits.
- Final targeted recovery and restored-output batch: **4 passed**, 20.26 s.
  In English/Japanese the first attempt was deliberately limited to one optimizer
  iteration, remained held, survived language/Help/review without another fit,
  and progressed only after explicit retry. The exception path retained only the
  bounded notice and supported explicit retry. Restored free-SD MML output kept
  its hold through stale settings and Word export, with zero automatic refits.

These are overlapping stages, not a single sum of independent tests. This work
does not requalify old saved inference flags by declaring all historical runs
valid, and it does not establish new numerical or interval accuracy. Public
release, full keyboard/narrow-screen acceptance and novice comprehension remain
separate from these local checks.

Chrome checked the Japanese Writing sample (960 ratings, 30 people, three
facets). It reached evidence review; the overclaim answer was rejected without
advancing the guide, and the bounded answer opened the final learning record.
The final record's dataframe clipped the end of its limitation at the observed
1200 px viewport. Its same three contents now appear as wrapping paragraphs.
Reloading the edited source left that development session outside the active
guide, so the final visual check restarted the guide with a second explicit fit;
this is not counted as one uninterrupted run or as evidence of hot-reload
continuity. The final source reached the learning record with the entire
limitation visible, then the explicit Finish button showed guide completion and
the ordinary result. The residual-PCA caution (first eigenvalue 3.74) remained.
Screenshots were inspected in conversation. The final paragraph-only change was
checked by browser and syntax validation; the completed numerical batches were
not rerun for it. Genuine termination failure and exception recovery were tested
by AppTest, not injected into this browser session. The task's tab and loopback
server were closed after checking.

<a id="report-review-notes-20261027"></a>
### Carry screen cautions into basic reports — 2026-10-27

**Finding.** The Japanese Writing sample reached Report & Export with a visible
residual-PCA warning: first eigenvalue 3.74 and stability requiring review.
The actual basic PDF and HTML bytes passed to the download control contained
neither `3.74` nor the specific PCA caution. They gave generic advice to add
residual checks. The screen asked users to carry cautions into their report
manually, so a self-contained exported document could lose known limitations.
This is an output-content gap, not evidence of a changed numerical estimate.

**Change.** Basic and detailed PDF, Word and HTML now contain a common
“Computation and interpretation checks” section before their tables. It reuses
`build_first_read_guide_rows(..., localize=False)` to retain the recorded
termination state plus non-OK checks, their evidence and follow-up. The actual
Writing example includes the PCA stability/loading sensitivity and 24 of 120
local-interaction review flags. The text distinguishes meeting a stopping rule
from model validity and explains that skipped optional checks support no claim
in their area; they are not prerequisites for every study. The report remains
English, with English/Japanese interface captions explaining what it includes.
The document cache version changes so a live session cannot reuse pre-fix bytes.
This section carries first-read findings, not every technical audit or a completed
study-specific interpretation. Existing calibration and inference holds remain.

**Verification.** Evidence is under
`validation/generated/report_save_20261027/`:

- The production AppTest route used the built-in Writing example (960 ratings,
  30 people, three facets), ran estimation once, and captured the actual PDF and
  HTML supplied to the download widget. The pre-change PDF had seven pages and
  did not contain the PCA finding. The fitted output was retained locally for
  subsequent output checks, avoiding new estimation.
- `tests.log`: **78 passed in 14.47 s**, covering run facts, all report formats,
  basic/detailed output, cache invalidation, bilingual controls, existing table
  layout and FAQ. New cases check PCA evidence, escaped text, skipped checks,
  unsuccessful/missing/truthy termination, fixed/free-SD MML and a JMLE boundary
  hold. Existing result and diagnostic values remain unmodified.
- `report_verification.json`: actual basic PDF, Word and HTML generated from the
  same saved Writing result all include stopping status, PCA evidence and
  interaction-review context. All four HTML result tables exactly match the
  pre-change tables; the analysis identity is unchanged and no refit occurs.
- PDF visual review found the Results tables heading orphaned on the preceding
  page. Keeping section headings with subsequent content fixed it. This final
  layout-only change was checked by rebuilding the PDF, confirming the heading
  and Table 1 share a page, and inspecting pages 1–2. Final pages 3–8 are byte-for-
  byte identical rendered PNGs to the six pages already visually reviewed.
  `visual_qa.json` records this reuse. Numerical tests were not rerun for pagination.
- Single local AppTest timings were 1.55 s for the fit action and 13.49 s for the
  first PDF/report action (with Matplotlib font-cache setup); selecting HTML took
  0.99 s. These are individual observations in that environment, not a speed
  benchmark or a measured regression. The final PDF contains eight pages,
  including its explicit review context, four tables, two figures, Q&A and
  references.

**Browser boundary.** Chrome reached the Japanese report screen and its PDF
control. The download click was performed, but completion of the saved file was
not established. Opening `chrome://downloads/` was rejected by the Browser Use
URL policy; no alternative browser surface or policy workaround was used.
The output inspection above independently captures the app's generated bytes
and does not establish that a native browser download succeeded. Public release,
real novice comprehension, full keyboard/narrow-screen acceptance and native
file-save completion remain unverified. No ConQuest/mfrmr study was rerun and no
inference qualification was expanded.

The final Chrome check used a fresh server/session because the earlier live
session retained the old translation caption after source reload. The final
Japanese caption explicitly states that the report includes stopping status and
first-read checks needing review, with the English report language and sharing
notice retained. Each of the two browser sessions used one explicit sample fit;
the AppTest capture used a separate fit. These are bounded UI/output checks, not
new simulation replications.

The final browser also switched among PDF, Word and HTML and exposed each
matching download control. Both task tabs and local servers were closed. This
does not alter the native file-save limitation above.

<a id="first-read-critical-review-20261027"></a>
### Critical review of shared first-read judgments — 2026-10-27

**Review-time verdict: fix these presentation/qualification inconsistencies before release.**
The authorized local fix and its checks are recorded [below](#first-read-qualification-fix-20261027); the following preserves the pre-fix findings.
The preceding change correctly carries the ordinary Writing sample's cautions
into its report. Sharing that function also shares its incorrect treatment of
some missing or withheld evidence. The review covered the recent document path,
its shared first-read builder, qualification helpers, actual diagnostic return
shapes and relevant tests. It assumes one analysis per Streamlit session; this
was not a whole-repository, security, concurrency or numerical-accuracy audit.
No application or test code was changed by this review.

1. **Unverified saved reliability becomes OK and supports ordering language.**
   `build_first_read_guide_rows` (lines 29726–29749 in the reviewed source) calls
   `qualified_reliability_table`, then compares the returned Person value against
   numeric thresholds without checking `ReliabilityStatus`. A three-person
   SE-based fixture has a valid recalculated reliability of 0.96. Giving its saved
   table 0.9999 makes qualification correctly return `unverified` and
   `stat.reliability.saved_values_mismatch`, while preserving the original number.
   First Read nevertheless returns OK, displays 1.00 and recommends using person
   ordering claims. The publication reliability table correctly withholds its
   precision cells. The newly shared report-note filter drops the OK row, leaving
   inconsistent summary and detailed evidence. Use the per-facet qualification
   before interpreting the value, retain its reason, and preserve genuinely valid
   neighboring facets. Recheck both this mismatch and valid known-basis controls.

2. **Unavailable PCA disappears; absent stability evidence can become OK.**
   The real PCA helpers return `None` and a reason for an observed residual matrix
   with only one item-combination column. Production diagnostics keep
   `pca_enabled=True` in that situation. First Read (lines 29776–29813) handles
   disabled PCA and a populated result, but no enabled/unavailable branch, so the
   Dimensionality row disappears from both the action plan and the report notes.
   The full readiness builder independently retains `NOT_ASSESSABLE` and the
   concrete reason. A second saved-result fixture with eigenvalue 1.5 and no
   stability table returns OK and “stable enough for scoped reporting”, because
   an empty status is merely tested as unequal to Review. Keep disabled,
   unavailable and computed evidence distinct. Missing audit evidence must not
   imply stability. Reuse the stored reason so a column-support failure directs
   the user to the data/diagnostic limitation, not just an already-enabled option.

3. **Zero bias flags can hide inference holds and skipped pairs.**
   First Read (lines 29824–29854) aggregates counts and a disconnectedness flag;
   it does not retain `ClaimStatus`, `PublicConclusionAllowed` or a skipped pair.
   The first probe used old-schema cells and correctly obtained Review, because
   their scalar-reference provenance was unavailable. That was not evidence for
   the hypothesized zero-flag defect. The follow-up selected two real, current-
   schema cells from the saved Writing output as a bounded presentation fixture.
   Both have no screening flag; the real audit still returns Review, Do not claim
   and `PublicConclusionAllowed=False`. First Read instead returns OK and the
   report's current-check section omits bias. Adding a separate skipped pair
   leaves exactly the same OK summary. This is not a new fitted study or a claim
   that the complete Writing family has no flags. Preserve each pair's missing
   evidence and the inference hold separately from descriptive flag counts.
   No-flags, not-computed and not-qualified are different states.

The report filter at lines 37428–37434 includes convergence and non-OK rows, so
these upstream defects reach the newly added notes. Keeping every row would not
repair an incorrect OK interpretation; the small fix belongs in the existing
shared first-read judgments. No new status framework, engine or model is needed.
The actual basic HTML confirms the omissions. Reliability-table withholding and
the later bias FAQ hold remain present: this review does **not** claim that all
output protection is bypassed or that bias inference was numerically requalified.

Evidence: `validation/generated/report_readiness_critical_20261027/` retains
source hashes, the standalone probes, JSON results and logs, plus three actual
basic HTML outputs. Figure generation was stubbed out in these semantic-output
checks; no new visual-layout claim is made. No estimator was called. The saved
ordinary Writing result remains a positive control with PCA eigenvalue 3.74 and
its 24/120 interaction review flags retained. The earlier 78-test success remains
valid for its scope; it did not establish these boundary cases. Re-running that
successful batch, simulations or external engines would not resolve the missing
qualification checks. Browser file-save completion, novice comprehension and
public deployment remain outside this review's evidence.


<a id="first-read-qualification-fix-20261027"></a>
### Preserve first-read qualification in screens and reports — 2026-10-27

The three preceding findings are fixed locally. The change uses the existing
qualification and diagnostic records; it does not grant new inference permission.
It stays within UX1's priority of correcting misleading current-result output.

- Reliability only uses rows qualified as `diagnostic` when interpreting values.
  Unverified, unassessable or withheld facets retain a readable reason (saved-value
  mismatch, incomplete inputs, unknown uncertainty basis, etc.). A valid neighboring
  facet remains usable. Missing tables are explicitly unassessed. Saved numbers
  are not replaced, and calibration holds still take precedence.
- Requested but unavailable PCA remains in the summary with its recorded reason.
  Disabled PCA remains skipped. Empty/non-finite eigenvalues cannot form a screen;
  a finite eigenvalue with no affirmative stability record cannot be OK. The final
  readiness builder uses the same affirmative status requirement. Exported reason
  payloads resolve in English independently of the sidebar, including embedded
  item-combination labels; legacy plain-string reasons remain as recorded.
- Bias retains its inference hold separately from flag counts. Unassessed pair
  names/reasons and cells without usable inference remain visible even when the
  computed cells have zero flags. This is a restriction on bias claims, not a
  new hard stop on all other qualified analysis results.
- The shared First Read path feeds both guided/advanced screens and existing basic
  and detailed PDF/Word/HTML notes. Document cache versioning prevents reuse of
  reports containing the old summaries. No new settings, model or status framework
  was added.

Validation and scope:

1. The targeted batch passed **164 tests in 45.24 s** under Python 3.14 / Streamlit
   1.66.0: `test_first_read_qualification`, `test_reliability_input_contract`,
   `test_report_run_facts`, `test_report_export_ux`, `test_final_readiness_evidence`,
   `test_calibration_termination_hold`, and `test_review_guidance`.
   The new 30 checks cover the three failures, valid controls, mixed valid/held
   facets, unknown stability, unavailable values, all six format/detail combinations,
   and the production First Read renderer in English/Japanese with fitting forbidden.
   Existing termination and method-specific holds also remain covered.
2. The initial new-test run had 25 passes and five assertion failures: four tests
   incorrectly rejected the valid wording “Person reliability is unavailable”,
   and the UI test read prose elements without the rendered detail DataFrame.
   Assertions were corrected to reject numeric endorsement and inspect the actual
   table as well; no qualification was relaxed to make these checks pass.
   The first shell attempt used a pre-existing broken `.venv` Python 3.12 link;
   tests ran with the available Python 3.14 used by the preceding work.
3. Reused the retained artificial Writing result
   `ana_0d4fb0699d783f6031c3416f` and its two existing figure images. Generated all
   three basic report formats with fitting prohibited. Existing four HTML tables
   are equal to the preceding output and a pickle digest confirms the complete
   saved result was not mutated. PCA 3.74 and its stability caution remain; the
   explicit bias inference hold is now also in the notes. Replaying the actual
   two-cell, zero-flag presentation probe now produces Review and retains the
   skipped `Rater x Task` reason. This is not an inference test of a reduced family.
4. Inspected all eight PDF pages as images. No clipping, overlapping text or
   orphaned table heading was seen. The first reuse harness accidentally repeated
   figure numbering already present in saved HTML captions; only that harness was
   corrected. The corrected figure page was inspected again, and extracted text
   outside that page was identical. This was not an application renderer defect.
   Format regression tests use small controlled inputs and omit unrelated figures
   and appendices; this saved-result check retains real figures and report tables.
5. `py_compile` and `git diff --check` passed. AST comparison against the turn's
   saved source finds changes only in `_render_reason`, `build_first_read_guide_rows`,
   `_render_publication_document_section` and `build_final_report_readiness`.
   Reliability and output-qualification modules are byte-for-byte unchanged.
   Earlier unrelated working-tree edits were preserved.

Evidence: `validation/generated/first_read_fix_20261027/` retains the pre-change
files/hashes, targeted and initial logs, the saved-output probe, generated reports,
PDF renders/text and checks. Existing numerical evidence is reused within its
original scope; no external ConQuest/mfrmr run, new numerical study, Tkinter
implementation, commit, push or deployment was performed. Browser-native save
completion, novice comprehension and public-environment acceptance are still
unverified. AppTest rendering and generated bytes do not establish those outcomes.


<a id="api-connectivity-review-20261027"></a>
### API connectivity: analysis ownership, persistence and reporting — 2026-10-27

**Review-stage verdict: four active-workflow handoff defects needed correction before release.**
The bounded follow-up below fixes these reproduced cases; this section preserves
the original evidence and scope.
Existing computation and qualification APIs are substantial, but not all consumers
use the same analysis identity, retained result or screening settings. Fix these
connections within UX1 before exposing more functions in the beginner interface.
This review changed no application or test-suite code and did not execute new
MFRM fits, external engines or numerical-validation studies.

Scope and method: one user's Streamlit session, including switching results,
restoring history, changing language and exporting. Statically inventoried the
entrypoint plus 78 package files (1,052 public-named top-level functions/classes).
The 138 definitions with no static name/attribute load are **candidates, not dead
APIs**: string-based core lookup, re-exports, tests, validation drivers and explicit
release exclusions change that classification. Read the main routes and followed
high-impact consumers through input/result identity, optional analyses, history,
DIF/bias screens and report/CSV generation. This is not a complete semantic review
of all 1,052 definitions, a security audit, or a concurrency/load assessment.

1. **High: additional analyses do not use the existing complete analysis identity.**
   `_model_choice_cache_key` (`streamlit_app.py:49655`) looks for `params['par']`.
   The retained real Writing fit's actual `params` keys are `theta`, `population`,
   `facets`, `steps`, `steps_mat`, `log_slopes` and `slopes`; `par` is absent.
   The key therefore ends in `0x0` and otherwise describes only model settings.
   Different data with the same settings receive the same key even though
   `build_result_analysis_identity` (line 40880) correctly distinguishes them.

   In the production design-network renderer (line 50066), first compute the
   retained 30-person design, then supply its 10-person subset under the same
   settings. The screen reuses 40 nodes / 332 edges without another calculation;
   a direct call to the same network API returns 20 nodes / 132 edges for the
   second design. Both records have different full analysis IDs. No MFRM fit is
   required for this observed-design diagnostic. The same key is used by model
   comparison (49716), severity/halo networks (49899/49989), assignment sensitivity
   (48971) and its export ownership checks (64371/65077). Those other consumers
   were traced in source, not independently rerun as numerical studies.

   Small correction: reuse the existing full analysis identity, plus the optional
   analysis's own settings/version, for cache and export ownership. Do not create
   another fingerprint scheme. Acceptance: changed data, estimates, anchors or
   relevant settings cannot reuse the old bundle; reopening the identical result
   still reuses it, without an automatic refit. Current cache tests distinguish
   population-SD options but do not test different retained data under equal model
   settings (`tests/test_information_criteria_units.py:107`).

2. **High: history restoration and general exports do not check ownership of prediction/refit results.**
   `restore_run_from_history` (34477) restores only `facets_mode_output` and its
   notice. The new-fit path clears prediction/refit session keys (39708), but
   restoration does not. `collect_download_frames` (64354–64365) accepts the current
   session's prediction/refit tables without comparing their source analysis with
   the result being exported. The collector's subsequent assignment check exists,
   but inherits finding 1's incomplete key.

   The probe creates an actual new-design prediction from the retained Writing
   fit, restores a distinct retained-data fixture through the real history API,
   and calls the real download collector. `new_design_predictions` still contains
   the previous analysis's table. A marked synthetic refit-summary table follows
   the same path; no costly refit study was started. This proves the collector
   boundary with `public_export_mode=False`; it is not a claim that every document
   format or the default privacy-filtered selection includes these tables.

   Small correction: retain source analysis ID and computation settings with each
   derived bundle and validate them when displaying/exporting. Restore matching
   derived results or explicitly mark them uncomputed; simply clearing one UI
   branch is insufficient for the shared collector. Acceptance: A → B → restore A
   never exports B's derived tables under A, and matching saved tables remain
   usable without recomputation. The existing prior-SD export helper already
   validates its evidence against the current result (27321); reuse that principle.

3. **Medium: classical DIF results exist only during the Run-button event.**
   `analyze_classical_dif` returns tables, settings, an audit and caveats, but
   `show_classical_dif_section` (59550) returns immediately whenever Run is not
   pressed. Its result is a local variable, with no retained bundle for later
   display, language changes, history or a report handoff.

   In AppTest, the production screen analyzed the retained artificial Writing
   rows with alternating artificial external groups and `min_group_n=5`. It
   showed four metrics and the result CSV download. Switching English to Japanese
   left the calculation count at one but removed all result metrics and that
   download. This tests result persistence, not the validity of the DIF statistics.
   The removed external difR-validation download is a separate, intentional
   release boundary and is not the cause of this failure.

   Small correction: retain the returned bundle with data/group/settings identity;
   let Run control calculation, not continued visibility. Language/navigation/save
   must reuse the result; changed scientific inputs must invalidate it. Preserve
   the DIF-specific audit and caveats when adding any report handoff, and keep
   external-group observed-score DIF distinct from fitted-facet bias.

4. **Medium: bias screen settings stop at the local screen.**
   `show_bias_section` passes the selected alpha, minimum cell count and practical
   size to screening/audit APIs (59696–59714). First Read (29872) and general CSV
   export (64433/64443) call those APIs with defaults. The local settings are not
   part of the result's shared reporting inputs.

   On the actual retained `Person x Rater` table, the production UI and its local
   CSV show 120/120 flags at `min_n=100`, while the shared report and general CSV
   show 24/120 under `min_n=5`. Both are correct for their different settings;
   their relationship is not explained by the current screen. The CSV's
   `MultiplicityFamily` field does record each criterion, but the user-facing
   summary does not identify that it reverted to defaults. Inference remains
   withheld throughout: this finding does not reopen bias/no-bias conclusions.

   Small correction: carry the pair-specific screening settings into shared
   summaries/exports, or explicitly designate the local change as a separate
   sensitivity view and identify the reporting settings in both places. The
   acceptance test must compare a non-default selection, local CSV, shared
   summary and general export; equality under defaults alone misses this break.

Other connections and deliberate boundaries:

| API family | Observed connection | Disposition |
| --- | --- | --- |
| Estimation → diagnostics → per-facet qualification → First Read → publication notes | Shared source and the preceding 164-test fix remain unchanged | Preserve this working path; do not duplicate calculations for Tkinter |
| Raw/fitted G/D → variance components → observed/selected D plan → CSV | Common engine and export frame; the recently added independent entry no longer requires an MFRM fit | Connected within its balanced-design and interval-hold limits; no new G/D study run here |
| Prior-SD evidence → current-result validation → export | `_add_mml_prior_sensitivity_export_frames` explicitly checks the evidence contract | Existing ownership pattern to reuse; no new MML inference approval |
| `simulation.planning.build_formula_plan_summary` / `evaluate_plan_structure` → prospective-design UI | Workload/cost/topology APIs exist; no current entrypoint call. Tests and package exports are present | Confirmed unconnected product capability, already scoped as D1. Keep behind the current UX1 priority; do not add it as another required novice screen |
| Individual rater feedback | Fit, category, exposure and network components exist; no single saved-result feedback-sheet contract | Previously planned integration candidate, not a newly demonstrated calculation bug; do not automatically connect conditional bias flags to judgments about a rater |
| CMLE/WLE and confirmatory research modules; Stan/Posterior viewer and external DIF-validation bundle | Research/compatibility references exist; public routes deliberately exclude specific legacy tools | Do not relabel them as accidental dead code or reconnect them automatically. The existing public-route exclusion test passed (1 test, 1.16 s) |
| `guided_goal_step_table`, `guided_progress_checklist` and similar old guidance helpers | Active-tree references found in definitions and unit tests, not current UI consumers | Maintenance candidates after checking compatibility promises; their tests are not evidence that the current user sees them. Do not restore redundant screens solely to create callers |

Evidence: `validation/generated/api_connectivity_review_20261027/` contains the
static candidate inventory, source hashes, three standalone probes and JSON/logs.
The first probe covers history/export, DIF persistence and default-vs-selected
bias calculations; the second exercises the actual network renderer; the third
checks the actual bias UI metric and captured local CSV against shared report
text. Its first two harness assertions incorrectly expected raw `Flag` and
`MultiplicityFamily` columns in a deliberately abbreviated displayed table;
inspection was corrected to the CSV bytes and displayed metric. Product code was
not changed to satisfy those assertions. Figure-download rendering alone was
stubbed in that last probe. The MFRM calibration entrypoint was forbidden; actual DIF screening and descriptive
network calculations ran.

Application, locale and preceding qualification/regression sources remain
identical to the prior fix. No full test suite or successful numerical studies
were repeated. This review does not establish browser-native file-save completion,
novice comprehension, public deployment behavior, every module's reachability or
all possible saved-result transitions. Existing working-tree changes are retained.


<a id="api-connectivity-fix-20261027"></a>
### Connect retained results and screening settings — 2026-10-27

**Outcome: the four reproduced handoff defects are corrected locally.** This is
an UX1 repair of existing analysis/report paths. It adds no estimation method,
mandatory novice screen, separate result-management framework or inference
qualification. The working directory remains unpublished.

| Broken connection | Change and observed result |
| --- | --- |
| Optional analyses → current fit | `_model_choice_cache_key` now uses `build_result_analysis_identity`, covering retained responses, configuration and fitted coordinates. Model comparison, network and assignment-sensitivity consumers share the corrected key. In the actual network renderer, switching the retained 30-person design (40 nodes/332 edges) to its 10-person subset hides the first result until explicit execution, then yields 20 nodes/132 edges, equal to direct computation. Returning to the first design reuses its result: two network calculations, no third calculation. |
| Session prediction/refit results → history/export | New prediction/simulation metadata begin with the same full identity. `_current_prediction_bundle` rejects foreign, legacy and missing ownership at the shared export boundary. Real three-row predictions and a synthetic refit-summary marker are included for their source analysis and excluded after an actual history restore to a different retained-data fixture. Settings JSON and generated runner defaults read retained, matching computation settings instead of transient widgets. No expensive refit stress test was run for this ownership check. |
| Classical DIF → continued viewing/save | The session retains group data, controls and the API's complete returned bundle, including settings, audit and caveats. Language changes, panel navigation and ordinary reruns preserve the table and CSV with one calculation. Changing alpha, Person/group associations or retained responses hides the old result until Run. Download clicks use Streamlit's non-rerunning mode. DIF remains an external-group observed-score screen, distinct from fitted-facet bias. |
| Bias controls → shared reporting | The callback stores thresholds in each pair's bundle. Shared screening, boundary/inference audits, First Read, general CSV, basic/detailed PDF/Word/HTML, settings JSON and the generated app-engine runner consume the retained settings. Explicit API overrides and older bundles' defaults still work. The real Person × Rater screen at minimum cell count 100 now reports 120/120 flags in local CSV, general CSV and report notes. The old default report had 24/120. Raw cell estimates are unchanged and public inference remains withheld. |

Settings are described in user-facing terms, including minimum ratings per cell
and the practical logit cutoff. Identical settings across pairs are grouped in
the common explanation instead of repeating every pair. A separate two-pair
fixture verifies differing settings are preserved across summaries, audits and
all six basic/detailed document-format combinations. Generated Python code was
compiled and its actual bias-export block executed against retained cells; that
check did not refit or validate the entire standalone runner.

Verification and evidence:

- `tests/test_api_connectivity.py` adds 18 bounded checks for ownership, DIF
  navigation/invalidation, per-pair controls, common exports, six document
  combinations and the generated runner handoff. Existing cache-version coverage
  now includes reports from `report-first-read-qualification-v1`.
- The connected regression selection passed **192 tests in 23.32 s**. It covers
  the new tests, information-criteria cache policy, conditional-bias metadata and
  mathematics, publication holds, decision stability, First Read, report export,
  optional panels, group assignment and intentional public-route exclusions.
- After the final explanatory wording change, the affected API-connectivity and
  First Read selection passed **48 tests in 5.83 s**. These overlap the 192 and
  are not an additional independent 48 cases. The final source hashes and both
  stage logs are recorded; this is staged verification, not a fresh full-suite run.
- `verify_saved_paths.py` and `saved_path_checks.json` capture the real network,
  history and bias UI paths using the previous artificial Writing output. Its
  numerical findings remain valid after the final wording-only change; its
  retained note text records the earlier wording. MFRM fitting was blocked in
  those probes. Existing regression tests include their own limited fit fixtures;
  they are not new numerical-validation studies.
- `check_reports.py` / `report_checks.json` contain final PDF/Word/HTML checks from
  that same retained result. All four existing report tables are unchanged; two
  already verified figures are reused. Final PDF text retains selected thresholds,
  120/120 counts, the PCA caution and inference hold. All eight rendered pages were
  inspected for clipping, overlap and section transitions. The saved result was
  unchanged by document generation. These files are local QA artifacts.

Development failures are retained separately from successful evidence. The first
new test run had 15 passes and one incomplete DIF fixture (missing score bounds),
then the corrected targeted UI tests passed. The saved-path script initially
needed the repository on `PYTHONPATH`; its next comparison mistakenly compared
rounded UI values with full-precision values. The corrected probe checks raw
results exactly and display values at their documented precision. The first
report assertion expected a second pair in a one-pair saved fixture; multiple
pairs are instead covered by the explicit two-pair tests. A wording edit briefly
introduced a test keyword typo and an unformatted fallback string; compilation
and output assertions caught both, and the final 48 checks and regenerated
reports passed. An editing-time indentation error and accidentally removed
model-choice label constant were corrected before the connected run; final AST
comparison confirms all pre-existing top-level definitions remain present.
Completed numerical calculations were reused through these harness/copy fixes.

Limits and disposition:

- This validates the reproduced local paths, not every API's semantic reachability,
  concurrent sessions, browser-native file saving or novice comprehension.
- DIF persistence is within the current Streamlit session. History snapshots are
  not expanded into portable archives of every optional analysis. Unknown/foreign
  optional evidence is excluded; this does not guarantee every non-default input
  or uploaded-file restoration transition in the prediction panel.
- The bias counts remain technical screening counts. No bias/no-bias conclusion,
  interval qualification or cross-engine numerical agreement is established here.
- D1 design planning, rater feedback, intentional research/legacy exclusions,
  Tkinter packaging and session-independent project files retain their existing
  disposition. No ConQuest/mfrmr rerun, stopped-study restart, push or deployment
  occurred. Previously valid numerical evidence is not reclassified by these UI
  tests. The next product priority remains the existing UX1/J0 acceptance scope.

Evidence directory: `validation/generated/api_connectivity_fix_20261027/`.


<a id="optional-analysis-lifecycle-20261027"></a>
### Optional-analysis lifecycle and generated API arguments — 2026-10-27

**Outcome: four further handoff defects were reproduced and corrected locally.**
This extends the preceding bounded review into its explicitly unverified
non-default prediction inputs and upload transitions. It also follows generated
API arguments through their actual consumer. The UX1/J0 scope, statistical
qualification and publication boundaries remain unchanged.

| Reproduced defect | Correction and evidence |
| --- | --- |
| Prediction/simulation controls disappear when their panel is hidden | Streamlit deletes widget keys. Returning with default controls invalidated a correctly owned prediction, fixed-parameter simulation and refit marker even though the user had not changed the inputs. The panel now retains controls and uploads separately from widget lifetime, scoped to the complete analysis identity. Same-analysis navigation reuses results; changed source/inputs invalidate them. |
| Explicit DIF file removal reuses the old group table | The preceding persistence fix fell back to the retained parsed group table even after the user removed the file. The new upload helper distinguishes cleanup from an explicit change callback. Removing/replacing a file invalidates old DIF input/result records. Without an upload, current fitted person data are used when available and identified by the existing caption; otherwise the screen waits for groups. This also avoids freezing a stale copy of fitted group data. |
| `bias_pairs` is accepted by the runner generator but ignored in execution | The generated `BIAS_PAIRS` constant had no runtime consumer. Explicit pairs now take precedence; `[]` means no bias execution and `None` uses the existing saved-mode fallback. Tests execute generated pair selection and the bias-export block, including the prior per-pair screening-setting contract. |
| Script input mappings/settings come from another analysis | The current-engine Python, portable Python and retained R generators read `facets_mode_output` without checking its owner. A foreign or unowned snapshot changed the requested script's Person/Rater/Score columns and maxit. All three now share `_script_input_context`, which accepts only an identity-matched result. Current matching mappings remain available. R is not reintroduced into public routes. |

A single small upload helper is used by new-design rows, supplemental person data
and DIF. It resets the file position before reading, keeps the retained filename
visible after navigation, offers a localized removal button and uses an
analysis-specific widget key. The file object is retained within session memory;
this adds no disk storage, portable project archive or background calculation.
Fitted group data are read from the current result rather than an obsolete copy.

Verification:

- Before modification, the first five checks produced **3 failures / 2 passes**:
  navigation lost the pasted design, explicit bias pairs were ignored, and an
  explicit empty list still ran the configured pairs. A separate DIF check
  reproduced the removed-file reuse. A later input-mapping check produced
  **2 failures / 1 pass**, distinguishing foreign/unowned and matching snapshots.
  These are expected reproductions, retained as `before_*` logs.
- The final connected selection passed **82 tests in 14.13 s**, including 17 new
  lifecycle/argument checks, the preceding API-connectivity checks, group-upload
  CSV/XLSX handling, optional tools, cache policy and public-route exclusions.
  Intermediate 24-, 72- and 7-test successes overlap this result and are not
  additional independent evidence. No full numerical suite was rerun.
- The generated Python scripts were compiled and their constants inspected.
  Tests execute pair selection plus bias export against recorded cell fixtures.
  R checks inspect generated mappings; they do not execute R or establish
  cross-engine equivalence. No full generated-runner refit was needed for these
  input-mapping and argument-consumption changes.
- `verify_saved_prediction.py` runs the production prediction panel with the
  previously saved artificial Writing fit (960 ratings, 30 persons). Actual
  APIs predict three supplied rows and make two fixed-parameter replicates.
  A synthetic refit-summary marker tests expensive-result retention without
  performing any calibration. After hiding the panel, switching English to
  Japanese and returning, all three bundles retain settings and values, and
  their CSV bytes are identical. Calculation counts stay at one each. Removing
  the retained design file invalidates only that prediction, without another
  calculation. The 2-replicate check concerns data flow, not simulation precision.
- File widgets are supplied with BytesIO objects and their change/cleanup events
  are simulated at the AppTest boundary. Both languages, repeated reads, explicit
  native-clear callbacks, retained-file removal, same-name replacement and
  a changed analysis are covered. This is not native browser upload/save acceptance.
- The AST comparison identifies changes only in the two UI sections and three
  script generators, plus the two small helpers. No pre-existing top-level
  definition is removed. Source hashes, before snapshots, test logs and the
  saved-result probe are in `validation/generated/optional_lifecycle_fix_20261027/`.
  The saved probe emits a nonfatal Streamlit log about restored widget values
  also having defaults; there is no AppTest exception or exposed error screen.

Disposition and limits: immutable history snapshots still represent their
creation-time outputs; no live settings are silently written back into old
snapshots. Optional results are retained for the current session/analysis, not
archived for all past analyses or across restarts. Other report and numerical
checks from the preceding fix remain applicable; no document layout or inference
rule changed here. Existing fitted-simulation CSV downloads remain a separate
save path, not a newly required novice report. Browser-native file operations,
novice comprehension and public deployment acceptance remain unverified. D1,
Tkinter packaging, research/legacy exposure and stopped numerical studies retain
their existing disposition; no external-engine run, push or deployment occurred.

<a id="desktop-foundation-20261027"></a>
## Shared basic analysis, durable project and Tkinter prototype — 2026-10-27

The [desktop prototype contract](desktop_prototype.md) records the implementation,
OS targets and acceptance boundaries. This turn moved the basic estimator,
diagnostics and identity dependency closure into `mfrm_app/engine.py`; the Web
entrypoint imports the same functions. Of the moved functions/classes, 147 retain
identical ASTs. The semantic PCA-label class no longer reads Streamlit state;
Web rendering retains English/Japanese translation, including active-locale
fallback. The publication measure-table implementation is shared separately.

The RSM/JMLE workflow now reaches typed JSON projects, reopened views and
self-contained HTML. Projects bind the original input, normalized input,
settings, numerical results, tables and saved qualification. Unsupported schema,
changed content, incomplete column mappings and unsupported objects are rejected.
The initial inverse-Hessian optimizer operator is deliberately omitted; this is
not an optimizer checkpoint. The desktop table uses the Web publication rule
that blanks withheld SE/interval cells. Saving a project/report no longer replaces
the persistent qualification explanation with a transient success message.

Evidence: `validation/generated/desktop_foundation_20261027/`. Tests were selected
for the changed dependency boundary and report/archive paths; counts overlap:

- `extraction_tests.log`: 170 passed, before the publication-table extraction.
- `numerical_regression.log`: 42 passed and one failed test targeting the former
  global namespace. Its v1-identity injection was moved to the defining engine
  namespace; `final_targeted_tests.log` then passed 53 tests (18 deselected).
  Allocation-guard monkeypatches were also moved, preserving their intended guard.
- `final_report_project_tests.log`: 118 passed after sharing the publication table.
- `final_project_tests.log`: final archive/API contracts, 27 passed. This includes
  fresh-process restoration with estimator and diagnostics calls forbidden,
  identical HTML, typed IDs/categories/nonfinite arrays, held intervals,
  atomic-save failure, unknown versions and missing column mappings.
- `final_locale_tests.log`: 31 passed for translated PCA fragments and first-read
  qualification after the Web renderer adjustment.
- `writing_roundtrip.json`: the saved 960-rating Writing result retained
  `ana_0d4fb0699d783f6031c3416f`; no re-estimation.
- `source_smoke/` and `bundled_final/`: real spawned-worker cancellation, input
  failure, restart, project and HTML creation. `bundled_reopened/` and
  `accepted_reopened/` reopen in separate executable processes; HTML bytes match.
- `gui/checks.json`: native Tk widget invocation reaches three tables, a plotted
  figure, project save, restore without fit, HTML save and persistent qualification.
  File-dialog selections were supplied programmatically. Window screenshots were
  inspected; this does not establish human usability or native file-dialog acceptance.

The first build failed at the sandbox-protected PyInstaller cache path. The cache
was relocated to the build output. A subsequent frozen run completed computation
and project save but failed at SVG export because `backend_svg` was absent; the
build now explicitly includes SVG and Tk backends. These were packaging failures,
not failed numerical experiments. Final source hashes and the retained `.app`
are recorded beside the logs; the final GUI/status and archive validation changes
were checked locally, and the accepted bundle reopened existing results without
re-running successful numerical work. The macOS copy passed `codesign --verify
--deep --strict` with its development ad-hoc signature, not Developer ID signing
or notarization.

Both Windows 11 x64 and macOS Apple Silicon are targets. Only macOS 27.0.1 arm64
was exercised here; Python is installed on that development machine. Windows,
minimum macOS version, a Python-free clean machine, real file-dialog operations,
read-only-input acceptance, signing/notarization and novice understanding remain
open. No new ConQuest/mfrmr comparison, stopped study, public push or deployment
was performed. The desktop prototype does not replace the Web UX1/J0 work.

<a id="public-web-review-20261027"></a>
## 公開Webの実画面と用語・利用条件の点検 — 2026-10-27

**公開版にも修正が必要。最優先はPathwayの誤解を招く判定説明と、図と説明の不一致である。** 商用利用と利用者責任の明示は、アプリ側の誤表示を残す理由にはしない。続いてUX1の既存範囲で、入力例→結果→保存の主経路を短くする。新モデル、停止中研究、mfrmrとの新しい比較実験を今回の必須作業に追加しない。

### 確認対象と事実

[公開Web](https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app/)をChromeで開き、休止画面から起動、日本語切替、ライセンスと版情報、作文サンプル960評定の入力確認・実行・ガイド完了、図の切替を操作した。公開画面の版は **0.2.18-beta / b483d97**。組込みサンプルの実行時間表示は **4.7秒**、RSM/JMLE、30反復で収束と表示された。これは1回のアプリ内計時で、起動待ちを含む利用者の総所要時間や速度保証ではない。私的データのアップロードは行っていない。

実画面のアクセシビリティ情報と画像、公開コミットのソース、今回の局所検証は`validation/generated/public_web_review_20261027/`に保存した。公開ソースは版表示に一致する`b483d97`をGit履歴から読み出した。統合版の未コミット変更が公開済みとは扱わない。

### 問題と採用する対応

| 優先度 | 確認できた問題・影響 | 対応と現在地 |
| --- | --- | --- |
| 最優先 | 公開Pathwayの実画面に「ZSTDの0が完全適合」「±2内はよく適合、外は有意な不適合」という説明がある。近似値・評定数・多重比較等の条件を読み飛ばして、色で合否を決めさせる。 | 開発版で実装済みのMnSq初期表示、ZSTDは任意の参照線、合否帯の撤去を公開候補に含める。公開未反映。 |
| 最優先 | 同じ公開Pathwayで、図は横軸Infit ZSTD・縦軸Measureなのに、別の読み方ガイド見出しは横軸Measure・縦軸Infit MnSqと記す。 | 開発版の修正を、表示・読み方・保存図の組として公開確認する。グラフが描けたことだけで合格にしない。 |
| 高 | サンプルの「1行が表すものを確認」には件数と列の役割表があるが、実際の評定行がない。「ステップ2/5」と進行状況の文章も重複する。906×768の取得画像では次へが画面外。 | 先頭数行を主に置き、列の役割はその近くで確認する。進行表示は一つにし、「次へ」を短くする。CSV説明自体は既存アップロードヘルプにもあるため「説明なし」とは判定しない。自分のCSV経路の実操作は未確認。 |
| 高 | サンプルの結果段階で、一般的な正誤二択が結果より先に出る。答えは当該結果を読まなくても選べる。ガイドは終了できるため分析全体の強制ゲートではないが、推奨経路の中心を占めている。 | 実際の結果と該当する保留・次の操作を先に置く。理解確認は任意の「学ぶ」に置き、残すなら実際の図表を読む課題にする。注意事項への同意を分析の妥当性確認と扱わない。未実装。 |
| 高 | 日本語でも最初の確認事項に英語のPCA文、P1/P2、技術的な相関値が並ぶ。ベータ版の説明にもReady／Remediation hold等の開発語と実装ライブラリ名がある。 | 利用者に必要な「何が確認を要するか」「どの結論を保留するか」「どこを開くか」を日本語で示し、数値の根拠・実装情報は詳細へ。実際の保留や根拠を削除しない。公開と開発側の差を確認して修正範囲を決める。 |
| 中 | 「図」→「図表」→「可視化パネル」と選択が重なり、その間に機能説明、再推定しないという説明、Evidence basisとlabel戦略等がある。図保存は別セクションの深い位置へ案内される。 | 図名の一つの選択欄、近接した「読み方」「保存」を基本とする。描画と保存で指標・対象・分析ID・保留を揃える。Plotly自体のPNGボタンは存在するため「保存できない」とは判定しない。公開の全保存形式の受入は未完了。 |
| 中 | 休止後の最初の訪問に英語の起動待ちが挟まる。 | Community Cloudの運用条件としてWeb／配布版の選択説明へ反映する。今回の起動成功を常時即時利用の保証にしない。 |

MnSqの期待値は1、ZSTDは0で、ZSTDは標本の大きさに敏感である。[Winstepsの一次資料](https://www.winsteps.com/winman/misfitdiagnosis.htm)もmean-squareを先に点検することを勧める。±2は近似検定の参照にはなり得るが、その内側は適合の証明ではなく、外側を自動的に実質的な不適合とする扱いにも条件が必要である。MnSqの特定の帯も普遍的な合格基準にはしない。今回、数値計算の正しさを新たに認定したわけではない。

### 用語と利用条件の修正・検証

- 配布版の操作・ダイアログ・ヘルプ・説明を「プロジェクト」に統一した。`.mfrm`のschema、内部API、内容は変更しない。
- 既存のMIT本文は変更しない。既存`LICENSE_NOTICE.md`とREADMEはすでに商用利用と利用者責任を明記していた。日本語の案内を補い、Webの日英「アプリ情報」と配布版ヘルプに、研究・教育・業務・商用で利用でき、データの取扱い、適用判断、結果の解釈・利用は利用者の責任で、無保証であることを揃えた。一般画面の常設警告や同意チェックは増やさない。
- 配布ビルドで漏れていたアプリ自身の`LICENSE`と`LICENSE_NOTICE.md`を同梱する。ビルド記録にも両ファイル、同梱locale、ビルドスクリプトと固定依存設定を含める。第三者ライセンス全件の配布条件を監査したことにはしない。
- 日英のWebアプリ情報をAppTestで確認し、Tk実部品で開く／保存する名称、ヘルプ、既存プロジェクトの往復を確認した。保存は元と全バイト一致、推定・診断再計算は0。初回は古い`.venv`が存在しないPython 3.12を参照して実行できず、利用可能な3.14へ切替。Tk確認はサンドボックス内で終了コード134となり、通常のネイティブ環境で成功した。これらを統計テストの失敗とは扱わない。
- macOSパッケージを更新し、同梱利用条件とソースhashを照合した。独立した実行ファイルで既存プロジェクトを再読込し、推定0・HTML全バイト一致、コピー先のアドホック署名を確認した。今回の文言変更のために成功済みの数値実験は再実行していない。

画面取得ツールには途中でScreenCaptureKitの取得失敗と更新前画像の返却があり、アプリ障害とは区別した。更新前の画像は実行前画面として保存名を修正し、結果確認には実際のアクセシビリティ情報を使った。公開CSVアップロード、各形式のダウンロード完了と保存ファイル、狭い画面全体、支援技術による操作、初心者の理解、Windows・クリーン端末は今回の合格範囲に含まない。公開更新・push・一般配布は実施していない。

<a id="literature-based-ux-20261027"></a>
## 先行研究と図の解釈・学習導線の照合 — 2026-10-27

**文献を付けるだけでは不十分である。文献が調べた条件と、本アプリで採用する判断を分けた。** 今回は前節のPathway・入力例・学習問題の問題に対応し、共通の図の説明も点検した。数値エンジンや判定資格を変更せず、未検証の閾値を新たに導入しない。主画面に参考文献一覧を常設せず、読み方ガイドと本記録で根拠を追えるようにする。

### 根拠と適用範囲

| 論点 | 確認した文献・証拠の強さ | 本アプリで採用すること／採用しないこと |
| --- | --- | --- |
| MnSqとZSTDの基準 | **Müller (2020)**。公開本文の方法・結果・考察を確認。二値Rasch項目を対象に、項目CMLE・Person WLEを用いた残差と条件付き残差を比較。Infit/Outfitの実験は150〜10,000人、10/15/20項目、各条件1,000反復。基準分布と誤検出率は推定方法・項目数等に依存した。 | MnSqはずれの大きさを読む入口にし、普遍的な検定や合格帯とは呼ばない。論文中の標本数別の境界値を、多相・多値のRSM/JMLEへ移植しない。今回のアプリの誤検出率を測定した証拠でもない。 |
| Wright→Pathwayの学習順 | **Boone (2008)**。全文確認。14人×12項目の教材を使った大学院授業の実践報告で、位置を読む図から適合度を加える図へ段階的に進む。無作為比較による学習効果の検証ではない。 | 既存のWright初期表示を維持する。縦軸の位置と横軸の適合度を明示する。元のPathwayで扱う精度表現と本アプリの一定サイズの点を混同しない。図の形が似ているだけで推定精度やMFRMの妥当性を保証しない。 |
| カテゴリー曲線の解釈 | **Linacre (2002)**。出版社掲載情報とPubMedの著者抄録・書誌情報を確認。カテゴリー度数、順序、測定との対応などを組み合わせる指針と2データセットへの適用を述べる。全文の全8指針を今回再検証したわけではない。 | 度数・閾値・採点基準の意味と併読する。ピークが見えないだけでカテゴリーの誤用や統合を決めない。既存の「Linacre 2004」という別題の論文への引用を、該当する2002年のカテゴリー研究へ修正した。 |
| 説明と対象の配置 | **Chandler & Sweller (1992)**。出版社の抄録・書誌を確認。NCプログラミングの図と説明、研究報告内の相互参照する文章を統合した2実験を報告する。全文・個別の効果量は今回未確認。 | 列の意味を実データの列に添え、図の近くから読み方を開く。本アプリの配置を直接比較した研究ではないので、改善効果は初学者の操作・理解で確かめる。説明を全削除する根拠にはしない。 |
| 単純化の条件 | **Chandler & Sweller (1996)**。出版社の抄録・書誌を確認。CAD/CAM学習では要素間の相互依存が高い条件で統合した教材の優位があり、低い条件では形式間の差がなかった。 | 「少ない説明・少ないクリックなら常に良い」と一般化しない。推論保留、列の役割、得点基準など、判断に必要な情報を残す。 |

二択問題を任意で開く形にする判断は、今回の画面で実データを読まなくても回答でき、結果の前に大きく表示されていたという観察に基づく。上記の認知負荷研究が「小テストは無益」と実証したわけではない。問題・フィードバック・学習状態は保持し、閉じたままを学習完了とは記録しない。

### 実装への反映

1. サンプル確認は、実際の先頭5評定と、同じ列に添えた役割を主にした。大きな件数カード、別の役割表、一般注意文の重複を外し、短い件数と単一の進行表示、「次へ」にした。元の960評定や分析設定は変えない。
2. 結果に先行していた学習問題・学習メモを、閉じた「学ぶ：結果の読み方」に収めた。実際の警告はその外に残る。既存の問題への誤答・正答・完了の処理、分析ID、退出・復帰は維持する。
3. Wright／Pathway／カテゴリー曲線の詳細説明をWebの日英・配布版・配布版HTMLで共有した。引用のある英語文を同じWeb画面で重ねて表示せず、日本語選択時は日本語ガイドを出す。Müller・Booneを既存の文献台帳と引用解決へ接続した。UI設計の研究を、利用者の分析手法の引用として報告書へ自動追加はしない。
4. Wrightの位置の重なりだけから天井・床効果を判断する説明、カテゴリーのピーク／閾値から誤用を断定する説明を修正した。ファセットの同定制約とカテゴリー境界を考慮する必要があるという実装上の条件を、授業実践の効果の主張と区別した。隣接するECDF説明の「密集しているから精度が高い」という誤った含意も除いた。
5. Pathwayの指標選択・参照線の説明を日英化し、点の詳細で利用可能なMnSq・近似ZSTD・評定数を確認できるようにした。一定の点の大きさがSEを表していると誤解されない説明を加えた。初期指標、参照線、座標、推定方法、保留規則は維持する。

### 文献

- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust them? *Journal of Statistical Distributions and Applications, 7*, Article 5. https://doi.org/10.1186/s40488-020-00108-7
- Boone, W. (2008). Teaching students about Rasch maps. *Rasch Measurement Transactions, 22*(2), 1163–1164. https://www.rasch.org/rmt/rmt222a.htm
- Linacre, J. M. (2002). Optimizing rating scale category effectiveness. *Journal of Applied Measurement, 3*(1), 85–106. https://pubmed.ncbi.nlm.nih.gov/11997586/
- Chandler, P., & Sweller, J. (1992). The split-attention effect as a factor in the design of instruction. *British Journal of Educational Psychology, 62*(2), 233–246. https://doi.org/10.1111/j.2044-8279.1992.tb01017.x
- Chandler, P., & Sweller, J. (1996). Cognitive load while learning to use a computer program. *Applied Cognitive Psychology, 10*(2), 151–170. https://doi.org/10.1002/(SICI)1099-0720(199604)10:2%3C151::AID-ACP380%3E3.0.CO;2-U

補助資料の[Winstepsの適合度解説](https://www.winsteps.com/winman/misfitdiagnosis.htm)は、指標を読む実務的な説明として扱う。論文の検証条件や本アプリの外部数値照合を代替しない。Smith et al. (1998)の抄録も調べたが、本文未確認の臨界値式は採用しなかった。

### 検証と未完了の範囲

証拠は`validation/generated/literature_ux_20261027/`。変更した関数は文献対応1、図の表示2、読み方1、初回案内3の計7。共通計算、保存形式、資格判定は前回パッケージと同一で、native図の関数ASTも変更していない。

- 初回案内の言語切替・退出／復帰・実際のサンプル分析から学習完了までの3テストが通過。閉じた問題欄を学習完了にしないこと、誤答・正答・分析IDの保持、実際の先頭5評定の一致を確認した。
- 図・文献・Help・読解案内の55テストが通過。続いて日英3図の共有説明・出典解決・閉じた読み方の6テストが通過。点の詳細ラベルの最終変更に対する実描画の3テストも通過（工程間に重複あり）。
- macOSパッケージを更新し、同梱ソースhashとアドホック署名を確認。既存の48評定プロジェクトを独立した実行ファイルで再読込し、推定0。HTMLには出典が追加され、3つのSVGは前回パッケージから全バイト不変だった。文言変更によりHTML全体が変わることは意図した差分である。
- ローカルの起動確認では、初回のサーバーbindがサンドボックス制約で失敗し、localhost限定の通常実行で起動した。Chromeで初期画面・日本語切替は確認したが、`noWindowsAvailable`、ScreenCaptureKitの取得失敗、古い要素番号が繰り返され、サンプル画面の実操作を最後まで確認できていない。ブラウザー上の応答と操作基盤の影響の切分けは未完了であり、AppTestの合格を実ブラウザー受入へ読み替えない。確認用サーバーは終了した。

公開版は前節のb483d97のままで、push・公開更新は行っていない。初心者による理解・操作、Windows・クリーン端末、全出力のネイティブ保存、新たなConQuest／mfrmr照合は未確認。図の入口と保存の統合、初回結果に残る英語の説明等は前節の優先課題として残る。今回の引用追加で、文献台帳全体や全診断の説明まで監査済みとはしない。

<a id="web-ui-rewrite-review-20261027"></a>
## Web UIのリライト範囲の点検 — 2026-10-27

**画面構成と画面状態の管理を段階的に書き直すことを推奨する。** 説明の短縮だけでは、同じ図への複数の入口、閲覧と保存の分断、表示設定を分析変更として扱う構造が残る。既存UX1の作文例で入力から報告保存まで進む目的に沿い、Streamlit内で主経路を組み直す。全計算の再実装、別フロントエンドへの移行、Webを配布版の現状のRSM/JMLE範囲へ縮小することは、この点検から必要とは判断しない。

今回は**評価と計画の精緻化**であり、新UIの実装・受入完了ではない。公開b483d97の実画面は前節の保存証拠を再利用し、現在の統合版ソースと照合した。開発版の改善済みの文言・任意の学習欄を、公開版にも反映済みとは扱わない。

### 確認した構造上の問題

| 優先 | 事実と確認箇所 | 利用者への影響・変更する範囲 |
| --- | --- | --- |
| 高 | `run_facets_mode`は小数桁、図の描画、保存図の作成、可視化設定の変更も`stale_reasons`へ追加する。通常の結果は非表示となり再分析ボタンを出す。復元済み結果には別の扱いがある。 | 「桁数や図の保存を変えたい」が分析のやり直しとして案内される。結果の取り違え防止は必要だが、入力・推定設定、追加診断、表示・出力設定を分ける。古い結果を現在の設定の結果として見せる修正はしない。これはコード経路の確認で、全ケースの実ブラウザー再現ではない。 |
| 高 | `resolve_analysis_depth_settings`の標準では図を表示できるが、保存図セットは無効。Full publicationは保存図だけでなく全ファセット対のバイアス検査やEB等も有効にする。`_render_downloads`は保存画面で図を生成するのに、この事前フラグで基本図への入力を空にする。 | 保存したくなった段階で、分析前の選択へ戻す構造。既存結果から作れる図はその場で生成・保存する。新たな計算が必要な診断は別の明示的操作とし、未計算をゼロ・問題なしとして扱わない。Full publicationという名前を報告の妥当性と混同させない。 |
| 高 | `_fit_summary_callout`は0.50–1.50外の要素比率で、0%をsuccess、10%以下をinfo、25%以下をwarning、それを超えるとerrorにする。日英の文言には既に「目安」「証明ではない」という限定がある。 | 文言の改善後も表示の強さが合否を印象づけるおそれがある。10%／25%をこの用途の警戒水準とする根拠は、今回のコードと参照文献記録では確認できない。目安外の件数と該当行を中立的に示し、推論保留・入力不備・処理失敗と区別する。分類の丸め規則や実際の保留を削除する意味ではない。 |
| 中 | `_render_guided_figures_section`から`show_visuals_section`へ入って図を選ぶ。診断側にもWright／図表がある。図セット・報告書は別の関数で図を組み立て、`_render_downloads`／`_publication_figure_payloads`には現在のPathway指標選択を渡す経路がない。 | 「図 → 図表 → 可視化パネル」と「レポートと出力 → 表・図を保存 → 図」を覚えさせる。図名の一つの選択欄と、その図の読み方・保存を近接させる。Plotly標準の画像保存は存在するので「保存不能」とはしない。Infit×Outfitの散布図はMeasure×InfitのPathwayとは別物であり、代替として扱わない。 |
| 中 | 「表示モード：コンパクト／全パネル表示」「設定の範囲：標準／詳細」「分析する内容：結果のみ／標準／図・出力を含む」が別に存在する。コンパクトは7セクション、全表示は16パネルの別分岐を持つ。 | 初心者に、表示量・操作項目・計算範囲の違いを先に学ばせる。表示は一つの主経路に揃え、専門的な設定と診断は該当箇所から開く。詳細機能を削除したり、見つけられない場所へ移したりしない。選択肢の数だけで認知負荷や誤操作率が実証されたとはしない。 |
| 中 | `_render_guided_essential_tabs`は選択した全セクションの手前で推奨確認先のprimaryボタンを出す。`main`はHome／Help／G・D画面でも入力部品を維持し、表示しない部品の値を失わないよう複数のキーを保持する。 | 次の確認先への案内は強みだが、保存等の主操作と競合し得る。推奨先は概要、実際の保留は関連する結果へ配置する。ページを単に分割すると設定・復帰を壊すため、部品の寿命と分析状態を分ける。現在のコードに保存結果の消失を確認したという指摘ではない。 |

適合度の色分けは、100要素のうち0／10／25／26要素だけInfitを1.51、残りとOutfitを1.0にした4ケースをAppTestで表示し、上表のコンポーネント種別と日本語文を確認した。数値は表示検査用の合成値で、分析や適合度の妥当性を検証したものではない。証拠は`validation/generated/web_ui_rewrite_review_20261027/fit_callout_probe.json`と同ディレクトリのソースhash・確認関数一覧。推定は実行していない。

### 残す強みと、書き直す境界

Pythonアプリの価値は、MFRMの計算に加え、**結果から確認すべき根拠をたどり、言える範囲を保って報告へ進めること**にある。目的別の案内、割付・カテゴリー・感度の点検、分析ID、資格判定、履歴、日英Help、根拠と報告文の対応を維持する。mfrmrの関数一覧の再現や、図だけを大きく見せる改変を完成条件にしない。

画面の外枠・ナビゲーション・設定編集・表示と保存の接続が書換え対象。既存の推定、結果の同一性・派生結果の所有元確認、出力資格、報告内容の生成は再利用する。`mfrm_app/basic_analysis.py`は現状の配布版の基本分析範囲であり、Webの全モデル・推定法を置き換える共通入口とみなさない。計算基盤の改善が必要な場合も、UI書換えを理由に統計的な契約を変更しない。

### 提案する主経路

| 画面 | 主に置く内容 | 必要時に開く内容 |
| --- | --- | --- |
| データ | CSVを開く／サンプル、実際の先頭数行、列に添えた役割、次へ | CSV見本、形式の説明、貼付け、ほかの既存データ源。読込エラーは該当箇所で具体的に説明する |
| 設定 | 列の役割、採点基準としての得点範囲、モデル、分析する | 推定法・アンカー・重み・追加診断等の詳細。観測最小・最大を、採点基準の可能な得点範囲と自動的に同一視しない |
| 結果 | 「概要／表／図／診断／報告」。選んだ内容を主に表示し、表・図の近くから保存する | 詳しい根拠、感度・追加分析、報告資料一式、履歴・比較。学習課題と理論はHelpまたは対象の「読み方」から開く |

「データ → 設定 → 結果」は入口の作業順であり、結果閲覧を一方向のウィザードに固定しない。結果内を自由に行き来でき、設定を編集しても直前の分析を識別して保持する。一般的な注意文や手順番号入りボタンで空間を埋めず、関連する入力と操作の配置で次を示す。数式・原著への案内は近接した読み方から開き、確認を要する値・保留理由は閉じたHelpだけに隠さない。

「報告」は研究・教育・業務に共通する入口にし、論文用の下書き・文献・再現資料はそこで選べるようにする。研究用の報告機能を消さず、商用利用者にも論文作成を必須経路として見せない。G・D研究は別の目的の入口を維持し、MFRM結果のタブへ混在させない。

### 研究と技術資料から言える範囲

- [Chandler & Sweller (1992)](https://doi.org/10.1111/j.2044-8279.1992.tb01017.x)は、離れた説明と図・文章を統合する教材を比較した研究である。図と読み方を近くに置く判断を支えるが、「3画面」「5つの結果メニュー」が最適だと実証した研究ではない。抄録・書誌の確認範囲と、1996年研究の要素相互依存による条件は前節のとおり。今回の具体的な画面案は、実装と観察からの設計上の推論である。
- [Boone (2008)](https://www.rasch.org/rmt/rmt222a.htm)の授業実践を参照し、図の初期表示はWrightで位置を読み、Pathwayで適合度を併読する順を維持する。経験者の直接選択は妨げない。「インパクト」は読み取れる関係で作り、根拠のない合否色や誇張した精度表現では作らない。
- [Müller (2020)](https://doi.org/10.1186/s40488-020-00108-7)は二値Raschの条件下で適合度の分布・判定を検討している。本アプリの多相・多値の適合度や「目安外25%」の警告水準を検証した研究ではない。読み方の限定と視覚表現を一致させる必要があり、注意書き一文で強い合否表示を相殺しない。
- [Streamlit公式のページ構成](https://docs.streamlit.io/develop/concepts/multipage-apps/page-and-navigation)は関数・ファイルによる画面分割を支援する。現在の課題だけから別フレームワークへの移行を必須化しない。`st.Page`／`st.navigation`は候補であり、3段階を必ず別URLにすることも要件にしない。
- [部品の状態保持に関する公式資料](https://docs.streamlit.io/develop/concepts/architecture/widget-behavior)を確認した。表示しない部品やページ移動で状態を失い得る。ローカルは1.66.0だが`requirements.txt`は`streamlit>=1.54`であり、1.59以降の`persist_state`等を無条件に使って最低対応版でも動くと仮定しない。新しい状態管理基盤を増設する前に、既存の`ux`・`guidance`・`help_navigation`と保存結果の保持を使う。

### 既存UX1内での実装順と受入

1. **誤解につながる表示を先行修正する。** 公開Pathwayの説明不一致と開発版の修正を公開候補へ対応づけ、適合度の件数表示と色の意味も揃える。UI全体の完成を誤表示の修正の前提にしない。
2. **最初の置換単位は「図を選ぶ → 読む → 同じ図を保存する」。** 一つの図選択欄へ集約し、表示した指標・要因・範囲・分析ID・保留を保存へ渡す。標準分析の既存結果だけで生成できる図は、設定へ戻ることや再推定を要求しない。画面を見ての操作と保存内容の一致までを一単位とする。
3. **設定編集と分析結果を分離し、3画面の外枠へ接続する。** 表示変更だけでは分析の同一性を変えず、出力の見た目が変わるときは出力キャッシュを適切に更新する。分析条件を変えた場合は編集中と分析済みを識別し、古い結果を新条件の結果として利用させない。追加診断は必要な処理だけを明示的に実行し、所有元を検査する。
4. **概要・表・診断・報告の重複を順に統合する。** コンパクトと全表示の別構成を、同じ結果の主画面とその詳細へ寄せる。高度な機能の移動先を確認してから旧入口を外す。利用者向けに新旧UI切替を常設し、第四の表示モードを増やす案は採らない。

現行テストには画面の文字列・分岐の存在だけを点検するものもある。例えば`test_run_staleness_tracks_report_scaling_and_visual_outputs`の合格は、その再分析誘導が適切である証拠ではない。書換え時はソース文字列の期待を機械的に維持せず、入力変更を誤って通さないこと、表示変更で推定しないこと、同じ結果の値・資格・保存が一致することを検証する。

既存の通常JMLE・ファセット境界・固定／自由SD MML、日英、Help／Home／学習の退出・復帰、履歴・比較、匿名化出力の証拠を参照し、変更で無効になる経路を選んで再確認する。既存CSV・Word／PDF／HTMLの生成確認を使い、ネイティブ保存の未確認部分を分ける。単純な表示変更のために成功済みの数値研究を一括再実行せず、計算や資格を変えた場合はその影響に必要な独立検証を行う。

実ブラウザーの狭い画面・キーボード・フォーカス・保存と、既定の5人＋2人の形成的レビュー、J0の待ち時間確認は既存手順を維持する。手順数・文字数の削減やAppTestの合格だけを理解の改善と数えない。公開候補の特定、必要チェック、公開後の版と代表操作の確認までが公開完了条件であり、本節の評価・計画更新は公開完了ではない。

<a id="web-figure-workflow-20261027"></a>
## Webの図の選択・読み方・保存を接続 — 2026-10-27

**リライトの第1段階を開発版へ実装した。** 「図」はWright／Pathway／カテゴリー曲線から詳細図まで10種類を一つの選択欄から開く。既存の図を描く関数へ表示・保存の処理を渡し、画面に出したPlotlyオブジェクトをそのままHTML／PNGに使う。標準分析や既存結果の閲覧では、保存図セットの事前設定や再分析を要求しない。計算済みデータがない場合は理由を表示し、空の図の保存ボタンは出さない。

- PathwayはMnSqを初期値とし、近似ZSTDへの切替と要因の絞込みを保存へ引き継ぐ。Helpからの復帰・別の図からの復帰でも指標と要因を保持する。図を直接ズーム・パンしたクライアント上の状態はサーバーのHTML／PNGには戻らないため、保存欄で範囲を説明し、その状態の画像はPlotlyのカメラボタンを案内する。
- HTMLにはPlotlyの実行コードを同梱し、図の閲覧に外部CDNを要求しない。分析ID、現在の推論保留の理由、3図の共有の読み方・出典、概要から渡された診断の確認事項を保存する。保存対象には表示した識別子や点の情報も含むことを保存欄に示す。これは表示中の図を手元へ保存する経路であり、既存の匿名化した共有用出力とは区別する。既存の匿名化出力処理は変更していない。
- 図より前に繰り返し表示していた大きな推奨先案内を、図の画面では外した。概要側の推奨先、上位の重大警告、図の推論保留表示は維持する。診断上の確認事項は図の近くから開け、HTMLにも残る。PNGには分析IDと保留の短い注記、確認事項がHTMLにあることを示す。
- `_fit_summary_callout`の10%／25%によるsuccess／warning／errorの切替を撤去し、目安外件数を中立的な文で示す。MnSqの数値、丸め前の分類、資格判定や実際の警告は変更しない。
- PNG実物で見つかった分析ID注記と軸名の重なりを修正した。さらに、保存済みの対称な48評定例では、ほぼ同じ評定者位置の丸め誤差が縦軸全体へ拡大されていた。Pathwayの縦軸表示幅を最低1 logitとし、微小差の過大な視覚表現を避けた。これは描画上の選択で、差の有意性・実質的重要性の閾値ではない。点の座標は変更していない。

### 検証と証拠の範囲

`validation/generated/web_figure_workflow_20261027/`へ変更前のソース、局所差分、工程別ログ、実ブラウザーのDOM／画像、保存物を保持した。

- 初回の新規8テストでは、日本語の2ケースがAppTest側のラベル解決で失敗した。表示名の辞書を部品生成時に確定する既存方式へ揃えて解消した。続く関連148項目では147が通過し、旧UIの関数呼出し文字列を要求する1項目が失敗した。この項目は、データがない図を保存できたように見せない動作確認へ置換し通過した。
- 図とWeb／配布版の座標・読み方に関する局所9項目、診断の確認事項保存を含む17項目が通過した。実際の通常JMLEの入力からHelp／Home／日英切替／CSV・Word・PDF・HTML出力／設定変更までの既存経路に、図の指標・要因の退出／復帰を追加した1項目も通過。分析後のこれらの操作で推定・再推定呼出しは0だった。
- 最終注記位置の変更に対する9項目が通過。縦軸の最小表示幅を加えた後は、影響する主要図・近似ZSTD・要因絞込み・微小差の5項目を再確認した。工程間の重複を合算して理解度や数値的妥当性の証拠数にしない。
- ローカルWebの日本語画面で作文960評定の標準分析から図へ進み、PathwayのRaterを選び、HTMLを実保存した。ファイル内の対象・MnSq・分析ID・読み方・7件の確認事項・外部script依存がないことを確認した。途中のサーバー再起動は翻訳キャッシュを更新するためで、公開環境の操作ではない。表示時間0.9秒はこのサンプルのアプリ内計時で、一般的な速度保証ではない。
- 最終レイアウトは保存済み48評定プロジェクト`ana_e43282ec8892d30f9da7793e`を推定0で読み込み、AppTestの実描画と保存HTMLの座標・配置・メタデータが一致することを確認し、PNGを独立生成した。ブラウザーの作文例でも、横軸名の下端534.58pxと注記の上端565.08pxが分かれていることをDOMと画像で確認した（注記位置修正後、縦軸幅の最終変更前の記録）。最終PNGも目視し、注記の非重複と微小差を拡大しない縦軸を確認した。

HTML保存イベントの操作ツール応答には約571秒かかった。これは原因を切り分けられていない操作基盤の経過時間であり、J0のアプリ性能値には採用しない。PNGはブラウザー上で作成とダウンロード操作まで進んだが、ファイルの到着を確認できなかった。Chromeのダウンロード管理画面へのアクセスはブラウザー操作ポリシーで拒否されたため、別のブラウザー操作で回避していない。**PNG生成の成功と、ネイティブのPNGダウンロード受入は区別し、後者を未完了とする。**

本実装はコンパクト表示の「図」の置換である。診断／全パネルの旧入口、既存の報告書用図セット、表示設定全般のstale判定、3画面の外枠、詳細図に残る英語は次の既存UX1課題として残す。公開Web更新、初心者の理解、狭い画面・キーボードの全面受入、新たなConQuest／mfrmr比較は実施していない。Streamlitの推定・資格判定、共有エンジン、配布版の実装・ビルドは今回変更していない。

<a id="web-mathematical-scope-20261027"></a>

## 図の数学的な意味と既存検証の適用範囲 — 2026-10-27

**前節のUI・保存検証を、アプリ全体の数学的妥当性の確認と扱わない。** 今回は現行ロードマップ、数理計画、既存の外部照合の保存結果、図と共有エンジンの実装を照合し、小さな代数対照と現行テストを実行した。推定器や図の挙動は変更していない。新しい外部ソフト実行、被覆率研究、停止中の正式研究の再開、公開は行っていない。

### 再利用した証拠と、今回直接確認したこと

- 10月23日の中規模研究の保存記録は、固定SD RSM、192データ・384推定、独立連続積分8点、TAM 4点という範囲である。保存されたTAM比較では全座標差は最大約1.03e−8、主対比差は約1.78e−8だった。24反復／条件の探索的な区間評価であり、全モデル・設計の性能保証ではない。`verification.json`も`scientific_inference_ready=false`を明示している。
- 当時の凍結ソースと現在のソースをASTで照合し、確率、線形予測子、JMLE／MMLの尤度・勾配、パラメータ展開、観測情報、条件付きSE、全体fitの10関数が同一だった。関数の同一性は、全依存関係・実行環境・通常アプリの最適化経路まで同じという証明ではない。
- 10月26日のmfrmr候補とのG研究比較は、完全均衡交差の正常4例・20分散成分で最大差2.66e−15。D研究はその時点の候補で内部関数欠落があり、10月21日の保存版との36計画・72係数の照合（最大差9.99e−16）と区別する。今回、候補の現在状態を再調査したり、mfrmrを再実行したりはしていない。
- 図のGPCM確率を、正規化前の指数を直接列挙する独立式と比較し、最大差1.11e−16。対数確率の有限差分から求めた情報量との最大差は9.69e−12未満だった。対象は同定条件を満たす2水準・3カテゴリーの指定点であり、一般的な数値精度保証ではない。
- 同定、条件付きSE、MMLの有限和による勾配・観測情報の既存69項目を実行し、**66通過・3失敗**。2失敗は共通エンジンへの移設後もテストが`streamlit_app`の依存関数だけを差し替えていたことによる。所有元`mfrm_app.engine`の該当依存を一時的に差し替える局所診断では、微分欠落／有限微分／ゼロ微分の3条件が通過し、共分散を含む分散18と分散0を確認した。残る1失敗は、EAPの説明欄が現行の推論保留文へ置き換わり、旧文言の完全一致期待と合わないものだった。`SE_Detail`の条件付き事後SDという説明は残る。テスト本体は未修正であり、69項目全通過とは報告しない。

証拠は既存の`validation/generated/web_figure_workflow_20261027/`内の`check_mathematical_scope.py`、`mathematical_scope.json`、`mathematical_scope.log`、`tests_mathematical_scope.log`に保存した。

### 数学的な対象・互換性について残った問題

**1. PCM／GPCMの「平均」の確率曲線は、確率の平均ではない。** `build_category_probability_curve_data`は閾値の算術平均と、GPCMでは識別力の算術平均を先に求め、その値から確率を計算する。平均パラメータの仮想モデルとしては有効だが、各水準の予測確率の平均とは一般に一致しない。

反例は閾値(−2, 2)と(−0.5, 0.5)、識別力0.5と2、調整後位置0。各閾値行の和は0、識別力の幾何平均は1で、現行の同定条件内である。中央カテゴリーの各水準の確率はともに0.576117なのに、現在の「Average across Task levels」は0.704611となった。約12.85ポイントの差は丸め誤差ではない。水準別表示を基本とし、平均パラメータの仮想曲線を残す場合はその名称・対象を明示する必要がある。単に「確率の平均」に差し替えると、情報量やカテゴリー診断の対象も変わるため、一緒に再定義する必要がある。

**2. 情報量曲線は実際の評定配置全体の情報量ではない。** 現行の`build_information_curve_data`は、上記の1曲線について`I(η)=a² Var(X|η)`と`1/√I(η)`を計算する。この式は対象の単一応答モデルでは成立するが、複数の課題・評価者による評定を合算していない。実際の評価計画全体のSEや、同時推定したパラメータの全不確実性として読ませてはいけない。Muraki（1993）の情報量の導出と、項目情報とテスト情報の区別を参照した。今回確認したのは[ETS研究報告版の本文pp. 3–5](https://files.eric.ed.gov/fulltext/ED390896.pdf)で、出版論文の同じ式番号まで照合したとは扱わない。

**3. Pythonの`whexact`はFACETSの同名設定と意味が異なる。** [FACETS公式資料](https://www.winsteps.com/facetman64/whexact.htm)では、`q=√(2/df)`として、Yesは`3(MnSq^(1/3)−1)/q+q/3`、Noは末尾の`q/3`を省く。Pythonの`zstd_from_mnsq_facets`ではFalseが前者、Trueが別の線形近似`(MnSq−1)/q`を返す。MnSq=0.975、df=20で公式式のYes／Noは0.025684／−0.079725、PythonのFalse／Trueは0.025684／−0.079057だった。FACETS実機を動かした比較ではなく、公式式との直接照合である。現行の`FitZSTDTransform`は計算の種類を記録するが、引数名からの互換性は成立しない。数値を黙って変更せず、APIの意味、移行方法、表示・保存の説明を整合させる必要がある。

また共有の図ヘルプは一律に「FACETS方式・±9」と説明するが、詳細設定ではengine／bothや上限変更も選べる。診断表が持つ実際の方式・変換・上限に基づく説明へ接続する必要がある。Infit MnSqの参照値1と、近似ZSTDの参照線±2は別物であり、いずれも合格の証明にしない。[Müller（2020）](https://doi.org/10.1186/s40488-020-00108-7)の二値Rasch研究は近似分布の問題を示す根拠であり、本アプリの多相・多値・重み付き設定の閾値保証ではない。

### 評価と範囲

現時点の評価は「限定条件の数式・数値照合はあるが、数学的な意味と表示がすべて整合したとは確認できていない」。Wright／Pathwayの座標と保存の一致だけで上の問題は解消しない。MMLの固定／自由SD、一般のSE・区間の被覆、非単位重みのfitの基準分布、JMLEの付随パラメータと極端反応の限界、モデル・割当仮定は既存の保留と検証範囲を維持する。

次の改修では、既存UX1の「数値・保留・解釈の不整合を先に直す」範囲として、平均曲線と情報量の対象、ZSTDの方式表示・API互換性、移設後のテスト参照を先に整える。大型の性能比較や推論保留解除を、この点検から自動的に開始しない。

<a id="math-ux-repair-conquest-20261027"></a>

## 数学的な対象・利用者向け表示の修正とConQuest実機照合 — 2026-10-27

**前節で残した曲線の意味、ZSTDの互換性、テスト参照を修正した。** UX1の「数値・保留・解釈の不整合を先に直す」範囲であり、MMLの推論保留や公開資格は変更しない。ConQuestを現在の環境で利用できるというユーザーの指摘を受け、既存の固定SD PCMデータ1件で実行環境と数値対応を確認した。新しいモデル・性能研究・一般的な同等性認定は追加しない。

### 画面・保存・共有API

- **カテゴリー曲線：** PCM／GPCMでは実際の先頭水準を初期表示する。従来の曲線は「平均パラメータの仮想曲線」と明記し、各水準の確率の平均と区別する。確率値の計算と既存の選択IDは維持し、CSVには`CurveBasis`を追加した。報告用の既定曲線も仮想曲線であることをタイトルへ残す。
- **カテゴリー診断への接続：** カテゴリー別の根拠表、尺度機能の要約、判断支援、再採点候補は仮想曲線を判定に使わず、実際の水準を調べる。観測件数・fitは従来どおり全体集計であり、水準別の観測fitを計算したと誤表示しない。新たに確認対象・優位になる範囲がない水準・ピークが低い水準・水準間の最小ピークを記録する。PCM／GPCMの根拠表から意味の定まらない単一の`CurvePeakTheta`等を除き、RSMの共有尺度では従来の単一曲線列も保持した。`CurveMinimumPeakProbability`等が後継列になる。元の診断目安を検証済みの合否閾値には昇格しない。
- **情報量：** 「1評定の情報量」と条件付きSEを表示する。他のパラメータを固定した1応答の`a² Var(X)`で、評価計画全体の情報量や推定した受検者位置の標本SEではない。詳細数値は閉じた欄へ移し、CSVにも`InformationBasis=single_response_conditional`を記録した。根拠は前節で確認した[MurakiのETS研究報告](https://files.eric.ed.gov/fulltext/ED390896.pdf)の項目情報とテスト情報の区別である。
- **ZSTD API：** `zstd_transform`として`wilson_hilferty`（補正項あり）、`wilson_hilferty_no_offset`、`linear_normal`を明示的に選べる。診断生成・自由度方式の再計算・カテゴリーfitへ渡す。旧`whexact=False/True`はそれぞれ従来のWH／線形近似を維持し、FACETSの同名Yes／Noとして扱わない。旧Trueと別の名前の併用はエラーにする。FACETS Yes／Noとの対応は[公式説明](https://www.winsteps.com/facetman64/whexact.htm)との式の照合であり、FACETS実機照合ではない。
- **実際の設定を図へ残す：** ZSTDのPathwayには自由度方式・変換・打切りを表示し、Web保存図と配布版の図に保持する。比較用の列を併記する`both`でも主表示の値は従来方式・打切りなしであることを区別する。設定がない旧結果では「未記録」とし、±9を推測して表示しない。従来方式のfit表にも方式のメタデータを追加したが、既定の数値は変えない。
- **内部文言：** 図の画面・保存HTMLの確認事項は利用者向けの表へ変換し、技術的な行番号・資格コードを表面に出さない。保留は理由と利用上の制約を普通の文で説明する。配布版HTMLの「prototype」と保存資格辞書の生表示を除き、図の詳細表と今回変更したfit設定の説明を日英で整えた。再現用JSON／CSVの安定した列名・理由コードまで削除する方針ではない。全パネル・全報告経路の内部文言除去完了とは扱わない。

### ConQuestの現在環境と、揃えた条件

`/Applications/ConQuest/ConQuest`は`arch -x86_64`で起動し、**ConQuest 5.47.5 Standard Version**を出力して正常終了した。アプリのライセンスや設定は変更していない。

既存の32人・256評定、3つの2水準ファセット、得点0–3、単位重み、完全データのPCMを再利用した。母集団の平均0・SD1.1を固定し、非ステップ要因は中心化、ステップ要因の位置と閾値の制約・符号・水準順を明示的に変換した。積分は両側とも−20～20・401点・間隔0.1の有限格子へ揃えた。

ConQuestで校正値を固定した実行と再推定を各1回行い、入力応答、モデル行列、固定分散、返却座標と最後の履歴の一致を確認した。Python側の再推定は同じ有限格子上の現行尤度・勾配を使う開発用L-BFGS-Bであり、通常アプリの既定GH／EM経路の性能を認定したものではない。

| 確認対象 | 結果 | 言える範囲 |
| --- | --- | --- |
| 同じ校正点の負の対数尤度 | 現行Pythonと独立有限和の差は固定点5.68e−14、ConQuest返却点0 | 正規化・同定・符号の対応。この値自体は2つのPython評価の差で、ネイティブ尤度は別に照合した |
| ConQuestの記録した負の対数尤度 | 固定299.533224、再推定299.533206 | 固定点は独立値との既定許容差内。再推定の表示丸め・記録時点と返却点での再評価は区別する |
| 同格子での再推定パラメータ | 最大差7.9413e−7 | この1条件の校正値の一致。ConQuestの出力は小数6桁 |
| Python仕上げ点の勾配 | 最大絶対値3.8791e−7 | この有限格子上の停留性。ConQuest丸め返却点では6.4831e−5 |
| 同じ点の勾配 | 独立式との差は最大3.39e−14 | 固定／返却の2点で座標微分も対応 |
| 表示用のカテゴリー確率 | 最大差3.55e−15 | ConQuest出力のモデル行列を独立式で評価し、8組合せ×401点×4カテゴリー×2校正点で照合。ConQuestが確率を直接出力したとの主張ではない |
| ConQuestのEAP／事後SD | 同じ校正値での決定的積分との差は最大約0.02312／0.00745 | 得点の完全一致は未確認。モンテカルロ安定性は評価していない |

ConQuestのEAP・事後SDは今回`p_nodes=2000, seed=2`のモンテカルロ計算であり、[公式マニュアルの説明](https://conquestmanual.acer.org/s4-00.html)とも対応する。積分法が異なることは確認できるが、観測された差の全量を乱数誤差と断定する反復実験は行っていない。校正の一致から得点や区間の一致まで推論しない。

### 検証・証拠の有効範囲

証拠は`validation/generated/math_ux_repair_20261027/`に保存した。`conquest/`には制御ファイル・入力・ネイティブ出力・履歴・起動情報・比較結果を保持し、入力・実行時ソース・ConQuestバイナリのhashも記録した。

- 最初の関連実行は90通過・3失敗。前節で未修正だったテスト差替え先と保留文期待は修正済みで通過した。今回の3失敗は、方式メタデータが存在しないという旧期待2件と、図の旧保留文期待1件。新しい動作の検証に更新し、該当テストを含む29項目が通過した。
- カテゴリー診断の追加修正では17項目が通過。実際のPCMで、平均パラメータ曲線に優位領域がある一方、一つの水準には存在しない反例を用い、対象水準が診断と再採点候補へ届くことを確認した。
- 最終の日英UI・曲線・同一図保存・配布版関連は60項目通過。工程間に重複があるため合計件数にはしない。AppTestでの検証であり、初心者の理解やブラウザーの全操作を代替しない。
- 前節の指摘に対応する組込みカテゴリー自己検証も通過。旧ZSTDの96条件は修正前と数値が完全一致した。確率・予測子・JMLE／MML尤度／勾配・パラメータ展開・情報行列・条件付きSEの9関数は修正前とASTが一致した。方式選択と記録を変更したfit関数は、この同一性の対象から除いた。
- 保存済み配布版プロジェクトから推定0でHTMLを生成し、今回指定した内部語が本文に残らないことを確認した。追加の検証スクリプトでは旧関数抽出時に既定定数を渡し忘れた`NameError`が1回発生し、検証側の名前空間を修正して完了した。アプリやConQuestの計算失敗ではない。

- 日英AppTestで保存した実際のPlotlyデータから、カテゴリー・情報量・ZSTD Pathwayの3 PNGを幅1100pxで作成し、図注記と軸名の非重複を目視した。カテゴリーの端のラベル欠けを修正し、最後の17項目も通過した。合成表示の確認であり、狭い画面・実利用者の受入ではない。最初のPNG生成は制限環境でChromeが起動直後に終了し、承認された実行環境で作成した。確認用データのZSTD列不足による1回の出力抽出失敗も記録し、該当列を補った。成功済みの情報量PNGは再生成していない。

新たなmfrmr実行、自由SD／GPCM等の外部再推定比較、一般的なSE・CIの性能保証、停止中研究の再開、公開Web更新、配布パッケージ再ビルドは行っていない。既存の保留を維持する。旧画面・報告用図セット・3画面構成への統合と、残る内部語・英語の整理は既存UX1の続行範囲である。

<a id="bias-em-model-comparison-20261027"></a>

## 探索的バイアス・EM・モデル比較の点検 — 2026-10-27

**3つとも既存の基盤を活用する。新しい推定器を増やすより、何を調べる機能なのか、比較に使える結果なのかを接続することが先である。** 今回は、未完了の再推定が情報量基準の比較へ入る不具合と、バイアス表の重み合計を件数と表示する不整合を修正した。以下の将来仕様を実装済みとは扱わず、UX1／既存数値計画の範囲・推論保留を維持する。

### 1. FACETSを参考にした探索的バイアス

問いは「全体として厳しい評価者か」ではなく、例えば「ある評価者が、特定の課題に限って全体のモデルから予測されるより低く採点しているか」である。FACETSは主分析の推定値を固定して、観測された組合せごとの局所的なずれを推定する。[公式Table 9](https://www.winsteps.com/facetman64/table9.htm)と[Table 13](https://www.winsteps.com/facetman64/table13.htm)にこの手順と得点単位／logit単位の差が示されている。

現行の`estimate_bias_interaction`も、個人位置・主効果・段階・傾きを固定し、各セルの予測子に加える1個のずれを条件付き尤度で推定する。ヒートマップ、観測−予測の平均得点差、logit差、条件付きSE・区間、局所LR、Holm補正の基盤がある。MMLでは個人のEAP値を代入しており、バイアスを含む周辺尤度を共同最大化したものではない。構造的な交互作用モデル、群別パラメータの共同推定、校正誤差を含む不確実性とは区別する。

- **実用上の価値：** 評価者×課題等の組合せを選び、差の方向・大きさ・観測件数から原データや採点基準を点検する。全体の厳しさの順位だけでは見つからない局所パターンを示せる。
- **意味の違い：** `FixedChiSq`は現行Pythonでは中心化した記述的な散らばりで、FACETS Table 13の「すべて0」を検定する統計量と同じではない。自由度・p値は未提供。キーに`Profile CI`があっても、現行区間は他の推定量を再最適化したprofile likelihood区間ではない。FACETSのt／自由度をそのまま再現する機能でもない。
- **推論の限界：** Holm補正は、元のp値の妥当性や同じ標本での探索・選択の問題を解決しない。疎なセル、極端反応、校正誤差、EAP代入、共通する個人・評価者による依存を考慮した誤検出率・検出力・区間被覆は、既存の独立式照合とは別の検証である。バイアスの検出を不公正の原因判定や評定の自動削除に結び付けない。
- **今回の表示修正：** 最初の表を「2要因の水準・評定数・観測−予測（点）・バイアスの大きさ（logit）」にした。件数は正の重みを持つ実評定数を優先し、旧記録では保存された件数を使用する。重み合計しかない記録は「重みの合計」と表示し、件数を推測しない。ヒートマップのhoverも同じ意味へ接続した。欠測セルを0へ置換せず、元の表・CSVの数値を保持する。Table 14への読み方リンクはセル単位のTable 13へ修正した。

既存の`validation/generated/g4_bias_reference_and_guards_20261004/summary.json`には、決定的な式・制御条件についてPython 20／20、独立R 47／47の通過が保存されている。今回はこの記録を再確認し、バイアス推定関数を変更していないことを照合した。新たなFACETS実機比較、mfrmr比較、独立反復研究を実施したとの主張にはしない。

### 2. ConQuestを参考にしたEMと数値診断

EMはモデルの種類やモデル選択基準ではなく、MMLを計算する方法である。能力分布を積分した尤度を扱うIRTの基盤は[Bock & Aitkin（1981）](https://www.cambridge.org/core/journals/psychometrika/article/abs/marginal-maximum-likelihood-estimation-of-item-parameters-application-of-an-em-algorithm/ABCDBCF23ECBBC837C2F3F38A11ACC7B)にある。今回原著は書誌・抄録を確認した範囲であり、全文を再点検したとは扱わない。現行アプリには`mfrm_em_mml`に加えて直接最適化・hybrid・autoがあり、Webの詳細設定へ接続済み。Eステップで求積点の事後重みを計算し、Mステップで期待対数尤度を制約付きで更新する。

「EMを導入すればConQuestと一致する」という判断はしない。[ConQuest公式コマンド説明](https://conquestmanual.acer.org/s4-00.html)は積分方法・点数、パラメータ変化、deviance変化、反復上限を別々に定め、いずれかの停止条件で終了すると説明している。本アプリでも終了フラグだけで停留性や積分の精度を証明できない。

既存の[クロスエンジン点検](mml_cross_language_development_check.md)には、更新前／後の履歴時点、返却点の勾配、GH次数の影響を分離した結果がある。前節のConQuest実機照合も、固定SD PCMの同一有限格子という限定条件である。今回EM・尤度・勾配・既定求積点を変更しておらず、これらの数値実験を再実行しない。

次の数値受入では、返却された座標での尤度、同じ積分則上の停留性、積分則を精密化したときの変動を分ける。ΔAICが小さいときは、尤度の数値誤差がその差を覆さないかを確認する必要がある。EM後の直接最適化は既存の仕上げ候補を評価する対象であり、一律に精度を保証する新しい既定値にはしない。自由SDで求積点を更新する処理についても、有限点での単調性と連続積分の尤度を混同しない。

### 3. 情報量基準による比較と今回の不具合

現行の比較基盤はRSM／PCM／GPCMを同じ入力から再推定する。基礎式は`AIC = −2ℓ + 2K`、`BIC = −2ℓ + K log(N)`だが、式の計算だけでは比較の妥当性を満たさない。

1. **比較する問い：** RSM→PCMは選んだ要因の水準ごとに段階構造を変える比較、PCM→GPCMはさらに選んだ要因の傾きを変える比較である。「複雑なモデルほど高度」という選択にしない。現在RSMでステップ要因を省略すると、複数水準を持つ最初の要因を自動選択する。要因順がRater／Criterionなら、課題差を調べたい場面でも評価者差を比較し得る。利用者が対象要因を選び、比較する仮定を確認できるUIが必要。**今回この選択UIは未実装**で、自動選択の公開推奨保留を解除しない。
2. **共通の評価対象：** 同じ観測行・得点定義・欠測処理で比較し、JMLEの同時尤度とMMLの周辺尤度を同じ順位表へ混ぜない。重み、アンカー、識別制約、母集団モデル、積分規約、自由パラメータ数を揃えるか、その違いを比較する候補モデルとして明示する。例えば固定SDと自由SDは理論上一律に比較禁止ではなく、同じ周辺尤度で追加パラメータを数えられることが条件。現行の自動再推定は母集団設定を継承し、非単位重み・アンカー・正則化・潜在回帰は対象外である。これらを理論上永久に比較不能とは説明しない。
3. **NとK：** 現行MML BICのNは独立な個人の周辺尤度の数であり、評定の総件数ではない。固定した構造パラメータ次元・個人あたり応答数が有界という適用範囲を保持する。EAP得点を個別の自由パラメータとしてKへ足さない。JMLEは個人パラメータが標本とともに増える問題があり、通常のIC解釈をそのまま認定しない。汎用の小標本AICc式はMFRM用の導出がなく、引き続き数値を出さない。
4. **数値と結論：** 終了・境界・停留性・積分誤差を確認する。AIC／BICの差や重みを真のモデルの確率と断定しない。AICとBICが異なる候補を支持することは、[ConQuestのIRTree例（2.15.5節）](https://conquestmanual.acer.org/s2-00.html)でも示される。LRTの参照分布は別の条件が必要で、GPCMの傾き制約・境界、JMLEの付随パラメータ等に対する現行の保留を維持する。ICの低さは絶対的な適合や公平性の証明ではない。

**再現した不具合：** 再推定が結果辞書を返せば、反復上限終了でも`refit_ok`となっていた。状態を分離した反例では、未完了PCMのAIC=194が、完了RSMの256・GPCMの258より低いため、PCMが内部の最上位・AIC重み約1になった。公開推奨は既に保留されていたが、差・重み・順位・LRの記録まで正しいわけではなかった。

修正後は最適化の成功フラグが明示的な真で、K・対数尤度・AIC・BICが有限の候補だけを差・重み・内部順位・LRへ含める。未完了・状態不明・非有限の候補は元の値と終了メッセージを保存し、比較へ入れない。2候補未満なら差・重み・順位を作らない。単なる終了確認を推論資格や停留性の合格へ昇格しない。

古い比較キャッシュは新しい除外規則を適用していないため再表示せず、利用者が再度「比較する」を押すまで再推定しない。他の任意分析が共有するキャッシュキーは変えない。画面には日英の状態・件数・パラメータ数・尤度・AIC／BICと差を示し、内部の状態コード・資格コード・AICc未提供コードを直接並べない。CSVには再現に必要な詳細を保持する。対象外理由の一部英語やバイアス詳細表等、残る内部向け表示をすべて除去したとは扱わない。

局所バイアス分析の尤度改善を主分析のMML尤度へ足してAICを作る接続はしない。交互作用モデルとして比較するには、同じ尤度内で必要なパラメータを共同推定し、識別・K・不確実性を定義する必要がある。探索で見つけた交互作用を同じデータで採用する場合の選択の影響も別に扱う。

### 4. Web／配布版への仕分けと受入

| 機能 | Streamlit Web | Tkinter配布版 |
| --- | --- | --- |
| 探索的バイアス | 結果から任意に開く。要因の組合せ、ヒートマップ、得点差・logit差・件数を中心にする。p値・条件付き区間・補正は詳細へ収める | 同じ共有計算と表示の意味を用い、選んだ組合せ・図・根拠・保留をプロジェクトに保存。初期版への接続は未実装 |
| EMによるMML | 「推定方法」の必要な詳細として提供。最初のCSV読込からEMを選ばせない。終了状態と精度の問題は結果に接続する | 同じMML APIを使用。長い反復の進行・中断・結果保持を確認してから追加。別エンジンを新設しない |
| モデル比較 | 調べたい要因と少数の候補を明示し、利用者の操作で実行。比較できない行を選択候補にしない。自動的な最適モデル宣言は保留 | 同一条件での比較履歴・設定・実行状態を保存。複数の比較や積分感度の反復をローカルで扱えることに価値を置く |

配布版の現在の初期範囲はRSM／JMLEの基本分析・保存・図であり、上表のMML・比較・バイアスがすでにTkinterへ接続済みという意味ではない。採用時には実行中断と保存・再読込を含めて受入する。Web／配布版で数理の正しさや推論資格を分けず、実行量・保存・作業の継続性を分ける。

独自性は他ソフトの表を増やすことより、**入力した評価配置から局所的な問題を見つけ、Wright／Pathway等の図と根拠に戻り、同じ分析条件で保存・報告できること**に置く。FACETSやConQuestと同名の欄があるだけで同等としない。条件付きバイアスと共同モデル、終了と精度、探索と検証を区別して渡すことを受入基準にする。

### 今回の検証と残る範囲

証拠は`validation/generated/bias_model_comparison_review_20261027/`に保持する。修正前スナップショット、未完了候補が最上位になる反例、実際の短反復テストの終了メッセージ、差分、検証ログを保存した。

- 初回の関連テストは58通過・5失敗。失敗した期待は、25／30反復の実データ生成フィクスチャがすべて反復上限に達しているのに比較・順位を要求していたものだった。反復上限を延長せず、元の有限ICを保持して比較から除外する期待へ修正した。完了した合成結果の差・重み・LR算術は別の明示的なテストで確認する。
- 修正と重み表示の追加テスト後は47通過・4失敗。4失敗は件数の整数表示「2／3」をテスト側が「2.000／3.000」と期待したためで、数値・列選択は通過していた。テストを実際の整数表示へ修正し、該当するJMLE／MML×件数情報3種類の6項目が通過した。工程は重複しており、件数を合算して単一実行の成功とはしない。
- 明示的な成功／失敗／状態不明、成功欄が文字列の旧記録、非有限値、1候補だけ完了、旧キャッシュ、日英の状態表示、重み合計しかない旧記録を確認した。表示変更による再推定0、元の結果・CSVの保持も確認している。
- 差分のAST照合で変更したトップレベル関数はモデル比較・2つのバイアス表示等の7関数に限定される。バイアス推定・尤度・勾配・推定エンジンへ新しい数値処理は加えていない。数値エンジンの再認定、ブラウザー全操作、初学者の理解、公開、配布版ビルド、外部ソフトの新規実行、停止中研究の再開は今回の実施範囲ではない。

<a id="model-comparison-facet-workflow-20261027"></a>

## モデル比較の対象要因・再利用・保存を接続 — 2026-10-27

**前節で残した、比較対象の要因を利用者が選ぶUIを実装した。** 対象は既存の任意ツール「モデル比較」であり、初心者の基本分析へ新しい必須設定を増やしていない。推定アルゴリズム、停止条件、求積既定値、モデル選択・LRT等の推論保留は変更しない。

### 操作と比較する仮定

RSMで複数の候補要因がある場合、「比較する要因」は未選択から始まり、「比較する」は選択後に有効になる。データで複数水準が観測された要因だけを候補にし、Person・単一水準・存在しない列を受け付けない。候補が1つだけなら、その要因を表示して使用する。元のPCM／GPCMに段階の要因が指定されていれば、その明示された設定を初期表示に使う。

GPCMの傾きの要因は、閉じた「比較条件」で選べる。元のGPCMの指定があれば維持し、それ以外は段階と同じ要因を初期値とする。画面には常に「RSMは共通の段階、PCMは選択した要因ごとの段階、GPCMはさらに指定した要因ごとの傾き」という実際の比較条件を表示する。段階と傾きが別の要因である設定を、黙って同一要因へ変更しない。

`compute_model_choice_comparison`にも`step_facet`／`slope_facet`を接続した。プログラムからの呼出しでも、複数候補のあるRSMで対象を省略した場合は、先頭要因を推測せず比較不可を返す。現行PCM／GPCMの設定と、唯一の候補は省略時にも利用できる。組込み自己検証と既存の呼出しテストも、比較する要因を指定する形へ更新した。旧APIで曖昧なRSM呼出しが自動実行された挙動は、意図的に維持しない。

### 別のモデルを混ぜず、再計算を必要な操作に限定する

「同じモデル名」だけでは既存結果を再利用しない。RSMは対象要因を変更しても共有段階モデルのままなので再利用できる。PCMは段階の要因、GPCMは段階・傾きの両方が元の分析と一致するときだけ再利用する。一致しなければ現在と同名の候補も再推定し、表では「再推定」と表示する。

比較のキャッシュを元の分析IDと段階・傾きの選択にひも付けた。他の任意分析が使う共有キャッシュキーは変更していない。要因の切替だけでは推定せず、「比較する」を押してから実行する。以前に計算した組合せへ戻ればその結果を再表示し、要因ごとの傾きの選択も保持する。別の分析・別の条件・古い方式のキャッシュは表示しない。

CSVにはモデルごとの`StepFacet`／`SlopeFacet`、全セクション共通の`ComparisonStepFacet`／`ComparisonSlopeFacet`／`SourceAnalysisID`を保存する。数値表・LR計算記録・利用上の制約が、同じ分析と選択条件へ戻れるようにした。現在の資格判定を表示・保存時に再適用し、キャッシュ内の古い推奨や資格で保留を解除できないようにした。これはWebセッション内の保持とCSV保存であり、Tkinterのプロジェクト保存へ接続済みという意味ではない。

### 主表を簡素化し、保存経路を揃える

主表はモデル・計算状況・AIC／BIC・それぞれの差に絞った。対象要因・件数・パラメータ数・対数尤度は「計算の詳細」、数理的な注意は「解釈上の注意」に収める。常時表示する保留通知は1つに統合し、値のない検定表や繰り返しの推奨保留文を並べない。非対応の重み・アンカー・正則化・潜在回帰等の説明も日英で表示し、内部の作業用語を出さない。非対応の設定を解除すれば科学的に適切になると誘導する説明は避けた。

ブラウザーで点検したところ、従来のグリッドには選択条件や保留を含まない独自のCSV保存があった。3モデルの表を標準の静的な表へ変更し、条件・終了状態・保留を含む専用CSV保存へ統一した。計算の詳細も同じ扱いにした。画面上の丸めとCSVの保存精度は従来どおり区別する。

### 検証の範囲

証拠は`validation/generated/model_comparison_facets_20261027/`へ保存した。

- 対象要因・再利用・キャッシュ・CSVを含む初回の関連テストは**76項目通過**。実際の小規模JMLE結果からTask／Raterをそれぞれ指定し、4つのPCM／GPCM再推定の設定まで届くこと、元の結果を変更しないことを確認した。これは経路の検証であり、モデル選択の精度や外部ソフトとの同等性の検証ではない。
- 主表と詳細の分離・非対応理由の日英表示後は**73項目通過・2項目選択対象外**。静的な表への変更後は、その影響を受けるUI・保存・非対応理由の**24項目通過・20項目選択対象外**。工程間に重複があり、合算しない。
- 未選択・不正な要因・単一水準での実行防止、段階と傾きの別指定、現在のPCM／GPCMと異なる仕様の再推定、要因を戻した際の再利用、日英切替、別分析への切替、古い／別条件のキャッシュ、保存した資格を確認した。表示・選択操作による自動再推定は0。
- ローカルChromeで独立した表示確認ページを使用した。入力は合成の状態制御であり、実データの数値受入ではない。初期未選択、選択後の実行、主表と折りたたみを確認した。確認用プロセスが古いimportを保持したため、ページ再読込だけで新しい表にならず、そのプロセスを再起動して最終画面を確認した。アプリの推定失敗ではない。
- 最終の1200×897表示で6列の主表、1つの保留通知、閉じた詳細を目視し、グリッド独自の保存ボタンがなくなったことをDOMでも確認した。専用CSVのボタンを操作したが、既定Downloadsに新しい該当ファイルを確認できず、ブラウザー経由の実保存完了・保存場所は未確認。CSVの内容と選択条件・保留の一致はAppTestによる検証と区別して記録する。
- 最終差分の点検で、不正な旧設定を要因設定より先に拒否する順序も維持した。`config`がNone・リスト・文字列の反例と正常な選択・再利用を含む15項目が通過。推定エンジンのファイルは着手前と同一で、同じ条件を渡した場合の尤度・勾配・EMの数値処理は変更していない。

公開Web更新、Tkinter実装・ビルド、EMの精度改善、外部ソフトの新しい実行、停止中研究の再開は行っていない。3画面の外枠や他の診断・出力の整理、初学者の理解、公開環境での受入は既存UX1の未完了範囲として保持する。

<a id="model-comparison-readable-export-20261027"></a>

## モデル比較の読みやすいCSV・実保存・画面復帰 — 2026-10-27

前節で未確認だったブラウザーからの実保存を、独立した検証用ブラウザーで確認した。同時に、通常の保存が数値・検定記録・内部の資格コードを混在させた横長の表になる問題を修正した。目的は既存UX1の「選んだ条件の結果を理解して保存できること」であり、モデル比較の推論資格の追加ではない。

### 保存内容と画面

- **表を保存（CSV）：** 日英の項目名・計算状況を使った3モデル×15列。モデルごとの段階・傾きの要因、件数、BICの標本数、パラメータ数、尤度、AIC／BICと差、推定方法、元の分析ID、解釈上の保留を含む。内部の状態コードや作業番号は出さない。日本語の読込みを考慮してUTF-8 BOMを付けるが、画面に文字コードの説明は追加しない。
- **計算記録（CSV）：** 閉じた「計算の詳細」の中に配置した。従来のモデル比較・LR計算記録・資格判定のセクションとコードを保持し、推定方法を追加する。通常表の保存とは用途を分ける。検証例では従来8行×40列、変更後は8行×41列。
- **画面の数値：** 小数3桁までで末尾のゼロを省き、値がない差は「—」とする。丸める前の表からCSVを作るため、保存値の精度は落とさない。英語の主表が不要なゼロで窮屈になった実画面を修正し、日英とも主表6列を維持した。
- **保存と復帰：** 両方のダウンロードで`on_click="ignore"`を使う。別の画面へ移って戻る、別の分析へ移って戻る、言語を変える操作で、対応する条件・結果を保持する。新しい比較を選択だけで開始しない。

### 検証と証拠の有効範囲

記録は`validation/generated/model_comparison_exports_20261027/`に保存した。既存のChrome for Testing 151.0.7922.34と対応ドライバーを、一時プロファイル・専用保存先で使用した。利用者のChrome 156のプロファイル・設定は変更していない。前節の通常Chromeで保存先を確認できなかった原因の解決や、公開環境の受入を意味しない。

1. **変更前の実保存：** 既存の計算記録をネイティブのダウンロード操作で取得し、8行×40列を確認した。保存時の推定は0だったが、画面処理は1回増えた。`baseline_download_check.json`に保存先・ハッシュ・処理回数を記録した。
2. **変更後の実保存：** 日英それぞれ通常表と計算記録の計4ファイルを取得した。数値列が元の計算記録と同じ精度であること、選択したTask／Rater、分析ID、推定方法、保留を照合した。4回とも保存操作前後の画面処理・再推定の増加は0。別画面との往復でも選択条件・結果を保持した。実ファイル、`final_browser_check.json`、成功ログを保存した。
3. **最終の表示変更：** 上記保存確認後は数値の表示文字列だけを変更した。日英とも738pxの表領域に6列が収まることをDOMと画像で確認した（`final_display_check.json`、`compact_ja.png`、`compact_en.png`）。CSV生成前の値は変えておらず、成功済み4ファイルのダウンロードは繰り返していない。ソースの段階差と各検証の有効範囲を`verification.json`に記録する。
4. **自動テスト：** 保存・表示・条件保持の初回11項目、表示変更後の最終11項目が通過し、各40項目は選択対象外。工程間で同じ項目を含むため合算しない。未完了PCMの差を欠測のまま保存する反例、元の結果・キャッシュを変更しないこと、画面・別分析からの復帰、日英のラベルと保留も確認した。
5. **確認手順の失敗：** 最初の自動操作で再描画中の要素を参照する失敗、次に選択変更の再描画前に実行を押す失敗があった。検証スクリプトが再描画・選択反映を待つよう修正し、失敗ログも保持した。アプリの推定失敗とは区別し、待機上限の延長はしていない。

ブラウザー入力と再推定応答は合成の状態制御であり、実データの推定・EM精度・モデル選択性能の証拠ではない。今回の実行コード変更は比較画面の表示・保存関数に限定し、推定エンジンは着手前と同一。既存のR／ConQuest照合を再実行せず、推論保留を維持した。

公開Web更新、Tkinterプロジェクトへの比較保存、配布物ビルド、3画面全体の受入、初学者の理解、通常ブラウザーの全保存環境は未完了として区別する。

<a id="information-criteria-validation-20261027"></a>

## AIC／BICの独立検証：計算の一致と比較可能性 — 2026-10-27

**結論：識別できる固定SD MMLの限定条件では、計算式・自由パラメータ数・独立ソフトとの値の一致を確認した。しかし、通常設定の積分誤差と、要因を分離できないデータを比較へ含める問題が残っていた。後者を修正し、前者を理由の一つとしてモデル選択の既存保留を維持する。** テストの成功やConQuestとの1例の一致によって、推論資格を変更しない。

### 問いと検証方法

今回の明示的な依頼はAIC・BICの検証である。G2のうち、式、K、標本単位、尤度の数値精度、同じデータのモデル間差を対象とした。停止中の正式研究やモデル選択率の大規模シミュレーションは再開していない。

1. 実装と既存テストから、制約後のK、MMLの積分対象、未完了候補・重み・アンカー・正則化等の扱いを確認した。
2. 保存済みの**合成データ32人×8評定=256行**を利用。各2水準のRater／Task／Criterion、4カテゴリー、固定母集団SD=1.1、Criterionを非中心化要因とした。PCMの段階とGPCMの段階・傾きはすべてCriterion。実調査データではない。
3. 現行の比較APIを通して、通常EM・15点求積・最大400反復・相対許容値1e−6でRSM／PCM／GPCMを推定した。現在結果の再利用と実際の候補再推定を含む。
4. アプリの確率・パラメータ展開関数を使わないPython実装で、同じ15点上の尤度を再計算。別にゼロ初期値・81点求積の有限差分最適化を行い、161点と適応積分で精度を確認した。積分範囲は標準正規尺度の±10／±12、事前に設定した参照値照合許容差は1e−7。局所勾配の上限1e−4と正の曲率も確認した。大域最大の証明ではない。
5. Rの別実装による連続積分で通常推定点・独立推定点の計6点を照合し、`stats::AIC`／`stats::BIC`を計算。ConQuest 5.47.5では同じデータ・固定SD・401点等間隔求積でRSMを新規実行し、既存の同条件PCM実行記録を再利用した。成功済みPCMは再実行していない。

再実行可能なPythonコードは`validation/information_criteria_check.py`、設定・入力識別・全数値・R／ConQuest出力は`validation/generated/information_criteria_validation_20261027/`に保存した。失敗を別試行で置換したり、結果に合わせて許容差を変更したりしていない。

### 式・K・標本単位

`AIC = −2ℓ + 2K`、`BIC = −2ℓ + K log(N)`を確認した。Kは表示される推定値の行数ではなく、制約後に独立に調整できるパラメータの数である。[Akaike（1974）の書誌・公開抄録](https://doi.org/10.1109/TAC.1974.1100705)、[Schwarz（1978）の公開要旨](https://doi.org/10.1214/aos/1176344136)、[R公式の定義・比較条件](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/AIC.html)と照合した。今回、原著PDF全文の取得はタイムアウト等で完了しておらず、全文を再査読したという記録にはしない。

検証例では、要因の位置の自由度はRaterの1＋Taskの1＋Criterionの2＝4。4カテゴリーの3段階には和0の制約があるため、RSMは2、PCMは2水準×2＝4、GPCMはさらに幾何平均1の傾きの自由度1が加わる。したがってKは**6／8／9**。MMLでは32人のEAPをKに加えない。自由SDの設定ではSDの1自由度を追加することも既存実装とテストで確認したが、今回の独立最適化は固定SDに限定した。

MMLの尤度は個人ごとの応答ベクトルについて積分し、その対数を足す。この標本のBICでは独立な個人の数**N=32**を使い、256評定を独立な標本数とはしない。これは固定した構造パラメータ次元・個人あたり応答数が有界で、独立性・正則性等の条件を満たす場合の単位の判断である。人数を使えばあらゆる割当や依存構造でBICが妥当になる、という結論ではない。

JMLEでは人数とともに個人パラメータも増えるため、通常のAIC／BICによる自動推奨は引き続き保留。JMLEの同時尤度とMMLの周辺尤度を同じ順位表へ混ぜない。AICcについても、他ソフトが出力していることだけを本アプリの小標本補正の根拠にせず、現行の未提供扱いを維持する。

### 数値結果

| モデル | K | 通常設定のAIC | 独立参照AIC | 通常設定のBIC | 独立参照BIC |
| --- | ---: | ---: | ---: | ---: | ---: |
| RSM | 6 | 620.140234 | 620.471061 | 628.934650 | 629.265476 |
| PCM | 8 | 614.730518 | 615.066412 | 626.456405 | 626.792299 |
| GPCM | 9 | 615.538373 | 615.972024 | 628.729996 | 629.163647 |

- 通常設定の返却点で、アプリと独立実装の**同一15点求積の尤度差は3モデルとも0**だった。式や保存値の取り違えではなく、積分と最適化の精度を分けて評価すべき差である。
- 独立PythonとRの連続積分・ICの最大差は**1.14e−12**。ConQuestのRSM／PCMのKも6／8で、AIC／BICの最大差は**3.95e−6**、印字の丸め内だった。ConQuestのPCM−RSM差はAIC=−5.40465、BIC=−2.47318。GPCMのConQuest実機照合を実施したという意味ではない。
- 81点で独立最適化した点を161点で点検した勾配最大絶対値は7.64e−6、1.81e−5、1.15e−5。数値Hessianの最小固有値は11.12、8.12、6.36で正、差分幅変更による相対差は最大5.46e−8だった。161点と適応積分のNLL差は最大2.40e−11。81点の目的関数自体には最大1.67e−6の差が残るため、表の独立参照ICには適応積分値を使った。
- 通常推定点を連続積分で評価し直すと、ICはRSM／PCM／GPCMで**+0.333306／+0.340824／+0.444987**変化した。その点から精密な参照推定点への移動による変化は**−0.002479／−0.004930／−0.011336**。今回の差の大部分は求積の誤差だった。
- この1例では、AIC・BICが最小なのはいずれもPCMだった。ただしGPCM−PCMのΔAICは通常設定の0.807855から参照値0.905612へ変わる。モデル間差の変化は最大約0.098であり、単なる全モデル共通の定数差として消せない。順位が変わらなかったことを、小さな差の比較が常に安定する根拠にはしない。

### 観測配置による識別不成立の反例と修正

同じ入力から、R1はT1、R2はT2の評定だけを残す128行の反例を構成した。全員が両方の組合せに回答していても、評価者効果と課題効果は分離できない。RSMの名目Kは6だったが、独立に構成したカテゴリー対数確率比のデザイン行列の階数は5。評価者の自由座標へ+0.7、課題へ−0.7を加えてもNLL=157.18061194156968は全く変わらなかった。

修正前はEMが正常終了し、このRSMを含む3モデルが差・重みの計算へ入った。公開推奨自体は保留されていたが、計算終了だけで技術比較へ含める条件としては不十分だった。

比較前に、保存された厳密な要因制約と実際の評定の組合せからMMLの要因位置の階数を確認し、既知の欠損を拒否するよう修正した。JMLEで既に検出された個人・要因の識別不成立も同じ入口で拒否する。古い比較キャッシュより先に判定し、日英で「要因の効果を分離できない」ことを説明する。元の推定値・K・保存済み結果は書き換えず、要因を自動削除したりKだけを減らしたりしない。修正後の同じ評定配置・制約では、再推定0で比較を停止した。

このガードは**必要条件の一部**である。段階・傾き・母集団SDを含む全面的な識別、有限な最大値の存在、境界、積分精度を認定しない。旧記録に制約情報がない場合は新しいMML階数検査を行えず、既存の推論保留を維持する。観測水準と保存設定の不一致、階数検査が対応規模（行×列500万）を超える場合も比較を保留する。資格記録にあった「パラメータ数は識別済み」という一律の文言も、座標の制約と観測配置の識別を区別する記述へ直した。

### 完了範囲と残る条件

数理・IC単位・自由SDの計数等の既存テストは**60項目通過**。ガード追加後の比較・保存・状態・識別の関連テストは**100項目通過**。最終の空データ／列欠落・資格説明変更後は影響する**17項目通過、10項目選択対象外**。工程間に重複があり合算しない。UIは日英AppTestによる局所確認で、新しいブラウザー受入や公開受入ではない。

数値検証はガード追加前のソースで実行し、着手時スナップショット・ハッシュを保存した。以後の実行コード差分は比較入口の拒否判定と資格記録の説明のみで、推定エンジンは不変。成功済み推定・外部照合は再実行せず、その有効範囲を`verification.json`に記録する。

現時点で追加確認が必要なのは、**利用する比較について、同じ観測・尤度・制約、識別と境界、返却点の停留性、積分を精密化したときのIC差の安定性を満たすこと**である。単に求積点数を一律に増やすことや、ConQuestと同じ式を使うことを資格判定の代用にしない。AICとBICは異なる候補を支持することもあり、その相違は[ConQuest公式例の2.15.5節](https://conquestmanual.acer.org/s2-00.html)にも見られる。両者の一致を「真のモデル」の証明にしない。

今回の1データと反例は、自由SD・他の配置・小標本でのモデル選択性能、JMLEの一般的妥当性、LRTの参照分布、GPCMのConQuest同等性を保証しない。mfrmr候補との新しい直接照合、公開Web更新、Tkinter接続・ビルド、推論保留解除、停止中研究の再開は行っていない。

<a id="ic-sensitivity-workflow-20261027"></a>

## AIC／BICの数値点検を既存の比較画面へ接続 — 2026-10-27

**前節の積分誤差を、比較した実際の推定点で点検・保存できるようにした。数値点検は再推定やモデル選択の認定を行わない。** 現行UX1／G2のうち、終了フラグと数値精度を区別し、同じ分析条件の記録を画面から保存する接続である。新しい比較モデルや正式研究は追加していない。

### 操作・計算・保存の対応

既存の閉じた「計算の詳細」に「数値精度を確認」を配置した。対応する比較を実行した際に、3候補の実際の推定座標・元データ・設定・終了状態を保持する。この段階では数値点検を実行しない。明示的なボタン操作で、各座標を固定したまま元の点数q、max(61, 2q+1)、max(121, 4q+1)のGauss–Hermite求積を評価する。元の尤度との再現差が1e−8を超えた場合や、入力・固定SD・制約が一致しない場合は結果を表示しない。

当面の対象は固定母集団SDのMML、単位重み、1,000評定以下、元の求積点数1〜121。潜在回帰・アンカー・正則化は対象外で、既存の比較入口の制限も適用する。上限はWebでの実行範囲であり、その範囲の全データの精度保証ではない。共通モジュール`mfrm_app/ic_sensitivity.py`に処理を置き、Web以外で異なる計算を作らない。Tkinterへの操作・保存の接続を完了したという意味ではない。

画面はモデル、確認点数、AIC・BICの変化、勾配の最大絶対値の4列。K・標本数・推定点を固定するため、この点検によるAICとBICの変化量は等しい。CSVは9行の評価値と6行のモデル対・基準別の差を含み、点数、K、BICの個人数、固定SD、終了状態、元の分析・候補・推定点・入力の識別情報、選択要因、日英の限界説明を保存する。再推定・連続積分照合・推論資格はいずれもFalseとし、通常の比較表・CSVの値は置き換えない。

点検のキャッシュをデータ・設定・座標に結び付け、別条件や資格の誤った保存状態から復活させない。画面復帰・言語変更・CSV保存で再評価せず、CSV保存による画面の再処理も防ぐ。旧バージョンの比較キャッシュは新しい座標記録を持たないため、明示的な「比較する」が必要で、自動再推定はしない。

### 既存の独立参照で確かめたこと

前節の32人・256評定の合成データと通常EMの保存推定点を復元し、点検モジュールを実行した。前回の入力・座標・連続積分結果とエンジンのハッシュを照合してから再利用し、推定、R、ConQuestは再実行していない。

| モデル | 15→121点でのIC変化 | 121点での勾配最大絶対値 | 保存済み連続積分とのNLL差 |
| --- | ---: | ---: | ---: |
| RSM | +0.333305829 | 0.151478 | −1.95e−9 |
| PCM | +0.340824024 | 0.182386 | +3.47e−9 |
| GPCM | +0.444986908 | 0.330696 | +8.27e−10 |

事前許容差1e−7以内で、前節の同じ推定点における積分感度を検出できた。モデル間の差の変化は最大0.111681。前節の約0.098は精密な参照点へ**再推定した後**の差であり、今回の座標固定の点検値と混同しない。勾配は制約後の自由座標で評価するため尺度に依存し、一律の合否閾値は設けない。

今回の121点はこの例の連続積分参照と一致したが、Webの点検自体には独立連続積分を含めていない。2つの求積結果の差が小さくても誤差上限ではない。また、精密な尤度を評価してもその値を最大化したことにはならない。この理由で、点検表を「補正済みAIC・BIC」として主表に採用せず、モデル選択の保留も解除しない。

### 検証の工程と限界

証拠は`validation/generated/ic_sensitivity_workflow_20261027/`に保存した。

- 初回の関連テストは**62項目通過**。独立した尤度・勾配の式との照合、入力や点の変更、異なる固定SD、元の尤度不一致、非対応範囲、日英の明示操作・保存・復帰を含む。
- 終了状態を文字列から成功へ読み替えない処理と、資格を誤って持つキャッシュの拒否を追加後、影響する数値点検・比較・識別の**51項目通過**。工程間に重複があり合算しない。エンジンと資格ポリシーは着手前と同一で、成功済み推定・外部照合を繰り返していない。
- 専用Chrome for Testingを一時プロファイル・保存先で実行。前節の保存推定点を復元した実際の比較表示から、日本語・英語の15行CSVをネイティブに保存し、値・識別情報・限界・Falseの資格を確認した。両保存で画面再処理・数値点検・再推定の増加は0。
- 同じブラウザー工程の画面移動で、検証スクリプトのXPathがラベルに一致せず失敗した。アプリの計算失敗ではない。成功済みダウンロードと失敗ログを保持し、ラベル選択だけ修正して、ダウンロードを省いた別工程で言語変更・画面往復を確認した。各新規セッションで準備の数値点検は1回、言語変更・復帰による追加点検と全工程の再推定は0。単一工程の一括成功とは扱わない。
- 最終の日英の数値点検表は、折りたたみ内の704pxの領域に4列が収まった。日本語画面の画像と英語の表示・DOMを確認した。これは専用の局所ページによる受入であり、通常ブラウザーの全環境や公開サイト、初心者の理解の検証ではない。

推定エンジン・既定設定・主表のIC・推論資格は変更していない。自由SD／JMLE等の数値点検、境界を含む全パラメータの識別、精密な尤度での再推定、モデル選択性能、Tkinterの画面・保存、公開・配布物更新は未完了であり、本点検の成功をそれらへ拡張しない。

<a id="conquest-ic-sensitivity-20261027"></a>

## ConQuestで積分感度と同じ求積での最適化を直接照合 — 2026-10-27

**RSM／PCMでは、ConQuestでも同じ推定点の15→121点の積分感度を再現した。同じ15点の尤度を厳しい条件で最大化した結果も、独立Python実装と一致した。** 前節までの「差の主因は積分近似」という説明を、ソフトウェア間の実行で直接確認した。この例では通常Pythonの終了点にもIC換算で約0.0012の改善余地があったが、約0.33〜0.34の積分による変化より小さい。一般的なソフト間同等性やモデル選択性能の認定ではない。

### 問い・比較条件・再利用

利用者の「ConQuestでも同じ結果か。異なるならアルゴリズムの違いか」という問いへの続行。G2の既存AIC／BIC点検の範囲で、積分による評価値の変化と、同じ近似尤度で推定値を更新する影響を区別した。32人・256評定、正規母集団SD=1.1固定、同じ要因制約・カテゴリー・単位重みの合成データを使用。RSM／PCMを対象とし、新しい入力条件や正式研究は追加しない。

1. 各モデルの通常Pythonの保存推定点をConQuestへ固定値として渡し、Gauss–Hermiteの15点と121点で尤度を評価した。推定点を変えないので、変化は求積条件の影響として切り分けられる。
2. 同じ点を初期値に、ConQuestで15点のまま再推定した。終了条件はパラメータ変化1e−9または逸脱度変化1e−11、反復上限2000。ConQuestの全既定設定による実行ではない。
3. アプリの確率・展開関数を呼ばない既存の独立Python式でも、ゼロ初期値・15点・L-BFGS-Bで推定した。異なる更新法が同じ近似尤度の解へ到達するかを確認する。大域最大の証明ではない。
4. 前回の401点ConQuestと連続積分の結果は再実行せず、保存出力・データ・ソース・ConQuest実行ファイルのハッシュを照合して再利用した。

ConQuest公式の[推定設定](https://conquestmanual.acer.org/s4-00.html#estimate)は、積分方法・点数とEMの終了条件を別に指定する。固定値の指定は[importの仕様](https://conquestmanual.acer.org/s4-00.html#import)に従った。通常PythonもConQuestもMML／EMの経路を持つので、EMの採用有無では今回の差を説明できない。EMと積分近似の関係は[公式技術説明](https://conquestmanual.acer.org/s3-00.html)を参照する。

### 同じ推定点での積分感度

| モデル | Pythonの15→121点のIC変化 | ConQuestの15→121点のIC変化 |
| --- | ---: | ---: |
| RSM | +0.333305829 | +0.33331 |
| PCM | +0.340824024 | +0.34082 |

ConQuestの固定点4実行と独立Pythonの逸脱度（−2対数尤度）の差は最大4.10e−6で、ConQuestの小数5桁出力の丸め内だった。事前照合許容差は1e−5で、結果に合わせて変更していない。

**全パラメータを固定したConQuest実行のKは0になる。** その出力のAIC／BICを、元の推定モデルのICと直接比較していない。まず逸脱度を照合し、元のモデルのK=6／8とN=32を用いて評価値と変化量を作った。これらは固定点での感度診断であり、精密な尤度を最大化したモデルのAIC／BICとは区別する。

### 同じ15点で十分に最適化した結果

| モデル・指標 | 通常Python・15点 | ConQuest・15点・厳しい終了条件 | 保存済みConQuest・401点 |
| --- | ---: | ---: | ---: |
| RSM・AIC | 620.140234 | 620.13904 | 620.47106 |
| RSM・BIC | 628.934650 | 628.93345 | 629.26548 |
| PCM・AIC | 614.730518 | 614.72935 | 615.06641 |
| PCM・BIC | 626.456405 | 626.45523 | 626.79230 |

15点でのConQuestと独立Python最適化の逸脱度差は最大4.84e−6、自由座標の最大差は8.11e−7。独立Pythonの勾配最大絶対値はRSMで3.47e−7、PCMで5.76e−6で、事前条件を満たした。ConQuestは46／45反復で逸脱度の終了条件に達した。表の通常Python結果からの改善は両モデルともIC換算で約0.0012であり、「通常の終了点から同じ15点の尤度をさらに改善できた」という証拠になる。

一方、ConQuestを15点で精密に終了させても、401点の結果への差は残る。終了条件だけを厳しくすれば積分誤差も消える、という説明は成立しない。この例では積分精度と終了点の確認を別に行う必要がある。RSM／PCMの2候補では、AIC・BICが小さいモデルはどの比較経路でもPCMだったが、選択性能や正しいモデルの判定を実証したとは扱わない。

### 証拠と未確認範囲

`validation/generated/conquest_ic_sensitivity_20261027/`の`protocol.json`、`check.py`、`review.py`、6実行の制御・出力、`summary.json`、`comparison.csv`、`pair_differences.csv`、`verification.json`に保存した。6実行とも正常終了し、実際の求積点数・個人数32・母集団平均0／分散1.21・ネイティブのKとIC算術を確認。モデル行列と採点済みデータも保存済み401点実行と一致した。固定座標の出力丸めと、最適化後の座標差・勾配は別に照合した。

今回の新規計算はConQuestの固定点4実行・15点再推定2実行と、独立Pythonの15点最適化2実行。アプリ本体・推定エンジンを変更しておらず、成功済みアプリテスト・401点推定・連続積分・R照合は再実行していない。GPCMのConQuest、自由SD、他のデータ配置、全既定設定、推論やモデル選択性能は未検証のまま。今回の結果を根拠とする既定点数の一律変更、推論保留解除、公開・配布物更新、停止中研究の再開は行わない。

<a id="em-stopping-record-20261027"></a>

## EMの停止条件を比較し、実際の判定値を保存 — 2026-10-27

**既存PythonのEMでも、同じ15点の尤度を厳しく最適化すると、ConQuestと一致した独立参照へ到達した。今回必要な修正は、計算法の置換ではなく、停止条件と実際の判定値を追える記録・説明の補完である。** 前節で残った約0.0012の改善余地を、積分精度とは分けて確認した。

### 終了条件の意味と追加計算の役割

PythonのEMは`abs(LL_new − LL_old) / (abs(LL_old) + 1e−10) < reltol`で終了する。今回のConQuestは絶対的な逸脱度変化1e−11またはパラメータ変化1e−9を終了条件に使う。逸脱度は−2LLなので、LLの大きさが約300の場合、ConQuestの逸脱度条件に対応する相対変化は概ね1.7e−14である。1e−6や1e−10という数値を、絶対許容値とそのまま比較しない。ConQuestの別の終了条件や初期値もあるため、この換算だけで実行経路全体を同一としない。

最初に、前節の通常Python保存点から、固定SD=1.1・同じ15点求積で、相対許容値を1e−10にしたRSM／PCMのEMを実行した。次に絶対・相対条件の対応を確認するため、1e−14を別条件として事前記録し実行した。Mステップや更新式は変更しない。独立参照とConQuestは前節の保存済み15点結果を再利用し、新しい反復研究は行っていない。

| モデル | EMの相対許容値 | 保存点からの反復数 | 独立15点参照とのNLL差 | 終了点の勾配最大絶対値 |
| --- | ---: | ---: | ---: | ---: |
| RSM | 1e−10 | 24 | +6.80e−8 | 1.01e−3 |
| PCM | 1e−10 | 24 | +7.24e−8 | 9.55e−4 |
| RSM | 1e−14 | 51 | +5.80e−12 | 8.40e−6 |
| PCM | 1e−14 | 64 | +5.85e−12 | 8.73e−6 |

事前の工学的な照合目標はNLL差1e−7、勾配最大絶対値1e−4。1e−10の条件はNLL差を満たしたが、勾配の目標を満たしていない。結果を削除したり合格へ読み替えたりしていない。絶対許容値との対応を調べた別条件1e−14は両方を満たし、自由座標の参照との差は最大5.79e−7だった。これは同じ有限15点の目的関数での結果であり、前節の約0.33〜0.34の積分近似による差を解消したことにはならない。

### 記録とHelpの修正

従来は終了メッセージと勾配等を保存していたが、実際に終了判定へ使った相対変化とMステップの実効許容値が記録にない。また、EMの尤度履歴には終了点を最後に重複して保存するため、末尾2値の差を取ると0となり、本来の停止判定値を復元できない。既存履歴は変更せず、実際の判定時点の値を別に保存する。

- `EMRelativeLogLikTolerance`：外側のEMで指定された相対許容値。
- `EMStoppingRelativeLogLikChange`：正常終了の判定に実際に使った相対変化。未完了・Mステップ失敗では欠測とし、最後の重複値から0を作らない。
- `EMMStepOptimizer`、`EMMStepFunctionTolerance`、`EMMStepGradientTolerance`、`EMMStepMaxIterations`：実際に渡した内部最適化の設定。通常L-BFGS-Bでは目的関数許容値1e−9、勾配許容値1e−5、反復上限50。呼出しに渡す元のftol=1e−7ではなく、共通設定関数による制限後の1e−9を記録する。方式にない許容値は欠測とする。

共通エンジンの結果から既存の収束表へ接続し、閉じた詳細画面・CSV・シリアライズ後の結果で保持する。古い結果に設定を推測して補わず、hybridの最終結果へEM準備段階の停止条件を混ぜない。主画面のボタンや警告は追加していない。日英の設定Helpと用語集から「小さい値ほど精密な推定になる」という保証的な説明を修正し、EMの相対尤度変化と積分点数を区別した。

推定値の更新、尤度、終了判定式、Mステップの設定、既定値は変更していない。数値処理の関数は不変で、変更したEM関数・収束表も記録追加を除くASTが着手前と一致する。この根拠で成功済みの数値計算と外部照合を再利用し、記録変更のために再推定を繰り返さない。

### 検証と適用限界

証拠は`validation/generated/em_precision_followup_20261027/`へ保存。新規計算4実行の全経路、事前目標、保存参照の識別、変更前後のソース、`numerical_source_identity.json`、テストログを含む。

初回の関連テストは30通過・11失敗。11件は既存テストが`streamlit_app.minimize`を差し替えても、移設済みの`mfrm_app.engine.minimize`へ不具合注入が届かない問題だった。参照先を直し、故障・非有限・目的関数悪化・制約違反・途中失敗を含む対象テストとHelpの21項目が通過した。期待する失敗条件を弱めていない。追加の日英詳細表示は2項目通過・4項目選択対象外。件数は工程間で重複し合算しない。

記録追加の検証には、固定SD／自由SDの小規模EM、正常終了の実際の相対変化、反復上限時の欠測、旧記録への補完防止、pickle再読込、CSVの丸め前の値、推論保留を含む。日英AppTestで保存済み記録の表示と表示時の推定0を確認した。新しいブラウザー実保存、公開・配布物更新、GPCMの外部照合や自由SDの精度検証は行っていない。1e−14を一般的な推奨値や既定値にせず、今回の限定条件での説明・検証結果として保持する。

<a id="ic-refinement-workflow-20261027"></a>

## 積分点を増やす任意の再推定を接続 — 2026-10-27

**同じ推定点での感度点検に続いて、細かい求積で推定値を更新する操作を追加した。保存済み32人・256評定では、独立した高精度参照とAIC／BICが最大8.30e−9の差で一致した。元の分析結果を置き換えず、モデル選択の保留も維持する。** 前節までの点検で、終了条件を厳しくするだけでは積分近似の差が残ること、点を固定した再評価だけでは精密な尤度の最大値にならないことが分かっていた。今回の追加は、この未接続の操作を少数の候補比較へ結び付けるもの。正式研究の再開やエンジンの置換ではない。

### 実装と判断の境界

既存の閉じた「計算の詳細」で「数値精度を確認」を実行した後、「積分点を増やして再推定」を任意で使う。初回表示・言語変更・画面復帰・保存で自動起動しない。主表と通常の比較CSVは元の結果を保持する。新しい表はモデル・計算状況・AIC・BIC・追加点検差・勾配の6列とし、専用CSVに座標・元の分析ID・選択した要因・終了状態・実効最適化設定・失敗理由・保留を残す。

`mfrm_app/ic_sensitivity.py`から既存の制約付き直接MMLを呼び出す。新しい尤度式・別エンジンは作らない。対象は従来の感度点検と同じ固定SD MML、単位重み、1,000評定以下、元の求積点数121以下。アンカー・潜在回帰・正則化を含めず、RSM／PCM／GPCMで同じ入力・SD・求積条件であることを検証する。3候補すべての元の尤度を再現してから最適化を開始する。

元の点数をqとし、max(121, 4q+1)点で元の推定座標を初期値として再推定、その2倍+1点で同じ終了座標を再評価する。通常の15点なら121点・243点となる。既存の方式選択に従い、今回RSM／PCMはL-BFGS-B、GPCMはSLSQPを使い、最大400反復・要求許容値1e−12と実効optionsを記録する。これらはこの任意操作の条件であり、通常分析の既定値ではない。SLSQPのftolをEMの相対尤度変化と同じ意味にはしない。

厳密な終了フラグ、座標と評価値の有限性、自由・従属の傾き制約、初期点に対する目的関数の非悪化を確認する。失敗・未完了の試行も保持し、いずれかが該当する組の比較差は欠測にする。境界付近ではその旨を示し、生の勾配だけで制約下の最適性を認定しない。「終了」は最適化の処理状態であって、積分精度・局所／大域最大・識別・選択性能の合格を意味しない。追加点検差が小さくても誤差上限として扱わず、連続積分との照合をアプリが実行したことにもしない。

### 数値の照合

今回の新規最適化は、保存済み合成入力に対するアプリの3実行のみ。独立Pythonの高精度最適化・連続積分、R、ConQuestの成功済み結果を再計算せず、保存した参照と入力・実装のハッシュを確認して使った。再推定後の座標では、別実装の式による243点の尤度も新たに照合した。

| モデル | 121点再推定後AIC | BIC（N=32） | 元のICからの変化 | 243点でのIC変化 | 243点の勾配最大絶対値 |
| --- | ---: | ---: | ---: | ---: | ---: |
| RSM（K=6） | 620.471060632 | 629.265476048 | +0.330826536 | +2.08e−9 | 2.85e−6 |
| PCM（K=8） | 615.066411756 | 626.792298979 | +0.335894041 | −8.29e−9 | 1.38e−5 |
| GPCM（K=9） | 615.972023637 | 629.163646762 | +0.433651001 | −3.22e−9 | 3.23e−8 |

3候補とも処理が終了し、制約内部で、事前に固定した照合目標を満たした。目標は同じ座標での独立NLL差1e−9、高精度参照とのIC差1e−5、印字済みConQuestとのIC差1e−4、追加求積のIC差1e−6、勾配最大絶対値1e−4。独立243点評価とのNLL差は出力精度で0、独立した保存済み高精度参照とのIC差は最大8.30e−9、RSM／PCMのConQuestとのIC差は最大3.96e−6だった。全3再推定と追加評価の処理時間は、この環境・この入力で18.14秒。別環境・大きな入力の所要時間の保証にはしない。

この入力では引き続きPCMのAIC／BICが最小である。GPCMは独立Python参照と照合したが、ネイティブConQuestのGPCMとは照合していない。固定SD・同じ尤度での数値一致を、自由SD・JMLE、モデル仮定の妥当性、モデル選択の性能へ広げない。

### 保存・画面・検証証拠

`validation/generated/ic_refinement_workflow_20261027/`に事前目標、3候補の全結果、座標、比較差、参照との照合、変更前後のソース、テストとブラウザー証拠を保存した。関連88テストが通過。独立式との照合、元の結果の非変更、実行前の3候補の整合、未完了・真偽値でない終了状態・非有限・座標次元・従属傾きの制約違反・目的関数悪化・例外・勾配次元の反例、境界表示、日英の操作・保存・キャッシュを含む。

専用Chrome for Testingで日英のCSV各1ファイルを実保存し、追加6列が704px幅に収まり、言語変更・画面復帰・保存による追加の再推定要求が0であることを確認した。ブラウザーでは今回保存した実際の数値結果を操作後に返す構成を使い、完了した3最適化は繰り返していない。数値工程とUI工程を組み合わせた確認であり、ブラウザー内で3最適化を再実行した通し試験とは区別する。最初のローカルサーバー起動はsandboxのbind制限で停止したが、許可された127.0.0.1限定の起動後に画面・保存の工程は通過した。

共通推定エンジンと資格判定は着手前と同一で、通常の設定・主表・推論保留は変更していない。Webはこの限定した任意操作を支え、Tkinterは同じ計算モジュールを使えるが、画面・プロジェクト保存への接続は未実装。公開・配布物更新、外部ソフトの新規実行、自由SDや他の配置の性能研究、停止中の正式研究は行わない。初学者本人の理解と公開環境での受入も未確認である。

<a id="ic-refinement-outcomes-20261027"></a>

### 再推定結果の区別と未完了理由の続行点検

前項の保存済み画面を点検すると、画面復帰後は実行ボタンが消え、再推定表の直前に見出しが残らなかった。通常の比較と同じAIC／BIC列を持つ別の表であるため、表を識別する「再推定」の短い見出しを残す。また、「未完了」だけでは終了条件の未達成と目的関数の悪化を区別できず、計算失敗の理由も技術CSVの生のエラー文字列にしかない場合があった。

計算結果に理由を記録し、終了条件の達成未確認、尤度の悪化、不正な推定座標、傾きの制約違反、尤度／勾配の評価失敗、その他の計算エラーを区別した。画面は6列を維持し、問題がある候補だけに日英の短い理由と「この候補との比較差は計算していない」旨を表示する。元の内部理由・エラーは技術記録として保持し、画面には表示しない。CSVへ同じ翻訳済みの計算状況と説明を追加し、画面から保存後の読み取りへ意味を維持する。

境界付近の説明には該当モデル名を入れ、範囲外の失敗と範囲内の境界接近を混ぜない。理由を持たない旧記録を読み直した場合は、不明な理由を推測したり再推定したりせず、詳しい理由が残っていないと明示する。正常な候補に新たな注意文は追加しない。

関連82項目が通過した。前項の88項目とは重複があり合算しない。既存の数値・比較・保存テストに加え、失敗種別の理由、全候補失敗、旧記録、境界、日英の画面とCSVの一致、内部エラーの非表示、画面復帰後の見出しを確認した。専用Chrome for Testingでは、保存済み正常結果と明示的に模擬した未完了・制約違反を表示し、日英2ファイルを実保存。6列は738pxの領域に収まり、言語変更・画面復帰・保存時の推定呼出しは0だった。この工程は保存結果の表示検証であり、実際に失敗するデータをアップロードして再推定する通し試験ではない。

証拠は前項の保存先内の`usability_followup/`に残した。理由の記録を除くモジュール全体のAST、共通エンジン、資格判定の同一性を確認し、数値・制御処理が変わっていないことを根拠に前項の高精度参照・ConQuest照合を再利用した。新しい最適化・外部ソフト実行は0。推定値・比較への算入条件・元の分析・推論保留は変更していない。内部ロードマップへ反映したが、公開・配布物更新や実利用者の理解の受入は行っていない。

<a id="history-ic-guard-20261027"></a>

## 分析履歴のΔAIC／ΔBICとBIC表示を整合 — 2026-10-27

**モデル比較の終了・数値点検とは別に、履歴が条件不一致や未完了の記録を一つの最小値に混ぜる経路が残っていた。履歴と保存結果の比較が使う共通関数を修正し、条件や終了を確認できない表では差を出さない。** 元の値は残し、差を省いた理由を日英の画面とCSVに結び付ける。

### 問い、反例と修正範囲

前節まではRSM／PCM／GPCMの候補比較を検証したが、`build_run_history_likelihood_table`は各指標の全行最小値から無条件に差を作り、その後に条件一致のフラグを付けていた。フラグがFalseでも差は表示・保存される。条件キーに固定SD・SD推定モード・求積点数・潜在回帰の実体がなく、元データを持たない同士も同じキーになり得た。終了状態も差の分母に使っていなかった。

保存済みの32人・256評定の結果を基に、片方の記録のAICを20小さくする模擬状態を作成した。これは数値性能の新規実験ではなく、表示・比較ロジックの反例である。固定SD変更、異なる入力、未完了、文字列の終了フラグ、元データ欠落の5例で、修正前にはいずれもΔAIC／ΔBICが20と0になった。入力が異なる例では条件一致Falseなのに差が出ており、他の4例では条件一致までTrueになっていた。

[R公式AIC文書](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/AIC.html)は、同じデータでの比較と尤度の最大化を前提にし、異なる尤度の定数にも注意を求めている。今回の修正はこれらを取り違えないための保護であり、終了フラグだけで最大化を認定しない。**履歴では同一条件の記録を保守的に照合する。異なるモデル・固定値・積分点数の比較が統計的に一律不可能だという意味ではない。** 適切な候補間のモデル比較と積分感度は既存の専用操作で扱い、モデル選択の利用保留は別に維持する。

### 実装

- 正規化後の実データ（得点・重み・対象者・要因）、モデル・推定法、得点範囲、制約と符号、母集団モデル、SD固定／推定、求積点数等を照合する。必要な記録が欠けている場合は、欠落値同士を一致として扱わない。固定SDの値は条件に含める一方、自由SDの推定値自体は条件へ入れず、推定モードと境界を使う。反復上限・許容値・通常EM／直接法の違いだけで異なる尤度と扱わない。
- 表全体で2件以上、条件既知・一致、厳密な終了フラグ、有限のAIC／BIC／逸脱度を確認できる場合だけ最小値との差を計算する。他の場合は差を欠測にし、異なる条件を混ぜた指標の平均等も出さない。正則化は履歴差の対応範囲外として、記録不足とは区別する。複数の条件が混在する履歴から任意の組を自動抽出せず、既存の保存結果選択で2件を選べる設計を維持する。
- 差を省く理由を、件数不足・記録不足・条件不一致・終了未確認・非有限・対応範囲外に分ける。画面の表示とCSVの説明を一致させ、条件一致を精度・モデル選択の資格へ昇格しない。
- 履歴の概要表も詳細表と同じ補正済みBICを使う。古いMMLの保存値は`SavedBIC`に保持する。対象者数の根拠が欠ける旧概要は、評定数で代替せずBICを欠測にする。今後の概要保存にはBICの標本数とその根拠を残し、正しく計算できる情報を失わせない。元の保存オブジェクトは書き換えない。
- 履歴CSVの保存を画面再処理なしに変更した。利用者向け説明から「理由コード」の保存を求める内部用語を除き、残る説明内容を示す表現へ変えた。

### 検証と限界

`validation/generated/history_ic_guard_20261027/`に変更前後の5反例、同一条件の正例、ソース、テスト、ブラウザー記録を保存した。修正後は5反例すべての差が欠測となり、既知・同一条件・終了確認の正例では20と0を保持した。設定・表示の保護であり、新しい最適化やConQuest実行による数値研究ではない。

初回の関連86項目、追加修正後の対象30項目、保存値・欠けた旧結果の5項目、正則化の説明1項目と既存の組込み自己テストが通過した。件数は重複があり合算しない。入力・重み・固定SD・自由SD・求積・母集団モデル・アンカー・符号・モデル・推定法、終了False／不明／文字列、非有限値、単独の記録、旧概要と新概要を含む。表示テストでは概要・詳細・CSVのBIC、差を省く理由、保留、保存値の非変更、表示時推定0を確認した。既存の小規模推定を含む回帰テストは実施したが、成功済みの独立高精度参照や外部ソフト照合は繰り返していない。

専用ブラウザーでは日英の4場面（未完了、固定SD不一致、条件一致、旧概要）とCSV各1ファイルを確認し、画面復帰・保存で追加推定0、CSV保存による再処理0だった。初回は旧形式を模擬する入力に現行のBIC標本数情報が残り、元値の保存に関する検証が停止した。製品側の誤りとは切り分け、記録を残して模擬入力を直し、未完了の保存工程と残る操作を確認した。最終の正則化用理由と説明文の変更は対象テスト・文言照合で確認し、影響しない4場面の保存は繰り返していない。

共通エンジン、AIC／BICの計算式、推論資格、前節までの精密な再推定機能は不変である。公開・配布物更新、停止中研究の再開、新しいConQuest／R実行、初心者本人の理解の受入は行っていない。履歴の横長の詳細表を含む全面的なUI整理も、この局所修正の完了範囲とは分ける。

<a id="beginner-ui-simplification-20261027"></a>

## 初心者向けStreamlit画面の説明と詳細表示を整理 — 2026-10-27

**説明の追加ではなく、操作・結果・必要な制限を先に見せる構成へ戻した。** 現行UX1と報告・学習要件の「入力→実行→解釈→保存」を対象とする。直近のIC機能だけを磨く作業へ限定せず、通常の入力画面と保存済み作文評価例の結果画面を実際に確認した。推定法の追加や推論資格の変更は行わない。

### 問題と修正

- 標準設定でもモデル・推定法・分析範囲の選択肢が常時並んでいた。「分析設定」にまとめ、条件の変更と各ヘルプは保持した。主な実行ボタンの直前には現在のモデル・推定法・分析範囲を短く表示する。入力確認と架空データの説明は短くし、実データの送信先・機密データの扱いは残した。
- 「最初に確認」には同じ根拠が状態表示、根拠欄、概要表、全項目表に重複していた。状態・実際の注意点・次に開く診断を主表示にし、全項目・利用範囲・凡例・CSVを一つの「確認の詳細」に集約した。他の結果ページは短い状態表示と選択した内容を中心にし、同じ誘導ボタンと長い説明を繰り返さない。
- 結果表は内部の計算方法列がInfit／Outfitより先に現れ、通常画面にも収束・構造的識別性の技術表が続いていた。初期表示の列をファセット・水準・推定値・Infit・Outfit・件数（存在する列）に絞り、数値の表示桁を指定した。列の選択表示だけを変更し、完全な表と資格の列は保持する。SE・CI等は列選択で確認できる。当初記載した「個別表」にはSE・CIのない表もあり、再点検で訂正した（下記の10月28日追補を参照）。収束の詳細と長い読み方は閉じた欄から開く。推論保留と全て最低点／最高点の受験者の制限は、短い警告として欄の外に残す。
- 実際に確認事項CSVを保存した後、結果への切替ボタンがStreamlit上部バーに隠れた。固定表示を内側の短いコンテナに指定していたため、スクロールを追従していなかった。実際のラッパーへ指定を移し、スクロール後のネイティブクリックと幅640pxの画面で操作を確認した。
- 同じ意味の総称警告と個別警告を重ねず、具体的なエラー・注意点を重大度ごとにまとめた。条件判定や警告内容、入力エラーによる実行拒否は変えていない。Helpの案内先も変更後の名称と配置に合わせた。

### 検証と証拠の範囲

`validation/generated/beginner_ui_simplification_20261027/`に変更前後の画面、差分、ログ、日英CSV、検証記録を保存した。既存の作文評価960評定の保存結果を再利用し、数値研究のための新しい推定を加えていない。日英の結果・図・保存画面、詳細の開閉、設定変更とHelp・言語切替を確認した。確認事項CSVは7項目の根拠・次の操作・状態を保持し、ネイティブ保存時の画面再実行は0、ブラウザー工程全体の推定は0で、分析IDも維持した。最終の結果表の表示変更後は影響する画面だけを再確認し、成功済みCSV保存は再実行していない。

初回の失敗には旧UIの期待値、Helpの旧名称、共通エンジン移設前の差替え先、実行中のテストファイル編集に伴うAppTestのソース行対応不整合が含まれた。原因を切り分けて修正後、関連バッチ166項目が通過し、残る3項目は更新後の対象18項目で通過した。結果の技術詳細を整理した後の26項目と、最終の境界表示2項目・表表示2項目も通過した。工程間には重複があり、件数を合算しない。入力不正の実行拒否、未知・未完了の推定、保留した区間の記録、極端得点の受験者、日英移動、Help、列の保持を対象とした。既存の小さな実推定を含む回帰テストと、保存結果だけを使うブラウザー検証は区別する。

共通エンジン・出力資格・IC計算のソース同一性を確認し、前節までの独立積分・ConQuestの証拠を再利用する。UIテストで数値性能やモデル選択の利用可否を新たに認定しない。公開版への反映、Tkinterの画面・配布物、実際の初心者による操作・理解の受入は今回の完了範囲に含めない。通常画面から技術情報を隠したことだけで、学習上の理解や適切な研究判断を保証しない。


<a id="beginner-ui-corrections-20261028"></a>

## 初心者向け結果表の資格表示と案内先を修正 — 2026-10-28

前節の簡素化には、情報を保持していることと利用者に伝わることを混同した不足があった。再点検で、全評定が最低点・最高点の受験者にも内部の有限値が通常の推定値として表示され、未収束時の主要ボタンが報告書作成へ進むことを再現した。列メニューには多数の内部名が並び、件数欄は実際の`N_Fit`と一致しない`N`を探していた。前回テストの警告の存在・元表の保持だけでは、該当行の識別と適切な移動先を保証できていなかった。

### 修正の範囲

- 基本表は表示用の列だけを渡し、「基本」「精度」で切り替える。列メニューから内部名が現れることを防ぎ、元の数値と判定列は別の「計算記録」に保持する。画面往復・言語切替で選択を維持する。
- Personの報告用推定値と既存の境界判定、分析全体の保留を表示へ反映する。全て最低点／最高点はその行に示し、推定不可・解釈保留の推定値、適合度と区間を基本／精度表で空欄にする。他の受験者とファセットまで一律に隠さない。元の保存結果・計算記録の値は変更しない。数値を文字列へ変えず、欠測はStreamlit標準の`placeholder`で「—」と表示する。
- 精度表は既存の方法・状態が確認できるSEと記録済み区間を使い、欠けた区間を画面側で自動生成しない。条件付き近似と事後分布の近似を区別し、根拠・水準が不明な区間は出さない。これは現在の資格に従う表示であり、被覆率や差の推論を新たに認定するものではない。
- `N_Fit`を表示する。計算定義まで照合すると重み合計であるため、重み列を指定した分析では「重み合計」と表示し、小数を整数へ丸めない。列ヘルプに評定数となる条件を残す。
- 収束の確認先を「計算の詳細」へ接続し、主要ボタンで開いたときは結果画面の先頭に開いて表示する。通常の結果閲覧では閉じた欄のままとする。日本語・英語の案内と確認事項の移動先も揃えた。

### 検証と残る範囲

証拠は`validation/generated/beginner_ui_corrections_20261028/`。初回は112項目通過・旧移動先を期待する1項目が失敗し、期待値を更新した。追補の表示・境界・終了状態・資格テスト68項目、最高点の境界1項目、重み付き件数等の表示3項目も通過した。工程間には重複がある。欠測の最終表示は対象テストと実画面で確認し、成功済みの独立数値比較を再実行しない。既存の終了状態回帰テストでは、意図的な終了失敗を作る反復上限1の小さなJMLEを1回実行した。ブラウザー表示時の推定0および新しい独立研究なしとは区別する。

保存済み960評定の結果で日英の基本／精度表示、列メニュー、画面往復を確認した。未収束を模擬した日英の操作では計算の詳細へ進み、報告書作成を呼ばない。既存の境界判定fixtureに与えた終端値を用いて、最低点行の欠測と状態、通常行の数値保持を確認した。これは新たに推定した境界実験ではない。800px幅の限定表示も確認したが、全操作・キーボード・実利用者の受入ではない。

検証の途中で、選択ボタンの属性とAppTestのプロトコルを誤って想定した検証コードを修正した。ブラウザーに読込済みの古いモジュールが残っていたため、検証サーバーを再起動し、最終のレンダラーを読み込んだ記録と実画面を照合した。欠測表示には[Streamlit 1.54の公式仕様](https://docs.streamlit.io/1.54.0/develop/api-reference/data/st.dataframe)と手元の1.66で確認した標準の`placeholder`を使い、依存版の引上げや数値の文字列化を避けた。最終表示の対象3項目が通過し、古いモジュールによる画面を最終修正の証拠にしない。失敗ログも保持する。

表内蔵のCSVボタンは検証用Chromeでファイル取得を確認できず、ダウンロード履歴も空だった。原因は特定しておらず、ブラウザーのポリシーによる拒否とは断定しない。画面へ渡す表の保留値・列・元表の非変更は確認したが、このCSVの実保存は未確認として残す。既存の出力画面の完全なCSV・計算エンジン・出力資格モジュール・報告用共通モジュールは変更していない。公開、Tkinterの配布物更新、ConQuest／Rの再実行、初心者本人の理解の受入は行っていない。


<a id="beginner-csv-recovery-20261028"></a>

## 表のCSV保存と、計算状況からの復帰を点検 — 2026-10-28

前節で未確認だった表内蔵CSVの保存は、ローカルの画面ありChromeで保存ファイルを読み戻して確認した。今回の画面なしChromeでは標準の`showSaveFilePicker`が`AbortError`を返していた。Streamlit 1.66の表はこのAPIを先に使い、中止時には終了するため、通常のダウンロード履歴や指定ダウンロードフォルダーだけを待っても保存を確認できない。アプリの保存失敗や組織ポリシーの拒否と同一視しない。

### 保存内容の照合と修正

- 保存APIが返したファイルハンドルから`getFile()`で内容を読み戻し、日本語の基本表と日英の精度表を照合した。既存の境界判定fixtureに与えた終端値を用いる3行の表示テストであり、新規に推定した統計実験ではない。推定不可のPerson行の推定値・適合度・SE・区間は空欄、行の状態は保持し、通常のPersonとRaterの値を維持した。技術列は基本／精度表のCSVへ混入しない。
- 元の結果との照合で、表示用に100倍した区間水準が同じ`CI_Level`名でCSVへ流れる不整合を確認した。表のデータは0.95のまま保持し、画面だけ百分率表示にした。[既存の最低依存版1.54の標準書式](https://docs.streamlit.io/1.54.0/develop/api-reference/data/st.column_config/st.column_config.numbercolumn)で対応でき、数値の文字列化や依存版の引上げは不要だった。影響する日英の精度表は修正後に再保存・読戻しし、既に通過した基本表は繰り返さなかった。
- 同じ既定ファイル名を使う2回目の保存要求が完了しなかった工程は残した。以後は検証用の別名を標準の保存ダイアログへ渡し、既存ファイルを上書きせず確認した。同名上書き操作の受入や、全ブラウザー・Windowsでの検証とは区別する。保存処理そのものを別実装へ置き換えたり、新たな保存ボタンを増やしたりしていない。

### 計算状況の表示

「次に確認」で計算の詳細を開いたとき、停止条件未達／終了記録不明／方法固有の推論保留を短く示す。入力を見直す場合は「入力と設定を確認」から既存の開始画面へ移動できる。設定・保存した結果は維持し、自動で再推定しない。技術的な表・終了メッセージ・勾配・識別性等は「計算記録と診断」にまとめ、必要時に開く。

同じ終了状態の警告を、ページ上部・操作の要約・計算画面で重ねない。計算画面が当該状態を表示している場合に限って重複を省き、他の画面へ戻ると上部の警告を戻す。アンカー・接続・残差等の別の問題は省かない。終了条件達成を、MMLやJMLE境界の推論保留解除へ読み替えない。

### 検証と範囲

証拠は`validation/generated/beginner_ui_save_recovery_20261028/`。表示・状態・復帰の対象30項目、要約の重複調整後11項目、上部警告を含めた調整後6項目、区間水準の修正後3項目が通過した。工程間で重複するため合算しない。False／不明の区別、日英表示、設定と結果の保持、推論保留、異なる警告の保持を確認した。数値エンジン・出力資格・共通の報告処理は不変であり、今回は推定やConQuest／Rの計算を実行していない。

日英のローカル画面で、短い状態表示→記録の開閉→入力への復帰を確認した。最終の重複抑制は日本語画面で再確認し、計算画面を離れると上部の警告が戻ることを確かめた。1440pxと800pxの限定表示を確認し、表示・保存の再推定は0。保存操作による画面の再実行も0だった。CSVの最終確認対象は日本語基本表、区間水準を修正した日本語・英語の精度表の3種類である。

CSVの証拠コピーは、保存済みファイルをブラウザーの標準APIでテキストとして読み戻したもの。保存元ファイル名・サイズ・更新時刻・APIの応答は別に記録した。初回の画面なし実行、同名保存で待機した工程、区間水準修正前のCSVも保持し、最終の合格証拠と区別する。公開・配布物更新、初学者本人の理解、全ての技術的な詳細画面の整理は完了範囲に含めない。
