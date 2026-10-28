# Release Checklist

Use this checklist before tagging or pushing a public beta release.

<a id="release-milestones-20261028"></a>

## 現行の公開マイルストーン — 2026-10-28

**品質を満たしてから公開する。公開日は固定しない。配布先はGitHub Releases。**
この節を今回の実行チェックリストの正本とする。下の英語の詳細項目は確認対象・証拠の参照一覧であり、
全項目を毎回ゼロから実行する指示ではない。チェックは成果物・対象ソース・環境・実施範囲の記録を根拠に付ける。
未完了の必須項目を「β」「近似」「未検証」の注記だけで完了にしない。

### M0：公開する範囲を固定する — 計画として完了

- [x] Web更新は、既存の分析・報告の保護と導線改善、固定SD RSM／PCMの指定した水準・差の近似区間を対象にする。
      この区間は数値受入v1を通った対象に限る。有限標本被覆率、自由SD、個人得点、同等性等の資格を追加しない。
- [x] 最初のTkinter配布は、既存の **RSM／JMLE・単一プロジェクト** に固定する。
      CSV入力、列・採点基準、別プロセスの分析と中止、表・3図、`.mfrm`保存／再読込、HTML・PNG・SVGが必須。
      Webと同じ計算・出力資格を使い、失敗・境界・識別不成立を成功にしない。
- [x] TkinterのMML対応は下記X1へ分ける。RSM／JMLE版の配布を「MML対応版の完成」と呼ばない。
      Word／PDF、全Web画面、複数プロジェクトの一括処理、計算途中の再開、自動更新は初回配布の必須条件にしない。
- [x] 対象はmacOS Apple SiliconとWindows 11 x64。OSごとに受入して公開し、両OSの同時公開は条件にしない。
      最低対応OS版・配布可能なOSは実証後に表示する。GitHubは配布先でありOS受入の代替ではない。
- [x] Developer ID署名・公証を当面使わないGitHub配布を前提にする。macOSはアドホック署名と
      Developer ID署名を区別し、未公証と初回起動上の制約を明示する。配布後の通常の起動経路をT4で確認する。
      Windowsも署名の有無を表示する。署名・公証の追加は別の配布改善であり、購入・登録をこの計画から自動実行しない。

### W1：Web公開候補を一つに閉じる — ローカル候補確定、2026-10-28

- [x] 既存の隔離候補・Python 3.12の確認、個人境界・履歴の修正、指定区間の影響範囲の証拠がある。
      [候補](docs/development_integration.md#web-candidate-20261028)、
      [境界・履歴](docs/development_integration.md#person-boundary-history-20261028)、
      [指定区間](docs/mml_quadrature_cross_engine_review.md#fixed-sd-interval-route-20261028)を再利用の起点とする。
- [x] 883ファイルの旧候補からの差分を閉じ、必要な未追跡モジュール・テスト・翻訳・Helpを含む候補を固定した。
      特に `mml_fixed_sd_review.py` と `structural_intervals.py` を取りこぼさない。
      研究生成物・原データ・ローカル計画・秘密設定を混ぜず、公開mainとの差を依存関係ごとに確認する。
- [x] 版表示・CITATION・変更履歴・READMEを候補 `0.2.19-beta` にそろえ、候補ソースと依存環境を特定した。
      未公開のためCITATIONの公開日は未設定。公開mainを基点とした隔離候補を使い、作業元の履歴・未コミット変更を保持する。
- [x] 有効な既存証拠と、新変更で無効になった経路を対応させた。隔離候補のimport／AppTest起動を確認した。
      **完了条件：不足ファイルがなく、レビュー対象と再実行範囲を特定した一つの候補がある。**

888ファイル、環境診断25項目・静的release-check・影響範囲81テストが通過。
[候補と証拠の対応](docs/development_integration.md#web-w1-candidate-20261028)を参照。
静的release-checkのReadyはCI成功や公開承認ではない。新しい区間の実ブラウザー操作・報告実物はW2、CIと公開はW3で確認する。

### W2：その候補の操作・出力を受け入れる — 未完了

- [ ] 通常例と保留例で、入力→推定→日英の結果→図表・報告保存→履歴復帰を確認する。
      既存の有効な操作記録は再利用し、変更がある経路に限定して新しい実ブラウザー確認を行う。
- [ ] 新しい指定区間について、RSM／PCMの提供可、求積等による保留、入力変更後の無効化、
      JSONの実アップロード／保存／再読込を確認する。一般表・個人得点等の保留は維持し、表示操作で再推定しない。
- [ ] 追加された区間表を含むWord／PDF／HTMLの実物で、数値・対象・注意書き・改ページを確認する。
      ファイルが生成できたことだけをレイアウトの合格としない。
- [ ] データ保護、終了失敗、境界、保存結果の取り違え、基本／詳細・日英切替の既存受入を候補へ対応させる。
      人による理解度評価は公開後。数値的な誤表示や操作不能を人による評価待ちで残さない。
      **完了条件：対象経路の既知の重大不具合がなく、各出力の利用条件が同じ候補で一致する。**

### W3：CIから公開先まで閉じる — 未完了、W1・W2に依存

- [ ] 同じ候補commitに必要なLinux／Python 3.11・3.12 CIを一組対応させる。
      現行workflowはintegrationブランチへのpushだけでは自動起動しない。PR／手動実行を選び、重複起動しない。
      必須CIの全体確認は省かず、通過済みのローカル数値研究を文書変更だけで再計算しない。
- [ ] 現在の公開への承認がそろった段階で、main反映・tag・GitHub prerelease・Web更新を実施する。
      今回の計画更新や過去のforce-push承認を、新しい公開・履歴置換の承認に読み替えない。
      公開時にCITATIONの実際の公開日とREADME／変更履歴の候補表示をそろえ、引用の影響範囲を確認する。
- [ ] 公開先の版・commitと、合成データの代表操作・実保存を確認する。不一致なら旧公開へ戻す経路を使う。
      **完了条件：GitHubの公開物と稼働Webが採用した候補に一致し、利用者の入口から保存まで通る。**

### T1〜T6：Tkinter初回配布に残る6つの作業群

基本経路・共通計算・schema 1・worker・図・macOSパッケージの実証は済んでいる。
[試作記録](docs/desktop_prototype.md)はこの証拠であり、一般配布済みの証拠ではない。
T1〜T4はmacOS先行、T5はWindows、T6は受入したOSごとに実施する。

**T1：現在のソースで配布候補を閉じる（未完了、W1の共通ソース確定を利用）**

- [ ] 10月28日の個人境界・出力保護を含め、旧パッケージと異なるソースを取り込む。
      版・OS／CPU・依存版・同梱資源をmanifestと実物で照合する。試作の再梱包を最新修正の反映と取り違えない。
- [ ] 最低対応OSを選んで検証範囲を固定する。アプリ・依存ライブラリの必要なライセンス表示と
      見本・翻訳・フォント等を同梱する。初版に不要な新依存を入れない。

**T2：実GUIの一巡を受け入れる（未完了、T1に依存）**

- [ ] 実ファイルダイアログでCSV入力→列・採点基準→分析→表・3図→HTML／PNG／SVG／プロジェクト保存を通す。
      日本語パス、ID保持、入力不備、通常のキーボード操作と表示範囲を確認する。
- [ ] 中止→再実行、終了失敗、極端得点・保留を確認する。未保存結果の扱い、入力変更と前回結果の区別、
      読み取り専用入力・書込不能先での説明と既存ファイル保全を確認する。表示・読込で再推定しない。

**T3：開発環境への依存を取り除く（未完了、T1に依存）**

- [ ] Python・開発用パッケージを利用しないクリーンな対象環境で、アプリの初回起動、worker、保存・再起動後の復元を確認する。
      開発端末でPATHを変えるだけの試験と、Python未導入環境の受入を区別する。VM・検証協力者等も使える。
- [ ] ネットワークなしで分析・3図・HTML・プロジェクトが利用できることを確認する。
      Web Helpだけはネット接続が必要と説明し、オフラインを全機能の保証へ広げない。

**T4：macOSのダウンロード後の起動を受け入れる（未完了、T1〜T3に依存）**

- [ ] ZIP等の実際に配るアーカイブから、展開・配置・起動・分析・保存を通す。
      公開前はdraft release／検証用artifactで取得操作を確認し、一般公開後のURLはT6で再確認する。
      GitHubからの取得で付く隔離属性を含め、初回のOS表示と必要な操作を確認する。公開前後で対象バイナリのhashを対応させる。
- [ ] Developer IDなし・未公証であること、対応OS、発行元・ファイル照合の方法を配布ページへ明記する。
      OS全体の保護解除や隔離属性の一括削除を通常の導入手順にしない。
      OSが許す個別の確認操作で起動できない環境を、対応済みと表示しない。
      **macOS配布判定：T1〜T4が通り、T6の公開物を同定できること。署名サービスの取得待ちを自動で追加しない。**

**T5：Windows版を作り、Windows 11で受け入れる（未完了）**

- [ ] GitHub ActionsのWindows runnerでビルドし、同梱ランタイムの起動・worker・保存／復元・HTML等を確認する。
      現行Linux CIだけでは不足する。デスクトップの変更でこの検査が対象になるようにする。
- [ ] Windows 11 x64の検証環境を確保し、T2・T3相当と実ダウンロード後の起動表示・日本語パスを確認する。
      検証協力者／VMを含み、開発者の所有は必須にしない。Windows Serverのrunner合格をWindows 11のGUI合格にしない。
      環境が得られるまではビルド検証用成果物とし、Windows 11対応の一般配布は完了にしない。

**T6：受入したOSのファイルをGitHub Releasesへ届ける（未完了、macOSはT4、WindowsはT5に依存）**

- [ ] 対象OS／CPU／最低版、モデル・出力範囲、署名状態、起動手順、既知の制約、ライセンス、SHA-256を添える。
      ランタイム入りアプリとGitHub自動生成のSource code ZIPを明確に区別する。Store配布・自動更新・専用installerを初版に要求しない。
- [ ] 新版で旧schema 1を開けることを確認する。既存プロジェクトの無断上書きを避け、旧版の配布物・復帰手順を残す。
- [ ] 公開承認後に、そのOSの受入済みバイナリを公開する。公開URLから再取得したhashと、起動・保存・復元を最終確認する。
      **完了条件：Pythonの導入なしで、配布ページから取得した同じファイルの一巡を再現できる。**

### X1：Tkinterで固定SD RSM／PCMの近似区間まで提供する後続版

初回JMLE版とは別の到達点。MMLを必要とする利用者に対しては、X1完了まではWeb／ローカルWebの対応範囲を案内する。

- [x] 固定SDの指定区間の共通計算、数値受入、既存保存点の照合は実装済み。有限標本被覆率は未認定。
- [ ] RSM／PCM・固定SDの設定とworker呼出しを共通APIへ接続する。
- [ ] MML結果・区間記録・出典を保持するプロジェクト仕様を追加し、旧schema 1の読込を維持する。
      数値ソースや環境の違う区間記録を自動承認しない。JSON区間記録だけを完全なプロジェクトと呼ばない。
- [ ] 共通資格を使うネイティブ表・区間図・保存報告を接続し、通常例・保留例・設定変更・復元を照合する。
- [ ] 対応内の大きめの入力で応答性・中止・資源を確認し、対象OSの配布物でT1〜T6の変更影響を受け入れる。
      全Web画面、自由SD、新しいG3研究、普遍的な被覆率認定を同じ版へ追加しない。

### 進め方を固定する

1. **W1はローカルで完了、次はW2→W3。** 待ち時間にT1・T5の既存資源確認は進められるが、Webへ届ける候補を追加機能で作り直さない。
2. 候補固定後の割込みは、数値・資格の誤表示、データ破損／漏えい、起動不能、必須操作の失敗など、
   この版の受入を壊す再現可能な問題に限る。修正後は無効になった証拠を更新し、該当しない改善は公開後へ置く。
3. 新しい必須作業は「この版のどの利用者の仕事が壊れるか」「既存証拠で不足する点」「完了条件への影響」を示して採否を決める。
   他ソフトの機能数、全保留の解除、文書の網羅性、テスト件数の増加を公開条件にしない。
4. 各工程の終了時に更新するのは、このチェックと既存の証拠リンク・残る阻害要因。
   新しい管理資料、同じ数値研究、細切れのパッケージビルドを追加しない。必須CIは候補に対して実行する。
5. 公開判断は**受入済みの機能・OS・実物**で行う。日付、完了率、成功テスト総数だけでは判断しない。
   人による理解度評価は公開後だが、初心者が使えるとの実証はその結果まで主張しない。

配布条件の参照（2026-10-28確認）：[Appleの起動時の説明](https://support.apple.com/en-gb/102445)、
[Developer ID](https://developer.apple.com/developer-id/)、
[Windowsの署名とSmartScreen](https://learn.microsoft.com/en-us/windows/apps/package-and-deploy/smartscreen-reputation)、
[GitHub runnerのOS](https://github.com/actions/runner-images)、
[PyInstallerの対象OS](https://pyinstaller.org/en/stable/operating-mode.html)。
署名はSmartScreen警告の即時消失を保証せず、GitHubから配布するだけでもOSの起動確認は省けない。

## 詳細な確認対象・既存証拠の参照

Actual-user evaluation is scheduled **after release** (2026-10-28 decision).
The user will arrange human evaluators. Recruitment and completion of that
review are not prerequisites for this Web release. Pre-release checks still
cover functioning journeys, output consistency and current inference holds;
automated or developer checks do not establish novice comprehension.

## Scope

- [ ] Confirm this repository is being used as an independent repo, not accidentally committed as ordinary files inside a private parent workspace.
- [ ] If the parent repo should reference this project, decide between a Git submodule and a parent-level ignore rule.
- [ ] Confirm the release label remains `standalone Python beta` unless the validation scope has changed.
- [ ] State whether this release retains inference holds or qualifies specific
      inferential outputs. A held beta update is not completion of a release
      intended to support those research conclusions. For each promoted output,
      require evidence for its actual estimator, target, uncertainty method,
      supported conditions and failure handling; UI/CI passes do not substitute.
      Review the [current causes and evidence](docs/mml_quadrature_cross_engine_review.md#inference-release-causes-20261028).
      Keep protections for extreme scores, unidentified designs and failed fits.
- [ ] Distinguish a numerically checked, explicitly conditional/asymptotic
      approximation from a finite-sample coverage or decision-error claim.
      Use the [target-specific SE contract](docs/mfrm_reporting_requirements.md#se-reference-contract-20261028);
      do not require every unfinished method's simulation for a bounded release,
      or excuse a known calculation failure by labelling it an approximation.
- [ ] Confirm `ROADMAP.md` remains the active product boundary and that no core
      feature invokes or requires an external estimation engine.
- [x] Confirm the repository-external Simulation engine loader and dead dispatch
      have been removed from `streamlit_app.py` and are covered by the
      entrypoint boundary test.

## Privacy

- [ ] Confirm no example or generated files contain confidential rating data.
- [ ] Confirm generated outputs under `validation/generated/` and retained source
      snapshots are excluded from the commit and deployment bundle; preserve the
      local research evidence.
- [ ] Confirm `.streamlit/secrets.toml` is not tracked.
- [ ] Confirm the in-app and README privacy warnings are still visible.
- [ ] Confirm README screenshots use only built-in, synthetic, or fully de-identified data.

## Deployment

- [ ] Review `DEPLOYMENT.md` before hosted deployment.
- [ ] Confirm hosted demos use synthetic, built-in, or fully de-identified data.
- [ ] Record the Streamlit Community Cloud Python version selected in Advanced settings.
- [ ] If a contextual Help adapter is marked `BROWSER_ACCEPTED`, confirm its
      stable evidence reference resolves to a reviewed record that passed
      `docs/help_browser_acceptance.md` against this exact deployed commit.

## Representative journeys

- [ ] Check all-minimum/all-maximum person scores: common report tables must
      blank their point, fit and uncertainty cells, retain raw audit values,
      and withhold person reliability while preserving valid other facets.
      Reopened projects and history snapshots must receive the same current
      qualification. Technical figure locations must identify this limitation.
- [ ] Opt into restorable history, navigate away and change inputs, then run.
      Confirm the preference survives, history remains reachable when outputs
      are stale, and restoring does not add a new fit or change the saved result.
      Download an actual report from the restored analysis and inspect its holds.
- [ ] Check G/D on supported complete balanced data and missing/repeated-cell
      counterexamples. Keep raw and zero-clipped variance estimates distinct,
      retain the observed design, and verify the dedicated CSV includes input
      identity, method, conditions and interval holds. Saving must not rerun the
      report or reset its selected tab. Browser download completion is separate
      from CSV generation tests; keep policy-blocked transfers explicitly open.
- [ ] Check observed design/category tables and manuscript drafts against retained
      rows, including weights, unused categories and stale saved summaries. Keep
      rating counts distinct from independent persons and likelihood-weight sums.
      Verify shared-person overlap counts without presenting them as a precision
      guarantee. Failed or missing optimizer termination must withhold publication
      SE/CI, reliability and inferential prose in all document formats.
- [ ] Check failed, absent and non-Boolean optimizer success in fresh and reopened
      analyses, including stale Ready/Converged tables. Technical CSVs, live
      intervals, scoring prose and document outputs must share the current hold;
      retained raw numbers must not carry affirmative screening instructions.
      Verify healthy explicit success preserves its existing conditional scope.
- [ ] Keep the separate robustness study's planned and eligible denominators,
      paired replication seeds and generating contrasts explicit. Continuous
      integration of the working likelihood does not validate its assumptions.
      Decisions use unrounded values and a justified, prespecified threshold.
- [ ] Retain the medium writing-study protocol, numerical source snapshot, all
      failures/timeouts and separate simulation denominators. Treat finite-rule
      stationarity, continuous-reference checks and Monte Carlo precision separately;
      neither a speed improvement nor this exploratory cohort releases MML inference.
- [ ] Check the simplified layout in both languages and at narrow widths:
      control labels and selections remain understandable, clean audits fold,
      problematic rows expose their reasons, and report options stay beside
      their save action. Check privacy and inference notices after navigation;
      removing repeated instructions does not remove qualification warnings.
- [ ] Check live and reopened forest plots against current inference holds;
      do not turn interval overlap into a contrast test. Verify reliability
      uses complete estimate/SE pairs, records calculation provenance, and
      withholds invalid or mismatched saved values in publications and exports.
      Keep residual Q–Q guidance appropriate for discrete responses. The bounded
      simulation does not release MML inference or replace browser/human review.
- [ ] Check current and restored step/measure/FACETS tables, EB figures, APA
      prose and scoring-quality outputs against the same inference hold. Confirm
      zero plug-in EB variance does not produce a zero-width confidence claim,
      and that an EAP hold label does not select a different reliability formula.
      Recheck document caches; UI tests do not qualify interval coverage.
- [ ] Try the two-school/common-recordings lesson in both languages. Verify
      the scenario's column defaults, shared IDs, reference exclusions and
      no-refit language switching. Confirm that readers distinguish joint
      calibration from fixed anchoring and constraint agreement from empirical
      reference stability; browser and human acceptance remain separate checks.
- [ ] Check that basic and detailed Word/PDF/HTML reports end with English Q&A
      and the cited references, and that the same guidance is readable in Japanese
      in the app. Inspect pagination, citation deduplication and cache invalidation;
      confirm that the explanations preserve current interpretation holds.
- [ ] Check basic/detailed report switching in both languages: explanatory
      Notes and holds remain visible, cached downloads match the selected scope,
      and table numbers, titles, borders and pagination survive Word/PDF export.
- [ ] Inspect nearby definitions and exported Notes for Strata and underfit:
      they must not equate Strata with actual groups or underfit with an
      exclusion decision. Actual-reader comprehension is evaluated after release.
- [ ] Inspect detailed methods/writing appendices: emphasis and list nesting
      survive export, headings stay with their text, and a calibration hold
      also restricts reliability, targeting and other dependent writing guidance.
      Retained technical values must not be labeled as ready for interpretation.
- [ ] Compare learning/reporting claims with the [current milestone checklist](ROADMAP.md#learning-reporting-milestones).
      Distinguish implemented outputs, accepted journeys and planned examples;
      documentation milestones do not qualify new SE/CI methods for inference.
- [ ] Review the [current UX1 task matrix](docs/ux_adversarial_audit.md#next-acceptance-slice--2026-10-14-planned)
      against the candidate: ordinary and held JMLE, fixed/free-SD MML, saved
      results and input recovery, in both languages.
- [ ] Record output consistency, presentation-only state preservation, browser
      operation and actual-user understanding separately. Keep unrun cases
      explicitly open; test counts do not complete the human review.
- [ ] Verify that report follow-up buttons retain the current analysis and its
      holds. Keep this functional check separate from post-release human review.
- [ ] Confirm that publication blanks and technical CSV values follow the same
      current hold, uncertainty basis and analysis identity. Inspect affected
      rendered tables, not only whether document generation succeeds.
- [ ] Record cold/repeated timing and memory for affected journeys when their
      computation or layout changes, following the [J0 comparison conditions](ROADMAP.md#performance-priority-20261014).
      Compare the same requested work and numerical accuracy; distinguish
      first use, reuse and changed-input recomputation. Retain failures and
      timeouts, and check full-bundle time when deferring individual formats.

## Statistical Validation

- [ ] Satisfy the candidate's required `make verify`/CI checks or documented
      equivalent evidence. Reuse unchanged source/environment results; do not
      rerun retained numerical studies solely because release documents changed.
- [ ] Confirm `tests/test_standalone_core_boundary.py` passes.
- [ ] Confirm every new unavailable/held computation exposes a stable ReasonCode.
- [ ] Confirm the advanced Facet Equivalence panel shows the current hold in
      English and Japanese for JMLE and fixed/free-SD MML, including reopened
      results with stale Ready labels, zero SEs or missing SEs. No TOST/BF/ROPE
      conclusions or interval-plot downloads may bypass that hold; verify the
      matching `facet.equivalence` qualification in exports.
- [ ] Confirm saved decisions reproduce from their exact AnalysisID, evidence
      fingerprints, prespecified sensitivity plan, and source records.
- [ ] Confirm native deterministic fixtures cover the supported model and design conditions changed in this release.
- [ ] Confirm the README preview image still matches the current public-beta UI after material layout changes.

The release-check payload, default Downloads archives, deterministic demo
archives, CLI, Make targets, built-in self-test registry, and normal pytest
selection are Python-native. Dormant compatibility tests carry the
`legacy_compat` marker and are excluded from `make apptest`, `make verify`, and
GitHub CI. They must be explicitly retained, deprecated, or retired before the
product-boundary roadmap item is complete; they are not release evidence.

## Commands

```bash
python -m pip install -r requirements-dev.txt
python -m py_compile streamlit_app.py
python streamlit_app.py --doctor
python streamlit_app.py --release-check
python streamlit_app.py --self-test
make apptest
python streamlit_app.py --benchmark-quick --benchmark-csv validation/generated/benchmark_smoke.csv
make clean
```

## GitHub Release

- [ ] Confirm the version, release notes, license notice, and documented feature
      scope match the candidate commit.
- [ ] Record successful required CI checks for that commit on Python 3.11 and
      3.12; earlier checked boxes are not evidence for a new candidate.
- [ ] Use the agreed publication workflow for this release, preserving repository
      history. A past history-replacement approval does not authorize a new one.
- [ ] Create the tag and prerelease only after the candidate checks pass and the
      current release is authorized.
- [ ] Record the published commit and verify the deployed app against it.

## Historical initial publication record

The completed items below describe the v0.1.x publication. They are retained
for provenance, not as a procedure or authorization for the next release.
In particular, the history replacement was specific to that publication.

- [x] Choose a license before public distribution.
- [ ] Confirm `LICENSE_NOTICE.md` still states commercial use is permitted and the software is provided as-is without warranty.
- [x] Add a remote with `git remote add origin <repo-url>`.
- [x] Because `origin/main` already has a separate history, choose a publication path: branch/PR, new repository, or explicit replacement.
- [x] Explicit replacement of `origin/main` approved by the project owner.
- [x] Push with `git push --force-with-lease origin main`.
- [x] Confirm GitHub Actions pass on Python 3.11 and 3.12.
- [x] Add a tag only after CI passes.
- [x] Create GitHub prerelease `v0.1.0-beta`.
- [x] Prepare `v0.1.1-beta` notes after cross-package validation protocol changes.
- [x] Prepare `v0.1.2-beta` notes after publication figure export and Simulation validation-template changes.

## Post-release human evaluation (not a release gate)

- [ ] The user arranges human evaluators using the existing
      [reporting comprehension task sheet](docs/ux_adversarial_audit.md#reporting-comprehension-protocol-20261016).
      Five novices and two experienced users remain a recruitment target, not
      confirmed participants. Record the evaluated release, assistance, initial
      errors, completion time and retests. Use the six comprehension tasks,
      including Strata versus actual groups and underfit versus exclusion.
- [ ] Prioritize observed obstacles and misinterpretations for a subsequent
      update. Keep human findings separate from developer/browser checks;
      do not claim novice usability has been established before this review.

## Retained research evidence

`retained_evidence` tests replay local, hash-bound study outputs. They remain
outside `make apptest` and CI because those output bundles are deliberately
not distributed. Restore the registered evidence before running
`python -m pytest -m retained_evidence tests`. Native unit tests generate
synthetic inputs in temporary directories and continue to run in CI.
A passing release gate does not certify retained research findings.

The seven timestamp-seal controller checks also require the retained trust root,
receipt and local verifier. They use the `retained_evidence` selection; their
separate environment/results must not be counted as a clean public-checkout pass.
