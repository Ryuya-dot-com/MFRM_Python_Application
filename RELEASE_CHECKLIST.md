# Release Checklist

Use this checklist before tagging or pushing a public beta release.

<a id="release-milestones-20261028"></a>

## 現行の公開マイルストーン — 2026-10-29

**品質を満たしてから公開する。公開日は固定しない。GitHubを公開案内の正本とし、版ごとのDropbox本体リンクも用いる。**
この節を今回の実行チェックリストの正本とする。下の英語の詳細項目は確認対象・証拠の参照一覧であり、
全項目を毎回ゼロから実行する指示ではない。チェックは成果物・対象ソース・環境・実施範囲の記録を根拠に付ける。
未完了の必須項目を「β」「近似」「未検証」の注記だけで完了にしない。

<a id="milestone-m0"></a>

### M0：公開する範囲を固定する — 計画として完了

- [x] Web更新は、既存の分析・報告の保護と導線改善、固定SD RSM／PCMの指定した水準・差の近似区間を対象にする。
      この区間は数値受入v1を通った対象に限る。有限標本被覆率、自由SD、個人得点、同等性等の資格を追加しない。
- [x] 最初のTkinter配布は、**RSM／PCM × JMLE／固定SD MML、バイアス、指定区間・共同交互作用、G/D** を一つのアプリに含める。
      CSV入力、列・採点基準、別プロセスの分析と中止、表・図、`.mfrm`保存／再読込、HTML・CSV・PNG・SVG、日英切替・論文用図が必須。
      Webと同じ計算・出力資格を使い、失敗・境界・識別不成立を成功にしない。
- [x] MML等を後続X1へ分ける旧案は採用せず、初回T0へ統合する。機能を削った候補の合格を初回全機能の完成にしない。
      Word／PDF、全Web画面、複数プロジェクトの一括処理、計算途中の再開、自動更新は初回配布の必須条件にしない。
- [x] 対象はmacOS Apple SiliconとWindows 11 x64。OSごとに受入して公開し、両OSの同時公開は条件にしない。
      最低対応OS版・配布可能なOSは実証後に表示する。GitHubは配布先でありOS受入の代替ではない。
- [x] Developer ID署名・公証を当面使わないGitHub配布を前提にする。macOSはアドホック署名と
      Developer ID署名を区別し、未公証と初回起動上の制約を明示する。配布後の通常の起動経路をT4で確認する。
      Windowsも署名の有無を表示する。署名・公証の追加は別の配布改善であり、購入・登録をこの計画から自動実行しない。

<a id="milestone-w1"></a>

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

<a id="milestone-w2"></a>

### W2：その候補の操作・出力を受け入れる — ローカル受入完了、2026-10-28

- [x] 通常例と保留例の入力→推定→日英の結果・図表→履歴復帰を確認した。
      新しい区間は合成180評定のRSM／PCMを実ブラウザーで操作。後から計算した区間が履歴から消える不具合を修正し、
      PCMの実復帰と、別runへの取り違え・非保存設定・再計算禁止のテストを通過した。
- [x] RSM／PCMの提供可、PCM求積5点での停止誤差・求積感度による保留、設定変更時の保護を確認した。
      PCMへRSMのJSONを実アップロードして拒否を確認。一般表・個人得点等の保留を維持し、表示時の再推定を禁止して検査した。
- [x] RSMの区間JSONを実保存・再アップロードし、AppTestから得た同じ分析の記録と一致した。
- [x] PCMの区間CSV／JSON／HTMLと通常Word／PDF／HTMLの6ファイルを、同じ候補の通常ボタンから実保存した。
      保存先の実ファイルとJSONの数値・資格・出典を照合し、PCM JSONを実アップロードして読み直した。
      別途、同じ分析の保持結果へ最適化・尤度評価を禁止してJSONを復元し、再推定を必要としないことを確認した。
- [x] RSM／PCMの指定区間を含む基本Word／PDF／HTMLを実生成し、数値・対象・注意書き・改ページを確認した。
      Word／PDFは各9ページ。共通ページの描画一致を使い、差のあるページも目視した。
      対象ファセット・単位・資格列の説明と折返しを修正した。生成・描画の成功をブラウザー実保存の代用にしない。
      追補で実保存したPCMのWord／PDF各8ページも、既に目視済みの同一描画ページを再利用して確認した。
      前者は差を含む4対象のAPI出力例、後者は3水準と実画面の分析情報を含む報告であり、報告全体が同一という意味ではない。
- [x] データ保護、終了失敗、境界、保存結果の取り違え、基本／詳細・日英切替の既存受入を候補へ対応させた。
      新差分は影響範囲35テストと実操作で確認。数値コア・資格条件はW1から同一であり、既存数値研究を再実行していない。
      人による理解度評価は公開後。数値的な誤表示や操作不能を人による評価待ちで残さない。
      **完了条件：対象経路の既知の重大不具合がなく、各出力の利用条件が同じ候補で一致する。**

[修正と受入範囲](docs/development_integration.md#web-w2-acceptance-20261028)を参照。
W2はこのローカル候補で完了。初回の「組織でブロックされました」は、同じ操作ツール管理下のタブでの手動クリックでも発生し、
ユーザーの手動確認では`DownloadRestrictions`は「なし」だった。その後、操作ツールの正式なダウンロード待受APIと
通常の保存ボタンを用いて6ファイルを取得した。製品コード・ブラウザー設定・保護設定は変更していない。
最初の遮断の内部原因は未特定であり、組織管理・Chrome全般・アプリの生成不具合と断定しない。
操作ツールが拒否した管理URLへの別経路アクセスや制限回避も行っていない。公開先での実保存はW3で別に確認する。

<a id="milestone-w3"></a>

### W3：CIから公開先まで閉じる — 必要CI通過、公開反映・公開先受入は未完了

追補候補`f5139bb7ec02fbfa2f1079ef6e381acff468d683`の[draft PR #5](https://github.com/Ryuya-dot-com/MFRM_Python_Application/pull/5)で、
[Linux／Python 3.11・3.12 CI](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/40687727467)の全工程が通過した。
各環境の通常テストは4,010件通過・10件skip・91件対象外。UX契約73件（通常テストと重複）、benchmark、demo報告、clean-checkも通過。
10件skipは外部R／immer不在によるもので、以前のローカル成功記録と区別して対応を残す。新たな推論資格を与えない。
CIが使うPRの合成mergeと候補のtreeは一致する。アプリ本体・翻訳・依存設定はW2受入時と同一。

初回候補`350b95323`の[CI](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/40679224500)は、
3.11で3,990件通過・8失敗・10準備エラー、3.12で3,988件通過・10失敗・10準備エラー。
両方とも10件skip・91件対象外。後段のbenchmark／demo報告／clean-checkは未実行。
固定一時パス、保存比較の丸め尺度、過去の最適化経路との過剰に厳しい一致、操作テストの既定3秒を修正した。
アプリ本体・数値基準・保存された生の参照値は変更していない。ローカルの影響範囲97項目に成功証拠を対応させ、
検証コード・テストの追補を含む候補で上記の必須CIを通過した。初回の失敗ログは保持する。
main・tag・GitHub release・稼働Webは変更していない。

**公開前追補：** 同じ候補の通常画面で、欠測セルを説明するG/D CSVとPathway図PNGを実保存した。
CSVの条件・区間保留、PNGの対象・軸・分析ID、保存前後の画面保持を確認し、古い「実保存未確認」を閉じた。
公開日・引用の確定、公開先の代表操作、失敗時の復帰は[公開・復帰手順](DEPLOYMENT.md#publication-and-hosted-acceptance)へ具体化した。
公開日未設定・未公開表示は維持する。[追補証拠](docs/development_integration.md#web-prepublication-20261028)を参照。
これらはローカルの文書追補であり、新しいCI・main反映・公開先受入ではない。下の未完了チェックは維持する。


**公開実行の承認（2026-10-28）：** ユーザーはmain反映・公開先受入を継続し、合格時のtag／GitHub prerelease公開を許可した。
公開日・引用を2026-10-28（UTCの公開日）へ同期する。数値・アプリ・依存設定は受入済み候補と同一。
引用の影響範囲をローカル確認し、文書・引用のPR追補だけはCIの重複起動を抑える。
マージcommitには抑止を付けず、最終mainのLinux 3.11／3.12 CIを一組実施する。
tag／prereleaseはそのCIと公開先受入の両方が通った後に作成する。公開結果はリリース記録へ残し、
結果記載だけの追加pushで同じ全体検査を再起動しない。下の未完了は実際に通った段階でローカル正本へ反映する。

- [x] 同じ候補commitに必要なLinux／Python 3.11・3.12 CIを一組対応させた。
      現行workflowはintegrationブランチへのpushだけでは自動起動しない。PR／手動実行を選び、重複起動しない。
      必須CIの全体確認は省かず、通過済みのローカル数値研究を文書変更だけで再計算しない。
- [ ] 現在の公開への承認がそろった段階で、main反映・tag・GitHub prerelease・Web更新を実施する。
      今回の計画更新や過去のforce-push承認を、新しい公開・履歴置換の承認に読み替えない。
      公開時にCITATIONの実際の公開日とREADME／変更履歴の候補表示をそろえ、引用の影響範囲を確認する。
- [ ] 公開先の版・commitと、合成データの代表操作・実保存を確認する。不一致なら旧公開へ戻す経路を使う。
      **完了条件：GitHubの公開物と稼働Webが採用した候補に一致し、利用者の入口から保存まで通る。**

### T1〜T6：Tkinter初回配布に残る6つの作業群

共通計算・schema 2（旧schema 1読込）・worker・各機能の表図を統合済み。既存のWindows凍結5分析・再推定なし復元の証拠を保持し、今回のmain統合による変更影響を確認する。
[試作記録](docs/desktop_prototype.md)はこの証拠であり、一般配布済みの証拠ではない。
T1〜T4はmacOS先行、T5はWindows、T6は受入したOSごとに実施する。

<a id="milestone-t1"></a>

**T1：現在のソースで配布候補を閉じる（未完了、W1の共通ソース確定を利用）**

事前点検：既存の`literature_ux_20261027`試作は、同梱manifestと現在のソースで9ファイルが異なる。
`Info.plist`の`CFBundleShortVersionString`も`0.0.0`であり、画面上のAPP_VERSIONだけでは配布物の版一致を満たさない。
同梱ライセンス一覧ではNumPy／pandas／SciPyと図用フォントを確認したが、Python／Tcl・Tkや他の同梱依存物までの網羅確認は残る。
同梱Mach-O 274ファイルのロードコマンドを確認すると、NumPy／SciPyにはmacOS 14.0を要求するバイナリが含まれる。
起動用実行ファイル単体の11.0を配布物全体の下限にしてはいけない。14.0での実動作を受け入れたという意味ではない。
記録は`validation/generated/web_candidate_20261028/w3_0_2_19/t1_macos_binary_inventory.json`に保持する。
この点検は既存配布物の読取りであり、新しいビルド・OS受入の成功ではない。

追補の同梱ファイル点検では、Python／Tcl・Tk、Matplotlib本体、PyArrow等のライセンス文書をファイル名一覧からは確認できなかった。
一覧だけで法的な適否は断定しない。次回ビルドの前に実際の同梱依存物と配布元のLICENSE／NOTICEを対応させ、不足を収める。
9ファイルの具体名と同梱版の記録は`prepublication/t1_existing_bundle_review.json`（上記W3証拠ディレクトリ内）。
最低OSは既存バイナリの要件14.0以上と実動作受入を区別し、未検証の14.0対応を表示しない。

- [ ] 10月28日の個人境界・出力保護を含め、旧パッケージと異なるソースを取り込む。
      版・OS／CPU・依存版・同梱資源をmanifestと実物で照合する。試作の再梱包を最新修正の反映と取り違えない。
- [ ] 最低対応OSを選んで検証範囲を固定する。アプリ・依存ライブラリの必要なライセンス表示と
      見本・翻訳・フォント等を同梱する。初版に不要な新依存を入れない。

<a id="milestone-t2"></a>

**T2：実GUIの一巡を受け入れる（未完了、T1に依存）**

- [ ] 実ファイルダイアログでCSV入力→列・採点基準→分析→各分析の表・図→HTML／CSV／PNG／SVG／プロジェクト保存を通す。
      日本語パス、ID保持、入力不備、通常のキーボード操作と表示範囲を確認する。
- [ ] 中止→再実行、終了失敗、極端得点・保留を確認する。未保存結果の扱い、入力変更と前回結果の区別、
      読み取り専用入力・書込不能先での説明と既存ファイル保全を確認する。表示・読込で再推定しない。

<a id="milestone-t3"></a>

**T3：開発環境への依存を取り除く（未完了、T1に依存）**

- [ ] Python・開発用パッケージを利用しないクリーンな対象環境で、アプリの初回起動、worker、保存・再起動後の復元を確認する。
      開発端末でPATHを変えるだけの試験と、Python未導入環境の受入を区別する。VM・検証協力者等も使える。
- [ ] ネットワークなしで各分析・図・HTML・プロジェクトが利用できることを確認する。
      Web Helpだけはネット接続が必要と説明し、オフラインを全機能の保証へ広げない。

<a id="milestone-t4"></a>

**T4：macOSのダウンロード後の起動を受け入れる（未完了、T1〜T3に依存）**

- [ ] ZIP等の実際に配るアーカイブから、展開・配置・起動・分析・保存を通す。
      公開前は許可されたdraft release／検証用artifact／Dropboxリンクで取得操作を確認し、一般公開後のURLはT6で再確認する。
      採用する配布元からの取得で付く隔離属性を含め、初回のOS表示と必要な操作を確認する。公開前後で対象バイナリのhashを対応させる。
- [ ] Developer IDなし・未公証であること、対応OS、発行元・ファイル照合の方法を配布ページへ明記する。
      OS全体の保護解除や隔離属性の一括削除を通常の導入手順にしない。
      OSが許す個別の確認操作で起動できない環境を、対応済みと表示しない。
      **macOS配布判定：T1〜T4が通り、T6の公開物を同定できること。署名サービスの取得待ちを自動で追加しない。**

<a id="milestone-t5"></a>

**T5：Windows版を作り、Windows 11で受け入れる（未完了）**

- [x] GitHub ActionsのWindows runnerでビルドし、同梱ランタイムの起動・worker・保存／復元・HTML等を確認する。
      [49eec16の検査](https://github.com/Ryuya-dot-com/MFRM_Python_Application/actions/runs/40867926816)でRSM／PCM×JMLE／固定SD MMLと独立G/Dの5分析、指定区間・共同交互作用と全保存形式、別プロセスの再推定なし復元を確認。main統合後の新実物は今回のCIで照合する。
- [ ] Windows 11 x64の検証環境を確保し、T2・T3相当と実ダウンロード後の起動表示・日本語パスを確認する。
      検証協力者／VMを含み、開発者の所有は必須にしない。Windows Serverのrunner合格をWindows 11のGUI合格にしない。
      Windows 11実機は利用可能。c55d082の基本一巡は利用者確認済みであり、同じ確認を重複要求しない。初回全機能の実物の受入完了までは検証用成果物とし、Windows 11対応の一般配布は完了にしない。

<a id="milestone-t6"></a>

**T6：受入したOSのファイルをGitHub Releasesへ届ける（未完了、macOSはT4、WindowsはT5に依存）**

- [ ] 対象OS／CPU／最低版、モデル・出力範囲、署名状態、起動手順、既知の制約、ライセンス、SHA-256を添える。
      ランタイム入りアプリとGitHub自動生成のSource code ZIPを明確に区別する。Store配布・自動更新・専用installerを初版に要求しない。
- [ ] 新版で旧schema 1を開けることを確認する。既存プロジェクトの無断上書きを避け、旧版の配布物・復帰手順を残す。
- [ ] 公開承認後に、そのOSの受入済みバイナリを公開する。公開URLから再取得したhashと、起動・保存・復元を最終確認する。
      **完了条件：Pythonの導入なしで、配布ページから取得した同じファイルの一巡を再現できる。**

<a id="milestone-x1"></a>

<a id="milestone-t0"></a>

### T0：初回全機能を共通計算・保存・表図へ接続する

旧X1を初回へ移管する。上の旧アンカーは参照先として残すが、後続版への延期を意味しない。

- [x] RSM／PCM・JMLE／固定SD MMLの設定とworkerを共通APIへ接続する。
- [x] バイアス、指定区間、共同交互作用、完全均衡1～2ランダム要因のG/Dを同じプロジェクトへ保持する。
- [x] schema 2に結果・依頼設定・利用条件・出典を保持し、旧schema 1の読込を維持する。表示・保存・再読込では再推定しない。
- [x] ネイティブ表・図・HTML／CSV／PNG／SVG、英語既定・日本語切替、説明付き／論文掲載用図を接続する。
- [ ] mainの最新修正を保持した全機能候補で、T1〜T6の残る配布条件を満たす。

数値・推論条件と適用設計は維持する。固定較正バイアスの有意性判断、G/D区間、多変量G/D、自由SD、普遍的な有限標本被覆率は、この統合によって提供済みにはならない。

### 進め方を固定する

1. **初回はT0の全機能を含める。** mainの最新修正と全機能ブランチを統合し、競合・変更影響を確認して同じ候補を配布へ進める。過去の成功証拠を再利用し、限定機能版の合格を初回配布の完了にしない。
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
