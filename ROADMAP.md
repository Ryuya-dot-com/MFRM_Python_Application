# MFRM Python Application Roadmap

## 現在の実行計画 — 2026-09-15

完成条件と優先順位の最終照合: **2026-10-09**。節の起点日と既存リンクは保持する。

**最終目的は、評定研究の計画から、根拠と限界を理解した分析・論文報告までを、一つのアプリで支えること。**
数学的正確性を最優先とし、利用者が「なぜこの分析を行い、何が分かり、何がまだ言えないか」を説明できる状態を目指す。

今後の着手順と現在地は本節を正本とする。下の実装履歴、旧M1–M4の時期、各文書の
日付付き「次の作業」は当時の記録として残す。統計的な資格は
[G0–G6](docs/statistical_engine_remediation_roadmap.html#phases)と対象別の検証計画が定め、
本節の優先順位や公開実績によって解除しない。

**2026-10-09の実行状態：ユーザーの指示により正式G3と後続MMLの起動を一時停止。**
G3は28,734終端・2試行進行中の状態でプロセスを停止し、途中のメモリと台帳を保持した。
現在の優先作業は、mfrmr開発版0.2.4.9000を参照したGPCMの数理・求積・計算効率の監査である。
固定1,281点でも探索範囲内の厳しい母数位置では連続積分との差が残り、単純な適応的求積への
置換も速度改善にはならなかった。根拠と再開前の条件は
[GPCM再監査](docs/mml_population_sd_uncertainty_audit.md#gpcmの数理求積とmfrmr開発版の照合--2026-10-09)に記録する。
凍結済み計画・ソース・既存試行は保持し、自動再開しない。推定対象や推定法を変える場合は、
既存の正式結果へ混在させず、計画と実装版の対応を整理してから再開の指示を受ける。

2026-10-02のローカル数理監査では、GPCMの条件付きSEに欠けていた識別力の二乗と、
デルタ法の計算不能な偏微分のゼロ扱いを修正した。自由SD・正則化ありMMLについては、
現行の構造Hessianから通常の構造SE・区間を作る経路を保留し、未計算のHessianと
階数欠損を区別した。尤度・同定座標・推定器・数値最適化は変更していない。
根拠と未完了事項は[不確実性監査](docs/mml_population_sd_uncertainty_audit.md#数学的な不確実性の再点検--2026-10-02ローカル修正)に記録する。
同日のZotero文献照合では、Murakiの情報量の原式を確認し、有限混合でMMLの観測情報と
EM M-stepの曲率を独立に照合した。`PosteriorMAPGrid` が求積質量を最大化していた問題を
事後密度の最大化へ修正し、EAPの事後SD・正規近似区間が校正値を固定していることを明記した。
書誌・要旨と本文の読取範囲は[文献を用いた追跡監査](docs/mml_population_sd_uncertainty_audit.md#zoteroの文献を用いた追跡監査--2026-10-02)に分けて記録する。
続くMyford・Wolfe文献照合では、バイアス対比較の尺度方向と共有主効果の相殺を修正した。
差のSEは校正値を固定した近似に限定し、根拠の未確認なWelch自由度・p値の生成を止めた。
詳細は[バイアス対比較の監査](docs/mml_population_sd_uncertainty_audit.md#バイアス対比較の符号と共有項--2026-10-02)に記録した。
これらはローカル修正であり、共同の不確実性・参照分布・誤検出率等の検証を要するG2・G3は保留を継続する。

### 出発点と、まだ残っている課題

| 領域 | 確認済みの到達点 | 次に埋める差 |
|---|---|---|
| 公開・運用 | 0.2.18-beta / `b483d97`を公開。候補treeのPython 3.11/3.12 CI、公開画面の日本語・合成例の分析・ホームからの復帰を確認 | 実利用者による理解の確認、正式なキーボード・フォーカス受入は未完了 |
| 数値計算 | 有限求積の尤度・勾配、二段階最適化、独立R・連続積分との限定した照合。固定SD PCMの現在の推定点を調べる画面を公開 | 停止・停留性・積分精度・推論の可否を、現在の分析に結び付けて表示。自由SDのSE・CI・被覆率の資格は未完了 |
| 他ソフトとの比較 | 現在の分析に対応した比較用出力と、限定PCMのオフラインTAM/ConQuestランナーを公開 | 外部結果の検証付き読込は未提供。重みの意味や対応モデルが違う比較、一般的同等性は未認定 |
| 結果・報告 | 残差PCA・バイアスへの導線、目的別の出力、英語APA原稿、MIT・Ryuya Komuro名義の引用を実装 | 初見の利用者が診断の限界と次の行動を理解し、原稿の主張を根拠まで追えるかを確認 |
| 研究計画 | 割付・評定負担・費用・接続性のPython計算基盤を統合 | 事前設計画面の接続は未着手。SE改善・必要人数の推奨を支える検証は含まない |
| 高速化・ベイズ拡張 | JAXのJ0–J4、Stan/NumPyroのB0–B4を計画済み | JAXの実装・性能比較、ベイズコードの公開配布は未着手 |

公開根拠と統合範囲は[統合台帳](docs/development_integration.md)、数値結果は
[比較検証記録](docs/mml_cross_language_development_check.md)、画面の確認範囲は
[UX点検記録](docs/ux_adversarial_audit.md)へたどる。テスト通過数は実装確認の記録であり、
被覆率や使いやすさの証拠数として合算しない。

### 大局：成果を積み上げる順序

| 段階 | 利用者にとっての成果 | 先に満たす条件・進め方 |
|---|---|---|
| 直近：数学的な受入条件を固定する | どのモデル・条件・出力まで根拠を揃えるかが明確になる | N2の対象・出力対応表と検証仕様を先に確定。共同不確実性を主課題とし、数値最適化・高速化とは別に判定する |
| 直近：現在の結果を正しく読める | 求積を変えた影響、未確認の事項、次の操作が分かる。自然な日本語でPCA・バイアス・原稿へ進める | 下記N1とUX1。通常画面を増築せず、既存の計算詳細と結果導線を使う |
| 次：数値と不確実性の根拠を固める | 公開してよい推定値・SE・区間の範囲を、条件別に説明できる | N2と既存G2–G4。自由SDの保留解除は、この根拠が揃ってから個別に判断 |
| 次：計画から報告をつなぐ | データ収集前に設計の弱点を調べ、分析後は判断の理由を論文・再現用出力へ残せる | D1とUX1。構造的な設計確認はN2の完了を待たずに進められる。精度・必要人数の提案は別途検証 |
| 条件が揃った拡張 | 同じ分析の待ち時間を短縮し、必要な人は外部ソフトやベイズ推定へ進める | E1、J0–J4、B0–B4。J0の測定とB0の仕様整理は可能。JAX採用は実測効果、ベイズ実装は既存G4等を前提とする |

暦上の期限だけで次段階へ移らない。2026-10-04以降の数理上の最初の区切りは
「N2の対象・出力対応表と検証仕様」とする。部分実装済みのN1は、N2に必要な数値の確認と
現在の分析への説明を支える範囲で仕上げる。UX1の代表操作の点検は独立して進める。
常時着手するのは数理上の主課題一つと、それを妨げないUX・文書課題一つを基本とする。
UX1は既存UX計画のU0–U6を代表操作で点検する作業であり、画面構成を置き換える別計画ではない。
小さな日本語・操作上の修正を自由SD研究の終了まで待たせず、新推定法を追加して主課題を先送りしない。

### 局所：次に実施する作業と完了条件

以下のIDは作業を参照するためのもの。統計資格のG番号や、実装の状態フラグを増やすものではない。

| 作業・優先度 | 最初の成果と既存実装の再利用 | 完了条件と依存関係 |
|---|---|---|
| **N1：求積感度と現在の分析の接続／N2を支援** | `mml_quadrature_sensitivity.py`、現在の入力保存・推定点確認を再利用し、固定点での再評価、次数ごとの再推定、再得点化を区別する対応表を作る。その後、限定した読み取り表示・保存を接続 | データ・全設定・制約・重み・座標・求積規則・実装版を照合し、元測定値から判定を再構成。古い結果、欠落、目的関数の取り違え、次数が近くても積分誤差が大きい例を拒否または未確認として保持。元の推定結果を置換せず、UIと出力の値・保留理由が一致する |
| **N2：自由SDの不確実性／最優先の研究課題** | [SD監査](docs/mml_population_sd_uncertainty_audit.md)と[v2計画](validation/MML_FREE_SD_STATIONARITY_V2_ROADMAP_20260811.md)を再利用。初期範囲、共同情報行列・プロファイルの問い、独立照合、求積誤差、失敗処理を固定する | まず単位重み・無罰則RSM/PCMの同定可能な内点条件に限定。数値検証→適切な不確実性計算→新規データでの回復・被覆率を別々に確認。境界・弱情報・欠測の扱い、全試行の分母、MCSE、反復数の根拠を結果を見る前に定める。合格対象以外と旧保存結果の保留は維持 |
| **UX1：理解できる診断と論文／直近** | 既存ホーム・サンプル・PCA・バイアス・Report & Exportを、同じ代表課題で通して読む。日本語の動的説明・エラー・ヘルプを文脈ごとに見直す | 日英、狭い画面、キーボードで実行・復帰・診断・英語原稿取得を確認。表示操作による再推定は0、結果の取り違えは0。少人数の初学者と経験者による形成的レビューで、無支援完了・迷い・解釈誤りを記録。実利用者未確認なら未完了のままにする |
| **D1：事前設計の入口／次** | `mfrm_app/simulation/planning.py`と保存済み`_working`画面を比較し、負担・費用の確認と明示操作での接続性評価を最小限接続 | 計画だけでも使え、分析済みデータを変更しない。切断、一本の橋への依存、過大な計画、評定者離脱の対応範囲を明記。構造的な改善案と精度改善を区別し、未計算のSEや必要人数を表示しない。Help・言語・ホーム復帰を確認 |
| **E1：外部比較の結果読込／次** | 公開済みの限定PCMランナーと`mml_comparison_report.py`を起点に、対応範囲内の保存結果だけを検証して読む | 固定校正のモデル対応を確認してから推定結果を比較。実行終了、校正値、実際の返却点の勾配、積分精度、得点MC変動を分離。数値許容差と実質的な同等性の幅は別に根拠を定め、結果を見て緩めない。版・入力・座標不一致や不足証拠では比較を成立させない。新モデルへの拡張より先に一往復を完了 |
| **J0→J1–J4：既存推定の高速化／測定から** | 公開版から残る改善を個別に照合。データ準備・推定・診断・描画・出力の時間とメモリを測り、支配的な計算核だけをJAX候補にする | 初回・反復と全分析時間を分け、同じ入力・停止条件・float64で比較。採用範囲の数値基準を先に固定し、推定値・勾配・利用可能なSE・判定を照合。総待ち時間が改善しない範囲はNumPy/SciPyを維持。自由SDの資格は継承しない |
| **B0→B1–B4：Stan/NumPyro配布／後続** | [拡張計画](docs/jax_numpyro_conquest_tam_plan.md)の保管コード監査とZotero参照記録を再利用し、共通RSM・同定・事前分布・出力仕様を先に確定 | 新規サンプリングは既存G4と計算範囲の条件を確認後に開始。選択エンジン一つの最小コード配布から検証し、ローカル実行・MCMC診断・結果読込を別々に判定。長時間サンプリングをStreamlitへ載せない |

**N1の最初の境界：** 既存Q感度モジュールは自由SD用で、公開中の推定点確認は固定SDの
限定PCM用である。同じ名前の診断として直結しない。固定SDではSD差を評価対象とせず、
自由SDでは`log(SD)`を含む共同勾配が必要。まず適用範囲と入力対応を確認し、現行の
読み取り専用導線を保つ。再最適化を含む候補を扱う場合も元結果と分け、計算コストと対応範囲を
確認するまで自動実行しない。候補を正式な結果へ反映する機能は、全得点・診断・共分散・
区間・原稿を一貫して再計算できるまで後続とする。

**N2の研究開始条件：** 観測済みの開発例を確認試験に再利用しない。保存入力から推定前に
問題を構成する処理、ソース・実行環境と結果の復元、計画がデータ生成前に存在した証拠を、
既存v2計画に従って整える。まず小規模な数値対照で手順を確認し、大規模な反復計算は
評価する問いとMC精度から必要量を定めてから行う。構造効果・母集団SD・個人得点の
不確実性を別の対象とし、校正値を固定した事後SDと校正を再推定する区間を混同しない。新しい結果で過去の失敗を置き換えない。

### 数値検証で必ず分ける問い

N1の最初の実装（2026-09-15、未公開）：限定した固定SD PCMで、現在の推定値を固定した
GH31/61/121の尤度・EAP・事後SD比較とJSON保存を接続した。入力・ソース・実行環境の
識別により古い確認結果を無効化し、数値差が小さくても精度の合格判定は発行しない。
画面内の確認は1,000回答行まで。詳細は[固定点比較の記録](docs/mml_cross_language_development_check.md#固定sdの同じ推定点での次数比較)。

N1の続き（同日、未公開）：保存済みの固定SD PCM 2例について、各次数・2開始点の計12経路を
オフラインで再推定し、同じ点の積分誤差と、連続積分上の再推定による変化を分けて保存した。
全経路で有限求積の最終勾配は1e-9以内だったが、16人の例のGH31では事後SDの同一点積分差が
最大約0.02071 logit残った。欠測例は連続積分参照の粗い側の誤差見積もりが事前基準を満たさず、
確認未達として保持した。停止・開始点間の一致を積分精度の保証に読み替えない。
画面への再推定接続、連続積分での再最適化、自由SD用の既存Q感度評価との接続は未完了。
詳細は[再推定と積分誤差の記録](docs/mml_cross_language_development_check.md#固定sdの再推定と積分誤差の分離)。

N1の参照誤差の点検（同日、未公開）：2例の元の点・全再推定点の計14点で、前回と同じ積分設定・
評価回数・数値を再現した。欠測例の事後SDの誤差見積もりは、全員分の総和と最悪のモーメントを
使う従来方式で最大1.414e-8、個人別では最大1.776e-10 logitだった。
未達の理由を集約方法の保守性に切り分け、同じ数値基準を個人別に適用した診断は14点で満たした。
総NLLには総和による誤差評価を維持し、過去の未達フラグや推論資格は変更していない。
詳細は[個人別誤差評価の記録](docs/mml_cross_language_development_check.md#固定sd参照の誤差集約の点検)。

N1・UX1の説明と手順（同日、未公開）：閉じた数値詳細に「読み方と次の確認」を追加し、
同一点比較・再推定・連続積分の違いと保存先を日英で案内した。既存の比較用ZIPには、
表示言語に合わせた[日本語手順](docs/fixed_sd_quadrature_check.ja.md)／[英語手順](docs/fixed_sd_quadrature_check.en.md)を同梱する。
JSONだけでは再実行できないこと、開発用Python確認には別途ソース一式が必要なことを明記した。
関連18テストと、実際のZIPから6経路を実行する入出力確認が通過した。推定関数は未変更。
N1・UX1の通常画面との接続（同日、未公開）：実際の入力・MML推定から求積確認、比較用ZIP、結果への復帰までを点検した。
比較資料の選択が消える問題、言語切替で入力画面の状態が失われる問題、翻訳前の選択肢をブラウザーが返して
詳細設定が標準設定へ戻る問題を修正した。最後の問題は実ブラウザーで発見し、翻訳後の選択値を送ることも回帰テストに追加した。
関連155テストが通過し、実ブラウザーから保存した日英ZIPの入力・推定値のハッシュも一致した。
「標準得点」という誤訳は「各ソフトの個人得点」へ修正した。計算式・推定処理・推論の資格は変更していない。
詳細は[通常画面での確認記録](docs/ux_adversarial_audit.md#通常画面から求積確認と比較資料への往復--2026-09-15ローカル実装)。
次は、優先度・信頼性の案内などに残る動的な英語文を整理し、PCA・バイアス・原稿への導線を通して読む。
初見の利用者による理解の確認と、アプリ全体のアクセシビリティ受入は引き続き未完了。

| 問い | 記録するもの | そこからは言えないこと |
|---|---|---|
| 指定した計算が終了したか | ソフト・版・停止理由・開始点・最後の点・返却点 | 最尤解に到達した、という保証 |
| 同じ有限求積目的関数で停留しているか | raw NLL、座標と尺度を明記した勾配、制約・再開始の差 | 連続積分の精度や大域最適性 |
| 積分近似は十分か | GH次数、固定格子の範囲と間隔、独立連続積分、構造値・EAP・事後SDの差 | ある次数や裾範囲が全データに十分、という保証 |
| 推定全体の結果が安定するか | 同じ点の積分誤差と、各規則で再推定した後の総変化を別記 | 近接する二次数の一致だけによる精度認定 |
| 不確実性と解釈を支えられるか | 対象別のSE・区間・回復・被覆率・失敗率と、研究目的に対応する根拠 | 停止やソフト間の一致だけによる推論・公平性の認定 |

mfrmr開発版は2026-09-15の`DESCRIPTION`で0.2.4.9000。内部ロードマップには高次数GHの
微小重みのゼロ化が未解決として残る。原因・fit/scoring両経路・積分への影響を別々に追い、
epsilonへの置換やQ増加だけで修正済みとしない。mfrmrの研究結果は対象と条件を対応づけて使い、
Pythonの被覆率の証拠へ転用しない。詳しくは[対象別の対応表](docs/jax_numpyro_conquest_tam_plan.md#mfrmr-024開発版から取り込む判断)。
TAM/ConQuestも同じ役割で使う。回答重みと受験者重み、有限事前分布の正規化、EAP/WLE、
校正の求積と得点の乱数計算を揃えない比較から同等性を主張しない。
ConQuest内部の追加再構成は、利用者の結果解釈に影響する未解決の問いがある場合に限定する。

### 統合アプリとしての完了判断と、見直し方

公開する**モデル・推定対象・用途の範囲ごと**に、次の四点を満たしたときに完成と判断する。
将来候補の全推定法を実装することは完成条件にしない。

1. **研究計画：** 想定する割付の弱点と負担を示し、構造の確認と精度の保証を区別できる。
2. **分析：** 公開する数値・SE・区間・診断に、その条件に適用できる根拠と未対応時の扱いがある。
3. **理解と報告：** 利用者が主要な限界と次の行動を説明できる。APA原稿は「問い→設計と理由→結果→答え」を追え、未入力の研究目的や考察を創作しない。各主張を元の分析・診断へたどれる。
4. **再現と運用：** 同じ入力・設定・依存環境で再実行し、記録した許容差内で再現できる。データ保護、日英操作、アクセシビリティ、公開環境での動作を確認する。

#### 完成の範囲と具体的な受入条件 — 2026-10-04

現行ソース・統合台帳・数理監査・検証計画を照合した。**G0（出力保護）・G1（同定座標）は完了、
G2以降の推論資格は未完了**。10月2日の数式修正はローカル実装であり、公開版への反映や
被覆率の確認を意味しない。この点検では公開環境の再確認や、新しい数値実験は行っていない。

「完成」は次の二つを区別する。

- **範囲を限定した研究用アプリの完成：** 対応するモデル・条件・出力を明示し、その範囲で
  計画→分析→診断→報告→再現を通せる。対応内の実用的な正例が推定・不確実性・報告まで
  到達し、対応外の例は理由付きで停止・保留される。すべてを保留している状態は完成と呼ばない。
- **現行の統計修復計画全体の完了：** G0–G6を通過する。現行計画のG6にはG5の構造的な
  非Person交互作用も前提として含まれる。限定版の完成をG5/G6の完了やbeta解除へ読み替えない。
  全計画の範囲を変更する場合は、その判断を別途記録する。

最初に範囲を閉じる対象は、既存N2どおり**単位重み・無罰則のRSM/PCM、同定可能で内点の条件**
とする。これは検証の開始範囲であり、すでに推論が適格という宣言ではない。
RSM/PCM/GPCM、JMLE/MML、母集団SD固定/自由、重み、罰則、アンカー・制約、データ設計、
点推定/SE/区間/検定を対応表にし、出力ごとに完成を判定する。たとえばRSM/PCMの成功を
GPCMのstep・識別力や、そのLRTへ継承しない。JMLEの記述用途が完成しても、
MMLと同じモデル選択の資格を与える必要はない。

| 対象 | どこまで実装・検証すれば完成か | 現時点の不足 |
|---|---|---|
| **数値の正しさ：N1・G1/G2** | 確率・尤度・勾配・同定を独立計算と照合。現在の返却点の停留性、同一点の積分差、再推定による差を区別して保存・説明できる。採用範囲と許容差を先に定め、失敗を保持する | 限定した固定SD PCMの画面とオフライン検証はある。再推定・自由SDへの接続と、対応範囲全体の確認は未完了 |
| **不確実性：N2・G2** | 自由SDと構造効果の結合を扱う共分散またはプロファイル等を実装し、独立照合と被覆率で確認。能力の条件付き事後SDと、校正不確実性を含む量を別に定義する。境界・弱情報・計算不能時の扱いを揃える | 共同情報・共分散の独立R照合は開発実装済み。自由SD・罰則ありの一般的な公開SE/CIは保留。求積精度と新規データによる被覆率の確認が必要 |
| **モデル比較：G2** | 同じ入力・尤度・母集団仕様・重み・制約での比較だけを認め、自由度と参照分布の根拠を持つ。適格範囲だけでICの解釈・推奨・LRTを公開する | Kと保護は実装済み。RSM–PCMの適格MML比較の資格は未完了。通常のGPCM LRTやJMLE自動推奨を一括解除しない |
| **バイアス診断：G3** | 尺度方向・差・SEの対象を一致させ、境界・最適化失敗・重みを記録。一次p値と多重比較の対象集合を固定し、採用する条件付き近似等の誤検出率・検出力・被覆率を確認する | 符号・共有項・単回計算・境界/失敗記録を修正済み。一次参照値を単位重み条件付きLR、一次補正をHolmへ統一。独立Rによる15ケースの計算対照は通過。校正後の誤検出率・検出力・被覆率と公開推論の資格が残る |
| **回復・性能の証拠：G4** | 水準別・真値別のbias/RMSE、SE、被覆率、失敗率を評価し、GPCMではstep・log-slope・slopeも対象にする。各公開用途を条件・独立反復数・seed・MCSE・全試行分母・再現物へ結ぶ | 水準/真値分位別bias、反復MCSE、同定後の真値、step/slopeの回復、共分散変換、全試行分母を実装。各反復でPersonを取り直す母集団回復、新規の凍結反復研究と推論資格への証拠接続は未完了 |
| **研究計画：D1** | 既存の割付・評定負担・費用・接続性計算を画面へ接続し、分析データがなくても計画を保存・再現できる。切断や脆弱な接続の理由を説明できる | 計算基盤は統合済み、画面接続は未着手。精度・必要人数の推奨は、この構造確認の完成条件には含めず、主張する際に別途検証する |
| **理解・報告・公開：UX1・G6** | 初学者と経験者が代表課題を通し、主要な限界を説明できる。日英・キーボード・狭い画面で操作でき、表・図・APA原稿・保存結果の数値と資格が一致する。clean checkout、旧結果、公開環境で再現を確認する | 画面と原稿出力は実装済み。実利用者の理解、正式なアクセシビリティ受入、出力別の根拠の結び付け、未公開修正のリリース確認が残る |

数値照合の許容差、被覆率・誤検出率等の許容範囲、反復数・MC精度、失敗処理は、
対象とする出力・設計ごとに**検証結果を見る前に**決める。観測された結果に合わせて基準を緩めない。
合格した範囲が狭ければその範囲を明示し、条件を広げる際に追加検証する。
現在の文書だけでは、未実施の研究を横断した全出力共通の数値基準は確定していない。
検証用コードの実装、独立な数理・数値照合、実際の反復研究、公開環境への反映を別の完了として記録する。

#### ここからの順序と後続拡張

1. **最初の成果物はN2の対象・出力対応表と検証仕様。** 固定SDで先に確認できるものと自由SDで
   共同推定が必要なものを分け、文献の式・同定・実装・独立対照を対応づける。
   データ生成前の計画の存在、推定前に保存入力から問題を構成する処理、凍結ソースからの復元は
   [既存v2計画](validation/MML_FREE_SD_STATIONARITY_V2_ROADMAP_20260811.md)の未完了条件として保持する。
2. **N2の共同不確実性とG2の比較資格を対象別に仕上げる。** 小規模対照で数式と実装を確かめてから、
   新規データで反復研究を行う。N1はその結果を現在の分析へ結び付ける範囲で進める。
3. **G3の残課題を修正し、G4の証拠を揃える。** G2–G4の正式な通過順は維持するが、
   G2/G3自身が必要とする回復・被覆率・誤検出率の試験を、G4開始まで延期しない。
   独立したUX1の改善・D1の画面接続は進められる。
4. **公開する範囲の受入を行い、全計画ではG5→G6へ進む。** 修正済みのローカルコードを
   公開済みと扱わず、証拠・移行・操作確認を伴ってリリースする。

E1（外部結果読込）は、対応する版・モデル・入力・座標・得点法を照合し、改変・不一致・
証拠不足を検出した一往復ができれば、その連携範囲を完成とする。内部の独立検証は必須だが、
全利用者向けの外部結果Viewerを単体分析アプリの完成条件にはしない。
G5を含む全計画では、まずRSM/PCMの二要因固定交互作用だけを扱い、条件付きscreenと共同推定を分離する。
JAX・GPU・Stan/NumPyro・多次元・高次交互作用は後続拡張であり、標準版の完成を無期限に延ばす条件にしない。
一方、対応規模での待ち時間・メモリ・操作性の確認は標準版にも必要で、採用技術と分けて判定する。

#### TAMとの数値照合の完了条件 — 2026-10-04

**RのTAMとの照合はN2・G2/G4の独立検証に含める。** E1の結果読込画面が完成するまで
待つ必要はない。既存のRランナーと保存入力を再利用し、Python MMLとTAMの多相MMLを
共通条件で比較する。既存の[PCM実機照合と後続監査](docs/mml_cross_language_development_check.md#pcm連続積分での再推定とtam実機照合)
には尤度・得点・積分範囲・停止／更新時点の確認があるが、SE出力の保存はSEの妥当性確認とは扱わない。

| 段階 | 揃える対象と残す証拠 |
|---|---|
| モデル対応 | 同じ回答とPersonの集約、カテゴリー、facetの符号、設計行列、同定・アンカー、母集団平均/分散、重みの意味を固定。最初は単位重み・無罰則RSM/PCM。固定したパラメータで反応確率が一致することを先に確認 |
| 尤度・積分 | 同一点での周辺尤度を比較し、求積点・重み・有限範囲・正規化・尤度定数を記録。共通の有限求積の一致と、独立連続積分への近さを別に判定 |
| 校正・停止 | それぞれの返却点を共通目的関数で再評価し、座標変換後の構造値・母集団SD・勾配・終了理由を保存。ソフト固有の停止や更新順序による差を、数学的対象の差と区別 |
| 個人得点 | 同じ最終校正値でEAP・条件付き事後SDを再計算して比較。ネイティブ保存値も残し、更新時点の差を隠さない。WLE/MLEとの単純な一致は要求しない |
| 不確実性 | TAM通常出力と `tam.se()` の対象パラメータ・固定する量・情報行列の定義を確認。座標変換のJacobianで共分散も変換し、比較可能なブロックを独立情報行列・プロファイル等と照合。自由SDとの結合を含むかを確認し、最終的に新規反復研究で被覆率と失敗率を評価 |

[TAM公式マニュアル](https://alexanderrobitzsch.r-universe.dev/TAM/doc/manual.html)の
`tam.se` の例にも、通常出力とのSEの差が示されている。列名が同じという理由だけで同じSEとは扱わない。
今回確認したオンライン文書は4.4-2、保存済み実機照合は4.3.25である。比較には実際のパッケージ版・
ソース・R環境を固定し、新しい版の結果で旧結果を置換しない。このロードマップ整理時点では
パッケージ更新・新規fitを行っていなかった。その後の実行は下の続行記録に区別して記載する。

完了には、対応する正例と境界・欠測・疎な設計の扱い、事前に定めた許容差、全実行の結果・失敗理由、
再実行可能な入力とスクリプトが必要。TAMとの数値一致に加え、解析的対照・独立積分・回復研究を使い、
不一致はモデル対応、積分、停止、得点時点、SE定義に切り分ける。ソフト間の一致だけでG2/G4を通過させない。

各作業の終了時に「解決した利用上の問題／根拠／未確認の範囲／次の一件」を既存台帳へ記録する。
進捗は機能数や総テスト数の割合ではなく、適用範囲が明確になった出力、減った解釈誤り、
実測した待ち時間、復元可能な分析で判断する。人による評価がない項目は未測定と記す。
文献が必要な数理判断では[既存のZotero参照記録](docs/jax_numpyro_conquest_tam_plan.md#文献仕様の根拠)を起点に原典を読む。
今回の整理は既存コード・実行記録・計画の照合であり、新たな文献レビューや数値資格の取得ではない。

#### G4へ向けた実装・検証の続行 — 2026-10-04

**2026-10-06の実行状況：** 利用者の明示許可に基づき、両T0のSHA-256とnonceだけを
指定先 `https://freetsa.org/tsr` へ送信し、元の照会・T0に対応する時刻署名を検証した。
保存容量の開始条件も満たし、15:17 JSTに凍結ソースから正式計算を開始した。
G3の固定対照はPython 20/20・独立R 47/47が通過し、反復計算を実行中。
G3 144,000反復、続いてMML 102,000反復と18異常入力対照を、最大2 workerで順次実行する。
MMLは待機中であり、両研究の全条件の結果・独立参照・推論資格は未確定。
開始前の承認・容量待ちは解消した。以下の2026-10-04記録は当日の状態として保持する。
実行記録は `validation/generated/g4_formal_execution_20261006`、署名と検証結果は各研究の保存先にある。
正式計算と独立した正規モデルの解析対照では、同一データ校正・独立標本校正・共同プロファイルを
分け、固定曲率のSEと残差の標本分散が異なることを16経路で確認した。
この対照の数値をMFRMの性能値として流用せず、凍結研究の反復数にも加算しない。
モデル比較の候補再推定が母集団SD・自由SD設定を継承しない不整合と、非単位重みの
拒否漏れを修正した。推定呼出し前の契約13検査と既存の出力保護11検査が通過。
続いてMMLのBICは、固定次元・各Personの応答数が有界な範囲で、独立なPerson数を使う
正式研究の定義に統一した。応答数Nを残したまま標本単位を明示し、保存済み結果の旧BICは
監査用に保持する。MFRM固有の小標本補正を導出していないAICcとその重みは利用不能とし、
N/Kだけでの自動推奨を止めた。保存済みの表のキャッシュ・CSV・APA出力まで含む
関連57検査が通過。正式研究で検証する推定経路との接続とG2の資格判定は未完了。

**既存計算の効率改善も並行する（2026-10-06の利用者指定）。** 数学的な対象が確定した部分を
実測し、尤度・導関数・情報行列・得点との一致を確認して改善する。最初の変更では、
grouped MMLの勾配を求める有限和の順序を変え、大きなPerson×求積点×パラメータ配列を
勾配だけの評価から省いた。27条件・各2評価点・情報行列の有無で旧実装との一致を確認し、
既存22検査も通過。N192・Q641/Q1281の18条件では、評価時間比の中央値3.50倍、
追跡した一時メモリのピークは中央値約80%減となった。同時に正式計算が動くMacでの
計算部分の測定であり、推定全体・アプリ全体の速度保証ではない。
さらに保存済み開発入力のN192・K4・Q641の9条件で、各版3回ずつ計54回の推定を比較した。
推定前の存在性監査、固定2段階最適化、共同情報行列・共分散まで含む時間比は
中央値2.20倍（1.52〜2.51倍）。27組すべてが事前の数値一致基準と既存の停留性検査を通過。
入力準備・アプリの報告生成は測定範囲外で、正式な標本性能検証の代替にはしない。
記録は `validation/generated/g4_grouped_fit_efficiency_20261006` に保存した。
実行中の正式研究は凍結版のまま維持する。求積点数・停止条件・許容差・失敗の扱いを
緩める変更とは分け、今後は構成処理・求積準備・最終診断を実測して優先順位を決める。

**追加の数値効率改善（2026-10-06）。** 有限logitの指数・正規化項を共有し、単純な縮約を
行列積へ置換、最適化途中だけ不要なEAP/SD計算を省略、標準正規求積則を最大8次数まで
読取専用で再利用する。共同情報行列・求積点数・2段階の停止設定・採否基準は維持した。
RSM/PCM/GPCM・通常/疎/不均衡・Q31/641/1281の27条件で3版各3回、計243推定を測定。
旧凍結版に対して中央値4.12倍、最初の改善版に対して1.88倍。162組の出力比較が通過した。
さらにG3原版のアプリソースを使った独立開発コホート72条件で、入力生成・両Q推定・
4診断・耐久保存までの時間比は中央値2.50倍（1.46〜4.90倍）。版・環境の各段階照合と
独立R参照の時間は含まず、正式研究全体の所要時間をこの比率で保証しない。
144組の推定出力と全72条件の採否・検出判定は一致したが、全数値へ一律の許容差を
適用した最初の比較は2条件の停止残差で不合格となった。この記録を残し、全条件・
両版の停止点・両次数で同一点の尤度/勾配/情報/事後量を追加照合すると全て通過した。
停止点の差に対する勾配差はHessianによる変化と一致し、正式基準を緩めてはいない。
この追加照合は事後診断であり、最初の不合格結果の書換えや標本性能の資格には使わない。
証拠は `g4_numerical_efficiency_20261006`、`g4_bias_numerical_efficiency_20261006`、
`g4_bias_numerical_equivalence_review_20261006`。両原版アプリに対する数理回帰検査は各26件、
補遺の改変拒否11件、試行版記録・再試行拒否はG3/MML合計6件が通過した。
元の2研究の封印を維持した実装補遺を `g4_efficiency_amendment_20261006` に凍結した
（T0 `7ab8190c672b1ada92d5e3448187f91aff50cb37db3af9b70eeb8255c25901b6`）。
利用者の明示許可後、補遺に対する08:54:24 UTCの時刻署名とnonceを検証し、17:58 JSTに
高速化版で正式計算を再開した。切替前の4,601開始済み試行（4,599終端済み）を保持し、
途中だった2試行は中断扱いで分母に残す。未着手試行にだけ新しい版を記録して続行する。
原ledgerの履歴ハッシュと新規試行の版記録順序を実物で確認した。
実行記録は `g4_formal_execution_efficiency_20261006`。G3を実行中、MMLは後続待機。
公開アプリとG2/G3/G4の推論資格は変更していない。

**アプリ側の高次数求積への接続（2026-10-06続行）。** アプリの共通求積処理に残っていた
NumPy `hermgauss` は、実機のQ641で487重み、Q1281で全重みが非有限となることを確認した。
作業版を、検証側と同じSciPyの標準正規求積に揃えた。標準則を最大8次数で再利用し、
返却配列の変更がキャッシュを汚染しないようにする。尤度・E-step・EAP/SD・事後要約/PVの
全経路でゼロの尾部重みを厳密なゼロとして扱い、不正な次数や事前確率質量を拒否する。
正規分布の既知のモーメント、低次数の互換性、高次数での尤度・勾配・共同情報・得点を
含む35検査が通過した。アプリ入口での3モデル×固定/自由SDの6比較も通過したが、
これはQ15・EM最大4反復の回帰確認であり、十分な収束や標本性能の資格を主張しない。
証拠は `g4_app_quadrature_entrypoint_20261006_v3`。既定の求積点数・推定法や実行中の
署名付き正式ソースは変更していない。公開環境への反映と共同推定経路の統合は未完了。
18:26 JSTの進捗確認ではG3が5,450開始済み・5,448終端済み、新版記録849件。
この確認は件数・履歴だけを対象とし、途中の誤検出率・被覆率等は参照していない。

**全体計測・入力監査の改善（2026-10-06続行）。** 署名済みG3実装を使う独立開発21条件で、
段階ごとのソース/環境照合・耐久保存・9条件の独立R参照まで計測した。Rの360照合は通過。
プロファイラ下では環境照合41.3秒、条件付き診断32.5秒、尤度/勾配29.4秒、共同情報7.8秒、
二重の入力構成18.8秒だった。Rを意図的に多く含む開発測定であり、正式研究の時間予測ではない。
入力監査の行別pandas検索を削減し、元indexが逆順・不連続・重複したときの値のずれ/例外も修正。
監査単体の18比較は既存表と完全一致し、中央値9.41倍。プロファイラを外した72条件の
二版比較では、各1,368保存イベントのpayloadが全て一致し、全体時間比は中央値1.039倍
（0.838〜1.596倍、合計325.1秒→305.2秒）。段階照合を含み、起動・独立Rは含まない。
併行計算の影響と小さい全体効果を踏まえ、正式計算は現在の署名済み版で続行する。
作業版では別の入力整合性修正として、無限大の得点/重みを監査・推定前検査・Python再現
スクリプトで明示的に拒否する。欠測・有限のゼロ/負重みの従来の扱いは維持する。
最終ソースで入力/再現スクリプト22件と既存数理/生成44件、計66検査が通過した。
最終版の監査を全72保存入力でも照合し、既存列・型・採否・理由の一致と新しい非有限フラグが
falseであることを確認した。証拠は `g4_pipeline_profile_20261006`、
`g4_input_audit_efficiency_20261006`、`g4_input_audit_math_20261006_v2`、
`g4_input_audit_review_20261006`。初回の表型不一致も記録を残し、型を維持する修正で解消した。
19:16 JSTにG3は7,179開始済み・7,177終端済み。空き44.24GB、計画上必要23.30GBを確認。
G3/MMLの署名済みソース108/296ファイルは不変。アプリの公開とG4の推論資格は未完了。
この後、以下の限定した共同推定のアプリ統合を進めた。

**共同推定のアプリ接続（2026-10-06続行、未公開・技術検証用）。**
`mml_engine="joint", estimate_population_sd=True` の明示指定で、RSM/PCMとbounded GPCMの
構造座標・log(母集団SD)を同時に推定する入口を追加した。単位重み・無罰則・通常の
同定制約に限定し、アンカー、潜在回帰、dummy facet等は推定前に拒否する。
`mfrm_app/mml_joint.py` の尤度・全共同情報・存在性の片側証明は署名済み実装と同一で、
Person内の依存、ゼロ求積重み、従属log-slopeの境界、2段階の固定SLSQP条件を保持する。
既定のEM経路とUIの選択肢は変更していない。

共同推定点と最終SDを通常のEAP・事後SD・plausible valuesへ渡し、共同情報・数値判定を
技術監査表として保存する。数値停留、求積精度、標本抽出下の推論資格は別々に記録する。
quick/fullダウンロードと再現runnerにも保留を渡す。runnerが自由SD指定を落とし、
推定後SDを初期値に使っていた問題を修正し、元の初期SDと推定指定を引き継ぐ。
条件付き近似SE・Person事後SDの既存ラベルは保ち、共同SE/CIの資格へ格上げしない。
最終ソースで関連113検査が通過した（486.74秒）。別途、既存エンジン分岐・自動
フォールバック・自由SD EMの既存2自己検査も通過した。検査結果・124ソース/資源・
生成runnerと入力/出力は `g4_app_joint_integration_20261006` に保存した。
これは実装の一致と保留の接続を示す開発検査で、標本性能の正式な合格ではない。

19:46 JST時点の正式G3は8,156開始・8,154終端（途中の性能指標は未参照）。
署名済みG3/MMLの108/296ファイルは不変で、正式計算への版切替は行わない。
G4の残件は全反復と事前指定のR/TAM照合、失敗を含む分母・MCSE・被覆率等の集計、
条件・出力ごとの合否判定と保留解除の証拠接続である。今回の実装一致はその代わりにならない。
次の効率化候補は環境照合内の重複メタデータ解析と、出力組立時の重複集計であり、
いずれも計測・同値性確認を先に行う。

**出力組立の重複処理削減（2026-10-06続行、未公開）。**
GPCM・不均衡割付の開発例では、core bundle 1回で確率曲線診断が461回、カテゴリー集計が
359回呼ばれていた。既存の4関数は本体・引数・数理式を変えず、2つの純粋な集計関数を
1回の出力生成中だけ再利用するようにした。キーは引数全体のバイト表現を使い、値・型・
index・設定の変化を検出する。返却表はコピーし、最大8件で、終了/例外時に必ず破棄する。
表示用の省略文字列やobject identityに基づく判定は使わない。

再利用機構6検査が通過。3モデル×3割付、JMLE、固定SD MML、192人条件を含む12入力で、
core/full計14組・1,036表の値・型・index・属性・CSVと返却コンテキストが完全一致した。
計測器なし・実行順を交互にした時間比は中央値1.647倍（1.006〜2.289倍）、
合計756.4秒→491.3秒。192人RSMでは101.35秒→100.72秒で実質的な効果は乏しい。
いずれもheadless実行・空のsession state・正式workerとの併行負荷がある開発測定であり、
公開ブラウザの待ち時間や正式反復の高速化とは区別する。証拠は `g4_export_efficiency_20261006`。
20:34 JSTの正式G3は9,966開始・9,964終端、署名済み108/296ファイルは不変。

**設定・結果・証拠の照合とUIの母集団SD（2026-10-06続行、未公開）。**
1,501要素の配列の777番目の変更を省略表示が隠す問題に加え、表の欠測と文字列`<NA>`、
数値と文字列、列名の型が旧ハッシュで同一化される例を確認した。型・全要素・shape・
index・カテゴリー定義を使うv2照合へ変更し、未対応の型を表示文字列へ置換しない。
AnalysisIDは保存済み設定ハッシュを信用せず、現在の設定、正規化した応答、推定座標、
展開パラメータ、得点・推定表・統計要約に結び付ける。感度分析の旧形式または不一致の
証拠はreadiness・CSV・JSONで再利用しない。感度分析画面のキャッシュと専用ZIPも同じ
AnalysisID照合で保護する。後から付加するSE/CIや実行時間は推定点の
同一性から分ける。自由SDの保留文言も特定の最適化エンジンに依存しない表記にした。

UIが推定後に`population_prior_sd`を初期値で上書きする不具合も修正した。入力値は
`ui_population_prior_sd`として保持し、推定SDを診断・求積へ渡す。24人・Q31の開発4例
（RSM EMと3モデルの共同推定）では旧処理により事後得点の再計算が最大0.0436〜0.0939
logitsずれ、修正版は得点・事後SDとも元の計算と完全一致した。これは実際のUI後処理を
抽出して実行したAPI検査であり、公開ブラウザを操作した検査ではない。
最終版の関連43検査が通過した。変更途中の共同推定・生成runnerを含む79/49検査のログも
ソース別に保持し、最終版だけの検査数に合算しない。証拠・ソース・回帰検査は
`g4_result_identity_review_20261006` に保存した。
21:08 JSTの正式G3は11,344開始・11,342終端、空き40.02GB。署名済み108/296ファイルは
不変で、途中の統計指標は参照していない。全反復・事前指定の独立R/TAM照合・分母/MCSEを
含む条件別判定は未完了であり、G4完了や推論資格の保留解除は宣言しない。

**MMLの最終証拠レビューへの接続（2026-10-06続行）。**
`validation/mml_population_final_review.py` を追加した。署名済み計算を変更せず、全予定RunIdの
終端を先に確認し、そろった後に凍結済み集計関数でSQLiteから再集計して保存済み要約と照合する。
18負例対照、全校正座標、Personの2基準×5区分、事前指定した独立R参照とTAMの両グリッドを
条件別にまとめる。欠落・重複・分母違い・保存内容の改変を拒否し、TAMの非ゼロ終了も参照未完了とする。
条件付き被覆・条件付き帰無棄却の登録基準は変更せず、既存の検定不能結果の補完範囲を併記する。
いずれの条件が通っても、このレビューだけでは公開・推論資格を発行しない。

関連54検査が通過した（12.87秒）。別の開発1例（RSM、24 Person、Q641/Q1281）では、実際の
独立R・連続積分・TAMの6照合が通過し、新レビューの判定が一致、参照ファイルの不変も確認した。
正式MML計画への読み取り専用起動では、時刻証明・ソース/環境を確認後、102,018試行が未実行で
あることから性能指標を読まずに終了した。正式計算の順序はG3→MMLを維持する。
証拠は `g4_mml_final_review_20261006`。21:34 JSTのG3は12,178開始・12,176終端、空き36.25GB。
署名済み108/296ソースは不変で、途中の誤検出率や被覆率は参照していない。

**G3の最終証拠レビューへの接続（2026-10-06続行）。**
`validation/bias_population_final_review.py` を追加した。全144,000試行の終端を先に確認し、
そろった後に凍結済み関数で各条件・4つの校正/Person基準を再集計する。独立R参照は保存済みの
数値から比較を作り直し、実行終了・入力/出力hash・台帳・全必要枠・保存要約との一致を確認する。
登録済みHolmの対象集合、失敗分母、部分的に観測されたfamilyの棄却、oracleだけの回復対象を
維持し、合否の基準や推論資格を変更しない。実行中はcontroller lockにより最終CLIを停止する。

最終ソースの関連36検査が通過（9.96秒）。生成開始前に保存された固定対照のPython側20件と
独立R側47件を再確認し、元ファイルは不変だった。一時コピーの5種類の改変もすべて拒否した。
最初のテストのSQL列名と確認用スクリプトのnamespace importを修正したログも保持する。
再推定・R参照の再実行はなく、固定対照を母集団反復の証拠に算入しない。
21:55 JSTの件数確認ではG3は144,000試行中12,768終端で、未完了のためendpoint payloadを
読まなかった。証拠は `g4_bias_final_review_20261006`。全反復・独立参照・条件別の科学的判定は
引き続き残件であり、G4完了や公開保留の解除は宣言しない。

**SE/区間の表示キャッシュの同一性（2026-10-06続行）。**
GPCM MMLのFACETS形式表で、キャッシュが推定パラメータとCI水準だけに依存していた点を修正した。
現在のAnalysisID、最適化座標、表示する元の表の全内容、CI水準を照合し、不十分な入力では
キャッシュを使わず計算する。関連47検査が通過した。数値計算・推定・デルタ法の式は変更していない。

固定SDの開発例では、48 Personの各応答パターンを別IDとして複製した96 Personの計算を比較した。
同じ推定点でHessianが2倍（残差最大1.31e-9）、SEが約0.707107倍になるのに、旧キャッシュキーは
同一だった。修正後の実際のUI分岐は新しいSE/区間を返し、各条件の再表示では再計算しなかった。
固定SDの数値をキャッシュした後に自由SDへ条件を変えた場合も、既存の保留が働きSE/区間は欠損となる。
これは計算上の人工対照で、応答の複製を独立な情報や検証反復には数えない。
確認スクリプトの列名と保存形式の不具合2件、同じ開発入力での再実行履歴も保持した。
証拠は `g4_fair_se_identity_20261006`。正式G3/MMLの計画・凍結ソース・資格保留は変更していない。

**Fair-average単独出力の資格表示（2026-10-06続行）。**
共分散監査の資格情報がFACETS形式表のSE/区間・単独CSVに引き継がれていなかった経路を修正した。
既存契約のTECHNICAL_ONLY/WITHHELD、理由、技術出力可否、確証的結論の可否を各SE/区間に付け、
正則化や欠落を合格表示に変えない。GPCMであっても非MMLに外部の共分散を渡す場合は適用外とする。
日英の画面説明も、返却点の観測情報・Person得点を固定した計算・正則化・被覆未検証を明記した。
旧い「MLEでのFisher情報の逆行列」という断定を使用せず、新注記を載せるためキャッシュをv3へ更新した。

関連50検査（65.15秒）と、実際の画面/CSV分岐の4条件×日英の8確認が通過した。
保存済みの開発例で数値が不変、自由SDは欠損、人工的な階数欠損は保留となり、各条件の再表示では
再計算しなかった。新たな推定・データからの観測情報の計算・R参照計算は行っていない。
階数欠損の対照では、保存行列から作った人工行列を1回逆行列化した。証拠は
`g4_fair_se_qualification_20261006`。22:36 JSTの正式G3は13,771開始・13,769終端、空き32.24GB。
署名済み108/296ソースは不変で、正式反復・独立参照・G4の条件別判定は未完了である。

**MML独立R参照の保存証拠の追跡（2026-10-06続行）。**
最終レビューに、有限求積・連続積分のR入力、終了コード、結果ファイルと実行台帳の照合を追加した。
凍結済み比較関数で数値差と判定も再構成する。欠落・異常終了を合格にせず、元の必要比較数を維持する。
保存済みRSM N24の開発対照では、Rの4比較が元の数値差・判定と一致し、TAMの2枠も従来の
ファイル確認を通過した。別の一時コピーでR結果を壊すと、旧レビューは通過し、修正版は拒否した。
関連55検査（12.26秒）が通過。新しい推定・R実行・連続積分はなく、真値点の有限求積評価のみ
元の比較どおりに再計算する。前段の6比較を新しい独立照合や母集団反復として加算しない。
証拠は `g4_mml_reference_artifact_review_20261006`。22:48 JSTの正式G3は14,076開始・14,074終端、
空き31.62GBで、108/296の凍結ソースは不変。全反復後の条件別判定とG4の完成は未完了である。

**MMLの被覆率と返却率の区別（2026-10-06続行）。**
最終レビューに、固定された校正パラメータごとの区間返却数、返却区間に限った被覆率、全予定試行に
対する被覆区間の返却率を明示する補足を追加した。区間がない試行を、返却区間の非被覆に置き換えない。
未返却部分の二値的な補完範囲も記述するが、実際には返していない区間の被覆率や新しい合否基準にしない。
関連74検査（16.26秒）が通過し、全生成失敗を表す人工台帳3条件・57校正対象でも、分母と保留が残った。
判定基準・正式集計・推定は変更せず、Personの事後区間とは区別する。証拠は
`g4_mml_coverage_denominator_review_20261006`。23:07 JSTのG3は14,476開始・14,474終端、
凍結108/296ソースは不変。途中の被覆率等は未参照で、G4の最終判定は引き続き未完了である。

上のロードマップ確認後、利用者の指定した通常・疎・不均衡な割付で開発対照を実行した。
夜間を含む数日規模の実行が許容されており、追加の規模確認は不要。

- **推定前の構成処理：** 単位重み・無罰則RSM/PCMを未推定データから構成。
  既存アプリとの尤度・勾配・digest一致と、不正入力の拒否を確認した。保存バイトからの
  registrar再構成も開発実装し、入力・設定・Pythonソース/環境を実物で確認する53テストが通過。
  外部TSAの署名付き時刻と改変/nonce違いの拒否を開発用文書で確認した。
  その後、高次数母集団研究のT0と単一試行の実行管理を実装した。T0は凍結済みだが、
  正式な時刻署名は外部送信の承認待ちで、新規資格データは未生成。資格フラグは変更していない。
- **共同不確実性：** 全共同情報と展開後の共分散をRの独立式で12ケース照合した。
  初回2ケースはR側の差分精度で不合格。その結果を保存し、同じ返却点・同じ許容差で
  全ケースに絶対刻みを適用した事後診断は12ケースで通過した。被覆率の根拠へは転用しない。
- **積分の開発対照：** 全36点と独立Rの24点照合が完了。Q241は全6条件・両初期値の12点で、
  スコア、従属水準を含むSE、初期値差の基準を通過した。展開後SEの最大相対差は約9.962×10⁻⁷。
  Q31〜Q181の未達14点も保存し、Q241での一致を他の設計や被覆率へ一般化しない。
- **TAM native再推定：** 同じ6データ×2固定gridの12fitについて、推定値・連続スコア・
  返却点で再構成した共同SEが基準を通過。保存済みEAPは1点で1×10⁻⁸基準をわずかに超えた。
  原判定を保ち、全12点をTAMで再得点する診断で、保存済み事後分布と返却点の時点差を確認した。
- **回復集計：** 制約でゼロになるブロックbiasを水準別biasへ変更し、独立反復単位のMCSE、
  中心化後の共分散によるCoverageSEを実装。共分散が対応しないPerson・JMLE・自由SD等の
  被覆率は評価不可とする。関係する77テストが通過した。
- **回復真値・失敗分母：** PCMのstep平均とGPCMの全体尺度を同定された真値に反映し、
  生成対数オッズが不変であることを確認。収束例ゼロ時の集計fallbackを削除し、全試行・失敗・
  除外の分母とMCSEをCSVにも保存する。固定したPerson値の回復を母集団SD回復とは呼ばない。
- **step・傾きの回復：** step・log-slope・slopeの回復行、従属座標のSE、正値を保つslope区間、
  真値四分位別biasを実装。固定パラメータを回復成功数から除外する。生成モデルv4は実際の
  step形状差を持ち、固定MMLのSDも同定後の生成尺度へ合わせる。関係する81テストが通過した。
  母集団回復・SE/区間・bias screenの凍結反復研究は未完了。
- **BIAS-004：** 固定校正の同一問題の反復を削除し、各セル1回の計算記録へ変更。関連61テスト
  通過。その後の境界・失敗・参照値方針の修正は次項。G3の推論資格は未完了。
- **BIAS-003とセルの数学的契約：** 条件付きスコアの単調性を使って有限内点解を確認し、
  無限MLE・探索範囲・失敗・重みを区別。RSM/PCM/GPCMで一次参照値を条件付きLR候補に揃え、
  各観測セル集合のHolm補正へ接続した。失敗・除外セルも比較数に保持し、BHは記述値に限定。
  根拠未確認のt参照と「固定校正＝プロファイル尤度」の説明を撤回した。
  独立Rの15ケース・多重比較8条件で照合済み。初回のR二階差分の丸め誤差13件を保持し、
  入力・許容差を変えず3刻み幅で再確認した。計算一致を被覆率検証と区別し、G3保留を
  画面・CSV・報告書・旧結果にも適用する。詳細は数理監査の2026-10-04続行記録を参照。
- **GPCMの全共同情報：** facet・step・log-slope・log(SD)の独立R照合は、3割付×Q31/Q61の
  6点で通過。ただし連続積分との精度基準は6点とも未達で、通常割付の展開後SE差は
  Q31で最大約41%、Q61でも約8.8%だった。全情報を求める独立連続参照と裾誤差評価を追加し、
  不合格を保持した。Q121/241/321の再推定と独立照合が完了し、Q321では3割付×2開始点の
  6点すべてが連続積分・共同情報・展開SEを含む基準を満たした（展開SE差の最大1.83e-7）。
  Q121には未達が残る。次数間の一致だけで精度合格としない。
  同じ有限求積でのソフト間一致を、連続尤度や区間被覆の資格へ読み替えない。
- **反復研究用の有限求積推定：** Person・セル・カテゴリの度数から再構成した同一尤度を、
  保存18点で元アプリ・独立Rと照合した。未使用カテゴリ等の構造方向を整数不等式で確認し、
  有限MLEの不存在を証明できる入力では推定を止める。Q321/Q641での初期値からの18推定も
  全数値基準を通過。構造方向が見つからないこと自体は有限MLEの存在証明ではない。
  開発用の計算であり、公開経路の変更・自由SD区間の保留解除はしていない。
- **母集団反復用の生成器：** 各反復でPersonを独立に抽出するRSM/PCM/GPCMのDGPを追加。
  能力・割付・応答の乱数を分離し、反復数や実行順で既存標本が変わらない128-bit seed方式を確認した。
  18テスト通過。未観測水準・同定不足を補充せず記録する。正式な回復・被覆率研究のデータは未生成。
- **母集団回復の集計：** nativeのCriterion位置を保持し、step・log-slope・slope・SDの
  共同共分散を展開する。真値分位は母集団の固定境界を使い、Person数が反復で変わる分位群では
  合計の比と反復単位のMCSEを用いる。区間条件付き被覆、既知真値全体での被覆、推定・区間の
  利用可能率、真値自体が生成できなかった反復を分ける。これは集計式の実装で、被覆率の実証ではない。
- **母集団反復研究の凍結：** 51条件×2,000反復と18異常入力対照、Q641/Q1281、
  事前指定した独立R/TAM照合、失敗分母・MCSE・受入幅・ソース/環境を
  `validation/generated/joint_mml_population_20261004_v1` に固定した。
  開発対照の全16比較と実行経路のテストは通過。外部時刻署名サービスへのハッシュ送信が
  自動承認審査で拒否されたため、指定送信先と内容の承認を依頼中。本試験は未開始である。
  読み取り専用の補足点検では、検定可能な反復での条件付き棄却率と、全予定反復あたりの
  観測棄却数を区別した。検定不能な結果の補完範囲を表示する13テストが通過。
  凍結済みの計算・受入基準は変更せず、条件付き合格を全予定反復の誤検出率保証へ読み替えない。
- **G3の有限標本対照：** 二項応答の全カウントを列挙し、既知校正と同一標本の切片校正を
  区別する開発計算を追加した。N=8・成功確率0.5では、現行のχ²近似による全標本あたり
  の棄却率は6.25%、推定不能を除いた条件付き率は16/254。同じ標本で切片を推定する
  対照では効果が吸収され、非ゼロ効果でも検出力は0。一般MFRMの率をこの例で代用せず、
  G3保留を維持する。
- **JMLEの付随パラメータ対照：** 1人2項目のRaschモデルで、人数だけを増やしても
  無制約JMLEの項目差が真値の2倍へ近づくことを直接導出した。40条件の二項全カウント列挙、
  独立Rの4組の行尤度/共同情報照合、16テストを完了し、推定差・SE・被覆を図示した。
  帰無点でも95% Wald区間の被覆は漸近的に約83.42%。一般MFRMアプリの性能値とは区別し、
  JMLEの自動モデル選択や推論の資格へ昇格させない。
- **G3の生成効果と残差：** 生成時のRater×Task shiftを、推定校正後の残差と区別する経路を
  実装。真の/推定した校正と、真のPerson/EAPの4基準を18開発データで確認し、432セルと
  144多重比較補正が独立Rとの許容差を通過した。2/6カテゴリーの追加照合も通過。
  失敗・未観測セルを非棄却としない反復集計を実装し、条件付き率と不明結果の補完範囲を分ける。
  全予定試行の完了まで集計しない経路、再試行禁止、校正失敗でも既知値の対照を保持する経路も実装。
  未補正p値のfamily誤検出率を5%とみなさず、正方向の検出・実質差・疎セルを分けて評価する。
  N192・6カテゴリー・不均衡の追加3開発データでも、72セル・24補正が独立R照合を通過した。
  R側で未観測・極端な水準・p値・多重比較を再構成する参照経路と、参照失敗を全予定分母へ
  残す経路を追加。負例で見つかった帰無点近傍のLRの桁落ちは、許容差を変えず同値式で修正した。
  固定閾値でのセルMLEと、閾値を推定するモデルの存在性を区別する対照を含む103テストが通過。
  実行管理と条件別の判定集約を追加し、関連35テストが通過。署名・隔離ソース・実行環境・容量・
  固定対照を確認してから生成し、中断試行の再実行を拒否する。条件やセルを平均して不合格を隠さず、
  独立R参照が欠けた条件は基準達成扱いにしない。
  72条件×2,000反復の仕様・seed・107ソース・環境を
  `validation/generated/bias_population_operating_20261004_v1` に別途凍結した。
  時刻署名は未取得で、正式データは未生成。開始時に22.6 GiBの空きを必要とする計画に対し
  現在は約13.5 GiBのため、保存先の指定も依頼中。新規の誤検出率・検出力・被覆率研究は未実施。

TAM 4.3.25の `tam.se()` は、今回の固定母集団の対照で対角曲率の逆平方根と一致し、
共同情報行列の逆行列から得るSEとは異なった。SEの列名だけで同等性を判断しない。
結果・失敗・出典・再現物は[数理監査の続行記録](docs/mml_population_sd_uncertainty_audit.md#推定前の問題構成とrによる共同情報の照合--2026-10-04-続行)にまとめた。
**G2/G3/G4はいずれも未完了。** 公開環境への反映や保留解除はしていない。

## 開発・数値検証の経過 — 2026-09-11〜15

<details>
<summary>過去の判断と数値結果を開く（現在の着手順は上の実行計画を参照）</summary>

Development-line consolidation, 2026-09-11: use
[`docs/development_integration.md`](docs/development_integration.md) for the active
checkout, preserved public baseline, and staged integration status. The original
scientific and product exit gates below remain in force.

Extension refinement, 2026-09-13: prioritize the numerical/inference contracts
of existing JMLE/MML and measure their performance before JAX acceleration.
Use the local mfrmr 0.2.4.9000 development evidence by target and condition;
it is not blanket validation of this Python application. Stan and NumPyro/JAX
are planned as selectable offline code packages, after a common-model audit
and the existing statistical gates. Their retained generators and Posterior
Viewer are not currently public routes. See the bounded work queue and B0–B4
conditions in [`docs/jax_numpyro_conquest_tam_plan.md`](docs/jax_numpyro_conquest_tam_plan.md).
The next UI work refines action-local feedback, visible selections, diagnostic
navigation and task-based export; its current-state evidence and acceptance
tasks are in [`docs/ux_adversarial_audit.md`](docs/ux_adversarial_audit.md#next-ui-refinement--2026-09-13-planned).
This revision updates internal planning only and does not reopen an estimator.

Numerical status, 2026-09-14: the observed PCM stress fixture has a checked
continuous-integral reference (SD approximately 1.6541712 from four starts).
Finite GH solutions still depend on order. TAM's range/density audit shows that
wider tails matter; nine fits with aligned stopping settings terminated. Earlier
strict stopping failures and stored-EAP timing discrepancies remain recorded.
ConQuest 5.47.5 runs through the user's Terminal; the agent's library-load failure
is route-specific. Matched fixed-grid estimates are close to TAM.

ConQuest scoring is a separate precision target: 31 fixed-calibration runs showed
EAP RMSE 0.014242 → 0.000857 and posterior-SD RMSE 0.007443 → 0.000428 logit as
p_nodes rose from 2,000 to 200,000 across ten seeds. A two-run control now isolates
the parameter-count display: removing the explicit mean anchor under CASES changes
-1 to 0 while calibration, deviance and scores remain identical. Preserve the old
failed metadata checks; use CASES alone for this fixed-zero-mean scoring condition.

The 1,543-row native GH replay confirms a one-row likelihood/coordinate timing
offset. Tiny final updates coexist with finite-GH log-SD gradients of approximately
1.827, -3.549, 0.374 and 0.03685 at Q31/61/121/181. Across all 1,539 adjacent
transitions, centering the previous GH posterior by its mean reproduces the next
variance within 1.111e-6. Unlike the TAM source audit, this is a post-hoc numerical
reconstruction, not proof of ConQuest's internal implementation.
Neither stopping nor a flat SD trace establishes likelihood stationarity.

The direct-start experiment is now complete: eight runs paired one iteration with
continued iteration at each order. All four first updates leave the finite-GH
stationary starts and increase their finite NLL; predicted centered variances match
within 3.423e-7. Low-order native endpoints depend strongly on starts. Q181 ends at
the same exported SD as before, but retains a log-SD gradient near 0.03681 and
same-point EAP/posterior-SD integration errors near 0.00118/0.00210 logit.
Returned checkpoints also differ from the minimum reported NLL over these histories.
Sixteen same-coordinate evaluations agree with independent R within 5.69e-13.

Structural replay now explains the observed update closely: old posterior weights,
new variance-scaled nodes and one native-order coordinate sweep reproduce 19 saved
transitions within 1.946e-6. Five sweeps do not reproduce the moving first updates.
Four GH61 controls confirm that keep-last changes exports without changing history;
keep-best does not mean the minimum reported NLL over the entire history.
The improving-iteration candidate hypothesis fits saved NLL values, but remains a
post-hoc description rather than a native implementation contract.

Independent-condition validation now covers four freshly simulated, unedited PCM
datasets (80 persons each; generating SD 1.0 or 2.5, two seeds each). Protocol and
engineering precision targets were fixed before generation. R and app arithmetic,
continuous-reference refinement and optimizer-gradient checks pass. Ordinary cases
meet the scoring target at GH121; neither wider-tail case meets it even at GH181.
Their GH181 EAP/posterior-SD errors reach 0.09846/0.11939 logit. One GH31 case reaches
different low-gradient solutions from two starts. A fixed grid spanning +/-12 at
0.1 spacing also misses the posterior-SD target in both wider-tail cases.

The separate range/density check is complete: halving spacing does not rescue
the +/-12 wider-tail grids; +/-16 meets the stated engineering targets at both
spacings, and +/-20 reduces the largest posterior-SD error to about 2.1e-9 logit.
Eight matched-start TAM fits and eight actual ConQuest runs are now reviewed.
Their free-coordinate differences reach 3.992e-5 and deterministic EAP differences
3.975e-5 logit. Native ConQuest MC scores remain a separate precision target.
ConQuest's prior-grid normalization changes the free-SD objective: its correction
is not constant in SD, even when calibration estimates are close to TAM.

Keep the new fixed-grid overall implementation flag false: one of 96 direct
optimizer calls returned ABNORMAL on restart, despite a small gradient and
unchanged coordinates. A TAM history-recorder failure was separately preserved
and corrected without changing estimation settings; all eight saved native fits
pass their stopping checks. Independent R arithmetic passes throughout.

The ABNORMAL restart is now reproduced with unchanged options. Native task capture
shows exhaustion of 50 line-search requests (36 actual objective calls). Independent
60/90-digit arithmetic finds a real trial improvement of 4.98e-14 that float64
reports as a 4.55e-13 worsening. A diagnostic Newton correction is at most 1.26e-7;
its true NLL improvement is 5.64e-13. This identifies rounding interference in this
case, without converting the raw failure into success or qualifying inference.
The existing termination/stationarity/quadrature contracts remain fail-closed;
their 34 tests pass. Successful controls also miss the requested 1e-8 gradient
tolerance, and one successful small-gradient fit fails the integration target.

The separately specified likelihood-difference arithmetic now passes 22 comparisons
at four anchors, including off-solution controls, against 60/90-digit references.
An exact fixed-grid probability-ratio identity recovers the small improvements;
the retained line-query/Newton differences have absolute errors below 2.65e-22.
Large changes, changes in SD, reversed comparisons and composition also pass.
This is a validation-only evaluator, not a new optimization run or stopping policy.

The fixed-anchor optimizer comparison now retains all 96 calls. With a near-solution
anchor and ftol=0, all eight primary/restart pairs reach requested gtol=1e-8;
their final gradients are at most 9.824e-9 and cross-start coordinate differences
at most 4.616e-9. Raw NLL and zero-anchor arms do not reach that gradient target.
Positive historical ftol also stops the near-anchor arm early. Six ABNORMAL calls
remain recorded across controls; execution/arithmetic success is a separate flag.
The +/-12 wider-tail cell still misses the integration target despite tighter
optimization. All 96 returned NLLs/gradients replay exactly; selected near-anchor
improvements agree with 60/90-digit arithmetic within 1.03e-20.

The two-stage procedure now passes on four fresh, unedited 80-person datasets,
with the protocol frozen before generation. Each run constructs its anchor from
its own preliminary endpoint. All 24 trajectories (three starting SDs, two grids)
pass the refinement diagnostics; all 48 refinement calls meet gtol=1e-8, while
none of the 48 preliminary calls do. Cross-start coordinates differ by at most
7.926e-10. All 96 calls terminate normally in this new check; historical failures
remain unchanged. Independent R and 60/90-digit comparisons pass. The wider-tail
datasets still fail integration targets on +/-12; +/-20 passes for these cases.

The moving-GH likelihood-ratio identity now passes 48 comparisons across two retained
datasets and orders 31/61/121/181. It includes theta changes through expm1(delta log-SD)
and no extra density/Jacobian term. Independent R, 60/90-digit NLL differences,
analytic log-SD derivatives and all-coordinate finite differences pass. At unchanged
wider-tail coordinates, GH181 still differs from continuous integration by 0.0890
in EAP and 0.1056 in posterior SD. Arithmetic agreement is not integration accuracy.
See [moving-GH arithmetic](docs/mml_cross_language_development_check.md#sdとともに動くgh点の尤度差).

Moving-GH two-stage fitting now passes all 48 trajectories on the four retained
datasets (three starting SDs, four orders). All 96 refinement calls reach gtol=1e-8;
the 96 preliminary calls terminate normally but miss it. Same-order cross-start
coordinates differ by at most 8.181e-10. Independent R and 60/90-digit checks pass.
The wider-tail fitted SDs remain nonmonotonic across orders: 2.527/3.054/2.543/2.768
and 2.198/2.473/2.895/2.621. GH181 posterior-SD integration errors remain 0.1046
and 0.0682. Finite-GH stationarity does not certify continuous-integral stationarity.
See [GH fitting and order sensitivity](docs/mml_cross_language_development_check.md#ghでの二段階推定と次数依存).

The optional two-stage API is now connected. It retains preliminary and refinement
runs separately, keeps public objective fields as raw NLL, and records fixed-anchor
differences and the last optimizer query separately. JointPolishOptions now accepts
ftol=0 with its positive default unchanged. Historical sources and hashes are preserved.
Eight PCM API trajectories match the prior solutions within 6.677e-13 in coordinates;
independent R and 60/90-digit checks pass. The unchanged Q-sensitivity assessment
accepts ordinary Q121/Q181 comparisons and rejects wider-tail comparisons. Passing
shifted values as raw NLL is rejected. All 83 related regression checks pass.
See [API integration and preserved evidence](docs/mml_cross_language_development_check.md#二段階推定apiへの接続と旧証拠の保全).

The development adapter now builds the moving-GH difference from the app's actual
observed rows and constrained parameter expansion. RSM/PCM checks include missing
rows, an empty person, zero/fractional response weights, opposite facet signs,
nonzero element/group anchors and person regression means. Independent 60/90-digit
scalar calculations and analytic-gradient checks pass; all 87 related tests pass.
A native RSM fit also completes the four-call two-stage adapter route. This extends
finite-rule arithmetic evidence, not weighted/missing-data inference qualification.
See [observed-row adapter evidence](docs/mml_cross_language_development_check.md#観測行を使う差分関数と開発用アダプター).

Weighted/missing PCM refitting now covers three starting SDs and Q31/61/121/181.
All 12 series pass finite-GH stationarity; all 48 PCM optimizer calls terminate
normally. Within-order coordinates differ by at most 7.426e-9 across starts.
Independent finite-rule arithmetic, continuous-reference refinement, prior scoring
for persons without information, and 60/90-digit differences pass. The Q31/61/121
solutions still miss the integration targets. Q181 meets them for this case, with
maximum EAP/posterior-SD errors 1.56e-6/3.25e-6 logit. The Q121/Q181 assessment
rejects all three pairs because structural coordinates change by about 3.46e-4.
This does not establish a universally sufficient GH order; the prior wider-tail
examples still fail at Q181. See [weighted PCM refits and integral checks](docs/mml_cross_language_development_check.md#欠測重み付きpcmの再推定と連続積分照合).

Continuous refitting from the four retained GH starts now reaches the same checked
local minimum: SD 2.27265279973, maximum computed gradient 2.11e-13 and cross-start
coordinate difference 1.85e-14. Two-step-size curvature checks are positive, and
independent R likelihood/moments agree. Relative to this refitted reference, GH181
has maximum structural/EAP/posterior-SD differences 5.36e-6/5.36e-6/3.24e-6.
GH121's total EAP difference is 5.55e-4, although its earlier same-coordinate
integration error was only 8.60e-5. Recalibration and scoring integration must both
be included in accuracy checks. See [continuous local reference and score effects](docs/mml_cross_language_development_check.md#欠測重み付きpcmの連続積分による局所解).

Independent data now cover two new 32-person response sets, each paired across
complete/unit-weight, missing/unit-weight and missing/response-weight designs.
All 24 GH series pass finite-rule stationarity; all 96 GH calls terminate normally,
and within-order coordinates agree across starts within 1.31e-8. The 12 continuous
refits pass local-reference checks. Independent R and 60/90-digit comparisons pass;
432 retained evaluation points and 101 artifacts/16 source hashes were audited.
All three ordinary designs meet total accuracy targets at Q121 and Q181; all three
wider-tail designs miss them at both orders. At Q181, wider-tail total EAP errors
range from 0.0114 to 0.0303 logit. Stable optimization does not establish adequate
integration. These are two paired numerical probes, not six independent replications.
See [independent observation designs and full results](docs/mml_cross_language_development_check.md#独立データにおける欠測と回答重みの対照).

Native TAM 4.3.25 / ConQuest 5.47.5 comparisons now cover the four complete/missing
unit-weight designs at ranges +/-12 and +/-20 (spacing 0.1). All 16 fixed-calibration
checks pass before 16 native refits. The native coordinate difference is at most
1.761e-5; deterministic EAP/SD recalculated at those coordinates differ by at most
1.703e-5/3.11e-6. Wider-tail +/-12 runs still miss continuous-reference accuracy;
both engines meet the total targets for all four designs at +/-20. Independent R,
96 initial/returned/final-iteration roles and 348 artifacts/22 sources were audited.
The normalized finite-prior NLL correction now accounts for each regression mean.
Stopping reasons and the two fixed TAM calls reaching their three-iteration cap
are retained, as is the first run's overstrict JSON-precision check and exact source.
See [native unit-weight comparisons and limitations](docs/mml_cross_language_development_check.md#独立した重み1条件のtamとconquest実機比較).

Fixed-calibration ConQuest scoring now covers those four designs at p_nodes
2,000 / 20,000 / 200,000 and ten seeds each, plus four exact-repeat controls.
All 124 executions preserve the audited calibration/data/design/likelihood exports;
the repeats also reproduce the previous baseline scores exactly. Across the four
designs, EAP RMSE against independent continuous R integration falls from
0.01097–0.03489 to 0.000830–0.002442 as p_nodes increases 100-fold. The largest
EAP error at 200,000 nodes is still 0.009442 logit. Finite-grid reference differences
are retained separately; 1,635 artifacts and 23 sources pass saved-output replay.
See [regression/missingness scoring variation](docs/mml_cross_language_development_check.md#回帰欠測条件でのconquest標準得点の乱数変動).

The weight-semantics audit now separates response powers inside the integral from
person weights outside it. All 24 fixed native calls completed, but both predeclared
ConQuest weighted-display hypotheses failed in the eight nonunit cases. Retained
post-hoc analysis finds mean-one weights in ConQuest's observed sufficient statistics
(maximum error 2.85e-14), while its fixed-calibration reported NLL matches the
unweighted normalized-prior sum (maximum error 1.81e-7). TAM's mean-one weighted
raw-prior NLL agrees within 2.28e-13. These are distinct quantities. Independent R,
252 artifacts and 29 source identities were audited, including the preserved failed
protocol/code. The inspected mfrmr development source is 0.2.4.9000 and applies
response powers, as does Python; no new mfrmr fit was performed.
See [weight definitions and rejected display hypotheses](docs/mml_cross_language_development_check.md#回答重みと受験者重みの区別および固定値対照).

The one-free-criterion control now covers 36 ConQuest calls. All first steps agree
with weighted scalar EM within 4.89e-7, and a retained post-hoc replay reproduces
99 nonduplicated full-run transitions within 5.79e-7. Keep-last returns agree with
independent weighted scalar roots within 4.89e-7 (1.05e-7 in the nonunit cases).
The default keep-best setting follows the same executed histories but returns
iteration-4 coordinates in all eight nonunit cases; root gaps are 0.00230–0.01063
logit and absolute weighted scores 0.0662–0.1344. Those returned points are not
the minimum displayed-deviance rows, so the internal selection rule remains open.
Another retained post-hoc finding aligns all 400 displayed history NLLs with the
preceding coordinates, not the coordinates written in the same row. Independent R
agrees within 5.41e-13; 547 execution artifacts/26 source identities were replayed,
and the history diagnosis preserves and checks 554 artifacts.
See [one-parameter weighted updates and returned points](docs/mml_cross_language_development_check.md#受験者重み付きの1パラメータ更新と返却点).

The older-history recheck now covers 24 retained runs: 16 GH runs including
duplicate start/return controls, and eight fixed-grid fits on four 80-person datasets.
All 5,139 comparable rows are compatible with the preceding coordinates after
propagating six-decimal structural/variance export rounding. Same-row coordinates
are excluded by conservative rounding bounds in 528 rows; the remaining rows do
not identify timing uniquely. All rows also satisfy the separately recorded
first-order sensitivity allowance. Maximum preceding-point discrepancies are
1.35e-5 for GH and 4.10e-5 for grids; retained independent R values agree within
4.55e-12. Four default-start first rows lack saved starting coordinates and are
excluded from timing comparisons. The previous stationarity and continuous-score
accuracy conclusions remain unchanged; a displayed NLL discrepancy at a returned
point is not by itself quadrature error. Original files are preserved, with 432
retained file identities and 25 new artifacts checked.
See [history timing and rounding recheck](docs/mml_cross_language_development_check.md#過去のconquest履歴の時点と丸め幅の再点検).

The regression/missingness history review now covers all eight retained unit-weight
ConQuest refits. All 1,335 rows are compatible with the preceding native coordinates,
including slope/variance rounding and person-specific prior normalization; 289 rows
exclude the same-row coordinates. Every returned point equals the last saved row.
Independent R agreement is within 5.12e-13 and two-step gradient differences within
1.87e-7. In the wider-tail +/-12 conditions, returned normalized-grid app-coordinate
gradient suprema remain 0.1522/0.2892; even raw-grid beta scores are 0.00568/0.00919.
The +/-20 conditions retain the previous engineering accuracy passes, without
certifying internal gradients at 1e-8 or SE/CI/coverage. Display timing agreement
does not resolve the range or stationarity problem. The audit preserves 348 old
artifacts and checks nine new artifacts with 24 source identities.
See [regression history, gradients and report requirements](docs/mml_cross_language_development_check.md#回帰欠測履歴の時点と返却点の勾配).

The four-question comparison record is now implemented in
`mfrm_app/mml_comparison_report.py`, with separate measurements, explicit targets,
missing prerequisites and computed status per question. No overall equivalence
or inference qualification is emitted. Rendering recomputes decisions instead of
trusting saved status labels. Eight retained comparisons now produce Japanese HTML,
CSV and JSON reports plus a ZIP. Four answers appear first; measurements and
point/objective/provenance details use native collapsed HTML sections.
The two wide +/-12 cases retain stationarity/integration target failures.
Integration meets the reporting targets in the other six cases, while their
internal stationarity remains unassessed because gradient rounding uncertainty is
not qualified. All eight calibration-equivalence and native-score-stability
answers remain unassessed: no equivalence margin or matched-calibration repeated
seed evidence is supplied. Reporting targets are explicitly post-hoc engineering
screens. Ten targeted tests pass; eight reports replay exactly; 34 ZIP entries
match local artifacts; desktop HTML was visually checked. No app estimator or
existing Streamlit UI was changed.
See [comparison report implementation and artifacts](docs/mml_cross_language_development_check.md#4つの問いを分けた比較レポートの実装).

The preparation path is now connected in Streamlit (2026-09-15): Report & Export
→ Check next → Compare with another program. It shows four unassessed questions,
one collapsed explanation and one ZIP download; other review details are hidden
while this resource is selected. Full current settings, response weights, anchors,
population design, returned coordinates and estimates are retained and bound by
fresh content hashes and an analysis identity, including when recorded run IDs are
stale. Missing comparison evidence never becomes a pass. The R worksheet now stops
before fitting: automatic design guesses, rounded response aggregation and the
rank-only/arbitrary-threshold guidance were removed. No external-result import or
TAM/ConQuest/mfrmr matching runner is supplied yet. The prior numerical-study
runtime was archived and its eight reports replayed exactly without changing old
evidence. See [Streamlit connection and checks](docs/mml_cross_language_development_check.md#streamlitへの比較準備の接続).

The first restricted runner is implemented (2026-09-15). Eligible exported PCM
fits now include offline Python/R scripts for TAM 4.3.25 and ConQuest 5.47.5.
Scope is three two-level facets, scores 0:3, unit weights, fixed zero-mean normal
prior, no anchors/penalty/regression, centered other facets and an uncentered step
facet. Duplicate cells, empty persons and combinations missing a score category
are rejected. Both fixed-calibration coding/design/anchor/NLL checks must pass
before optional structural refits. Native files, code and input hashes are retained;
the app itself never launches an external engine or imports a pass.
Twelve retained native calls passed these implementation checks. Complete and
one-response-missing examples have maximum TAM/ConQuest calibration differences
7.67e-7/1.26e-6 logits, while single-seed ConQuest EAP differences remain about
0.0231 logits. Range and spacing probes are fixed-point checks, not universal
quadrature qualification. See [runner scope and findings](docs/mml_cross_language_development_check.md#限定pcmランナーと実行確認).

The fixed-SD stopping/order audit is now retained for the same two examples:
GH31/61/121, two starts, 36 optimization calls and an independent adaptive integral.
All finite-GH endpoints have gradient suprema below 9.22e-10; the largest cross-start
coordinate difference is 8.55e-11. EM continuation still stops with gradients
0.00788/0.00879. Polishing GH31 improves its own NLL but slightly worsens continuous
NLL in both examples: stopping and integration errors partly cancelled.
At polished GH121, absolute continuous-NLL differences are 4.15e-9/3.28e-9, and
calibration differences from retained TAM/ConQuest outputs are at most
3.05e-8/1.29e-6 logits. These observed checks do not qualify Q121 as a default;
continuous gradients still differ from finite-GH gradients. See the
[stopping/order audit](docs/mml_cross_language_development_check.md#固定sdのem停止とgh次数感度の切り分け).

The opt-in fixed-SD development API now uses the shared MML likelihood and stable
likelihood-difference implementation. Disabled calls do no work. Enabled calls
retain the original result and return separate candidate parameters, labeled EAP/
posterior-SD scores, the exact fixed rule, gradients and all three solver records.
Unsupported models/constraints/weights and inconsistent inputs stop before fitting.
The two previous endpoints were reproduced within 3.65e-15 logits; three additional
SD/data checks confirm the need to keep integration separate. On the cyclic test
fixture with SD1.6, GH121 has a stationary finite objective but continuous NLL/SD
differences of 0.00704/0.000775, despite nearly identical EAPs. These candidates are
not fitted-result replacements and do not inherit SEs, CIs or inference readiness.
See the [API checks and counterexample](docs/mml_cross_language_development_check.md#固定sd仕上げ候補を返す開発用api).

The fixed-SD MML view now starts with one short status and a closed numerical-details
panel. Supported PCM results have an explicit returned-point check: it reevaluates
the same fixed-GH likelihood and gradient without refitting. Stopping, gradient and
integration remain separate; a small gradient does not supply an integration pass.
Fresh input hashes invalidate stale checks, while language changes retain them.
All 168 relevant tests passed; Japanese/English desktop and mobile component checks
were completed. Extracting the shared input serializer preserved all 110 assets in
six export comparisons. See the [UI scope and evidence](docs/mml_cross_language_development_check.md#返却点診断のコンパクトな画面接続).

Japanese copy across the main workflow has been revised in 672 locale entries,
including setup, returned-point checks, residual PCA, bias and report/export.
The wording separates a stopping rule from gradient and integration checks;
the old blanket reporting-order recommendation has been removed. Numerical code,
decision thresholds and defaults are unchanged. See the
[Japanese UI review](docs/ux_adversarial_audit.md#日本語の見直し--2026-09-15ローカル実装).

Next connect saved diagnostics and the existing Q-sensitivity assessment using the
same current-input/returned-point checks. Keep application
of candidate coordinates to a full result deferred until every derived table and
uncertainty quantity can be recomputed consistently. Do not automatically increase
Q, mark scientific equivalence or import native-score stability from these cases.
Before displaying measured external evidence for a current fit, verify the data,
full model/constraints/weights, coordinates, runtime and settings, then recompute
the four answers from retained measurements and justified targets. A file hash or
native convergence flag alone is insufficient. Further native-internal reconstruction needs a concrete
unresolved comparison question. Reevaluate the actual returned point
even when the native summary says convergence was reached. Keep nonunit response-weight
native equivalence unestablished; do not substitute person weights for powers inside
each person's integral. Retain
calibration quadrature settings and posterior scoring method/p_nodes/seed separately
in comparison reports. Native calibration agreement does not certify standard-score
precision, and these ten-seed summaries do not establish a universal Monte Carlo
rate or a sufficient default count. Keep finite-objective stationarity separate
from the reported stopping reason. Global optimization, SE/CI/coverage,
default-estimator changes and inference qualification remain unestablished.
Reuse the existing Q-sensitivity assessment and its raw-NLL reconstruction checks.
Finite-order agreement remains separate from continuous-reference and score accuracy.
These observed numerical probes are not SE/CI/coverage qualification.
ftol=0 still permits termination at a computed zero/nonpositive reduction and is
not gradient-only stopping. Use range and density checks together.
Always retain initial, returned and last coordinates and reevaluate each objective.
Further ConQuest internal reconstruction is a secondary diagnostic, not the
definition of correctness. Qualify SE/CI/coverage separately. No default order or
inference gate changes.
See the [numerical evidence ledger](docs/mml_cross_language_development_check.md)
and [anchor control / GH history audit](docs/mml_cross_language_development_check.md#conquest平均アンカー対照とgh全履歴監査).
The [direct-start results](docs/mml_cross_language_development_check.md#conquest直接gh解からの開始値対照)
separate first-step movement, start sensitivity and continuous-integral accuracy.
See [structural replay and checkpoint controls](docs/mml_cross_language_development_check.md#conquest構造更新と返却点の対照)
for the update reconstruction, its limits and the next validation priority.
See [independent ordinary and wider-tail conditions](docs/mml_cross_language_development_check.md#独立した通常条件と広い裾での求積検証)
for the frozen protocol, numerical targets, full results and remaining scope.
See [range/density and matched-start external checks](docs/mml_cross_language_development_check.md#積分範囲と密度の対照およびtamconquestの実測)
for the new evidence, objective-normalization distinction and retained failures.
See [optimizer resolution and status contract](docs/mml_cross_language_development_check.md#異常終了の再現と数値精度の判定契約)
for the bounded high-precision diagnosis and its negative controls.
See [stable likelihood differences](docs/mml_cross_language_development_check.md#尤度差の安定計算と高精度照合)
for the exact ratio identity, arithmetic scope and stopping-scale caveat.
See [fixed-anchor optimizer controls](docs/mml_cross_language_development_check.md#固定基準点を用いた最適化対照)
for the 96 retained runs, gradient attainment and unresolved integration error.
See [two-stage checks on fresh data](docs/mml_cross_language_development_check.md#保存済み解を使わない二段階検証)
for self-constructed anchors, independent arithmetic checks and the API boundary.
The [fixed-calibration scoring comparison](docs/mml_cross_language_development_check.md#conquest校正固定での得点mc精度)
separates posterior Monte Carlo draws from MML quadrature order.

</details>

- Status: active
- Last updated: 2026-09-15 (current execution plan and released baseline reconciled)
- Product boundary: standalone Python; optional offline interchange within validated scope

## Product thesis

The application will become a design-to-decision workbench for many-facet
rating studies. Its primary value is not the number of estimators it exposes.
Its value is the ability to show:

1. where a rating design is fragile before data collection;
2. which model and data conditions support a reported result;
3. whether a conclusion survives reasonable sensitivity checks; and
4. what evidence-based action a researcher can take next; and
5. which intended interpretations and uses remain outside the evidence that
   the application can assess.

The intended workflow is:

```text
Plan -> Check -> Estimate -> Diagnose -> Stress -> Decide -> Archive
```

The active statistical-engine remediation lane is governed separately by the
self-contained
[`docs/statistical_engine_remediation_roadmap.html`](docs/statistical_engine_remediation_roadmap.html).
As of 2026-08-04, G0 output protection and G1 exact identified
parameterization are complete; S2 inference reconstruction is active.
Until its release-qualification gate is complete, identification, inference,
bias-screening, and validation work follows that document's dependency order
and temporary output-safety boundaries. New estimator scope must not bypass
those gates.

Every step must remain executable in Python. The Streamlit interface is the
guided research surface; the statistical and evidence contracts belong in
testable modules under `mfrm_app/`.

`StudyContext` is a first-class input to interpretation, not decorative report
metadata. The same numerical analysis may support different bounded
interpretations for different populations, constructs, purposes, and stakes;
those differences change interpretation and decision identities without
silently changing the fitted AnalysisID.

## Execution roadmap: August 2026 to July 2027

The dates and milestone labels below preserve the original planning context.
Use the September 15 execution plan above for current sequencing and completion;
the applicable scientific and product exit gates below remain in force.

This section is the operational view of the detailed delivery tracks below.
Work is selected by milestone exit gate, not by the number of UI features that
can be added. The default priority order is: prevent unsafe interpretation,
shorten the primary journey, modularize the owning contract, then expand
statistical scope only when validation can travel with it.

### Baseline — consolidated on 2026-08-04

The current baseline includes the standalone-Python runtime and release
boundary, AnalysisIdentity and EvidenceRecord foundations, bilingual Help and
terminology contracts, privacy-safe problem routes, contextual Help return,
phase-aware UX presentation policy, native release checks, and a clean native
test profile. The adversarial UX baseline is recorded in
[`docs/ux_adversarial_audit.md`](docs/ux_adversarial_audit.md).

Baseline verification:

- `make apptest`: 1,068 passed, 62 compatibility tests deselected;
- `--self-test`: 56 passed;
- `--doctor` and `--release-check`: passed;
- default sample and real-data privacy paths: AppTest verified; and
- worktree and roadmap: consolidated into a reviewable commit.

### Repository-only CMLE research lane — started 2026-08-09

The pure `mfrm_app.cmle` structural-calibration core now evaluates exact conditional
likelihoods for additive RSM/PCM structures, conditions out Person parameters,
and fails closed on deficient conditional-information rank. A scaled moment
dynamic program returns the exact gradient and conditional information, so
rank and structural covariance no longer depend on a numerical Hessian. Exact
response enumeration, finite-difference checks, matched
`immer::immer_cml()` complete/incomplete-design tests, and the reused
20-dataset stress run are recorded in
[`docs/cmle_phase0_design.md`](docs/cmle_phase0_design.md) and
[`validation/cmle_phase0_20260809/RESULTS.md`](validation/cmle_phase0_20260809/RESULTS.md).
The first fresh-process scaling envelope is retained in
[`validation/cmle_scaling_20260809/RESULTS.md`](validation/cmle_scaling_20260809/RESULTS.md):
all 14 eligible cases were inference-ready, work and 512 MiB memory guards
blocked two oversized audits before second-moment evaluation, and exact-Newton
polishing resolved two above-threshold BFGS precision-loss endings while
preserving raw optimizer evidence.

This is not yet a public estimator route. Public integration remains blocked
until the remaining in-app/multi-platform scaling gate, anchors,
post-calibration Person scoring, downstream diagnostic
qualification, localization, saved-result schema migration/public disclosure,
and simulation operating characteristics pass the gate in the design record.
The current-schema private archive identity is complete but does not promote
the UI. Exact CMLE will
not be offered for estimated-slope GPCM, and a future composite conditional
method cannot silently replace it.

### Decision-stability and design-stress increment — 2026-08-09

Fit and conditional-bias decisions now use a Streamlit-free numerical
contract in `mfrm_app.decision_stability`. Raw finite values determine the
classification; display rounding is presentation only; exact endpoint
semantics are explicit; and numerical/display proximity is retained as
evidence. The first-read surface, fit scatter, table colouring, casebook,
bias inference audit, one-click results bundle, and full table download now
share this contract. Missing fit statistics remain missing/review evidence,
not a stable decision.
The result overview also exposes one compact interpretation snapshot before
the goal router instead of requiring users to discover the First-read section.

Anchor audits now report unique anchored/unanchored counts and descriptive
coverage share, including a stacked coverage figure. No universal percentage
cutoff was introduced, and non-finite group-anchor values are rejected rather
than silently coerced to zero. The deterministic stress record in
`validation/decision_stability_20260809/RESULTS.md` covers MNSQ and bias
threshold neighbours, planned missingness, empty/low facet-pair cells,
zero-count categories, and anchor coverage from 0% to 100%. It found six
raw-versus-rounded fit classification mismatches and six bias-decision
mismatches in deliberately constructed boundary cases; sparse negative
controls exposed up to three empty and 115 low-count facet-pair cells.

This increment is a presentation/numerical sensitivity gate, not a completed
operating-characteristic study. M3 must extend the same contract to remaining
decision surfaces (reliability, PCA, category, DIF/equivalence), add repeated
data-generating mechanisms for false-positive, power, bias, RMSE, and interval
coverage, and retain conclusion-level sensitivity across failed and
unsupported specifications.

The 2026-08-11 assignment-stress increment now has two deliberately separate
engines. Exchangeable equal-context Person-Rater blocks retain a connected
degree-preserving 2-switch dose path. Unequal Task/Criterion block signatures
use a score-free binary MILP endpoint that exactly preserves every
Rater-by-context row margin and locks observed spanning-tree overlap witnesses.
The observed assignment is audited as a feasible baseline before optimization;
only a certified zero-gap HiGHS optimum with post-solve invariant checks enters
the existing same-method paired JMLE/MML refit contract. This closes the first
unequal-context product gap without pretending that an integer endpoint defines
a continuous assignment propensity. M3 remains open for frozen feasibility and
operating-characteristic matrices, a mathematically qualified nested MILP dose
construction (or an explicit decision to remain endpoint-only), anchor-stratum
interpretation, and an offline large-design runner. FACETS remains a qualified
same-estimand JMLE comparator, not a reference for marginal/conditional
estimator ranking or for values reconstructed from its two-decimal display.

The fresh four-vector assignment preflight is now operationally complete. The
earlier crash diagnosis was too broad: interactive FACETS 4.5.0 and visible
`BATCH=NO` work, whereas every tested hidden-window route (`BATCH=YES` or a
hidden Windows startup state) null-dereferences in the shared Xojo GUI runtime.
The qualified adapter therefore launches visible `BATCH=NO`, waits for stable
nonempty report and per-facet Scorefiles, and closes the owned windows with
`WM_CLOSE`. A separate 260-character auxiliary-output stall established a
second boundary; shared code now enforces a conservative 220-character FACETS
path budget before launch.

The clean short-path supplement passed 12/12 FACETS/Python JMLE calibrations
with zero retries and was joined without changing any original completion
marker or failed FACETS artifact. Original and supplemental native Python JMLE
rows replayed within `4.44e-16`. All operational gates and all four registered
n=4 sign-stability directions passed, with R 4.5.1 reconstruction within
`1.67e-16`. This authorizes freezing a separate confirmation design; it is not
itself confirmation, an estimator ranking, or general FACETS compatibility.
The retained dump and vendor packet document the hidden-window failure mode,
not a general inability of FACETS 4.5.0 to run on this host.

The MML generator boundary is now executable rather than only prose. The
default fixed-Person lane holds posterior/EAP measures as non-truth reference
coordinates. A separate rank-preserving population lane draws normal order
statistics at the fitted fixed/free MML population SD and maps them to strict
source-EAP ranks, sharing one generated vector across every assignment scenario
but using an independent response-uniform stream. This answers a conditional
rank-preserving population-realization question; it is not an unconditional
new-person sample or a generated assignment mechanism. The latter remains an
M3 repository-only research requirement, together with tie handling,
covariate-stratified latent-regression generation, and frozen operating
characteristics before any broader claim.

The first M3 repeated-simulation scaffold is now retained under
`validation/operating_characteristics_20260809/` and governed by
`validation/OPERATING_CHARACTERISTICS_PROTOCOL_20260809.md`. It deterministically
pairs null and +0.60-logit local-bias conditions, adds sparse/missing and paired
correct/contaminated-anchor controls, separates attempted/converged/eligible
denominators, and reports Monte Carlo uncertainty, recovery, conditional-Wald
coverage, and raw/display conclusion sensitivity. The retained two-replicate
Python JMLE run verifies the orchestration and negative controls only. M3 is
still open: study-depth precision, remaining decision surfaces, CMLE
Person/anchor/local-bias scope, and multi-
platform evidence are required before any public simulation-evidence button is
enabled.
The cross-engine input gate now retains exact generated ratings/truth/anchors,
per-run data and fit-input identities, and byte-level file hashes. Its R bridge
passed 16/16 import checks on the smoke bundle and records the isolated mfrmr
0.2.3 source snapshot plus TAM/immer/sirt function identities. This closes the
RNG/seed ambiguity but not the pending estimator-adapter or agreement gates.
The mfrmr JMLE adapter is now complete for the smoke matrix. Its strict mode
was inference-ready for 12 identified RunIds, enforced all hard anchors
exactly, and agreed closely with Python on jointly included facet and focal
bias estimates; its matched-control mode demonstrated that optimizer code 0
can coexist with unacceptable terminal gradients. mfrmr rejected four sparse
rank-deficient designs before optimization, whereas Python returned numerical
fits and withheld only the sparse focal decision. This difference is now a
named readiness-policy gap.

The matched unanchored Python-CMLE/`immer_cml` adapter is also complete for
the smoke matrix. Both engines returned all 12 structurally eligible RunIds;
96 free coordinates differed by at most 0.0000000563 logits and conditional
log likelihoods by at most 0.00000000000477. Four sparse RunIds failed the
independently reproduced Python/R conditional-information rank gate (5/12,
nullity 7), with rank/nullity agreement on 16/16 RunIds. This strong
implementation agreement also exposed a consequential stopping-rule boundary:
all 12 `immer` fits reported optimizer code 0, but terminal-gradient readiness
was 12/12 at `1e-4`, 7/12 at the primary `1e-5`, and 0/12 at `1e-6`.
Anchor-requesting runs remain unanchored parity references, not anchor
validation. The +0.60 local-interaction generator shifted balanced-small
additive coordinates by 0.0612 logits on average and at most 0.172, retained as
omitted-interaction leakage rather than a CMLE bias decision. CMLE local-bias
estimation and Person scoring remain unsupported.

The sirt adapter is now complete as a deliberately non-parity MML sensitivity
lane. It maps Task x Criterion to virtual PCM items, estimates a common normal
Person population, and compares 30- versus 61-point quadrature without pooling
raw and mean-centered item locations. Thirty of 32 returned fits passed the
analysis gate; the same sparse +0.60 RunId hit the 1,200-iteration cap under
both grids and remains visible but excluded. Balanced-small Python-JMLE/sirt
Rater differences were small descriptively, while sparse differences exposed
the distinct assumption-based MML identification route. Exact CMLE and strict
mfrmr rejected all four sparse designs structurally, whereas primary sirt was
eligible for three. Paired contaminated anchors transmitted the injected
+0.25-logit location shift essentially unchanged. These findings do not imply
package parity, anchor robustness, or performance. Study-depth precision,
PCM operating characteristics, CMLE anchor/Person/local-bias scope, and
multi-platform evidence remain open gates.

The TAM adapter is also complete as an additive-response-surface MML
sensitivity lane. It preserves Criterion + Rater + Task + common-step RSM
structure while integrating an estimated normal Person population. All 32
q21/q61 fits returned and passed the retained loop/progress convergence audit;
30 were eligible. Two q21 sparse fits stored variance as `0.0010000001` at the
configured `0.001` lower boundary and were excluded through a `1e-8` numerical
band. All 16 primary q61 fits remained eligible. This boundary would be missed
by naive exact comparison despite identical rounded display.

Balanced-small Python-JMLE/TAM-MML Rater MAE was 0.01383 logits, and TAM/sirt
MML Rater MAE was 0.001513. Sparse primary TAM was eligible for four RunIds
where exact CMLE and strict mfrmr rejected all four structurally. Quadrature
sensitivity reached 0.02224 logits for Rater, 0.3814 for Person EAP, and 0.6990
for cumulative response-surface difficulty. Under contaminated anchors, TAM's
sum-zero constraint shifted fixed Raters by +0.25 and redistributed about
-0.25 to unanchored Raters. Derived last-level SE coverage and cross-estimator
information-criterion comparisons remain withheld. The first-read projection
is retained but no public TAM/cross-engine button is enabled.

The prospective M3 precision plan and 20-replicate Python pilot are now also
retained. Exact manifest-extension auditing preserved the 16 smoke rows inside
the 160-row pilot. Runtime and failure accounting scaled, but the pilot blocks
study expansion: 156/160 optimizer flags converged and 124/160 app decisions
were eligible, while 0/160 passed the separately labelled `1e-4` terminal-
gradient sup-norm sensitivity. Sparse null/alternative conditions each yielded
4/20 eligible decisions and exceed the registered 2,500-attempt planning cap.
The row-level floating audit found five raw/display fit-label disagreements,
and paired anchor contamination transmitted +0.25 exactly to fixed Raters.
The prospectively frozen strict-Python-JMLE follow-up has now separated two
questions that the pilot had conflated. The selected L-BFGS-B precision polish
passed every frozen numerical gate in 160/160 runs at raw gradient sup norm
`<=1e-4` (126/160 at `1e-5`; 43/160 at `1e-6`) without a raw-versus-displayed
gate disagreement. However, a post-fit movement audit found flat-direction
movement up to 57.02 logits, and the exact eta free-coordinate design audit
showed nullity 7 with eight Person-Rater components in every one of the 40
sparse runs. Thus numerical qualification did not authorize unconditional
JMLE integration.

The application now performs that eta-rank audit before JMLE optimization and
separates `Converged` from `InferenceReady`; rank-deficient fits withhold the
conditional bias result while retaining machine-readable rank, nullity,
connectivity, and null-space evidence. The integrated guard matched the frozen
audit contract in 160/160 retained RunIds. It does not change the optimizer or
estimates and never auto-switches to MML. M3 next requires a connected-design
extreme-score policy and downstream precision-result reconstruction for the
rank-full conditions, plus prospective redesign, anchoring, or an explicitly
qualified alternative estimand for the sparse conditions. The guard covers
only the Person/facet eta block; threshold/slope identification remains under
the separate engine gates. The frozen v1 continuous-outcome budget conflict
also requires a prospective amendment; no public evidence button is enabled.

A second output-specific gate now covers JMLE Person score boundaries. Across
all 8,320 retained Person/run rows it reproduced the frozen all-minimum/all-
maximum classification exactly. The 120 structurally identified runs contained
four boundary Persons, and every >=1-logit movement in that rank-full subset was
confined to those Person coordinates; maximum non-theta movement was 0.000293
logits. The app retains the optimizer/constraint `Estimate` only as technical
evidence, sets `ReportableEstimate` unavailable, and marks Person inference
unready at the boundary. Classification uses integer score patterns rather
than a parameter-magnitude cutoff, so floating-point display cannot change it.
No finite correction, estimator switch, or precision-polish core integration
has been selected. A repository-only fixed-calibration Warm WLE comparison has
now completed the first prospective correction experiment for rank-full
designs. Python matched TAM on 16/16 frozen RSM/PCM/GPCM Person fixtures within
`1.68e-11` for theta and `1.18e-11` for SE. Replaying the 120 rank-full Stage-B2
runs compared 7,600 Persons while holding facets and steps fixed. The four
exact-extreme patterns moved by 23.72--25.19 logits; median absolute interior
movement was about 0.009 logits. Person-fit zones changed 51 times, and `.3g`
display rounding disagreed with the raw fit decision 54 times across JMLE and
WLE surfaces. The focal Holm and combined strong-bias decisions were stable,
but one practical `|bias| >= 0.50` decision crossed the threshold. The 40
rank-deficient sparse runs remained withheld. WLE therefore remains an
explicit alternative Person estimand, not an automatic JMLE repair. The next
gate is a user-facing estimator-policy design plus calibration-uncertainty and
bootstrap sensitivity; precision polish, WLE, JMLE boundary output, and MML
must not be silently pooled behind one result label.
The first repository-only exact-CMLE-to-WLE fit-sample bridge now also passes
its frozen RSM/PCM contracts: category surfaces, direct generic WLE scores,
and SEs matched exactly for 14 Persons including two exact extremes; parameter
ordering was invariant and missing identity failed closed. This establishes
the candidate two-stage architecture, but not public readiness. Its first
asymptotic calibration-sensitivity gate also passed: 2,000 covariance draws per
RSM/PCM model returned all 32,000 draw-Person scores, with model-median
calibration-draw SD `0.0616`--`0.0801` logits, a maximum `0.518` logits, and a
maximum calibration-SD/conditional-SE ratio `0.427`. Exact extremes were most
sensitive. The retained quantiles are explicitly not confidence intervals and
the quadrature sum is not a total inferential SE because same-sample CMLE/WLE
dependence is omitted. A six-card first-read projection keeps research-ready,
caution, and withheld states distinct. Conditional bootstrap qualification,
new-Person and unseen-unit design construction, and the public Streamlit
estimator-policy surface remain prospective gates. A later five-state
repository contract now implements the six-card one-click view model, but it
does not establish user comprehension or promote the UI. Native hard facet
anchors are now implemented as fixed likelihood offsets and passed a
prospectively amended 32-fit engineering smoke. That evidence does not promote
the UI.
The bootstrap design and machine-readable plan were frozen before
implementation. The initial repository-only core now implements fixed-score
conditional response-pattern sampling, joint plug-in response generation,
fail-closed refitting, and complete attempt ledgers. Six unit tests pass,
including exact enumeration and RSM/PCM total preservation. The frozen
200-replicate RSM/PCM matrix completed 800/800 full-rank, inference-ready,
WLE-available refits; the Wilson 95% lower bound is nevertheless `0.981` per
cell. Fixed-score totals never changed. Median Person bootstrap SD was `0.0602`
logits under fixed-score pattern resampling and `0.615` under joint plug-in
resampling; maximum SD was `0.744`, and 441 Person-replicates changed exact-
extreme status in the joint lane. Sampler-contract status passed, but promotion
is withheld. The variance difference between lanes is not a decomposition and
a separate post-pilot repeated-truth coverage protocol is still required.

The next prospectively registered downstream-diagnostic gate is complete at a
repository research level. The untrimmed fixed-calibration CMLE-WLE Person
Infit/Outfit kernel matched `sirt::pcm.fit` for 32 RSM/PCM Person rows and 162
administered observations, including four exact-extreme Persons and planned
missingness; maximum Infit and Outfit differences were `4.89e-15` and
`4.00e-15`. TAM's default Person Outfit remains a separately labelled estimand
because its JML fit code trims unusually large squared-standardized residuals
before aggregation. The separate 800-replicate fit extension returned Person
fit for every attempt and preserved exact persisted WLE identity with the
original pilot. It retained 2,372 Infit transitions, 2,454 Outfit transitions,
24 display-rounding boundary statistics, and nine raw/display classification
mismatches; classifications always used unrounded MnSq. The first extension
run failed the unchanged exact-zero identity gate because it compared a parsed
CSV with an unpersisted in-memory float. That failed evidence remains
immutable; a registered comparison-method amendment required a full corrected
rerun using symmetric `%.17g` persisted readback, which passed with zero WLE
difference and without changing the original or failed-run inventories.

The subsequent post-result threshold audit confirms that rounding did not
generate the retained result: raw canonical classes/transitions reproduced
exactly, while a counterfactual three-decimal classification disagreed on nine
replicate-statistic classes and five either-statistic transition indicators;
six decimals reproduced raw decisions exactly. Dependence on the substantive
threshold triplet is much larger: over a frozen 125-triplet grid, aggregate
either-statistic transitions ranged from 1,975 to 3,675. A transparently
registered class-decomposition addendum exposed another limitation of the
binary transition projection. Varying only the noisy/distorting boundary from
1.90 to 2.10 left Infit/Outfit/either transition totals fixed at
2,372/2,454/2,554, yet moved 37 Infit and 51 Outfit replicate labels from
`distorting` to `noisy`. No threshold or ordinal severity score was selected.

A prospectively registered fixed-calibration known-truth fixed-calibration pilot now
separates rule behavior from calibration estimation. All 160 runs returned and
retained 624,000 responses. At the canonical raw either-upper rule, clean rates
were 3.4%--3.7% for 24 observations and 10.1%--10.8% for six observations;
affected detection ranged from 3.8% to 46.6% across frozen mechanisms. The
post-result registered recovery addendum retained 163 exact-extreme WLE rows,
found non-extreme RMSE of 0.263--0.720 logits, and found 2,333/32,000 WLE versus
generating-theta flag disagreements. It also corrects the nominal 125-grid
interpretation: the four rules have only 5/5/25/5 effective threshold input
configurations. Twenty replicates remain design debugging, not validation.

This closes the formula-translation and descriptive bootstrap-transition
tasks, not the downstream-diagnostic qualification milestone. ZSTD, p-values,
threshold optimality, confirmatory repeated-truth precision, estimated-CMLE
calibration uncertainty, sparse connectivity beyond the initial randomized-
anchor topology stress, drift distributions, differential/common-anchor provenance,
matched repeated mfrmr/TAM/immer/sirt fits, new-Person/unseen-unit
prediction, saved-result schema migration/public disclosure, and confirmatory
human comprehension remain open. Current-schema private archive identity and a
non-public bilingual comprehension instrument/scorer are complete; structural
language parity is not human comprehension, translation equivalence, or
accessibility evidence. A balanced private cognitive-interview schedule,
answer-free packets, moderator guide, and fail-closed blank record schema are
also complete, but provide neither ethics approval nor recruitment authority.
Frozen instrument content identity and directional confirmatory-gate mechanics
are also complete. They separate any comprehension error from an unsafe error
direction, retain invalid slot denominators, and prohibit pooling changed
instrument content. They do not register a confirmatory sample size or provide
human evidence. No public CMLE, WLE, or bootstrap fit button is enabled.

The native hard-anchor hard-anchor smoke reused 24,000 fixed response rows for 32
RSM fits spanning dense/sparse, clean/random, zero/one/three correct, and three
uniformly contaminated anchors. All fits returned inference-ready and all
engineering gates passed. A post-result registered diagnostic then confirmed
that a common +0.25 Rater-anchor shift is an origin transformation: WLE moves
by +0.25, while conditional likelihood, probabilities, MnSq, and raw fit flags
remain invariant within numerical tolerance. Consequently, later UI work must
warn that fit cannot validate anchor origin and must expose anchor provenance,
drift, and common-versus-differential contamination.

The larger prospectively registered differential milestone is now complete:
320/320 fits returned inference-ready over 120,000 retained response rows, four
dense/sparse clean/random conditions, and eight anchor scenarios. It found no
universal anchor share or monotonic recovery benefit. A zero-mean differential
error changed paired raw Person-fit decisions by 2.4%--15.8% depending on
condition/group, while the common-origin error changed none; 3-decimal
display-driven MnSq decisions also disagreed with raw decisions 12 times among
64,000 Person-scenario rows. Full rank and acceptable gradients therefore do
not certify anchor validity.

The prospectively registered randomized-anchor-content and bridge-topology
milestone is now complete. Its graph-component rule matched exact prefit rank
in 350/350 cases, but the minimally connected chain showed why rank alone is
not enough: 287 cases passed prefit, 286 returned, and only 274 were inference-
ready; all 13 shortfalls occurred with one Person per edge. Partial anchors did
not create missing contrasts, and two anchors covered the two-component design
in only 7/10 randomized assignments.

The prospectively registered finite-domain remediation is now complete. It
changed the one non-finite exception into a typed returned-but-nonready result,
raising returnability to 287/287 while leaving readiness at 274/287. Baseline
ready states and 82,200 Person raw/rounded decisions were unchanged. This does
not solve likelihood existence: a post-remediation binary perfect-separation
check shows that small gradients and full finite-point rank can coexist with a
maximum at infinity.

The registered conditional-separation/finite-MLE research milestone is now
complete. Exact exhaustive support enumeration matched all 93 frozen controls
but exposed its intended scalability limit. The separately frozen fixed-score
support oracle then matched 1,184/1,184 exhaustive maxima, 93/93 fixture
statuses, 123/123 completed retained statuses, and six high-cap sentinels. It
resolved the same-byte 350-case stress as 63 structural, 13 boundary, and 274
interior cases with no unavailable/tolerance-unstable result. All 274 current
ready fits were interior. The eligible-case P95/max was 0.521/0.568 seconds;
the 350-case total passed the prospectively fixed 120-second gate at 109.216
seconds. This is an internal integration candidate, not a button. The next
prospective milestone is also complete: default repository-only `fit_cmle`
readiness now requires an oracle interior. All 287 same-byte fits returned;
274 interior cases remained ready and 13 boundary cases were non-ready with a
typed reason, while estimates, likelihood, WLE, MnSq, and raw/rounded decisions
remained numerically invariant. The next milestone is a structured pre-
optimizer boundary result, followed by drift distributions, matched
mfrmr/TAM/immer/sirt boundary replay, and comprehension-tested staged
explanations.

The structured pre-optimizer milestone is now complete. A stable four-stage
Japanese/English result matched 93/93 registered fixtures and 350/350 retained
cases. It stopped 63 structural and 13 boundary cases before optimization and
optimized only the 274 interior cases, with zero forbidden optimizer calls and
binary64-level numerical invariance. P95 time was 0.068 seconds for a design
block, 0.306 for a boundary stop, and 1.084 for a ready interior path. Person
scoring remains an explicit available/not-run stage, and message presence is
not comprehension evidence.

The first matched cross-engine boundary milestone is also complete for 13
unanchored additive-RSM controls. Python and matched immer CMLE agreed on five
oracle-interior cases within `1.78e-9` for free coordinates and `5.33e-15` for
conditional log likelihood. Yet on seven oracle-boundary cases, immer returned
code zero in 3/7 and finite coefficient/SE vectors in 4/7 at the primary cap;
the frozen mfrmr 0.2.3 JML snapshot reported convergence in 7/7. These are not
contradictions: the support oracle answers finite conditional-optimum
existence, while mfrmr JML, TAM MML, and sirt MML have different Person
treatment and estimands. The public workflow must therefore gate before
optimization and cannot promote finite output or convergence text to a green
CMLE result.

The separately registered PCM and hard-anchor cross-engine milestone is now
also complete. Python exact PCM CMLE and matched immer `W`/`b_const` coordinates
agreed across all six oracle-interior cases within `2.45e-9`, with conditional
log-likelihood difference at most `3.55e-14`; all six separation/unused-support
cases remained oracle-blocked at every retained iteration cap. This covers
Rater and Criterion main-effect anchors, contaminated anchors, and a
differential two-facet anchor, but not step anchors. mfrmr JML and TAM/sirt MML
remain descriptive sensitivity lanes, and exact enforcement cannot validate
anchor content.

The downstream scoring/rounding card milestone is now complete at the
Streamlit-free contract level. Five registered workflow states returned the
same six bilingual cards. Only two ready calibrations reached fixed-calibration
WLE and Person MnSq; three input/design/existence stops made zero downstream
attempts. Fourteen ready Person handoffs reproduced the direct scoring path
exactly, and the registered `1.5004` Infit probe remained raw-classified as
`noisy` despite a displayed `1.500`. This proves orchestration and raw-value
control, not user comprehension. The next instrument-readiness increment is now
complete: 10 non-public bilingual previews, 100 participant tasks, a separate
100-row key, six non-compensatory dangerous-misconception domains, and 30
synthetic fail-closed scoring packets passed their automated contract. Human
participants remain zero, so cognitive interviews, an instrument-frozen pilot,
and a freshly registered item-by-language confirmatory gate remain mandatory.
The private interview operations kit now balances all case pairs, positions,
and target experience strata and rejects PII columns and record drift; it does
not authorize or complete those interviews. High-dimensional sparse PCM, a
clean released-mfrmr reproduction, private-archive schema migration, real-
browser/assistive-technology QA, and an independently governed public
disclosure contract also remain before any public button.

The confirmatory arithmetic increment is now complete without human data. Six
directional danger domains are assessed separately by language, giving 12
primary cells. A zero-error cell at `n=24` misses the strict one-sided 95%
Wilson upper-bound rule (`0.101310`), while `n=25` passes it narrowly
(`0.097654`) but fails the 12-cell Bonferroni sensitivity (`0.217782`). These
values expose the pilot/confirmatory distinction; they do not select 25 as an
adequate sample size. A future protocol must prospectively justify minimum n,
anticipated error, power/precision, attrition, clustering, and accessibility
before any responses.

The registered no-human-data sensitivity surface now sharpens that requirement.
Exact binomial pass probability under the Wilson rule is sawtoothed: with true
dangerous-error probability `0.05`, the 90% single-cell first crossing is
`n=224`, whereas the first value sustained through the registered `n=5000`
search is `n=260`. Conservative 12-cell union-bound projections are `n=422`
and `n=456`. Independent attrition assurance gives 231 planned slots for 200
valid responses at 90% homogeneous retention and 95% assurance; a separate
design-effect heuristic gives nominal `n=145` for effective `n=100`, mean
cluster size 10, and ICC 0.05. At 75 slots per language, each of three blocked
danger mechanisms is exposed only 25 times. These diagnostics deliberately do
not collapse into a selected n. A future human-study protocol must choose its
estimand, familywise rule, anticipated rates, clustering model, retention
assumptions, mechanism-specific minima, and accessibility strata prospectively.

An explicit dependence/MNAR stress now demonstrates why those choices cannot be
reduced to a design-effect multiplier. In 64,000 registered synthetic
replicates, participant dependence and clustering changed the all-12 gate even
when marginal truth stayed fixed. More seriously, at 500 planned slots per
language, dangerous-response under-retention transformed complete-data truth
near 0.10 into an observed valid-record rate of `0.060315` and yielded a
`0.4035` false-reassuring pooled-gate rate. A rotated blocked-mechanism hotspot
yielded `0.4955`; the combined mechanism/dependence/cluster/MNAR stress yielded
`0.9675`. Increasing n can therefore make a biased or incorrectly pooled
conclusion more precise. Before human recruitment, the protocol must define
invalidity estimands and denominator accounting, collect reasons for
invalidity without PII, govern mechanism-specific stopping/protection, and
predefine cluster-aware analysis. These simulations do not estimate realistic
correlations, missingness, language effects, or user behavior and do not select
a sample size.

The private protocol-preflight increment is now complete at the repository
contract level. It converts the remaining work into 13 explicit blocking
decisions and binds decision content, invalidity codes, denominator schema, and
partial-identification sensitivity to one protocol SHA-256. Every scheduled
slot must remain in a categorical denominator ledger; direct identifiers,
free-text substitutions, missing/duplicate slots, unregistered reasons, and
post-outcome exclusions fail closed. At valid n=500 and observed 5% error, an
invalid/valid ratio of 5% yields an all-invalid-dangerous Wilson upper bound of
`0.118434`, despite the observed-valid analysis appearing acceptable. The
deterministic private ZIP therefore leads with `BLOCKED`, not a green research
status. This prepares evidence structure but does not resolve the estimand,
sample size, multiplicity, cluster model, mechanism protection, language or
accessibility objectives, ethics, consent, privacy, retention, or incident
response. Those decisions require prospective external ownership before any
recruitment manifest can be frozen.

The private decision-workbench increment now makes those 13 decisions
reviewable without allowing the software to choose them. Each has three
non-ranked options, yielding 39 rows with explicit strengths, risks, evidence,
and qualitative impacts. The registered dependency graph has 32 edges and no
cycles; sample-size decision SCI-02 is downstream of seven scientific choices,
so n cannot be entered first and used to back-justify the estimand or analysis.
The read-only bilingual HTML is network-free and contains no form controls.
Selection references remain unverified even when a synthetic worksheet is
complete, so the workbench cannot authorize recruitment. The next genuine
milestone is external deliberation and prospective evidence for SCI-01 and its
downstream scientific decisions—not another numerical optimization over the
existing simulations.

The next SCI-01 support increment is complete at the private software-contract
level. Three formal target-population definitions are projected over the same
72 partial-identification scenarios. Thirty-six scenarios demonstrate that an
observed-valid pass can coexist with an all-invalid-dangerous proxy failure and
a non-robust dual conclusion, without any displayed-rounding involvement. This
exercise exposed a material data-contract gap: `V+I` cannot be silently
substituted for scheduled slots, and invalid records cannot all be called
composite failures without a prospective code map. Consequently Option B is
shown only as a diagnostic proxy. The repository still does not select SCI-01.
The next genuine milestone remains external scientific deliberation: define
the scheduled population, adjudicate invalidity codes with participant-access
and withdrawal implications, choose the allowed claim, and attest the choice
before outcome inspection. Only then can SCI-04/05 and later SCI-02 sample-size
work proceed. Static offline HTML contract tests are complete; public button
workflow, keyboard/screen-reader/browser acceptance, and comprehension remain
withheld.

The invalidity-adjudication support layer now makes the next dependency
explicit without resolving it. Fifteen reason codes are separated into
consent, withdrawal, eligibility, technology, record/response integrity,
missingness, administration, accessibility, protocol, and amendment families.
Only `none` is fixed; 14 rows require external review. In the registered
diagnostic grid, 36/72 source scenarios cross the decision boundary as the
fraction of invalid outcomes classified as composite events varies. The
opposite movement of the zero-fraction and all-invalid-dangerous curves shows
why the fraction cannot be optimized after results. The next milestone is a
prospective, owner-attested code map and approved withdrawal/data-use policy,
followed by SCI-05 missingness/MNAR rules. No sample-size work should consume a
chosen composite risk until those dependencies are frozen.

Product sequencing is now explicitly UI-first before another private research
workflow is added. The first setup increment moves the Guided-defaults Run
action from the technical sidebar into the main workspace after mapping and
readiness. It uses goal-oriented bilingual wording and preserves the expert
sidebar route only for Advanced controls. This satisfies the single-primary-
action rule for the setup phase at the Streamlit contract/AppTest level. It
does not yet satisfy real-browser keyboard, focus, narrow-screen, screen-reader,
or first-time comprehension acceptance; those remain the governing next gate.
The second increment replaces the ten-choice source list with four task-level
source classes and an always-visible example-scenario selector. Legacy source
IDs, restored workspaces, immediate privacy severity, and last-sample memory
are preserved. The third increment removes technical controls from the Guided
setup surface and restores explicit standard values when switching back from
Advanced controls. Weight columns are not silently applied or offered as
facets. The remaining work is rendered-browser keyboard, focus, narrow-screen,
screen-reader, and first-time comprehension acceptance, followed by
privacy-safe task telemetry if the product decision permits it.

### M1 — one authoritative workflow shell (0-30 days)

Outcome: the default journey exposes one primary next action per phase while
retaining every scientific caveat and expert detail.

Deliverables:

- define one widget-independent workflow-shell record for source, setup, run,
  first-read, evidence detail, and archive phases;
- remove remaining duplicate orientation surfaces from the default result
  route while keeping stable Help targets and compatibility projections;
- keep the authoritative result navigation visible during long-page scrolling,
  with keyboard focus, Help return, reduced motion, and narrow-viewport browser
  acceptance;
- measure initial and post-run render topology so new visible controls require
  an explicit information-architecture decision; and
- add a clean-checkout CI job that runs compile, doctor, release-check,
  self-test, native AppTests, and privacy/localization gates.

Exit gate:

- one primary action is visually exposed in every default phase;
- the authoritative result selector remains operable while a long selected
  section scrolls, without covering required content on a narrow viewport;
- the sample journey and real-data setup are keyboard-completable without a
  hidden required control;
- view/language/help transitions trigger no estimation and preserve
  AnalysisID; and
- no raw response row is visible by default.

Progress on 2026-08-04:

- introduced a widget-independent six-phase workflow-shell contract with one
  named primary-action owner per phase;
- removed the duplicate top-level first-read overview from the default result
  route while retaining its full detail in the First read section;
- captured initial and post-run AppTest control budgets so additions to the
  rendered topology require an explicit information-architecture change; and
- converted the Python workflow into an explicit clean-checkout matrix with
  named workflow/privacy/localization contracts and a final worktree check.

Navigation refinement on 2026-08-04:

- removed the sidebar keyboard-shortcut cheat sheet and its maintenance/cache
  commands; the app defines no custom hotkeys;
- placed Essential and All-panels result selectors in one sticky navigation
  dock that remains available during long-page scrolling; and
- kept the narrow layout compact with a touch-scrollable one-line section
  switcher instead of a multi-row fixed overlay.

First-run route refinement on 2026-08-04:

- added the pure `mfrm_app.guidance` five-node catalog and exhaustive reducer;
  it imports neither Streamlit nor pandas and has no claim-readiness field or
  estimator callback;
- replaced the initial configuration-heavy page with one optional three-route
  decision: learn with a sample, start with user data, or continue without the
  guide;
- connected data-role review, production-path estimation, a bounded formative
  interpretation, and a non-scientific learning checkpoint into one
  skippable/resumable/restartable sample route; and
- kept guide completion independent from AnalysisID and preserved the fitted
  identity across formative answers and completion.

The browser accessibility matrix is now executable and fail-closed through
`mfrm_app.accessibility_acceptance`: eight tasks, two locales, seven targeted
profiles, 52 cases, and 482 required evidence rows. Static/AppTest evidence
cannot promote it. The generated blank bundle is `NOT_ACCEPTED`, stale catalog
fingerprints are rejected, and the evaluator exits nonzero until every required
browser record is `PASS`. Executing the matrix—including actual screen-reader
observation—remains the next M1 acceptance item.

Increment verification:

- `make apptest`: 1,117 passed, 62 compatibility tests deselected;
- complete five-step sample AppTest: production fit, overclaim retry,
  learning-only completion, and stable AnalysisID passed;
- guide state/isolation AppTests: focused landing, own-data route, Japanese
  rerender, exit/resume, and existing-workspace restoration passed;
- representative initial/post-run AppTests: two scenarios passed;
- All-panels Help return: panel state and AnalysisID preserved without refit;
- `--self-test`: 56 passed;
- `--doctor` and `--release-check`: passed; and
- benchmark and demo-report generation smokes: passed.

### M2 — scalable input and modular presentation (31-90 days)

Outcome: data entry and result navigation can grow without extending the
70,000-line entrypoint or lengthening one global control list.

Deliverables:

- replace the flat source list with source class -> scenario/source detail,
  preserving stable IDs and migrating existing session/help state;
- extract source ingestion, setup/readiness, workflow shell, and result-router
  adapters into focused `mfrm_app/` modules;
- keep computation, evidence, and presentation policies independently tested;
- add rendered-route locale coverage for Japanese and English, including
  errors and advanced settings; and
- establish privacy-safe UX telemetry schemas for phase reached, stable reason
  code, rerun count, and time to first evidence, without input values or IDs.

Exit gate:

- `streamlit_app.py` is smaller after extraction and owns no new pure business
  rule;
- adding a sample scenario changes a registry, not the global page topology;
- old stable source IDs and Help returns continue to work; and
- no telemetry event can contain uploaded values, labels, paths, free text, or
  fingerprints.

### M3 — coherent evidence and design resilience (3-6 months)

Outcome: the app explains tensions between diagnostics and identifies fragile
rating designs without turning either into an automatic verdict.

Deliverables:

- normalize fit, dimensionality, local dependence, category, bias,
  reliability, and design evidence into one coherence view;
- implement deterministic leave-one-level and concentrated-missingness design
  perturbations with stable reason codes;
- provide constrained repair candidates with burden and non-claim boundaries;
- generate a human-readable Analysis Brief plus deterministic JSON sidecar;
  and
- introduce a conclusion-level sensitivity ledger that retains failed and
  unsupported specifications in its denominator.

Exit gate:

- contradictory fixtures produce an explicit tension state;
- disconnected, bridge-dependent, sparse, and imbalanced fixtures reproduce
  known behavior;
- every recommended repair states what it addresses and what it cannot prove;
  and
- no report claim cites stale, unavailable, or differently identified
  evidence.

### M4 — calibrated native validation and beta decision (6-12 months)

Outcome: public claims are backed by frozen Python-native validation evidence,
and the project can make an explicit beta-exit decision.

Deliverables:

- freeze ADEMP-style RSM, PCM, and bounded-GPCM scenario contracts;
- report bias, RMSE, interval coverage, recovery, flag rates, numerical
  failures, and Monte Carlo uncertainty for applicable conditions;
- establish checked runtime and peak-memory baselines for supported Python and
  dependency versions;
- complete the retain/deprecate/remove decision for dormant compatibility
  generators, viewers, fixtures, and locale copy; and
- publish a machine-readable validation and release-scope summary used by the
  UI, documentation, and release gate.

Exit gate:

- every advertised diagnostic or estimator claim links to a named validation
  condition or is explicitly labelled uncalibrated;
- clean-environment release gates pass on the supported version matrix;
- compatibility-only code cannot re-enter the public journey or release gate;
  and
- the public-beta label is retained or removed by a documented evidence-based
  decision, never by schedule alone.

### Operating scorecard

Review these measures at each monthly roadmap checkpoint:

| Dimension | Measure | Direction |
|---|---|---|
| Journey clarity | Visible primary actions per phase | Exactly 1 |
| Safety | Raw response rows visible by default | 0 |
| State integrity | Presentation-only transitions that refit or change AnalysisID | 0 |
| Evidence integrity | User-facing claims without current EvidenceIDs and boundaries | 0 |
| Accessibility | Sample and real-data setup tasks passing the browser matrix | Increase to 100% |
| Architecture | New pure rules added to `streamlit_app.py` | 0 |
| Reliability | Native test, doctor, self-test, and release-check failures | 0 |
| Performance | Material regressions without reviewed benchmark evidence | 0 |

Items outside the current milestone may proceed only when they do not delay its
exit gate and do not create another workflow, evidence, or state authority.

## Product boundaries

### In scope

- Native Python RSM, PCM, and bounded GPCM workflows within their documented
  identification and estimation limits.
- Native JMLE and MML estimation, uncertainty audits, prediction, simulation,
  and parameter-recovery evaluation.
- Planned optional JAX acceleration of the same JMLE/MML likelihoods, subject
  to numerical equivalence and cold/warm end-to-end performance checks.
- Planned NumPyro Bayesian RSM/PCM research after the acceleration evaluation,
  governed by separate estimand, prior, diagnostic, and release gates.
- Scoped offline ConQuest/TAM exports, including the released restricted PCM
  runner; verified return ingestion remains planned. The standalone default
  runtime remains independent of external installations.
- Repository-only native exact CMLE research for eligible RSM/PCM structures;
  public workflow exposure remains governed by its separate integration gate.
- Prospective rating-design diagnostics and post-collection design audits.
- Facet-aware design, severity/leniency, halo, local-dependence, and residual
  network diagnostics.
- Residual PCA, PCA stability, DIMTEST, model-fit, category, bias/interaction,
  covariance, and reliability evidence.
- Sensitivity analysis that records when a substantive conclusion changes.
- Guided and research-oriented interfaces backed by the same statistical
  contracts.
- A structured study-context, rating-design, and scoring-process record that
  distinguishes what the application computed from what a study intends to
  infer or do.
- Reproducible Python-native reports, tables, figures, settings, fingerprints,
  and claim boundaries.

### Explicitly out of scope

- EGA or automatic discovery of a confirmatory Q matrix.
- Automatic calls from the hosted/default product journey to TAM, Shiny,
  FACETS, ConQuest, `mfrmr`, `Rscript`, `rpy2`, or other external estimators.
- Automatic external-job submission or unvalidated external-result ingestion.
  The released optional PCM helper runs locally within its explicit scope;
  a general external-job service is not planned.
- Unqualified cross-engine parameter comparisons or posterior ingestion that
  bypass identity, privacy, numerical, or estimator-specific validation.
- A remote R worker or a hosted external-solver queue.
- Confirmatory multidimensional MFRM unless it is later implemented and
  validated as a native Python model through a separate scope decision.
- Claims that a single diagnostic proves model truth, unidimensionality,
  fairness, rater quality, or suitability for a high-stakes decision.

The 2026-09-11 extension decision supersedes the blanket prohibition on
ConQuest/TAM interchange. Selected archived export and comparison contracts
will be reused through the extension plan, while other legacy surfaces retain
their inventory and retain/deprecate/remove decision. External installations
must not become dependencies of the core workflow. External comparison
evidence may qualify its specific adapter or numerical scope, without making
ConQuest or R/TAM availability a requirement for the standalone release checks.

Python-native report archives may remain as application-specific download
packaging. They are not cross-engine interchange formats, and the application
will not promise that another engine can execute or return them.

## Starting point

This roadmap extends existing Python work rather than proposing a replacement
application.

| Area | Current foundation | Principal gap addressed here |
|---|---|---|
| Estimation | Native RSM/PCM/bounded-GPCM JMLE and MML paths; repository-only exact RSM/PCM CMLE core | Conclusion-level sensitivity, conditional-design qualification, and public CMLE integration gate |
| Dimensionality screens | Residual PCA, leave-one-column/bootstrap stability, DIMTEST | Coherent interpretation when the screens disagree |
| Rating networks | Design, severity/leniency, and halo networks | Perturbation resilience, effective-N evidence, and repair candidates |
| Uncertainty | Covariance/SE audits and prior-SD sensitivity | One versioned sensitivity ledger across supported decisions |
| Simulation | Synthetic data and parameter-recovery utilities | Immutable scenario contracts, formal runner, failure accounting, and MCSE |
| Evidence and reports | Readiness, claim-to-evidence, publication, and visual outputs | One normalized evidence source for UI, prose, and exports |
| Architecture | A growing `mfrm_app/` helper package plus a large entrypoint | Tested computation modules with compatibility wrappers |
| UX | Guided essential view, research detail, bilingual scaffolding | A shorter decision journey and evidence-linked actions |

The adversarial baseline, delivered first corrections, and task-level UX
acceptance measures are recorded in
[`docs/ux_adversarial_audit.md`](docs/ux_adversarial_audit.md). Future UX work
should reduce competing orientation surfaces against that baseline rather than
adding another local guide.

The current network functions are therefore modularized and extended, not
reimplemented. Likewise, existing PCA, DIMTEST, sensitivity, and reporting
outputs become inputs to a common evidence contract rather than parallel new
versions.

## North-star outcomes

### Defensible design

Before collecting data, a user can compare candidate assignments and identify
disconnection, single-link dependence, overloaded raters, weak overlap, and
likely failure under missing ratings or rater withdrawal.

### Explainable evidence

After estimation, every surfaced conclusion can be traced to a named evidence
record with its scope, applicability, warning state, sensitivity state, and
recommended next action.

### Robust conclusions

The application reports whether important conclusions persist across
prespecified estimator, scale, missingness, category, and leave-one-level-out
checks. It does not silently convert a fragile point result into a firm claim.

### Reproducible decisions

The application preserves the data fingerprint, resolved settings, seeds,
numerical diagnostics, evidence state, claim limits, and user-visible action
record needed to reproduce the analysis in the same Python application.

## Evidence vocabulary

Scientific status uses two different state systems, which will be kept
separate from workflow and learning progress.

### Computation state

- `AVAILABLE`: the evidence was computed and passed its applicability checks.
- `CAUTION`: the evidence was computed but has a documented limitation.
- `HOLD`: a failed prerequisite prevents the intended interpretation.
- `NOT_ASSESSABLE`: the data or model do not support the calculation.

### Conclusion stability

- `STABLE`: the conclusion persists across all required sensitivity checks.
- `CONDITIONALLY_STABLE`: it persists within an explicitly bounded subset of
  supported conditions.
- `SENSITIVE`: at least one required, credible specification changes it.
- `NOT_ASSESSED`: the required sensitivity analysis was not run.

These states are not pass/fail labels for a person, rater, institution, or
assessment. A conclusion may be computationally available and still be
sensitive.

### Workflow, learning, and claim readiness

- `workflow_complete`: the required interface actions for a bounded route were
  performed, including a route that ends with a documented failed analysis.
- `learning_complete`: the learner supplied the defined evidence of
  understanding for a guide node; merely opening a card is insufficient where
  a formative check is specified.
- `claim_ready`: every required scientific prerequisite for a particular
  ClaimID is satisfied for the current StudyContextID, AnalysisID, and
  InterpretationID.

These flags are independent. A learner may complete a node while the analysis
remains `HOLD`; a failed or unassessable run may be archived as a reproducible
no-conclusion record; and neither guide completion nor an acknowledgement may
make a substantive claim ready.

## Implementation record

### Help and audience-language shadow foundation — 2026-07-24

Delivered locally, without changing the visible Streamlit route:

- Added a pure three-layer terminology registry for task-facing, method, and
  technical vocabulary. Twenty initial review-tracked concepts retain exact
  reproducibility labels without exposing them through the Standard
  projection.
- Added immutable Help topic, section, link, target, and verified-context
  contracts plus a pure navigation reducer. Opening, changing, invalidating,
  localizing, and returning from Help carries presentation identity only and
  cannot mutate an analysis, evidence conclusion, or learning state. Route
  construction is reducer-only, including protection against forged
  `dataclasses.replace()` transitions; closing Help clears route history.
- Registered 22 task-centered Help topics, a reserved localized safe fallback,
  42 source/return targets, and 74 Help links, including exact mappings for 18
  existing popover keys. No invented topic alias is retained. All terminology
  and problem references resolve, every return focus is owned by its target,
  and the graph has no orphan nodes.
- Separated the fallback from Help-home entries and placed all 42 intended
  surface, focus, and panel identities in an independent provisional catalog.
  Consumer typos can no longer create their own valid targets; the remaining
  activation slices must bind and test each catalog entry against the real
  Streamlit surface.
- Added 11 privacy-safe user-problem records for input, mapping, estimation,
  resource, and residual-structure failures. A rendered notice can contain
  only a registered problem code, a bounded occurrence phase, and a freshly
  generated random support reference; exception text, paths, uploaded values,
  and analysis identifiers are not retained. Generic application failures use
  their own troubleshooting topic rather than an estimation-failure route.
- Added matched English and Japanese seed copy for 100 terminology keys, 34
  problem keys, 437 Help-topic keys, and 24 Help-navigation keys, including
  explicit boundaries for fit, residual structure, differential interaction,
  privacy, estimator choice, and reporting claims. Every reporting example is
  structurally paired with a warning to verify its conditions first.
- Added contract, graph, reducer, locale, exposure, privacy, and standalone-core
  tests. The seed Help copy contains no active external-runtime handoff and no
  default implementation vocabulary.
- Marked terminology and Help content `in_review`; SLA-method and UX review
  remain explicit gates before this copy replaces the visible application.
  The current ClaimBoundary/reference IDs and surface-target catalog are
  provisional contract identifiers, not completed bibliographic records or
  evidence that a live UI target exists. Those records and the remaining
  adapter bindings are separate activation gates across R2 and R3.

This was deliberately delivered first as a shadow slice. The bounded R2
activation below now consumes it for the global Help entry and first-failure
paths. The first R3 vertical slice now also activates the fit-scatter link;
legacy Help sections and every other chart popover remain compatibility
surfaces until their exact registered adapters pass the later gates.

### Persistent Help and safe-problem activation — 2026-07-24

Delivered locally:

- Added a persistent bilingual Help launcher to the application shell and a
  full main-content Help surface available before data are present. It uses
  stable topic/link IDs, renders only the Standard audience layer by default,
  checks lifecycle applicability, and excludes the system fallback from the
  22-topic picker.
- Kept the source data and settings controls alive in the sidebar while Help
  owns the main content. Paste drafts, delimiters, uploads, selected source,
  and keyed analysis settings therefore remain in Streamlit's normal widget
  lifecycle; the Help route never consumes a run or refit trigger.
- Added an in-Help topic selector, persistent English/Japanese switching, an
  explicit return control, localized return status, and a safe generic route
  for unknown, removed, or inapplicable destinations. Locale changes retain
  route identity rather than translating IDs.
- Activated privacy-safe problem notices on paste/upload, wide-to-long,
  simulation-threshold, and covariate-preview parsing, plus primary estimation
  failures. The visible notice contains registered bilingual copy, a random
  support reference, and exact registered Help actions; raw exception text,
  paths, input values, and identifiers remain absent from the public UI even
  if a legacy technical-error environment variable is present.
- Replaced visible “technical details” failure gates on these paths with
  expanded, task-facing next-action panels connected to Help.
- Added AppTests for no-data open → locale change → return, source-state
  preservation, zero estimator/refit calls, unknown-link recovery, denied and
  allowed lifecycle states, pending-trigger preservation, and sentinel
  exception redaction through both direct and simulation-threshold
  problem-to-Help paths. The zero-call guard instruments both the cached core
  namespace and direct module-level estimation/refit functions.

This is a bounded activation, not a claim that all 42 provisional targets are
live. The global start target and first problem routes are connected; legacy
popover, guided, result, report, and download adapters still require exact
surface/focus verification. Topic copy, ClaimBoundary IDs, and reference IDs
remain `in_review` until SLA-method and UX review are recorded.

### First contextual Help return slice — 2026-07-24

Delivered locally as a deliberately narrow R3 slice:

- Activated only `fit_scatter` among the 18 registered core popover aliases.
  Its button resolves the exact registered
  `link.popover.fit_scatter` -> `help.results.fit#fit-scatter` edge and the
  exact `target.results.figure.fit_scatter` /
  `focus.results.figure.fit_scatter` return pair. The other 17 aliases and the
  two legacy posterior popovers have no detailed-Help action in this slice.
- Required a current fitted context built from the retained result's
  AnalysisID and its stored sample/real provenance. A legacy output with no
  stored provenance is not inferred from the current sidebar selection, and
  missing, stale, mismatched, or malformed context fails closed.
- Added separate allowlists for presentation projection and target-owned
  scientific view state. Return restores the Essential or Full fit-details
  panel plus only `fit_df_method_method` and `fit_df_method_cap` through the
  return intent. Native MML prior-SD sensitivity draft controls and misfit
  ranking display controls use separately validated, AnalysisID-scoped durable
  session state, so Streamlit widget cleanup cannot reset them while Help owns
  the main surface or carry them into a different fitted analysis. Analysis
  settings, data, evidence, and guide state remain untouched.
- Removed the legacy imported-FACETS fit-table comparator from the default Fit
  Details journey. Its compatibility function remains isolated and is not an
  active Help dependency or standalone-core workflow; no external comparison
  upload needs to be retained by this return adapter.
- Bound the return intent to the exact link, target, focus, projection,
  AnalysisID, provenance, schema version, and validated source-view values.
  Missing or altered intent opens the reserved localized fallback instead of
  guessing a destination or reporting a successful return.
- Rechecks the current fitted phase, AnalysisID, and sample/real provenance
  before restoring any panel or source value. Browsing a static Help topic can
  still return through an intact origin intent, while an A-to-B analysis or
  provenance change fails before any A-owned value is written into B.
- Suppressed analysis triggers for the one return rerun without consuming
  pending user-initiated run/quick-start flags. The real-app AppTest verifies
  zero estimator/refit calls, the same AnalysisID, bilingual route retention,
  recovery of a `both` / `12.5` fit-display convention after Streamlit removes
  the hidden widgets, retention of native MML-sensitivity and misfit-ranking
  drafts, and consumption of the intent at the real figure marker.
- Made a new Help journey an explicit cleanup boundary for an abandoned
  `awaiting_render` intent, its one-shot trigger suppression, and an old return
  notice. A target that was not rendered therefore cannot poison a later
  global Help journey.
- Anchored the real localized Fit scatter `H3` with the registry-owned marker
  ID and, only after a valid intent reaches that rendered heading, emits one
  static focus/scroll request. The progressive-enhancement script contains no
  AnalysisID, provenance, locale copy, or input data; validates the exact
  connected `H3`; waits two animation frames; avoids stealing focus after a
  user move; and guards duplicate execution. The visible localized return
  confirmation is a non-live caption, so an empty ARIA marker and duplicate
  success-alert announcement are not introduced.
- Encoded adapter maturity as `APPTEST_ONLY` or `BROWSER_ACCEPTED`. Exactly one
  AppTest-only pilot is allowed; every other active adapter requires a stable
  retained browser-evidence reference. `fit_scatter` remains the sole pilot,
  has no browser evidence ID, and is not browser accepted. The executable
  matrix and promotion procedure are in
  [`docs/help_browser_acceptance.md`](docs/help_browser_acceptance.md).

This is not R3 completion. AppTest verifies the native heading anchor, exact
static controller body, JavaScript opt-in flag, one-shot emission, localized
caption, and fail-closed absence for invalid or stale intents; it does not run
the browser script. The exact `document.activeElement`, viewport placement,
visible focus outline, subsequent Tab order, CSP behavior, zoom, reduced-motion
behavior, and single screen-reader announcement still require the browser
acceptance pass. A remaining popover can become the next sole pilot only after
the current pilot passes those gates and its reviewed evidence reference is
recorded; defining the checks is not evidence that they passed.

### Evidence-contract foundation — 2026-07-24

Delivered:

- Added `mfrm_app/evidence.py` with versioned contracts for deterministic
  analysis identity and lineage, evidence records, sensitivity comparisons,
  sensitivity decisions, workflow dispositions, and stable reason codes.
- Fixed the distinction between computation state and conclusion stability in
  executable validation rules rather than relying on display strings.
- Added order-independent sensitivity synthesis: a conclusion is `STABLE` only
  when every required credible variant succeeds, is comparable, and leaves the
  conclusion unchanged.
- Bound each prespecified sensitivity label to an embedded, self-validating
  baseline/variant AnalysisIdentity ledger, evaluate a closed versioned rule
  against archived metrics with fail-closed three-valued logic, and
  content-fingerprint the complete aggregate outcome.
- Added strict bundle restoration: sensitivity and workflow decisions are
  accepted from JSON only when their linked plan, evidence, and source records
  reproduce the saved decision exactly; unknown fields, type-coerced arrays,
  registry-version mismatches, and unresolved EvidenceIDs are rejected.
- Added a static core-boundary test that prevents `mfrm_app/` from importing
  Streamlit, invoking subprocess/dynamic external engines, or embedding new
  TAM/sirt/mirt/mfrmr/CmdStan handoffs.
- Removed the legacy parity-fixture exporter from the CLI, `make verify`, and
  GitHub CI, and removed cross-package, Stan, and Posterior compatibility tests
  from the standalone `--self-test` registry.
- Updated contribution and release guidance so new statistical computation is
  modular, Python-native, reason-coded, and validated with native fixtures.
- Added `mfrm_app/**` and this roadmap to CI path filters so changes to the
  modular core cannot bypass the Python test workflow.
- Removed the repository-external `Simulation/engines/mfrm_python.py` loader,
  its import-time execution path, and the unreachable shared-engine dispatch.
  An entrypoint AST guard now prevents that runtime loader from returning.
- Added `mfrm_app/readiness.py` and routed final-report readiness through it.
  The original five display columns remain compatible, while every row now
  carries a deterministic AnalysisID/EvidenceID, explicit computation and
  stability states, stable ReasonCode, scope, source-artifact links, and a
  validated full EvidenceRecord export.
- Separated the different meanings previously hidden behind final-readiness
  `Review`: non-convergence is `HOLD`, an executed diagnostic requiring review
  is `CAUTION`, and skipped PCA, bias, or strict-marginal evidence is
  `NOT_ASSESSABLE`.
- Added `mfrm_app/mml_prior_sensitivity.py` and connected fixed-prior-SD MML
  refits to exact baseline/variant identities, a prespecified executable plan,
  sensitivity records, and a fail-closed SensitivityDecision. The executed
  plan and JSON contract are now included in the Python-native sensitivity
  archive, and final readiness consumes the same decision.
- Closed prior-SD edge cases: free-population-SD fits and baseline-only grids
  are rejected as inapplicable/unassessable, `max_values` fixes the plan before
  refitting, and run failure, non-convergence, missing rule metrics, and
  non-comparable measures cannot produce a false stable result.
- Hardened the sensitivity comparison against false stability: outcome flags
  require exact booleans; measure and population-term keys must be complete,
  finite, and duplicate-free; rank evidence must be evaluable for every
  eligible facet; and each returned refit must attest the planned fixed-SD,
  engine, quadrature, tolerance, iteration, no-plausible-value controls, and
  unchanged population/anchor/constraint semantics.
- Revalidate saved MML bundles by reconstructing their identities, evidence,
  plan, records, and synthesized decision before readiness or export consumes
  them; a mutated same-AnalysisID payload is ignored rather than displayed.
- Unified the MML sensitivity baseline with the final-readiness AnalysisID and
  reject stored sensitivity evidence whose baseline identity differs from the
  current fitted result.
- Added reconstructable AnalysisIdentity, EvidenceRecord, sensitivity-contract,
  settings, and result-table sidecars to the standard Python-native ZIPs. Their
  source-artifact links are tested against the actual archive contents rather
  than surviving only as in-memory DataFrame attributes.
- Reduced `--release-check` to the standalone Python contract: its JSON has
  only release metadata, native readiness, and native limitation rows, and it
  no longer generates or reports external inventories, Stan handoffs, or
  `mfrmr` migration coverage.
- Removed external simulation inventories/templates, R/Julia reproducibility
  bundles, TAM/sirt/mirt handoffs, Stan/Posterior packages, and `mfrmr`
  migration tables from the default Downloads and deterministic demo archives.
  Python engine/local-batch/self-contained-JMLE assets and the Evidence
  sidecars remain first-class native outputs.
- Removed Posterior Viewer and advanced Stan generation from the sidebar,
  Report, Help, and tutorial routes. Public Yardstick and rating-scale recode
  downloads now enter Python-only generator branches directly instead of
  constructing legacy R/Julia assets and filtering them afterward.
- Split dormant compatibility tests behind the `legacy_compat` marker. The
  standard `make apptest`, `make verify`, and GitHub workflow now test only the
  standalone Python release surface; frozen compatibility fixtures are not
  release evidence.
- Added demo-specific handoff filenames so the generated checklist and
  manuscript handoff point to the files that the deterministic demo actually
  writes.
- Applied the demo artifact profile to every generated table, Markdown/HTML
  report, visual evidence map, and binder README, and added the previously
  referenced `export_privacy_manifest.csv` as an actual synthetic-demo
  artifact. Recursive archive verification now finds no stale app-download
  names in the demo packages.

Remaining before the product-boundary completion gate:

- Record and execute the explicit retain/deprecate/remove decision for dormant
  generator definitions, unreachable locale copy, and frozen compatibility
  fixtures. Public Help, Posterior, Report, Downloads, demo, CLI, self-test, and
  CI routes are already isolated; none of the dormant surfaces may re-enter
  them.

## Delivery sequence

The order below is dependency-driven. Later outcomes should not be built on
unversioned tables or UI-only logic.

### Product boundary and modular contracts

Objective: make the standalone-Python boundary enforceable and give new work a
stable home outside the large Streamlit entrypoint.

Deliverables:

- Record this roadmap as the active product plan and retain completed release
  plans as historical documents.
- Define versioned Python contracts for analysis identity, evidence records,
  sensitivity results, decision states, and skip/hold reasons.
- Move computational network and decision helpers into focused `mfrm_app/`
  modules while retaining compatibility wrappers where tests or public imports
  require them.
- Separate UI rendering, computation, and export-frame construction.
- Update contributor guidance so new statistical work is modular by default.
- Add a boundary test that rejects runtime imports, subprocess calls, and UI
  actions that invoke external estimation engines in the core workflow.
- Classify existing legacy external-comparison surfaces as compatibility-only
  so they cannot be selected as prerequisites by new features, appear in the
  default workflow, or gate a release.
- Produce a separate compatibility decision for legacy code generators,
  posterior-result viewers, cross-engine fixtures, and handoff downloads. Do
  not remove them incidentally during statistical refactoring.

Completion gate:

- Existing public fixtures and result columns remain compatible or have an
  explicit migration note.
- All tests pass after each extraction, with no statistical change hidden in a
  modularization commit.
- Every unavailable computation returns a stable reason code rather than an
  empty result or silent exception.
- Core estimation, diagnostics, sensitivity, and reporting execute without an
  external engine installed.

### Resilient rating design

Objective: progress from a static connectivity summary to an analysis of how a
design fails and how it can be strengthened.

Deliverables:

- Maintain distinct graph views for the global multipartite design, direct
  rater-person links, rater-task/criterion coverage, and projected rater
  overlap. Do not compare density across unlike graph definitions.
- Add weighted edge strength, minimum shared cases, bridges, articulation
  points, component structure, workload imbalance, and coverage summaries.
- Add deterministic perturbation audits for:
  - leave-one-rater-out;
  - leave-one-task/item/criterion-level-out where applicable;
  - removal of singleton or weakest links;
  - prespecified random rating loss; and
  - concentrated missingness in one facet region.
- Report the proportion and identity of perturbations that disconnect the
  design or materially reduce overlap.
- Generate constrained repair candidates, such as additional shared ratings,
  while displaying burden and the exact diagnostic each candidate repairs.
- Keep connectivity and precision distinct: a repaired connection is not
  automatically an adequately precise design.

Completion gate:

- Complete, chain, star, bridge-dependent, disconnected, and imbalanced
  synthetic fixtures return their known graph properties.
- Leave-one-level-out results are deterministic and use the same validated
  observation set as the fitted analysis.
- Every proposed repair reconnects the targeted synthetic design when a repair
  is feasible, and no repair is labelled as sufficient evidence of precision.
- Exported graph definitions, thresholds, removed observations, and
  perturbation seeds reproduce the displayed result.

### Coherent diagnostic evidence

Objective: combine diagnostics without treating them as a vote or collapsing
them into an unsupported universal score.

Deliverables:

- Normalize fit, covariance, residual PCA, PCA stability, DIMTEST,
  local-dependence, category, bias/interaction, reliability, design-network,
  severity-network, and halo-network outputs into evidence records.
- Record for each item of evidence:
  - analysis identity and scope;
  - data and configuration fingerprints;
  - applicability prerequisites;
  - computation state and reason code;
  - observed statistic and uncertainty when available;
  - interpretation boundary;
  - recommended next inspection; and
  - source table or figure.
- Build an evidence-coherence view that highlights agreement, tension, and
  unresolved prerequisites. It must not use majority voting.
- Encode explicit interpretation rules for common tensions, including:
  - a residual PCA signal with unstable bootstrap loadings;
  - a DIMTEST result that conflicts with residual PCA;
  - acceptable global fit with concentrated local dependence;
  - high reliability with fragile network connectivity; and
  - a rater severity signal that disappears after residual adjustment.
- Use the same evidence records in guided UI, research UI, readiness checks,
  narrative drafts, and exports.

Completion gate:

- Contradictory synthetic evidence produces a documented tension state, not a
  false all-clear or automatic dimensionality claim.
- No final-readiness statement cites an unavailable, stale, or differently
  fingerprinted diagnostic.
- Every user-facing recommendation links to its evidence and caveat.
- Locale and export tests confirm that status meanings are consistent in
  Japanese and English.

### Conclusion robustness

Objective: evaluate the stability of the decisions users actually make, not
only the stability of individual coefficients.

Deliverables:

- Define a versioned sensitivity specification before running a grid.
- Support applicable checks across:
  - JMLE/MML and supported MML engine choices;
  - fixed population prior SD settings;
  - RSM/PCM/GPCM only when the data and comparison question permit it;
  - category handling and declared score support;
  - influential facet levels;
  - reasonable missingness perturbations; and
  - documented fit-statistic conventions.
- Track conclusion-level outcomes such as:
  - facet severity direction and rank bands;
  - misfit and bias-review flags;
  - reliability interpretation;
  - category-functioning warnings;
  - dimensionality-screen wording; and
  - final claim availability.
- Distinguish a numerical change from a substantive conclusion change.
- Preserve failed specifications in the denominator and explain their failure;
  do not summarize only successful refits.
- Produce a compact stability map plus a full specification ledger.

Completion gate:

- Fixed fixtures demonstrate each stability state and prevent stale-base-run
  reuse.
- The specification ledger reproduces every displayed comparison from frozen
  inputs and seeds.
- A failed or unsupported refit cannot strengthen the final claim.
- Rank, flag, and narrative changes use documented thresholds rather than
  post-hoc visual judgment.

### Calibrated simulation and stress testing

Objective: show under which data-generating and design conditions the app's
diagnostics and decisions behave as intended.

Deliverables:

- Complete the modular condition/scenario system for immutable, serializable,
  fingerprinted simulation specifications.
- Expand native DGP coverage for supported RSM, PCM, and bounded GPCM use cases
  without implying support for an unfitted model family.
- Add prespecified stressors for:
  - disconnected and bridge-dependent assignments;
  - unequal rater workload and restricted coverage;
  - severity spread and extreme raters;
  - halo and supported local interactions;
  - rater drift across occasions;
  - missing completely at random and concentrated design loss;
  - sparse and unused score categories; and
  - non-normal person distributions where the estimator permits evaluation.
- Evaluate bias, RMSE, interval coverage, rank recovery, flag rates,
  disconnection rates, conclusion reversals, numerical failure, and Monte
  Carlo uncertainty.
- Maintain an ADEMP-style scenario ledger that separates development checks,
  calibration runs, and frozen validation simulation runs.
- Add checkpoint/resume only as a Python-native computation feature, not as an
  external-engine exchange contract.

Completion gate:

- Repeating a scenario with the same specification reproduces task identities,
  seeds, truth, and summaries.
- Scenario additions do not change the random streams of existing frozen
  scenarios.
- Rates retain all requested replications in their denominator and report
  numerical failures separately.
- Proportions include Monte Carlo uncertainty; small development runs are not
  presented as performance evidence.
- Each public diagnostic claim is backed by a named set of frozen simulation
  conditions or is explicitly labelled uncalibrated.

### Guided decisions and research transparency

Objective: serve first-time and expert users without maintaining two different
statistical implementations.

This work is an information-architecture consolidation, not another tutorial
surface layered onto the existing onboarding, Tutorial, Guided, and Help
systems. The application already contains substantial explanatory material.
The long-term problem is to make one route and one evidence state authoritative
while rendering different levels of detail for different users.

#### P0 research, interpretation, and learning contracts

The interface will not be redesigned until the following contracts are frozen
in testable, widget-independent form. The detailed implementation contract is
maintained in
[`docs/study_context_guidance_contract.md`](docs/study_context_guidance_contract.md).

| Contract | Minimum responsibility | Boundary |
|---|---|---|
| `StudyContext` | Purpose, stakes, target population/domain, construct and subconstructs, performance mode, intended score interpretation and use, intended decision, and unsupported inferences | Does not assert validity, fairness, or fitness for use |
| `RatingDesignRecord` | Object of measurement, nuisance facets, crossed/nested/anchored/fixed/random roles, planned and observed assignment, common performances, occasions/order, double scoring, adjudication, rubric version, and missingness plan | Graph connectedness alone does not establish comparability or generalizability |
| `ScoringProcessRecord` | Rater selection, training, qualification, monitoring, rescoring, adjudication, and benchmark process | Statistical flags do not authorize punitive personnel action |
| `AnalysisSpec` | Canonical data mappings, model, estimator, identification, diagnostic, sensitivity, and seed settings independent of Streamlit widget visibility | View, guide, language, and theme are excluded |
| `GuidanceState` | Stable route/node IDs, progress, sample context, bindings, invalidation, and lifecycle events | Cannot mutate computation or conclusion state |
| `HelpTopic` / `HelpLink` | Stable concepts, audience vocabulary, applicability, claim boundaries, references, exact contextual targets, and return/focus behavior | Help navigation cannot mutate science/learning or present stale analysis context as current |
| `EvidenceIssue` | One prioritized question with all relevant EvidenceIDs, tensions, applicability, boundaries, and next actions | Never hides other `HOLD`, `CAUTION`, or required-unassessed evidence |
| Interpretation and decision records | Stable claim, limitation, action, authorship, and review-state identities linked to study and analysis context | An analysis record may document that no interpretable conclusion is available |

Planned absence and unexpected missingness are represented separately. At a
minimum the application distinguishes planned incomplete assignment,
test-taker nonresponse, not-reached responses, technical failure, rater skip
or unratable response, criterion-level absence, administrative exclusion, and
adjudication-only ratings. An unexplained `NA` is not silently treated as a
planned absent design cell.

For analytic rubrics, the study context must state whether criteria are
intended as one construct, distinct subconstructs, or a reported composite.
Treating criterion as a facet does not by itself justify an overall person
measure. If the intended claim requires a multidimensional model outside the
native scope, the claim is `HOLD` or `NOT_ASSESSABLE`, not rescued by a clean
residual screen.

All currently visible explanatory text is part of the evidence surface. Help,
tutorial, chart guides, popovers, warnings, captions, and exported prose are
therefore audited in G0 against
[`docs/interpretation_copy_audit.md`](docs/interpretation_copy_audit.md).
Later Evidence-linked copy cannot coexist with stronger legacy statements
elsewhere in the application.

The file-level remediation sequence, registry design, concrete bilingual copy,
Help lifecycle behavior, migration slices, and test changes are specified in
[`docs/help_and_terminology_remediation_plan.md`](docs/help_and_terminology_remediation_plan.md).

#### Non-negotiable UX principles

- The guide is optional, skippable from every step, restartable, and resumable
  within the session.
- Skipping the guide never deletes uploaded data, fitted results, the current
  AnalysisID, EvidenceRecords, or a draft decision note.
- Acknowledging a guide card never changes `HOLD`, `NOT_ASSESSABLE`,
  `CAUTION`, or conclusion-stability states. Only new data, settings,
  computation, or evidence may change an analysis state.
- Guided and research surfaces use the same estimator results,
  AnalysisIdentity, EvidenceRecords, sensitivity records, and reason codes.
- Tutorial progress is not analysis progress, and guide completion is not a
  certificate of MFRM competence or evidence that a real-data analysis is
  valid.
- The guide teaches the standalone Python workflow. TAM, R, Julia, external
  posterior ingestion, and cross-engine ZIP compatibility do not appear in the
  guide or become prerequisites for it.
- The sample route is allowed to teach operation and interpretation, but a
  successful sample run must explicitly state that it does not validate the
  user's data or design.
- The guide remains a thin projection over the production computation path; it
  does not maintain canned results or a separate tutorial estimator.

#### Entry routes

The landing surface offers three clear choices without asking users to label
themselves by an assumed ability level.

| Entry route | Behaviour | Guidance level |
|---|---|---|
| Learn with a sample | Continue from the landing surface into the deterministic guide | Full first-run guide |
| Start with my data | Open the real-data workflow | Contextual readiness and evidence guidance |
| Continue without the guide | Open the normal Standard view | No tutorial sequence |

The landing surface itself is the `welcome` node. It records the selected
route and moves forward in one action; users are not asked to make the same
sample-versus-data choice on a second welcome screen. A future pre-collection
route may begin from `Plan`, but it must reuse the same StudyContext and
RatingDesignRecord rather than become a fourth onboarding system.

The fixed tutorial uses a built-in sample because an arbitrary uploaded data
set may be disconnected, sparse, malformed, or non-convergent. Real data
therefore receives evidence-dependent coaching rather than being forced
through a linear tutorial. After the basic sample, an optional single-warning
practice scenario may teach one `CAUTION` or `NOT_ASSESSABLE` branch at a time.
It must not combine many artificial problems into one confusing exercise.

#### Five-step guided first-run route

Stable internal node IDs, rather than translated labels, drive routing.

| Node ID | User-facing purpose | Learning/workflow completion | Scientific effect | Primary output |
|---|---|---|---|---|
| `welcome` | Choose the sample or leave the guide | A landing route is selected | None | Bound route |
| `data_check` | Confirm roles, mappings, design, and readiness | The required concepts and current preflight outcome were reviewed | `HOLD` may remain and directs repair | Data fingerprint, design record, and readiness state |
| `estimate` | Run the bound specification and inspect its numerical outcome | The current run outcome, including failure, was reviewed | Only a successful current fit creates usable fitted evidence | Analysis identity or a versioned failure record |
| `evidence_review` | Examine the priority evidence issue and its claim boundary | The learner distinguishes one supported from one unsupported statement in the sample exercise | Required evidence may remain `HOLD`, `CAUTION`, or `NOT_ASSESSABLE` | Evidence-linked boundary and next action |
| `archive` | Save a bounded interpretation/decision draft | A supported claim or explicit no-conclusion statement, a limitation, and a next action are recorded | Archiving does not certify or improve claim readiness | Human-readable Analysis Brief plus JSON sidecar |

The visual sequence is:

```text
Landing/welcome -> Data check -> Estimate -> Read evidence -> Save record
                                      |
                                      +-- HOLD: repair, or archive why no
                                          interpretable conclusion is available
```

The guide is a five-node learning projection over the full research workflow,
not a replacement vocabulary. The mapping is shown throughout the guide so
the final transition does not create a second onboarding cliff.

| Guide node | Long-term research stages |
|---|---|
| `welcome` | Plan |
| `data_check` | Plan -> Check |
| `estimate` | Estimate |
| `evidence_review` | Diagnose -> Stress -> interpretation substep |
| `archive` | Decide -> Archive |

The canonical seven-stage workflow remains stable. `Interpret` is an explicit
substep between `Stress` and the use decision inside the `Decide` stage; its
InterpretationRecord and the subsequent DecisionRecord remain distinct
identities even when the navigation groups them on one screen.

Each guide screen has at most one primary card, three headline metrics, and one
state-changing primary call to action. `Exit guide` is persistent escape
chrome, not the screen's secondary action budget. Non-mutating Help and detail
disclosures may remain available without competing visually with the primary
action. Long tables, model theory, all-panel tours, and publication options
remain in contextual detail or Help.

#### Skip, exit, resume, and restart contract

The guide state has four lifecycle values: `NOT_STARTED`, `ACTIVE`, `SKIPPED`,
and `COMPLETED`. It also records the active node, last valid node, completed
nodes, entry route, bound data fingerprint, bound AnalysisID, reviewed
EvidenceID, and the node at which the guide was skipped.

- Before starting, the low-emphasis action is `Continue without the guide`.
- During a guide, the action is `Exit guide`, because users are skipping an
  explanatory sequence rather than statistical prerequisites.
- Exiting before estimation opens the current data/check surface.
- Exiting after estimation opens Current focus without refitting.
- Exiting while reading evidence preserves the same EvidenceRecord and detail
  target in the normal workspace.
- Resuming uses the first valid incomplete node, not blindly the previously
  displayed node.
- Restarting resets guide progress only; it does not clear input data, results,
  AnalysisIdentity, EvidenceRecords, or decision drafts.
- A skipped guide does not repeatedly reopen or display a blocking prompt in
  the same session. A quiet resume action remains in Help or the sidebar.
- Starting or resuming the sample guide while a real-data analysis exists uses
  an isolated tutorial context. It snapshots and restores the real-data
  `AnalysisSpec`, fitted identity, evidence, and draft notes; the sample may
  never overwrite them or silently become the current research analysis.
- Guide state is session-only in the initial implementation. Persistent browser
  preference is deferred until observed user friction justifies the privacy,
  accessibility, and maintenance cost.

No confirmation dialog is required for exiting because the operation is
non-destructive. The UI states that data and current results are retained.

#### Standard and Detailed views

Completing or skipping the tutorial leads to the Standard view, not directly
to every expert panel. This prevents an onboarding cliff in which a simple
guide ends at the existing high-density interface. The user-facing names avoid
using `Guided` for both the guide and a workspace, and avoid implying that only
a `Research` workspace is suitable for research.

| Stable view ID | User-facing projection | Purpose |
|---|---|---|
| `standard` | Standard view | One current evidence issue, its boundary, and next action |
| `detailed` | Detailed view | Resolved settings, full diagnostics, sensitivity ledgers, simulations, and exports |

The view selector is a presentation projection only. It never applies
estimation defaults, writes hidden widget values, changes a draft or fitted
`AnalysisSpec`, changes AnalysisID, invalidates a fit, triggers estimation, or
starts a hidden diagnostic. `AnalysisSpec` is canonical state outside
Streamlit widgets so temporarily undisplayed advanced controls retain their
values. The interface distinguishes fitted settings from an edited draft and
offers an explicit return to the fitted specification. Analysis depth remains
an Estimate setting because it changes which diagnostics and artifacts are
computed.

The canonical seven-stage navigation remains:

```text
Plan -> Check -> Estimate -> Diagnose -> Stress -> Decide -> Archive
```

Existing surfaces migrate without duplicating computation:

| Existing surface | Long-term home |
|---|---|
| Start | Plan / Check |
| First Read | Diagnose / Decide |
| Results | Estimate / Diagnose |
| Diagnostics | Diagnose / Stress |
| Figures | The relevant result or diagnostic panel |
| Report & Export | Decide / Archive |
| Learn | Persistent Help, not an analysis stage |

The Standard view starts with one Current focus card backed by an
`EvidenceIssue`, not a single convenient result. It prioritizes stale results,
non-convergence, disconnected designs, required-but-unassessable evidence,
material diagnostic cautions, evidence tensions, unassessed sensitivity,
reporting limits, and then the standard next result. It also shows the total
number of `HOLD`, `CAUTION`, and required-unassessed issues and why the current
issue was selected. Additional items remain available without competing with
the highest-priority action.

#### Information density and terminology

During the sample guide, the tutorial binds a documented sample
`AnalysisSpec`; it does not simulate a choice among display, control, and
computation presets. Standard-view sidebar content is limited to data, column
mapping, a resolved estimation summary, and the run action. Optimizer,
tolerance, anchors, visual preferences, custom diagnostic options, simulation,
and publication controls stay collapsed or appear in Detailed view. Hiding a
control never resets its canonical value.

Guide labels use plain language followed by the canonical technical term
where useful, for example:

- Model discrepancy (`Fit`);
- Residual-structure check (`Residual PCA`);
- Element differentiation (`Reliability / separation`); and
- Differential interaction screen (`facet interaction`; legacy `bias` labels
  are retained only where required for a clearly marked method reference).

Rater flags are framed as review prompts for the design, observation count,
rating examples, rubric, and scoring process. They do not label a rater as
accurate, fair, competent, or unsuitable, and they do not automatically
recommend exclusion, retraining, or personnel action.

The default run action uses standalone product language such as `Run MFRM
estimation`, not a label that implies runtime dependence on FACETS or another
external product.

#### Audience language and Help continuity

Implementation vocabulary is not automatically user vocabulary. Internal
enums, record types, IDs, cache/schema terms, and support codes remain in the
core contract and machine-readable sidecar; they are not the default wording
of the guide, Standard view, or Help landing page.

| Presentation layer | Intended language | Examples |
|---|---|---|
| Guide and Standard view | Plain task, status, boundary, and action language | “Study purpose,” “analysis settings,” “why this result is limited,” “stop and repair,” “analysis reference” |
| Detailed view and method Help | Canonical MFRM/SLA terminology with a short definition and applicability boundary | Residual PCA, Infit/Outfit, separation, differential interaction, JMLE/MML |
| Technical details, JSON, and privacy-safe support | Exact implementation names and stable identifiers | `StudyContextID`, `AnalysisID`, `EvidenceID`, `ReasonCode`, schema version |

For example, the default interface says “current evidence question,” not
`EvidenceIssue`; “interpretation draft,” not `InterpretationRecord`; and “why
this is limited,” not `ReasonCode`. A stable analysis reference may be shown
and copied in Detailed/technical details, but the user is not instructed to
reason in terms of internal class or enum names. Status codes such as `HOLD`
are paired with a localized action label such as “Stop and repair”; the code is
secondary technical detail.

Help is one connected content system, not a second tutorial or a collection of
unrelated tables. A versioned `HelpTopic` registry owns stable topic and
concept IDs, locale keys, audience layer, prerequisites, applicability,
ClaimBoundaryIDs, references, related workflow targets, and review metadata.
Every guide node, EvidenceIssue, chart popover, warning, and advanced control
links through a `HelpLink` to that registry.

Persistent static Help is reachable before data entry, during mapping and
preflight, after a failed first run, in every `HOLD`/`NOT_ASSESSABLE` state, and
after a successful fit. It never depends on a fitted result merely to render.
Only the optional “For this analysis” block requires a verified current
analysis binding.

The Help path has three depths:

```text
Inline boundary -> focused Help topic -> method/reference detail
       ^                    |
       +------ return ------+
```

- Inline Help answers what is shown, how to read it, and the most important
  limit without leaving the current task.
- `Open detailed help` is an actionable control, not a caption telling the
  user to find another page manually. It opens the exact topic/section,
  preserves the originating stage/card and current analysis context, announces
  the new heading after rerun, and provides a return action.
- Opening Help, following a HelpLink, searching, or returning changes no
  AnalysisID, setting, evidence, claim readiness, or learning-completion state.
- A missing/unknown HelpTopic fails a contract test and renders a safe generic
  fallback with a support reference; it never silently removes the Help
  control.
- Topic/target routing uses stable IDs only, never translated labels, screen
  names, or Japanese/English substring matching. The glossary has one
  canonical source for inline definitions, Help search, and exports.
- A claim-critical limitation remains inline. Help may explain it but cannot
  be the only place where it appears.
- Standard and Detailed views use the same topic and ClaimBoundaryID. Detailed
  view may reveal method/reference depth, but it does not receive a different
  substantive conclusion.

The Help home supports search and browsing by users' questions rather than
internal object names: research question; design and missingness; model and
estimation; diagnostics; interpretation and fairness; reporting and review;
reproducibility, privacy, and ethics; troubleshooting; and glossary and
references. Each focused topic states the question, prerequisites, what is
computed, what it can and cannot show, the next action, safe/unsafe reporting
language, references, related screens, and return target. Dense developer or
reviewer inventories and machine filenames move to optional technical details
or exported support/reviewer artifacts, not the guided Help path.

#### Evidence-linked teaching and claim boundaries

The evidence-review card is rendered directly from an EvidenceRecord. It shows
the computation state, what the evidence supports, what it does not establish,
the recommended action, the next inspection, and the source artifact. No
parallel tutorial copy is allowed to silently contradict the evidence ledger.

All guide-facing statements are bounded by the current data, model, and
computed scope. In particular:

- convergence means that the numerical stopping rule was met; it does not
  establish fit, validity, fairness, or conclusion stability;
- fit statistics are screening evidence, not universal pass/fail rules or
  labels for a person or rater;
- residual PCA can identify residual-structure signals but cannot prove
  unidimensionality when no large component is detected;
- orderly categories do not prove scale validity, and a disordered threshold
  does not mechanically require category collapse;
- rater reliability can indicate distinguishable severity differences; it is
  not evidence that raters agree or are interchangeable;
- dropping missing rows does not establish an ignorable missingness mechanism;
  and
- no flagged interaction within the computed pairs and detection conditions is
  not proof that bias or unfairness is absent.

Normal-state language therefore uses formulations such as `no high-priority
flag was detected under the configured screen`, not `good`, `normal`, or `no
problem`. Each card retains the applicable method and limitation note.

Before real-data upload, the application presents deployment-specific data
governance guidance covering processing location, retention/cache limits,
de-identification, consent or secondary-use authority, access control,
small-cell and intersectional re-identification risk, and data that must not be
uploaded to a public deployment. An acknowledgement records that guidance was
shown; it does not authorize processing or unlock a high-stakes claim.

#### Guidance architecture

Add `mfrm_app/guidance.py` as a pure-Python module that imports neither
Streamlit nor pandas. It owns stable node, goal, stage, completion-rule, and
target IDs plus a deterministic state reducer. It stores i18n keys, not display
strings. Streamlit adapters translate only at render time.

The authoritative model separates:

- learning state: `LOCKED`, `AVAILABLE`, `ACTIVE`, `COMPLETE`;
- workflow completion: whether the bounded route reached a terminal record;
- computation state: `AVAILABLE`, `CAUTION`, `HOLD`, `NOT_ASSESSABLE`;
- conclusion stability: `STABLE`, `CONDITIONALLY_STABLE`, `SENSITIVE`,
  `NOT_ASSESSED`; and
- claim readiness: a ClaimID-specific derivation from required evidence,
  context, applicability, review, and stability rather than a user toggle.

User acknowledgement may complete a learning node but cannot alter the latter
three scientific systems. Existing onboarding, first-run route, Guided, and
Help functions remain compatibility wrappers until their AppTests are
migrated. Routing by translated label or English substring is replaced by
explicit target IDs.

The first implementation consolidates route identity, prerequisites,
completion rules, targets, and locale keys. Long conceptual Help text remains
separate content projected from the same concept IDs; forcing every paragraph
into one catalog would create another hard-to-maintain monolith.

#### Invalidation and identity rules

| Change | Guidance effect | Analysis effect |
|---|---|---|
| Data or required mapping changes | Re-evaluate from `data_check` | New or stale AnalysisID |
| Model or estimator setting changes | Re-evaluate from `estimate` | New or stale AnalysisID |
| Analysis depth or computed diagnostic changes | Re-evaluate affected evidence | New evidence/analysis identity under the current contract |
| New fitted AnalysisID | Clear analysis-dependent review and archive completion | Bind to new result |
| Intended purpose, population, interpretation, use, or stakes change | Re-evaluate interpretation and decision nodes | Preserve AnalysisID unless the change also alters `AnalysisSpec`; create a new StudyContextID/InterpretationID/DecisionID |
| Planned or observed rating-design/scoring-process facts change | Re-evaluate `data_check` and affected claims | Preserve or invalidate AnalysisID according to whether data mappings/specification changed; always version the design/process record |
| View, guide, or language changes | Preserve guide validity | Preserve AnalysisID and the fitted/draft `AnalysisSpec` |
| Evidence card is opened | Record learning progress only | Preserve evidence state |
| Decision note changes | Update the decision artifact | Preserve AnalysisID |
| Pure display-theme change | Preserve guidance | Preserve statistical AnalysisID; update only artifact identity if needed |

A result whose data or settings no longer match the current context cannot
complete `estimate`, satisfy evidence review for the current context, or
produce a claim-ready archive. The guide returns to the first invalidated node
and explains what changed. It may still preserve an explicitly stale or
no-conclusion history record that cannot be mistaken for the current result.

#### Native analysis record

Add a versioned Python-native analysis-record contract after the study,
evidence, identity, and guidance contracts are stable. The primary user
artifact is a human-readable Markdown or HTML `MFRM_Analysis_Brief`; a
machine-readable `MFRM_Analysis_Record.json` is its deterministic sidecar.
Both derive from the same structured record and reference the same stable
StudyContextID, AnalysisID, EvidenceIDs, InterpretationID, DecisionID, and
ClaimIDs.

The structured record contains resolved settings, study/design/scoring-process
context, EvidenceRecords and unresolved evidence issues, supported and
unsupported claim entries, an explicit no-conclusion state when applicable,
the next action, authorship/review status, application/schema versions, and a
privacy declaration. Raw response rows and direct person/rater identifiers are
absent by default. Display language is artifact metadata; canonical ClaimIDs
and user-authored text are stored separately so translation does not silently
change a scientific identity.

The scientific record does not include tutorial progress, skipped-node
history, view preference, or external-engine compatibility metadata. Guide
state is not scientific evidence. Generating either artifact is distinct from
claiming that a browser download completed successfully, and the default
status is `draft_for_review`, not `validated` or `approved`.

#### Delivery slices

The UI migration is deliberately incremental. Accessibility and bilingual
meaning are acceptance criteria for every visible slice, not a final polish.

1. **G0 — Freeze contracts, behaviour, and interpretive risk.** Approve the
   StudyContext/Guidance/Identity contract, inventory every interpretation
   surface, classify unsafe, stale, or internal-only copy, map the current Help
   link graph and orphans, and add baseline tests for the built-in sample,
   AnalysisID, EvidenceRecords, stale detection, modes, locales, and exports
   without redesigning the UI.
2. **G1 — Canonical state and guidance core.** Add widget-independent
   `AnalysisSpec`, study/design/scoring-process records, the pure guidance
   catalog and reducer, the public-vocabulary map, `HelpTopic`/`HelpLink`
   registries, lifecycle and invalidation rules, locale-key parity tests, and
   compatibility wrappers. Replace G0 `STOP` copy before exposing a new
   evidence-linked route.
3. **G2 — Optional sample guide.** Replace the overlapping landing
   onboarding/Tutorial presentation with the five-step sample route, including
   skip, exit, resume, restart, formative interpretation, `HOLD`, stale, and
   non-destructive sample-context branches. Every contextual Help action opens
   the exact topic and can return to the originating node/card.
4. **G3 — View projection.** Map the existing Essential/Guided and
   Full/Advanced surfaces to Standard and Detailed views. Switching view must
   preserve fitted and draft specifications and invoke no estimator or hidden
   diagnostic.
5. **G4 — Evidence-to-decision artifacts.** Render the priority EvidenceIssue
   with its tensions and aggregate counts, collect bounded
   claim/limit/action or no-conclusion notes, and generate the human-readable
   Analysis Brief plus its privacy-safe JSON sidecar.
6. **G5 — Research and education companion.** Add an unknown-case transfer
   exercise, data-governance preflight, Methods/Results/Limitations reporting
   matrix, reviewer checklist, and versioned sample manifest without turning
   learning progress into scientific evidence or the application into an LMS.
7. **G6 — Consolidate and harden.** Move long theory to Help, remove only
   proven-unreachable duplicate guidance, complete the copy-audit remediation,
   bilingual AppTests, browser accessibility checks, and then retire
   compatibility wrappers.

Statistical changes, state-model extraction, visible UI changes, and legacy
code removal remain separate review units. Old routes are not removed until
new end-to-end gates pass.

The first G1/G2 vertical slice shipped on 2026-08-04: the pure guide catalog
and reducer, focused bilingual five-step sample journey, direct skip/exit,
resume/restart, production-estimator fit binding, formative overclaim guard,
and learning-only completion are implemented and AppTest-covered. G1 is not
complete: canonical `AnalysisSpec` and the remaining study/design/scoring
records still precede a full claim of canonical state. G2 is also not complete
until every node has exact contextual Help return, `HOLD`/stale/failure paths,
a complete existing-real-fit sample roundtrip, and browser keyboard/focus
acceptance.

#### Validation and evaluation

Unit tests cover unique IDs, valid routes and targets, completion/invalidation,
locale-key parity, one-time quick-start consumption, and the rule that user
acknowledgement cannot clear a machine `HOLD`. They also verify that view,
language, Help, learning answers, and guide lifecycle events cannot mutate the
fitted or draft `AnalysisSpec`, AnalysisID, EvidenceRecords, or claim readiness.

AppTest covers at least:

- Japanese and English sample completion;
- immediate skip by an experienced user;
- exit after estimation without refitting or losing AnalysisID;
- stale settings returning to the correct step;
- a `HOLD` path with a named repair and an archiveable no-conclusion record but
  no claim-ready continuation;
- a `CAUTION` path that carries its boundary into the decision note;
- language change at evidence review with the same AnalysisID and EvidenceID;
- resume and restart without deleting data or results;
- sample start/resume without overwriting an existing real-data context;
- Standard/Detailed projections displaying the same evidence state while
  retaining hidden advanced values and making zero estimator calls;
- independent `workflow_complete`, `learning_complete`, and `claim_ready`
  states;
- planned-absence and unexpected-missingness reason paths; and
- an evidence-tension fixture that cannot collapse into an all-clear message.
- persistent static Help access from landing/no-data, mapping/preflight,
  first-run failure, `HOLD`/`NOT_ASSESSABLE`, and post-fit states;
- one-action contextual Help routing to the exact topic and back to the same
  guide/EvidenceIssue/control without changing scientific or learning state;
- unknown, unavailable, and stale Help targets producing a visible safe
  fallback rather than a silent no-op; and
- absence of internal enum/class/file-pipeline language from the guide and
  Standard Help path, with optional technical identifiers still available in
  the declared technical layer.

Browser-level checks begin with the first visible slice and target WCAG 2.2 AA
in the supported Streamlit/browser scope. They cover keyboard-only completion
and exit, rerun focus/status announcement, logical headings and reading order,
320 CSS-pixel layouts and 200--400% zoom, long Japanese text, reduced motion,
touch-target size, non-colour state recognition, text/table alternatives for
figures, and keyboard-safe popovers/dialogs. AppTest alone is not treated as an
accessibility audit.

Initial evaluation uses deterministic tests and moderated observation rather
than production telemetry. Success is measured by whether a user can run the
sample, recover from a blocker, find Current focus, identify one supported and
one unsupported claim, and create the analysis record. A second, unfamiliar
single-problem SLA scenario measures transfer: the user must identify what is
supported, what is not supported, what to inspect next, and when supervisory
review is needed. Satisfaction or guide completion alone is not learning
evidence. If telemetry is ever
considered, it is opt-in, coarse, and excludes raw ratings, file/column names,
identifiers, fingerprints, AnalysisID, free text, and resolved settings.

#### Explicit non-goals for the initial guide

- A forced wizard or a separate top-level Tutorial tab.
- DOM-selector coach marks, automatic tutorial progression or selector-driven
  tour scrolling, hover-only instruction, or a video-first experience.
- A tour of every diagnostic, figure, export, and advanced setting.
- A second estimation or result pipeline for tutorial data.
- Universal threshold teaching, automatic validity claims, completion badges,
  or qualification certificates.
- Automatic high-stakes authorization, rater personnel decisions, or claims
  that a completed guide establishes research competence.
- A learning-management system, graded student accounts, or hidden instructor
  surveillance. Education worksheets and rubrics remain separable companion
  artifacts.
- Persistent browser tracking before a demonstrated need.
- TAM, R, Julia, external posterior, or cross-engine ZIP workflows.

Exact-source focus/scroll restoration after an explicit Help return is not
tutorial progression. It is permitted only for a registry-owned real heading,
must yield to user-moved focus, uses no smooth motion, and remains governed by
the contextual Help browser-acceptance runbook.

#### Completion gate

- Standard and Detailed views consume identical current EvidenceRecords and
  never display conflicting states for the same AnalysisID.
- The complete sample journey, immediate skip, exit at every step, resume,
  restart, stale invalidation, `HOLD`, and `CAUTION` paths pass in both locales.
- Skipping or switching view preserves data, fitted results,
  AnalysisIdentity, EvidenceRecords, and draft decision notes and does not
  trigger re-estimation.
- Every visible status has text, an interpretation boundary, a next action,
  and a method or limitation note; no critical meaning depends on colour.
- A user can identify one supported and one unsupported claim and can return
  from stale evidence to the setting that invalidated it without restarting
  the application.
- The Analysis Brief and JSON sidecar are schema-consistent, linked to the
  current study/analysis/interpretation/decision identities, privacy-safe by
  default, and contain no tutorial-progress or external-engine compatibility
  state.
- A failed, `HOLD`, or required-`NOT_ASSESSABLE` analysis can preserve a
  reviewable no-conclusion record without being presented as claim-ready.
- The default journey contains one primary action per guide screen, remains
  keyboard and narrow-screen usable, and does not eagerly render every heavy
  diagnostic panel.
- Every visible Help action reaches a registered bilingual topic in one action,
  retains a usable return target, and exposes claim-critical limits inline.
- Guide and Standard surfaces do not require users to understand internal
  class, enum, schema, manifest, runner, or raw ID names.

### Public-beta hardening

Objective: make release claims depend on native validation evidence and
operational reliability.

Deliverables:

- Replace external-parity language as a release gate with native numerical
  consistency, simulation calibration, and explicitly bounded method claims.
- Establish checked benchmark baselines for built-in scenarios and review
  material regressions in runtime or peak memory.
- Test supported Python versions in clean environments with pinned dependency
  bounds and deterministic fixtures.
- Exercise upload limits, malformed inputs, cancellation, stale state, cache
  isolation, privacy redaction, and report generation.
- Maintain a machine-readable release-scope summary covering supported models,
  estimators, diagnostics, sensitivities, and uncalibrated areas.
- Keep the public-beta label until the documented native validation matrix and
  release checklist are complete.

Completion gate:

- Compile, doctor, release-check, self-test, unit/integration, AppTest, privacy,
  localization, and benchmark gates all pass from a clean checkout.
- The release contains no confidential, machine-specific, or generated test
  data.
- README, in-app help, claim-to-evidence output, and release-scope summary agree
  on supported and unsupported functionality.
- No public claim depends on an external-engine run or an external result
  archive.

## Cross-cutting engineering rules

- `main` is the canonical public integration target. Experimental checkouts are
  donor workspaces only; pure-Python changes are selected by contract and test,
  not merged wholesale with rejected external-handoff work.
- New statistical computation belongs in `mfrm_app/`; Streamlit rendering
  calls it through a narrow contract.
- Canonical `AnalysisSpec`, StudyContext, rating-design, scoring-process,
  evidence, interpretation, decision, and guidance state live outside
  Streamlit widgets. Views are projections over those records, not alternate
  sources of truth.
- Refactoring and statistical changes are separate review units.
- Public tables are versioned contracts. Column removal or semantic change
  requires migration documentation and tests.
- Fixed seeds are part of the analysis identity, not hidden implementation
  details.
- Cache keys include all data, model, diagnostic, sensitivity, and code inputs
  that can change a result.
- Raw uploaded data and direct identifiers are not written to disk or included
  in reports unless the user explicitly requests an allowed export.
- Real-data entry states the deployment-specific processing and retention
  boundary before upload. Education or support artifacts exclude raw rows,
  column names, identifiers, fingerprints, free text, and person/rater labels
  by default.
- Every warning has a stable reason code, user-facing explanation, and
  recommended action.
- Every automatic recommendation states what it repairs and what it does not
  establish.
- No threshold is promoted from a descriptive convention to a universal
  acceptance rule without calibration and documentation.
- No help text, chart guide, popover, caption, success message, export prose,
  or manuscript template may express a stronger claim than the governing
  EvidenceRecord and StudyContext permit.
- Statistical rater or subgroup flags are investigation prompts. Automated
  exclusion, retraining, punishment, fairness verdicts, or high-stakes use are
  outside the core decision contract.

## Prioritized work inventory

One serial queue would postpone user, interpretation, and accessibility risk
until after most statistical work. Delivery therefore proceeds through three
parallel tracks with a shared P0 gate. Sequence is strict within a track, not
between tracks. A visible feature may ship only when its applicable scientific,
experience, and research-practice gates all pass.

### Shared P0 gate

| Outcome | Prerequisite | Principal evidence of completion |
|---|---|---|
| Canonical scope, evidence, and identity contracts | None | Boundary and contract tests |
| StudyContext, rating-design, scoring-process, guidance, and `AnalysisSpec` contracts | Scope/identity decisions | Approved schemas, transition/invalidation matrix, and no widget ownership |
| Full interpretation-copy inventory and risk classification | Evidence vocabulary | Every locale/help/chart/export surface has an owner and gate |
| Audience terminology and Help connection graph | Guidance/claim-boundary decisions | Surface-exposure matrix plus zero-orphan topic/link audit |
| Priority user tasks and SLA research scenarios | Study-context contract | Scenario matrix covering planning, first analysis, reproduction, teaching, and review |

### Track A — Scientific evidence

| Sequence | Outcome | Prerequisite | Principal evidence of completion |
|---:|---|---|---|
| A1 | Modular network and decision computation | Evidence contracts | Compatibility fixtures |
| A2 | Design perturbation and fragility audit | Modular network layer | Known-graph tests |
| A3 | Evidence-coherence computation | Normalized evidence records | Conflict fixtures |
| A4 | Conclusion-level sensitivity ledger | Analysis identity contract | Stability-state fixtures |
| A5 | Design repair candidates | Fragility audit | Feasibility and burden tests |
| A6 | Modular simulation condition system | Analysis identity contract | Serialization/fingerprint tests |
| A7 | Facet-aware stress scenarios | Simulation conditions | Recovery and flag-rate summaries |
| A8 | Native validation matrix | Frozen scenarios | Calibration reports with MCSE |

### Track B — Guided experience

| Sequence | Outcome | Prerequisite | Principal evidence of completion |
|---:|---|---|---|
| B1 | Canonical guidance state, audience vocabulary, Help registry, and safe-copy remediation | Shared P0 | Reducer, terminology/help-link/copy contracts, and locale tests |
| B2 | Optional five-node sample route with contextual Help return | Guidance and isolated tutorial context | Bilingual skip/exit/resume/stale/formative/Help AppTest journeys |
| B3 | Standard/Detailed view projections | Widget-independent `AnalysisSpec` | State-consistency and zero-computation view-switch tests |
| B4 | Current EvidenceIssue and recovery actions | Evidence-coherence computation | Conflict, aggregate-count, and repair-routing fixtures |
| B5 | Human-readable Analysis Brief and JSON sidecar | Interpretation/decision identities | Schema, privacy, language, and deterministic-generation tests |
| B6 | Accessibility and legacy-route consolidation | Proven replacement journeys | Browser matrix and unreachable-route evidence |

### Track C — SLA research and learning practice

| Sequence | Outcome | Prerequisite | Principal evidence of completion |
|---:|---|---|---|
| C1 | Study/rating/scoring records and RQ-to-design preflight | Shared P0 | Holistic/analytic and crossed/nested scenario tests |
| C2 | Missingness provenance and data-governance preflight | C1 | Reason-code, redaction, and hosted/local boundary tests |
| C3 | Claim-to-evidence and Methods/Results/Limitations matrix | Evidence and interpretation contracts | Reviewer fixtures and bounded prose tests |
| C4 | Formative sample, unfamiliar transfer case, and education companion | B2 and C3 | Learner evidence, versioned worksheet/rubric, and no-state-mutation tests |
| C5 | Privacy-safe support/reviewer bundle | B5 and C3 | Blind handoff and disclosure tests |

Public-beta hardening is a cross-track release gate. It does not require every
future item in all tracks, but every advertised claim and journey must have its
scientific validation, safe interpretation, accessibility, privacy, and clean-
environment evidence complete.

## Research validation matrix

The validation program will cross the conditions below selectively rather than
claiming that one giant factorial experiment is always feasible.

| Dimension | Representative conditions |
|---|---|
| Response model | RSM, PCM, bounded GPCM within supported scope |
| Estimator | JMLE; supported MML engine choices; fixed/free population SD where implemented |
| Performance mode | Synthetic writing, speaking, and human-scored short response within declared applicability |
| Rubric structure | Holistic; analytic criteria intended as one construct; analytic subconstructs/composites that trigger bounded or unsupported claims |
| Assignment | Fully crossed, balanced incomplete, chain, bridge-dependent, disconnected |
| Facet relation | Crossed, nested, anchored, fixed/random role declarations, and partial confounding |
| Scoring operation | Single, double, adjudicated, rescored, common-performance, and monitoring designs |
| Rater structure | Few/many levels, balanced/imbalanced workload, narrow/wide severity, session/order/drift, training cohort |
| Score support | Well populated, sparse middle category, unused declared category, extreme concentration |
| Missingness | None, planned incomplete assignment, nonresponse/not-reached, rater skip/unratable, technical/administrative loss, rater/task/subgroup-concentrated loss |
| Dependence | None, supported halo/local interaction, criterion/performance clustering, occasion drift |
| Differential interaction | Rater x task/criterion/proficiency/subgroup with overlap, multiplicity, sparse cells, and small-cell suppression |
| Person distribution | Reference normal and documented non-normal sensitivity conditions |
| Sample information | Person count, observations per level, and pairwise overlap varied separately |
| Resampling unit | Person, performance, task, rater, and justified clusters rather than rating-row bootstrap by default |

For each applicable cell, the report must identify estimand, requested and
successful replication counts, numerical failures, Monte Carlo uncertainty,
and the exact claim the result can or cannot support.

## Release decision rules

A feature is complete only when all applicable items below are satisfied:

- The statistical question and non-goals are documented.
- Applicability checks fail closed with stable reason codes.
- Deterministic unit fixtures and at least one end-to-end scenario exist.
- Sensitivity to missingness or sparse overlap is tested when relevant.
- UI, export, narrative, and evidence-ledger outputs use the same result.
- Japanese and English interfaces preserve the same statistical meaning.
- Privacy and cache behavior are reviewed.
- Performance is measured against a checked baseline.
- README/help/claim boundaries are updated.
- Current Help, tutorial, chart guides, popovers, captions, success messages,
  exported prose, and manuscript templates do not contradict or exceed the
  governing EvidenceRecord and StudyContext claim boundary.
- Every contextual Help control reaches a registered bilingual topic/section
  and returns to its source; stale or mismatched analysis context cannot appear
  as a current dynamic explanation.
- Guide and Standard surfaces use the approved audience vocabulary. Raw
  enum/class/schema/backend/environment/runner terms and IDs appear only in
  declared technical details or machine-readable artifacts.
- Production errors reveal no raw exception text, backend symbol, environment
  override, uploaded value, or identifier; they provide a stable user reason,
  reversible action, and privacy-safe support reference.
- Optional guidance can be skipped, exited, resumed, and restarted without
  changing or deleting the current analysis state.
- Workflow completion, learning completion, and ClaimID-specific readiness are
  tested as independent states; acknowledgement cannot unlock scientific or
  high-stakes use.
- Pure presentation changes, including language and view, do not alter
  AnalysisID or trigger re-estimation.
- Standard/Detailed switching preserves all hidden draft and fitted settings
  and invokes no estimator or newly requested diagnostic.
- A `HOLD` or required `NOT_ASSESSABLE` state can be archived honestly as a
  no-conclusion record but cannot produce an unbounded supported claim.
- Real-data upload and exported/support artifacts pass the applicable
  data-governance, privacy, redaction, and small-cell checks.
- Statistical flags about a person, rater, or subgroup never become an
  automated punitive or fairness decision.
- The full relevant test gate passes without suppressing unrelated failures.

## Deferred questions

The following require explicit product decisions and are not silently assumed
by this roadmap:

- Which legacy external-comparison and advanced-model handoff surfaces beyond
  the selected ConQuest/TAM plan should be deprecated or retained.
- Whether the large Streamlit entrypoint should ultimately become a thin shell
  or retain some stable compatibility wrappers.
- What native evidence would be required before considering a confirmatory
  multidimensional Python model.
- Whether guide-completion preference should persist beyond one session; this
  requires observed user need and a separate privacy/accessibility decision.
- Whether an explicit, privacy-safe learning checkpoint should be exported and
  imported for multi-session teaching before any browser persistence is added.
- Whether `Interpret` should later become an eighth top-level navigation stage;
  the InterpretationRecord remains separate even while it is grouped inside
  `Decide`.
- Which domain profiles beyond the initial synthetic SLA cases deserve
  versioned applicability rules rather than generic prose.
- When the legacy seven-section navigation and compatibility wrappers can be
  removed after Standard/Detailed view migration.
- What calibrated criteria, if any, justify removing the public-beta label.

Until those decisions are made, they cannot block the standalone core or be
used to expand its claims.

## Known assignment mechanism qualification (2026-08-11)

The repository now has a qualified simulation primitive for
`Pr(G | fixed degrees, connected) proportional to exp(gamma*T(G))`. A 90-state
exact oracle verified direction, normalization, finite-chain recovery, degree
margins, connectivity, and score-free materialization. This closes the first of
four deliberately separate layers:

1. assignment-kernel correctness — passed for the frozen finite-state oracle;
2. large-design mixing and dose calibration — next, with multiple chains and
   fixed Person/Rater exposure;
3. response-generation and within-estimator bias sensitivity — requires a new
   frozen PCM study;
4. FACETS/Python workbench qualification — FACETS versus Python only on the
   JMLE estimand, with MML/CMLE reported as different-estimand sensitivities.

The mechanism remains repository-only until all downstream gates pass. It is
not a fitted MAR/MNAR model, and its known `gamma` must not be presented as an
empirical estimate. See
`validation/KNOWN_ASSIGNMENT_MECHANISM_PILOT_20260811.md`.

The 80×4 large-design phase subsequently retained two failed finite-MCMC
settings, qualified a four-Rater exact-DP sampler, and selected `|gamma|=0.8`
without observing responses. A separate fixed-Person response screening then
completed 120/120 attempts and 30/30 FACETS/Python JMLE calibration pairs. It
identified a candidate, domain-selective MML Rater recovery sensitivity, while
all JMLE/CMLE and non-Rater screening intervals included zero. Therefore layers
1–4 are operationally connected for this narrow DGM, but public UI promotion
remains blocked.

The next evidence step is a fresh, prospectively fixed confirmatory design with
multiple independent Person vectors or population draws. It should register
the MML Rater-RMSE contrast and severity compression/expansion estimand before
generation, retain gamma dose-response, and treat FACETS solely as the JMLE
engine calibration. Until then, no robustness, MAR/MNAR, or estimator-ranking
claim is allowed. See
`validation/KNOWN_ASSIGNMENT_LARGE_DESIGN_CALIBRATION_20260811.md` and
`validation/KNOWN_ASSIGNMENT_RESPONSE_SCREENING10_20260811.md`.

The computational prerequisite for that step is now closed. A prospectively
gated dense-margin exact DP reproduced the 90-state oracle and all five
retained 80-by-4 recursive log normalizers (maximum difference `2.84e-14`).
The largest partition build took 0.202 seconds versus the retained 92.68-second
same-machine recursive diagnostic at `gamma=+0.8`. The route remains bounded
by an explicit dense-cell cap and does not itself qualify any estimator. Before
the confirmatory design is frozen, use it in a fresh multi-Person-vector
preflight and keep those observations out of the confirmatory sample. See
`validation/KNOWN_ASSIGNMENT_DENSE_DP_QUALIFICATION_20260811.md`.

The registered four-vector preflight completed all 48 attempt units and the
separate visible-mode FACETS supplement passed 12/12 same-input calibrations.
The derivative aggregate passed every operational gate, and PF1--PF4 all met
their registered n=4 direction/count rules. R 4.5.1 independently reconstructed
the values. The next valid step is a separately frozen confirmation using new
Person vectors; none of these four vectors or earlier screening observations
may be pooled into it. See
`validation/KNOWN_ASSIGNMENT_MULTIVECTOR_PREFLIGHT_STATUS_20260811.md`.
