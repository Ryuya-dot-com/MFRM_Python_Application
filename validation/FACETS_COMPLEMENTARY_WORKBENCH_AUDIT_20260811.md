# FACETS 4.5.0 補完ワークベンチ監査

## 2026-08-11 運用訂正とmulti-vector追補

後続診断により、Windows上のFACETS 4.5.0呼出し契約を訂正した。hidden
`BATCH=YES`だけでなく、`BATCH=NO`をWindows `STARTUPINFO`でhiddenにした経路も
`XojoGUIFramework64.dll`内で `0xC0000005` となる。一方、direct/Pythonのvisible
`BATCH=NO`は完走した。現在のqualified adapterはwindowをvisibleに保ち、reportと
全facet Scorefileの非空・安定を確認してから `WM_CLOSE` する。強制終了は成功では
ない。また、補助Scorefileが260文字へ達すると停止したため、共通呼出しは220文字の
保守的path budgetを起動前に強制する。

この訂正経路で、既知の確率的割付multi-vector preflightの同一12入力を別保存し、
FACETS/Python JMLE校正12/12、zero retry、facet内Spearman最小1.0を得た。元の48
attemptと失敗artifactは変更せず、native Python JMLE replayが最大 `4.44e-16` で
一致した後だけderivative aggregateを作成した。登録済みPF1--PF4は全て4/4で登録
方向を示し、R 4.5.1が16値と4要約を最大 `1.67e-16` で再構成した。これは別確認
計画へ進むためのn=4 preflightであり、確認、一般的バイアス定数、推定法ランキング
ではない。通常2桁のFACETS fit表示は未丸め計算値として使っていない。

実施日: 2026-08-11  
対象: Streamlit MFRM application / `C:\Facets\Facets.exe` 4.5.0  
結論の強度: アダプタ資格確認・小標本スモーク。推定法の優劣を一般化する検証ではない。

## 結論

本アプリの妥当な長期的位置づけは、FACETSの代替品ではなく、**FACETSで外部照合された補完的MFRMワークベンチ**である。

現在、登録した加法RSM/JMLEに加え、四カテゴリPCM/JMLE（Criterionを単一step facetとする範囲）もFACETS 4.5.0との外部数値照合が成立した。PCMではcomplete/planned-connectedの独立100反復まで進み、200/200ペアが直接一致・校ブレーションgateを通過した。同時に、MML（母集団SD固定・自由推定）とexact CMLEを同じ生成データへ通し、推定法の優劣ではなく、Person母集団仮定と観測設計が有限標本結果をどう変えるかを監査できる基盤ができた。

ただし「FACETSで検証済み」と製品全体へ一般化してはならない。外部照合は、登録済み符号・中心化・アンカー規約、適格な収束、構造的に識別された加法RSMと限定PCMに限られる。group anchor、複数scale/model、混合category support、fitの生値一致、biasの安定したfalse-positive/power一致は未検証または非対応である。FACETSの丸め表示とPythonの生値も同一精度として扱わない。

## 今回確立した証拠

### FACETS対Python JMLE

- FACETS 4.5.0は16/16反復で返却・収束した。
- 構造的に識別された12/16反復、108主ファセット母数を直接比較対象とした。
- 尺度調整後のFACETS対Python JMLEのMAEは **0.00173366 logits**、最大絶対差は **0.0104228 logits** だった。
- ファセット内順位のSpearman相関は、直接比較した全群で1.0だった。
- 疎な一人一評定者条件4反復は、観測デザインの構造的nullityが7であり直接比較から除外した。FACETSの `Subset connection O.K.` は、このアプリの線形予測子全体がfull rankであることを保証しない。

### アンカー制約の発見と修正

初回照合では、FACETSは指定アンカーで原点を固定した後に未アンカーレベルを自由推定する一方、Python実装は未アンカーレベルだけを再中心化していた。この過剰制約により、アンカー条件のRater MAEは最大約0.23 logits、最大差0.56 logitsまで増大した。

`build_facet_constraint` を修正し、hard/group anchorがfacet原点を固定した場合は未アンカーレベルを独立座標として扱うようにした。修正後、正しいアンカー条件でFACETS対PythonのRater MAEは約0.00077–0.00103 logits、アンカードリフト条件でも約0.00086 logitsとなった。これは単なる丸め差ではなく、外部照合が実装上の識別制約バグを検出した例である。

### ネイティブPython MML

同じ16反復を以下の2モードへ通した。

- `PYTHON_MML_FIXED_SD1_Q31`: 生成時のPerson SD=1を固定
- `PYTHON_MML_FREE_SD_Q31`: Person SDをEMで自由推定

32/32適合が返却・収束・主ファセット有限性の適格条件を満たした。完全デザインでは自由SD平均が約1.05–1.15だった一方、疎な一人一評定者デザインでは平均 **0.485** へ縮小した。疎条件のRater MAEは固定SDで **0.868**、自由SDで **0.680** logitsだった。

これは自由SDが「より正しい」と示す結果ではない。疎条件ではPerson分散とRater差の情報分離が弱く、正規母集団仮定と分散推定が結果を大きく規定している。MMLでは推定が返っても、設計だけで識別されたJMLEと同じ意味ではない。そのため疎条件は `normal_population_distribution` による仮定依存感度分析として記録した。

意図的なアンカードリフト条件では、固定／自由MMLのRater MAEは約0.277／0.280 logitsだった。アンカーfacetを再中心化せず絶対尺度で評価しているため、アンカー汚染が見かけ上消えていない。

### 独立エンジン

- Python exact CMLEとR `immer` の自由座標は、12共同返却反復で平均絶対差 `9.36e-9`、最大差 `5.55e-8` だった。
- ネイティブfree-SD MML Q31とTAM MML Q61は、非アンカー完全条件と疎条件で多くの主ファセットが約0.0002–0.0013 logitsの平均差だった。ただしアンカー伝達規約が異なる条件では差が拡大したため、同一推定量の証明とは扱わない。
- sirtとのRater順位は完全条件でよく一致したが、疎条件では差が増えた。
- `mfrmr` 0.2.3 strict JMLEは非アンカー完全条件でPython JMLEとの平均差が約0.000013 logitsだった。アンカー条件では制約意味の差が拡大した。使用したmfrmrは隔離インストール済みsnapshotで、現在のsource treeを回収できていないため、補助証拠に限定する。

このMML/R独立エンジンlaneの反復数は各セル2である。上記はパイプライン、尺度、制約、明白な故障を確認する値であり、バイアス、被覆率、検出力、推定法優越性の安定した推定値ではない。FACETS/JMLE/Table 8/13 laneは、後述の各セル20反復pilotまで進めた。

### 20反復pilotの実行会計

- 事前登録済み8条件×20反復のimmutable bundleを用いた。FACETSは160/160返却、158/160収束。118反復を直接JMLE比較へ入れ、構造非識別の疎40反復と非収束2反復を除外した。
- 直接比較した1,053主ファセット母数のFACETS対現在Python JMLE weighted MAEは **0.001535 logits**、最大絶対差は **0.012015 logits**、条件×facet群の最小Spearman相関は **0.999833** だった。
- 現在Python完全診断refitは160/160返却、159/160収束、119/160 analysis-eligible、記録fit時間337.9秒だった。FACETSの主・Table 8補助のdual passはwall time 328.8秒だった。
- Windows refitのCSV byteは改行・文字列表現により旧macOS bundleと異なるが、manifest 160行、rating 185,184行、truth 1,600行、anchor 160行の値、および160/160のDataId/FitInputIdは完全一致した。FACETSとPython step adapterには旧immutable byteを直接渡し、現在Python比較結果は別ディレクトリ・別hashとして記録した。
- Python pilotは全成果物を書いた後、最後のCP932コンソール表示だけで `UnicodeEncodeError` になった。UTF-8成果物とbridge検証は完了しており、CLIは表示不能文字を安全にescapeして計算完了を偽失敗にしないよう修正した。

### Table 8 category/threshold

- FACETS 4.5.0のTable 8を160/160反復で解析し、640/640カテゴリ行のTotal countがimmutable生成bundleから再集計した度数と一致した。
- `Umean=0,1,6` はmeasureを高精度で出す一方、固定幅のcategory Outfit欄が `1.` や `.` へ切れるため、その文字列から値を推測しない。主passのthreshold/measureと、同一データを `Umean=0,1,2` で再報告した補助passの一桁category Outfitを結合するdual-pass契約にした。
- category Outfitは表示値±0.05として扱い、記述的0.5–1.5帯に633行がpass、5行が `boundary_uncertain`、2行がflagだった。17行はlow-countだった。これはカテゴリ機能の証明ではない。
- 構造的・収束適格な118反復の351 RSM thresholdで、FACETS対Python JMLEのweighted MAEは **0.004789 logits**、最大絶対差は **0.020472 logits** だった。大標本・アンカー条件の条件別MAEは約0.00162–0.00264、小標本では約0.00993–0.01018だった。
- 事前計画にthresholdの等価許容幅を置いていないため、これを「完全一致」へ事後的に昇格させない。現在の表現は、登録加法RSMに限定したformula/numerical alignmentである。

### Table 13 Rater×Task bias

- 三つの並べ替え出力を重複集計せず、明示的に `arranged by N` とされた一つをcanonical Table 13として160/160反復で解析した。
- 独立生成データから期待した2,400セルのうちFACETSは2,292セルを報告した。未報告108セルはすべて疎条件のextreme/unmeasurable候補である。監査ledgerに残し、negative findingへの変換を禁止した。
- FACETS表示p値へrun内Holm補正、`|bias| >= 0.50` の実用閾値、`n < 5` のsparse no-claim、表示丸め区間を適用した。全報告セルでは1,492 no-flag、792 sparse no-claim、8 flagだった。
- FACETS/Python双方のrun適格性を要求した1,404セルで、strong判定は **1,404/1,404一致**し、8 flagを含む非退化なmappingになった。bias sizeの平均絶対差は **0.001412 logits**、最大絶対差は **0.009200 logits** だった。
- これはRater×Task、加法RSM、20反復pilotに限定した外部decision mappingであり、confirmatoryなfalse-positive/power一致ではない。

### 外部照合がbias結論を変えた例

アンカー制約修正前の旧Python pilotでは、0.6-logit条件のfocal strong flagがanchor driftで18/20、clean anchorで11/20だった。現在のabsolute-origin制約ではそれぞれ3/20、3/20となり、FACETSの適格focal判断とも一致した。非アンカー小標本条件は旧新とも1/20だった。

したがって、旧結果を「アンカー条件ではbias検出力が高い」と解釈していれば誤りだった可能性が高い。原因は推定法の性能ではなく、アンカーで原点を固定した後にも未アンカーレベルを再中心化した過剰制約である。この例は、simulation rateを見る前に外部ソフトで制約意味を照合すべき理由を具体的に示す。

## FACETSの丸めに対する契約

FACETS 4.5.0の出力を実機で確認した。

- `Umean=0,1,6` とtab Scorefileにより、Measure、S.E.、Displacementは6桁で保持する。
- Infit/Outfit MnSqとZSTDはScorefile、Residualfileとも小数点以下2桁の表示値である。
- Table 8の固定幅category Outfitは `Umean=...,6` で文字が切れるため値を推測せず、同一入力・仕様の `Umean=0,1,2` 補助passだけから取得する。threshold/measureは主passを正本とする。
- FACETS表示値をfitの生値として扱わない。表示値 `x` の可能な生値を保守的に `[x-0.005, x+0.005]` とする。
- 閾値の片側に区間全体がある場合だけpass/flagを確定し、閾値と重なる場合は `boundary_uncertain` とする。
- Python側は丸め前の生値を保存し、別にFACETS 2桁表示への射影を持つ。Pythonの生値を丸め済み値で置換しない。
- 丸め後に一致しても、fit式の同一性を証明したとは結論しない。丸め済みResidualfile成分からFACETSの「生fit」を再構成することもしない。

FACETSではmeasurement・reliability・agreementはTable 7、rating scale/categoryはTable 8である。アプリ内に残っていた「Table 8-style agreement」という対応付けは修正した。

## 現在のアプリとFACETSの差

| 領域 | 現在の位置づけ | 解釈 |
|---|---|---|
| 加法RSM JMLE | 限定範囲で外部照合済み | FACETSとの直接比較対象 |
| PCM JMLE | 四カテゴリ・単一step facet・complete/planned-connectedで外部照合済み | 任意のmissingness、複数scale、mixed supportへ一般化しない |
| GPCM/discrimination | アプリ独自拡張 | FACETS互換性を主張しない |
| MML | 固定/free SD、EM/direct/hybrid/auto、EAP、latent regression | 本アプリの主要な補完価値 |
| exact CMLE | repository core、immer照合 | FACETS equalityではなく異なるestimand |
| 情報的割付感度 | 等密度・等Rater露出でability–severity割付を摂動し、MMLのRater圧縮とPerson SD移動をfresh-100で確認 | 実データのMAR/MNAR診断・補正ではない |
| hard anchor | 加法RSMのelement anchorを外部照合 | group anchorと全推定法の意味一致は未完 |
| fit | Pythonは生値、FACETS textは2桁表示 | 区間・三値判断で比較 |
| rating scale | Table 8 parser、度数監査、RSM/限定PCM threshold照合済み | 複数scale・mixed support・生category fitは未検証 |
| local bias | canonical Table 13、Holm・実用閾値・sparse・欠落セル監査を実装 | 退化no-flag mappingのみ。false-positive/power一致は未検証 |
| 複数model/scale | 一分析一model family | FACETSの複数model statementには未対応 |
| capacity | Streamlitのpreflight上限あり | FACETS級の大規模処理能力を暗示しない |
| 再現性 | 同一byte bundle、hash、runtime identity、失敗台帳 | 補完ワークベンチとして強い領域 |

## メタ認知的なリスク点検

1. **推定法名を比較単位にしない。** JMLE/MML/CMLEというラベルだけでなく、尤度、Personの扱い、母集団分布、quadrature、中心化、anchor、extreme score、SEの定義を一組のestimand contractとして比較する。
2. **収束を妥当性と同一視しない。** 今回、疎条件のMMLは全適合が収束したが、自由SDは0.49付近へ縮小した。収束は計算成功であり、識別・頑健性の証明ではない。
3. **FACETSを生成器と審判に兼用しない。** 主研究は独立Python生成器を使い、全エンジンへ同じbyteを渡す。FACETS `Simul=` は二次的parametric bootstrapだけに用いる。
4. **煙試験を性能研究へ昇格させない。** 各セル2反復はパーサ、尺度、制約、失敗処理の資格確認にすぎない。
5. **モデルが正しい場合だけを調べない。** 歪んだ／混合Person分布、PCMをRSMでfit、rater slope、severity-dependent missingness、person×rater依存を登録因子に含める。
6. **UIの便利さで証拠階層を隠さない。** 「FACETS matched」「formula aligned」「analogous」「workbench addition」「unsupported」「unvalidated」を画面・export・レポートで明示する。
7. **対話処理と出版用診断を分ける。** 32 MML完全診断に約7分45秒を要した。大規模MMLでは点推定、標準診断、出版用不確実性を深度別にし、cache/background jobを使う必要がある。
8. **欠測率と割付機序を同一視しない。** 等しい50%観測・Rater露出・接続性でも、ability–severity割付はMMLのRater回収を変えた。実データでは真のThetaが未知なので、推定Thetaとの相関だけからMNARを宣言しない。

## PythonかRか

全てをどちらか一方へ寄せるより、役割を固定したハイブリッドが最適である。

- Python: 独立生成器、Streamlit製品、主要JMLE/MML/CMLE実装、実験orchestration、hash・schema・失敗台帳。
- FACETS 4.5.0: 外部JMLE標準、Table 7/8/13の参照、Windows batch golden master。
- R: TAM、sirt、immer、mfrmrによる独立感度分析。Python実装と同じコード経路を共有しないこと自体が価値になる。

全Python化は独立検証を弱め、全R化は製品runtimeとの距離を増やす。Rを「正解生成器」ではなく独立な反証可能性の層として使う。

## 長期ロードマップ

### P0: 「FACETS-validated」の境界を固定

- 現在の加法RSM golden bundle、FACETS spec/report/Scorefile、入力hash、実行ファイルhashをread-only evidenceとしてversion化する。
- **完了:** Table 8 parser、カテゴリ度数・平均measure・threshold・順序・dual-pass表示精度を対応付けた。
- **完了:** Table 13 parser、canonical arrangement、Holm・実用閾値・sparse・欠落セルの三値／no-claim mappingを追加した。ただし十分な反復研究まではbiasのfalse-positive/power「FACETS一致」を主張しない。
- **完了:** 実行時に構造的rank/nullityを計算し、外部ledgerがない場合に直接比較をfail closedする。
- **完了:** JMLE/MML/CMLEへestimand classとlikelihood basisを付け、cross-basis likelihood/AICと推定法ランキングを作らないsame-row bridgeを構成した。
- **完了:** FACETS report障害とPython統計証拠を分離し、exit-0/missing-reportだけの有限retry、cross-process lock、全try保持、依存関係hashを実装した。
- Streamlitの結果画面とexportへ互換性状態、エンジンidentity、丸め境界を表示する。

### P1: 20反復/セルのpipeline pilot（FACETS/JMLE core完了）

- **完了:** 事前登録済みbalanced coreとstress contrastの160反復を現在Python JMLEとFACETS dual passへ通した。
- **完了:** runtime、failure、構造非識別、Table 8丸め境界、Table 13欠落・非退化判定を記録した。
- **完了:** 歪み／混合／heavy-tail Person分布、complete/planned-connected missingness、anchor contamination、PCM-as-RSMを登録因子として実行した。
- **完了:** fresh-100でfree-SD MMLのshape感度を確認し、別のfresh-100 H6Rでplanned-connectedのJMLE threshold RMSE増加を独立再現した。
- **完了:** 等密度・等Rater露出のlatent ability–Rater severity割付を20反復でscreeningし、fresh-100で自由/固定SD MMLのRater RMSE増加、Rater圧縮傾き、自由MML Person SD低下を確認した。
- **未完:** rater slope、response-level person×rater依存、割付強度のfactorial化、既知propensity下の補正法を追加する。
- MML完全診断は深度別にし、全反復では点推定と最小収束情報を保存する。対話UIへ長時間診断を同期実行で持ち込まない。

### P2: 500反復/セルの登録研究

- 因子表、estimand、除外規則、失敗会計、decision threshold、Monte Carlo SEを実行前に凍結する。
- bias、RMSE、coverage、local-bias false-positive/power、fit判断不一致を条件別に推定する。
- 次の優先因子は、確認済みability–severity割付の強度・positivity・Rater数/Person数をfactorial化し、rater slopeとresponse-level person×rater依存を追加することである。今回の一つの構成機序を、欠測機構一般の主張へ拡張しない。
- 推定法の「勝者」ではなく、条件×損失関数ごとのtrade-offを示す。
- 登録研究と外部reviewが終わるまで、アプリが自動で推定法を推奨しない。

## 主要成果物

- 最終FACETS Table 7/8/13 pilot結果: `validation/facets_450_pilot20_table8_table13_v2_20260811/`
- 現在Python 20反復refit: `validation/operating_characteristics_pilot20_anchor_fix_20260811/`
- 同一byteの入力・R/Python結果: `validation/operating_characteristics_anchor_fix_smoke_20260811/`
- FACETS adapter: `validation/operating_characteristics_facets.py`
- native MML adapter: `validation/operating_characteristics_python_mml.py`
- native JMLE step adapter: `validation/operating_characteristics_python_table8.py`
- 互換性マトリクス: `validation/facets_compatibility_matrix.json`
- 事前計画: `validation/facets_complementary_workbench_plan_20260810.json`
- Table 8/13事前計画・事後追記: `validation/facets_table8_table13_parser_plan_20260811.json`, `validation/facets_table8_table13_parser_amendment_20260811.json`
- 20反復pilot事後追記: `validation/facets_pilot20_amendment_20260811.json`
- PCM資格確認: `validation/FACETS_PCM_QUALIFICATION_20260811.md`
- same-row estimand bridge: `validation/ESTIMAND_BRIDGE_PILOT_20260811.md`
- Person-shape fresh-100: `validation/ESTIMAND_DISTRIBUTION_CONFIRMATORY100_20260811.md`
- H6独立再現: `validation/H6_DESIGN_REPLICATION100_20260811.md`
- 情報的割付screening: `validation/INFORMATIVE_ASSIGNMENT_SCREENING20_20260811.md`
- 情報的割付fresh-100確認: `validation/INFORMATIVE_ASSIGNMENT_CONFIRMATORY100_20260811.md`

この監査後の最も重要な未完了事項は、確認済み情報的割付効果のfactorial一般化と実データ用positivity/selection-model契約、rater slope、複数scale/mixed support、fit判断のoperating characteristics、そしてUI/exportでの証拠階層表示である。現時点で証拠に見合う最も強い表現は、「限定した加法RSM/JMLEと四カテゴリPCM/JMLEをFACETS 4.5.0へ外部照合でき、丸め・構造非識別・report I/Oを失敗台帳に保持しながら、JMLE/MML/CMLEのestimand感度、観測密度感度、情報的割付感度を同一生成データで監査できる補完ワークベンチ基盤が成立した」である。
