# Estimand／情報的割付：Streamlit製品統合点検

実施日: 2026-08-11  
対象: MFRM Streamlit application  
証拠版: `informative_assignment_confirmatory100_20260811`

## 結論

FACETS補完ワークベンチとして、次の順序を製品に固定した。

1. 推定値を見る前に、現在の推定法が対象とするestimandと尤度基盤を示す。
2. 得点・推定能力・推定重症度・残差・fitを使わず、Person×Rater割付を記述する。
3. 観測割付だけではMCAR／MAR／MNAR／情報的割付を同定できないと明示する。
4. 同一密度の反実仮想感度分析は、同一文脈blockなら厳格なV1 dose経路、異なる文脈blockなら
   Rater×文脈margin厳密保存MILP endpointとして、ゲートを通った現在の適合結果に限り
   利用者が明示的に実行できる。未実行・不適格・部分失敗を区別する。
5. 外部fresh-100の効果量を、利用者データの補正式・prior・MNAR判定に転用しない。

この順序により、「接続している」「曝露が均等である」「欠測率が低い」という記述的
事実を、「割付が無視可能である」「推定にバイアスがない」という結論へ昇格させない。

## Estimand契約

画面と一括出力は、以下をランキングではなく契約表として並べる。

| Method | Estimand | Personの扱い | FACETSとの関係 |
|---|---|---|---|
| JMLE | fixed-Person joint calibration | Personをincidental fixed parameterとして同時推定 | 適格化済みRSM/PCM JMLE範囲の外部比較対象 |
| MML | Gaussian-population marginal calibration | Personを母集団分布について積分し、EAP/posteriorを報告 | FACETS JMLEはmarginal estimandの正解基準ではない |
| exact CMLE | Person-total conditional facet calibration | Person nuisance parameterを条件付き尤度から除去 | 補完的conditional estimator。公開selectorには未掲載 |

joint／marginal／conditionalの間で尤度、deviance、AIC、BICを比較しない。FACETSの
二桁表示値からfitを再計算せず、境界判定に必要な精度を表示値から創作しない。

## Outcome-blind Design audit

監査単位は、尤度に含まれた行から作る固有Person–Raterペアである。同じペアに複数の
Task／Criterion行があっても1割付と数える。算出するのは以下の記述量である。

- 固有割付数、割付セル密度
- Person当たりRater数、複数Raterに評定されたPerson割合
- Rater当たりPerson数、曝露CV、曝露Gini
- Raterペアの共有Person数、Jaccard overlap
- 直接Rater-overlap graphの成分数と孤立Rater

Score、Estimate、StdResidual等を変更しても結果が変わらない不変性テストを置いた。
監査は一般MFRMの全体接続性を置換しない。Task等を介して全体が接続していても、ここで
報告する直接Rater overlapは非接続になり得る。

## Fail-closed開始条件

記述監査を作れることと、感度分析を開始できることを分離した。V1開始条件は次の通り。

- Rater roleが設定または名称キーワードで確認されていること
- 元の適合が `InferenceReady` で、JMLE/MMLかつ加法的RSM/PCMであること
- regularization、latent regressionを用いていないこと
- Person座標とRater重症度座標が有限かつ報告可能であること
- frequency weightで圧縮されていない尤度行であること
- 各Person–Rater blockが同一のTask／Criterion等context signatureを持つこと
- 正確なswitch探索の計算上限内（500 assignment edges、100,000 rows）であること
- 直接Rater-overlap graphが1成分であること
- 接続性を保つ次数保存2-switchが少なくとも1つ存在すること

V1は全条件を満たす場合だけ選択する。common context signatureで停止した場合は、同じ
estimator/readiness/role/weight/connectedness契約を保ったまま、一般化MILPの変数・制約・非零数
上限とoverlap witness全域木を別に監査する。一般化側も不適格なら実行ボタンを出さない。
完全割付のように目的を改善できない設計はsolver後の不変条件で停止する。最初のfacetを
便宜的候補にしただけの場合、非接続、境界Person、非収束でも停止する。eligibleでも
実行済み・MNAR同定済みとは表示しない。

## 実装したV1 parametric runner

`mfrm_app/assignment_sensitivity.py` は観測Scoreを入力結果へ持ち出さず、次を行う。

1. Person次数、RaterのPerson曝露、Raterの応答行曝露、context、直接overlap接続性を
   不変条件として監査する。
2. 適合Person順位とRater重症度順位の整列を強める、決定論的greedy接続2-switch列を作る。
   これは局所経路であり、大域的worst caseの達成を主張しない。
3. 観測割付、要求した中間dose、整列endpointの各snapshotについて、元の適合パラメータから
   既知の適合済みPersonに
   条件づけたcategory確率を得る。MMLでは報告済みposterior/EAP座標を用い、新規Personを
   母集団から生成する研究とは区別する。
4. 反復内では `_SourceRow` ごとの同じ一様乱数を全doseで共有し、新しい応答を生成する。
5. 現在のmodel、method、facet signs、anchor、MML SD設定を固定して各シナリオを再適合する。
6. 各doseと観測の両側が `InferenceReady` かつmetric取得可能なペアだけで、Rater reference
   centered RMSE／MAE、recovery slope、Spearman、MML population SDの「dose－観測」差を出す。

利用者が要求した0／25／50／75／100%は、switch数ではなくdirection-adjusted objectiveの
始点–終点間進捗として定義する。離散経路上で同じswitchに到達する要求doseは1回だけ適合し、
要求dose集合と実際に達成したdoseの両方を保存する。doseは割付propensityの推定値でも、
大域的worst-case scaleでもない。

尤度、deviance、AIC、BICは出力せず、methodも切り替えない。JMLEのpopulation SDのような
非該当metricは、推論可能ペア数とfinite差分数を分けて0件として残す。全scenario refit、
収束、推論可能、完成ペアの分母を別々に記録する。

差分、RMSE、slope、相関はPythonの未丸め数値から計算し、CSVにも未丸め値を渡す。画面の
表示桁やFACETSの既定二桁表示を計算入力には戻さない。特にfitや閾値判断をFACETSの二桁値
から再構成する経路は、このrunnerにも存在しない。

このrunnerは、適合済み座標を既知真値ではない参照座標として用いる局所的parametric stressである。観測割付機構を
同定せず、fresh-100の既知真値効果を補正式として移植せず、整列方向以外の全MNAR機構を
覆うものでもない。

## Unequal-context MILP endpoint

`mfrm_app/assignment_context_milp.py` は、Person–Rater blockごとのTask／Criterion構成が
異なるためV1の2-switchを適用できない場合を扱う。単純なmin-cost flowでは
Person–target-Rater重複、複数contextの同時margin、直接overlap接続性を同時には表せないため、
binary MILPを採用した。

- 各元blockを分割せず、ちょうど1つのtarget Raterへ割り当てる。
- 各PersonのRater次数、各RaterのPerson数・応答行数を厳密保存する。
- 各Rater×Task×Criterion等の全context cell行数を厳密保存する。
- 元のdirect-overlap graphから決定論的全域木を作り、各辺の共有Personを両Raterへ固定する。
- solver前に元割付を全制約へ代入し、可行baselineであることを監査する。
- HiGHSが大域的optimalを返し、MIP gap、constraint residual、integrality residual、全不変条件、
  edge変更、目的改善が通った場合だけendpointを返す。

「大域的」は、このbinary変数、margin equality、witness lockの宣言済み可行領域の内部だけを
意味する。接続制約を外したworst case、実際の割付propensity、因果効果を意味しない。また、
整数解間を補間しても可行な割付にならないため、現段階では観測0とendpoint 1だけを公開し、
中間doseを作らない。score生成・common uniform・同一method再適合・未丸めmetric・
cross-basis比較禁止はV1と同じ契約を再利用する。

## MML Person generatorの分離

同じ「MML simulation」という名称で異なる問いを混ぜないため、応答generatorを二つに分けた。

1. `source_fitted`（既定）: sourceのposterior/EAP Person measureを固定する。これは局所的な
   fitted-Person条件付き感度であり、EAPを既知真値とは扱わない。
2. `mml_population_rank_preserving`: source MMLの固定または推定population SDから正規乱数を
   Person数だけ生成し、そのorder statisticsを厳密なsource EAP順位へ割り当てる。同一replicateの
   全assignment scenarioで同じ生成Person座標を使い、応答用common uniformとは独立した
   SeedSequence streamを使う。

後者はpopulationの距離・裾・実現値をreplicateごとに変えながら、counterfactualを構築した
順位–割付関係を保持する。そのため、無条件の新Person標本ではなく、source EAP順位に条件づけた
rank-preserving fitted-assignment stressである。JMLE、EAP tie、非正population SD、latent regression
ではfail closedとする。生成Person座標はreplicate内では既知だが、Rater referenceはsource fitの
ままで既知真値ではない。個人別drawはprivate exportだけに残し、公開bundleにはgenerator契約と
replicate集約だけを残す。

## Fresh-100証拠の製品内での役割

保持済みassessmentから、primary 1件とHolm-confirmed secondary 4件をfull precisionで
製品証拠registerへ写した。primaryのfree-SD MML Rater RMSE差は
`+0.09863414655507183`、95% CI
`[+0.0843219376607637, +0.11294635544937996]` である。

このregisterの用途は「なぜ同一密度感度分析が必要か」を説明することだけである。
現在のデータへの効果量移植、補正、prior設定、推定法順位づけには使わない。テストは
registerの5推定値を
`informative_assignment_confirmatory100_assessment_20260811.json` と直接照合する。

## 製品面

- GuidedのData quality: compactなestimand／割付境界
- AdvancedのData panel: 完全な契約、監査、曝露、overlap、感度計画、外部証拠
- Report tables: estimation summaryより前にestimand契約
- one-click ZIP／Excel: estimand／監査／外部証拠に加え、runner gateとblock profileを追加
  - `estimator_estimand_contract`
  - `assignment_design_summary`
  - `assignment_rater_exposure`
  - `assignment_rater_overlap`
  - `fixed_density_sensitivity_plan`
  - `informative_assignment_validation_evidence`
  - `fixed_density_assignment_runner_gates`
  - `fixed_density_assignment_block_profiles`
  - `assignment_context_milp_solver_audit`
  - `assignment_context_margin_audit`
  - `assignment_context_connectivity_witnesses`（Person IDを含むため公開bundleでは除外）
  - `assignment_generator_contract`
  - `assignment_generator_gates`
  - `assignment_generation_summary`
  - `assignment_generation_draws`（Person IDを含むため公開bundleでは除外）
- 実行後の完全・非公開Downloads bundle: dose table、dose×replicate contrasts、dose curve summary、
  endpoint互換表、completion denominator、Rater recovery、switch trajectory／ledger、全doseの
  invariants、dose別assignment map、MILP solver／context margin／witness監査
- inline ZIP: 公開用privacy filterを適用し、Person単位map／block profile／switch ledger／witnessを除外
- config export: outcome-blind scope、mechanism未同定、要求dose、dose定義、実行状態、
  seed／replicate／maxit、current-method lock、cross-basis禁止を明記

## 検証

- assignment runner focused（2-switch、MILP core、source固定／rank-preserving generator、multi-dose JMLE、MILP endpoint JMLE／固定SD・自由SD MML、negative control、失敗分母、privacy、UI、i18n）: 61 passed
- 周辺MML fixed/free SD・prior sensitivity・custom simulation回帰: 54 passed
- assignment／privacy／app smoke／end-to-endを含む広域選択回帰: 76 passed、12件は既知の相対AppTest path解決だけで失敗
- Streamlit runner harness: 実行ボタンから2 paired replicates、downloadまでexception 0
- Streamlit main absolute-path initial render: exception 0

- 新規・i18n統合: 20 passed
- estimand／FACETS丸め境界／情報的割付を含む選択回帰: 48 passed
- より広いGuided・export・qualification回帰: 120 passed
- FACETS／JMLE／MML／CMLE証拠再構築を含む選択回帰: 77 passed
- `py_compile`: pass
- locale JSON parse／key parity: pass
- Streamlit AppTest絶対パス初期描画: exception 0
- Advanced既定分析: 推定・結果再描画はexception 0。受入スクリプトの最後に
  SafeSessionStateを通常dictとして `.get()` した検査コードだけが失敗した

相対パスで `AppTest.from_file("streamlit_app.py")` を呼ぶ既存11ケースは、現在の
Streamlitが呼出元の `tests/` に対して解決するため `tests/streamlit_app.py` を探して
失敗した。絶対パスでは初期描画が通るため、本変更の製品例外とは区別した。

## 次工程

次は、(a) context-margin MILPで数学的に意味のあるnested dose集合を事前定義・検証できるかの
研究（未成立ならendpointのまま維持）、(b) 現在のrank-preserving MML generatorに加えて、
割付機構自体を生成する事前登録済みsuperpopulation研究lane、(c) anchor群別の感度解釈、
(d) 20,000行超のoffline shard runner、(e) MILPとgeneratorの多様な可行・不可行設計に対する
frozen operating-characteristic試験である。
FACETS比較は引き続き同一estimandのJMLE較正に限定し、このrunnerのMML差をFACETS不一致とは
解釈しない。
