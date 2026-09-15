# JAX・NumPyro・ConQuest・TAMの拡張計画

更新日: 2026-09-13。状態: **内部計画を精緻化。JAX/NumPyroの実装・性能検証は未着手**。
ユーザー指定の優先順位は、既存JMLE・MMLの数理検証と高速化、その後にベイズ推定の拡張。
今回、StanとNumPyro/JAXをローカルで実行するコードのダウンロードを将来の配布方針に加えた。
公開版 `0.2.17-beta` / `530b92468` の機能追加や公開資格の変更ではない。

## 今回の判断と着手順

| 順序 | 解決する問題 | 次に残す成果・完了条件 |
|---|---|---|
| 1: 現行結果の意味を固定 | 計算終了と、信頼できる推論が混同され得る | mfrmr開発版の知見を下記の対象別に照合し、Pythonの推論保留を表・図・論文・保存出力まで追跡する |
| 2: 待ち時間と迷いを測る | 大きな処理と細かなUI操作の両方が待ち時間を生む | 統合台帳の段階2とJ0を実施。推定・診断・出力・再描画を分け、[UXの次期作業](ux_adversarial_audit.md#next-ui-refinement--2026-09-13-planned)と同じ代表操作で測る |
| 3: 外部ベイズ配布の仕様を固める | 保管中のStanコードを現在のモデルと誤認するおそれ | 下記B0の静的監査、モデル対応表、最小パッケージ仕様。数値計算を増やさず着手できる |
| 4: 検証した範囲だけ実装・公開 | コード生成の成功だけでは正しい事後分布を保証できない | J1–J4は既存推定の改善として進める。新MCMC/PPCのB1以降は既存G4と計算予算の条件を満たしてから開始し、B4で配布・読込を個別に判定する |

Hosted Streamlitは設定とファイルの受け渡しを担い、長時間のStan/NumPyroサンプリングは
研究者のローカル環境で実行する。タイムアウトや再起動、切断、計算資源の制約を避ける
設計上の選択であり、「Streamlitは一律に何秒でタイムアウトする」という仕様の断定ではない。
NumPyroなら必ず速い、GPUがあれば完走する、という約束もしない。
新しいホスト上のジョブ管理や別フロントエンドは、この配布方針に必要ない。

## 解決したい問題

反復シミュレーション、感度分析、不確実性の監査では、同じモデルを何度も計算する。
利用者が結果を待つ時間を短縮し、研究に必要な反復数を実行しやすくすることが目的である。
現時点では、どの処理が総時間を支配するかを新しい統合版で測定していない。
JAXの導入だけで高速化や統計的な妥当性が達成されたとは扱わない。

## 計算基盤と照合先の役割

| 対象 | このアプリでの役割 | 最初に実施すること |
|---|---|---|
| NumPy/SciPyの既存実装 | 現行動作の基準。既存のJMLE・MMLと出力資格を保持 | 公開版の性能改善を照合し、同じ条件で時間と数値を測定 |
| JAX | 同じ統計モデル・推定対象を保つ計算基盤の候補 | 支配的な尤度・勾配計算を限定的に移植し、float64 CPUで比較 |
| NumPyro | JAXを使う後続のオフライン・ベイズ推定基盤 | Stanと共通のモデル・事前分布・出力仕様をB0で固定し、検証後に選択式コード配布を検討 |
| ConQuest | ローカルの外部照合と研究者の既存分析手順への受け渡し | 保存済みの書き出し契約を統合し、出力の読取・パラメータ対応を検証 |
| TAM | 多相MMLのスクリプトによる独立照合 | 既存のRアダプターを再利用し、共通条件を明示した照合を追加 |
| mfrmr開発版 | 共通GH求積と異なる最適化経路を持つR側の比較先 | 0.2.4.9000の感度確認機能・TAM/ConQuest比較記録を再利用し、インストール済み旧版と区別 |

NumPyroはJAXによるJMLE・MML高速化の必須依存にはしない。
JAXは自動微分・コンパイルの計算基盤、NumPyroはNUTS/HMCなどを備えた
確率的プログラミング基盤であり、役割を分ける。
[JAXの自動微分](https://docs.jax.dev/en/latest/notebooks/autodiff_cookbook.html)、
[NumPyro公式説明](https://num.pyro.ai/en/stable/getting_started.html)

## 実機で確認した出発点

これは2026-09-11に確認したローカル環境であり、サポート対象や推奨版の宣言ではない。

- Mac: arm64、Python 3.14.3。現在のPython環境には `jax`、`jaxlib`、`numpyro` は未導入。
- R: 4.6.1、TAM: 4.3.25。パッケージの存在と版を確認した。新しいモデル推定は実行していない。
- ConQuest: `/Applications/ConQuest/ConQuest` とローカルのマニュアルが存在。
  実行ファイルは x86_64。SHA-256 は
  `61d0b87f379f1578466b789866366c5cc633d31a6c3501e872861d44ff02da48`。
  これは保存済みのConQuest 5.47.5検証記録の実行ファイルと一致する。
  今回の起動、ライセンス利用可能性、arm64ホスト上の互換実行は未確認。
- 過去のR検証記録にはR 4.5.1のものもある。現在の環境での再実行を、過去の凍結証拠の置換には使わない。

最初のJAX検証対象はCPUとする。公式のインストール表はmacOS arm64 CPUを対象に含むが、
Apple GPUは実験的な位置づけであり、Mac GPUを性能目標の前提にしない。
Python・JAX・jaxlibの互換性は導入時に確認し、実際に検証した組合せを固定する。
[JAXの対応環境](https://docs.jax.dev/en/latest/installation.html)

## JAXによる高速化の順序

| 段階 | 対象と問い | 完了条件 |
|---|---|---|
| J0: 基準測定 | 公開版の改善を反映した後、初回分析と反復分析の時間を何が支配するか | データ準備・推定・診断・出力を分けた時間、メモリ、条件、版を記録 |
| J1: 小さな計算核 | 同じ引数・同定座標で尤度と勾配を再現できるか | 小規模な無アンカー・単位重みの加法的RSMから開始し、確率・尤度・勾配・制約の照合を通す |
| J2: 推定への接続 | 実際のJMLE・固定分散MMLで総所要時間を減らせるか | 既存最適化器と停止条件を維持した比較。初回と反復を分け、推定値・SE・終了状態・丸め前判定の差を記録 |
| J3: 適用範囲の拡張 | PCM、自由分散MML、複数反復でも数値と計算資源を管理できるか | 各範囲の独立した数値監査。自由分散MMLは停留性・Q31/Q61・独立R照合のゲートを通す |
| J4: 任意の計算基盤として公開 | 利用者が既存実装と選択・再現できるか | 任意依存、設定保存、版・dtype・deviceの記録、日英表示、clean checkoutの試験 |

JMLEとMMLは別々に効果を測る。初回のコンパイル時間が支配的なら、小規模な一回の分析は
NumPy/SciPyを既定に残す。反復分析で速いことを、すべての入力で速いという主張にはしない。
改善目標と数値許容差はJ0の基準測定後、JAX側の比較結果を見る前に固定する。

最初に調べる既存処理は `streamlit_app.py` の
`mfrm_loglik_jmle_value_grad`、`mfrm_loglik_mml_value_grad`、
`_m_step_expected_ll_value_grad`、`mfrm_em_mml` とカテゴリー確率の関数群。
共通の同定座標、カテゴリー支持、ファセット符号を再利用し、画面全体のJAX化は行わない。
まず既存SciPy最適化器に計算結果を渡す小さな接続で評価し、転送コストも含めて判断する。
最適化器やEM方式自体の置換は、計算核の高速化と別の変更として評価する。

自由分散MMLでは、現在の `mfrm_app/mml_stationarity.py` が `log(sigma)` の微分に
中心差分を使っている。JAXの自動微分は候補だが、sigmaに依存する求積点も含めた
同じ有限求積目的関数を微分し、既存の中心差分と独立R実装で検査する。
自動微分への置換だけで `InferenceReady` を昇格させない。

### 数値と性能の比較規則

- JAX開始時に64ビットを明示し、実際の配列dtypeも検査する。float32での結果を既存float64の同等品として扱わない。
- 学習・コンパイル前の初回時間、コンパイル後の反復時間、CPU/device間の転送、全分析の時間を分ける。
  非同期処理の完了を待って計測する。JIT呼出しが返るまでの時間だけでは比較しない。
- 比較中はデータ、初期値、乱数、反復上限、収束基準、求積、正則化、アンカーと重みを揃える。
  同一入力を保存し、乱数基盤の変更でデータ生成条件まで変えない。
- 欠測、疎な割付、未使用カテゴリー、境界、同定不足、非有限値、極端反応を含める。
  成功例だけを性能比較の分母にしない。
- SE・Hessian・勾配と、適合度などの丸め前の判定も照合する。数値差が判定を変える場合は個別に報告する。
- `backend`、実装版、依存版、device、dtypeを計算記録に含める。backendを変えた再推定は別のAnalysisIDとして残す。
  明示的にJAXを要求した実行を、無記録でNumPyに切り替えない。

[JAXの64ビット設定](https://docs.jax.dev/en/latest/notebooks/Common_Gotchas_in_JAX.html#double-64bit-precision)、
[JAXの性能計測](https://docs.jax.dev/en/latest/benchmarking.html)

## ConQuest・TAMとの関係

2026-09-11にmfrmr 0.2.4.9000のソース・既存検証とTAM実機関数、ConQuest公式仕様を
確認した。[MML求積の比較レビュー](mml_quadrature_cross_engine_review.md)を比較条件の基準にする。
点数だけで求積を同一視しない。ConQuestの既定GHと保存比較の明示的 `quadrature`、
TAMの点の密度と裾範囲、mfrmrの高次数での微小重み、推定と得点計算の積分を区別する。
既存のQ31/Q61開発確認を、一般的な積分精度の保証には使わない。

両者はアプリの代わりに裏で自動実行する必須エンジンとはせず、
計算結果を独立に照合し、利用者の既存手順へ受け渡す役割にする。
ConQuestは公式にMML・JML・ベイズMCMCを提供しており、MML専用ソフトとは扱わない。
ただし、今回の初期連携対象は保存済み契約の範囲内のRSM/PCM JML・MMLである。
TAMでは `tam.mml.mfr` と関連する多相MMLの機能を主な照合先とする。
[ConQuest公式説明](https://www.acer.org/au/conquest)、
[TAM公式マニュアル](https://alexanderrobitzsch.r-universe.dev/TAM/doc/manual.html)

| 比較する組合せ | 比較の目的 | 揃える条件・留意点 |
|---|---|---|
| NumPy/SciPy JMLE と JAX JMLE | 実装の数値互換と速度 | 同一のjoint尤度・同定・初期値・停止条件 |
| NumPy/SciPy MML と JAX MML | 実装の数値互換と速度 | 同一のmarginal尤度、母集団モデル、求積点・重みと分散の扱い |
| Python MML と TAM/ConQuest MML | 独立実装による校正 | モデル行列、符号、閾値、母集団平均・分散、求積、重み、アンカー、尤度定数を明示 |
| Python JMLE と ConQuest JML | 同じ推定対象での校正 | 同定、極端反応、制約、収束と報告精度を確認。FACETSとの既存照合も保持 |
| 将来のNumPyro と他の推定法 | 仮定や推定対象による感度 | ベイズ事後分布、MML、JMLE、CMLEを同じ推定対象とは扱わない |

個人得点はEAP・WLE・MLEを区別し、同じ種類だけを直接比較する。
既存のConQuest書き出しはMML結果にWLE報告を要求するため、Python側のEAPとそのまま比較しない。
ConQuestが持つ他の得点出力を利用する場合も、新しい対応表とfixtureを用意する。
joint・marginal・conditional尤度やベイズの予測評価量を横断してAIC/BIC順位を作らない。
尤度が完全には揃わない場合は、明示した共通条件での反応確率・構造パラメータの比較に限定する。
[ConQuestの推定・個人得点指定](https://conquestmanual.acer.org/s4-00.html)

実装は次の順に進める。

1. `archive/working-20260911` に保存済みのConQuest書き出し・manifest検証・privacy処理を統合する。
   `validation/operating_characteristics_tam.R` と独立した
   `validation/known_assignment_mml_crossfit.R` を再利用する。
2. ローカルで生成した出力を、元のbundleとの対応を検証して読み取る。
   実行終了、ファイル整合性、収束、パラメータ対応、数値一致を別々の状態として残す。
3. 同じ入力を使う小規模RSMから照合し、PCM、分散の扱い、アンカー等を範囲ごとに拡張する。
   ソフト間の一般的な同等性を一度の一致から宣言しない。
4. 必要性が確認された後に、別のローカルCLIでの実行補助を検討する。
   `/Applications/ConQuest/ConQuest` はこのMacの候補パスであり、製品に固定しない。
   実行ファイル・版・ライセンス、タイムアウト、ログ、分離した出力先を確認する。
   Hosted Streamlitから外部推定器を起動する機能は初期範囲に含めない。

公開モードの受け渡しはデータを含まないテンプレートとし、実データを含む実行用bundleは
非公開の保存範囲を維持する。ConQuestやR/TAMの未導入はPythonの通常分析を妨げない。

## mfrmr 0.2.4開発版から取り込む判断

2026-09-13に指定ルートの `development/` を再確認した。DESCRIPTIONは **0.2.4.9000 / development**、
HEADは `6dfba8258403cb0ab2f86979e78b31519753e82f`。多数の未コミット変更があるため、
今回読んだ作業ツリーをそのコミットの内容と同一視しない。通常インストール版の使用や
Rパッケージへの編集、新しい大規模シミュレーションは実施していない。
統制文書は[内部作業順序](../../mfrmr/development/inst/validation/internal-roadmap-0.2.3.md)と
[公開方向](../../mfrmr/development/ROADMAP.md)。9月7日の求積修正記録に残る旧「次の作業」より、
9月9–12日の内部作業順序を優先して読む。

| 対象 | mfrmrで確認した根拠と限界 | Python側の具体的な作業・受入条件 |
|---|---|---|
| 固定母集団RSM/PCMの構造効果・SE | 20,000データの研究と別の30,000データの確認は、指定されたq61・単位重み・設計条件の根拠 | 既存の独立尤度・勾配・共分散照合に対象と条件を対応づける。研究数を合算して全条件の反復数やPythonのcoverage証拠にしない |
| MMLの求積感度 | `mml_quadrature_sensitivity()` は同一データの再最適化と得点変化を示す。点数や安定判定を自動決定しない | `mfrm_app/mml_quadrature_sensitivity.py` と既存R照合を再利用。固定推定点での再評価、再最適化、Personの再得点化を分ける。尤度、構造効果、確率、EAP、事後SD、終了・資格状態を記録 |
| 高次数GHの微小重み | 9月12日内部記録ではQ61/Q121/Q181で4/36/74個のゼロ、未解決。生成重みと保存校正の正重み要件の不整合 | 既存記録を維持し、共通生成器・fit/scoring両経路で原因と影響を確認。任意のepsilon置換や、点数増加だけでの解決扱いをしない。Pythonへの影響は独立に判定 |
| 母集団分散とPerson得点 | 母集団の本確認は未実行。固定校正EAPと、校正を再推定するPerson区間は別の対象 | 自由SDの共同停留性・情報行列・coverageを別々に判定。EAP事後SDを構造パラメータのSEや反復標本の信頼区間に読み替えない |
| アンカー変更・Fair Score・DRF | 固定/再推定した参照、同一データ上の差の共分散、群能力差だけのDRF nullを区別。FairZの小規模確認はcoverage承認ではない | fit → 診断 → 表/図 → APA → exportで対象・参照・不確実性・保留理由を保持。現在のPythonにない機能は照合課題として残し、一括移植しない |
| GPCM・JML・CMLE | 全GPCMと推定法固有のJMLの根拠は未解決部分がある。単純なTAM共通モデルの一致では閉じない | slope owner、尺度制約、raw/corrected、極端反応処理を固定。PythonのCMLE研究経路も別資格のまま。mfrmrを無条件の正解や実行時依存にしない |

今回のソース点検を識別するSHA-256（過去の実行証拠は上書きしない）:

- `R/mfrm_core.R`: `be2355ddc2c62c6f4a83bd6f4c77b924b61c81464bc9e257e3fb0d5ce92fd4bc`
- `R/api-quadrature-sensitivity.R`: `cb06e51f8618ecaeb070527814da900c215907cc1167b6ebfe647d6fafbf927f`
- `R/core-fixed-calibration.R`: `26b5dab0b1656a41a082cb280d1ede5e1ac1e632ee50b3579e0d5b281032e5e6`

詳細な積分式・TAMの裾範囲・ConQuestの方式・分散座標の対応は、既存の
[求積レビュー](mml_quadrature_cross_engine_review.md)を再利用する。
潜在Person効果をサンプリングするNUTSは、その積分に同じGH重みを使わないが、
このことは頻度主義MMLの求積問題の修正でも、Monte Carlo誤差の解消でもない。

## Stan出力の棚卸しとNumPyro/JAXの配布仕様

### 現在残っているものと、再利用前の修正点

以下は `530b92468` の `streamlit_app.py` と公開境界テストを読んだ結果。
**Stan生成器・Posterior Viewerは保管中のコードであり、公開UIの提供機能ではない。**
`tests/test_standalone_core_boundary.py` が公開経路からの呼出しを遮断している。
`tests/test_validation_contract.py` 等の旧互換テストや合成CSV fixtureの存在も、
現行Stanコンパイル・実サンプリング・数理資格の通過を意味しない。

| 保管コード | 再利用できる部分 | 配布前の課題 |
|---|---|---|
| `build_generic_mfrm_stan_code()` | RSM/PCMの隣接カテゴリーlogit、`log_lik`、`y_rep` | RSMの `ordered[C-1] tau` は閾値の順序制約を追加する。最初のfacetだけのsoft centeringも現行の厳密な同定座標と異なる |
| `build_generic_mfrm_stan_data_export()` | prepared rows、カテゴリー符号、facet変数・ID対応表 | 非単位重みをReviewにするだけでデータは出せる。現在のアンカー・交互作用・母集団設定等を完全に表現する契約ではない。未対応条件は生成前に理由付きで停止させる |
| `stan_reproducibility_package_assets()` / `bayesian_stan_runner_templates()` | モデル・データhash、seed、chain別CSV、診断記録、Python/R/Juliaランナー | 最初は選んだエンジンのPythonランナー一つに絞る。古いファイル名・環境変数・既定事前分布を一つの実行例で照合。Uto系scaffoldの自動同梱を外す |
| `parse_stan_run_manifest_upload()` / `_posterior_load_*()` | 読取・ハッシュ照合・chain/drawを保つ設計 | 旧Stan専用manifestと、新しいエンジン共通仕様を版付きで区別。古いArviZ InferenceData前提と現行DataTreeの互換を確認し、読めない値を正常値で埋めない |

順序付きの得点カテゴリーと、ステップパラメータの大小順は別である。
例えばRSMで `eta=0, tau=(1,-1)` のlogitsは `(0,-1,0)`、確率は約
`(0.422319,0.155362,0.422319)` となり正規化できるが、旧 `ordered` 宣言では表せない。
現行RSMに対応する配布仕様ではこの追加制約を除く。順序制約モデルを提供する場合は
別モデルとして扱い、NumPyroに旧制約を無断で移さない。
またproper priorで事後分布が定義できることと、尤度だけで座標を識別できることは別である。

### B0で固定する共通モデル仕様

最初の候補は、単位重み・無アンカー・一つの観測尺度・加法的RSM・固定Person分布
`N(0,1)`。PCMは次段階とし、GPCM、交互作用、潜在回帰、多次元、driftは含めない。
これは旧Stanの未知 `sigma_theta` とhalf-Cauchy既定事前をそのまま配布する計画ではない。
未知分散を扱う場合は別のモデル版とし、事前分布と識別条件を追加検証する。

- 現在の同定・展開処理を起点に、位置/切片、facet、ステップ、Personの対応表を作る。
  すべてのブロックを機械的に中心化して必要な位置パラメータを失わない。
  Person母集団平均0を、実現したPerson標本の平均を厳密に0にする制約と混同しない。
- StanとNumPyroで同じ自由座標の基底、符号、支持、事前分布の**同時分布**を定義する。
  既存の「最後の要素を負の和で復元する」座標に独立Normalを置くと、最後の要素の
  周辺分散は他と異なる。直交基底等を使う場合も、展開後の共分散とscaleの意味を照合する。
  変換や未知scaleがある場合は密度の正規化・Jacobianを含めて比較する。
- 行ごとの基準カテゴリーlogitを0に置き、`log_prob[k] = log_prob[k-1] + eta - tau[k-1]`
  を用いる同じ確率関数を独立に照合する。共通保存形式は0始まり、Stan境界では1始まりとし、
  元の評点ラベル・欠測行・Person/facet/step索引の対応を固定する。
- 範囲外モデルや設定を別モデルへ黙って置き換えない。疎な設計、未使用カテゴリー、
  極端反応について、データの情報と事前依存を区別する。切断・同定不足には明示した扱いを持つ。

厳密なsum-to-zero変換と事前分散の関係は
[Stan公式の制約変換](https://mc-stan.org/docs/reference-manual/transforms.html)を参照する。
既定の事前scaleは文献から一律に借りず、評点確率の事前予測と感度確認から決める。

### 最小ダウンロードと結果の受け渡し

配布先は既存 **Report & Export → Files** 内の任意の「外部でベイズ推定する」入口とする。
Stan / NumPyro (JAX) はそこで初めて選ぶ。初回分析の推定法セレクターは増やさない。
一つのZIPには選択エンジン用のモデル、Pythonランナー、データ作成手順、環境定義、
短いREADME、実行設定とmanifest仕様を入れる。多数の言語・拡張モデルを一括提示しない。
公開配布はデータなしテンプレートか明示した合成例を既定とし、実データ入りの再実行bundleは
既存の非公開出力範囲で選択する。符号化済み評定やPersonの事後drawも匿名公開可能とは扱わない。

| 出力 | 共通契約と確認事項 |
|---|---|
| モデル・設定・manifest | モデル/スキーマ版、基底と制約、尤度/事前、input/model/runner hash、元AnalysisID、独立したRunID、エンジン/依存版、dtype/device、seed、chain数、warmup/draw数、終了状態、ファイル一覧 |
| `posterior.nc` | `posterior`, `sample_stats` と明示した座標・`chain × draw`。同じ保存拡張子でもスキーマ一致を検証する。元のCmdStan CSVはローカル原記録として保持 |
| `summary.csv` / 診断 | 事後平均・中央値・SD・区間方式/水準・MCSE・rank R-hat・bulk/tail ESS。Stan/NumPyroのサンプラー統計を対応づけ、欠けた診断はUnavailableとする |
| 予測出力 | 実際に生成した場合だけ事前/事後予測とpointwise log-likelihoodを保存。対象行、観測Personへの予測、新Personへの予測を区別し、容量とデータ範囲を示す |

NumPyroはfloat64 CPUから始め、配列生成前に64ビット設定を有効化して実測dtypeを検査する。
複数chainを保持し、`group_by_chain=True` と必要な `extra_fields` を実行時に指定する。
`potential_energy` は運動エネルギーを含む `energy` の代用品ではなく、E-BFMIには後者が必要。
木の深さとleapfrog数も単純な列名変更では同一視しない。
初回コンパイル、warmup、sampling、診断/保存、最大メモリ、失敗数、ESS/秒を比較する。
等しいseedでStanとNumPyroのdrawが一致することは求めず、独立chainのMCSEを含めて比較する。
[NumPyro MCMC](https://num.pyro.ai/en/stable/mcmc.html)、
[JAX dtype](https://docs.jax.dev/en/latest/default_dtypes.html)

現在の `arviz_base.from_numpyro()` はDataTreeを返し、log-likelihood生成の既定はFalseである。
旧Viewerの `arviz.from_netcdf()` / InferenceData前提と接続するには、実際の依存版を固定した
往復試験が必要。NetCDFを出せることだけで互換性を認めない。
[ArviZ公式変換仕様](https://python.arviz.org/projects/base/en/latest/api/generated/arviz_base.from_numpyro.html)

結果の読込は、ファイル整合、実行完了、MCMC診断、予測・感度確認、報告可能範囲を別々に示す。
hash一致は同一ファイルの確認であり、推定結果の正しさの証明ではない。
外部結果は元の頻度主義fitを上書きせず、異なるRunIDとして保持する。
要約CSVだけの読込ではR-hat・ESS・発散等を再計算したことにしない。

`log_lik` があるだけでLOO/WAICを自動提供しない。まず予測の単位と問いを決める。
評定行を一つ外す評価と、Person全体を外して新Personを予測する評価は別であり、後者では
Person効果の積分や実際のgroup refitが必要になり得る。条件付き行尤度の単純な足し算だけで
新Personへの予測評価が完成したと扱わない。JMLE/MMLのAIC/BICとの混合順位も作らない。

### 検証と配布の段階

| 段階 | 具体的な成果 | 受入条件 |
|---|---|---|
| B0: 仕様・静的監査（棚卸し済み、仕様確定が次） | 上記のモデル対応、保管Stanの問題一覧、配布/読込仕様、対応外条件 | 閾値逆転例、位置変換、展開後prior covariance、索引/欠測の照合例を定義。生成可能・実行可能・推論資格を区別 |
| B1: ローカルRSM実行 | 選択エンジン別の最小bundleと同じ合成入力 | G4と計算予算を確認後に開始。Stanのcompile、NumPyroのdtype、独立確率/対数尤度/勾配、固定点のlog joint（定数差は説明）、複数chain、異常中断と出力の往復を確認 |
| B2: ベイズ推定の検証 | 事前予測・事後予測、事前感度、独立エンジン比較、SBCの限定した研究計画と実行記録 | rank R-hat・ESS・MCSE・発散・energy等を判定。許容差、対象、失敗分母、反復数とMC精度は確認結果を見る前に固定。SBCと固定真値での頻度主義coverageを別の問いとして報告 |
| B3: PCM等の追加 | PCMのstep facet・ladder別制約の照合。未知分散はさらに別のモデル版 | RSMの成功を継承せず対象ごとにB1/B2を実施。Utoの一般化・多次元・driftモデルは原論文の式・事前・同定との対応監査後の別研究課題 |
| B4: 配布とViewer | 選択エンジンのダウンロード、対応する結果読込、日英案内・再現例 | 配布は単独で判定できるが、読込未提供なら「戻して閲覧できる」と案内しない。旧境界テストを削るだけで再公開しない。契約を満たした範囲の正負テスト、対応環境、公開/非公開出力、UI受入を確認 |

## 文献・仕様の根拠

Zotero local APIを2026-09-13に検索し、以下の書誌・DOIを照合した。
キーはZotero item keyでありBibTeX citekeyではない。添付PDFや索引全文の取得、
ライブラリへの追加・修正は行っていない。文献モデルの忠実再現を認定したレビューではなく、
B0/B3で読むべき原論文と、検証の問いの対応を固定するための参照である。

| Zotero item key | 文献 | この計画での役割 |
|---|---|---|
| `IR3LFDRM` | Uto & Ueno (2020), [A generalized many-facet Rasch model and its Bayesian estimation using Hamiltonian Monte Carlo](https://doi.org/10.1007/s41237-020-00115-7) | 一般化モデルとHMCの原仕様。加法的RSMの配布と混同しない |
| `38TX837G` | Uto (2021), [A multidimensional generalized many-facet Rasch model for rubric-based performance assessment](https://doi.org/10.1007/s41237-021-00144-w) | 多次元モデルの後続範囲 |
| `4DTXAU8H` | Uto (2023), [A Bayesian many-facet Rasch model with Markov modeling for rater severity drift](https://doi.org/10.3758/s13428-022-01997-z) | 時間を意味する設計とdriftモデル。刊行元でオンライン公開2022-10-25、巻号年2023・55巻3910–3928頁を確認。Zoteroの2022は保持し、巻号付き引用は2023に統一する計画 |
| `LAA54HWF` | Gabry et al. (2019), [Visualization in Bayesian workflow](https://doi.org/10.1111/rssa.12378) | 事前/事後予測を含む確認手順の参照 |
| `GF6R9HF7` | Vehtari et al. (2017), [Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC](https://doi.org/10.1007/s11222-016-9696-4) | 予測評価の対象と診断。重複するZotero項目を別の証拠数にしない |
| `PKQMUBH7` | Morris et al. (2019), [Using simulation studies to evaluate statistical methods](https://doi.org/10.1002/sim.8086) | 問い、生成条件、推定対象、指標、MC精度、失敗の分母を事前に決める |
| `NBA428KJ` | Wind et al. (2023), [Does sparseness matter?](https://doi.org/10.1177/01466216231182148) | 疎な評定設計を検証条件に含める際の文献候補 |

NumPyro・SBC・rank-normalized R-hatの直接対応項目は今回のZotero検索では見つからなかった。
実装APIは上記公式文書、MCMC診断は
[Vehtari et al. (2021)の著者所属機関公開版](https://aaltodoc.aalto.fi/items/332dd844-1b47-4876-a7ca-f2b54a33e4b1)、
SBCは[Talts et al., Validating Bayesian Inference Algorithms with Simulation-Based Calibration](https://arxiv.org/abs/1804.06788)
の書誌・要旨で補った。SBCの具体的プロトコルと事前分布はB2の着手前に原典の方法を精読して固定する。
UI参考投稿の読取範囲と適用案は[UX点検記録](ux_adversarial_audit.md#next-ui-refinement--2026-09-13-planned)に集約する。

## 統合計画への適用

- [開発統合台帳](development_integration.md)の段階2で公開版の改善を照合してからJ0を開始する。
  外部照合の準備は段階5に組み込む。
- JAXの数値互換検証は現行モデルの改善として進める。G0/G1を維持し、G2以降や
  [自由分散MML v2の資格](../validation/MML_FREE_SD_STATIONARITY_V2_ROADMAP_20260811.md)を省略しない。
- NumPyroの研究・公開は[統計エンジンの修正計画](statistical_engine_remediation_roadmap.html)の
  新モデル・MCMCのゲートを維持する。今回の将来計画への追加は、公開機能の検証完了を意味しない。
- 2026-09-13の今回の作業は内部文書の更新。依存関係のインストール、既存推定コードの変更、
  MCMCやConQuest/TAMの新規推定、今回の変更のコミット・push・デプロイは実施していない。
