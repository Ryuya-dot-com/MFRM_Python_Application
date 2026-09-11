# 自由分散MMLの母集団SD・不確実性監査（2026-09-11）

## 問題と修正範囲

統合ブランチの自由分散MMLは、他の推定パラメータを固定して母集団SDの
対数尤度を数値微分し、その値を「Profile SE」として表示していた。
さらに、その値から正規近似の95%信頼区間を作成していた。
これはプロファイル尤度の曲率ではなく、「やや楽観的」と誤差の程度を
保証することもできない。保存済み結果の再表示と論文用出力にも影響する。

今回、点推定アルゴリズムは維持し、不確実性の計算・表示・出力契約を修正した。

| 対象 | 修正後 |
| --- | --- |
| 推定SD | 技術的な点推定値として保持 |
| 他パラメータ固定の曲率 | `population_sd_conditional_curvature_scale` に技術診断値として保持 |
| `population_sd_se` / `population_sd_ci` | 新規推定では `None`。通常のSE・信頼区間を生成しない |
| 日英UI・方法付録・論文用案内 | SE・信頼区間の保留と根拠を表示。古い保存結果の数値も推論用に再表示しない |
| 出力資格表 | `mml.population_sd_uncertainty` を `WITHHELD`、`PublicConclusionAllowed=False` とする |
| 数値微分 | 差分点がSDの上下限に触れる場合、非有限値や非正の情報量の場合は診断値も利用不可 |

技術出力に残る旧設定値や曲率診断値は、公開結論への利用を認めるものではない。
保存された資格フラグによって現在の保留を解除することもできない。

## 数学的根拠

SDを σ、他の自由パラメータを β とする。プロファイル対数尤度は

\[
\ell_p(\sigma)=\max_\beta\ell(\beta,\sigma)
\]

であり、σを動かすごとにβも再最適化する。
[Rのプロファイル実装](https://svn.r-project.org/R/trunk/src/library/stats/R/glm-profile.R)
でも、対象係数を固定した各点で残りの係数を再推定している。

内点の停留解で観測情報行列が正定値である等の正則条件の下では、局所的な
プロファイル情報はSchur補行列

\[
I_p=I_{\sigma\sigma}
 - I_{\sigma\beta}I_{\beta\beta}^{-1}I_{\beta\sigma}
\]

となる。旧コードの値は \(1/\sqrt{I_{\sigma\sigma}}\) に対応し、
他のパラメータとの結合項を含まない。両者の一致は一般には保証されない。
また、局所曲率からのWald区間と、尤度比を反転したプロファイル尤度区間は別物である。

解析的な反例として、\(x=\sigma-1.5\) とおき、

\[
\ell(\beta,\sigma)=-\tfrac12(\beta^2+2\rho\beta x+x^2),
\qquad \rho=0.99
\]

を考える。βを0に固定すると曲率スケールは1。βを \(-\rho x\) に再最適化すると
\(1/\sqrt{1-\rho^2}\approx7.0888\) となる。これはアプリの実データで誤差が
7倍だったという主張ではなく、誤差を「わずか」とする一般的な根拠がないことの証明である。
現在のEM解が上記の正則条件を満たすとの認定も行っていない。

## 検証と残る限界

実行コマンド：

```bash
python3 -B -m pytest -q -p no:cacheprovider \
  tests/test_mml_free_population_sd.py tests/test_mml_population_sd_help.py \
  tests/test_output_qualification.py tests/test_i18n_parity.py \
  tests/test_identified_parameterization.py tests/test_mml_stationarity.py \
  tests/test_mml_engine_v2.py tests/test_readiness_report.py \
  tests/test_final_readiness_evidence.py
python3 -B -c 'import streamlit_app as app; app._self_test_gradient_checks()'
python3 -B streamlit_app.py --doctor
```

解析的反例を情報行列の逆行列と照合し、上下限、平坦・凸・非有限な尤度、
新規推定の出力、古い結果の日英表示、結果バンドル・全ダウンロード・APA表・
方法付録・論文用案内の保護を確認する。
既存の勾配検証はJMLE/MMLのRSM・PCM・GPCM、およびMMLのM-stepを対象とする。

上記の関連118テスト、勾配検証、起動診断は通過した。最後のヘルプ文言修正後も
日英表示と翻訳整合性を再確認した。
さらに、変更前の `a6582da4a` の `mfrm_estimate` を同じ数値環境で実行し、
12行の既存勾配検証データ、積分点9、最大50反復を用いて、RSM・PCM・GPCMの
固定SD／自由SDの6条件を比較した。最適化座標、目的関数、使用SD、SD更新履歴は
変更後と完全一致した。この小規模比較は変更範囲の確認であり、回復性能や
被覆率の資格検証ではない。

旧EMの尤度変化量による停止は、全パラメータの停留性を保証しない。
SDの更新に伴う積分点の再配置も、有限個の積分点による目的関数の単調改善を
自動的には保証しない。今回、推定器の変更や正規近似区間の再導入は行っていない。

次は、既存の[自由分散MML v2検証計画](../validation/MML_FREE_SD_STATIONARITY_V2_ROADMAP_20260811.md)
に従い、同定された座標での共同最適化、停留性、積分点数への安定性、
他の推定パラメータを考慮した不確実性、独立データでの被覆率を検証する。
本修正と回帰テストの通過をG2の完了やMML全体の妥当性の証明に読み替えない。

## 共同停留性・Q31/Q61の追跡確認

既に観測済みの開発用データ（24 Person、192評定、seed 20260811）で、
RSMとPCMを比較した。生成モデルはいずれもRSMであり、PCMの回復試験ではない。
`tests/test_mml_free_sd_quadrature_adapter.py` のデータと開発用の数値基準を再利用し、
積分点数を31／61に変更した。データは今後の資格検証から除外する。

| 診断 | RSM | PCM |
| --- | ---: | ---: |
| EMの停止報告 | 成功 | 成功 |
| EM停止時の共同勾配の最大絶対値 | 0.141788 | 0.161810 |
| Q31共同再最適化後の共同勾配の最大絶対値 | 1.42e-8 | 2.26e-7 |
| Q61共同再最適化後の共同勾配の最大絶対値 | 3.44e-9 | 1.40e-7 |
| Q31推定SD | 1.300951 | 1.317754 |
| Q61推定SD | 1.306430 | 1.324306 |
| Q31／Q61の構造パラメータ最大絶対差 | 0.004064 | 0.003270 |
| Q31／Q61のlog(SD)絶対差 | 0.004203 | 0.004960 |

共同勾配は、構造パラメータとlog(SD)を座標とする合計負対数尤度の微分である。
SDは内点で、再最適化後の射影勾配は通常の勾配と一致する。
これらの値は、EMの停止と有限Qでの停留性を区別する必要を示す。
Q31／Q61は既存の開発用基準を満たしたが、その基準には構造パラメータ差0.1、
log(SD)差0.1という許容幅がある。積分が高精度・厳密だとの主張には使わない。
独立Rでの今回の再照合、より高精度な積分との照合、被覆率評価は未実施である。

追跡中に、アプリの `InferenceReady` が自由分散MMLでも単にEMの成功フラグから
作られていることを確認した。現在は、`Converged` の停止記録を保持したまま、
自由分散MMLの `InferenceReady=False` と理由コードを返す。
古い保存結果には共通の要約関数を通して同じ保留を適用し、原記録は変更しない。
First Read・結果案内・最終報告・APA表・ダウンロード・感度分析の開始判定を修正した。
構造パラメータだけの `GradientNorm` がSDの微分を含まないことも明記した。
仮定監査で自由SDを「固定SD」と説明していた分岐も修正した。

再実行用の[開発プローブ](../validation/mml_free_sd_development_probe.py)と
[実行記録](../validation/mml_free_sd_development_probe_20260911.json)に、
使用ソースのSHA256、実行環境、共同最適化・再スタート・差分勾配・情報行列診断・
Q31／Q61の相互目的関数評価を保存した。保存先は新しいファイルを指定する。

```bash
python3 -B validation/mml_free_sd_development_probe.py \
  --output /tmp/mfrm_free_sd_probe_new.json
```

このプローブには開発依存関係（pytestを含む）が必要。
出力の `qualification_eligible` と `scientific_inference_ready` は常にfalseであり、
科学的な資格を発行する登録器ではない。

今回の変更後に、関連180テストと起動診断を通過した。新規推定、旧保存結果の
再表示・出力、日英の収束画面、First Read、最終報告、自由SDを基にした感度分析の
保留、および既存JMLE・固定SDの感度分析を確認した。実行記録のソースハッシュを
照合し、プローブが既存の出力ファイルを上書きしないことも確認した。
