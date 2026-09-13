# FACETS 4.5.0 対照 PCM/JMLE 資格確認

実施日: 2026-08-11  
対象: Streamlit MFRM application / `C:\Facets\Facets.exe` 4.5.0

## 結論

本アプリの次の範囲を、FACETS 4.5.0 と外部照合済みとする。

- 4カテゴリ（0, 1, 2, 3）の PCM/JMLE
- `Criterion` を唯一の step facet とする
- Criterion の各水準が、和を 0 とする独立の3閾値ベクトルを持つ
- FACETS 側は `Models=?,?,?,#,R3`
- Person は正方向の能力、Rater・Task・Criterion は正方向の厳しさ／困難度

これは FACETS の PCM 全般の互換性宣言ではない。名前付き General/Specific scale、カテゴリ数混在、再コード、閾値アンカー、欠測・疎デザイン、複数の `#` facet、raw fit parity は未検証である。

## 検証を段階化した理由

最初から多数反復を実行すると、構文、尺度表の対応付け、符号、識別制約の誤りを「安定したシミュレーション結果」と誤認しうる。このため、結果を見る前に合否基準を固定し、次の順で進めた。

1. RSM生成データをPCMで当てる構文・parser負の対照
2. 閾値共有条件と異質条件による4実行の既知真値スモーク
3. 親4実行を完全保存した20反復×2条件のパイプラインパイロット

事前計画は `facets_pcm_qualification_plan_20260811.json`、重みと方向判定の補足は `facets_pcm_qualification_amendment_20260811.json`、20反復ゲートは `facets_pcm_pipeline_pilot_plan_20260811.json` に保存した。

## 実機結果

### 構文・parserプローブ

- FACETS 4.5.0 と Python はともに収束
- Criterion 2水準を FACETS Table 8.1 / 8.2 に一対一で対応
- カテゴリ度数 8/8 が入力データと完全一致
- 主効果9個: FACETS対Python MAE 0.004396、最大差 0.010443 logits
- 閾値6個: MAE 0.009799、最大差 0.015142 logits

この段階は RSM 生成データを PCM で当てた負の対照であり、PCM真値回収の証拠には数えない。

### 4実行の既知真値スモーク

- 4/4 が両エンジンで収束・推論準備完了
- 展開した PCM 隣接カテゴリ設計の nullity は 4/4 で 0
- Table 8 度数 32/32 が完全一致
- 主効果36個: weighted MAE 0.001325、最大差 0.003059 logits
- 閾値24個: weighted MAE 0.003247、最大差 0.005929 logits
- FACETS対Python の facet内 Spearman 最小値は実質1.0
- 異質条件の基準間閾値距離が共有条件より大きい方向判定は 4/4

### 20反復パイプラインパイロット

- 試行40/40、比較適格40/40、実行別の直接一致ゲート40/40
- 親スモークの既存4実行を4/4完全保存
- 展開設計 nullity 0: 40/40
- FACETS Table 8 の Criterion×Category 度数: 320/320完全一致
- 主効果360個: weighted MAE **0.001230**、最大差 **0.003059** logits
- 閾値240個: weighted MAE **0.003128**、最大差 **0.006815** logits
- 実行別主効果 MAE の中央値 0.001251、95%点 0.001525、最大 0.001655
- 実行別閾値 MAE の中央値 0.003088、95%点 0.003784、最大 0.004157
- facet内 Spearman の全実行最小値は実質1.0
- 異質性の方向判定40/40

既知真値に対する記述的 RMSE は、Rater 0.0590–0.0615、Task 0.0436–0.0454、Criterion 0.0349–0.0383、閾値 0.0869–0.0998 logits だった。FACETS と Python の値が近いことは実装一致を補強するが、20反復は bias、coverage、false-positive、power、推定法ランキングの根拠にはしない。

### 欠測・誤指定境界パイロット

続く12生成run・24 fitの境界パイロットでは、完全データに加えて、観測数を半減してもPerson–Rater連結性を保つリング型計画欠測と、同じ観測数でPerson–Raterを2成分へ分断する負の対照を比較した。

- rank-fullのPCM/RSM fit 16/16が比較適格かつ直接一致ゲート合格
- Table 8度数144/144完全一致
- 主効果144個のweighted MAE 0.002036、最大差0.011621 logits
- 閾値72個のweighted MAE 0.003472、最大差0.009261 logits
- 異質PCM生成時のPCM対RSM log-likelihood改善方向4/4合格
- 非連結fit 8/8で独立設計とアプリがnullity 1を検出し、`InferenceReady=false`
- 同じ非連結fitでFACETSは8/8収束し、8/8で `Subset connection O.K.` と表示

最後の結果から、FACETSの収束・subset文言と、アプリの全制約展開後の予測子rankは同じ保証ではないと判断した。詳細は `FACETS_PCM_BOUNDARY_PILOT_20260811.md` に分離した。

### 同一行 estimand bridge

境界パイロットのrank-fullな8生成runを一切再生成せず、FACETS/Python JMLE、Python MML（Person SD固定・自由）、native exact CMLEへ渡した。新規fitはMML 32/32、CMLE 16/16が収束・推論適格で、同一尤度基底内のRSM対PCM改善方向は16/16合格した。

ただし、これは推定法ランキングではない。JMLEはPerson固定効果のjoint likelihood、MMLは正規Person母集団を仮定するmarginal likelihood、CMLEはPerson総得点で条件化するconditional likelihoodである。したがって、尤度・AICは各基底内でのみ比較し、方法間のパラメータ移動はdifferent-estimand sensitivityとして保持した。2反復のRMSEは次段階の仮説生成に限る。詳細は `ESTIMAND_BRIDGE_PILOT_20260811.md` に分離した。

### Person分布形状 preflight

次に、Person平均0・実現母SD 0.8を一致させ、正規・右歪み・対称二峰混合・t(3)重尾の形状だけを動かした。正しい異質閾値PCM、完全／連結計画欠測の8データセットで、FACETS/Python JMLE pair 8/8、固定／自由SD MML 16/16、exact CMLE 8/8が適格だった。FACETS対Pythonの主効果最大差は0.010392 logits、閾値最大差は0.006188 logitsだった。

32 attemptのresume再実行は全件をhash一致でskipし、再fit 0・marker変更0だった。Windows/Dropboxのdirectory finalizationと260文字path境界、初回aggregateの分布metadata欠落を検出して個別amendmentへ残した。解釈対象は `estimand_distribution_preflight_v3_20260811/aggregate_corrected/` のみである。自由SD 0.687–0.811や極端Person 0–1は1反復の候補仮説であり、biasや分布頑健性の証拠ではない。詳細は `ESTIMAND_DISTRIBUTION_PREFLIGHT_20260811.md` に分離した。

### Person分布形状 20反復 screening

登録済みscreeningは160データセット・640/640 attemptを完了した。FACETS/Python JMLE 160/160、MML 320/320、exact CMLE 160/160が適格で、640 markerのresume replayは再fit 0・変更0だった。FACETS対Pythonの主効果最大差は0.013211 logits、閾値最大差は0.009030 logits、facet内Spearman最小は1.0だった。

平均・実現SDを一致させた非正規対正規contrastでは、全体RMSEのscreening信号は主効果0/90・閾値0/30だった。一方、計画欠測対完全ではRater 20/20、Task 20/20、閾値20/20でRMSEが増えた。自由SD MMLは右歪み・重尾で正規より低く、対称二峰混合のC02局所閾値には全5モードで同方向のbias再配分が見られた。CMLE極端Personは12,800 Person-run中20件で、完全4・計画欠測16だった。

この結果は「分布形状は局所推定、観測設計は全体精度へ作用する」という次段階の仮説である。20反復から選んだcontrastを同じ20反復込みで確証してはならない。詳細は `ESTIMAND_DISTRIBUTION_SCREENING20_20260811.md` に分離した。

### Person分布形状 fresh 100反復 confirmatory

replicate 21–120だけを用いる固定N=100のconfirmatoryを実行した。800データセット・3,200 markerを完了し、resumeは3,200/3,200をhash一致でskipして再fit 0だった。MMLは1,600/1,600、exact CMLEは800/800で適格だった。

Holm補正した6仮説のうち、planned-connected設計における自由SD MMLの right_skew − normal は `-0.08983`、heavy_tail_t3 − normal は `-0.06502` で、両方向をfresh dataで確証した。後者は95%半幅 `0.02030` が登録目標 `0.02000` を僅かに超えたため、方向確証と精度missを分離して報告する。screeningでJMLE・自由SD MML・exact CMLEに共通して見えたsymmetric-mixture C02 step 2の負方向は、fresh 100では全3法が小さな正方向となり再現しなかった。

FACETS/Python JMLEは適格799/799 pairでdirect agreementを通過し、main最大差0.013272 logits、閾値最大差0.015546 logits、facet内Spearman最小1.0だった。ただしFACETSが1回だけexit 0のままU6 reportを生成せず、3,199/3,200 success、799/800 pairとなった。したがって厳格なconfirmatory batch/workbench資格は不合格で、100 pairを要求したJMLE閾値RMSE設計効果も99 pairのためinconclusiveである。

元markerを変えない別planのserial診断は3/3成功し、一過性report I/O障害という分類を支持したが、診断runをconfirmatoryへ補充していない。R 4.5.1によるpaired t、CI、Holm、論理gateの再計算は最大差 `1.78e-15` でPythonと一致した。詳細は `ESTIMAND_DISTRIBUTION_CONFIRMATORY100_20260811.md` に分離した。

## FACETS の表示精度契約

FACETS の通常表示は fit 等を小数点以下2桁程度に丸める。さらに Table 8 は固定幅のため、`Umean=0,1,6` で category Outfit が `1.` や `.` に切れる場合がある。

- measure と threshold の直接比較は `Umean=0,1,6` の主パスを使用した
- category Outfit が必要な場合だけ、同一データを `Umean=0,1,2` で再報告する補助パスを使用した
- 主・補助パスの Scorefile は分離し、高精度主出力の上書きを防止した
- 表示値から未報告の raw fit 精度を推測しない
- 表示値による閾値判定は丸め区間を持つ三値判定（pass / flag / boundary uncertain）として扱う

## 大局的な位置づけ

FACETS は照合基準であり、アプリの目的は FACETS の複製ではない。長期的な強みは、同じ生成データと同じ判断契約の下で JMLE、MML（固定／自由母集団SD）、exact CMLE などを比較し、設計・分布仮定・識別性・欠測・アンカー誤指定が推定バイアスと意思決定をどう変えるかを観察できる点にある。

計画欠測・非連結対照、RSM誤指定、同一行estimand bridge、Person分布形状screening、fresh 100 confirmatoryまでを分離して確認した。右歪み／重尾に対する自由SD MMLの分布仮定感度は確証された一方、局所mixture閾値patternは再現しなかった。次の高価値段階は、失敗runを後から補充することではなく、外部engine運用とPython evidence persistenceを分離したversioned runnerで観測設計効果を独立再検証することである。

1. Python JMLE artifactをFACETS report成功と独立に永続化する
2. exit-0/missing-reportに限る回数制限retryまたはFACETS直列化を事前登録する
3. replicate 121以降のfresh dataでH6-only観測設計効果を必要に応じて再検証する
4. UIでは推定法ランキングではなく、estimand・分布仮定・観測設計・運用資格を別表示する

## 主要成果物

- `facets_pcm_syntax_probe_20260811/`
- `facets_pcm_known_truth_smoke_20260811/`
- `facets_pcm_pipeline_pilot20_20260811/`
- `facets_pcm_boundary_pilot_20260811/`
- `estimand_bridge_pilot_20260811/`
- `estimand_distribution_preflight_v3_20260811/aggregate_corrected/`
- `estimand_distribution_screening20_20260811/screening_analysis_v2/`
- `estimand_distribution_confirmatory100_20260811/confirmatory_analysis/`
- `ESTIMAND_DISTRIBUTION_CONFIRMATORY100_20260811.md`
- `estimand_distribution_confirmatory100_assessment_20260811.json`
- `facets_pcm_probe.py`
- `facets_pcm_known_truth_smoke.py`
- `facets_pcm_boundary_pilot.py`
- `estimand_bridge_pilot.py`
- `estimand_distribution_study.py`
- `estimand_distribution_confirmatory_analysis.py`
- `estimand_distribution_confirmatory_verify.R`
- `facets_compatibility_matrix.json`
- `tests/test_facets_pcm_probe.py`
- `tests/test_facets_pcm_known_truth_smoke.py`
- `tests/test_facets_pcm_boundary_pilot.py`
- `tests/test_estimand_bridge_pilot.py`
- `tests/test_estimand_distribution_study.py`
- `tests/test_estimand_distribution_confirmatory_analysis.py`

FACETS モデル構文の参照: <https://www.winsteps.com/facetman64/models.htm>  
FACETS Rasch/PCM の参照: <https://www.winsteps.com/facetman64/raschmodels.htm>  
FACETS Rating scale の参照: <https://winsteps.com/facetman/ratingscale.htm>  
FACETS Table 8 の参照: <https://www.winsteps.com/facetman64/table8_1ratingscale.htm>
