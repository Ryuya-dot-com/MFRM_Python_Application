# Person分布形状 20反復 screening

実施日: 2026-08-11  
対象: FACETS 4.5.0 JMLE、Python JMLE、固定／自由SD MML、native exact CMLE  
条件: 4分布形状 × 完全／連結計画欠測 × 20 paired replicate

## 結論

160データセット・640 attemptの20反復screeningが全件成功した。FACETS/Python JMLE 160 pair、MML 320 fit、exact CMLE 160 fitがすべて推論適格で、resume replayも640/640 skip・再fit 0・completion marker変更0だった。

大局的には、二つの作用層が分かれた。

1. **観測設計**は、ほぼ全estimandで主効果・閾値の全体RMSEを一貫して悪化させた。
2. **Person分布形状**は全体RMSEを安定して変えなかったが、自由SD MMLと特定level・PCM閾値の局所biasへ作用した。

これは推定法ランキングではない。JMLE、MML、CMLEは異なるestimandであり、尤度基底間の数値比較は行っていない。また20反復はscreeningで、リポジトリの100適格反復という研究深度未満である。

## FACETS対Python JMLE

| 指標 | 160 pairの最悪値 |
|---|---:|
| 主効果 weighted MAE | 0.004545 logits |
| 主効果最大差 | 0.013211 logits |
| PCM閾値 weighted MAE | 0.005310 logits |
| PCM閾値最大差 | 0.009030 logits |
| facet内Spearman最小 | 1.0 |

160/160が比較適格かつ直接一致ゲート合格だった。FACETS measure/thresholdは登録済み高精度主パス、Table 8は別の補助パスで取得した。表示fit値から未報告raw精度を推定しない。

## 全体RMSE: 分布形状より観測設計

非正規−正規のpaired RMSE contrastでは、主効果0/90、PCM閾値0/30のscreening区間が0を除外しなかった。平均contrastの最大絶対値は、Rater 0.00882、Task 0.00572、Criterion 0.00867、閾値0.02448 logitsだった。これは分布不変性の証明ではなく、20反復で安定した全体RMSE信号を検出しなかったという意味である。

一方、計画欠測−完全のRMSE contrastは次のとおりだった。

| 対象 | 区間が0を除外 | 平均RMSE増加範囲 |
|---|---:|---:|
| Rater | 20/20 | 0.03044–0.06642 |
| Task | 20/20 | 0.01226–0.02599 |
| Criterion | 14/20 | 0.00618–0.01951 |
| PCM閾値 | 20/20 | 0.01865–0.08472 |

Person数・Rater数・生成パラメータが同じでも、観測を半減した影響は推定法を越えて残った。長期的には「推定法を選べば設計不足を補える」というUIメッセージを避け、design sensitivityを先に提示すべきである。

## 自由SD MML

全実現Person SDは0.8に固定した。

| 分布 | 完全の平均推定SD | 計画欠測の平均推定SD |
|---|---:|---:|
| 正規 | 0.7993 | 0.8286 |
| 右歪み | 0.7373 | 0.7424 |
| 対称二峰混合 | 0.8135 | 0.8417 |
| t(3)重尾 | 0.7699 | 0.7631 |

正規との差は、右歪みで完全−0.0619・計画欠測−0.0862、重尾で完全−0.0294・計画欠測−0.0655となり、4つのscreening区間が0を除外した。二峰混合との差は両設計で0を含んだ。

計画欠測−完全では、正規だけが+0.0293（screening区間+0.0025–+0.0562）で0を除外した。他の3分布は0を含んだ。1反復preflightでは計画欠測でSDが一様に低く見えたが、20反復では再現しなかった。単一seedの方向をbiasと呼ばなかった判断が正しかった。

## level別・閾値別bias

facet内平均signed errorは、mean alignmentとsum-to-zero制約によって正負が相殺される。このため初回分析を診断用として残し、level別の20反復平均誤差へ修正した。

主効果の非正規−正規contrastでは、未調整270区間中12が0を除外した。

- FACETS JMLE、Python JMLE、exact CMLE: 0件
- 固定SD MML: 9件
- 自由SD MML: 3件

特に右歪み×計画欠測では、Criterion C01/C02が固定SD MMLで約−0.0141/+0.0141、自由SD MMLで約−0.0109/+0.0109動いた。これはMMLの正規母集団仮定に固有の局所感度候補である。ただし多数の未調整screening区間から選ばれたため、新しいデータで確認する必要がある。

PCM閾値では、対称二峰混合のCriterion C02に全5モードで同方向の局所信号があった。

| 設計・閾値 | 非正規−正規contrastの5モード範囲 |
|---|---:|
| 完全・C02 category 2 | −0.0804 – −0.0748 |
| 計画欠測・C02 category 1 | +0.0817 – +0.1002 |
| 計画欠測・C02 category 2 | −0.1238 – −0.1145 |

全体threshold RMSEが安定していても、制約された閾値ベクトル内でbiasが再配分されうる。したがってworkbenchは平均RMSEだけでなく、level/step別ledgerを第一級出力にすべきである。

## exact CMLEの情報損失

CMLEは160/160 fitが適格だったが、12,800 Person-run中20件が極端総得点だった。

- 完全: 4件
- 計画欠測: 16件
- 正規: 1件
- 右歪み: 13件
- 対称二峰混合: 0件
- t(3)重尾: 6件

全体率は0.156%と小さい。ただし欠測と形状で偏りがあり、条件付き推定がPerson分布をモデル化しないことと、観測総得点が常に同じ情報量を持つことは別である。

## 計算・運用

- 8 shard wall time: 約26.5分
- FACETS/Python JMLE pair: 合計592秒
- 固定SD MML: 合計2,619秒
- 自由SD MML: 合計2,865秒
- exact CMLE: 合計523秒

MMLが計算量を支配したが、attempt単位checkpointと8 shardにより現実的に完了した。640 completion markerはresume replayで不変だった。

## 次の研究判断

100反復へ単純延長して最良条件を再利用してはならない。今回の20反復で候補を選んだため、確証段階はscreeningのreplicateを含めず、**新しい100 replicate**で行うべきである。

候補は次の4群である。

1. 自由SD MMLの右歪み／重尾対正規contrast
2. 対称二峰混合のC02 category 1/2局所閾値contrast
3. 計画欠測によるRater・Task・閾値RMSE増加
4. 右歪み／重尾と計画欠測によるCMLE極端Person率

次回は生成前に対象contrast、必要なMonte Carlo半幅、多重性、fresh-seed manifestを固定する。それまでは自動推定法推薦をUIへ追加しない。

## 主要成果物

- `estimand_distribution_screening20_20260811/retained_input/`
- `estimand_distribution_screening20_20260811/attempts/`
- `estimand_distribution_screening20_20260811/aggregate/`
- `estimand_distribution_screening20_20260811/screening_analysis_v2/`
- `estimand_distribution_screening_analysis_amendment2_20260811.json`
- `estimand_distribution_screening20_assessment_20260811.json`
- `estimand_distribution_screening_analysis.py`
- `tests/test_estimand_distribution_screening_analysis.py`
