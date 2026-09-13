# Person分布形状 estimand preflight

実施日: 2026-08-11  
実機: FACETS 4.5.0 (`C:\Facets\Facets.exe`)、Python JMLE/MML/exact CMLE

## 結論

Person平均0・母標準偏差0.8を各replicateで一致させたまま、分布形状だけを正規・右歪み・対称二峰混合・t(3)重尾へ変える検証基盤が成立した。完全データとPerson–Rater連結性を保つ計画欠測を組み合わせた8データセットに対し、32/32 attemptが成功した。

- FACETS/Python JMLE parity: 8/8比較適格・直接一致ゲート合格
- Python MML: 固定SD 8/8、自由SD 8/8推論適格
- native exact CMLE: 8/8条件付きrank full・有限解・推論適格
- 制約: MML/CMLE 24/24合格
- resume replay: 32/32 skip、再fit 0、completion marker変更0
- 訂正aggregateのmetadata完全性ゲート合格

これは運用preflightであり、推定法の順位、bias、分布頑健性を示す研究ではない。

## 分離した研究軸

応答モデルは異質なcriterion別閾値を持つ正しいPCMだけに固定した。RSM/PCM誤指定を同時に動かさないため、今回観察される差は主としてPerson分布形状、観測設計、estimandの組合せに限定できる。

各実現Personベクトルは平均0・母SD 0.8へ正規化した。したがって固定SD MMLは4分布で同じ一次・二次モーメントを使い、非正規条件では形状だけが正規母集団仮定に反する。CMLEも「分布に無関係」と短絡せず、極端総得点による情報喪失を別に保持した。

## FACETS照合

8契約の最悪値は次のとおりだった。

| 指標 | 最悪値 |
|---|---:|
| 主効果 weighted MAE | 0.003114 logits |
| 主効果最大差 | 0.010392 logits |
| PCM閾値 weighted MAE | 0.003634 logits |
| PCM閾値最大差 | 0.006188 logits |
| facet内Spearman最小 | 実質1.0 |

FACETSは高精度主パスとTable 8用補助パスを分離した。表示値から未出力のraw精度を復元せず、通常の丸められたfit値は丸め区間を持つ判断材料として扱う契約を維持した。

## 記述的な候補仮説

自由SD MMLの推定値は次のとおりだった。

| Person分布 | 完全 | 連結計画欠測 |
|---|---:|---:|
| 正規 | 0.7951 | 0.7604 |
| 右歪み | 0.7441 | 0.6984 |
| 対称二峰混合 | 0.8110 | 0.7833 |
| t(3)重尾 | 0.7536 | 0.6869 |

全条件の生成時実現SDは0.8である。計画欠測および一部の非正規形状で自由SD推定が低下したことは、20反復のpaired screeningで確認すべき仮説である。1反復なのでbiasとは呼ばない。

exact CMLEの極端Personは7条件で0、重尾×計画欠測で1だった。これはCMLEの分布非依存性と、条件付き推論に利用できる総得点情報が別問題であることを示す小さな実例だが、頻度評価には反復が必要である。

## 実行基盤から得た教訓

統計以前に二つのWindows運用境界を検出した。

1. Dropbox配下では完了ディレクトリ全体の`os.replace`が拒否されたため、世代別作業物をattempt配下へ保持し、安定表と小さなcompletion markerを最後に確定する方式へ変更した。
2. 長いRunIdを保存先へ展開するとFACETS作業パスが265文字になったため、ディレクトリ名をimmutableな5桁AttemptOrdinalへ短縮し、完全なAttemptIdはmanifestとmarker内に保持した。

どちらも推定結果を見る前にamendmentを登録し、旧ディレクトリは不採用のまま新規生成した。

最初のv3 aggregateではMML recoveryに`PersonDistribution`が伝播せず4分布を混合した。fit成果物は変更せず、immutable manifestをRunIdでmany-to-one結合する事後metadata amendmentを登録し、`aggregate_corrected/`を新規作成した。元の`aggregate/`は診断証拠として残すが、解釈には使用しない。

## 計算予算

- FACETS/Python JMLE pair 8 attempt: 合計20.1秒
- 固定SD MML 8 fit: 合計98.7秒
- 自由SD MML 8 fit: 合計109.2秒
- exact CMLE 8 fit: 合計16.7秒
- 2 shard並列の観測wall time: 約132.8秒

20反復screeningでは160データセット・640 attemptとなる。MMLが計算量を支配するため、4 shard程度、attempt単位resume、completion hash監査を前提とする。ただし20反復は仮説screeningであり、既存の研究深度契約では100適格反復未満である。

## 次段階の停止条件

20反復bundleはv3の合格を受けて別ディレクトリへimmutableに準備できる。ただし、次を守る。

1. preflightのRunId、分布、設計、fit controlを変更しない
2. 20反復ではpaired contrastのMonte Carlo SD/SEを報告し、勝者を選ばない
3. 失敗・極端Personを分母から消さず、condition別に明示する
4. 100反復へ進む前に、対象contrastと精度目標を結果盲検のamendmentで固定する

## 主要成果物

- `estimand_distribution_screening_plan_20260811.json`
- `estimand_distribution_preflight_infrastructure_amendment_20260811.json`
- `estimand_distribution_preflight_infrastructure_amendment2_20260811.json`
- `estimand_distribution_preflight_aggregation_amendment_20260811.json`
- `estimand_distribution_study.py`
- `estimand_distribution_preflight_v3_20260811/retained_input/`
- `estimand_distribution_preflight_v3_20260811/attempts/`
- `estimand_distribution_preflight_v3_20260811/aggregate_corrected/`
- `estimand_distribution_preflight_assessment_20260811.json`
- `tests/test_estimand_distribution_study.py`
