# 既知の確率的割付 multi-vector preflight 最終状態

実施日: 2026-08-11  
現在の判定: **運用ゲート PASS・非確認的 advancement signal PASS**

## 大局的な位置づけ

4本の新しい独立Personベクトルについて、固定周辺度数を保つ既知の確率的Person–Rater割付を `gamma=-0.8, 0, +0.8` で生成した。12データセット、48推定attempt、seed、推定法、指標、方向判定はデータ生成前に固定した。以前の固定Person screeningも、この4本のpreflightも、将来の確認試験にはpoolしない。

このpreflightの役割は、FACETSとの同一estimand校正を維持しながら、JMLE、固定SD MML、自由SD MML、exact CMLEが割付ストレスにどう反応するかを点検することである。推定法の順位付けや普遍的バイアスの推定ではない。

## 元の48 attempt

- completion marker: 48/48
- independently persisted Python JMLE: 12/12 ready
- fixed-SD MML: 12/12 ready
- free-SD MML: 12/12 ready
- exact CMLE: 12/12 ready
- 再開監査: 48/48 hash検証、再fit 0

最初のFACETS呼出しは12/12で `0xC0000005` を返した。この失敗、完了マーカー、成果物は削除・置換・上書きしていない。

## FACETS障害の訂正診断

当初の「FACETS 4.5.0ランタイム全体が使用不能」という判断は広すぎた。

- 引数なしの対話GUI: 起動成功
- direct CMD `BATCH=YES`: `0xC0000005`
- Python `BATCH=YES`: `0xC0000005`
- `BATCH=NO`でもWindows側でwindowをhiddenにした場合: `0xC0000005`
- visible `BATCH=NO`: 正常完了

したがって、Pythonそのものや統計仕様ではなく、FACETS/Xojoのhidden-window初期化経路が必要原因だった。可視ランチャーは、reportと全facet Scorefileが非空かつ安定するまで待ち、FACETS所有windowへ `WM_CLOSE` を送る。強制終了は成功として扱わない。

さらに、最初の補完実行ではprimary passが完了した一方、auxiliary `Umean=0,1,2` passがTable 5後に停止した。primaryの `scores.4.txt` は257文字、auxiliaryの `scores_u2.4.txt` はちょうど260文字だった。同一spec/dataを最長98文字の短いパスで実行すると5.3秒でreportと4 Scorefileが揃った。このため、共通FACETS呼出し層に220文字の保守的な事前path gateを追加した。

## 別保存したFACETS補完校正

短い新規ディレクトリで同一12入力をゼロから実行した。

- calibration ready: 12/12
- direct agreement pass: 12/12
- retry: 0
- 主効果 weighted MAEの最大値: 0.004645 logits
- 主効果 absolute differenceの最大値: 0.012604 logits
- PCM閾値 weighted MAEの最大値: 0.005354 logits
- PCM閾値 absolute differenceの最大値: 0.009033 logits
- facet内Spearmanの最小値: 1.0

FACETSの通常fit表示は小数第2位までの表示値である。未報告のraw fitを逆算・補間していない。measure/threshold校正には登録済みの `Umean=6` primary passを用い、Table 8表示層とは分離した。

## 非破壊evidence join

元の48 markerと補完12 markerを `AttemptId + RunId + RunInputSHA256` で一対一結合した。元の科学行はそのまま保持し、補完FACETS行は `CalibrationOnly=True`、`IncludedInStudy=False` の別表に隔離した。

- 元のPython JMLEと補完Python JMLEの主効果差: 0
- 閾値差の最大値: `4.44e-16`
- join audit: 全gate PASS
- 元の失敗artifact変更: なし

## 登録済み4指標

| 指標 | n=4平均 | 範囲 | 登録方向 | 方向一致 |
|---|---:|---:|---|---:|
| PF1 free-MML symmetric-stress Rater RMSE | +0.09533 | +0.04510 ～ +0.12764 | 正 | 4/4 |
| PF2 fixed-MML symmetric-stress Rater RMSE | +0.09216 | +0.03879 ～ +0.11956 | 正 | 4/4 |
| PF3 free-MML symmetric-stress SD shift | -0.07288 | -0.11150 ～ -0.03169 | 負 | 4/4 |
| PF4 free-MML direction-aligned Rater slope | +0.36545 | +0.17550 ～ +0.47308 | 正 | 4/4 |

PF1は「平均が正、かつ少なくとも3/4が正」という登録済み前進条件を満たした。全指標とも4/4で登録方向だったが、n=4の設計安定性確認であり、p値も確認的結論も生成していない。

R 4.5.1は、集約済み要約を再利用せず、native-precision recovery/run行から16個の値と4要約を再計算した。Pythonとの差は最大 `1.67e-16` で、方向数とadvancement flagも一致した。このR検証はPython endpoint読後の実装照合であり、独立確認試験ではない。

## 次の長期ゲート

新しいPersonベクトルとseedを用いる別確認計画を、生成前にfreezeできる。ただし次を維持する。

1. この4ベクトルと過去screeningを確認試験へpoolしない。
2. FACETSは同一estimandのJMLE校正であり、MML/CMLEのgold standardではない。
3. joint・marginal・conditional likelihoodのAIC/BICを横断比較しない。
4. estimator rankingを作らない。
5. FACETSの2桁fit表示をraw値として扱わない。
6. visible-window依存と220文字path budgetを運用前提として明示する。

## 主要artifact

- [元study identity](known_assignment_multivector_preflight4_20260811/study_identity.json)
- [補完校正assessment](known_assignment_multivector_preflight4_20260811/v/assessment.json)
- [join audit](known_assignment_multivector_preflight4_20260811/visible_join_audit_20260811.json)
- [derivative aggregate assessment](known_assignment_multivector_preflight4_20260811/aggregate_visible_join_20260811/assessment.json)
- [PF1–PF4 summary](known_assignment_multivector_preflight4_20260811/aggregate_visible_join_20260811/registered_preflight_summary.csv)
- [R verification](known_assignment_multivector_preflight4_20260811/r_verification_20260811/assessment.json)
- [path/launch diagnosis](FACETS_450_NATIVE_CRASH_DIAGNOSIS_20260811.md)
