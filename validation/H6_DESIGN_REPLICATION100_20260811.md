# H6 観測設計感度の独立100反復再現（2026-08-11）

## 結論

元の fresh-100 確認試験で99ペアのため inconclusive とした H6 を、既存データを補完・置換・併合せず、replicate 121–220 の新規固定100反復で独立再現した。登録した主要量は、normal Person・異質閾値PCMにおける、未丸め Python JMLE の6閾値RMSEについての `planned_connected - complete` である。

主要差は `+0.075464` logits、95% CI `[+0.065075, +0.085853]`、片側 `p=2.36e-26` だった。100/100ペアが有限で、95%区間半幅 `0.010389` は登録目標 `0.015` 以下である。したがって方向の独立再現と精度要件はいずれも合格した。100反復中93反復が正、7反復が負であり、これは各標本で必ず悪化するという主張ではない。

## 証拠の分離

- screening: replicate 1–20
- original confirmatory: replicate 21–120（H6は99ペアのため、現在も inconclusive）
- independent replication: replicate 121–220
- 過去120反復を主要検定に併合しない
- 元の欠測行を診断再実行で置換しない
- fixed N=100、optional stoppingなし
- JMLE・MML・CMLEの順位付け、異なる尤度基盤間のAIC比較を作らない

今回の合格は元の H6 判定を遡及変更しない。`H6R` という別の登録済み再現結果である。

## FACETS 4.5.0 副次校ブレーション

耐障害ペア層では Python JMLE を先に独立保存し、その後に FACETS を実行した。200データセットすべてで Python statistical evidence、FACETS calibration、pair full qualification が合格した。FACETS は200 try、再試行0、校ブレーション失敗0だった。

| 指標 | 200 run 中の最悪値 |
|---|---:|
| main weighted MAE | 0.004504 |
| main max absolute difference | 0.012475 logits |
| threshold weighted MAE | 0.005920 |
| threshold max absolute difference | 0.010601 logits |
| within-facet Spearman minimum | 1.0 |

FACETS 4.5.0 の表示精度契約は維持した。measure/threshold は `Umean=0,1,6`、Table 8 count/fit は別の `Umean=0,1,2` pass を使う。小数第2位に丸められた表示fitから未報告のraw fitを再構成していない。主要統計はFACETS表示値ではなく、校正済みの未丸め Python JMLE 閾値から計算した。

## 実行完全性

- preexecution: 200 datasets / 200 attempts、replicate 121–220、normalのみ、両設計、全件nullity 0
- 完全デザイン: 1,920 ratings/run
- planned-connected: 960 ratings/run
- 8 shard wall time: 673秒（先行1件スモークを除く）
- completion markers: 200/200
- Python evidence ready: 200/200
- FACETS calibration ready: 200/200
- operational failure: 0
- completion-marker set SHA-256: `9a038bbd808aa1d334cff30b6693df3fab45c370b066ee85c1e38fc8ca0cd6e8`
- resume: 200 valid skip、refit 0
- R 4.5.1 再計算: 数値・論理一致、Pythonとの差の最大 `1.78e-15`

## 事前監査が捕捉した運用境界

推定前後の都合のよい修正を避けるため、準備・実行の各版を保持した。

1. v1 は正しい行数を異なる `RunId` 順で比較して監査が停止した。fitは0件。
2. v2 のcomplete 1件ではPython証拠を保存したが、FACETS作業基底が254文字となり、3回ともexit 0のままU6 reportを作らなかった。plannedとのペア差・主要統計は観測していない。
3. v3 は短い作業領域を導入したが、監査のパス長表示が相対パス基準だったため、fit前に停止した。
4. v4 は絶対パスで最悪想定220文字を確認し、200/200を完了した。

これは単なる実装上の逸話ではない。FACETSとの比較ワークベンチでは、統計モデルの一致と、古いデスクトップ実行系のI/O信頼性を別々の品質軸にする必要がある。Python evidence と FACETS calibration を分離したことにより、外部report障害を統計的欠測へ変換せず、同時にFACETS batch gateを甘くしない構造になった。

## 大局的な解釈

独立再現されたのは「この登録DGMで、観測数を半減させたring-connected設計がJMLEの閾値回収RMSEを平均的に増加させる」という感度である。欠測一般の因果効果、任意の疎設計、MAR/MNAR、実データでのバイアスを証明したわけではない。

長期的には、FACETSとの数値一致だけを製品価値にせず、次の二層を維持すべきである。

1. 同一estimandの外部校正: FACETS対Python JMLEを、モデル・識別・表示精度を固定して検証する。
2. 補完ワークベンチ: JMLE/MML/CMLEを順位付けせず、仮定・estimand・観測設計を変えたときに結論がどう動くかを同じ入力で可視化する。

今回の結果は第1層のFACETS batch信頼性を強化し、第2層では「観測設計感度を主要な軸として明示する」根拠を与える。ただし、MMLやCMLEがJMLEより優れているという証拠には使わない。

## 主要成果物

- plan: `h6_design_replication_plan_20260811.json`
- execution registration: `h6_design_replication_execution_registration_20260811.json`
- runner: `h6_design_replication.py`
- resilient pair component: `facets_resilient_pair.py`
- qualified study: `h6_design_replication100_v4_20260811/`
- primary results: `h6_design_replication100_v4_20260811/aggregate/primary_results.csv`
- paired contrasts: `h6_design_replication100_v4_20260811/aggregate/paired_threshold_rmse.csv`
- FACETS/Python aggregate: `h6_design_replication100_v4_20260811/aggregate/`
- R verifier: `h6_design_replication_verify.R`
- assessment: `h6_design_replication100_assessment_20260811.json`

## Claim limit

本結果は、登録したnormal-Person・異質閾値PCM DGMの観測設計感度1件を独立再現し、この200-pair FACETS 4.5.0校ブレーションbatchを合格とする。FACETS全般の互換性、raw fit parity、任意の疎設計、一般的な欠測バイアス、推定法ランキングは主張しない。
