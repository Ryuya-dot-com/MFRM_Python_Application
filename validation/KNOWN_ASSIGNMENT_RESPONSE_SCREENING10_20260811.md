# 既知の確率的割付機構：response screening 10

実施日: 2026-08-11  
最終資格判定: **PASS（screeningのみ）**  
FACETS: 4.5.0 (`C:\Facets\Facets.exe`)  
独立再構築: base R 4.5.1

## 問い

Person/Rater次数、Rater exposure、Task×Criterion context、連結性、固定Person truth、潜在的な完全応答を保ったまま、既知の能力×Rater severity割付依存を変えたとき、各estimandの回復誤差はどう動くかを点検した。

割付機構は次式で、4-Rater固定周辺度数に対してexact dynamic programmingから独立抽出した。

\[
P(G\mid d_P,d_R,\text{connected})\propto
\exp\{\gamma\sum_{(p,r)\in G}z(\theta_p)z(\delta_r)\}.
\]

`gamma` は実データから推定した値ではなく、事前に資格化したsimulation doseである。

## 凍結設計

- 80 Persons、4 Raters、3 Tasks、2 Criteria、4 categories
- 各Personは2 Raters、各Raterは40 Persons／240 response rows
- 1 datasetは960 rows、constrained PCM nullityは全て0
- Person truthは平均0、母集団SD 0.8の1つの凍結vector
- `gamma=-0.8, 0, +0.8`、各10独立assignment graphs
- 各replicateで3 gammaは同じ完全潜在応答行とrow-uniform streamを共有
- 正しく特定されたheterogeneous-threshold PCM
- FACETS/Python JMLE、固定SD Q31 MML、自由SD Q31 MML、native exact CMLE

実現した平均割付相関は `-0.4329`, `-0.0047`, `+0.3980` だった。割付生成はScoreを入力にせず、全30 graphで固定次数、連結性、行・context exposureを保持した。

これはPerson-populationを反復抽出するstudyではない。Monte Carlo変動は、1つの固定Person vectorに条件付けたassignment drawとresponse drawによる。

## 実行資格

- 登録datasets: 30/30
- 登録attempts: 120/120 completion、native evidence ready
- FACETS/Python JMLE: 30/30 calibration ready
- 固定SD MML: 30/30
- 自由SD MML: 30/30
- exact CMLE: 30/30
- MML/CMLE constraints: 90/90 pass
- FACETS tries: 30、retry: 0
- 5 estimator modes×4 recovery domains×2 metrics×2 stress doses×10 pairs: 800/800

FACETS対Python JMLEの最大差は次の範囲だった。

| JMLE calibration metric | Worst case |
|---|---:|
| main weighted MAE | 0.00404 |
| main absolute difference | 0.01088 |
| threshold weighted MAE | 0.00570 |
| threshold absolute difference | 0.00933 |
| minimum within-facet Spearman | 1.00000 |
| independently repeated Python main/threshold difference | 0 / 0 |

これは同一fixed-person joint estimandの外部engine校正である。FACETSをMMLやCMLEの正解基準にはしていない。

## 主なRater recovery結果

以下は各推定法内でのRater RMSEであり、行間の順位比較を意図しない。

| Estimator lane | gamma -0.8 | gamma 0 | gamma +0.8 |
|---|---:|---:|---:|
| FACETS JMLE | 0.1106 | 0.0764 | 0.0769 |
| Python JMLE | 0.1070 | 0.0788 | 0.0801 |
| MML fixed SD 0.8 | 0.1503 | 0.0596 | 0.1027 |
| MML free SD | 0.1540 | 0.0615 | 0.1092 |
| exact CMLE | 0.0809 | 0.0662 | 0.0662 |

stress minus neutralのpaired screening contrastは次のとおりだった。

| Estimator lane | gamma | Mean Rater-RMSE contrast | Screening 95% interval |
|---|---:|---:|---:|
| FACETS JMLE | -0.8 | +0.0342 | [-0.0163, +0.0847] |
| FACETS JMLE | +0.8 | +0.0004 | [-0.0370, +0.0378] |
| Python JMLE | -0.8 | +0.0282 | [-0.0254, +0.0818] |
| Python JMLE | +0.8 | +0.0012 | [-0.0369, +0.0394] |
| MML fixed SD 0.8 | -0.8 | **+0.0907** | **[+0.0423, +0.1391]** |
| MML fixed SD 0.8 | +0.8 | **+0.0430** | **[+0.0090, +0.0771]** |
| MML free SD | -0.8 | **+0.0925** | **[+0.0412, +0.1438]** |
| MML free SD | +0.8 | **+0.0477** | **[+0.0087, +0.0867]** |
| exact CMLE | -0.8 | +0.0147 | [-0.0319, +0.0613] |
| exact CMLE | +0.8 | 約0 | [-0.0351, +0.0351] |

同じ選択性はRater MAEにも現れた。全80のRMSE/MAE contrast groupのうちscreening intervalが0を含まなかったのは、固定／自由SD MMLのRater domain 8行だけだった。Task、Criterion、PCM threshold、JMLE、CMLEのintervalは全て0を含んだ。

これはp値による確証ではない。n=10のscreening intervalを、普遍的な頑健性や推定法優越性の判定へ変換しない。

## biasの方向

MMLのRater level truth errorは、割付方向に応じた明瞭な候補patternを示した。

| MML lane / gamma | R01 truth error | R04 truth error | pattern |
|---|---:|---:|---|
| fixed SD / -0.8 | -0.1627 | +0.1791 | severity spreadの拡大 |
| fixed SD / +0.8 | +0.0756 | -0.1357 | severity spreadの圧縮 |
| free SD / -0.8 | -0.1670 | +0.1833 | severity spreadの拡大 |
| free SD / +0.8 | +0.0866 | -0.1450 | severity spreadの圧縮 |

正のgammaでは高能力者がよりsevereなRaterへ割り当てられ、Rater間の観測応答差が相殺されやすい。負のgammaでは高能力者がよりlenientなRaterへ割り当てられ、差が増幅されやすい。Normal-person MMLは共通のPerson分布とRater assignmentの独立性を暗黙に使うため、この既知の依存をRater severityへ配分し、圧縮／拡大として現した可能性が高い。

この説明は本DGMに対する機構的解釈であり、MML一般の欠陥やJMLE/CMLE一般の優越を意味しない。実データではgammaもmissingness mechanismも未同定である。

## 自由MMLのPerson SD

| gamma | Mean estimated SD | Screening 95% interval |
|---:|---:|---:|
| -0.8 | 0.7229 | [0.6863, 0.7595] |
| 0 | 0.8154 | [0.7802, 0.8507] |
| +0.8 | 0.7622 | [0.7251, 0.7994] |

真の固定vector SDは全条件で0.8である。割付依存はRater recoveryだけでなく、marginal population spreadの推定にも移動候補を生じさせた。

## FACETS表示精度

FACETS measure/thresholdは登録済み`Umean=6` primary passを用いた。通常表示で丸められるfitは別の補助証拠であり、Rater RMSE、truth error、interval、JMLE parityのraw入力に使っていない。小数第2位fitから未報告のraw residual／fitを復元していない。

## 独立R再構築

base R 4.5.1で、`recovery.csv`から5 estimator lanes×2 stress dosesのRater RMSE contrast 10行を独立に再構築した。Pythonとの差はmean、SDとも最大 `1.39e-17` でPASSした。これはRater RMSE contrastだけの数値再構築であり、confirmatory claimを追加しない。

## 保持された運用失敗と分析修正

1. 最初のstudyはmanifest互換列`Seed`欠落により30 JMLE attemptsが全てengine前に失敗した。失敗markerを保持し、ratings/truth/assignmentを変えず`Seed=UniformSeed`を追加したv2を作成した。
2. v2 wrapperの最初の起動はrepository import path不足でattempt前に停止した。`sys.path`起動修正だけをv3登録した。
3. 最初のaggregateはMML/CMLE threshold行のoptional `Replicate` metadataが空でgroupbyから落ち、構造ゲートにFAILした。数値式を変えず、RunIdから登録manifest contextを再付与した`aggregate_v2`が全ゲートにPASSした。

これらは「最終的に動いた」だけでなく、外部engine、統計推定、runner互換性、集計metadataを別の品質層として扱う必要性を示す。

## 結論と次段階

このアプリの補完的価値は、FACETS JMLEを外部校正したうえで、異なるPerson処理が既知の割付依存にどう反応するかを同じDGM上で分離できる点にある。今回のscreeningは、Normal-person MMLのRater recoveryとpopulation SDが割付依存に選択的に反応する候補を示した。

次の確証studyを行うなら、今回の固定Person vectorをそのまま反復するだけでは不十分である。少なくとも複数の事前固定Person vectors／population draws、gamma dose-response、独立assignment draws、固定サンプルサイズ、MML Rater RMSEとseverity compression slopeの事前登録が必要である。FACETS比較は引き続きJMLE同士に限定する。

主要artifact:

- 計画: [known_assignment_response_screening_plan_20260811.json](known_assignment_response_screening_plan_20260811.json)
- 最終判定: [assessment.json](known_assignment_response_screening_v2_20260811/aggregate_v2/assessment.json)
- paired contrast: [gamma_contrast_summary.csv](known_assignment_response_screening_v2_20260811/aggregate_v2/gamma_contrast_summary.csv)
- R再構築: [assessment.txt](known_assignment_response_r_verification_20260811/assessment.txt)
- Seed remediation: [known_assignment_response_seed_remediation_20260811.json](known_assignment_response_seed_remediation_20260811.json)
- aggregate metadata remediation: [known_assignment_response_analysis_metadata_amendment_20260811.json](known_assignment_response_analysis_metadata_amendment_20260811.json)
