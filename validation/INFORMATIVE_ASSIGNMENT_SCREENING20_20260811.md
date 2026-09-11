# 情報的割付 20反復スクリーニング

実行日: 2026-08-11  
対象: FACETS 4.5.0 JMLE、Python JMLE、固定SD MML、自由SD MML、exact CMLE  
設計: 3観測設計 × 20反復、60データセット、240推定attempt

## 結論

全240推定が成功し、補正後集計は全ゲートに合格した。等しい疎データ密度のまま、潜在能力が低いPersonを寛大なRater、高いPersonを厳しいRaterへ割り付けると、**MMLのRater厳しさがゼロ方向へ圧縮される候補シグナル**が現れた。一方、JMLEとexact CMLEのRater RMSE、ならびに全推定法のPCM閾値RMSEには20反復で同じシグナルが見られなかった。

これは推定法ランキングではない。20反復はscreeningであり、4つの区間は未調整である。さらに、人工的に構成した割付機序であって、実データのMAR/MNAR判定や補正法の検証ではない。次の100反復を独立に事前登録して初めて確認段階へ進める。

## 何を固定し、何を変えたか

complete response matrixから、同一反復内で次の3設計を切り出した。

| 設計 | 行数 | Rater/Person | 役割 |
|---|---:|---:|---|
| complete | 1,920 | 4 | データ密度の文脈 |
| planned_connected | 960 | 2 | 潜在能力と独立な疎設計 |
| ability-severity aligned | 960 | 2 | 能力四分位とRater厳しさを整列 |

2つの疎設計では、各Raterが40 Person・240行を担当し、Person–Raterグラフは1成分、制約後PCM nullityは0である。生成済みスコア、Personベクトル、Rater露出、Task/Criterionの観測数を揃えた。変えたのは「誰がどのRaterに割り付けられるか」だけである。割付は生成後のスコアではなく潜在Thetaと真のRater厳しさを使用した。

この設計により、先のH6Rで確認した**観測密度低下**と、今回の**潜在特性に依存する割付**を分離して考えられる。

## 実行・校正・独立検証

- 60/60 FACETS–Python JMLE pairが校正可能
- FACETS try 60、retry 0
- MML 120/120、exact CMLE 60/60が推論可能
- 制約監査 180/180合格
- completion marker 240/240、resume replayは240/240 skip・再fit 0
- R 4.5.1が原始truth-error行から32集計セル、640 paired metric行、16 Rater-levelセル、20自由SDペアを再構築
- RとPythonの最大絶対差は `2.78e-17`

FACETSとのJMLE直接校正は、主効果weighted MAE最大0.004732、主効果最大差0.014063、閾値weighted MAE最大0.006434、閾値最大差0.010917 logits、facet内Spearman最小1.0だった。FACETSの通常の小数点以下2桁表示からfitやtruth errorの隠れた精度を復元していない。登録済みUmean=6 measureと別取得のTable 8を用途別に扱い、科学統計はPythonのネイティブ精度行から計算した。

## 32セルのscreening結果

`aligned − planned` のpaired contrastを推定法内で計算した。区間が0を除外したのは32セル中4セルで、すべてMMLのRater回収だった。

| 推定 | 指標 | 平均差 (logits) | 95% screening区間 |
|---|---|---:|---:|
| 自由SD MML | Rater RMSE | +0.07575 | +0.04747 ～ +0.10404 |
| 自由SD MML | Rater MAE | +0.07239 | +0.04506 ～ +0.09972 |
| 固定SD MML | Rater RMSE | +0.04109 | +0.01624 ～ +0.06595 |
| 固定SD MML | Rater MAE | +0.03706 | +0.01503 ～ +0.05908 |

RMSEとMAEは同じ現象を異なる損失で見ているため、4件を4つの独立発見とは数えない。

## 「圧縮」と呼ぶ根拠

Raterレベル別の平均truth-error差は次のとおりだった。正値はaligned設計で推定値がより正方向へ、負値はより負方向へ動いたことを表す。

| 推定 | R01 (-0.45) | R02 (-0.15) | R03 (+0.15) | R04 (+0.45) |
|---|---:|---:|---:|---:|
| exact CMLE | +0.0034 | +0.0051 | +0.0048 | -0.0132 |
| JMLE | +0.0034 | +0.0034 | +0.0070 | -0.0137 |
| 固定SD MML | +0.1705 | +0.0917 | -0.0860 | -0.1762 |
| 自由SD MML | +0.2189 | +0.1161 | -0.1079 | -0.2271 |

MMLでは負の厳しさが上へ、正の厳しさが下へ動き、両端が0へ近づいた。4レベルの真の厳しさに対する記述的傾きは、固定SD MMLで−0.406、自由SD MMLで−0.521だった。JMLEは−0.016、exact CMLEは−0.017である。この傾きは登録済みのlevel別出力を要約した診断値であり、今回の推論エンドポイントではない。

## 自由SD MML

真の実現Person SDは両設計とも0.8である。それでも自由SD MMLの平均推定値は、plannedの0.79933からalignedの0.67455へ低下した。paired差は−0.12478、95% screening区間は−0.16836～−0.08120だった。

Rater厳しさの圧縮とPerson分散の縮小が同時に生じたことは、正常Person分布とRater厳しさの分離が情報的割付の下で曖昧になる可能性と整合する。ただし、これは観察されたパターンに対する機序仮説であり、因果的に同定された説明ではない。

## シグナルが見られなかった領域

主要なRMSE差は次のとおりで、screening区間はいずれも0を含んだ。

| 推定・対象 | 平均差 | 95% screening区間 |
|---|---:|---:|
| JMLE Rater | −0.00781 | −0.04299 ～ +0.02738 |
| exact CMLE Rater | −0.00597 | −0.03558 ～ +0.02364 |
| JMLE PCM閾値 | −0.01045 | −0.05298 ～ +0.03208 |
| 固定SD MML PCM閾値 | +0.00449 | −0.02516 ～ +0.03414 |
| 自由SD MML PCM閾値 | +0.00036 | −0.03482 ～ +0.03555 |
| exact CMLE PCM閾値 | +0.00204 | −0.03251 ～ +0.03659 |

したがって、先のH6Rで確認されたJMLE閾値RMSEの悪化を、単に「割付が能力依存だから」と説明する証拠は得られなかった。H6Rはcomplete対plannedの密度差、今回は等密度のplanned対alignedであり、異なる摂動を検査している。

区間が0を含むことは、JMLE/CMLEの不変性や優越性を証明しない。20反復で候補シグナルを検出しなかったという限定的な結果である。

## メタ認知的な位置づけ

今回の最も重要な知見は、「欠測率」だけでは推定の難しさを表せないことである。同じ50%観測、同じRater露出、同じ接続性でも、観測される組合せが潜在能力とRater厳しさに依存すると、正常Person母集団を積分するMMLはRater差とPerson分散を別の形で配分しうる。

一方、JMLEはPersonを固定効果として推定し、exact CMLEはPerson総得点で条件づける。このestimand/尤度基底の違いは今回のパターンに対する妥当な説明候補だが、推定法の一般的優劣を意味しない。FACETSはJMLEの外部基準として機能し、MML/CMLEに対する「正解ラベル」ではない。

長期的には、workbenchの表示順も次のようにすべきである。

1. 観測設計・接続性・割付依存性を先に診断する。
2. estimandとPerson分布仮定を明示する。
3. 同一推定法内で設計感度を示す。
4. FACETS parityはJMLE校正として別枠に置く。
5. 異なる尤度基底のAIC・log-likelihood・順位表を作らない。

## 集計補正の監査

最初のaggregateは、MML/CMLEの閾値行に`Design`と`Replicate`が伝播せず、期待32セルに対して26セルしか形成できなかったためqualification falseだった。これはfit失敗ではない。

結果を見る前にmetadata-only amendmentを固定し、manifestから登録済みメタデータをjoinした。推定値・truth error・completion marker・fitは変更せず、元aggregateも保存した。補正後は32/32セル、各20ペアとなり全ゲートに合格した。この区別により、解析コードの欠陥を科学的失敗として数えることも、逆に黙って修正して成功扱いすることも避けた。

## 次の判断

screening反復221～240はここで閉じる。確認試験は反復241～340の新しい100反復を用い、自由SD MMLのRater RMSE増加を単一primary endpointとする。固定SD MMLのRater RMSE、自由MMLのPerson SD低下、事前定義したRater厳しさ圧縮傾きをsecondary/mechanistic endpointとして固定する。反復1～240は再利用しない。

主な証跡:

- `informative_assignment_screening_plan_20260811.json`
- `informative_assignment_screening20_20260811/preexecution_audit.json`
- `informative_assignment_screening20_20260811/aggregate_corrected/`
- `informative_assignment_screening_verify.R`
- `informative_assignment_screening20_20260811/aggregate_corrected/r_verification.csv`
- `informative_assignment_screening20_assessment_20260811.json`
