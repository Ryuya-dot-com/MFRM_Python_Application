# 情報的割付 100反復確認試験

実行日: 2026-08-11  
fresh replicate: 241～340（screening 1～240は不使用）  
対象: FACETS 4.5.0 JMLE、Python JMLE、固定SD MML、自由SD MML、exact CMLE

## 結論

情報的割付に対する**MMLのRater厳しさ圧縮**は、独立した100反復で確認された。

単一primary endpointである自由SD MMLのRater RMSE差 `aligned − planned` は **+0.09863 logits**、95%区間 **+0.08432～+0.11295**、片側p値 `7.70e-25` だった。登録した95%区間半幅目標0.015に対し実測0.01431で、方向確認とMonte Carlo精度の両方に合格した。

secondary 4件もprimary gate通過後、Holm family-wise 0.05で全件確認された。固定SD MMLでもRater RMSEは+0.06851増え、自由SD MMLの推定Person SDは−0.11839低下した。自由/固定SD MMLの厳しさ圧縮傾きもそれぞれ−0.489、−0.384だった。

この結果は、FACETSよりMMLが「悪い」というランキングではない。異なるestimandを同じ尤度/AICで競わせる根拠にもならない。成立した主張は限定的である。**同じ50%観測、同じRater露出、同じ接続性、同じ生成responseでも、潜在能力とRater厳しさが結び付く割付は、正常Person分布を仮定するMMLのRater回収とPerson分散推定を変えうる。**

## 事前登録と分母

20反復screeningを報告した後、結果生成前に次を固定した。

- fresh 100反復、241～340
- 2設計 × 100反復 = 200データセット
- FACETS/Python JMLE 200 pair、MML 400 fit、exact CMLE 200 fit
- 合計800 attempt、補充・置換・optional stoppingなし
- primary 1件、secondary 4件の方向とHolm補正
- primary 95%区間半幅目標0.015
- joint/marginal/conditional間のlikelihood・AIC・BIC比較とランキングを禁止

800/800 attemptがexecution/statistical readinessに合格した。FACETS 200/200、MML 400/400、exact CMLE 200/200、制約600/600である。FACETS tryは200、retryは0だった。resume replayは800/800 skip、再fit 0、completion marker set不変だった。

## Primary

| Endpoint | N | 平均差 | 95%区間 | 片側p | 半幅 | 判定 |
|---|---:|---:|---:|---:|---:|---|
| 自由SD MML Rater RMSE、aligned − planned | 100 | +0.09863 | +0.08432～+0.11295 | 7.70e-25 | 0.01431 | 方向確認・精度合格 |

絶対水準では、自由SD MMLのRater RMSE平均はplanned 0.07969、aligned 0.17832 logitsだった。情報的割付で約2.24倍になったが、この比は補助的な記述であり、登録判定はpaired差に基づく。

## Secondary family

| ID | 登録contrast | 平均 | 95%区間 | Holm p | 判定 |
|---|---|---:|---:|---:|---|
| IA2 | 固定SD MML Rater RMSE差 | +0.06851 | +0.05656～+0.08045 | 5.51e-20 | 確認 |
| IA3 | 自由SD MML Person SD差 | −0.11839 | −0.13394～−0.10284 | 1.89e-27 | 確認 |
| IA4 | 自由SD MML 圧縮傾き | −0.48932 | −0.53039～−0.44825 | 3.33e-42 | 確認 |
| IA5 | 固定SD MML 圧縮傾き | −0.38384 | −0.41825～−0.34942 | 6.16e-40 | 確認 |

自由SD MMLの推定Person SD平均はplanned 0.80631、aligned 0.68792だった。生成時の実現SDは両設計とも0.8である。正しいSD=0.8へ固定してもRater RMSE効果が+0.06851残ったため、この現象は自由SD推定の失敗だけでは説明できない。一方、自由化するとRater RMSE差が+0.09863へ大きくなるため、分散推定が感度を増幅することとは整合する。

## Rater厳しさ圧縮

Rater-level truth-error差 `aligned − planned` は次のように再現した。

| 推定 | R01 (−0.45) | R02 (−0.15) | R03 (+0.15) | R04 (+0.45) |
|---|---:|---:|---:|---:|
| 固定SD MML | +0.1647 | +0.0785 | −0.0752 | −0.1679 |
| 自由SD MML | +0.2100 | +0.1009 | −0.0979 | −0.2130 |

寛大側の負値は上方へ、厳しい側の正値は下方へ動く。これは平均位置の一様なずれではなく、Rater間の広がりをゼロ方向へ縮めるパターンである。事前定義した傾きが自由SDで−0.489、固定SDで−0.384となったことが、この機序的な読みを支持する。

## Screeningからの再現性

| 指標 | 20反復screening | fresh 100確認 |
|---|---:|---:|
| 自由SD MML Rater RMSE差 | +0.07575 | +0.09863 |
| 固定SD MML Rater RMSE差 | +0.04109 | +0.06851 |
| 自由SD MML Person SD差 | −0.12478 | −0.11839 |
| 自由SD MML 圧縮傾き | −0.52060 | −0.48932 |
| 固定SD MML 圧縮傾き | −0.40590 | −0.38384 |

5指標すべてがfresh dataで同方向となり、4つのsecondaryはHolm後にも確認された。screening効果量をconfirmationへ混ぜず、選択と検証を分離した。

## FACETS校正と表示精度

| 指標 | 200 pairの最悪値 |
|---|---:|
| 主効果weighted MAE | 0.004551 logits |
| 主効果最大絶対差 | 0.013946 logits |
| PCM閾値weighted MAE | 0.008656 logits |
| PCM閾値最大絶対差 | 0.014610 logits |
| facet内Spearman最小 | 1.0 |

FACETS/Python JMLEは200/200でdirect calibration readyだった。FACETSの通常の小数点以下2桁表示からfitやtruth errorの隠れた精度を復元していない。登録済みUmean=6 measureを直接校正に用い、確認endpointはPythonのネイティブ精度truth-error行から作った。

FACETSはここでJMLEの外部校正器である。MMLやCMLEは別estimandなので、FACETSとの差をそのまま誤差や優劣とみなさない。

## Rによる独立再計算

R 4.5.1がaggregate済みendpoint値を読むだけでなく、`recovery.csv`と`run_ledger.csv`から次を再構築した。

- 5 endpoint × 100ペア = 500 contrast
- Rater RMSE
- 自由MMLのPerson SD差
- 4 Raterレベルからの圧縮傾き
- 片側t検定、95%区間、primary精度判定
- secondary 4件のHolm補正とprimary gate

Pythonとの最大絶対差は `3.55e-15`、数値・論理とも全件一致した。

## JMLE/CMLEの文脈

登録primary/secondaryではない記述的Rater RMSE差は、JMLE +0.00441（95%区間−0.00673～+0.01554）、exact CMLE +0.00367（−0.00630～+0.01364）だった。今回、MMLと同程度のRater圧縮は見られなかった。

ただし、これはJMLE/CMLEの一般的な頑健性や優越性を証明しない。JMLEはPersonを固定効果として扱い、exact CMLEはPerson総得点で条件づけ、MMLは正常Person母集団を積分する。異なる問いに答えているためである。

JMLE閾値RMSE差は記述的に−0.01413（−0.02717～−0.00110）だったが、登録確認endpointではなく、多重性昇格もしない。「情報的割付がJMLE閾値を改善する」という主張にはしない。exact CMLEの極端Personはplanned 4、aligned 1で、全200 fitは適格だった。

## メタ認知的な解釈

今回確認できたのは、単純な欠測率や接続性だけでは不十分だということである。

1. 2設計はともに960行、2 Rater/Person、40 Person/Rater、240行/Raterだった。
2. 同一反復では完全response matrixとPersonベクトルが同じだった。
3. Person–Raterグラフは1成分、制約後nullityは0だった。
4. それでも、能力四分位とRater厳しさを整列するとMMLだけに大きなRater圧縮が再現した。

最も整合的な説明候補は、情報的割付の下で、正常Person母集団の分散とRater厳しさの分離が曖昧になることである。しかし、これは既知の生成機序を使ったsimulation内のモデルベース説明であり、実データの欠測機序を同定したものではない。

現実のデータでは真のThetaが観測されない。推定Thetaと割付の相関をそのままMNAR診断に使うと循環論法になりうる。割付時点の外生共変量、運用ルール、Rater availability、割付確率などが必要である。

## Workbenchへの長期的含意

FACETS補完workbenchは、推定ボタンを増やすだけでは不十分である。推奨する表示順は次のとおりである。

1. **Design audit**: 観測密度、Rater露出、Person–Rater接続性、割付バランスを先に表示する。
2. **Estimand cards**: fixed-person JMLE、normal-population MML、person-total CMLEの問いと仮定を明示する。
3. **Fixed-density stress test**: 行数を固定したまま割付依存性を摂動し、同一推定法内の変化を示す。
4. **Local profile**: global RMSEだけでなくRaterレベル別の移動と推定Person SDを示す。
5. **FACETS calibration panel**: JMLE parityと表示精度契約を独立パネルに置く。
6. **No leaderboard**: 異なる尤度基底のlog-likelihood/AIC/BICや総合順位を作らない。
7. **Correction gate**: MAR/MNAR補正は、外部割付変数・明示したselection model・positivity監査がある場合にだけ別研究として提供する。

この順序なら、ユーザーは「どの推定法が勝ったか」ではなく、「どの仮定と設計摂動に結果が敏感か」を判断できる。ここがFACETSを置き換えずに補完する本アプリの最も価値ある方向である。

## 次の検証

今回の効果を一般化してはならない。次は登録済みfactorial operating-characteristics studyとして、少なくとも次を変える。

- 能力–厳しさ割付の強度
- Rater/Person数とRater厳しさの分散
- Person数とPerson分布形状
- PCM閾値の異質性
- 1～4 Rater/Personとoverlap構造
- 既知の割付確率とpositivity

補正法を評価するなら、既知propensity下のIPW、selection model、joint model等を、未補正推定とは別フェーズで比較する。FACETS parityは引き続きJMLE校正として保持する。

主な証跡:

- `informative_assignment_confirmatory_plan_20260811.json`
- `informative_assignment_confirmatory_execution_registration_20260811.json`
- `informative_assignment_confirmatory100_20260811/preexecution_audit.json`
- `informative_assignment_confirmatory100_20260811/aggregate/`
- `informative_assignment_confirmatory_verify.R`
- `informative_assignment_confirmatory100_20260811/aggregate/r_verification.csv`
- `informative_assignment_confirmatory100_assessment_20260811.json`
