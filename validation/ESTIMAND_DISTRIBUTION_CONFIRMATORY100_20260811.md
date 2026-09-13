# Person分布形状 confirmatory 100（2026-08-11）

## 結論

screening 20から選んだ6仮説を、重複しないreplicate 21–120のfresh dataで固定標本100として検証した。確証されたのは、連結計画欠測下の自由SD MMLで、正規Personに比べて右歪みPersonとt(3)重尾Personの推定母集団SDが低くなる2仮説である。対称二峰混合のC02 category 2にscreeningで見えた局所的な負方向は、Python JMLE・自由SD MML・exact CMLEの3法すべてで再現せず、点推定も小さな正方向へ反転した。

計画欠測が正常分布条件のJMLE閾値RMSEを増やすH6は、平均差 `+0.06594`、95%区間 `[0.05772, 0.07416]` と明瞭だった。しかし、FACETSが1回だけexit 0のままU6 reportを生成しなかったため、登録済み100-pair gateに対して99 pairとなった。したがってH6は、統計量が強くてもconfirmatoryには **inconclusive** とする。

厳格なFACETS補完ワークベンチ資格も、3,200 attempt中1件の失敗により不合格である。一方、適格だった799 FACETS/Python JMLE pairは799/799でdirect agreementを通過した。数値校正と外部engineのbatch I/O信頼性は、別の品質レイヤーとして扱う必要がある。

## 事前登録と独立性

- screening: replicate 1–20
- confirmatory: replicate 21–120
- confirmatory標本数: 固定100、optional stoppingなし
- DGM: 正しく指定した4カテゴリ・異質閾値PCM
- Person: 各実現ベクトルを平均0、母SD 0.8へ厳密に標準化
- 形状: normal / right_skew / symmetric_mixture / heavy_tail_t3
- 設計: complete / planned_connected
- 推定: FACETS 4.5 JMLE、Python JMLE、固定／自由SD Q31 MML、native exact CMLE
- 検定: 6つの片側paired t検定にHolm補正、familywise alpha 0.05
- pair gate: 各仮説に有限な100 pairが必須、補完・仮説差替え・方向反転なし

planはconfirmatory入力生成・fit・aggregate・分析より前に固定した。screeningの20反復はprimary検定へ一切混入していない。全800データセットはrank-fullで、confirmatory seedはscreening seedと非重複だった。

## Primary family

| ID | 登録contrast | N | 平均差 | 95%区間 | Holm p | 精度 | confirmatory判定 |
|---|---|---:|---:|---:|---:|---|---|
| H1 | free-SD: right_skew − normal、planned | 100 | -0.08983 | [-0.10380, -0.07585] | 3.24e-22 | pass (HW 0.01398 ≤ 0.02) | 方向確証 |
| H2 | free-SD: heavy_tail_t3 − normal、planned | 100 | -0.06502 | [-0.08532, -0.04472] | 1.28e-08 | **miss** (HW 0.02030 > 0.02) | 方向確証、精度missを併記 |
| H3 | JMLE: mixture − normal、planned C02 step 2 | 100 | +0.01600 | [-0.02206, +0.05407] | 1.0 | pass | 負方向は再現せず |
| H4 | free-MML: mixture − normal、planned C02 step 2 | 100 | +0.01044 | [-0.02794, +0.04882] | 1.0 | pass | 負方向は再現せず |
| H5 | exact-CMLE: mixture − normal、planned C02 step 2 | 100 | +0.01501 | [-0.02283, +0.05285] | 1.0 | pass | 負方向は再現せず |
| H6 | JMLE threshold RMSE: planned − complete、normal | 99 | +0.06594 | [+0.05772, +0.07416] | 1.88e-28（記述のみ） | pass | **100-pair gate不合格、inconclusive** |

H2は方向検定を通過したが、95%半幅が登録目標を `0.0003005` 上回った。方向確証と精度資格を混同せず、「方向は確証、目標精度は未達」と報告する。

## Screeningとの対照

| Endpoint | screening 20 | confirmatory | 判断 |
|---|---:|---:|---|
| H1 right-skew free-SD | -0.08620 | -0.08983 | 大きさも方向も再現 |
| H2 heavy-tail free-SD | -0.06552 | -0.06502 | 大きさも方向も再現 |
| H3 mixture JMLE local threshold | -0.11452 | +0.01600 | 非再現・符号反転 |
| H4 mixture free-MML local threshold | -0.12294 | +0.01044 | 非再現・符号反転 |
| H5 mixture exact-CMLE local threshold | -0.11534 | +0.01501 | 非再現・符号反転 |
| H6 design threshold RMSE | +0.08092 | +0.06594 (N=99) | 大きさ・方向は整合、登録gateで未確証 |

screeningで5モードに共通して見えた局所threshold patternさえ、fresh 100では再現しなかった。この結果は、多モード同方向という見た目だけでは、選択後推論と多重な局所探索の影響を克服できないことを示す。アプリの警告や自動判断へscreening patternを直結させてはならない。

一方、自由SD MMLのplanned条件平均は、normal `0.79925`、right_skew `0.70942`、heavy_tail_t3 `0.73422`だった。生成時の実現SDは全て0.8である。これは「MMLが劣る」というランキングではなく、正規Person母集団を仮定して自由SDを推定するestimandが、同じ第1・第2モーメントでも分布形状に感応することの確証である。

## FACETS 4.5.0校正と運用失敗

3,200 completion markerは全て作成され、3,199 attemptが成功した。成功した799 FACETS/Python pairの校正は次のとおりだった。

- direct agreement: 799/799
- main weighted MAEのrun最大: 0.004989
- main最大絶対差: 0.013272 logits
- threshold weighted MAEのrun最大: 0.009313
- threshold最大絶対差: 0.015546 logits
- facet内Spearman最小: 1.0

唯一の失敗は `complete__normal::rep-00053::FACETS_PYTHON_JMLE_PCM` である。FACETSはexit 0を返したが、登録済みU6 reportを作らず、stdout/stderrも空だった。これは数値的不一致や統計的非収束ではない。

primary分析後、元markerを変更せず、confirmatoryへ戻さない診断planを別途登録した。同じ入力を孤立serialで3回再生すると3/3成功し、3/3でdirect agreementを通過した。よって一過性report I/O障害という分類は支持されるが、並列性が原因とは断定できない。診断値でH6を100へ補充せず、workbench資格も不合格のままとした。

FACETS表示精度契約も維持した。measure/thresholdは `Umean=0,1,6`、Table 8 fit/countは別の `Umean=0,1,2` passを使う。小数第2位程度に丸められた表示fitから未報告raw値を逆算していない。primaryのJMLE thresholdは、FACETSと校正されたunrounded Python JMLEを用いた。

## 再開・同一性・cross-engine検証

- 8 shardの壁時計時間: 8,111.1秒（135.2分）
- completion: 3,200/3,200
- success: 3,199/3,200
- resume: 3,200 valid skip、再fit 0、failure 0
- marker結合SHA-256: resume前後とも `725e6f98ec24a52812bbc6689009fcf2a1f18f441b9842aaa40bca7274e75dd4`
- R 4.5.1再計算: paired mean、SD/SE、t、95% CI、片側p、Holm、100-pair論理が全一致
- R/Python最大数値差: `1.78e-15`

CMLEは800/800が推論適格で、64,000 Person-run中113件（0.1766%）が極端総得点、111/800 runで少なくとも1件だった。これはsecondary descriptive outputであり、推定法の優劣や失敗率比較には使わない。

## 大局的・長期的判断

1. **確証された製品価値はestimand sensitivityである。** FACETSと同じJMLEを再現するだけでなく、正規母集団MMLの自由SDが分布形状にどう反応するかを、同じ生成データで可視化できる。
2. **方法間ランキングは禁止したままにする。** JMLEはPerson固定効果joint likelihood、MMLは母集団分布を置くmarginal likelihood、exact CMLEはPerson総得点で条件化する。likelihood/AIC/BICを基底横断で比較しない。
3. **局所screening signalは自動ルールにしない。** H3–H5の非再現は、screeningを明示的にscreeningへ留めた判断の正しさを裏付ける。
4. **観測設計は依然として高価値な感度軸である。** H6のpoint estimateはscreeningと整合したが、規則どおり未確証である。必要ならreplicate 121以降の新規H6-only studyを別planで行い、失敗runの置換はしない。
5. **外部engine運用を統計から分離する。** 将来のversioned runnerでは、Python JMLE artifactをFACETS report成功と独立に永続化し、exit-0/missing-reportだけを対象とする回数制限retryまたはFACETS直列化を事前登録する。今回のstudyへ遡及適用しない。
6. **資格は二層で表示する。** 「799/799 numerical calibration pass」と「800-run batch operational qualification fail」を同時に示し、片方で片方を隠さない。

## 主要成果物

- plan: `estimand_distribution_confirmatory_plan_20260811.json`
- preexecution audit: `estimand_distribution_confirmatory100_preexecution_audit_20260811.json`
- retained study: `estimand_distribution_confirmatory100_20260811/`
- aggregate: `estimand_distribution_confirmatory100_20260811/aggregate/`
- primary analysis: `estimand_distribution_confirmatory100_20260811/confirmatory_analysis/`
- assessment: `estimand_distribution_confirmatory100_assessment_20260811.json`
- FACETS I/O diagnostic: `estimand_distribution_confirmatory100_20260811/facets_io_diagnostic/`
- Python analysis: `estimand_distribution_confirmatory_analysis.py`
- R verification: `estimand_distribution_confirmatory_verify.R`

## Claim limit

confirmatory表現は、この登録DGM・設計におけるH1とH2に限る。H6はinconclusive、H3–H5は非再現、厳格workbench資格は不合格である。FACETSの全面互換、一般的分布頑健性、fit raw parity、推定法ランキング、基底横断likelihood比較は主張しない。
