# 既知割付機構：80×4大設計の校正

実施日: 2026-08-11  
最終判定: **PASS（exact-DP lane）**  
選択されたresponse-study dose: **`gamma = 0, ±0.8`**

## なぜ小規模oracleだけで終えなかったか

6 Persons×3 Ratersの90状態oracleは遷移核の定常分布を検証できるが、80 Persons×4 Ratersでの有限連鎖混合を保証しない。そこでScoreも推定結果も存在しない段階で、固定Person次数2、固定Rater次数40、160 edgesの設計を校正した。

Person座標は1つの凍結normal draw、Rater座標は真のseverity `[-0.45, -0.15, 0.15, 0.45]` とした。候補は `gamma=-0.8,-0.4,0,+0.4,+0.8`。dose選択は割付相関、混合、連結性、周辺次数だけで行い、後のbiasを見て選ばない。

## 保持された2つの失格

### 1. edge-pair Metropolis 2-switch

- 4 chains/gamma、各3,000標本、burn-in 10,000、thin 5
- 全次数・連結・受理率・split-R-hatは合格
- 最小Statistic ESSは **13.79**（登録下限300）
- 判定: **FAIL、dose選択なし**

### 2. Person-pair heat-bath Curveball

90状態oracleで同じ目標分布を再確認した後、開始状態、gamma、標本数、burn-in、thin、ESS下限を変えずにblock Gibbsへ変更した。

- 最大split-R-hat 1.0171、最小movement rate 0.3141
- 最小Statistic ESSは **33.34**（登録下限300）
- 判定: **FAIL、dose選択なし**

R-hatと受理／移動率だけでは、割付統計量の長い自己相関を除外できなかった。閾値は緩めず、両失格bundleを保持した。

## exact dynamic-programming remediation

4 Raters・既知周辺次数という限定を利用し、`(Person位置, 残余Rater次数)` を状態とする後向きDPで分配関数を計算した。そのfull conditionalから各PersonのRater集合を逐次抽出し、非連結グラフのみを棄却する。したがって受理標本は固定次数・連結条件付き分布からの独立標本である。

gammaごとの分配関数構築が約100秒かかるため、最初のmonolithic実行2回は300秒上限で出力前に停止した。科学条件を変えず5 gamma shardへ分割する実行amendmentを先に記録し、その後のみaggregateした。

| gamma | 独立標本 | 平均割付相関 | split-R-hat | 連結棄却 |
|---:|---:|---:|---:|---:|
| -0.8 | 4,000 | -0.40782 | 1.00021 | 0 |
| -0.4 | 4,000 | -0.24529 | 0.99932 | 0 |
| 0.0 | 4,000 | +0.00062 | 1.00025 | 0 |
| +0.4 | 4,000 | +0.24642 | 1.00032 | 0 |
| +0.8 | 4,000 | +0.40824 | 0.99935 | 0 |

- 全gammaでDP statesは985,681／上限2,000,000。
- 20 batchesの最小ESS診断は576.41／下限300。
- 全20,000標本生成attemptで連結棄却は0。
- 全gammaでPerson/Rater次数、direct overlap連結性、Rater×Task×Criterion度数、Score-free materializationが合格。
- 割付相関はgammaに対して厳密に単調。
- `|gamma|=0.8` の相関絶対値は登録範囲0.20–0.60内で、正負対称誤差は0.00043／上限0.08。

事前ルール「全ゲートを通る最大の候補dose」により、`|gamma|=0.8` を選択した。機械可読結果は [assessment.json](known_assignment_large_design_dp_20260811/assessment.json)、実行分割の理由は [known_assignment_large_design_dp_execution_amendment_20260811.json](known_assignment_large_design_dp_execution_amendment_20260811.json) にある。

## 解釈と長期設計

この結果はMCMC一般の失敗を意味しない。登録した有限連鎖設定がこの統計量に不十分だったことを意味する。4-Rater response studyではexact-DPがより強い選択肢だが、Rater数が増えると状態空間が急増する。DP capを超える設計では、より長いMCMC、別のblock kernel、または近似誤差を明示した別laneが必要である。

`gamma=0.8` は現実の割付機構の推定値ではない。中程度に強い既知simulation stress doseである。実データUIには昇格させず、次の凍結PCM response studyでFACETS 4.5 JMLE／Python JMLE／MML／exact CMLEのwithin-estimator design sensitivityを調べる。

FACETSはこの割付校正には関与しない。次段階でもFACETSとPythonの直接parityはJMLE同士に限定し、通常の小数第2位fit表示をraw計算へ戻さない。

## 証拠管理上の補足

最初のoracle bundleとMetropolis大設計bundleは、manifestにimport module hashが欠けていた。後続追加前のruntime moduleを再構築したSHA-256と影響範囲を [known_assignment_dependency_disclosure_20260811.json](known_assignment_dependency_disclosure_20260811.json) に開示した。heat-bath／exact-DP以降はmodule dependency hashをmanifestへ含めている。
