# 既知の割当機構: 密配列 exact-DP 資格化

実施日: 2026-08-11  
判定: **PASS（実装経路の資格化のみ）**

## 目的

既に資格化された固定周辺度数の割当分布

\[
P(G\mid d_P,d_R)\propto\exp\{\gamma T(G)\}
\]

について、再帰・メモ化DPと同じ分配関数を、4-Raterの小さな周辺度数空間に
特化した密配列DPで計算できるかを検証した。これは割当モデル、`gamma`、
estimand、応答生成、推定法を変更するものではない。

計画・入力・実装のhashと許容差は実行前に
[`known_assignment_dense_dp_plan_20260811.json`](known_assignment_dense_dp_plan_20260811.json)
へ固定した。

## 三重照合

1. 6 Person × 3 Rater の全90状態oracle
2. 既存の再帰・メモ化DP
3. 80 Person × 4 Rater の既存資格化デザイン

小規模oracleでは各 `gamma=-0.8, 0, +0.8` から12,000個の独立標本を生成した。
oracle外確率はすべて0、分配関数の再帰DPとの差は最大
`8.88e-16`、全変動距離は最大`0.03483`で、事前登録上限`0.06`を通過した。

80×4では、5条件すべてで再帰DPとの差が最大`2.84e-14`だった。固定Person
次数2、固定Rater次数40、連結性、最終統計量の再計算も全条件で一致した。

| gamma | recursive log Z | dense log Z | absolute difference | partition seconds |
|---:|---:|---:|---:|---:|
| -0.8 | 165.941831190697 | 165.941831190697 | 0 | 0.2014 |
| -0.4 | 144.575905306945 | 144.575905306945 | 2.84e-14 | 0.1982 |
| 0.0 | 136.342612878623 | 136.342612878623 | 0 | 0.1973 |
| +0.4 | 144.575905306945 | 144.575905306945 | 0 | 0.1980 |
| +0.8 | 165.941831190697 | 165.941831190697 | 0 | 0.1990 |

## 性能境界

80×4の配列は`41×41×41`、全81層で5,582,601セル、44,660,808 bytesである。
最大構築時間は0.202秒だった。既存の同一マシン上の `gamma=+0.8` 再帰DP診断
92.68秒に対し約466倍で、事前登録した5倍以上を通過した。

この倍率は一般的なハードウェア性能主張ではない。密配列は明示したセル上限を
超えるとfail closedする。Rater数や周辺度数が大きい問題へ無条件には拡張しない。

## 大局的な意味

従来は1つのPersonベクトル・1つの`gamma`について分配関数を作るだけで約100秒
かかり、複数のPerson母集団drawを扱う確認研究が実務上困難だった。今回の変更で、
「1つの固定Personベクトルに条件づけたscreening」を繰り返すのではなく、独立した
Personベクトルを研究上の変動源として扱える。

ただし、計算可能性が科学的妥当性を与えるわけではない。次段階では新しいPerson
ベクトルを用いた非確認的preflightを行い、応答生成、FACETS/Python JMLE校正、
MML/CMLE制約、endpoint符号のPerson-vector依存性を確認してから、確認研究の標本数と
primary endpointを固定する。

## FACETS精度境界

この資格化には応答、Score、推定結果、FACETS出力を一切使用していない。
したがってFACETSの通常表示で小数第2位に丸められるfit値も入力されていない。
後続研究でも、FACETSは同一estimandのJMLE外部校正に限定し、通常表示fitから未報告の
raw residualやfitを復元しない。

主なartifact:

- [assessment.json](known_assignment_dense_dp_20260811/assessment.json)
- [small_oracle_equivalence.csv](known_assignment_dense_dp_20260811/small_oracle_equivalence.csv)
- [large_design_equivalence.csv](known_assignment_dense_dp_20260811/large_design_equivalence.csv)
- [manifest.json](known_assignment_dense_dp_20260811/manifest.json)

