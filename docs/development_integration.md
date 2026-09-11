# 開発系統の統合台帳

更新日: 2026-09-11。方針: 一つのアプリへ段階的に統合する。

## 現在の作業場所と公開版

- 今後の開発場所: `MFRM_Python_Application`。
- 開発ブランチ: `integration/unified-app`。
- 公開ブランチ: `main`。確認時の GitHub HEAD は
  `d8558fc6b4c05d96e4c13a0d067212fde6423b59`、公開表示は `0.2.16-beta`。
- 開発側の `0.2.15-beta` は継承した表示であり、公開版との新旧比較には使わない。
  統合版のバージョン変更は公開準備時に行う。
- `MFRM_Python_Application_working` は原本保全用。今後の変更は開発ブランチに集約する。
  独自の画面・診断・外部連携をすべて移したわけではないため、削除しない。
- GitHub への push、公開 `main` への merge、デプロイは今回行っていない。

## 保全した履歴

| 参照 | 内容 |
|---|---|
| `archive/public-main-20260911` | 確認した公開版 `d8558fc` |
| `archive/local-main-20260911` | 整理前のローカル `main`、`e729f19` |
| `archive/analysis-20260911` | 分析系統の未コミット実装・テスト・検証プロトコルを含む `3cec72d` |
| `archive/working-20260911` | `_working` の未コミット変更も含む `b4dd534` |

二系統の共通祖先は `ce4dd03`。公開版と分析系統の単純な merge は11ファイル、
`_working` と分析系統の merge は15ファイルで競合した。これはファイル数であり、
解消すべき論点や競合ブロックの数ではない。いずれも自動統合済みとは扱わない。

退避先はローカルの `backups/consolidation_20260911/`。
各系統の `history.bundle`、`files.tar.gz`、`manifest.json`、`tracked.patch`、
`status.z` を保存した。manifest は各ファイルのサイズと SHA-256 を保持し、
tar 内の全ファイルを読み直して一致を確認した。

| 原本 | 退避したファイル数 | 合計バイト数 |
|---|---:|---:|
| 分析系統 | 95,371 | 1,603,599,549 |
| `_working` | 149 | 5,859,219 |

対象は Git の追跡ファイルと非無視の未追跡ファイル。仮想環境、キャッシュ、
無視済み秘密設定などを含むディレクトリ全体のバックアップではない。
原本ファイルは移動・削除していない。

分析系統の source snapshot は560ファイルを保存している。
残る94,811ファイルは主に検証結果であり、元の `validation/` と tar に保持する。
その既存結果ディレクトリと vendor packet ZIP のみをローカル
`.git/info/exclude` に列挙し、ソース変更と実験結果を分けた。
この除外は公開リポジトリには伝播しない。

## 統合状況と順序

| 段階 | 状態 | 内容と完了条件 |
|---|---|---|
| 0: 保全・作業系統の整理 | 完了 | 公開版を固定、原本と未コミット成果を検証付きで退避、開発ブランチを一本化 |
| 1: 事前設計の計算基盤 | 完了 | `mfrm_app/simulation/` 14ファイルと既存テスト6ファイルを取り込み、172テスト通過 |
| 2: 公開版の改善を照合 | 未着手 | 公開版の遅延読み込み、Standard実行の計算量制御、設定保存、ダウンロード生成を分析系統へ適用。統計的な出力保護・AnalysisID・日英画面の回帰を確認 |
| 3: 事前設計の画面 | 未着手 | `_working` の `mfrm_app/ui/design_planner.py` を共通の入口・Help・状態管理へ接続。画面テストと実ブラウザ確認を通す |
| 4: 適合度計算と表示 | 未着手 | `mfrm_app/diagnostics/` の式・極端反応・raw/display方針を現在の判定規則と照合。既存推定・raw値・丸め前分類を変えないことを検証 |
| 5: 外部検証への書き出し | 未着手 | `mfrm_app/external_validation/` の privacy・ConQuest 契約を共通ダウンロードへ統合。公開テンプレートと非公開データの境界を維持し、アプリから外部推定器を実行しない |
| 6: リリース資格 | 未着手 | clean checkout でのCI、検証成果の復元手順、統計ゲート、アクセシビリティ、公開する機能・版番号を確定 |

段階1は `b4dd534` からファイル内容を変更せず移植した。既存の分析コード、
画面、推定法セレクターには接続変更を加えていない。事前設計基盤の内容は、
評定負担・費用、決定的割付、接続性、構造的な追加評定案、科学条件と乱数の契約。
得点生成、反復推定、標準誤差改善、必要サンプルサイズ推奨は含まない。
画面と統計推論の完成を、計算基盤のテスト成功から推定しない。

元の設計説明は `archive/working-20260911:docs/sample_size_design_mvp.md`、
外部連携の説明は同ブランチの `docs/conquest_external_handoff.md` を参照する。
未統合の実装はそのブランチで比較・復元できる。

## 統計・検証の現在地

- 出力保護 G0 と同定座標 G1 は完了記録があり、G2 以降の資格を統合で解除しない。
- exact CMLE は研究用のまま。事前設計モジュールの追加で公開セレクターに出さない。
- `kac200` は200独立Personベクトル・600データセットの確認試験まで完了。
  v1登録基準では方向と精度の判定を満たしたが、追加監査では Python MML の
  最適化停留性の資格が未完了とされた。過去の結果をv2で置換しない。
- 自由分散MML v2 は開発段階。独立した新規データでの資格検証が必要。
- 実ブラウザのアクセシビリティ受入は未完了。

詳細は [MML v2](../validation/MML_FREE_SD_STATIONARITY_V2_ROADMAP_20260811.md)、
[統計エンジンの修正計画](statistical_engine_remediation_roadmap.html)、
[ブラウザ受入](browser_accessibility_acceptance.md) に従う。
過去の README / ROADMAP にある「次の確認試験」という記述より、対応する
最新の実行・監査記録を優先する。本台帳は開発統合状況の入口であり、統計ゲートの代替ではない。

## 再確認用コマンド

```bash
git branch --show-current
git log -1 --oneline main
git diff --stat archive/public-main-20260911 integration/unified-app
git diff --stat archive/working-20260911 integration/unified-app
python3 -B -m pytest -q -p no:cacheprovider \
  tests/test_design_spec_cost.py tests/test_rating_design_assignments.py \
  tests/test_rating_design_topology.py tests/test_design_planning.py \
  tests/test_simulation_conditions.py tests/test_structural_transition_preflight.py
python3 streamlit_app.py --doctor
```

既存の研究テストにはローカルに保全した `validation/` の結果を読むものがある。
source snapshot だけの新規 checkout ですべての研究テストが通るとは主張しない。
必要な凍結fixtureの選別・配布・復元は段階6で扱う。
