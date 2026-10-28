"""Short offline instructions and CSV checks for the desktop entrypoint."""
from __future__ import annotations

import csv
import io
from pathlib import Path

import numpy as np
import pandas as pd

from mfrm_app import engine
from mfrm_app.io_tables import TABLE_FILE_HARD_LIMIT_BYTES

WEB_APP_URL = "https://mfrmpythonapplication-xodghmcghngx8cz5zmhfeg.streamlit.app"
EXAMPLE_PATH = Path(__file__).resolve().parents[1] / "examples" / "desktop_writing.csv"

CSV_GUIDE = """1行 = 1人を、1人の評価者が、1つの課題・観点について採点した記録です。
同じ人に複数の評定があれば、複数行に分けます（長形式）。

Person,Rater,Task,Criterion,Score
E01,R1,Task1,Content,4
E01,R1,Task1,Organization,3
E01,R2,Task1,Content,5
E02,R1,Task1,Content,2

列の意味
• Person：受験者・学習者など、測定する対象者のID。
• Rater：評価者のID。同じ評価者には同じIDを使います。
• Task / Criterion：課題や評価観点。実際の設計にある要因だけを用意します。
• Score：1回の評定の整数得点（例：0〜5）。平均点・合計点を入力しません。

準備のポイント
• 1行目を列名にし、Excelでは「CSV UTF-8（コンマ区切り）」で保存します。
• 列名は日本語でも構いません。読込後に、各列の役割を選びます。
• 評価者ごとに列が並ぶ表は、上のように評価者IDと得点を行に分けてください。
• 未採点は0点にせず、その行を作らないか得点を空欄にします。
• 同じ対象者・評価者・課題・観点の重複行を確認します。再採点なら採点回等の区別が必要です。
• 最低点・最高点には、今回観測されなかった点も含め、採点基準上の範囲を入力します。
• 評価者や課題の間に共通の対象者・評定がないと、比較できない場合があります。

空欄や文字の得点は計算対象から外れます。分析前に対象外の行数を表示します。
小数・無限大・範囲外の得点は修正してから分析します。
この初期版は各評定を等しく扱います。重み列を追加しても重み付き分析にはなりません。

「見本CSVを保存」は、Web版と同じ作文の人工データを保存します。
見本は操作・図の読み方の練習用であり、実在の学習者や評価者のデータではありません。
"""

SETTINGS_HELP = """列
対象者：受験者や学習者のIDの列。
得点：1回の評定の整数得点の列。
評価者・課題など：分析に含める要因の列をクリックして選択します。複数選べます。
同じ列を対象者・得点・要因に重ねて指定することはできません。

採点基準
ルーブリックの最低点と最高点を入力します。
例えば採点基準が0〜5点なら、今回0点や5点が出ていなくても最低点0、最高点5です。
「データ内の得点」は入力した得点の一覧で、採点基準を自動決定する値ではありません。

入力の問題があれば、分析前に該当箇所を表示します。
空欄等による対象外の行、同じ組合せの重複がある場合だけ、続行を確認します。
"""

QUICK_HELP = """データ → 設定 → 結果
「CSVを開く」または「サンプル」でデータを読み込み、「次へ」で設定に進みます。
列と採点基準を指定し、「分析」を押すと結果に切り替わります。
結果では「図」「推定値」「得点」「データ概要」を選べます。各図の説明は「読み方」にあります。

保存
「プロジェクト」は入力データ・識別子と結果を含み、再起動後も再推定せず開けます。
「レポート」は図と表を含むHTMLです。図はPNG／SVGでも保存できます。
データや設定に戻っても結果は保持します。設定変更後の結果は前回の分析であることを表示します。

この版でできる分析
RSM（共通の評定尺度）をJMLEで推定します。評価者等の効果は加法的に扱います。
推論の保留が表示された場合は、図が描けても差や妥当性を結論にしないでください。

数式・理論を調べるには
下の「Web版の理論・ヘルプ」を開き、Web版のHelp（ヘルプ）または
Tutorial（チュートリアル）で、MFRM、推定方法、適合度、各図の説明を参照します。
リンクはWeb版の入口を開きます。デスクトップのデータや結果を送信する機能ではありません。
ネット接続が必要です。公開Webとこの開発版では、表示や対応機能が異なる場合があります。

利用条件
研究・教育・業務・商用に利用できます（MITライセンス）。
データの取扱い、分析の適用判断、結果の解釈・利用は利用者の責任です。無保証で提供します。
再配布時は著作権・ライセンス表示を保持してください。第三者の構成要素には各ライセンスが適用されます。
"""


def writing_example():
    return read_desktop_csv(EXAMPLE_PATH)


def read_desktop_csv(path):
    """Read UTF-8 long-form CSV without converting identifier strings to numbers."""
    path = Path(path)
    with path.open("rb") as stream:
        raw = stream.read(TABLE_FILE_HARD_LIMIT_BYTES + 1)
    if len(raw) > TABLE_FILE_HARD_LIMIT_BYTES:
        raise ValueError("CSVが読込上限（100 MiB）を超えています。")
    try:
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("文字コードを読み取れません。Excelで「CSV UTF-8（コンマ区切り）」として保存し直してください。") from exc
    rows = csv.reader(io.StringIO(text), strict=True)
    try:
        header = next(rows)
        if len(header) < 3 or any(not name.strip() for name in header):
            raise ValueError("1行目に対象者・評価者等・得点の列名が必要です。CSVの形式と見本を確認してください。")
        if len(set(header)) != len(header):
            raise ValueError("同じ列名が重複しています。それぞれ異なる名前にしてください。")
        for number, row in enumerate(rows, start=2):
            if row and len(row) != len(header):
                raise ValueError(f"CSVの{number}行目の列数が1行目と違います。区切りや引用符を確認してください。")
    except (csv.Error, StopIteration) as exc:
        raise ValueError("CSVを読み取れません。空のファイルや引用符の閉じ忘れがないか確認してください。") from exc
    frame = pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)
    if frame.empty:
        raise ValueError("列名だけでなく、評定データの行も必要です。")
    return frame


def input_review(data, settings):
    """Use the engine's own row audit; explain exclusions before starting a worker."""
    person, score = settings["person_col"], settings["score_col"]
    facets = list(settings["facet_cols"])
    columns = [person, *facets, score]
    if not person or not score or not facets:
        raise ValueError("対象者と得点の列、および評価者・課題等を1つ以上選択してください。")
    if len(set(columns)) != len(columns) or any(c not in data for c in columns):
        raise ValueError("対象者・得点・ファセットには、それぞれ別の列を選択してください。")
    minimum, maximum = settings["rating_min"], settings["rating_max"]
    if minimum >= maximum:
        raise ValueError("採点基準の最低点より大きい最高点を入力してください。")
    scores = pd.to_numeric(data[score], errors="coerce")
    present = scores.notna()
    if (present & ~np.isfinite(scores)).any():
        raise ValueError("得点に無限大があります。有限の整数得点に修正してください。")
    fractional = present & ((scores - np.round(scores)).abs() > np.sqrt(np.finfo(float).eps))
    if fractional.any():
        raise ValueError("得点に小数があります。平均点ではなく、1回ごとの整数の評定を入力してください。")
    audit = engine.build_response_data_audit(data, person_col=person, facet_cols=facets, score_col=score)
    excluded = len(audit["excluded_rows"])
    retained = data.loc[audit["rows"]["IncludedInLikelihood"].to_numpy(dtype=bool)]
    kept_scores = pd.to_numeric(retained[score], errors="coerce")
    if kept_scores.nunique() < 2:
        raise ValueError("計算対象となる得点が2種類以上必要です。列の対応や空欄を確認してください。")
    if ((kept_scores < minimum) | (kept_scores > maximum)).any():
        raise ValueError(f"得点に指定範囲（{minimum}〜{maximum}）の外の値があります。得点列と採点基準を確認してください。")
    duplicates = int(retained.duplicated([person, *facets], keep=False).sum())
    messages = []
    if excluded:
        messages.append(f"空欄や数値でない得点などにより、{len(data):,}行のうち{excluded:,}行が対象外になります。")
    if duplicates:
        messages.append(f"同じ対象者・ファセットの組合せが重なる行が{duplicates:,}行あります。重複コピーや採点回の区別を確認してください。")
    return {"included": len(retained), "excluded": excluded, "duplicate_rows": duplicates, "messages": messages}
