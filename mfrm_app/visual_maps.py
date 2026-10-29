"""Shared map coordinates and native, offline figures. Never refits a model."""
from __future__ import annotations

import numpy as np
import pandas as pd
import textwrap

from mfrm_app import output_qualification

COLORS = ("#2176ae", "#d88627", "#8064a2", "#27988d", "#a34d72", "#6e7d35")
MAP_GUIDES = {
    "wright": (
        "Wright map",
        "左は対象者の測定値の分布、右は評価者・課題などの位置です。共通の縦軸（logit）で順序と間隔を読みます。"
        "通常の設定では上ほど対象者は高い測定値、評価者は厳しく、課題は難しくなります。\n\n"
        "ファセットの原点は同定制約に依存します。位置の重なりだけで天井・床効果や妥当性を結論にせず、"
        "得点分布、カテゴリー境界、測定精度を合わせて確認します。保留されたSEは利用しません。\n\n"
        "位置の違いを読んでから、Pathway mapで適合度も確認できます (Boone, 2008)。\n"
        "Boone (2008), Teaching students about Rasch maps: https://www.rasch.org/rmt/rmt222a.htm",
    ),
    "pathway": (
        "Pathway map",
        "縦軸は測定値、横軸は選択した適合度です。点の大きさは推定誤差を表しません。\n\n"
        "Infit MnSq：初期表示。1を基準に、予測に対する評定のばらつきの大きさを見ます。"
        "大きい値は予測以上、小さい値は予測以下のばらつきを示します。固定の合格帯は設けません。\n\n"
        "近似Infit ZSTD：0を基準とする標準化値で、±2は参照線です。評定数・自由度・推定方法・"
        "多重性に依存するため、範囲内を適合の証明、範囲外を自動的な不適合としません。"
        "変換方法と打切りの有無は図の注記で確認できます。\n\n"
        "MnSq、評定数、元の評定と内容を合わせて点検します。MnSq自体の基準分布も普遍的ではありません。"
        "データの種類と推定法に合う解釈が必要です (Müller, 2020)。\n"
        "Müller (2020), Item fit statistics for Rasch analysis: can we trust them?: https://doi.org/10.1186/s40488-020-00108-7",
    ),
    "category": (
        "確率曲線",
        "横軸は対象者の測定値から評価者・課題等の効果を差し引いた位置、縦軸は各得点の予測確率です。"
        "中間カテゴリーが最も選ばれやすい範囲の有無を、実際の度数、閾値の順序、採点基準の意味と合わせて確認します。"
        "ピークが見えないことだけでカテゴリーの誤用や統合を結論にしません (Linacre, 2002b)。\n"
        "Linacre (2002), Optimizing rating scale category effectiveness: https://pubmed.ncbi.nlm.nih.gov/11997586/",
    ),
}

MAP_GUIDES_EN = {
    "wright": (
        "Wright map",
        "Read the ordering and spacing of person locations (left) and facet locations (right) on the shared logit axis. "
        "Under the usual coding, higher locations mean higher person measures, greater rater severity or greater task difficulty. "
        "Facet origins depend on identification constraints. Overlap alone does not establish targeting, floor/ceiling effects or validity; "
        "inspect score distributions, category boundaries and qualified precision estimates together.\n\n"
        "After reading locations, use the Pathway map to inspect fit (Boone, 2008).\n"
        "Boone (2008), Teaching students about Rasch maps: https://www.rasch.org/rmt/rmt222a.htm",
    ),
    "pathway": (
        "Pathway map",
        "The vertical axis shows measures; the horizontal axis shows the selected fit statistic. Point size does not encode uncertainty.\n\n"
        "Infit MnSq (default): compare response variation with the model reference of 1. Larger values indicate more variation, "
        "smaller values less. No fixed acceptance band is imposed. Approximate Infit ZSTD uses 0 and reference lines at ±2; "
        "counts, degrees of freedom, estimation and multiple comparisons affect interpretation. The figure note identifies the transformation and any cap.\n\n"
        "Inspect MnSq, counts and the original ratings and content together. MnSq and ZSTD are not universally calibrated tests. "
        "Interpretation must account for the data and estimation method (Müller, 2020).\n"
        "Müller (2020), Item fit statistics for Rasch analysis: can we trust them?: https://doi.org/10.1186/s40488-020-00108-7",
    ),
    "category": (
        "Category probability curves",
        "The horizontal axis is the person location adjusted for rater, task and other facet effects; the vertical axis is category probability. "
        "Inspect whether intermediate categories become most probable alongside observed frequencies, threshold order and rubric meaning. "
        "An absent peak alone does not establish category misuse or justify collapsing categories (Linacre, 2002b).\n"
        "Linacre (2002), Optimizing rating scale category effectiveness: https://pubmed.ncbi.nlm.nih.gov/11997586/",
    ),
}


def wright_map_values(person_table, facet_table):
    """Normalize finite coordinates identically for the Web and desktop maps."""
    person_table = person_table if isinstance(person_table, pd.DataFrame) else pd.DataFrame()
    persons = pd.to_numeric(person_table.get("Estimate", pd.Series(dtype=float)), errors="coerce")
    persons = persons[np.isfinite(persons)]
    facets = facet_table.copy() if isinstance(facet_table, pd.DataFrame) else pd.DataFrame()
    if "Estimate" in facets:
        facets["Estimate"] = pd.to_numeric(facets["Estimate"], errors="coerce")
        facets = facets.loc[np.isfinite(facets["Estimate"])].copy()
    return persons, facets


def pathway_map_values(diagnostics, statistic="InfitZSTD"):
    if statistic not in {"Infit", "InfitZSTD"}:
        raise ValueError("Choose Infit or InfitZSTD")
    measures = diagnostics.get("measures")
    required = {"Estimate", statistic, "Facet", "Level"}
    if not isinstance(measures, pd.DataFrame) or not required.issubset(measures):
        return pd.DataFrame(columns=sorted(required))
    frame = measures.copy()
    for name in ("Estimate", statistic):
        frame[name] = pd.to_numeric(frame[name], errors="coerce")
    return frame.loc[np.isfinite(frame["Estimate"]) & np.isfinite(frame[statistic])].copy()


def map_caution(result, lang="en"):
    reason = output_qualification.calibration_inference_hold_reason(result)
    if not reason:
        boundary = (result.get("person_boundary") or {}) if isinstance(result, dict) else {}
        if isinstance(boundary, dict) and boundary.get("status") == "extreme_persons_present":
            return ("最低点・最高点だけの人の位置は、有限のJMLEではなく計算の終端値です。個人の位置・適合度・分布は技術的な点検用として扱ってください。"
                    if lang == "ja" else "Extreme-score person locations are terminal calculation values, not finite JMLEs. Person locations, fit and distributions are for technical review only.")
        return ""
    if str(reason).startswith("stat.mml_"):
        return ("MMLの数値精度と区間推定の検証が十分ではありません。表示値だけから差の有無や順位を結論づけないでください。"
                if lang == "ja" else "MML numerical accuracy and interval estimates are not sufficiently validated. Do not use these values alone to conclude differences or rankings.")
    if str(reason).startswith("stat.jmle."):
        return ("効果を分けて推定できない、または有限の推定値が得られない可能性があります。評定配置と最低点・最高点への偏りを確認してください。"
                if lang == "ja" else "Effects may not be separately estimable, or finite estimates may not exist. Check the rating design and extreme-score patterns.")
    return ("計算の正常終了を確認できていません。結果の解釈に進む前に、分析の終了状況を確認してください。"
            if lang == "ja" else "Successful estimation has not been confirmed. Check the estimation status before interpreting results.")


def fit_convention_note(diagnostics, lang="en"):
    """Describe the plotted primary values, not the optional comparison columns."""
    ja = lang == "ja"
    frame = diagnostics.get("measures")
    if not isinstance(frame, pd.DataFrame) or "FitDfMethod" not in frame:
        frame = diagnostics.get("fit")
    row = frame.iloc[0] if isinstance(frame, pd.DataFrame) and not frame.empty else {}
    method = row.get("FitDfMethod", diagnostics.get("fit_df_method"))
    wm = method in {"facets", "facets_wright_masters"}
    known = wm or method in {"engine", "both", "engine_primary_facets_available", "engine_homogeneous_variance_approximation"}
    basis = ("Wright–Mastersの自由度" if wm else "分散・重みの総和による自由度") if ja else (
        "Wright–Masters degrees of freedom" if wm else "degrees of freedom from variance/weight totals")
    if not known:
        basis = "自由度の方式は未記録" if ja else "degrees-of-freedom method not recorded"
    transform = row.get("FitZSTDTransform", diagnostics.get("fit_zstd_transform"))
    labels = {"Wilson-Hilferty": "Wilson–Hilferty", "linear normal approximation": "線形正規近似" if ja else "linear normal approximation",
              "Wilson-Hilferty without offset": "補正項なしのWilson–Hilferty" if ja else "Wilson–Hilferty without offset"}
    transform = labels.get(transform, "変換方法は未記録" if ja else "transformation not recorded")
    cap = row.get("FitZSTDCap", diagnostics.get("facets_zstd_cap")) if wm else None
    capped = isinstance(cap, (int, float, np.number)) and np.isfinite(cap) and cap > 0
    cap_recorded = known and (not wm or "FitZSTDCap" in row or "facets_zstd_cap" in diagnostics)
    if capped:
        cap_text = f"±{cap:g}で打切り" if ja else f"capped at ±{cap:g}"
    elif cap_recorded:
        cap_text = "打切りなし" if ja else "uncapped"
    else:
        cap_text = "打切りの記録なし" if ja else "cap not recorded"
    return " · ".join((basis, transform, cap_text))


def pathway_axis_range(frame, statistic):
    x = frame[statistic]
    if statistic == "Infit":
        return (0., max(2., float(x.max()) + max(.1, float(x.max() - x.min()) * .12)))
    return (min(-2.8, float(x.min()) - .6), max(2.8, float(x.max()) + .6))


def _font_family():
    # User-provided Japanese labels should survive native drawing and SVG export.
    from matplotlib import font_manager
    available = {font.name for font in font_manager.fontManager.ttflist}
    return next((name for name in ("Hiragino Sans", "Yu Gothic", "Meiryo", "Noto Sans CJK JP") if name in available), "DejaVu Sans")


def style_figure(figure):
    from matplotlib.text import Text
    family = _font_family()
    for artist in figure.findobj(Text):
        artist.set_fontfamily(family)
    figure.set_facecolor("white")
    for axes in figure.axes:
        axes.set_facecolor("#f8fafc")
        axes.spines[["top", "right"]].set_visible(False)
        axes.tick_params(labelsize=9, colors="#334155")
        axes.grid(axis="y", color="#dce3eb", linewidth=.6, zorder=0)
    return figure


def _record_points(figure, artist, rows):
    if not hasattr(figure, "_mfrm_hover_points"):
        figure._mfrm_hover_points = []
    figure._mfrm_hover_points.append((artist, rows.reset_index(drop=True)))


def _extreme_note(result):
    persons = result.get("facets", {}).get("person", pd.DataFrame())
    if "ExtremeScorePattern" not in persons:
        return ""
    count = int(persons["ExtremeScorePattern"].fillna(False).astype(bool).sum())
    return (f" {count} extreme-score persons: displayed locations are approximations; finite JMLEs do not exist."
            if count else "")


def wright_figure(result):
    from matplotlib.figure import Figure

    facets = result.get("facets", {})
    person_table = facets.get("person", pd.DataFrame())
    persons, elements = wright_map_values(person_table, facets.get("others", pd.DataFrame()))
    figure = Figure(figsize=(9, 4.8), layout="constrained")
    left, right = figure.subplots(1, 2, sharey=True, gridspec_kw={"width_ratios": [1, 2.5], "wspace": .08})
    if persons.empty:
        left.text(.5, .5, "No finite person locations", ha="center", transform=left.transAxes)
    else:
        values = persons.to_numpy()
        # Near-identical estimates otherwise produce subpixel bars on the shared logit axis.
        bins = (np.linspace(float(np.median(values)) - .5, float(np.median(values)) + .5, 6)
                if np.ptp(values) < 1e-9 else min(20, max(5, int(np.sqrt(len(persons))))))
        left.hist(values, bins=bins, orientation="horizontal",
                  color="#309b92", edgecolor="white", linewidth=.7)
    left.invert_xaxis()
    left.set(xlabel="Persons (count)", ylabel="Measure (logits)", title=f"Persons · n = {len(persons)}")
    names = list(elements["Facet"].drop_duplicates()) if "Facet" in elements else []
    coordinates = np.r_[persons.to_numpy(), elements.get("Estimate", pd.Series(dtype=float)).to_numpy()]
    span = max(1., float(np.ptp(coordinates))) if len(coordinates) else 1.
    label_limits = list(coordinates)
    for index, name in enumerate(names):
        rows = elements.loc[elements["Facet"].eq(name)]
        artist = right.scatter(np.full(len(rows), index), rows["Estimate"], color=COLORS[index % len(COLORS)],
                               s=65, edgecolors="white", linewidths=.8, zorder=3, label=str(name))
        _record_points(figure, artist, rows)
        # Label text may move; the scientific point coordinates never do.
        if len(rows) <= 10 and len(names) <= 4:
            ordered = rows.sort_values("Estimate")
            yvals = ordered["Estimate"].to_numpy()
            previous = -np.inf
            for (_, row), value in zip(ordered.iterrows(), yvals):
                text_y = max(float(value), previous + span * .055)
                previous = text_y
                label_limits.append(text_y)
                label = str(row.get("Level", ""))
                right.annotate(label if len(label) < 19 else label[:16] + "…", (index, value),
                               xytext=(index + .13, text_y), fontsize=8, color="#334155",
                               arrowprops=dict(arrowstyle="-", color="#94a3b8", lw=.5))
    right.set_xticks(range(len(names)), [str(name) for name in names])
    right.set_xlim(-.4, max(.6, len(names) - .15))
    if label_limits:
        right.set_ylim(min(label_limits) - .08 * span, max(label_limits) + .08 * span)
    right.set(title="Facet locations")
    for axes in (left, right):
        axes.axhline(0, color="#8593a4", linewidth=.8, linestyle="--", zorder=1)
    figure.suptitle("Wright map", fontsize=13, fontweight="bold", color="#17344b")
    figure.supxlabel(textwrap.fill(map_caution(result) + _extreme_note(result), 120), fontsize=8, color="#526274")
    return style_figure(figure)


def pathway_figure(result, diagnostics, *, facet=None, statistic="Infit"):
    from matplotlib.figure import Figure

    frame = pathway_map_values(diagnostics, statistic)
    if facet is not None:
        frame = frame.loc[frame["Facet"].eq(facet)]
    figure = Figure(figsize=(9, 4.8), layout="constrained")
    axes = figure.subplots()
    if statistic == "InfitZSTD":
        for x in (-2, 2):
            axes.axvline(x, linestyle="--", linewidth=.9, color="#8a98a7")
    axes.axvline(1 if statistic == "Infit" else 0, linestyle=":", linewidth=.8, color="#9aa8b5")
    order = list(dict.fromkeys([*result.get("config", {}).get("facet_names", []), "Person", *frame["Facet"].unique()]))
    for index, name in enumerate(order):
        rows = frame.loc[frame["Facet"].eq(name)]
        if rows.empty:
            continue
        artist = axes.scatter(rows[statistic], rows["Estimate"], s=62,
                              color=COLORS[index % len(COLORS)], alpha=.88,
                              edgecolors="white", linewidths=.7, zorder=3, label=str(name))
        _record_points(figure, artist, rows)
        if len(frame) <= 12:
            for number, (_, row) in enumerate(rows.iterrows()):
                label = str(row["Level"])
                axes.annotate(label if len(label) < 23 else label[:20] + "…",
                              (row[statistic], row["Estimate"]),
                              xytext=(7, 7 if number % 2 == 0 else -13), textcoords="offset points", fontsize=8)
    if frame.empty:
        axes.text(.5, .5, "No finite measure / fit pairs", ha="center", transform=axes.transAxes)
    else:
        axes.set_xlim(*pathway_axis_range(frame, statistic))
        axes.margins(y=.22)
        axes.legend(loc="best", frameon=False, fontsize=9)
    axes.set(xlabel="Infit MnSq" if statistic == "Infit" else "Approximate Infit ZSTD", ylabel="Measure (logits)")
    if statistic == "InfitZSTD":
        axes.text(.02, .98, "±2: reference only", va="top", fontsize=9, color="#526274", transform=axes.transAxes)
    title = "Pathway map" + (f" · {facet}" if facet else "")
    figure.suptitle(title, fontsize=13, fontweight="bold", color="#17344b")
    note = fit_convention_note(diagnostics) + "\n" if statistic == "InfitZSTD" else ""
    figure.supxlabel(note + textwrap.fill(map_caution(result) + _extreme_note(result), 120), fontsize=8, color="#526274")
    return style_figure(figure)


def attach_map_hover(canvas, figure, on_text):
    """Show full identifiers on hover, including when direct labels are crowded."""
    def inspect(event):
        for artist, frame in getattr(figure, "_mfrm_hover_points", []):
            if event.inaxes is not artist.axes:
                continue
            hit, info = artist.contains(event)
            if hit and len(info.get("ind", [])):
                row = frame.iloc[int(info["ind"][0])]
                text = f"{row['Facet']} / {row['Level']}  |  測定値 {row['Estimate']:.3f} logit"
                for key, label in (("InfitZSTD", "近似Infit ZSTD"), ("Infit", "Infit MnSq"), ("N_Fit", "評定数")):
                    if key in row and pd.notna(row[key]):
                        text += f"  |  {label} {row[key]:.3g}"
                on_text(text)
                return
        on_text("")
    return canvas.mpl_connect("motion_notify_event", inspect)
