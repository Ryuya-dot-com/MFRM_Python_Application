"""Native figures and offline reports from retained, qualified analysis records."""
import html
import json
import textwrap

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from scipy.stats import norm

from mfrm_app import basic_analysis as basic, desktop_analysis as desktop, engine, output_qualification, project
from mfrm_app import structural_intervals as intervals, visual_maps as maps
from mfrm_app import joint_interaction as joint

KIND_ALIASES = {
    'fit_pathway': 'pathway',
    'ccc': 'category',
    'thresholds': 'step',
    'tcc': 'expected',
    'residual_pca': 'pca',
}

DIAGNOSTIC_GUIDES = {
    'ja': {
        'facet': '各相水準の推定位置（logit）です。条件付き／固定SD MMLの区間バーは探索的スクリーニング用（screening_only, SupportsFormalInference=FALSE）であり、固定アンカー・極端得点・識別不能時は点位置のみを示します。',
        'rater_severity': '相の推定位置（logit）の順位図です。正方向の意味は横軸に示します。区間バーは探索的スクリーニング用（screening_only, SupportsFormalInference=FALSE）であり、正式な差の検定ではありません。',
        'step': '中心化された段階閾値オフセット（logit）のラダーです。閾値の逆転は中間カテゴリーの最頻領域が狭いことを示しますが、自動的にカテゴリー統合を決める規則ではありません。',
        'bubble': '縦軸は推定値（logit）、横軸は適合度、円の面積はスクリーニングSE（保留時は評定数）です。0.5・1.5（またはZSTDの±2）は点検用の参照線であり、固定の合否基準ではありません。',
        'expected': '正味予測子（logit）に対する期待得点曲線 E[X | eta] です。採点基準の最低点・最高点と中央値を破線で示します。',
        'information': '単一評定あたりのカテゴリーFisher情報量 I(eta) と条件付きSE = 1/sqrt(I(eta)) です。較正パラメータを固定した条件付き指標（SupportsFormalInference=FALSE）です。',
        'fair_average': '観測平均得点と、相手相の厳しさ・難しさを調整した公正平均（FairM／FairZ）の比較です。45度破線からの縦方向のずれが文脈調整量を示します。',
        'unexpected': '観測得点と期待得点の乖離が大きい評定（標準化残差 StdResidual）の上位一覧です。±2・±3は点検の目安であり、入力誤りや採点基準の解釈差を原データで確認します。',
        'residual_qq': '対象者ごとに集約した標準化残差（sum(z)/sqrt(n)）の正規Q-Qプロットです。極端得点境界の対象者は除外しています。離散評定のため45度線からの軽微な乖離は生じ得ます。',
        'pca': '標準化残差の主成分分析（PCA）における上位コントラストの固有値です。第1コントラスト固有値が2.0前後以上の場合は、負荷量の正負両極と欠測の重なりを併せて確認してください。',
        'pca_loadings': '標準化残差PCAの第1コントラスト負荷量です。正負両極（±0.4参照線）に位置する相組合せの内容的な共通性を点検します。',
        'bias_ranked': '主効果較正を固定した条件付きバイアス交互作用セルの記述的シフト量（logit）の順位図です。非内部解セルは除外し、検定・区間は保留としています。',
    },
    'en': {
        'facet': 'Ranked facet-level locations (logits). Displayed interval bars are exploratory screening bands (screening_only; SupportsFormalInference=FALSE); fixed, extreme-score and unidentifiable levels show point locations only.',
        'rater_severity': 'Ranked facet locations (logits); the axis shows the direction of increasing estimates. Displayed interval bars are exploratory screening bands (screening_only; SupportsFormalInference=FALSE), not formal significance tests.',
        'step': 'Centered step threshold ladder (logits). Disordered crossings indicate narrow intermediate modal zones, not an automatic category-collapsing rule.',
        'bubble': 'Diagnostic bubble chart of measure (logits) vs fit statistic, with bubble area proportional to screening ModelSE (or observation count when SE is withheld). Reference lines at 0.5/1.5 (or ±2 ZSTD) are inspection prompts, not fixed cutoffs.',
        'expected': 'Expected score curve E[X | eta] across the net logit predictor, bounded by rubric minimum, midpoint and maximum reference lines.',
        'information': 'Single-observation category Fisher information I(eta) and conditional SE = 1/sqrt(I(eta)), holding calibration parameters fixed (screening_only; SupportsFormalInference=FALSE).',
        'fair_average': 'Fair average (FairM / FairZ) versus observed average score across facet levels. Vertical offsets from the 45-degree identity line show context adjustment for partner-facet locations.',
        'unexpected': 'Observations with the largest standardized residuals (StdResidual). Reference lines at ±2 and ±3 guide data and rubric inspection rather than automatic deletion.',
        'residual_qq': 'Normal Q-Q plot of person-aggregated standardized residuals (sum(z)/sqrt(n)), excluding JMLE extreme-score boundary persons. Discrete ratings produce approximate normality only.',
        'pca': 'Standardized residual PCA scree plot. Inspect substantive contrast poles alongside missingness overlap when the first contrast eigenvalue approaches or exceeds 2.0.',
        'pca_loadings': 'First-contrast loadings from standardized residual PCA. Inspect substantive patterns contrasting opposite loading poles (±0.4 guides).',
        'bias_ranked': 'Ranked descriptive conditional eta shifts (logits) for interior bias cells with main calibration held fixed. Non-interior cells are omitted and formal tests/intervals remain withheld.',
    },
}


def request_settings(analysis):
    """The saved calculation request, independent of view choices."""
    if analysis['kind'] == 'gd-analysis-v2':
        return dict(workflow='gd', gd_settings=analysis['gd_analysis']['settings'])
    request = dict(workflow='mfrm', **analysis['fit_settings'],
        bias_requests=[e['settings'] for e in analysis.get('bias_analyses', [])],
        interval_requests=[e['settings'] for e in analysis.get('interval_analyses', [])],
        gd_settings=analysis.get('gd_analysis', {}).get('settings'))
    if 'joint_analysis' in analysis:
        request['joint_settings'] = analysis['joint_analysis']['settings']
    return request


current_tables = desktop.current_tables


def figure_choices(analysis, *, language=None):
    """Each choice has explicit model/level/pair provenance, including in reports."""
    choices = []
    if 'result' in analysis:
        config = analysis['result']['config']
        levels = config['facet_levels'][config['step_facet']] if config['model'] == 'PCM' else [None]
        for kind, title in [('wright', 'Wright map'), ('category', 'カテゴリー確率')]:
            for level in levels:
                suffix = f" · {config['step_facet']} = {level}" if level is not None else ''
                choices.append((title + suffix, dict(kind=kind, level=level)))
        choices.insert(1, ('Pathway map', dict(kind='pathway')))
        if 'Rater' in config['facet_names']:
            label = '評価者の寛大さ' if config.get('facet_signs', {}).get('Rater', -1) > 0 else '評価者の厳しさ'
            choices.append((label, dict(kind='rater_severity')))
        choices.extend([
            ('相水準の推定値', dict(kind='facet')),
            ('段階閾値ラダー', dict(kind='step')),
            ('適合度バブルチャート', dict(kind='bubble')),
            ('期待得点曲線', dict(kind='expected')),
            ('情報量・条件付きSE曲線', dict(kind='information')),
            ('公正平均 vs 観測平均', dict(kind='fair_average')),
            ('予期しない応答', dict(kind='unexpected')),
            ('標準化残差 Q-Q', dict(kind='residual_qq')),
        ])
    for i, entry in enumerate(analysis.get('bias_analyses', [])):
        s = entry['settings']
        choices.append((f"バイアス · {s['facet_a']} × {s['facet_b']}", dict(kind='bias', index=i)))
        choices.append((f"バイアス順位 · {s['facet_a']} × {s['facet_b']}", dict(kind='bias_ranked', index=i)))
    for i, entry in enumerate(analysis.get('interval_analyses', [])):
        s = entry['settings']
        suffix = ' / ' + ' − '.join(s['difference']) if s['difference'] else ''
        choices.append((f"指定区間 {i+1} · {s['facet']}{suffix}", dict(kind='interval', index=i)))
    if 'gd_analysis' in analysis:
        choices.append(('D-Study · G / Phi', dict(kind='gd')))
    if 'joint_analysis' in analysis:
        s = analysis['joint_analysis']['settings']
        choices.append((f"共同交互作用 · {s['facet_a']} × {s['facet_b']}", dict(kind='joint')))
    if analysis.get('diagnostics', {}).get('pca_enabled'):
        choices.append(('残差PCA · 固有値スクリー', dict(kind='pca')))
        choices.append(('残差PCA · 第1コントラスト負荷量', dict(kind='pca_loadings')))
    if language is not None:
        from mfrm_app.desktop_i18n import text
        return [(text(title, language), opts) for title, opts in choices]
    return choices


def _steps(result, level):
    config = result['config']
    if config['model'] == 'RSM':
        return np.asarray(result['params']['steps']), 'RSM · common steps'
    levels = list(config['facet_levels'][config['step_facet']])
    if level not in levels:
        raise ValueError('Choose a retained PCM step-facet level')
    return np.asarray(result['params']['steps_mat'])[levels.index(level)], f"PCM · {config['step_facet']} = {level}"


def _empty(axes, reason):
    axes.text(.5, .5, textwrap.fill(str(reason), 70), ha='center', va='center', transform=axes.transAxes)
    axes.set_axis_off()


def _step_series_list(result, level):
    config = result['config']
    if config['model'] == 'RSM' or level is not None:
        steps, label = _steps(result, level)
        if not np.isfinite(steps).all():
            raise ValueError('Category steps must be finite')
        return [(label, steps)]
    levels = list(config['facet_levels'][config['step_facet']])
    mat = np.asarray(result['params']['steps_mat'], dtype=float)
    out = [(f"{config['step_facet']} = {lv}", mat[i]) for i, lv in enumerate(levels) if np.isfinite(mat[i]).all()]
    if not out:
        raise ValueError('Category steps must be finite')
    return out


def _require_mfrm(analysis, kind):
    if 'result' not in analysis:
        raise ValueError(f"Figure kind '{kind}' requires a fitted MFRM result")
    return analysis['result'], analysis.get('diagnostics') or {}


def _validate_facet_choice(facet, result, extra_facets=()):
    if facet is None:
        return
    known = {'Person', *(result.get('config', {}).get('facet_names') or []), *map(str, extra_facets)}
    if str(facet) not in known:
        raise ValueError(f'Unknown facet: {facet}')


def _facet_estimates_figure(analysis, *, facet=None, rater_default=False, top_n=40):
    result, diagnostics = _require_mfrm(analysis, 'rater_severity' if rater_default else 'facet')
    fnames = list(result.get('config', {}).get('facet_names') or [])
    if rater_default and facet is None:
        facet = 'Rater' if 'Rater' in fnames else (fnames[0] if fnames else None)
    raw = diagnostics.get('measures')
    extra = raw['Facet'].astype(str).unique() if isinstance(raw, pd.DataFrame) and 'Facet' in raw.columns else ()
    _validate_facet_choice(facet, result, extra)
    top_n = max(1, int(top_n))
    title_base = 'Facet-level estimates'
    xlabel = 'Measure (logits)'
    if facet is not None:
        positive = result['config'].get('facet_signs', {}).get(facet, -1) > 0 or facet == 'Person'
        xlabel += ' · higher values ' + ('raise' if positive else 'lower') + ' expected scores'
        if rater_default and facet.casefold() == 'rater':
            title_base = 'Rater leniency profile' if positive else 'Rater severity profile'
        title_base += f' · {facet}'

    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    hold_reason = output_qualification.calibration_inference_hold_reason(result)
    if hold_reason == output_qualification.OPTIMIZER_TERMINATION_HOLD_REASON:
        _empty(ax, 'Withheld: ' + maps.map_caution(result))
        fig.suptitle(title_base)
        fig.supxlabel(textwrap.fill(maps.map_caution(result), 110), fontsize=8)
        return maps.style_figure(fig)

    q = engine.qualify_measures_for_reporting(
        raw if isinstance(raw, pd.DataFrame) else pd.DataFrame(),
        result,
        obs_df=diagnostics.get('obs'),
    )
    if q.empty or not {'Facet', 'Level', 'Estimate'}.issubset(q.columns):
        _empty(ax, 'No facet-level estimates available')
        fig.suptitle(title_base)
        return maps.style_figure(fig)

    q['Estimate'] = pd.to_numeric(q['Estimate'], errors='coerce')
    rows = q.loc[q['Facet'].astype(str).eq(str(facet)) if facet is not None else ~q['Facet'].astype(str).eq('Person')].copy()
    if facet == 'Person':
        for col in ('FiniteJMLEEstimate', 'PersonInferenceReady'):
            if col in rows.columns:
                rows = rows.loc[~rows[col].eq(False).fillna(False)]
    if 'UnboundedFacetEstimate' in rows.columns:
        rows = rows.loc[~rows['UnboundedFacetEstimate'].eq(True).fillna(False)]
    rows = rows.loc[np.isfinite(rows['Estimate'])].copy()
    if rows.empty:
        _empty(ax, 'No finite estimates available for the selected facet')
        fig.suptitle(title_base)
        return maps.style_figure(fig)

    truncated = False
    if len(rows) > top_n:
        rows = rows.assign(_abs=rows['Estimate'].abs()).sort_values('_abs', ascending=False).head(top_n).drop(columns=['_abs'])
        truncated = True
    facet_order = {name: idx for idx, name in enumerate([*fnames, 'Person'])}
    rows = rows.assign(_f_ord=rows['Facet'].astype(str).map(lambda f: facet_order.get(f, 999)))
    rows = rows.sort_values(['_f_ord', 'Estimate'], kind='mergesort').drop(columns=['_f_ord']).reset_index(drop=True)

    ax.axvline(0, color='#8593a4', linestyle='--', linewidth=.8, zorder=1)
    multi_facet = rows['Facet'].astype(str).nunique() > 1
    ci_ok = (
        rows.get('CIUse', pd.Series('', index=rows.index)).astype(str).eq('screening_only')
        & np.isfinite(pd.to_numeric(rows.get('EffectiveCI_Lower', np.nan), errors='coerce'))
        & np.isfinite(pd.to_numeric(rows.get('EffectiveCI_Upper', np.nan), errors='coerce'))
    )
    rows['_y'] = np.arange(len(rows))
    for idx, fname in enumerate(dict.fromkeys(rows['Facet'].astype(str))):
        sub = rows.loc[rows['Facet'].astype(str).eq(fname)]
        color = maps.COLORS[idx % len(maps.COLORS)]
        sub_ci = sub.loc[ci_ok.loc[sub.index]]
        sub_pt = sub.loc[~ci_ok.loc[sub.index]]
        if not sub_ci.empty:
            x = sub_ci['Estimate'].to_numpy(dtype=float)
            lo = pd.to_numeric(sub_ci['EffectiveCI_Lower'], errors='coerce').to_numpy(dtype=float)
            hi = pd.to_numeric(sub_ci['EffectiveCI_Upper'], errors='coerce').to_numpy(dtype=float)
            ax.errorbar(
                x, sub_ci['_y'].to_numpy(), xerr=np.vstack([np.maximum(0., x - lo), np.maximum(0., hi - x)]),
                fmt='o', color=color, capsize=3, markersize=5.5, linewidth=1.1, label=fname if multi_facet else 'Screening interval',
            )
        if not sub_pt.empty:
            ax.scatter(
                sub_pt['Estimate'].to_numpy(dtype=float), sub_pt['_y'].to_numpy(),
                marker='D', s=42, color=color, edgecolors='#1e293b', linewidths=.7, zorder=4,
                label=(f'{fname} (point only)' if multi_facet else 'Point only (fixed / withheld SE)'),
            )
        artist = ax.scatter(sub['Estimate'].to_numpy(dtype=float), sub['_y'].to_numpy(), s=1, alpha=0)
        maps._record_points(fig, artist, sub)

    ylabels = [
        (f"{r.Facet}: {r.Level}" if multi_facet else str(r.Level))[:28]
        for r in rows.itertuples(index=False)
    ]
    ax.set_yticks(rows['_y'].to_numpy(), ylabels)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    if multi_facet or (ci_ok.any() and (~ci_ok).any()):
        ax.legend(fontsize=8, loc='best', frameon=False)
    fig.suptitle(title_base)
    trunc_note = f' Showing top {top_n} levels by |Estimate|.' if truncated else ''
    caution = maps.map_caution(result)
    band_note = (
        'Exploratory screening bands (screening_only; SupportsFormalInference=FALSE).'
        if ci_ok.any()
        else 'Uncertainty intervals withheld; showing point estimates only.'
    )
    fig.supxlabel(textwrap.fill((band_note + trunc_note + (' ' + caution if caution else '')).strip(), 115), fontsize=8)
    return maps.style_figure(fig)


def _step_ladder_figure(analysis, *, level=None):
    result, _ = _require_mfrm(analysis, 'step')
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    if output_qualification.calibration_inference_hold_reason(result) in {
        output_qualification.OPTIMIZER_TERMINATION_HOLD_REASON,
        output_qualification.CATEGORY_UNSUPPORTED_HOLD_REASON,
        output_qualification.UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON,
        output_qualification.UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON,
    }:
        _empty(ax, 'Withheld: ' + maps.map_caution(result))
        fig.suptitle('Step threshold ladder')
        fig.supxlabel(textwrap.fill(maps.map_caution(result), 110), fontsize=8)
        return maps.style_figure(fig)

    series_list = _step_series_list(result, level)
    ax.axhline(0, color='#8593a4', linestyle='--', linewidth=.8, zorder=1)
    rmin = int(result['prep']['rating_min'])
    disordered_groups = 0
    max_k = 0
    for idx, (lbl, steps) in enumerate(series_list):
        steps = np.asarray(steps, dtype=float)
        xs = np.arange(1, len(steps) + 1)
        max_k = max(max_k, len(steps))
        color = maps.COLORS[idx % len(maps.COLORS)]
        ax.plot(xs, steps, marker='o', linewidth=1.5, markersize=6, color=color, label=lbl)
        diffs = np.diff(steps)
        if np.any(diffs < -1e-12):
            disordered_groups += 1
            for k_idx in np.where(diffs < -1e-12)[0]:
                ax.plot(
                    [xs[k_idx], xs[k_idx + 1]], [steps[k_idx], steps[k_idx + 1]],
                    color='#dc2626', linestyle='--', linewidth=2.0, marker='D', markersize=6.5,
                    label='Disordered step' if disordered_groups == 1 and k_idx == np.where(diffs < -1e-12)[0][0] else None,
                )
    ax.set_xticks(range(1, max_k + 1), [f'{rmin + k - 1}->{rmin + k}' for k in range(1, max_k + 1)])
    ax.set(xlabel='Adjacent category transition', ylabel='Centered step threshold (logits)')
    if len(series_list) > 1 or disordered_groups:
        ax.legend(fontsize=8, ncols=min(3, len(series_list) + 1), frameon=False)
    disorder_txt = f' · {disordered_groups} disordered ladder(s)' if disordered_groups else ''
    fig.suptitle('Step threshold ladder' + disorder_txt)
    caution = maps.map_caution(result)
    note = 'Centered step offsets (logits). Step disorder indicates narrow intermediate modal zones, not an automatic category-merging rule.'
    fig.supxlabel(textwrap.fill((note + (' ' + caution if caution else '')).strip(), 115), fontsize=8)
    return maps.style_figure(fig)


def _bubble_figure(analysis, *, facet=None, statistic='Infit', top_n=60):
    result, diagnostics = _require_mfrm(analysis, 'bubble')
    if statistic not in {'Infit', 'Outfit', 'InfitZSTD', 'OutfitZSTD'}:
        raise ValueError('Choose Infit, Outfit, InfitZSTD, or OutfitZSTD')
    raw = diagnostics.get('measures')
    extra = raw['Facet'].astype(str).unique() if isinstance(raw, pd.DataFrame) and 'Facet' in raw.columns else ()
    _validate_facet_choice(facet, result, extra)
    top_n = max(1, int(top_n))
    q = engine.qualify_measures_for_reporting(
        raw if isinstance(raw, pd.DataFrame) else pd.DataFrame(),
        result,
        obs_df=diagnostics.get('obs'),
    )
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    is_zstd = statistic in {'InfitZSTD', 'OutfitZSTD'}
    rows = q.loc[q['Facet'].astype(str).eq(str(facet)) if facet is not None else ~q['Facet'].astype(str).eq('Person')].copy() if not q.empty else pd.DataFrame()
    if facet == 'Person' and not rows.empty:
        for col in ('FiniteJMLEEstimate', 'PersonInferenceReady'):
            if col in rows.columns:
                rows = rows.loc[~rows[col].eq(False).fillna(False)]
    if not rows.empty and 'UnboundedFacetEstimate' in rows.columns:
        rows = rows.loc[~rows['UnboundedFacetEstimate'].eq(True).fillna(False)]
    guard_detail = ''
    if not rows.empty:
        f_statuses = set(rows.get('FitStatus', pd.Series(dtype=str)).astype(str))
        if is_zstd and 'zstd_unavailable_boundary_proxy' in f_statuses:
            guard_detail += ' ZSTD omitted where JMLE extreme-score boundary proxies are included.'
        if is_zstd and 'zstd_unavailable_low_df' in f_statuses:
            guard_detail += ' ZSTD omitted where degrees of freedom < 1.'
        if 'unavailable' in f_statuses:
            guard_detail += ' Fit omitted where variance/residuals/weights are unavailable.'
        rows['Estimate'] = pd.to_numeric(rows['Estimate'], errors='coerce')
        rows[statistic] = pd.to_numeric(rows[statistic], errors='coerce')
        rows = rows.loc[np.isfinite(rows['Estimate']) & np.isfinite(rows[statistic])].copy()

    if rows.empty:
        _empty(ax, 'No finite measure / fit pairs' + (f' ({guard_detail.strip()})' if guard_detail else ''))
        fig.suptitle(f'Diagnostic bubble chart · {statistic}')
        fig.supxlabel(textwrap.fill((maps.map_caution(result) + guard_detail).strip() or 'Bubble chart unavailable.', 115), fontsize=8)
        return maps.style_figure(fig)

    ref_x = 0.0 if is_zstd else 1.0
    if len(rows) > top_n:
        rows = rows.assign(_dev=(rows[statistic] - ref_x).abs()).sort_values('_dev', ascending=False).head(top_n).drop(columns=['_dev'])

    se_vals = pd.to_numeric(rows.get('ModelSE', np.nan), errors='coerce')
    se_ok = rows.get('SEUse', pd.Series('', index=rows.index)).astype(str).eq('screening_only') & np.isfinite(se_vals) & (se_vals > 0)
    if se_ok.any():
        raw_s = se_vals.where(se_ok, float(se_vals.loc[se_ok].min())).to_numpy(dtype=float)
        size_note = 'Bubble area proportional to screening ModelSE (SupportsFormalInference=FALSE).'
    else:
        n_vals = pd.to_numeric(rows.get('N_Fit', 1.0), errors='coerce').fillna(1.0).clip(lower=1.0)
        raw_s = n_vals.to_numpy(dtype=float)
        size_note = 'Bubble area proportional to observation count N_Fit (SE withheld).'
    s_ptp = float(np.ptp(raw_s))
    areas = np.full(len(rows), 130.0) if s_ptp < 1e-12 else 65.0 + 235.0 * (raw_s - float(raw_s.min())) / s_ptp
    rows = rows.assign(_area=areas)

    ax.axhline(0, color='#8593a4', linestyle='--', linewidth=.75, zorder=1)
    ax.axvline(ref_x, color='#9aa8b5', linestyle=':', linewidth=.85, zorder=1)
    for gx in ((-2.0, 2.0) if is_zstd else (0.5, 1.5)):
        ax.axvline(gx, color='#8a98a7', linestyle='--', linewidth=.85, zorder=1)

    fnames = list(dict.fromkeys([*(result.get('config', {}).get('facet_names') or []), *rows['Facet'].astype(str)]))
    for idx, fname in enumerate(fnames):
        sub = rows.loc[rows['Facet'].astype(str).eq(fname)]
        if sub.empty:
            continue
        artist = ax.scatter(
            sub[statistic], sub['Estimate'], s=sub['_area'],
            color=maps.COLORS[idx % len(maps.COLORS)], alpha=.68,
            edgecolors='#1e293b', linewidths=.8, zorder=3, label=fname,
        )
        maps._record_points(fig, artist, sub)
        if len(rows) <= 16:
            for _, r in sub.iterrows():
                ax.annotate(str(r['Level'])[:18], (r[statistic], r['Estimate']), xytext=(5, 4), textcoords='offset points', fontsize=8)
    stat_base = statistic[:-4] if is_zstd else statistic
    ax.set(xlabel=f'Approximate {stat_base} ZSTD' if is_zstd else f'{stat_base} MnSq', ylabel='Measure (logits)')
    ax.set_xlim(*maps.pathway_axis_range(rows, statistic))
    ax.margins(y=.18)
    ax.legend(loc='best', frameon=False, fontsize=8)
    fig.suptitle(f'Diagnostic bubble chart · {statistic}' + (f' · {facet}' if facet else ''))
    caution = maps.map_caution(result)
    fig.supxlabel(textwrap.fill(' '.join(p for p in (size_note, guard_detail.strip(), caution) if p), 115), fontsize=8)
    return maps.style_figure(fig)


def _expected_score_figure(analysis, *, level=None):
    result, _ = _require_mfrm(analysis, 'expected')
    if output_qualification.calibration_inference_hold_reason(result) in {
        output_qualification.OPTIMIZER_TERMINATION_HOLD_REASON,
        output_qualification.CATEGORY_UNSUPPORTED_HOLD_REASON,
        output_qualification.UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON,
        output_qualification.UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON,
    }:
        fig = Figure(figsize=(9, 4.8), layout='constrained')
        ax = fig.subplots()
        _empty(ax, 'Withheld: ' + maps.map_caution(result))
        fig.suptitle('Expected score curve')
        fig.supxlabel(textwrap.fill(maps.map_caution(result), 110), fontsize=8)
        return maps.style_figure(fig)
    series_list = _step_series_list(result, level)
    rmin = int(result['prep']['rating_min'])
    rmax = int(result['prep']['rating_max'])
    all_steps = np.concatenate([s for _, s in series_list])
    eta = np.linspace(min(-4., float(all_steps.min()) - 2.), max(4., float(all_steps.max()) + 2.), 201)
    scores = np.arange(rmax - rmin + 1, dtype=float)

    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    for bound in (rmin, rmax):
        ax.axhline(bound, color='#94a3b8', linestyle='--', linewidth=.75)
    ax.axhline(0.5 * (rmin + rmax), color='#cbd5e1', linestyle=':', linewidth=.75)
    for idx, (lbl, steps) in enumerate(series_list):
        probs = engine.category_prob_rsm(eta, np.r_[0., np.cumsum(steps)])
        expected = rmin + probs @ scores
        ax.plot(eta, expected, linewidth=1.8, color=maps.COLORS[idx % len(maps.COLORS)], label=lbl)
    ax.set(xlabel='Net predictor (logits)', ylabel='Expected score', ylim=(rmin - 0.15, rmax + 0.15))
    if len(series_list) > 1:
        ax.legend(fontsize=8, ncols=min(3, len(series_list)), frameon=False)
    title_lbl = series_list[0][0] if len(series_list) == 1 else f"PCM · {result['config']['step_facet']} levels"
    fig.suptitle(f'{title_lbl} · expected score curve')
    fig.supxlabel(textwrap.fill(maps.map_caution(result) or 'Model-implied expected score E[X | eta] on the original rating scale.', 115), fontsize=8)
    return maps.style_figure(fig)


def _information_figure(analysis, *, level=None):
    result, _ = _require_mfrm(analysis, 'information')
    if output_qualification.calibration_inference_hold_reason(result) in {
        output_qualification.OPTIMIZER_TERMINATION_HOLD_REASON,
        output_qualification.CATEGORY_UNSUPPORTED_HOLD_REASON,
        output_qualification.UNOBSERVED_BOUNDARY_CATEGORY_HOLD_REASON,
        output_qualification.UNOBSERVED_OWNER_ENDPOINT_HOLD_REASON,
    }:
        fig = Figure(figsize=(9.2, 4.6), layout='constrained')
        ax = fig.subplots()
        _empty(ax, 'Withheld: ' + maps.map_caution(result))
        fig.suptitle('Information and conditional SE')
        fig.supxlabel(textwrap.fill(maps.map_caution(result), 110), fontsize=8)
        return maps.style_figure(fig)
    series_list = _step_series_list(result, level)
    all_steps = np.concatenate([s for _, s in series_list])
    eta = np.linspace(min(-4., float(all_steps.min()) - 2.), max(4., float(all_steps.max()) + 2.), 201)

    fig = Figure(figsize=(9.2, 4.6), layout='constrained')
    left, right = fig.subplots(1, 2)
    sem_interior_max = 1.5
    for idx, (lbl, steps) in enumerate(series_list):
        probs = engine.category_prob_rsm(eta, np.r_[0., np.cumsum(steps)])
        ks = np.arange(probs.shape[1], dtype=float)
        mean_k = probs @ ks
        info = np.maximum((probs * (ks - mean_k[:, None]) ** 2).sum(axis=1), 1e-12)
        sem = 1.0 / np.sqrt(info)
        sem_interior_max = max(sem_interior_max, float(np.percentile(sem, 90)) * 1.35)
        color = maps.COLORS[idx % len(maps.COLORS)]
        left.plot(eta, info, linewidth=1.7, color=color, label=lbl)
        right.plot(eta, sem, linewidth=1.7, color=color, label=lbl)
    left.set(xlabel='Net predictor (logits)', ylabel='Fisher information I(eta)', title='Category information', ylim=(0, None))
    right.set(xlabel='Net predictor (logits)', ylabel='Conditional SE (logits)', title='Conditional SE = 1 / sqrt(I)', ylim=(0, min(8.0, sem_interior_max)))
    if len(series_list) > 1:
        left.legend(fontsize=8, frameon=False)
    title_lbl = series_list[0][0] if len(series_list) == 1 else f"PCM · {result['config']['step_facet']} levels"
    fig.suptitle(f'{title_lbl} · information and conditional SE')
    caution = maps.map_caution(result)
    note = 'Single-observation Fisher information and conditional SE holding calibration parameters fixed (screening_only; SupportsFormalInference=FALSE).'
    fig.supxlabel(textwrap.fill((note + (' ' + caution if caution else '')).strip(), 115), fontsize=8)
    return maps.style_figure(fig)


def _fair_average_figure(analysis, *, facet=None, statistic='FairM', top_n=60):
    result, _ = _require_mfrm(analysis, 'fair_average')
    report_tables = analysis.get('report_tables') or {}
    _validate_facet_choice(facet, result, report_tables.keys())
    metric = 'FairZ' if statistic == 'FairZ' else 'FairM'
    target_facets = [str(facet)] if facet is not None else [f for f in result['config']['facet_names'] if f in report_tables]
    frames = []
    extreme_persons = engine._boundary_person_set(result) if facet == 'Person' else set()
    b_audit = output_qualification.audit_boundary_and_support(result)
    unbounded_levels = set(b_audit['unbounded_facet_levels']) | set(b_audit['unbounded_step_owners'])
    for fname in target_facets:
        tbl = report_tables.get(fname)
        if isinstance(tbl, pd.DataFrame) and not tbl.empty and 'ObservedAverage' in tbl.columns:
            sub = tbl.copy()
            sub['Facet'] = fname
            if fname == 'Person' and extreme_persons and 'Level' in sub.columns:
                sub = sub.loc[~sub['Level'].astype(str).isin(extreme_persons)]
            if unbounded_levels and 'Level' in sub.columns:
                sub = sub.loc[~sub['Level'].astype(str).map(lambda lv, f=fname: (f, lv) in unbounded_levels)]
            frames.append(sub)
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    if not frames:
        _empty(ax, 'Fair-average table is unavailable')
        fig.suptitle('Fair vs observed average')
        return maps.style_figure(fig)
    df = pd.concat(frames, ignore_index=True)
    df['ObservedAverage'] = pd.to_numeric(df['ObservedAverage'], errors='coerce')
    fallback_note = ''
    if metric not in df.columns or not np.isfinite(pd.to_numeric(df[metric], errors='coerce')).any():
        if metric == 'FairM' and 'FairZ' in df.columns and np.isfinite(pd.to_numeric(df['FairZ'], errors='coerce')).any():
            metric = 'FairZ'
            fallback_note = 'FairM unavailable (unbounded person mean); displaying zero-reference FairZ.'
    df[metric] = pd.to_numeric(df.get(metric, np.nan), errors='coerce')
    df = df.loc[np.isfinite(df['ObservedAverage']) & np.isfinite(df[metric])].copy()
    if df.empty:
        _empty(ax, 'No finite observed/fair average pairs available')
        fig.suptitle(f'Fair ({metric}) vs observed average')
        return maps.style_figure(fig)
    if len(df) > max(1, int(top_n)):
        df = df.assign(_gap=(df[metric] - df['ObservedAverage']).abs()).sort_values('_gap', ascending=False).head(int(top_n)).drop(columns=['_gap'])

    lo = min(float(df['ObservedAverage'].min()), float(df[metric].min()))
    hi = max(float(df['ObservedAverage'].max()), float(df[metric].max()))
    pad = max(0.1, (hi - lo) * 0.08)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color='#8593a4', linestyle='--', linewidth=.9, label='Identity (Fair = Observed)')
    for _, r in df.iterrows():
        ax.plot([r['ObservedAverage'], r['ObservedAverage']], [r['ObservedAverage'], r[metric]], color='#cbd5e1', linewidth=.8, zorder=2)
    for idx, fname in enumerate(dict.fromkeys(df['Facet'].astype(str))):
        sub = df.loc[df['Facet'].astype(str).eq(fname)]
        artist = ax.scatter(
            sub['ObservedAverage'], sub[metric], s=58,
            color=maps.COLORS[idx % len(maps.COLORS)], edgecolors='white', linewidths=.7, zorder=3, label=fname,
        )
        maps._record_points(fig, artist, sub)
        if len(df) <= 16:
            for _, r in sub.iterrows():
                ax.annotate(str(r['Level'])[:18], (r['ObservedAverage'], r[metric]), xytext=(5, 4), textcoords='offset points', fontsize=8)
    ax.set(xlim=(lo - pad, hi + pad), ylim=(lo - pad, hi + pad), xlabel='Observed average score', ylabel=f'Fair average ({metric})')
    ax.legend(fontsize=8, loc='best', frameon=False)
    fig.suptitle(f'Fair ({metric}) vs observed average' + (f' · {facet}' if facet else ''))
    caution = maps.map_caution(result)
    note = 'Vertical segments show score-scale adjustment for partner-facet context.'
    fig.supxlabel(textwrap.fill(' '.join(p for p in (note, fallback_note, caution) if p), 115), fontsize=8)
    return maps.style_figure(fig)


def _valid_obs_residuals(analysis, *, exclude_extreme_persons=False):
    result, diagnostics = _require_mfrm(analysis, 'unexpected')
    obs = diagnostics.get('obs')
    if not isinstance(obs, pd.DataFrame) or obs.empty:
        obs = analysis.get('scorefile')
    if not isinstance(obs, pd.DataFrame) or obs.empty or 'StdResidual' not in obs.columns:
        return result, pd.DataFrame(), 0
    df = obs.copy()
    total_n = len(df)
    std_r = pd.to_numeric(df['StdResidual'], errors='coerce')
    var_ok = pd.to_numeric(df.get('Var', 1.0), errors='coerce') > 1e-10
    proxy = df.get('FitBoundaryProxy', pd.Series(False, index=df.index)).fillna(False).astype(bool)
    mask = np.isfinite(std_r) & var_ok & (~proxy)
    if exclude_extreme_persons and str(result.get('config', {}).get('method', '')).upper() == 'JMLE' and 'Person' in df.columns:
        ext_p = engine._boundary_person_set(result)
        if ext_p:
            mask = mask & (~df['Person'].astype(str).isin(ext_p))
    df = df.loc[mask].assign(StdResidual=std_r.loc[mask])
    return result, df, total_n - len(df)


def _unexpected_figure(analysis, *, top_n=25):
    result, df, n_excluded = _valid_obs_residuals(analysis, exclude_extreme_persons=False)
    fig = Figure(figsize=(9, 5.0), layout='constrained')
    ax = fig.subplots()
    if df.empty:
        _empty(ax, 'Standardized residuals unavailable')
        fig.suptitle('Unexpected responses')
        return maps.style_figure(fig)
    abs_z = df['StdResidual'].abs()
    n_ge2 = int((abs_z >= 2.0).sum())
    n_ge3 = int((abs_z >= 3.0).sum())
    cap = max(1, min(int(top_n), 25))
    flagged = df.loc[abs_z >= 2.0].copy()
    if flagged.empty:
        shown = df.assign(_az=abs_z).sort_values('_az', ascending=False).head(min(cap, 15)).drop(columns=['_az'])
        rule_note = f'0 observations with |StdResidual| >= 2.0; showing {len(shown)} largest absolute residuals.'
    else:
        shown = flagged.assign(_az=flagged['StdResidual'].abs()).sort_values('_az', ascending=False).head(cap).drop(columns=['_az'])
        rule_note = f'Flagged |StdResidual| >= 2.0: {n_ge2} (|StdResidual| >= 3.0: {n_ge3}); showing top {len(shown)}.'
    shown = shown.sort_values('StdResidual', kind='mergesort').reset_index(drop=True)
    id_cols = [c for c in ['Person', *(result.get('config', {}).get('facet_names') or [])] if c in shown.columns]
    labels = []
    for _, r in shown.iterrows():
        who = '·'.join(str(r[c]) for c in id_cols)
        obs_v = r.get('Observed', r.get('Score', np.nan))
        exp_v = r.get('Expected', np.nan)
        suffix = f' ({obs_v:g} vs {exp_v:.2f})' if np.isfinite(obs_v) and np.isfinite(exp_v) else ''
        labels.append((who[:22] + ('…' if len(who) > 22 else '')) + suffix)
    ys = np.arange(len(shown))
    zs = shown['StdResidual'].to_numpy(dtype=float)
    ax.axvline(0, color='#64748b', linewidth=.8)
    for gx, ls in ((-2, '--'), (2, '--'), (-3, ':'), (3, ':')):
        ax.axvline(gx, color='#8a98a7', linestyle=ls, linewidth=.85)
    for y_i, z_i in zip(ys, zs):
        ax.plot([0, z_i], [y_i, y_i], color='#cbd5e1', linewidth=1.1, zorder=2)
    pos = zs >= 0
    if np.any(pos):
        ax.scatter(zs[pos], ys[pos], color='#d88627', s=52, zorder=3, label='Higher than expected (z > 0)')
    if np.any(~pos):
        ax.scatter(zs[~pos], ys[~pos], color='#2176ae', s=52, zorder=3, label='Lower than expected (z < 0)')
    ax.set_yticks(ys, labels, fontsize=8)
    ax.set_xlabel('Standardized residual (StdResidual)')
    lim = max(3.3, float(np.max(np.abs(zs))) * 1.12)
    ax.set_xlim(-lim, lim)
    ax.legend(fontsize=8, loc='best', frameon=False)
    fig.suptitle('Unexpected responses · standardized residuals')
    excl_note = f' Excluded {n_excluded} zero-variance/boundary row(s).' if n_excluded else ''
    fig.supxlabel(textwrap.fill((rule_note + excl_note + ' ' + maps.map_caution(result)).strip(), 115), fontsize=8)
    return maps.style_figure(fig)


def _residual_qq_figure(analysis):
    result, df, n_excluded = _valid_obs_residuals(analysis, exclude_extreme_persons=True)
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    if df.empty or 'Person' not in df.columns:
        _empty(ax, 'Person standardized residuals unavailable for Q-Q plot')
        fig.suptitle('Normal Q-Q of person standardized residuals')
        return maps.style_figure(fig)
    grouped = df.groupby(df['Person'].astype(str), sort=False)['StdResidual'].agg(
        lambda s: float(s.sum() / np.sqrt(max(1, len(s))))
    )
    grouped = grouped[np.isfinite(grouped)]
    if len(grouped) < 2:
        _empty(ax, 'At least 2 interior persons are required for a residual Q-Q plot')
        fig.suptitle('Normal Q-Q of person standardized residuals')
        return maps.style_figure(fig)
    ordered = grouped.sort_values(kind='mergesort')
    n = len(ordered)
    a = 0.375 if n <= 10 else 0.5
    theo = norm.ppf((np.arange(1, n + 1) - a) / (n + 1.0 - 2.0 * a))
    samp = ordered.to_numpy(dtype=float)
    lo = min(float(theo.min()), float(samp.min()), -2.5)
    hi = max(float(theo.max()), float(samp.max()), 2.5)
    ax.plot([lo, hi], [lo, hi], color='#8593a4', linestyle='--', linewidth=.9, label='Identity N(0, 1)')
    for gy in (-2.0, 2.0):
        ax.axhline(gy, color='#cbd5e1', linestyle=':', linewidth=.8)
    q_rows = pd.DataFrame({'Person': ordered.index, 'TheoreticalQuantile': theo, 'SampleQuantile': samp})
    artist = ax.scatter(theo, samp, color='#2176ae', s=46, alpha=.85, edgecolors='white', linewidths=.6, zorder=3, label=f'Persons (n = {n})')
    maps._record_points(fig, artist, q_rows)
    ax.set(xlabel='Theoretical N(0, 1) quantile', ylabel='Person residual sum(z) / sqrt(n)', xlim=(lo - .15, hi + .15), ylim=(lo - .15, hi + .15))
    ax.legend(fontsize=8, loc='upper left', frameon=False)
    fig.suptitle('Normal Q-Q of person standardized residuals')
    excl_txt = f' Excluded {n_excluded} boundary/zero-variance observation(s).' if n_excluded else ''
    note = f'Aggregated across ratings per interior person (n = {n}). Discrete polytomous ratings yield approximate normality.{excl_txt}'
    fig.supxlabel(textwrap.fill((note + ' ' + maps.map_caution(result)).strip(), 115), fontsize=8)
    return maps.style_figure(fig)


def _pca_loadings_figure(analysis, *, facet=None, top_n=30):
    result, diagnostics = _require_mfrm(analysis, 'pca_loadings')
    _validate_facet_choice(facet, result)
    if facet is not None:
        pca = (diagnostics.get('pca_by_facet') or {}).get(str(facet))
        reason = (diagnostics.get('pca_by_facet_reasons') or {}).get(str(facet), '')
    else:
        pca = diagnostics.get('pca')
        reason = diagnostics.get('pca_reason', '')
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    loadings = pca.get('loadings') if isinstance(pca, dict) else None
    if not isinstance(loadings, pd.DataFrame) or loadings.empty or 'PC1' not in loadings.columns:
        msg = str(reason) if reason else 'Residual PCA loadings unavailable: set compute_pca=True with at least 2 valid residual columns.'
        _empty(ax, msg)
        fig.suptitle('Standardized residual PCA · Contrast 1 loadings')
        return maps.style_figure(fig)
    s = pd.to_numeric(loadings['PC1'], errors='coerce').dropna()
    if s.empty:
        _empty(ax, 'No finite Contrast 1 loadings available')
        fig.suptitle('Standardized residual PCA · Contrast 1 loadings')
        return maps.style_figure(fig)
    if len(s) > max(1, int(top_n)):
        s = s.loc[s.abs().sort_values(ascending=False).head(int(top_n)).index]
    s = s.sort_values(kind='mergesort')
    ys = np.arange(len(s))
    vals = s.to_numpy(dtype=float)
    ax.axvline(0, color='#64748b', linewidth=.8)
    for gx in (-0.4, 0.4):
        ax.axvline(gx, color='#8a98a7', linestyle='--', linewidth=.85)
    for y_i, v_i in zip(ys, vals):
        ax.plot([0, v_i], [y_i, y_i], color='#cbd5e1', linewidth=1.2, zorder=2)
    pos = vals >= 0
    if np.any(pos):
        ax.scatter(vals[pos], ys[pos], color='#2176ae', s=54, zorder=3, label='Positive pole (PC1 >= 0)')
    if np.any(~pos):
        ax.scatter(vals[~pos], ys[~pos], color='#d88627', s=54, zorder=3, label='Negative pole (PC1 < 0)')
    ax.set_yticks(ys, [str(idx)[:26] for idx in s.index], fontsize=8)
    ax.set(xlabel='Contrast 1 loading (PC1)', xlim=(-1.05, 1.05))
    ax.legend(fontsize=8, loc='best', frameon=False)
    ev1 = float(pca['eigenvalues'][0]) if len(pca.get('eigenvalues', ())) else np.nan
    ev_txt = f' (EV1 = {ev1:.2f})' if np.isfinite(ev1) else ''
    fig.suptitle('Standardized residual PCA · Contrast 1 loadings' + (f' · {facet}' if facet else '') + ev_txt)
    fig.supxlabel('Inspect substantive content contrasting opposite poles (±0.4 guides) alongside missingness overlap.', fontsize=8)
    return maps.style_figure(fig)


def _bias_ranked_figure(analysis, *, index=0, top_n=30):
    from mfrm_app import bias_analysis as bias
    result, _ = _require_mfrm(analysis, 'bias_ranked')
    entries = analysis.get('bias_analyses') or []
    if not entries:
        raise ValueError('No bias_analyses recorded in this project')
    entry = entries[index]
    s = entry['settings']
    label = s['facet_a'] + ' × ' + s['facet_b']
    table = bias._bias_table_for_publication({label: entry['result']}, result['config']['model'])
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    if table is None or table.empty:
        _empty(ax, entry['result'].get('_skip_reason', 'No computed cells'))
        fig.suptitle(label + ' · ranked conditional shifts')
        fig.supxlabel('Tests / intervals withheld.', fontsize=8)
        return maps.style_figure(fig)
    interior = table.loc[table['Calculation'].eq('Interior solution')].copy()
    interior['Bias Size'] = pd.to_numeric(interior['Bias Size'], errors='coerce')
    interior = interior.loc[np.isfinite(interior['Bias Size'])].copy()
    omitted = len(table) - len(interior)
    if interior.empty:
        _empty(ax, 'No interior conditional bias estimates available')
        fig.suptitle(label + ' · ranked conditional shifts')
        fig.supxlabel(f'Omitted {omitted} non-interior cell(s). Tests / intervals withheld.', fontsize=8)
        return maps.style_figure(fig)
    if len(interior) > max(1, int(top_n)):
        interior = interior.assign(_ab=interior['Bias Size'].abs()).sort_values('_ab', ascending=False).head(int(top_n)).drop(columns=['_ab'])
    interior = interior.sort_values('Bias Size', kind='mergesort').reset_index(drop=True)
    ys = np.arange(len(interior))
    vals = interior['Bias Size'].to_numpy(dtype=float)
    ax.axvline(0, color='#64748b', linewidth=.8)
    for gx in (-0.5, 0.5):
        ax.axvline(gx, color='#8a98a7', linestyle='--', linewidth=.85)
    for y_i, v_i in zip(ys, vals):
        ax.plot([0, v_i], [y_i, y_i], color='#cbd5e1', linewidth=1.2, zorder=2)
    pos = vals >= 0
    if np.any(pos):
        ax.scatter(vals[pos], ys[pos], color='#d88627', s=52, zorder=3, label='Positive shift')
    if np.any(~pos):
        ax.scatter(vals[~pos], ys[~pos], color='#2176ae', s=52, zorder=3, label='Negative shift')
    cell_labels = [
        f"{r['Level A']} × {r['Level B']} (n={int(r.get('Ratings', r.get('Observed Count', 0)))})"[:30]
        for _, r in interior.iterrows()
    ]
    ax.set_yticks(ys, cell_labels, fontsize=8)
    lim = max(0.75, float(np.max(np.abs(vals))) * 1.15)
    ax.set(xlabel='Conditional eta shift (logits)', xlim=(-lim, lim))
    ax.legend(fontsize=8, loc='best', frameon=False)
    fig.suptitle(label + ' · ranked conditional shifts')
    omit_txt = f' Omitted {omitted} non-interior cell(s).' if omitted else ''
    fig.supxlabel(f'Descriptive cell estimates with main calibration held fixed; tests / intervals withheld.{omit_txt}', fontsize=8)
    return maps.style_figure(fig)


def analysis_figure(analysis, kind, *, level=None, index=0, facet=None, statistic='Infit', presentation='explained', gd_swap=False, top_n=40):
    from mfrm_app.figure_presentation import matplotlib_figure
    project.require_reportable_input(analysis)
    canonical = KIND_ALIASES.get(str(kind), str(kind))
    fig = _analysis_figure(analysis, canonical, level=level, index=index, facet=facet, statistic=statistic, gd_swap=gd_swap, top_n=top_n)
    return matplotlib_figure(fig, presentation)


def _analysis_figure(analysis, kind, *, level=None, index=0, facet=None, statistic='Infit', gd_swap=False, top_n=40):
    """Draw retained values only; unavailable intervals never become error bars."""
    if kind == 'pathway':
        result, diagnostics = _require_mfrm(analysis, kind)
        return maps.pathway_figure(result, diagnostics, facet=facet, statistic=statistic)
    if kind == 'facet':
        return _facet_estimates_figure(analysis, facet=facet, rater_default=False, top_n=top_n)
    if kind == 'rater_severity':
        return _facet_estimates_figure(analysis, facet=facet, rater_default=True, top_n=top_n)
    if kind == 'step':
        return _step_ladder_figure(analysis, level=level)
    if kind == 'bubble':
        return _bubble_figure(analysis, facet=facet, statistic=statistic, top_n=top_n)
    if kind == 'expected':
        return _expected_score_figure(analysis, level=level)
    if kind == 'information':
        return _information_figure(analysis, level=level)
    if kind == 'fair_average':
        return _fair_average_figure(analysis, facet=facet, statistic=statistic, top_n=top_n)
    if kind == 'unexpected':
        return _unexpected_figure(analysis, top_n=top_n)
    if kind == 'residual_qq':
        return _residual_qq_figure(analysis)
    if kind == 'pca_loadings':
        return _pca_loadings_figure(analysis, facet=facet, top_n=top_n)
    if kind == 'bias_ranked':
        return _bias_ranked_figure(analysis, index=index, top_n=top_n)
    if kind in {'wright', 'category'}:
        result, _ = _require_mfrm(analysis, kind)
        steps, label = _steps(result, level)
        if not np.isfinite(steps).all():
            raise ValueError('Category steps must be finite')
        if kind == 'wright':
            fig = maps.wright_figure(result)
            ax = fig.axes[1]
            for step in steps:
                ax.axhline(float(step), color='#c07125', linestyle=':', linewidth=.8)
            lo, hi = ax.get_ylim()
            ax.set_ylim(min(lo, float(steps.min())-.2), max(hi, float(steps.max())+.2))
            ax.set_title('Facet locations / step offsets\n' + label, fontsize=10)
            return fig
        return maps.category_probability_figure(result, steps, title=label + ' · category probabilities', figsize=(9, 4.8))
    fig = Figure(figsize=(9, 4.8), layout='constrained')
    ax = fig.subplots()
    if kind == 'bias':
        from mfrm_app import bias_analysis as bias
        entry = analysis['bias_analyses'][index]
        s = entry['settings']; label = s['facet_a'] + ' × ' + s['facet_b']
        table = bias._bias_table_for_publication({label: entry['result']}, analysis['result']['config']['model'])
        if table is None:
            _empty(ax, entry['result'].get('_skip_reason', 'No computed cells'))
        else:
            # Non-interior cells remain visible as gaps; their reason is in the table.
            values = table.assign(value=table['Bias Size'].where(table.Calculation.eq('Interior solution')))
            pivot = values.pivot(index='Level A', columns='Level B', values='value')
            array = np.ma.masked_invalid(pivot.to_numpy(dtype=float))
            from matplotlib import colormaps
            limit = max(1., float(np.max(np.abs(array)))) if array.count() else 1.
            cmap = colormaps['RdBu_r'].with_extremes(bad='#cbd5e1')
            artist = ax.imshow(array, cmap=cmap, vmin=-limit, vmax=limit, aspect='auto')
            if array.size <= 100:
                for (y, x), value in np.ndenumerate(array.data):
                    if not np.ma.getmaskarray(array)[y, x]:
                        ax.text(x, y, f'{value:.2f}', ha='center', va='center',
                                color='white' if abs(value) > .65*limit else '#172330', fontsize=9)
            ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=30, ha='right')
            ax.set_yticks(range(len(pivot.index)), pivot.index)
            ax.set(xlabel=s['facet_b'], ylabel=s['facet_a'])
            fig.colorbar(artist, ax=ax, label='Conditional eta shift (logits)')
        fig.suptitle(label + ' · descriptive cell estimates')
        fig.supxlabel('Tests / intervals withheld. Gaps: no interior estimate; see cell status in table.', fontsize=8)
    elif kind == 'interval':
        entry = analysis['interval_analyses'][index]
        try:
            if entry['status'] != 'computed':
                raise ValueError(entry['reason'])
            table = intervals.read_structural_intervals(analysis['result'], entry['record'])
            available = table.loc[table.CIEligible.eq(True)]
            if available.empty:
                _empty(ax, 'Intervals withheld: ' + '; '.join(table.Reason.unique()))
            else:
                x = available.Estimate.to_numpy()
                ax.errorbar(x, range(len(x)), xerr=np.array([x-available.Lower, available.Upper-x]), fmt='o', capsize=4)
                ax.set_yticks(range(len(x)), available.Target)
                ax.set_xlabel('Logits'); ax.invert_yaxis()
            fig.suptitle(f"{100*entry['settings']['level']:g}% intervals · {entry['settings']['facet']}")
            fig.supxlabel(textwrap.fill(intervals.CAUTION, 110), fontsize=8)
        except ValueError as exc:
            _empty(ax, 'Withheld: ' + str(exc))
    elif kind == 'joint':
        entry = analysis['joint_analysis']
        try:
            if entry['status'] != 'computed':
                raise ValueError(entry['reason'])
            table = joint.read_joint_interaction(analysis['result'], entry['record'])
            available = table[table.CIEligible.eq(True)]
            if available.empty:
                _empty(ax, 'Withheld: '+', '.join(sorted(set(table.Reason))))
            else:
                x = available.Estimate.to_numpy(dtype=float)
                ax.errorbar(x, range(len(x)), xerr=np.array([x-available.Lower, available.Upper-x]), fmt='o', capsize=4)
                ax.set_yticks(range(len(x)), available.Target)
                ax.axvline(0, color='gray', linewidth=.7)
                ax.set_xlabel('Additional predictor shift (logits)'); ax.invert_yaxis()
            fig.suptitle(f"Joint interaction · {100*entry['settings']['level']:g}% pointwise intervals")
            fig.supxlabel(textwrap.fill(joint.CAUTION, 110), fontsize=8)
        except ValueError as exc:
            _empty(ax, 'Withheld: '+str(exc))
    elif kind == 'gd':
        entry = analysis['gd_analysis']; result = entry['result']
        if not result['available']:
            _empty(ax, result['reason'])
        else:
            table = result['d_study']; facets = list(entry['settings']['random_facets'])
            if gd_swap and len(facets) == 2:
                facets = [facets[1], facets[0]]
            if len(facets) == 2:
                unique_levels = sorted(int(v) for v in table[facets[1]].dropna().unique())
                if len(unique_levels) > 6:
                    obs_val = int(result['design']['observed_n'].get(facets[1], unique_levels[0]))
                    sel_val = int(result['selected_plan'].get(facets[1], obs_val))
                    mandatory = sorted({unique_levels[0], unique_levels[-1], obs_val, sel_val} & set(unique_levels))
                    step = max(1, (len(unique_levels) - 1) // 5)
                    candidates = [unique_levels[idx] for idx in range(0, len(unique_levels), step)] + unique_levels
                    fill = []
                    for val in candidates:
                        if val not in mandatory and val not in fill:
                            fill.append(val)
                        if len(mandatory) + len(fill) >= 6:
                            break
                    plot_levels = sorted((mandatory + fill)[:max(6, len(mandatory))])
                    plot_table = table.loc[table[facets[1]].isin(plot_levels)]
                else:
                    plot_table = table
                groups = plot_table.groupby(facets[1], sort=True)
            else:
                groups = [(None, table)]
            for group, rows in groups:
                rows = rows.sort_values(facets[0])
                suffix = f' · {facets[1]}={group}' if group is not None else ''
                line, = ax.plot(rows[facets[0]], rows.G, marker='o', label='G' + suffix)
                ax.plot(rows[facets[0]], rows.Phi, marker='x', linestyle='--', color=line.get_color(), label='Phi' + suffix)
            selected = result['selected_plan']
            ax.scatter([selected[facets[0]]]*2, [selected['G'], selected['Phi']], marker='*', s=130, color='black', label='Selected plan', zorder=5)
            ax.set(xlabel=facets[0] + ' count', ylabel='Coefficient', ylim=(0, 1.03))
            ax.legend(fontsize=8, ncols=2)
        fc = entry['settings'].get('fixed_condition')
        cond_suffix = ''
        if fc:
            cond_suffix = f" · {fc['facet']}={fc['level']}" if fc['mode'] == 'level' else f" · {fc['facet']} ({fc['mode']})"
        fig.suptitle('D-Study · original observed scores' + cond_suffix)
        fig.supxlabel('G: relative decisions; Phi: absolute decisions. Point projections; intervals withheld.', fontsize=8)
    elif kind == 'pca':
        pca = analysis.get('diagnostics', {}).get('pca')
        if not isinstance(pca, dict) or len(pca.get('eigenvalues', [])) == 0:
            _empty(ax, 'Residual PCA unavailable: fewer than 2 valid facet-combination residual columns.')
        else:
            ev = np.asarray(pca['eigenvalues'], dtype=float)
            k = min(15, len(ev))
            xs = np.arange(1, k + 1)
            ys = ev[:k]
            ax.bar(xs, ys, color='#2563eb', alpha=0.82, width=0.62, label='Eigenvalue')
            ax.plot(xs, ys, color='#0f172a', marker='o', linewidth=1.4, markersize=4.5)
            ax.axhline(2.0, color='#b45309', linestyle='--', linewidth=1.1, label='Inspection prompt (EV = 2.0)')
            ax.axhline(1.0, color='#64748b', linestyle=':', linewidth=0.9, label='Unit variance (EV = 1.0)')
            ax.set_xticks(xs, [f'PC{i}' for i in xs], rotation=0 if k <= 10 else 30)
            ax.set(
                xlabel='Residual principal contrast',
                ylabel='Eigenvalue (standardized residual units)',
                ylim=(0, max(2.4, float(ys.max()) * 1.15)),
            )
            ax.legend(fontsize=8, loc='upper right')
        fig.suptitle('Standardized residual PCA · scree plot')
        fig.supxlabel(
            'Eigenvalue sum equals residual column count. Inspect substantive contrast poles alongside missingness overlap.',
            fontsize=8,
        )
    else:
        raise ValueError('Unknown desktop figure')
    return maps.style_figure(fig)


def report_html(analysis, *, pathway_facet=None, pathway_statistic='Infit', facet=None, statistic=None,
                presentation=None, language=None, gd_swap=None):
    """Offline HTML with every requested table and explicitly labelled figure."""
    from mfrm_app.figure_presentation import caption_html
    from mfrm_app.desktop_i18n import text
    if facet is not None and pathway_facet is not None and facet != pathway_facet:
        raise ValueError('Specify either facet or pathway_facet, not both')
    if statistic is not None and pathway_statistic != 'Infit' and statistic != pathway_statistic:
        raise ValueError('Specify either statistic or pathway_statistic, not both')
    pathway_facet = facet if facet is not None else pathway_facet
    pathway_statistic = statistic if statistic is not None else pathway_statistic
    view = analysis.get('view_settings', {})
    presentation = presentation or view.get('presentation', 'explained')
    language = language or view.get('language', 'en')
    gd_swap = (view.get('gd_swap') == '1') if gd_swap is None else bool(gd_swap)
    tables = current_tables(analysis)
    escape = html.escape
    kind = ('G/D-Study' if analysis['kind'] == 'gd-analysis-v2' else
            analysis['fit_settings']['model'] + ' / ' + analysis['fit_settings']['method'])
    parts = [f'<!doctype html><html lang="{language}"><meta charset="utf-8"><title>MFRM analysis</title>',
        '<style>body{font:16px sans-serif;margin:2em;max-width:1200px}table{border-collapse:collapse;font-size:12px}',
        'td,th{padding:.4em;border-bottom:1px solid #bbb}section{overflow-x:auto}svg{max-width:100%;height:auto}pre{white-space:pre-wrap}</style>',
        '<h1>' + escape(kind) + '</h1><p>Analysis ID: ' + escape(analysis['analysis_id']) + '</p>',
        '<p>' + escape(text('プロジェクトには入力データと識別子を含みます。表示・報告では再推定しません。', language)) + '</p>',
        '<h2>Settings</h2><pre>' + escape(json.dumps(request_settings(analysis), ensure_ascii=False, indent=2)) + '</pre>']
    if 'result' in analysis:
        result = analysis['result']
        current = basic.qualification_snapshot(result)
        original = analysis.get('qualification_at_save')
        for note in (
            basic.reporting_facts.optimizer_termination_text(result),
            basic.reporting_facts.recorded_method_details(result),
            f"Run qualification: state={current['run_state']}; SupportsFormalInference=FALSE; "
            f"formal_hold_reason={current['calibration_hold_reason'] or 'none'}; "
            f"screening_block_reason={current['screening_block_reason'] or 'none'}.",
            current['calibration_hold_message'],
        ):
            if note:
                parts.append('<p>' + escape(note) + '</p>')
        if isinstance(original, dict) and (
            original.get('calibration_hold_reason') != current['calibration_hold_reason']
            or original.get('calibration_hold_message') != current['calibration_hold_message']
        ):
            saved_reason = str(original.get('calibration_hold_reason') or 'none')
            saved_msg = str(original.get('calibration_hold_message') or '')
            parts.append(
                '<p>' + escape(
                    f"Saved qualification at archive creation (retained): reason={saved_reason}"
                    + (f" — {saved_msg}" if saved_msg else "")
                ) + '</p>'
            )
        parts.append('<p>' + escape(text('個人得点と一般表の利用条件を、指定区間の利用条件から判断しないでください。自動的な同等性の結論は保留です。', language)) + '</p>')
        if 'joint_analysis' in analysis:
            parts.append('<p>' + escape(text('通常の表・Wright map等は元の主効果モデルです。共同交互作用は別モデルとして再推定した結果を専用の表・図に示します。', language)) + '</p>')
    for title, table, note in tables:
        parts.extend([
            '<h2>' + escape(text(title, language)) + '</h2><section>',
            table.to_html(index=False, escape=True),
            '</section><p>' + escape(text(note, language)) + '</p>',
        ])
    for i, (title, options) in enumerate(figure_choices(analysis)):
        options = dict(options)
        if options['kind'] == 'pathway':
            options.update(facet=pathway_facet, statistic=pathway_statistic)
        elif options['kind'] == 'gd':
            options.update(gd_swap=gd_swap)
        fig = analysis_figure(analysis, **options, presentation=presentation)
        svg = basic.figure_bytes(analysis, 'figure_' + str(i), figure=fig).decode('utf-8')
        parts.extend(['<h2>' + escape(text(title, language)) + '</h2>', svg[svg.index('<svg'):]])
        if presentation == 'manuscript':
            parts.append(caption_html(fig._mfrm_caption))
        if options['kind'] in maps.MAP_GUIDES:
            parts.append('<p>' + escape((maps.MAP_GUIDES if language == 'ja' else maps.MAP_GUIDES_EN)[options['kind']][1]) + '</p>')
        elif options['kind'] in DIAGNOSTIC_GUIDES.get(language, {}):
            parts.append('<p>' + escape(DIAGNOSTIC_GUIDES[language][options['kind']]) + '</p>')
    parts.append('<h2>Source</h2><pre>' + escape(json.dumps(analysis['source'], ensure_ascii=False, indent=2)) + '</pre></html>')
    return '\n'.join(parts)
