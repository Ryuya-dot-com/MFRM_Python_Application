"""Native figures and offline reports from retained, qualified analysis records."""
import html
import json
import textwrap

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from mfrm_app import basic_analysis as basic, desktop_analysis as desktop, engine, project
from mfrm_app import structural_intervals as intervals, visual_maps as maps
from mfrm_app import joint_interaction as joint


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


def current_tables(analysis):
    return ([] if analysis['kind'] == 'gd-analysis-v2' else basic.current_tables(analysis)) + desktop.current_additional_tables(analysis)


def figure_choices(analysis):
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
    for i, entry in enumerate(analysis.get('bias_analyses', [])):
        s = entry['settings']
        choices.append((f"バイアス · {s['facet_a']} × {s['facet_b']}", dict(kind='bias', index=i)))
    for i, entry in enumerate(analysis.get('interval_analyses', [])):
        s = entry['settings']
        suffix = ' / ' + ' − '.join(s['difference']) if s['difference'] else ''
        choices.append((f"指定区間 {i+1} · {s['facet']}{suffix}", dict(kind='interval', index=i)))
    if 'gd_analysis' in analysis:
        choices.append(('D-Study · G / Phi', dict(kind='gd')))
    if 'joint_analysis' in analysis:
        s = analysis['joint_analysis']['settings']
        choices.append((f"共同交互作用 · {s['facet_a']} × {s['facet_b']}", dict(kind='joint')))
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


def analysis_figure(analysis, kind, *, level=None, index=0, facet=None, statistic='Infit', presentation='explained', gd_swap=False):
    from mfrm_app.figure_presentation import matplotlib_figure
    fig = _analysis_figure(analysis, kind, level=level, index=index, facet=facet, statistic=statistic, gd_swap=gd_swap)
    return matplotlib_figure(fig, presentation)


def _analysis_figure(analysis, kind, *, level=None, index=0, facet=None, statistic='Infit', gd_swap=False):
    """Draw retained values only; unavailable intervals never become error bars."""
    if kind == 'pathway':
        return maps.pathway_figure(analysis['result'], analysis['diagnostics'], facet=facet, statistic=statistic)
    if kind in {'wright', 'category'}:
        result = analysis['result']
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
        fig = Figure(figsize=(9, 4.8), layout='constrained')
        ax = fig.subplots()
        eta = np.linspace(min(-4., float(steps.min())-2), max(4., float(steps.max())+2), 201)
        if result['config']['model'] == 'PCM':
            cumulative = np.r_[0., np.cumsum(steps)][None, :]
            probabilities = engine.category_prob_pcm(eta, cumulative, np.zeros(len(eta), dtype=int))
        else:
            probabilities = engine.category_prob_rsm(eta, np.r_[0., np.cumsum(steps)])
        for i, values in enumerate(probabilities.T):
            ax.plot(eta, values, label=str(i + result['prep']['rating_min']))
        ax.set(xlabel='Net predictor (logits)', ylabel='Category probability', ylim=(0, 1))
        ax.legend(title='Score', ncols=min(6, probabilities.shape[1]))
        fig.suptitle(label + ' · category probabilities')
        fig.supxlabel(textwrap.fill(maps.map_caution(result), 110), fontsize=8)
        return maps.style_figure(fig)
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
    else:
        raise ValueError('Unknown desktop figure')
    return maps.style_figure(fig)


def report_html(analysis, *, pathway_facet=None, pathway_statistic='Infit', presentation=None, language=None, gd_swap=None):
    """Offline HTML with every requested table and explicitly labelled figure."""
    from mfrm_app.figure_presentation import caption_html
    from mfrm_app.desktop_i18n import text
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
        for note in (basic.reporting_facts.optimizer_termination_text(result), basic.reporting_facts.recorded_method_details(result),
                     basic.qualification_snapshot(result)['calibration_hold_message']):
            parts.append('<p>' + escape(note) + '</p>')
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
    parts.append('<h2>Source</h2><pre>' + escape(json.dumps(analysis['source'], ensure_ascii=False, indent=2)) + '</pre></html>')
    return '\n'.join(parts)
