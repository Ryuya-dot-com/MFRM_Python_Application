"""Separate manuscript titles/notes from graphics without changing plotted data.

APA figure components: https://apastyle.apa.org/style-grammar-guidelines/tables-figures/figures
Panel labels, legends, data labels and unavailable-result messages stay in the image.
"""
import html

MODES = ('explained', 'manuscript')


def mode(value):
    if value not in MODES:
        raise ValueError('Unknown figure presentation mode')
    return value


def matplotlib_figure(fig, presentation='explained'):
    """Apply to a newly constructed figure; retain a separate caption in both modes."""
    mode(presentation)
    title = getattr(fig, '_suptitle', None)
    note = getattr(fig, '_supxlabel', None)
    titles = [title.get_text()] if title and title.get_text() else []
    notes = [note.get_text()] if note and note.get_text() else []
    axes = [ax for ax in fig.axes if ax.get_label() != '<colorbar>']
    # Titles identifying panels of a multi-panel map remain essential labels.
    if len(axes) == 1 and axes[0].get_title():
        titles.append(axes[0].get_title())
        if presentation == 'manuscript':
            axes[0].set_title('')
    fig._mfrm_caption = dict(title='; '.join(titles), notes=notes)
    if presentation == 'manuscript':
        if title: title.set_visible(False)
        if note: note.set_visible(False)
    return fig


def plotly_figure(fig, presentation='explained'):
    """Copy before removing only overall titles and outside-plot explanatory notes."""
    import plotly.graph_objects as go
    mode(presentation)
    out = go.Figure(fig)
    meta = dict(out.layout.meta) if isinstance(out.layout.meta, dict) else {}
    prior = meta.get('mfrm_presentation', {})
    title = out.layout.title.text or prior.get('title', '')
    notes, keep = list(prior.get('notes', [])), []
    for annotation in out.layout.annotations or ():
        outside = (annotation.xref == 'paper' and annotation.yref == 'paper' and isinstance(annotation.y, (int, float))
                   and (annotation.y < 0 or annotation.y > 1 or (annotation.yshift or 0) < -30))
        explanatory = annotation.name == 'mfrm_note' or (outside and not annotation.showarrow)
        if explanatory:
            if annotation.text not in notes: notes.append(annotation.text)
        if presentation != 'manuscript' or not explanatory: keep.append(annotation)
    meta['mfrm_presentation'] = dict(mode=presentation, title=title, notes=notes)
    out.update_layout(meta=meta)
    if presentation == 'manuscript':
        out.layout.title = None
        out.layout.annotations = tuple(keep)
        out.update_layout(margin_t=40, margin_b=65)
    return out


def caption_text(caption, analysis_id=''):
    parts = [caption.get('title', '')]
    parts.extend('Note. ' + note for note in caption.get('notes', []) if note)
    if analysis_id: parts.append('Analysis ID: ' + analysis_id)
    return '\n\n'.join(part for part in parts if part)


def caption_html(caption):
    title = caption.get('title', '')
    parts = ['<p><em>' + html.escape(title) + '</em></p>'] if title else []
    parts.extend('<p><em>Note.</em> ' + html.escape(note) + '</p>' for note in caption.get('notes', []) if note)
    return '\n'.join(parts)
