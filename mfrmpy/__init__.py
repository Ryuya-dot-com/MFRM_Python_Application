"""Public Python interface to the shared MFRM application calculations.

Use :func:`analyze` for RSM/PCM or observed-score G/D studies, :func:`tables`
for qualified results, and :func:`save_project` / :func:`load_project` for
portable .mfrm archives. Raw retained results do not override output limits.
"""
from importlib.resources import files
from pathlib import Path

import pandas as pd

from mfrm_app import basic_analysis, desktop_analysis
from mfrm_app.desktop_analysis import run_analysis as analyze
from mfrm_app.project import ProjectError, load_project, save_project

__version__ = "0.1.0a2"

__all__ = [
    "__version__", "analyze", "tables", "load_example", "save_project",
    "load_project", "ProjectError", "plot", "report_html",
]


def tables(analysis):
    """Return (title, DataFrame, interpretation note) tuples without refitting.

    Current output-specific checks apply to retained results. Unavailable
    uncertainty stays unavailable; inspect each table's Use/Reason fields and
    accompanying note. Do not substitute raw archived SEs for withheld output.
    """
    basic = [] if analysis.get("kind") == "gd-analysis-v2" else basic_analysis.current_tables(analysis)
    return basic + desktop_analysis.current_additional_tables(analysis)


def load_example():
    """Return a fresh copy of the synthetic 0–5 writing ratings as a DataFrame."""
    resource = files(__package__).joinpath("data", "writing.csv")
    if not resource.is_file():
        resource = Path(__file__).resolve().parents[1] / "examples" / "desktop_writing.csv"
    with resource.open("r", encoding="utf-8") as handle:
        return pd.read_csv(handle)


def plot(analysis, kind="wright", *, level=None, index=0, facet=None,
         statistic="Infit", presentation="explained", gd_swap=False):
    """Return a Matplotlib Figure; install the ``figures`` extra first.

    Kinds and output restrictions follow the desktop figures. PCM Wright maps
    and category curves require ``level``. No estimator is rerun to draw results.
    """
    from mfrm_app.desktop_outputs import analysis_figure
    return analysis_figure(analysis, kind, level=level, index=index, facet=facet,
                           statistic=statistic, presentation=presentation, gd_swap=gd_swap)


def report_html(analysis, *, pathway_facet=None, pathway_statistic="Infit",
                presentation=None, language=None, gd_swap=None):
    """Return an offline HTML report; requires the ``figures`` extra.

    The caller chooses where to save it. Reports and .mfrm archives can contain
    the original ratings and identifiers; handle them as research data.
    """
    from mfrm_app.desktop_outputs import report_html as render
    return render(analysis, pathway_facet=pathway_facet,
                  pathway_statistic=pathway_statistic, presentation=presentation,
                  language=language, gd_swap=gd_swap)
