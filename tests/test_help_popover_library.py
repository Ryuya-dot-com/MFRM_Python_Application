"""Contract tests for the just-in-time help popover library.

Each of the 18 topics in `_HELP_POPOVER_LIBRARY` powers a "❓ How to
read" popover next to a specific diagnostic chart. Empty / missing /
mis-keyed entries degrade the UX silently, so these tests pin the
shape every entry must satisfy.
"""

from __future__ import annotations

import pytest

import streamlit_app as app


REQUIRED_FIELDS = ("title", "what", "how")
OPTIONAL_FIELDS = ("watch",)


def test_library_has_at_least_eighteen_topics():
    assert len(app._HELP_POPOVER_LIBRARY) >= 18


@pytest.mark.parametrize("topic_key", list(app._HELP_POPOVER_LIBRARY.keys()))
def test_topic_has_required_fields(topic_key: str):
    entry = app._HELP_POPOVER_LIBRARY[topic_key]
    for field in REQUIRED_FIELDS:
        assert field in entry, f"{topic_key!r} missing {field!r}"
        assert isinstance(entry[field], str) and entry[field].strip(), (
            f"{topic_key!r}: {field!r} is empty"
        )


@pytest.mark.parametrize("topic_key", list(app._HELP_POPOVER_LIBRARY.keys()))
def test_what_field_is_substantive(topic_key: str):
    entry = app._HELP_POPOVER_LIBRARY[topic_key]
    assert len(entry["what"]) >= 30, (
        f"{topic_key!r}: 'what' too short ({len(entry['what'])} chars)"
    )


@pytest.mark.parametrize("topic_key", list(app._HELP_POPOVER_LIBRARY.keys()))
def test_how_field_is_substantive(topic_key: str):
    entry = app._HELP_POPOVER_LIBRARY[topic_key]
    assert len(entry["how"]) >= 50, (
        f"{topic_key!r}: 'how' too short ({len(entry['how'])} chars)"
    )


def test_v026_topics_are_present():
    """v0.2.6-beta added 5 new topics; ensure they stayed in the library."""
    v026_topics = {
        "misfit_ranking", "facet_distribution", "rater_agreement",
        "threshold_map", "zstd_distribution",
    }
    assert v026_topics.issubset(app._HELP_POPOVER_LIBRARY.keys())


def test_core_visual_topics_are_present():
    """Core diagnostic topics must exist so every chart can be wired up."""
    core_topics = {
        "wright_map", "category_probability", "pathway_map", "scree",
        "fit_scatter", "bias_heatmap", "forest_measures", "qq_residuals",
        "ecdf_measures", "posterior_trace", "posterior_rhat_ess",
        "coverage_heatmap", "category_usage",
    }
    assert core_topics.issubset(app._HELP_POPOVER_LIBRARY.keys())


def test_render_helper_handles_unknown_key_without_raising():
    """`render_help_popover` must no-op on unknown keys (not crash the page)."""
    # The helper calls `st.popover`, which fails outside a ScriptRunContext.
    # We only verify the look-up path: unknown key → early return.
    try:
        app.render_help_popover("__definitely_not_a_topic__")
    except KeyError:
        pytest.fail("render_help_popover raised KeyError on unknown topic key")
    except Exception:
        # Streamlit's popover may raise outside a ScriptRunContext — that
        # is fine; we only care that unknown topic keys do not raise.
        pass
