"""Regression tests for `render_readiness_panel`.

The pre-run readiness panel previously referenced undefined ``icon`` and
``status`` names while formatting its banner, and its only caller wrapped the
call in a bare ``except Exception: pass``. The result: the banner raised a
``NameError`` on every render and was silently swallowed, so users never saw
the pre-estimation readiness gate at all. These tests render the panel for each
severity and assert the banner actually appears with its status label.
"""

from __future__ import annotations

import streamlit_app as app


class _ReadinessRenderRecorder:
    """Stand-in for the module-level ``st`` that records banner text.

    Unknown attributes (e.g. ``session_state``) delegate to the real Streamlit
    module so ``t()`` keeps resolving translations exactly as in a live run.
    """

    def __init__(self, real):
        self._real = real
        self.messages: list[str] = []

    def __getattr__(self, name):
        return getattr(self._real, name)

    def success(self, msg="", *a, **k):
        self.messages.append(str(msg))

    def warning(self, msg="", *a, **k):
        self.messages.append(str(msg))

    def error(self, msg="", *a, **k):
        self.messages.append(str(msg))

    def caption(self, *a, **k):
        return None

    def markdown(self, *a, **k):
        return None

    def dataframe(self, *a, **k):
        return None

    def expander(self, *a, **k):
        class _NullCtx:
            def __enter__(self_inner):
                return None

            def __exit__(self_inner, *exc):
                return False

        return _NullCtx()


def _render_with_recorder(report: dict) -> list[str]:
    recorder = _ReadinessRenderRecorder(app.st)
    saved = app.st
    app.st = recorder
    try:
        app.render_readiness_panel(report)
    finally:
        app.st = saved
    return recorder.messages


def _report(overall: str) -> dict:
    return {
        "overall": overall,
        "checks": [{
            "name": "n_obs",
            "severity": overall,
            "headline": "render regression",
            "detail": "synthetic check",
        }],
        "n_warnings": 1,
        "n_issues": 1,
    }


def test_panel_renders_each_severity_with_label():
    for overall, label in (("ok", "OK"), ("warning", "CAUTION"), ("issue", "ISSUE")):
        messages = _render_with_recorder(_report(overall))
        assert any(f"[{label}]" in m for m in messages), (
            f"readiness panel did not render the '{overall}' banner"
        )


def test_panel_does_not_raise_nameerror():
    # Exercises the exact code path that the caller used to swallow.
    for overall in ("ok", "warning", "issue"):
        _render_with_recorder(_report(overall))


def test_registered_self_test_passes():
    # The in-app --self-test guard for the same regression must pass.
    app._self_test_readiness_panel_renders()
