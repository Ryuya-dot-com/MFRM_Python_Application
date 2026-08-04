"""Pure presentation policy for the Streamlit research workflow.

The app has many statistical surfaces, but the page should still answer one
question at a time.  This module keeps the small amount of state that controls
progressive disclosure independent from Streamlit and from the estimator.

It deliberately does *not* decide scientific readiness.  It only describes
which setup surfaces should be prominent for the current data source and
workflow phase; readiness remains governed by the evidence contracts.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class DataSourceKind(str, Enum):
    """Stable data-origin categories used by the presentation layer."""

    SAMPLE = "sample"
    SIMULATION = "simulation"
    PASTE = "paste"
    UPLOAD = "upload"
    UNKNOWN = "unknown"


class WorkspacePhase(str, Enum):
    """Coarse journey phase; separate from statistical readiness."""

    CHOOSE_DATA = "choose_data"
    REVIEW_AND_RUN = "review_and_run"
    READ_RESULTS = "read_results"


@dataclass(frozen=True, slots=True)
class WorkspacePresentation:
    """Progressive-disclosure policy for one Streamlit rerun."""

    phase: WorkspacePhase
    collapse_input: bool
    show_pre_run_checks: bool


def classify_data_source(option_id: object) -> DataSourceKind:
    """Classify the stable data-source widget value without display labels."""

    value = str(option_id or "").strip()
    if value.startswith("scenario:"):
        return DataSourceKind.SAMPLE
    if value == "simulate":
        return DataSourceKind.SIMULATION
    if value == "paste":
        return DataSourceKind.PASTE
    if value == "upload":
        return DataSourceKind.UPLOAD
    return DataSourceKind.UNKNOWN


def data_source_may_contain_user_data(source: DataSourceKind | str) -> bool:
    """Return whether the source can carry a user's real rating records.

    Unknown values fail closed: an unrecognized future source receives the
    stronger privacy treatment until it is explicitly classified.
    """

    try:
        kind = source if isinstance(source, DataSourceKind) else DataSourceKind(source)
    except (TypeError, ValueError):
        kind = DataSourceKind.UNKNOWN
    return kind not in {DataSourceKind.SAMPLE, DataSourceKind.SIMULATION}


def workspace_presentation(
    *,
    has_data: bool,
    has_result: bool,
) -> WorkspacePresentation:
    """Return the page-density policy for the current workflow phase."""

    if has_result:
        return WorkspacePresentation(
            phase=WorkspacePhase.READ_RESULTS,
            collapse_input=True,
            show_pre_run_checks=False,
        )
    if has_data:
        return WorkspacePresentation(
            phase=WorkspacePhase.REVIEW_AND_RUN,
            collapse_input=False,
            show_pre_run_checks=True,
        )
    return WorkspacePresentation(
        phase=WorkspacePhase.CHOOSE_DATA,
        collapse_input=False,
        show_pre_run_checks=False,
    )
