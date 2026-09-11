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


DATA_SOURCE_CLASS_IDS = ("sample", "simulate", "paste", "upload")


class WorkspacePhase(str, Enum):
    """Legacy three-phase projection retained for compatibility."""

    CHOOSE_DATA = "choose_data"
    REVIEW_AND_RUN = "review_and_run"
    READ_RESULTS = "read_results"


class WorkflowPhase(str, Enum):
    """Authoritative journey phases, independent from Streamlit widgets."""

    SOURCE = "source"
    SETUP = "setup"
    RUN = "run"
    FIRST_READ = "first_read"
    EVIDENCE_DETAIL = "evidence_detail"
    ARCHIVE = "archive"


class WorkflowSurface(str, Enum):
    """Stable shell surfaces that may own the primary next action."""

    SOURCE_WORKSPACE = "source_workspace"
    SETUP_WORKSPACE = "setup_workspace"
    RUN_STATUS = "run_status"
    RESULT_GOAL_ROUTER = "result_goal_router"
    RESULT_SECTION = "result_section"


@dataclass(frozen=True, slots=True)
class WorkflowShell:
    """Normalized page topology for one application rerun.

    ``primary_action_surface`` is singular by construction.  Result sections
    can still contain domain controls, but only the named shell surface owns
    orientation and advertises the workflow's primary next action.
    """

    phase: WorkflowPhase
    primary_action_surface: WorkflowSurface
    supporting_surfaces: tuple[WorkflowSurface, ...] = ()
    selected_result_section: str | None = None

    @property
    def collapse_setup(self) -> bool:
        """Return whether setup should yield visual priority to later phases."""

        return self.phase in {
            WorkflowPhase.RUN,
            WorkflowPhase.FIRST_READ,
            WorkflowPhase.EVIDENCE_DETAIL,
            WorkflowPhase.ARCHIVE,
        }

    @property
    def show_pre_run_checks(self) -> bool:
        """Return whether readiness checks belong on the visible setup surface."""

        return self.phase is WorkflowPhase.SETUP

    @property
    def show_result_router(self) -> bool:
        """Return whether the result goal router owns shell orientation."""

        return self.primary_action_surface is WorkflowSurface.RESULT_GOAL_ROUTER


@dataclass(frozen=True, slots=True)
class WorkspacePresentation:
    """Progressive-disclosure policy for one Streamlit rerun."""

    phase: WorkspacePhase
    collapse_input: bool
    show_pre_run_checks: bool


RESULT_SECTION_IDS = (
    "start",
    "first_read",
    "results",
    "diagnostics",
    "figures",
    "report_export",
    "learn",
)
_EVIDENCE_DETAIL_SECTION_IDS = frozenset({
    "results",
    "diagnostics",
    "figures",
    "learn",
})


def resolve_workflow_shell(
    *,
    has_data: bool,
    has_result: bool,
    estimation_running: bool = False,
    selected_result_section: object = None,
) -> WorkflowShell:
    """Resolve the single workflow-shell contract for the current rerun.

    Display labels and widget keys do not enter this policy. Unknown result
    sections fail back to ``start`` so a stale or future session value cannot
    bypass first-read orientation.
    """

    section = str(selected_result_section or "start").strip()
    if section not in RESULT_SECTION_IDS:
        section = "start"
    if section == "report_export":
        result_phase = WorkflowPhase.ARCHIVE
    elif section in _EVIDENCE_DETAIL_SECTION_IDS:
        result_phase = WorkflowPhase.EVIDENCE_DETAIL
    else:
        result_phase = WorkflowPhase.FIRST_READ

    if estimation_running:
        return WorkflowShell(
            phase=WorkflowPhase.RUN,
            primary_action_surface=WorkflowSurface.RUN_STATUS,
        )
    if has_result:
        return WorkflowShell(
            phase=result_phase,
            primary_action_surface=WorkflowSurface.RESULT_GOAL_ROUTER,
            supporting_surfaces=(WorkflowSurface.RESULT_SECTION,),
            selected_result_section=section,
        )
    if has_data:
        return WorkflowShell(
            phase=WorkflowPhase.SETUP,
            primary_action_surface=WorkflowSurface.SETUP_WORKSPACE,
        )
    return WorkflowShell(
        phase=WorkflowPhase.SOURCE,
        primary_action_surface=WorkflowSurface.SOURCE_WORKSPACE,
    )


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


def data_source_class_id(option_id: object) -> str | None:
    """Project a stable legacy source option onto the two-level chooser."""

    kind = classify_data_source(option_id)
    if kind is DataSourceKind.SAMPLE:
        return "sample"
    if kind is DataSourceKind.SIMULATION:
        return "simulate"
    if kind is DataSourceKind.PASTE:
        return "paste"
    if kind is DataSourceKind.UPLOAD:
        return "upload"
    return None


def data_source_option_id(
    source_class_id: object,
    *,
    scenario_key: object = None,
) -> str | None:
    """Compose the unchanged analysis-facing option ID from chooser state."""

    source_class = str(source_class_id or "").strip()
    if source_class == "sample":
        scenario = str(scenario_key or "").strip()
        return f"scenario:{scenario}" if scenario else None
    if source_class in {"simulate", "paste", "upload"}:
        return source_class
    return None


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
    """Project the authoritative shell onto the legacy setup policy."""

    shell = resolve_workflow_shell(has_data=has_data, has_result=has_result)
    if shell.phase in {
        WorkflowPhase.FIRST_READ,
        WorkflowPhase.EVIDENCE_DETAIL,
        WorkflowPhase.ARCHIVE,
    }:
        return WorkspacePresentation(
            phase=WorkspacePhase.READ_RESULTS,
            collapse_input=shell.collapse_setup,
            show_pre_run_checks=shell.show_pre_run_checks,
        )
    if shell.phase is WorkflowPhase.SETUP:
        return WorkspacePresentation(
            phase=WorkspacePhase.REVIEW_AND_RUN,
            collapse_input=shell.collapse_setup,
            show_pre_run_checks=shell.show_pre_run_checks,
        )
    return WorkspacePresentation(
        phase=WorkspacePhase.CHOOSE_DATA,
        collapse_input=shell.collapse_setup,
        show_pre_run_checks=shell.show_pre_run_checks,
    )
