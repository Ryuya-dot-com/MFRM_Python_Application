from pathlib import Path

import pytest

from validation.facets_visible_launcher import _output_state, invoke_facets_batch_no
from validation.operating_characteristics_facets import (
    _facets_command,
    _facets_path_budget_violations,
    _facets_readiness_paths,
)


def test_output_state_requires_every_nonempty_file(tmp_path: Path):
    report = tmp_path / "report.txt"
    score = tmp_path / "score.txt"
    assert _output_state((report, score)) is None
    report.write_text("report", encoding="utf-8")
    assert _output_state((report, score)) is None
    score.write_text("", encoding="utf-8")
    assert _output_state((report, score)) is None
    score.write_text("score", encoding="utf-8")
    state = _output_state((report, score))
    assert state is not None
    assert [row[1] for row in state] == [6, 5]


def test_visible_launcher_requires_readiness_paths(tmp_path: Path):
    with pytest.raises(ValueError, match="readiness path"):
        invoke_facets_batch_no(
            ["unused"], cwd=tmp_path, readiness_paths=(), timeout_seconds=1
        )


def test_visible_launcher_rejects_nonpositive_timing(tmp_path: Path):
    with pytest.raises(ValueError, match="positive"):
        invoke_facets_batch_no(
            ["unused"],
            cwd=tmp_path,
            readiness_paths=(tmp_path / "report",),
            timeout_seconds=0,
        )


def test_facets_command_supports_explicit_visible_mode(tmp_path: Path):
    command = _facets_command(
        Path("Facets.exe"),
        tmp_path / "analysis.txt",
        tmp_path / "report.txt",
        ("Umean=0,1,2",),
        batch_value="NO",
    )
    assert command[1] == "BATCH=NO"
    assert command[-1] == "Umean=0,1,2"


def test_readiness_paths_include_all_overridden_scorefiles(tmp_path: Path):
    spec = tmp_path / "analysis.txt"
    spec.write_text(
        "Facets=4\nScorefile=original.txt\n", encoding="utf-8"
    )
    report = tmp_path / "report.txt"
    score_base = tmp_path / "override.txt"
    paths = _facets_readiness_paths(
        spec, report, (f"Scorefile={score_base}",)
    )
    assert paths == (
        report.resolve(),
        tmp_path / "override.1.txt",
        tmp_path / "override.2.txt",
        tmp_path / "override.3.txt",
        tmp_path / "override.4.txt",
    )


def test_facets_path_budget_reports_only_overlong_paths(tmp_path: Path):
    short = tmp_path / "short.txt"
    long = tmp_path / ("x" * 80) / "score.4.txt"
    budget = len(str(short.resolve()))
    assert _facets_path_budget_violations(
        (short, long), maximum_characters=budget
    ) == ((long.resolve(), len(str(long.resolve()))),)


def test_facets_path_budget_rejects_nonpositive_limit(tmp_path: Path):
    with pytest.raises(ValueError, match="positive"):
        _facets_path_budget_violations((tmp_path,), maximum_characters=0)
