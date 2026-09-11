#!/usr/bin/env python3
"""Create frozen fixed-calibration Warm-WLE parity fixtures and Python output."""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

import sys

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from mfrm_app.person_scoring import score_fixed_calibration_persons


@dataclass(frozen=True)
class FixtureCase:
    name: str
    intercepts: np.ndarray
    slopes: np.ndarray
    responses: np.ndarray


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixture_cases() -> list[FixtureCase]:
    binary_intercepts = np.zeros((5, 2), dtype=float)
    binary_slopes = np.tile(np.array([0.0, 1.0]), (5, 1))
    binary_responses = np.array(
        [
            [0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0],
            [1, np.nan, 0, 1, 0],
            [1, 1, 1, 1, 1],
        ],
        dtype=float,
    )

    rsm_difficulty = np.array([-1.1, -0.6, -0.1, 0.3, 0.8, 1.2])
    rsm_step_cumulative = np.array([0.0, -0.8, -0.7, 0.2])
    categories4 = np.arange(4, dtype=float)
    rsm_intercepts = -(
        rsm_difficulty[:, None] * categories4[None, :]
        + rsm_step_cumulative[None, :]
    )
    rsm_slopes = np.tile(categories4, (len(rsm_difficulty), 1))
    rsm_responses = np.array(
        [
            [0, 0, 0, 0, 0, 0],
            [0, 1, 2, 3, np.nan, 2],
            [3, 3, 2, 1, 0, np.nan],
            [3, 3, 3, 3, 3, 3],
        ],
        dtype=float,
    )

    pcm_steps = np.array(
        [
            [-1.2, -0.2, 0.7],
            [-0.8, 0.3, 1.0],
            [-0.4, -0.1, 0.6],
            [-1.0, 0.5, 1.3],
            [-0.6, 0.1, 0.9],
        ]
    )
    pcm_intercepts = np.column_stack(
        [np.zeros(len(pcm_steps)), -np.cumsum(pcm_steps, axis=1)]
    )
    pcm_slopes = np.tile(categories4, (len(pcm_steps), 1))
    pcm_responses = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 1, 2, 3, np.nan],
            [3, 2, np.nan, 1, 0],
            [3, 3, 3, 3, 3],
        ],
        dtype=float,
    )

    gpcm_steps = np.array(
        [
            [-1.0, -0.1, 0.8],
            [-0.7, 0.2, 1.1],
            [-1.3, 0.0, 0.5],
            [-0.5, 0.4, 1.4],
            [-0.9, -0.2, 0.9],
        ]
    )
    gpcm_discrimination = np.array([0.65, 0.9, 1.15, 1.4, 1.8])
    gpcm_intercepts = np.column_stack(
        [
            np.zeros(len(gpcm_steps)),
            -gpcm_discrimination[:, None] * np.cumsum(gpcm_steps, axis=1),
        ]
    )
    gpcm_slopes = gpcm_discrimination[:, None] * categories4[None, :]
    gpcm_responses = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 1, 2, 3, np.nan],
            [3, np.nan, 2, 1, 0],
            [3, 3, 3, 3, 3],
        ],
        dtype=float,
    )

    return [
        FixtureCase(
            "binary_rasch_extremes", binary_intercepts, binary_slopes, binary_responses
        ),
        FixtureCase(
            "rsm_four_categories_missing", rsm_intercepts, rsm_slopes, rsm_responses
        ),
        FixtureCase(
            "pcm_item_specific_intercepts_missing",
            pcm_intercepts,
            pcm_slopes,
            pcm_responses,
        ),
        FixtureCase(
            "gpcm_varying_discrimination_missing",
            gpcm_intercepts,
            gpcm_slopes,
            gpcm_responses,
        ),
    ]


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    calibration_rows: list[dict[str, object]] = []
    response_rows: list[dict[str, object]] = []
    python_rows: list[pd.DataFrame] = []
    for case in fixture_cases():
        n_items, n_categories = case.intercepts.shape
        if case.slopes.shape != case.intercepts.shape:
            raise AssertionError(f"Mismatched calibration arrays: {case.name}")
        if case.responses.shape[1] != n_items:
            raise AssertionError(f"Mismatched response items: {case.name}")
        for item in range(n_items):
            for category in range(n_categories):
                calibration_rows.append(
                    {
                        "Case": case.name,
                        "Item": f"I{item + 1}",
                        "Category": category,
                        "Intercept": case.intercepts[item, category],
                        "Slope": case.slopes[item, category],
                    }
                )

        persons: list[str] = []
        observed: list[int] = []
        intercepts: list[np.ndarray] = []
        slopes: list[np.ndarray] = []
        for person_index, response_row in enumerate(case.responses):
            person = f"P{person_index + 1}"
            for item, response in enumerate(response_row):
                response_rows.append(
                    {
                        "Case": case.name,
                        "Person": person,
                        "Item": f"I{item + 1}",
                        "Observed": response,
                    }
                )
                if np.isfinite(response):
                    persons.append(person)
                    observed.append(int(response))
                    intercepts.append(case.intercepts[item])
                    slopes.append(case.slopes[item])
        scored = score_fixed_calibration_persons(
            persons,
            observed,
            np.asarray(intercepts),
            np.asarray(slopes),
            person_levels=[f"P{index + 1}" for index in range(len(case.responses))],
        )
        scored.insert(0, "Case", case.name)
        python_rows.append(scored)

    return (
        pd.DataFrame(calibration_rows),
        pd.DataFrame(response_rows),
        pd.concat(python_rows, ignore_index=True),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output_dir = args.output.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    calibration, responses, python_output = build_tables()
    paths = {
        "calibration": output_dir / "fixture_calibration.csv",
        "responses": output_dir / "fixture_responses.csv",
        "python_output": output_dir / "python_wle.csv",
    }
    calibration.to_csv(paths["calibration"], index=False, lineterminator="\n", float_format="%.17g")
    responses.to_csv(paths["responses"], index=False, lineterminator="\n", float_format="%.17g")
    python_output.to_csv(paths["python_output"], index=False, lineterminator="\n", float_format="%.17g")

    tracked_inputs = {
        "plan": ROOT / "validation/fixed_calibration_wle_plan_20260809.json",
        "amendment": ROOT / "validation/fixed_calibration_wle_plan_amendment_20260809.json",
        "python_core": ROOT / "mfrm_app/person_scoring.py",
        "fixture_generator": Path(__file__).resolve(),
        **paths,
    }
    manifest = {
        "schema_version": "mfrm-fixed-calibration-wle-fixture-v1",
        "case_names": [case.name for case in fixture_cases()],
        "persons": int(python_output.shape[0]),
        "all_python_status_ok": bool(python_output["Status"].eq("ok").all()),
        "maximum_python_adjusted_score_residual": float(
            python_output["AdjustedScoreResidual"].max()
        ),
        "files": {
            name: {"path": str(path.relative_to(ROOT)), "sha256": sha256_file(path)}
            for name, path in tracked_inputs.items()
        },
    }
    (output_dir / "fixture_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
