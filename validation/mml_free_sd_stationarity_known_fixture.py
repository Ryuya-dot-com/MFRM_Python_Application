#!/usr/bin/env python3
"""Replay one known KAC case as excluded posthoc engineering evidence.

The selected run was inspected before this script and therefore cannot be used
to tune or qualify prospective stationarity thresholds.  Its sole purpose is
to preserve a reproducible regression case for the legacy-EM plateau followed
by finite-Q joint structural/log-SD polishing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd
import scipy


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import streamlit_app as app  # noqa: E402
from mfrm_app.mml_stationarity import (  # noqa: E402
    JointPolishOptions,
    audit_joint_gradient,
    make_joint_free_sd_functions,
)
from validation.mml_free_sd_stationarity_adapter import (  # noqa: E402
    prepare_app_free_sd_problem,
    run_app_free_sd_stationarity_v2,
)
from validation.operating_characteristics_python_mml import build_mml_kwargs  # noqa: E402


SCHEMA_VERSION = "mml-free-sd-stationarity-known-fixture-v2"
CLASSIFICATION = "POSTHOC_ENGINEERING_KNOWN_EXCLUDED_FROM_QUALIFICATION"
RUN_ID = "gamma_pos_0p0__vector-05"
DEFAULT_STUDY = ROOT / "validation" / "kac200_20260811"
DEFAULT_OUTPUT = ROOT / "validation" / "mml_free_sd_stationarity_known_fixture_v2_20260811"
PLAN_PATH = ROOT / "validation" / "mml_free_sd_stationarity_known_fixture_plan_20260811.json"
SOURCE_FILES = (
    ROOT / "streamlit_app.py",
    ROOT / "mfrm_app" / "mml_stationarity.py",
    ROOT / "mfrm_app" / "mml_engine_v2.py",
    ROOT / "validation" / "mml_free_sd_stationarity_adapter.py",
    ROOT / "validation" / "operating_characteristics_python_mml.py",
    ROOT / "validation" / "known_assignment_mml_crossfit.R",
    ROOT / "requirements.txt",
    ROOT / "requirements-dev.txt",
    PLAN_PATH,
    Path(__file__).resolve(),
)
EXPECTED_INPUT_SHA256 = {
    "study_identity.json": "a75869852519ca4947c25efc599bdb3a7c4d84532eea5623160c0ae51ad1a806",
    "retained_input/generated_ratings.csv": "6fb16fc5291a437062c8689e2cb098b4fc4ded30b765efbeffc9954642fec353",
    "retained_input/manifest.csv": "c19d44fd3fea390f2bd0565ba8a7678ce3a4280a78314b37fd21be31780833c0",
    "r_mml/crossfit_runs.csv": "a27f6c962812d844e1ddfad734e83033e9b5561d04961bf5cacec13974a6b009",
    "r_mml/verification_identity.json": "4509aea77577f7b73f052362274bf4443da6b130dbaa4658f1631e777aa10070",
}
EXPECTED_STRICT_AUDIT_SHA256 = {
    "audit_identity.json": "148ea5c5cbc72a2123f6f2007f4816c7a38faa5fe5035ad216cf88c974e1d116",
    "assessment.json": "25f8d091ac45ad65e9f4bbcec4be5955e5aef41bb12bfe2b77db2647f5e8fef9",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        return numeric if np.isfinite(numeric) else None
    return value


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(_json_value(value), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _load_fixture(study: Path) -> tuple[pd.DataFrame, pd.Series]:
    retained = study / "retained_input"
    ratings_path = retained / "generated_ratings.csv"
    manifest_path = retained / "manifest.csv"
    if not ratings_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError("frozen KAC retained input is unavailable")
    ratings = pd.read_csv(ratings_path)
    ratings = ratings.loc[ratings["RunId"].astype(str) == RUN_ID].copy()
    manifest = pd.read_csv(manifest_path)
    row = manifest.loc[manifest["RunId"].astype(str) == RUN_ID]
    if len(row) != 1 or len(ratings) != 960:
        raise ValueError("known fixture must resolve to one manifest row and 960 ratings")
    return ratings, row.iloc[0]


def validate_frozen_inputs(study: Path) -> dict[str, str]:
    observed: dict[str, str] = {}
    for name, expected in EXPECTED_INPUT_SHA256.items():
        path = study / name
        if not path.is_file():
            raise FileNotFoundError(f"expected frozen input is missing: {path}")
        actual = sha256_file(path)
        if actual != expected:
            raise ValueError(f"frozen input hash mismatch: {name}")
        observed[name] = actual
    strict_root = ROOT / "validation" / "kac200_posthoc_strict_audit_20260811"
    for name, expected in EXPECTED_STRICT_AUDIT_SHA256.items():
        path = strict_root / name
        if not path.is_file():
            raise FileNotFoundError(f"expected strict-audit evidence is missing: {path}")
        actual = sha256_file(path)
        if actual != expected:
            raise ValueError(f"strict-audit evidence hash mismatch: {name}")
        observed[f"strict_audit/{name}"] = actual
    return observed


def _fit_legacy(ratings: pd.DataFrame, manifest_row: pd.Series) -> dict[str, Any]:
    kwargs = build_mml_kwargs("PYTHON_MML_FREE_SD_Q31", int(manifest_row["Categories"]))
    kwargs.update(
        {
            "model": "PCM",
            "step_facet": "Criterion",
            "estimate_population_sd": True,
            "population_prior_sd": 0.8,
        }
    )
    return app.mfrm_estimate(
        data=ratings[["Person", "Rater", "Task", "Criterion", "Score"]].copy(),
        **kwargs,
    )


def replay(study: Path = DEFAULT_STUDY) -> dict[str, Any]:
    ratings, manifest_row = _load_fixture(study)
    legacy = _fit_legacy(ratings, manifest_row)
    problem = prepare_app_free_sd_problem(legacy, quadrature_points=31)
    start = np.concatenate(
        [np.asarray(problem.structural_start, dtype=float), [np.log(problem.sigma_start)]]
    )
    objective, value_gradient = make_joint_free_sd_functions(
        problem.value,
        problem.value_gradient,
        log_sigma_relative_step=1e-5,
    )
    legacy_gradient = audit_joint_gradient(
        objective,
        value_gradient,
        start,
        relative_step=1e-5,
    )
    _, run = run_app_free_sd_stationarity_v2(
        legacy,
        quadrature_points=31,
        options=JointPolishOptions(
            maxiter=500,
            gtol=1e-8,
            ftol=1e-15,
            maxls=50,
            log_sigma_relative_step=1e-5,
        ),
    )
    summary = legacy["summary"].iloc[0]
    r_crossfit_path = study / "r_mml" / "crossfit_runs.csv"
    r_rows = pd.read_csv(r_crossfit_path)
    r_row = r_rows.loc[r_rows["RunId"].astype(str) == RUN_ID]
    if len(r_row) != 1:
        raise ValueError("known fixture lacks its single retained R crossfit row")
    r = r_row.iloc[0]
    return {
        "schema_version": SCHEMA_VERSION,
        "classification": CLASSIFICATION,
        "qualification_eligible": False,
        "reason_excluded": (
            "Run and numerical behavior were known before this replay; values cannot "
            "select or justify prospective thresholds."
        ),
        "run_id": RUN_ID,
        "input_identity": {
            "study_identity_sha256": sha256_file(study / "study_identity.json"),
            "generated_ratings_sha256": sha256_file(
                study / "retained_input" / "generated_ratings.csv"
            ),
            "manifest_sha256": sha256_file(study / "retained_input" / "manifest.csv"),
            "r_crossfit_runs_sha256": sha256_file(r_crossfit_path),
            "selected_rows": int(len(ratings)),
        },
        "runtime": {
            "python": platform.python_version(),
            "python_executable_sha256": sha256_file(Path(sys.executable)),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "blas_lapack": np.__config__.CONFIG,
            "requirements_sha256": {
                "requirements.txt": sha256_file(ROOT / "requirements.txt"),
                "requirements-dev.txt": sha256_file(ROOT / "requirements-dev.txt"),
            },
            "app_version": str(legacy["config"].get("app_version", "")),
        },
        "legacy_em": {
            "converged": bool(summary.get("Converged", False)),
            "inference_ready_legacy_contract": bool(summary.get("InferenceReady", False)),
            "iterations": int(summary.get("Iterations", 0)),
            "log_likelihood": float(summary["LogLik"]),
            "estimated_population_sd": float(summary["EstimatedPopulationSD"]),
            "reported_structural_gradient_l2norm": float(summary["GradientNorm"]),
            "joint_gradient_supnorm": legacy_gradient.analytical_supnorm,
            "joint_gradient_fd_supnorm": legacy_gradient.finite_difference_supnorm,
            "joint_gradient_analytic_fd_max_abs_difference": (
                legacy_gradient.maximum_absolute_difference
            ),
        },
        "joint_stationarity_v2": run.to_dict(),
        "retained_independent_r": {
            "r_log_likelihood_q31": float(r["RLogLikQ31"]),
            "r_gradient_supnorm_q31": float(r["RGradientSupNormQ31"]),
            "r_sigma_q31": float(r["RSigma"]),
            "r_optimized_log_likelihood_q61": float(r["ROptimizedLogLikQ61"]),
            "r_gradient_supnorm_q61": float(r["RGradientSupNormQ61"]),
            "r_sigma_q61": float(r["RQ61Sigma"]),
        },
        "interpretation": {
            "allowed": (
                "Regression evidence that the legacy likelihood plateau and finite-Q "
                "joint stationarity are distinct numerical conditions."
            ),
            "forbidden": (
                "Prospective qualification, threshold selection, endpoint re-analysis, "
                "or revision of the frozen kac200 scientific conclusion."
            ),
        },
    }


def _markdown(payload: dict[str, Any]) -> str:
    legacy = payload["legacy_em"]
    run = payload["joint_stationarity_v2"]
    primary = run["primary_polish"]
    restart = run["restart_polish"]
    return f"""# Free-SD MML stationarity known fixture

Status: `{CLASSIFICATION}`. Qualification eligible: **No**.

This is a posthoc engineering regression case. The run and its behavior were
known before replay, so none of these values may define or justify prospective
acceptance thresholds.

## Observed numerical separation

- Legacy EM: log likelihood `{legacy['log_likelihood']:.15g}`, sigma
  `{legacy['estimated_population_sd']:.15g}`, joint score sup norm
  `{legacy['joint_gradient_supnorm']:.9g}`.
- First joint polish: objective improvement
  `{primary['objective_improvement']:.9g}`, maximum coordinate displacement
  `{primary['maximum_coordinate_displacement']:.9g}`.
- Restart polish: improvement per observation
  `{run['restart_improvement_per_observation']:.9g}`, displacement
  `{run['restart_displacement']:.9g}`.
- Final joint projected score sup norm
  `{run['final_projected_gradient_supnorm']:.9g}`; optimizer termination alone
  is not interpreted as scientific readiness.

Allowed use: regression testing of the numerical pathway. Forbidden use:
threshold selection, prospective qualification, endpoint re-analysis, or
revision of the frozen `kac200` scientific conclusion.
"""


def publish(output: Path, study: Path = DEFAULT_STUDY) -> dict[str, Any]:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    frozen_before = validate_frozen_inputs(study)
    source_before = {path.relative_to(ROOT).as_posix(): sha256_file(path) for path in SOURCE_FILES}
    payload = replay(study)
    if validate_frozen_inputs(study) != frozen_before:
        raise RuntimeError("frozen evidence changed during replay")
    source_after = {path.relative_to(ROOT).as_posix(): sha256_file(path) for path in SOURCE_FILES}
    if source_after != source_before:
        raise RuntimeError("source changed during replay")
    identity = {
        "schema_version": SCHEMA_VERSION,
        "classification": CLASSIFICATION,
        "source_sha256": source_before,
        "frozen_evidence_sha256": frozen_before,
    }
    output.mkdir(parents=True, exist_ok=False)
    stage_replay = output / f".replay.json.stage.{os.getpid()}"
    stage_readme = output / f".README.md.stage.{os.getpid()}"
    stage_identity = output / f".identity.json.stage.{os.getpid()}"
    committed = [output / "replay.json", output / "README.md", output / "identity.json"]
    try:
        _write_json(stage_replay, payload)
        stage_readme.write_text(_markdown(payload), encoding="utf-8")
        identity["artifact_sha256"] = {
            "replay.json": sha256_file(stage_replay),
            "README.md": sha256_file(stage_readme),
        }
        _write_json(stage_identity, identity)
        os.replace(stage_replay, committed[0])
        os.replace(stage_readme, committed[1])
        os.replace(stage_identity, committed[2])
        observed_files = sorted(path.name for path in output.iterdir())
        if observed_files != ["README.md", "identity.json", "replay.json"]:
            raise RuntimeError("published fixture file set differs from the contract")
        reread = json.loads(committed[2].read_text(encoding="utf-8"))
        for name, expected in reread["artifact_sha256"].items():
            if sha256_file(output / name) != expected:
                raise RuntimeError(f"published fixture hash mismatch: {name}")
        return payload
    except Exception:
        for path in (stage_replay, stage_readme, stage_identity, *reversed(committed)):
            try:
                path.unlink(missing_ok=True)
            except Exception:
                pass
        try:
            output.rmdir()
        except Exception:
            pass
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study", type=Path, default=DEFAULT_STUDY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = publish(args.output.resolve(), args.study.resolve())
    print(json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
