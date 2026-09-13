#!/usr/bin/env python3
"""Cross-engine RSM/PCM stress pilot for the standalone Python MFRM app.

The script has three deliberately separate stages:

1. ``generate`` creates deterministic long-format rating data and truth tables.
2. ``fit-python`` fits only the standalone Python engine.
3. ``summarize`` combines the Python and R-engine CSV outputs and aligns every
   structural estimate on a common cumulative-difficulty surface.

The companion ``cross_engine_stress.R`` fits mfrmr, TAM, immer, and sirt to the
same generated CSV rows.  This is a deterministic stress pilot, not a
high-replication coverage or estimator-ranking study.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


SCHEMA_VERSION = "cross-engine-stress-v1"
N_CAT = 4


@dataclass(frozen=True)
class Cell:
    model: str
    scenario: str
    replicate: int
    run_jml: bool
    run_mml: bool
    run_sirt: bool
    run_q61: bool = False


def cells() -> list[Cell]:
    out: list[Cell] = []
    for model in ("RSM", "PCM"):
        for rep in (1, 2):
            out.append(Cell(model, "balanced", rep, True, True, model == "PCM"))
            out.append(Cell(model, "sparse_connected", rep, True, True, model == "PCM"))
        out.append(Cell(model, "forced_extremes", 1, True, False, False))
        out.append(Cell(model, "disconnected", 1, True, True, model == "PCM"))
    for model in ("RSM", "PCM"):
        for rep in (1, 2):
            out.append(Cell(model, "high_variance", rep, False, True, model == "PCM", True))
    out.extend(
        [
            Cell("PCM", "rare_categories", 1, True, True, True, True),
            Cell("PCM", "local_dependence", 1, False, True, True),
            Cell("PCM", "score_mnar", 1, False, True, True),
            Cell("PCM", "small_sample", 1, True, True, True),
        ]
    )
    return out


def scenario_settings(name: str) -> dict[str, Any]:
    base = {
        "n_person": 80,
        "n_rater": 4,
        "n_criterion": 4,
        "theta_sd": 1.0,
        "assignment": "crossed",
        "forced_extreme_fraction": 0.0,
        "local_dependence_sd": 0.0,
        "missing_mechanism": "none",
        "rare_steps": False,
        "expected_connected": True,
    }
    changes = {
        "balanced": {},
        "sparse_connected": {"n_rater": 5, "assignment": "cycle_two"},
        "forced_extremes": {"forced_extreme_fraction": 0.20},
        "disconnected": {
            "n_person": 60,
            "assignment": "two_components",
            "expected_connected": False,
        },
        "high_variance": {"n_person": 100, "theta_sd": 2.2},
        "rare_categories": {"n_person": 100, "rare_steps": True},
        "local_dependence": {"n_person": 100, "local_dependence_sd": 0.9},
        "score_mnar": {"n_person": 100, "missing_mechanism": "score_mnar"},
        "small_sample": {"n_person": 28, "n_rater": 3, "n_criterion": 2},
    }
    if name not in changes:
        raise KeyError(name)
    base.update(changes[name])
    return base


def dataset_seed(cell: Cell) -> int:
    scenario_index = sorted({c.scenario for c in cells()}).index(cell.scenario) + 1
    model_offset = 0 if cell.model == "RSM" else 100_000
    return 2_026_080_900 + model_offset + scenario_index * 1_000 + cell.replicate * 17


def centered_sequence(n: int, scale: float) -> np.ndarray:
    values = np.linspace(-scale, scale, n, dtype=float)
    return values - values.mean()


def step_matrix(model: str, n_criterion: int, rare: bool) -> np.ndarray:
    base = np.array([-4.0, 0.0, 4.0] if rare else [-1.1, 0.0, 1.1], dtype=float)
    if model == "RSM":
        return base[None, :]
    rows: list[np.ndarray] = []
    for i in range(n_criterion):
        perturb = np.array(
            [-0.18 * math.cos(i + 0.3), 0.15 * math.sin(i + 0.7), 0.12 * math.cos(i + 1.1)]
        )
        row = base + perturb
        rows.append(row - row.mean())
    return np.vstack(rows)


def assignment_pairs(n_person: int, n_rater: int, assignment: str) -> list[tuple[int, int]]:
    if assignment == "crossed":
        return [(p, r) for p in range(n_person) for r in range(n_rater)]
    if assignment == "cycle_two":
        return [(p, r) for p in range(n_person) for r in (p % n_rater, (p + 1) % n_rater)]
    if assignment == "two_components":
        half = n_person // 2
        return [
            (p, r)
            for p in range(n_person)
            for r in ((0, 1) if p < half else (2, 3))
        ]
    raise ValueError(f"Unknown assignment: {assignment}")


def softmax_sample(logits: np.ndarray, rng: np.random.Generator) -> int:
    logits = logits - np.max(logits)
    probs = np.exp(logits)
    probs /= probs.sum()
    return int(rng.choice(len(probs), p=probs))


def rater_graph_components(data: pd.DataFrame) -> int:
    levels = sorted(data["Rater"].astype(str).unique())
    neighbors = {level: set() for level in levels}
    for _, sub in data[["Person", "Rater"]].drop_duplicates().groupby("Person"):
        raters = sorted(sub["Rater"].astype(str).unique())
        for left in raters:
            neighbors[left].update(r for r in raters if r != left)
    remaining = set(levels)
    components = 0
    while remaining:
        components += 1
        stack = [remaining.pop()]
        while stack:
            found = neighbors[stack.pop()] & remaining
            remaining -= found
            stack.extend(found)
    return components


def generate_one(cell: Cell) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    settings = scenario_settings(cell.scenario)
    seed = dataset_seed(cell)
    rng = np.random.default_rng(seed)
    n_person = int(settings["n_person"])
    n_rater = int(settings["n_rater"])
    n_criterion = int(settings["n_criterion"])
    theta = rng.normal(0.0, float(settings["theta_sd"]), n_person)
    rater = centered_sequence(n_rater, 0.65)
    criterion = centered_sequence(n_criterion, 0.45)
    steps = step_matrix(cell.model, n_criterion, bool(settings["rare_steps"]))
    rows: list[dict[str, Any]] = []
    for p, r in assignment_pairs(n_person, n_rater, str(settings["assignment"])):
        cluster = rng.normal(0.0, float(settings["local_dependence_sd"]))
        for c in range(n_criterion):
            step = steps[0 if cell.model == "RSM" else c]
            cumulative = np.concatenate(([0.0], np.cumsum(step)))
            k = np.arange(N_CAT, dtype=float)
            eta = theta[p] + cluster - rater[r] - criterion[c]
            score = softmax_sample(k * eta - cumulative, rng)
            rows.append(
                {
                    "Person": f"P{p + 1:03d}",
                    "Rater": f"R{r + 1:02d}",
                    "Criterion": f"C{c + 1:02d}",
                    "Score": score,
                }
            )
    data = pd.DataFrame(rows)
    forced_n = int(math.floor(n_person * float(settings["forced_extreme_fraction"]) / 2.0))
    if forced_n:
        low = {f"P{i + 1:03d}" for i in range(forced_n)}
        high = {f"P{i + 1 + forced_n:03d}" for i in range(forced_n)}
        data.loc[data["Person"].isin(low), "Score"] = 0
        data.loc[data["Person"].isin(high), "Score"] = N_CAT - 1
    before_missing = len(data)
    if settings["missing_mechanism"] == "score_mnar":
        probability = np.where(data["Score"].to_numpy() >= 2, 0.42, 0.06)
        keep = rng.random(len(data)) >= probability
        data = data.loc[keep].reset_index(drop=True)
    dataset_id = f"{cell.model}-{cell.scenario}-r{cell.replicate:02d}"
    data.insert(0, "DatasetId", dataset_id)

    truth_rows: list[dict[str, Any]] = []
    for r in range(n_rater):
        for c in range(n_criterion):
            step = steps[0 if cell.model == "RSM" else c]
            cumulative = np.cumsum(step)
            for category in range(1, N_CAT):
                truth_rows.append(
                    {
                        "DatasetId": dataset_id,
                        "Rater": f"R{r + 1:02d}",
                        "Criterion": f"C{c + 1:02d}",
                        "Category": category,
                        "Truth": category * (rater[r] + criterion[c]) + cumulative[category - 1],
                    }
                )
    truth_surface = pd.DataFrame(truth_rows)
    person_truth = pd.DataFrame(
        {
            "DatasetId": dataset_id,
            "Person": [f"P{i + 1:03d}" for i in range(n_person)],
            "ThetaTruth": theta,
        }
    )
    grouped = data.groupby("Person")["Score"].agg(["sum", "count"])
    observed_extreme = (grouped["sum"] == 0) | (grouped["sum"] == grouped["count"] * (N_CAT - 1))
    person_truth["ObservedExtreme"] = person_truth["Person"].map(observed_extreme).fillna(True).astype(bool)
    person_truth["ObservedRatings"] = person_truth["Person"].map(grouped["count"]).fillna(0).astype(int)
    support_index = pd.MultiIndex.from_product(
        [
            sorted(data["Rater"].astype(str).unique()),
            sorted(data["Criterion"].astype(str).unique()),
            range(N_CAT),
        ],
        names=["Rater", "Criterion", "Score"],
    )
    support_counts = data.groupby(["Rater", "Criterion", "Score"]).size().reindex(support_index, fill_value=0)
    top_counts = support_counts.xs(N_CAT - 1, level="Score")
    criterion_support = data.groupby(["Criterion", "Score"]).size().unstack(fill_value=0)

    meta = {
        "SchemaVersion": SCHEMA_VERSION,
        "DatasetId": dataset_id,
        "Model": cell.model,
        "Scenario": cell.scenario,
        "Replicate": cell.replicate,
        "Seed": seed,
        "Persons": n_person,
        "Raters": n_rater,
        "Criteria": n_criterion,
        "Categories": N_CAT,
        "ThetaSDTruth": float(settings["theta_sd"]),
        "Assignment": settings["assignment"],
        "ForcedExtremeFraction": float(settings["forced_extreme_fraction"]),
        "LocalDependenceSD": float(settings["local_dependence_sd"]),
        "MissingMechanism": settings["missing_mechanism"],
        "PlannedRatings": before_missing,
        "ObservedRatings": len(data),
        "ObservedFraction": len(data) / before_missing,
        "ObservedExtremePersons": int(person_truth["ObservedExtreme"].sum()),
        "RaterGraphComponents": rater_graph_components(data),
        "MinimumPseudoitemCategoryCount": int(support_counts.min()),
        "PseudoitemCategoryZeroCells": int((support_counts == 0).sum()),
        "PseudoitemsMissingTopCategory": int((top_counts == 0).sum()),
        "CriterionCategorySupportComplete": bool(
            criterion_support.reindex(columns=range(N_CAT), fill_value=0).gt(0).all().all()
        ),
        "ExpectedConnected": bool(settings["expected_connected"]),
        "RunJML": cell.run_jml,
        "RunMML": cell.run_mml,
        "RunSirt": cell.run_sirt,
        "RunQ61": cell.run_q61,
    }
    return data, truth_surface, person_truth, meta


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_generation(output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    data_frames: list[pd.DataFrame] = []
    surfaces: list[pd.DataFrame] = []
    persons: list[pd.DataFrame] = []
    manifest: list[dict[str, Any]] = []
    for cell in cells():
        data, surface, person, meta = generate_one(cell)
        data_frames.append(data)
        surfaces.append(surface)
        persons.append(person)
        manifest.append(meta)
    pd.concat(data_frames, ignore_index=True).to_csv(output / "simulated_ratings.csv", index=False)
    pd.concat(surfaces, ignore_index=True).to_csv(output / "truth_surface.csv", index=False)
    pd.concat(persons, ignore_index=True).to_csv(output / "truth_persons.csv", index=False)
    pd.DataFrame(manifest).to_csv(output / "manifest.csv", index=False)
    all_data = pd.concat(data_frames, ignore_index=True)
    support = (
        all_data.groupby(["DatasetId", "Rater", "Criterion", "Score"])
        .size()
        .rename("N")
        .reset_index()
    )
    support.to_csv(output / "observed_category_support.csv", index=False)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256_file(Path(__file__).resolve()),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "python": sys.version.replace("\n", " "),
    }
    (output / "generation_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def as_number(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if np.isfinite(number) else float("nan")


def python_surface(result: dict[str, Any], dataset_id: str, mode: str) -> pd.DataFrame:
    facets = result["facets"]["others"].copy()
    rater = facets.loc[facets["Facet"] == "Rater"].set_index("Level")["Estimate"].astype(float)
    criterion = facets.loc[facets["Facet"] == "Criterion"].set_index("Level")["Estimate"].astype(float)
    steps = result["steps"].copy()
    model = str(result["config"]["model"])
    rows: list[dict[str, Any]] = []
    for rr, rvalue in rater.items():
        for cc, cvalue in criterion.items():
            if model == "RSM":
                step = steps["Estimate"].astype(float).to_numpy()
            else:
                step = steps.loc[steps["StepFacet"].astype(str) == str(cc), "Estimate"].astype(float).to_numpy()
            for category, cumulative in enumerate(np.cumsum(step), start=1):
                rows.append(
                    {
                        "DatasetId": dataset_id,
                        "Engine": "PythonApp",
                        "Mode": mode,
                        "Rater": str(rr),
                        "Criterion": str(cc),
                        "Category": category,
                        "Estimate": category * (float(rvalue) + float(cvalue)) + float(cumulative),
                    }
                )
    return pd.DataFrame(rows)


def python_persons(result: dict[str, Any], dataset_id: str, mode: str) -> pd.DataFrame:
    table = result["facets"]["person"].copy()
    return pd.DataFrame(
        {
            "DatasetId": dataset_id,
            "Engine": "PythonApp",
            "Mode": mode,
            "Person": table["Person"].astype(str),
            "Estimate": pd.to_numeric(table["Estimate"], errors="coerce"),
            "ParameterStatus": "not_reported",
        }
    )


def app_git_identity(repo: Path) -> dict[str, str]:
    def run(*args: str) -> str:
        return subprocess.check_output(args, cwd=repo, text=True).strip()

    try:
        return {
            "git_head": run("git", "rev-parse", "HEAD"),
            "git_branch": run("git", "branch", "--show-current"),
            "git_status": run("git", "status", "--short"),
        }
    except Exception as exc:  # pragma: no cover - provenance fallback
        return {"git_error": str(exc)}


def fit_python(input_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[1]
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    import streamlit_app as app  # pylint: disable=import-outside-toplevel

    manifest = pd.read_csv(input_dir / "manifest.csv")
    ratings = pd.read_csv(input_dir / "simulated_ratings.csv")
    runs: list[dict[str, Any]] = []
    surfaces: list[pd.DataFrame] = []
    persons: list[pd.DataFrame] = []

    def execute(row: pd.Series, data: pd.DataFrame, mode: str, estimator: str, **kwargs: Any) -> None:
        started = time.perf_counter()
        caught: list[str] = []
        try:
            with warnings.catch_warnings(record=True) as warning_records:
                warnings.simplefilter("always")
                result = app.mfrm_estimate(
                    data=data,
                    person_col="Person",
                    facet_cols=["Rater", "Criterion"],
                    score_col="Score",
                    rating_min=0,
                    rating_max=N_CAT - 1,
                    model=str(row.Model),
                    step_facet="Criterion" if str(row.Model) == "PCM" else None,
                    min_obs_per_element=1,
                    min_obs_per_category=1,
                    **kwargs,
                )
                caught = [str(item.message) for item in warning_records]
            summary = result["summary"].iloc[0]
            config = result["config"]
            elapsed = time.perf_counter() - started
            runs.append(
                {
                    "DatasetId": row.DatasetId,
                    "Engine": "PythonApp",
                    "Mode": mode,
                    "Estimator": estimator,
                    "Model": row.Model,
                    "Scenario": row.Scenario,
                    "Replicate": int(row.Replicate),
                    "Quadrature": config.get("quad_points"),
                    "PopulationSDMode": "estimated" if config.get("estimate_population_sd") else ("fixed" if estimator == "MML" else "not_applicable"),
                    "PopulationSD": config.get("estimated_population_sd") if config.get("estimate_population_sd") else (config.get("population_prior_sd") if estimator == "MML" else np.nan),
                    "FitReturned": True,
                    "Converged": bool(summary.get("Converged", False)),
                    "ConvergenceBasis": "scipy_optimizer_success",
                    "InferenceReady": np.nan,
                    "LogLik": as_number(summary.get("LogLik")),
                    "Npar": as_number(summary.get("KParams")),
                    "Iterations": as_number(summary.get("Iterations")),
                    "GradientNorm": as_number(summary.get("GradientNorm")),
                    "ElapsedSeconds": elapsed,
                    "Error": "",
                    "Warnings": " | ".join(caught),
                    "ExpectedConnected": bool(row.ExpectedConnected),
                }
            )
            surfaces.append(python_surface(result, str(row.DatasetId), mode))
            persons.append(python_persons(result, str(row.DatasetId), mode))
        except Exception as exc:  # pylint: disable=broad-except
            runs.append(
                {
                    "DatasetId": row.DatasetId,
                    "Engine": "PythonApp",
                    "Mode": mode,
                    "Estimator": estimator,
                    "Model": row.Model,
                    "Scenario": row.Scenario,
                    "Replicate": int(row.Replicate),
                    "Quadrature": kwargs.get("quad_points"),
                    "PopulationSDMode": "unknown",
                    "PopulationSD": np.nan,
                    "FitReturned": False,
                    "Converged": False,
                    "ConvergenceBasis": "exception",
                    "InferenceReady": False,
                    "LogLik": np.nan,
                    "Npar": np.nan,
                    "Iterations": np.nan,
                    "GradientNorm": np.nan,
                    "ElapsedSeconds": time.perf_counter() - started,
                    "Error": f"{type(exc).__name__}: {exc}",
                    "Warnings": " | ".join(caught),
                    "ExpectedConnected": bool(row.ExpectedConnected),
                }
            )

    for row in manifest.itertuples(index=False):
        data = ratings.loc[ratings["DatasetId"] == row.DatasetId, ["Person", "Rater", "Criterion", "Score"]].copy()
        if bool(row.RunJML):
            execute(
                row,
                data,
                "APP_JML_RAW",
                "JML",
                method="JMLE",
                noncenter_facet="Person",
                maxit=500,
                reltol=1e-8,
            )
        if bool(row.RunMML):
            execute(
                row,
                data,
                "APP_MML_FREE_Q31",
                "MML",
                method="MML",
                noncenter_facet="Criterion",
                mml_engine="EM",
                quad_points=31,
                estimate_population_sd=True,
                population_prior_sd=1.0,
                maxit=350,
                reltol=1e-7,
            )
            if bool(row.RunQ61):
                execute(
                    row,
                    data,
                    "APP_MML_FREE_Q61",
                    "MML",
                    method="MML",
                    noncenter_facet="Criterion",
                    mml_engine="EM",
                    quad_points=61,
                    estimate_population_sd=True,
                    population_prior_sd=1.0,
                    maxit=350,
                    reltol=1e-7,
                )
            if str(row.Scenario) == "balanced":
                execute(
                    row,
                    data,
                    "APP_MML_FIXED_Q15_PERSON",
                    "MML",
                    method="MML",
                    noncenter_facet="Person",
                    mml_engine="EM",
                    quad_points=15,
                    estimate_population_sd=False,
                    population_prior_sd=1.0,
                    maxit=300,
                    reltol=1e-7,
                )
                execute(
                    row,
                    data,
                    "APP_MML_FIXED_Q31_PERSON",
                    "MML",
                    method="MML",
                    noncenter_facet="Person",
                    mml_engine="EM",
                    quad_points=31,
                    estimate_population_sd=False,
                    population_prior_sd=1.0,
                    maxit=300,
                    reltol=1e-7,
                )

    pd.DataFrame(runs).to_csv(output_dir / "fit_runs_python.csv", index=False)
    pd.concat(surfaces, ignore_index=True).to_csv(output_dir / "surfaces_python.csv", index=False)
    pd.concat(persons, ignore_index=True).to_csv(output_dir / "persons_python.csv", index=False)
    identity = {
        "schema_version": SCHEMA_VERSION,
        "app_file": str((repo / "streamlit_app.py").resolve()),
        "app_sha256": sha256_file(repo / "streamlit_app.py"),
        **app_git_identity(repo),
    }
    (output_dir / "python_runtime_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def correlation(left: pd.Series, right: pd.Series) -> float:
    x = pd.to_numeric(left, errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(right, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2 or np.std(x[ok]) == 0 or np.std(y[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(x[ok], y[ok])[0, 1])


def combine_csv(paths: Iterable[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in paths if path.exists()]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def summarize(input_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(input_dir / "manifest.csv")
    truth_surface = pd.read_csv(input_dir / "truth_surface.csv")
    truth_person = pd.read_csv(input_dir / "truth_persons.csv")
    runs = combine_csv([output_dir / "fit_runs_python.csv", output_dir / "fit_runs_r.csv"])
    surfaces = combine_csv([output_dir / "surfaces_python.csv", output_dir / "surfaces_r.csv"])
    persons = combine_csv([output_dir / "persons_python.csv", output_dir / "persons_r.csv"])
    if runs.empty or surfaces.empty:
        raise RuntimeError("Both Python and R engine outputs are required before summarize.")
    run_keys = ["DatasetId", "Engine", "Mode"]
    if runs.duplicated(run_keys).any():
        duplicate = runs.loc[runs.duplicated(run_keys, keep=False), run_keys]
        raise RuntimeError(f"Duplicate run identities detected: {duplicate.to_dict('records')[:5]}")
    returned_keys = set(map(tuple, runs.loc[runs["FitReturned"].astype(bool), run_keys].to_numpy()))
    surface_keys = set(map(tuple, surfaces[run_keys].drop_duplicates().to_numpy()))
    person_keys = set(map(tuple, persons[run_keys].drop_duplicates().to_numpy()))
    missing_surface = sorted(returned_keys - surface_keys)
    missing_person = sorted(returned_keys - person_keys)
    if missing_surface or missing_person:
        raise RuntimeError(
            f"Incomplete returned fits: missing surfaces={missing_surface[:5]}, "
            f"missing persons={missing_person[:5]}"
        )
    expected_surface_n = truth_surface.groupby("DatasetId").size().to_dict()
    actual_surface_n = surfaces.groupby(run_keys).size()
    wrong_surface_n = [
        (identity, int(count), int(expected_surface_n[identity[0]]))
        for identity, count in actual_surface_n.items()
        if int(count) != int(expected_surface_n[identity[0]])
    ]
    if wrong_surface_n:
        raise RuntimeError(f"Incomplete surface grids: {wrong_surface_n[:5]}")

    surface_metrics: list[dict[str, Any]] = []
    aligned_frames: list[pd.DataFrame] = []
    keys = ["DatasetId", "Rater", "Criterion", "Category"]
    for (dataset_id, engine, mode), estimate in surfaces.groupby(["DatasetId", "Engine", "Mode"], sort=False):
        joined = estimate.merge(truth_surface, on=keys, how="inner", validate="one_to_one")
        expected_n = len(truth_surface.loc[truth_surface["DatasetId"] == dataset_id])
        if len(joined) != expected_n:
            continue
        category = joined["Category"].to_numpy(dtype=float)
        delta = joined["Estimate"].to_numpy(dtype=float) - joined["Truth"].to_numpy(dtype=float)
        finite = np.isfinite(category) & np.isfinite(delta)
        if not finite.any():
            continue
        shift = float(np.sum(category[finite] * delta[finite]) / np.sum(category[finite] ** 2))
        joined["LocationShift"] = shift
        joined["EstimateAligned"] = joined["Estimate"] - joined["Category"] * shift
        joined["ErrorAligned"] = joined["EstimateAligned"] - joined["Truth"]
        joined["Engine"] = engine
        joined["Mode"] = mode
        aligned_frames.append(joined)
        err = joined["ErrorAligned"].to_numpy(dtype=float)
        surface_metrics.append(
            {
                "DatasetId": dataset_id,
                "Engine": engine,
                "Mode": mode,
                "SurfaceN": len(joined),
                "LocationShift": shift,
                "SurfaceBias": float(np.mean(err)),
                "SurfaceRMSE": float(np.sqrt(np.mean(err**2))),
                "SurfaceMAE": float(np.mean(np.abs(err))),
                "SurfaceMaxAbs": float(np.max(np.abs(err))),
                "SurfaceCorrelation": correlation(joined["EstimateAligned"], joined["Truth"]),
            }
        )
    surface_metrics_df = pd.DataFrame(surface_metrics)
    aligned = pd.concat(aligned_frames, ignore_index=True)

    enriched = runs.merge(surface_metrics_df, on=run_keys, how="left")
    enriched = enriched.merge(
        manifest[["DatasetId", "ThetaSDTruth", "ObservedExtremePersons", "RaterGraphComponents"]],
        on="DatasetId",
        how="left",
    )

    person_metrics: list[dict[str, Any]] = []
    shift_lookup = surface_metrics_df.set_index(run_keys)["LocationShift"].to_dict()
    estimator_lookup = runs.set_index(run_keys)["Estimator"].to_dict()
    for (dataset_id, engine, mode), estimate in persons.groupby(run_keys, sort=False):
        joined = estimate.merge(truth_person, on=["DatasetId", "Person"], how="inner")
        estimator = estimator_lookup.get((dataset_id, engine, mode), "unknown")
        if estimator == "JML":
            joined = joined.loc[~joined["ObservedExtreme"].astype(bool)].copy()
        shift = shift_lookup.get((dataset_id, engine, mode), np.nan)
        joined["EstimateAligned"] = pd.to_numeric(joined["Estimate"], errors="coerce") - shift
        joined["ErrorAligned"] = joined["EstimateAligned"] - joined["ThetaTruth"]
        finite = np.isfinite(joined["ErrorAligned"].to_numpy(dtype=float))
        joined = joined.loc[finite]
        if joined.empty:
            continue
        err = joined["ErrorAligned"].to_numpy(dtype=float)
        person_metrics.append(
            {
                "DatasetId": dataset_id,
                "Engine": engine,
                "Mode": mode,
                "PersonN": len(joined),
                "PersonBias": float(np.mean(err)),
                "PersonRMSE": float(np.sqrt(np.mean(err**2))),
                "PersonCorrelation": correlation(joined["EstimateAligned"], joined["ThetaTruth"]),
            }
        )
    person_metrics_df = pd.DataFrame(person_metrics)
    enriched = enriched.merge(person_metrics_df, on=run_keys, how="left")

    summary = (
        enriched.groupby(["Estimator", "Model", "Engine", "Mode"], dropna=False)
        .agg(
            Attempts=("DatasetId", "size"),
            FitsReturned=("FitReturned", "sum"),
            Converged=("Converged", "sum"),
            MedianSurfaceRMSE=("SurfaceRMSE", "median"),
            MaxSurfaceRMSE=("SurfaceRMSE", "max"),
            MedianSurfaceCorrelation=("SurfaceCorrelation", "median"),
            MedianPersonRMSE=("PersonRMSE", "median"),
            MedianPersonCorrelation=("PersonCorrelation", "median"),
            MedianSeconds=("ElapsedSeconds", "median"),
        )
        .reset_index()
    )
    scenario_summary = (
        enriched.groupby(["Estimator", "Model", "Scenario", "Engine", "Mode"], dropna=False)
        .agg(
            Attempts=("DatasetId", "size"),
            FitsReturned=("FitReturned", "sum"),
            Converged=("Converged", "sum"),
            MedianSurfaceRMSE=("SurfaceRMSE", "median"),
            MaxSurfaceRMSE=("SurfaceRMSE", "max"),
            MedianPersonRMSE=("PersonRMSE", "median"),
            MedianPersonCorrelation=("PersonCorrelation", "median"),
        )
        .reset_index()
    )
    connected_summary = (
        enriched.loc[enriched["ExpectedConnected"].astype(bool)]
        .groupby(["Estimator", "Model", "Engine", "Mode"], dropna=False)
        .agg(
            Attempts=("DatasetId", "size"),
            FitsReturned=("FitReturned", "sum"),
            Converged=("Converged", "sum"),
            MedianSurfaceRMSE=("SurfaceRMSE", "median"),
            MaxSurfaceRMSE=("SurfaceRMSE", "max"),
            MedianPersonRMSE=("PersonRMSE", "median"),
            MedianPersonCorrelation=("PersonCorrelation", "median"),
        )
        .reset_index()
    )

    pair_spec = [
        ("APP_MML_FREE_Q31", "APP_MML_FREE_Q61"),
        ("MFRMR_MML_FREE_Q31", "MFRMR_MML_FREE_Q61"),
        ("TAM_MML_Q21", "TAM_MML_Q61"),
        ("SIRT_MML_Q30", "SIRT_MML_Q61"),
    ]
    sensitivity_rows: list[dict[str, Any]] = []
    for low_mode, high_mode in pair_spec:
        low = aligned.loc[aligned["Mode"] == low_mode].copy()
        high = aligned.loc[aligned["Mode"] == high_mode].copy()
        for dataset_id in sorted(set(low["DatasetId"]) & set(high["DatasetId"])):
            left = low.loc[low["DatasetId"] == dataset_id, keys + ["EstimateAligned"]]
            right = high.loc[high["DatasetId"] == dataset_id, keys + ["EstimateAligned"]]
            joined = left.merge(right, on=keys, suffixes=("Low", "High"))
            diff = joined["EstimateAlignedHigh"] - joined["EstimateAlignedLow"]
            low_run = enriched.loc[(enriched["DatasetId"] == dataset_id) & (enriched["Mode"] == low_mode)].iloc[0]
            high_run = enriched.loc[(enriched["DatasetId"] == dataset_id) & (enriched["Mode"] == high_mode)].iloc[0]
            sensitivity_rows.append(
                {
                    "DatasetId": dataset_id,
                    "Engine": low_run.Engine,
                    "LowMode": low_mode,
                    "HighMode": high_mode,
                    "SurfaceRMSChange": float(np.sqrt(np.mean(np.asarray(diff, dtype=float) ** 2))),
                    "SurfaceMaxChange": float(np.max(np.abs(diff))),
                    "LogLikChange": as_number(high_run.LogLik) - as_number(low_run.LogLik),
                    "PopulationSDChange": as_number(high_run.PopulationSD) - as_number(low_run.PopulationSD),
                }
            )
    sensitivity = pd.DataFrame(sensitivity_rows)

    negative = enriched.loc[~enriched["ExpectedConnected"].astype(bool)].copy()
    negative["RejectedBeforeOrDuringFit"] = ~negative["FitReturned"].astype(bool)
    negative["FalseReadySignal"] = negative["FitReturned"].astype(bool) & negative["Converged"].astype(bool)

    enriched.to_csv(output_dir / "fit_metrics.csv", index=False)
    aligned.to_csv(output_dir / "surface_estimates_aligned.csv", index=False)
    person_metrics_df.to_csv(output_dir / "person_metrics.csv", index=False)
    summary.to_csv(output_dir / "engine_summary.csv", index=False)
    scenario_summary.to_csv(output_dir / "scenario_summary.csv", index=False)
    connected_summary.to_csv(output_dir / "engine_summary_connected.csv", index=False)
    sensitivity.to_csv(output_dir / "quadrature_sensitivity.csv", index=False)
    negative.to_csv(output_dir / "disconnected_negative_controls.csv", index=False)
    pd.DataFrame(
        [
            {
                "SchemaVersion": SCHEMA_VERSION,
                "Datasets": int(manifest["DatasetId"].nunique()),
                "Attempts": int(len(runs)),
                "FitsReturned": int(runs["FitReturned"].astype(bool).sum()),
                "FailedFits": int((~runs["FitReturned"].astype(bool)).sum()),
                "SurfaceRows": int(len(surfaces)),
                "PersonRows": int(len(persons)),
                "DuplicateRunKeys": 0,
                "ReturnedFitsMissingSurface": 0,
                "ReturnedFitsMissingPersons": 0,
                "IncompleteSurfaceGrids": 0,
            }
        ]
    ).to_csv(output_dir / "run_completeness.csv", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("generate", "fit-python"):
        command = sub.add_parser(name)
        command.add_argument("--input", type=Path)
        command.add_argument("--output", type=Path, required=True)
    command = sub.add_parser("summarize")
    command.add_argument("--input", type=Path, required=True)
    command.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "generate":
        write_generation(args.output)
    elif args.command == "fit-python":
        fit_python(args.input or args.output, args.output)
    elif args.command == "summarize":
        summarize(args.input, args.output)


if __name__ == "__main__":
    main()
