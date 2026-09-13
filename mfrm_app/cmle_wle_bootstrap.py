"""Repository-only bootstrap research core for exact CMLE plus Warm WLE.

Two bootstrap lanes remain deliberately distinct. ``fixed_score_conditional_pattern``
samples response patterns conditional on every observed Person total. It measures
pattern/calibration coupling at fixed totals, not total Person-score uncertainty.
``joint_plugin_parametric`` samples from fitted CMLE category kernels plus fitted
Warm-WLE Person locations. It is a plug-in sensitivity distribution and does not
become a coverage-qualified confidence interval without repeated-truth evidence.
"""

from __future__ import annotations

from itertools import product
import time
from typing import Callable

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from mfrm_app.cmle import fit_cmle, prepare_cmle_design
from mfrm_app.cmle_person_scoring import (
    _ready_cmle_inputs,
    score_cmle_persons_wle,
)
from mfrm_app.cmle_wle_fit import compute_cmle_wle_person_fit


FIXED_SCORE_LANE = "fixed_score_conditional_pattern"
JOINT_PLUGIN_LANE = "joint_plugin_parametric"
SUPPORTED_BOOTSTRAP_LANES = frozenset({FIXED_SCORE_LANE, JOINT_PLUGIN_LANE})
BOOTSTRAP_SCOPE = {
    FIXED_SCORE_LANE: (
        "conditional response-pattern/calibration coupling at fixed observed "
        "Person totals; not total Person-score uncertainty"
    ),
    JOINT_PLUGIN_LANE: (
        "joint plug-in CMLE-WLE sampling sensitivity; not a coverage-qualified "
        "confidence interval"
    ),
}


def _explicit_seed(value: object) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError("seed must be an explicit integer, not boolean.")
    if isinstance(value, (int, np.integer)):
        integer = int(value)
    elif isinstance(value, (float, np.floating)):
        numeric = float(value)
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise ValueError("seed must be an explicit non-negative integer.")
        integer = int(numeric)
    else:
        raise ValueError("seed must be an explicit integer.")
    if integer < 0:
        raise ValueError("seed must be an explicit non-negative integer.")
    return integer


def _positive_integer(value: object, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a positive integer.")
    if isinstance(value, (int, np.integer)):
        integer = int(value)
    elif isinstance(value, (float, np.floating)):
        numeric = float(value)
        if not np.isfinite(numeric) or not numeric.is_integer():
            raise ValueError(f"{name} must be a positive integer.")
        integer = int(numeric)
    else:
        raise ValueError(f"{name} must be a positive integer.")
    if integer < 1:
        raise ValueError(f"{name} must be a positive integer.")
    return integer


def _row_category_kernels(cmle_result: dict[str, object]) -> tuple[object, np.ndarray]:
    design, parameters = _ready_cmle_inputs(cmle_result)
    kernels = design.row_offset + np.einsum(
        "rkp,p->rk", design.row_design, parameters, optimize=True
    )
    if kernels.shape != (len(design.data), design.n_categories):
        raise ValueError("CMLE category-kernel dimensions are inconsistent.")
    if not np.isfinite(kernels).all():
        raise ValueError("CMLE category kernels must be finite.")
    return design, kernels


def _suffix_log_coefficients(kernels: np.ndarray, target: int) -> np.ndarray:
    """Return suffix log normalizers for a fixed integer score total."""
    kernels = np.asarray(kernels, dtype=float)
    if kernels.ndim != 2 or kernels.shape[0] < 1 or kernels.shape[1] < 2:
        raise ValueError("Conditional kernels must be rows by at least two categories.")
    if not np.isfinite(kernels).all():
        raise ValueError("Conditional kernels must be finite.")
    rows, categories = kernels.shape
    maximum_score = rows * (categories - 1)
    if target < 0 or target > maximum_score:
        raise ValueError("Conditional Person total is impossible for the administered rows.")
    suffix = np.full((rows + 1, target + 1), -np.inf, dtype=float)
    suffix[rows, 0] = 0.0
    for row in range(rows - 1, -1, -1):
        for score in range(target + 1):
            candidates = [
                kernels[row, category] + suffix[row + 1, score - category]
                for category in range(min(categories - 1, score) + 1)
                if np.isfinite(suffix[row + 1, score - category])
            ]
            if candidates:
                suffix[row, score] = float(logsumexp(candidates))
    if not np.isfinite(suffix[0, target]):
        raise ValueError("Conditional Person total has zero fitted probability.")
    return suffix


def _sample_fixed_score_pattern(
    kernels: np.ndarray,
    target: int,
    rng: np.random.Generator,
) -> np.ndarray:
    kernels = np.asarray(kernels, dtype=float)
    suffix = _suffix_log_coefficients(kernels, int(target))
    rows, categories = kernels.shape
    sampled = np.zeros(rows, dtype=int)
    remaining = int(target)
    for row in range(rows):
        available: list[int] = []
        log_weights: list[float] = []
        for category in range(min(categories - 1, remaining) + 1):
            tail = suffix[row + 1, remaining - category]
            if np.isfinite(tail):
                available.append(category)
                log_weights.append(float(kernels[row, category] + tail))
        if not available:
            raise RuntimeError("Conditional sampler reached an impossible suffix state.")
        centered = np.asarray(log_weights) - float(logsumexp(log_weights))
        probabilities = np.exp(centered)
        probabilities /= probabilities.sum()
        chosen = int(rng.choice(len(available), p=probabilities))
        sampled[row] = available[chosen]
        remaining -= sampled[row]
    if remaining != 0 or int(sampled.sum()) != int(target):
        raise RuntimeError("Conditional sampler did not preserve the Person total.")
    return sampled


def _enumerated_fixed_score_distribution(
    kernels: np.ndarray,
    target: int,
    *,
    maximum_patterns: int = 100_000,
) -> pd.DataFrame:
    """Enumerate a small conditional distribution for validation only."""
    kernels = np.asarray(kernels, dtype=float)
    if kernels.ndim != 2:
        raise ValueError("kernels must be a two-dimensional array.")
    pattern_count = kernels.shape[1] ** kernels.shape[0]
    if pattern_count > int(maximum_patterns):
        raise ValueError("Conditional enumeration exceeds maximum_patterns.")
    patterns = [
        values
        for values in product(range(kernels.shape[1]), repeat=kernels.shape[0])
        if sum(values) == int(target)
    ]
    if not patterns:
        raise ValueError("No response pattern has the requested conditional total.")
    log_weights = np.asarray(
        [sum(kernels[row, value] for row, value in enumerate(values)) for values in patterns],
        dtype=float,
    )
    probabilities = np.exp(log_weights - float(logsumexp(log_weights)))
    return pd.DataFrame({"Pattern": patterns, "Probability": probabilities})


def generate_cmle_wle_bootstrap_sample(
    cmle_result: dict[str, object],
    *,
    lane: str,
    seed: int,
) -> dict[str, object]:
    """Generate one score-replaced fit-sample data frame for a named lane."""
    lane = str(lane)
    if lane not in SUPPORTED_BOOTSTRAP_LANES:
        raise ValueError(
            "lane must be fixed_score_conditional_pattern or joint_plugin_parametric."
        )
    seed_value = _explicit_seed(seed)
    design, kernels = _row_category_kernels(cmle_result)
    rng = np.random.default_rng(seed_value)
    persons = design.data[design.person_col].astype(str).to_numpy(dtype=object)
    original = (
        pd.to_numeric(design.data[design.score_col], errors="raise").to_numpy(dtype=int)
        - int(design.rating_min)
    )
    sampled = np.full(len(design.data), -1, dtype=int)

    baseline_wle: pd.DataFrame | None = None
    theta_by_person: dict[str, float] = {}
    if lane == JOINT_PLUGIN_LANE:
        baseline_wle = score_cmle_persons_wle(cmle_result)
        if not baseline_wle["WLEPersonScoringReady"].all():
            raise ValueError("All fit-sample WLE scores must be ready for joint sampling.")
        theta_by_person = dict(
            zip(
                baseline_wle["Person"].astype(str),
                pd.to_numeric(baseline_wle["Estimate"], errors="raise").astype(float),
            )
        )
        if set(theta_by_person) != set(persons):
            raise ValueError("WLE Person identity differs from the fitted CMLE design.")

    for person in pd.unique(persons):
        row_indices = np.flatnonzero(persons == person)
        if lane == FIXED_SCORE_LANE:
            target = int(original[row_indices].sum())
            sampled[row_indices] = _sample_fixed_score_pattern(
                kernels[row_indices], target, rng
            )
        else:
            theta = float(theta_by_person[str(person)])
            logits = kernels[row_indices] + theta * np.arange(
                design.n_categories, dtype=float
            )[None, :]
            logits -= logsumexp(logits, axis=1, keepdims=True)
            probabilities = np.exp(logits)
            for local_index, row_index in enumerate(row_indices):
                sampled[row_index] = int(
                    rng.choice(design.n_categories, p=probabilities[local_index])
                )

    if np.any(sampled < 0) or np.any(sampled >= design.n_categories):
        raise RuntimeError("Bootstrap sampler returned an invalid category.")
    original_totals = pd.Series(original, index=persons).groupby(level=0, sort=False).sum()
    sampled_totals = pd.Series(sampled, index=persons).groupby(level=0, sort=False).sum()
    total_mismatches = int((original_totals != sampled_totals).sum())
    if lane == FIXED_SCORE_LANE and total_mismatches:
        raise RuntimeError("Fixed-score bootstrap changed at least one Person total.")

    frame = design.data.copy()
    frame[design.score_col] = sampled + int(design.rating_min)
    frame["__bootstrap_original_internal_score__"] = original
    frame["__bootstrap_lane__"] = lane
    return {
        "data": frame,
        "audit": pd.DataFrame(
            [
                {
                    "Lane": lane,
                    "Seed": seed_value,
                    "Rows": len(frame),
                    "Persons": len(original_totals),
                    "ChangedResponses": int(np.sum(sampled != original)),
                    "ChangedPersonTotals": total_mismatches,
                    "TotalPreservationRequired": lane == FIXED_SCORE_LANE,
                    "TotalPreservationPassed": (
                        total_mismatches == 0 if lane == FIXED_SCORE_LANE else np.nan
                    ),
                    "Scope": BOOTSTRAP_SCOPE[lane],
                }
            ]
        ),
    }


def _refit_arguments(cmle_result: dict[str, object]) -> dict[str, object]:
    design, _ = _ready_cmle_inputs(cmle_result)
    config = cmle_result["config"]
    return {
        "person_col": design.person_col,
        "facet_cols": list(design.facet_cols),
        "score_col": design.score_col,
        "rating_min": design.rating_min,
        "rating_max": design.rating_max,
        "model": design.model,
        "step_facet": design.step_facet,
        "weight_col": config.get("weight_col"),
        "response_unit_col": design.response_unit_col,
        "positive_facets": [
            facet for facet, sign in design.facet_signs.items() if int(sign) > 0
        ],
        "hard_anchors": design.hard_anchors.copy(),
        "allow_duplicate_units": bool(design.audit.get("duplicate_rows", 0)),
        "rank_audit_max_parameters": int(config["rank_audit_max_parameters"]),
        "rank_audit_max_work": int(config["rank_audit_max_work"]),
        "rank_audit_max_bytes": int(config["rank_audit_max_bytes"]),
        "newton_polish_maxiter": int(config["newton_polish_maxiter"]),
    }


def run_cmle_wle_bootstrap(
    cmle_result: dict[str, object],
    *,
    lane: str,
    n_replicates: int,
    seed: int,
    gtol: float = 1e-8,
    maxiter: int = 500,
    maximum_replicates: int = 10_000,
    maximum_row_replicate_work: int = 5_000_000,
    maximum_estimated_output_bytes: int = 512 * 1024 * 1024,
    progress_callback: Callable[[dict[str, object]], None] | None = None,
) -> dict[str, object]:
    """Generate, refit, and score a fully accounted bootstrap lane."""
    lane = str(lane)
    if lane not in SUPPORTED_BOOTSTRAP_LANES:
        raise ValueError("Unsupported CMLE-WLE bootstrap lane.")
    replicates = _positive_integer(n_replicates, "n_replicates")
    maximum_replicates = _positive_integer(maximum_replicates, "maximum_replicates")
    maximum_row_replicate_work = _positive_integer(
        maximum_row_replicate_work, "maximum_row_replicate_work"
    )
    maximum_estimated_output_bytes = _positive_integer(
        maximum_estimated_output_bytes, "maximum_estimated_output_bytes"
    )
    if replicates > maximum_replicates:
        raise ValueError("n_replicates exceeds maximum_replicates.")
    seed_value = _explicit_seed(seed)
    design, _ = _ready_cmle_inputs(cmle_result)
    row_work = int(len(design.data) * replicates)
    if row_work > maximum_row_replicate_work:
        raise ValueError("Requested row-replicate work exceeds the configured maximum.")
    estimated_output_bytes = int(
        replicates * (len(design.data) * 24 + design.audit["persons_total"] * 1536 + 2048)
    )
    if estimated_output_bytes > maximum_estimated_output_bytes:
        raise ValueError("Estimated bootstrap output exceeds the configured memory maximum.")
    if not np.isfinite(float(gtol)) or float(gtol) <= 0:
        raise ValueError("gtol must be finite and positive.")
    maxiter = _positive_integer(maxiter, "maxiter")

    baseline = score_cmle_persons_wle(cmle_result).copy()
    baseline_fit = compute_cmle_wle_person_fit(cmle_result)["persons"].copy()
    baseline_fit_columns = [
        "Person",
        "Infit",
        "Outfit",
        "InfitClass",
        "OutfitClass",
        "WorstFitClass",
        "InfitDisplay",
        "OutfitDisplay",
        "InfitBoundaryStatus",
        "OutfitBoundaryStatus",
        "InfitDecisionStable",
        "OutfitDecisionStable",
        "InfitDisplayDecisionConsistent",
        "OutfitDisplayDecisionConsistent",
        "PersonFitReady",
        "PersonFitStatus",
    ]
    baseline = baseline.merge(
        baseline_fit[baseline_fit_columns],
        on="Person",
        how="left",
        validate="one_to_one",
    )
    if not bool(baseline["PersonFitReady"].all()):
        raise ValueError("Baseline CMLE-WLE Person fit must be ready before bootstrap.")
    baseline_estimate = baseline.set_index("Person")["Estimate"].astype(float)
    baseline_infit = baseline.set_index("Person")["Infit"].astype(float)
    baseline_outfit = baseline.set_index("Person")["Outfit"].astype(float)
    baseline_infit_class = baseline.set_index("Person")["InfitClass"].astype(str)
    baseline_outfit_class = baseline.set_index("Person")["OutfitClass"].astype(str)
    refit_arguments = _refit_arguments(cmle_result)
    children = np.random.SeedSequence(seed_value).spawn(replicates)
    child_seeds = [int(child.generate_state(1, dtype=np.uint64)[0]) for child in children]
    ledger_rows: list[dict[str, object]] = []
    person_parts: list[pd.DataFrame] = []
    generator_parts: list[pd.DataFrame] = []
    started = time.monotonic()
    for replicate, child_seed in enumerate(child_seeds, start=1):
        row: dict[str, object] = {
            "Lane": lane,
            "Model": design.model,
            "Replicate": replicate,
            "SeedIdentity": child_seed,
            "CMLEConverged": False,
            "CMLEInferenceReady": False,
            "Rank": np.nan,
            "Nullity": np.nan,
            "WLEAvailable": False,
            "PersonFitAvailable": False,
            "PersonsFitReady": 0,
            "PersonsFitTotal": int(len(baseline)),
            "FailureStage": "",
            "FailureReason": "",
        }
        try:
            generated = generate_cmle_wle_bootstrap_sample(
                cmle_result, lane=lane, seed=child_seed
            )
            generator_audit = generated["audit"].copy()
            generator_audit.insert(2, "Replicate", replicate)
            generator_parts.append(generator_audit)
            preflight_arguments = dict(refit_arguments)
            preflight_arguments.pop("newton_polish_maxiter")
            prepared = prepare_cmle_design(
                generated["data"],
                **preflight_arguments,
            )
            row["Rank"] = int(prepared.audit["conditional_rank"])
            row["Nullity"] = int(prepared.audit["conditional_nullity"])
            if not bool(prepared.audit["eligible"]):
                blocking = [
                    issue["Message"]
                    for issue in prepared.audit["issues"]
                    if issue["Severity"] == "Block"
                ]
                row["FailureStage"] = "cmle_eligibility"
                row["FailureReason"] = " ".join(blocking)
                row["ElapsedSeconds"] = float(time.monotonic() - started)
                ledger_rows.append(row)
                if progress_callback is not None:
                    progress_callback(dict(row))
                continue
            refit = fit_cmle(
                generated["data"],
                **refit_arguments,
                maxiter=maxiter,
                gtol=float(gtol),
            )
            fit_row = refit["summary"].iloc[0]
            row["CMLEConverged"] = bool(fit_row["Converged"])
            row["CMLEInferenceReady"] = bool(fit_row["InferenceReady"])
            row["Rank"] = int(fit_row["InformationRank"])
            row["Nullity"] = int(fit_row["InformationNullity"])
            if not row["CMLEInferenceReady"]:
                row["FailureStage"] = "cmle_inference_readiness"
                row["FailureReason"] = str(fit_row["ReadinessReasons"])
            else:
                scored = score_cmle_persons_wle(refit).copy()
                row["WLEAvailable"] = bool(scored["WLEPersonScoringReady"].all())
                fit_columns = baseline_fit_columns
                if not row["WLEAvailable"]:
                    row["FailureStage"] = "wle_scoring"
                    row["FailureReason"] = "At least one Person WLE score was unavailable."
                    for column in fit_columns:
                        if column != "Person":
                            scored[column] = np.nan
                else:
                    try:
                        fit_result = compute_cmle_wle_person_fit(refit)
                        fit_persons = fit_result["persons"].copy()
                        scored = scored.merge(
                            fit_persons[fit_columns],
                            on="Person",
                            how="left",
                            validate="one_to_one",
                        )
                        row["PersonsFitReady"] = int(
                            scored["PersonFitReady"].fillna(False).sum()
                        )
                        row["PersonsFitTotal"] = int(len(scored))
                        row["PersonFitAvailable"] = bool(
                            row["PersonsFitReady"] == row["PersonsFitTotal"]
                            and np.isfinite(scored[["Infit", "Outfit"]]).all().all()
                        )
                        if not row["PersonFitAvailable"]:
                            row["FailureStage"] = "person_fit"
                            row["FailureReason"] = (
                                "At least one Person Infit/Outfit was unavailable."
                            )
                    except Exception as fit_exc:
                        row["FailureStage"] = "person_fit"
                        row["FailureReason"] = f"{type(fit_exc).__name__}: {fit_exc}"
                        for column in fit_columns:
                            if column != "Person":
                                scored[column] = np.nan
                scored.insert(0, "Replicate", replicate)
                scored.insert(0, "Model", design.model)
                scored.insert(0, "Lane", lane)
                scored["BaselineEstimate"] = (
                    scored["Person"].astype(str).map(baseline_estimate)
                )
                scored["BootstrapMinusBaseline"] = (
                    scored["Estimate"] - scored["BaselineEstimate"]
                )
                scored["BaselineInfit"] = scored["Person"].astype(str).map(
                    baseline_infit
                )
                scored["BaselineOutfit"] = scored["Person"].astype(str).map(
                    baseline_outfit
                )
                scored["BaselineInfitClass"] = scored["Person"].astype(str).map(
                    baseline_infit_class
                )
                scored["BaselineOutfitClass"] = scored["Person"].astype(str).map(
                    baseline_outfit_class
                )
                fit_ready = scored.get(
                    "PersonFitReady", pd.Series(False, index=scored.index)
                ).fillna(False).astype(bool)
                scored["InfitZoneTransition"] = np.where(
                    fit_ready,
                    scored["BaselineInfitClass"].astype(str)
                    + "->"
                    + scored["InfitClass"].astype(str),
                    "unavailable",
                )
                scored["OutfitZoneTransition"] = np.where(
                    fit_ready,
                    scored["BaselineOutfitClass"].astype(str)
                    + "->"
                    + scored["OutfitClass"].astype(str),
                    "unavailable",
                )
                scored["InfitZoneChanged"] = fit_ready & scored[
                    "InfitClass"
                ].astype(str).ne(scored["BaselineInfitClass"].astype(str))
                scored["OutfitZoneChanged"] = fit_ready & scored[
                    "OutfitClass"
                ].astype(str).ne(scored["BaselineOutfitClass"].astype(str))
                scored["AnyFitZoneTransition"] = (
                    scored["InfitZoneChanged"] | scored["OutfitZoneChanged"]
                )
                display_consistent = (
                    scored.get(
                        "InfitDisplayDecisionConsistent",
                        pd.Series(False, index=scored.index),
                    )
                    .fillna(False)
                    .astype(bool)
                    & scored.get(
                        "OutfitDisplayDecisionConsistent",
                        pd.Series(False, index=scored.index),
                    )
                    .fillna(False)
                    .astype(bool)
                )
                scored["RawDisplayMismatch"] = fit_ready & ~display_consistent
                person_parts.append(scored)
        except Exception as exc:  # fail closed while preserving the replicate
            row["FailureStage"] = "bootstrap_or_refit"
            row["FailureReason"] = f"{type(exc).__name__}: {exc}"
        row["ElapsedSeconds"] = float(time.monotonic() - started)
        ledger_rows.append(row)
        if progress_callback is not None:
            progress_callback(dict(row))

    ledger = pd.DataFrame(ledger_rows)
    person_draws = (
        pd.concat(person_parts, ignore_index=True) if person_parts else pd.DataFrame()
    )
    generator_audit = (
        pd.concat(generator_parts, ignore_index=True)
        if generator_parts else pd.DataFrame()
    )
    summary = pd.DataFrame(
        [
            {
                "Lane": lane,
                "Model": design.model,
                "AttemptedReplicates": replicates,
                "GeneratedReplicates": len(generator_audit),
                "CMLEConvergedReplicates": int(ledger["CMLEConverged"].sum()),
                "CMLEInferenceReadyReplicates": int(
                    ledger["CMLEInferenceReady"].sum()
                ),
                "WLEAvailableReplicates": int(ledger["WLEAvailable"].sum()),
                "PersonFitAvailableReplicates": int(
                    ledger["PersonFitAvailable"].sum()
                ),
                "FailedReplicates": int((~ledger["WLEAvailable"]).sum()),
                "InferenceReadyShare": float(ledger["CMLEInferenceReady"].mean()),
                "WLEAvailableShare": float(ledger["WLEAvailable"].mean()),
                "PersonFitAvailableShare": float(
                    ledger["PersonFitAvailable"].mean()
                ),
                "RowReplicateWork": row_work,
                "EstimatedOutputBytes": estimated_output_bytes,
                "ElapsedSeconds": float(time.monotonic() - started),
                "Scope": BOOTSTRAP_SCOPE[lane],
                "IntervalCoverageQualified": False,
                "TotalInferentialSEQualified": False,
                "FitZSTDQualified": False,
                "FitPValueQualified": False,
                "StreamlitIntegrationAuthorized": False,
            }
        ]
    )
    return {
        "summary": summary,
        "ledger": ledger,
        "generator_audit": generator_audit,
        "person_draws": person_draws,
        "baseline_persons": baseline,
    }


__all__ = [
    "FIXED_SCORE_LANE",
    "JOINT_PLUGIN_LANE",
    "SUPPORTED_BOOTSTRAP_LANES",
    "generate_cmle_wle_bootstrap_sample",
    "run_cmle_wle_bootstrap",
]
