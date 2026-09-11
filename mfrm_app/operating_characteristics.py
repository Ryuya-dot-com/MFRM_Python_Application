"""Streamlit-free contracts for Monte Carlo operating-characteristic studies.

The functions in this module do not fit an MFRM.  They make the orchestration
and aggregation around a fitting engine explicit and reproducible:

* condition and replicate identities are deterministic;
* common-random-number pairing is recorded instead of inferred;
* attempted, failed, unavailable, and eligible decisions remain distinct;
* false-positive rate and power are never pooled across different truths;
* Monte Carlo uncertainty accompanies every rate; and
* small pilot runs are labelled as such rather than promoted to validation.

Engine-specific simulation and fitting live in validation runners.  Keeping
this layer independent of Streamlit lets Python, mfrmr, TAM, immer, and sirt
emit the same input schema in later cross-engine studies.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


SCHEMA_VERSION = "mfrm-operating-characteristics-v1"
PRECISION_PLAN_SCHEMA_VERSION = "mfrm-operating-characteristics-precision-plan-v1"
DEFAULT_NOMINAL_CI = 0.95
DEFAULT_MIN_STUDY_REPLICATES = 100


def _json_scalar(value: Any) -> Any:
    """Return a stable JSON-compatible representation of a scalar."""

    if value is None or value is pd.NA:
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        if not np.isfinite(numeric):
            return str(numeric)
        return numeric
    if isinstance(value, (list, tuple)):
        return [_json_scalar(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _json_scalar(value[key]) for key in sorted(value, key=str)}
    return str(value)


def canonical_condition_id(settings: Mapping[str, Any], *, prefix: str = "condition") -> str:
    """Build a compact content-derived condition identifier."""

    payload = json.dumps(
        {str(key): _json_scalar(settings[key]) for key in sorted(settings, key=str)},
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    )
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
    safe_prefix = "".join(char if char.isalnum() or char in "-_" else "-" for char in str(prefix))
    return f"{safe_prefix or 'condition'}-{digest}"


def deterministic_replicate_seed(
    base_seed: int,
    seed_group: str,
    replicate: int,
) -> int:
    """Return a stable uint32 seed independent of Python's hash randomization."""

    if int(replicate) < 1:
        raise ValueError("replicate must be >= 1")
    payload = f"{int(base_seed)}|{str(seed_group)}|{int(replicate)}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big", signed=False)


def build_replicate_manifest(
    conditions: Iterable[Mapping[str, Any]],
    *,
    replicates: int,
    base_seed: int,
    max_runs: int = 10_000,
) -> pd.DataFrame:
    """Expand condition dictionaries to one auditable row per attempted run.

    Each condition must have a unique ``ConditionId`` or one is derived from
    its content.  ``SeedGroup`` may intentionally be shared by conditions to
    implement common random numbers (for example, null and alternative bias
    conditions with the same latent draws).  Shared seeds are therefore valid,
    but duplicate ``(ConditionId, Replicate)`` identities are not.
    """

    reps = int(replicates)
    if reps < 1:
        raise ValueError("replicates must be >= 1")
    condition_rows: list[dict[str, Any]] = []
    for raw in conditions:
        row = {str(key): value for key, value in dict(raw).items()}
        condition_id = str(row.get("ConditionId", "")).strip()
        if not condition_id:
            condition_id = canonical_condition_id(row)
        row["ConditionId"] = condition_id
        row["SeedGroup"] = str(row.get("SeedGroup", condition_id))
        condition_rows.append(row)
    if not condition_rows:
        return pd.DataFrame()
    ids = [row["ConditionId"] for row in condition_rows]
    if len(ids) != len(set(ids)):
        raise ValueError("ConditionId values must be unique")
    attempted = len(condition_rows) * reps
    if attempted > int(max_runs):
        raise ValueError(
            f"requested {attempted} runs exceeds max_runs={int(max_runs)}; "
            "raise the explicit budget only after reviewing runtime and storage"
        )

    rows: list[dict[str, Any]] = []
    for condition_order, condition in enumerate(condition_rows, start=1):
        for replicate in range(1, reps + 1):
            seed = deterministic_replicate_seed(
                int(base_seed),
                str(condition["SeedGroup"]),
                replicate,
            )
            rows.append({
                "SchemaVersion": SCHEMA_VERSION,
                "ConditionOrder": condition_order,
                **condition,
                "Replicate": replicate,
                "Seed": seed,
                "BaseSeed": int(base_seed),
                "SeedCoupling": (
                    "common_random_numbers"
                    if sum(row["SeedGroup"] == condition["SeedGroup"] for row in condition_rows) > 1
                    else "independent_condition_stream"
                ),
                "RunId": f"{condition['ConditionId']}::rep-{replicate:05d}",
            })
    manifest = pd.DataFrame(rows)
    if manifest.duplicated(["ConditionId", "Replicate"]).any():
        raise RuntimeError("replicate manifest contains duplicate run identities")
    return manifest


def audit_manifest_extension(
    prior: pd.DataFrame,
    expanded: pd.DataFrame,
) -> dict[str, Any]:
    """Verify that a larger manifest is a seed-preserving replicate extension.

    The manifest is condition-major, so an expanded manifest is not a literal
    row prefix of the smaller one.  Extension instead means that every prior
    ``RunId`` and every value on that row are preserved, condition definitions
    are unchanged, and newly added replicate indices continue after the
    previous per-condition maximum without gaps.
    """

    if not isinstance(prior, pd.DataFrame) or prior.empty:
        raise ValueError("prior manifest must be a non-empty DataFrame")
    if not isinstance(expanded, pd.DataFrame) or expanded.empty:
        raise ValueError("expanded manifest must be a non-empty DataFrame")
    required = {"RunId", "ConditionId", "Replicate", "Seed", "BaseSeed"}
    for label, frame in (("prior", prior), ("expanded", expanded)):
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{label} manifest missing required columns: {sorted(missing)}")
        if frame["RunId"].astype(str).duplicated().any():
            raise ValueError(f"{label} manifest contains duplicate RunId values")
        if frame.duplicated(["ConditionId", "Replicate"]).any():
            raise ValueError(
                f"{label} manifest contains duplicate ConditionId/Replicate values"
            )
    if list(prior.columns) != list(expanded.columns):
        raise ValueError("manifest columns or column order changed during extension")

    prior_ids = prior["RunId"].astype(str)
    expanded_ids = set(expanded["RunId"].astype(str))
    missing_ids = sorted(set(prior_ids) - expanded_ids)
    if missing_ids:
        raise ValueError(f"expanded manifest dropped prior RunIds: {missing_ids[:3]}")

    expanded_by_id = expanded.assign(
        _RunIdKey=expanded["RunId"].astype(str)
    ).set_index("_RunIdKey", drop=True)
    retained = expanded_by_id.loc[prior_ids.tolist(), prior.columns].reset_index(drop=True)
    prior_fingerprint = frame_fingerprint(prior.reset_index(drop=True))
    retained_fingerprint = frame_fingerprint(retained)
    if prior_fingerprint != retained_fingerprint:
        raise ValueError("expanded manifest changed one or more prior row values")

    prior_conditions = set(prior["ConditionId"].astype(str))
    expanded_conditions = set(expanded["ConditionId"].astype(str))
    if prior_conditions != expanded_conditions:
        raise ValueError("condition set changed during manifest extension")

    new_rows = expanded.loc[~expanded["RunId"].astype(str).isin(set(prior_ids))].copy()
    if new_rows.empty:
        raise ValueError("expanded manifest adds no new replicate rows")
    run_specific = {"Replicate", "Seed", "RunId"}
    condition_columns = [column for column in prior.columns if column not in run_specific]
    for condition_id in sorted(prior_conditions):
        before = prior.loc[prior["ConditionId"].astype(str).eq(condition_id)].copy()
        after = expanded.loc[expanded["ConditionId"].astype(str).eq(condition_id)].copy()
        before_reps = sorted(pd.to_numeric(before["Replicate"], errors="raise").astype(int))
        after_reps = sorted(pd.to_numeric(after["Replicate"], errors="raise").astype(int))
        if before_reps != list(range(1, max(before_reps) + 1)):
            raise ValueError(f"prior replicate sequence has gaps for {condition_id}")
        if after_reps != list(range(1, max(after_reps) + 1)):
            raise ValueError(f"expanded replicate sequence has gaps for {condition_id}")
        if max(after_reps) <= max(before_reps):
            raise ValueError(f"expanded manifest did not extend {condition_id}")
        added_reps = [replicate for replicate in after_reps if replicate not in before_reps]
        if added_reps != list(range(max(before_reps) + 1, max(after_reps) + 1)):
            raise ValueError(f"new replicate indices do not append after prior rows for {condition_id}")
        before_definition = frame_fingerprint(before[condition_columns].head(1).reset_index(drop=True))
        after_definitions = after[condition_columns].drop_duplicates().reset_index(drop=True)
        if len(after_definitions) != 1 or frame_fingerprint(after_definitions) != before_definition:
            raise ValueError(f"condition definition changed during extension for {condition_id}")

    return {
        "Passed": True,
        "PriorRuns": int(len(prior)),
        "ExpandedRuns": int(len(expanded)),
        "AddedRuns": int(len(new_rows)),
        "Conditions": int(len(prior_conditions)),
        "PriorManifestSHA256": prior_fingerprint,
        "RetainedRowsSHA256": retained_fingerprint,
        "ExpandedManifestSHA256": frame_fingerprint(expanded.reset_index(drop=True)),
        "Boundary": (
            "This proves deterministic seed/condition extension only; it does not "
            "authorize pooling pilot outcomes into confirmatory performance estimates."
        ),
    }


def load_precision_plan(path: str | Path) -> dict[str, Any]:
    """Load and validate the registered operating-characteristics precision plan."""

    plan_path = Path(path)
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if not isinstance(plan, dict):
        raise ValueError("precision plan must be a JSON object")
    if plan.get("schema_version") != PRECISION_PLAN_SCHEMA_VERSION:
        raise ValueError(
            f"precision plan schema must be {PRECISION_PLAN_SCHEMA_VERSION!r}"
        )
    profiles = plan.get("profiles")
    if not isinstance(profiles, dict):
        raise ValueError("precision plan profiles must be an object")
    expected_profiles = {"smoke": 2, "pilot": 20, "study": 100}
    for profile, expected in expected_profiles.items():
        entry = profiles.get(profile)
        if not isinstance(entry, dict):
            raise ValueError(f"precision plan missing {profile} profile")
        if int(entry.get("replicates_per_condition", -1)) != expected:
            raise ValueError(
                f"precision plan {profile} profile must retain {expected} replicates per condition"
            )

    decisions = plan.get("authoritative_decisions")
    if not isinstance(decisions, dict):
        raise ValueError("precision plan authoritative_decisions must be an object")
    if decisions.get("classification_values") != "finite_unrounded":
        raise ValueError("primary decisions must use finite unrounded values")
    if decisions.get("primary_bias_rule") != "DecisionStrongRaw":
        raise ValueError("primary bias rule must remain DecisionStrongRaw")
    if not bool(decisions.get("boundary_audit_required", False)):
        raise ValueError("floating-point boundary auditing must remain required")

    precision = plan.get("confirmatory_precision")
    if not isinstance(precision, dict) or not precision:
        raise ValueError("precision plan confirmatory_precision must be a non-empty object")
    for estimand, entry in precision.items():
        if not isinstance(entry, dict):
            raise ValueError(f"precision entry {estimand} must be an object")
        rate = float(entry.get("reference_rate", np.nan))
        minimum = int(entry.get("minimum_eligible_per_condition", 0))
        target_mcse = float(entry.get("target_mcse", np.nan))
        target_half = float(entry.get("target_wilson_half_width_95", np.nan))
        if not (0.0 < rate < 1.0 and minimum > 0):
            raise ValueError(f"invalid rate or eligible minimum for {estimand}")
        if not (np.isfinite(target_mcse) and target_mcse > 0):
            raise ValueError(f"invalid MCSE target for {estimand}")
        if not (np.isfinite(target_half) and target_half > 0):
            raise ValueError(f"invalid Wilson half-width target for {estimand}")
        achieved_mcse = math.sqrt(rate * (1.0 - rate) / minimum)
        expected_events = int(round(rate * minimum))
        lower, upper = wilson_interval(expected_events, minimum)
        achieved_half = (upper - lower) / 2.0
        if achieved_mcse > target_mcse + 1e-15:
            raise ValueError(
                f"{estimand} minimum does not attain its registered MCSE target"
            )
        if achieved_half > target_half + 1e-15:
            raise ValueError(
                f"{estimand} minimum does not attain its registered Wilson target"
            )

    budget = plan.get("budget_rule")
    if not isinstance(budget, dict):
        raise ValueError("precision plan budget_rule must be an object")
    if budget.get("eligibility_rate_bound") != "two_sided_wilson_lower_95":
        raise ValueError("attempt inflation must use the registered Wilson lower bound")
    if not bool(budget.get("pilot_excluded_from_confirmatory_performance", False)):
        raise ValueError("pilot outcomes must remain excluded from confirmatory performance")
    if int(budget.get("maximum_attempts_per_condition", 0)) < max(
        int(entry["minimum_eligible_per_condition"]) for entry in precision.values()
    ):
        raise ValueError("maximum attempt cap is below a registered eligible target")
    return plan


def plan_confirmatory_attempts(
    eligible: int,
    attempts: int,
    *,
    target_eligible: int,
    maximum_attempts: int,
) -> dict[str, Any]:
    """Inflate a fixed confirmatory budget from the pilot eligibility lower bound.

    Only eligibility/failure accounting is used here.  Decision rates, power,
    recovery, and coverage outcomes must not be consulted when fixing the
    confirmatory manifest.
    """

    n = int(attempts)
    x = int(eligible)
    target = int(target_eligible)
    cap = int(maximum_attempts)
    if n <= 0 or x < 0 or x > n:
        raise ValueError("eligible and attempts must describe a valid binomial count")
    if target <= 0 or cap < target:
        raise ValueError("target and maximum attempts must be positive with maximum >= target")
    lower, upper = wilson_interval(x, n)
    observed = x / n
    required = math.ceil(target / lower) if np.isfinite(lower) and lower > 0 else None
    if required is None:
        status = "blocked_for_redesign_zero_lower_bound"
    elif required > cap:
        status = "blocked_for_redesign_attempt_cap"
    else:
        status = "fixed_attempt_budget_available"
    return {
        "PilotAttempts": n,
        "PilotEligible": x,
        "ObservedEligibilityRate": observed,
        "EligibilityWilsonLower95": lower,
        "EligibilityWilsonUpper95": upper,
        "TargetEligible": target,
        "RequiredFixedAttempts": required,
        "MaximumAttempts": cap,
        "Status": status,
        "Boundary": (
            "Use pilot eligibility and failure accounting only. Freeze a fixed "
            "confirmatory manifest before inspecting confirmatory decisions."
        ),
    }


def frame_fingerprint(frame: pd.DataFrame) -> str:
    """Return a stable SHA-256 fingerprint for a tabular study manifest."""

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    records = [
        {str(column): _json_scalar(row[column]) for column in frame.columns}
        for _, row in frame.reset_index(drop=True).iterrows()
    ]
    payload = json.dumps(records, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _coerce_bool(series: pd.Series) -> pd.Series:
    """Coerce common bool encodings while retaining unavailable values."""

    if pd.api.types.is_bool_dtype(series.dtype):
        return series.astype("boolean")
    lowered = series.astype("string").str.strip().str.lower()
    mapped = lowered.map({
        "true": True,
        "1": True,
        "yes": True,
        "y": True,
        "false": False,
        "0": False,
        "no": False,
        "n": False,
    })
    return mapped.astype("boolean")


def wilson_interval(successes: int, trials: int, *, z: float = 1.959963984540054) -> tuple[float, float]:
    """Wilson score interval for a binomial rate."""

    n = int(trials)
    x = int(successes)
    if n <= 0 or x < 0 or x > n:
        return float("nan"), float("nan")
    p = x / n
    z2 = float(z) ** 2
    denom = 1.0 + z2 / n
    centre = (p + z2 / (2.0 * n)) / denom
    half = float(z) * math.sqrt((p * (1.0 - p) + z2 / (4.0 * n)) / n) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def monte_carlo_rate_summary(successes: int, trials: int) -> dict[str, float | int]:
    """Summarize one Bernoulli operating characteristic with Monte Carlo error."""

    n = int(trials)
    x = int(successes)
    if n <= 0:
        return {
            "Events": 0,
            "Trials": 0,
            "Rate": np.nan,
            "MonteCarloSE": np.nan,
            "WilsonLower95": np.nan,
            "WilsonUpper95": np.nan,
        }
    rate = x / n
    lower, upper = wilson_interval(x, n)
    return {
        "Events": x,
        "Trials": n,
        "Rate": rate,
        "MonteCarloSE": math.sqrt(rate * (1.0 - rate) / n),
        "WilsonLower95": lower,
        "WilsonUpper95": upper,
    }


def evidence_tier(eligible_replicates: int, *, min_study_replicates: int = DEFAULT_MIN_STUDY_REPLICATES) -> str:
    """Label Monte Carlo depth without making a model-validity claim."""

    n = int(eligible_replicates)
    minimum = max(1, int(min_study_replicates))
    if n <= 0:
        return "no eligible decisions"
    if n < minimum:
        return f"pilot only (<{minimum} eligible replicates)"
    if n < 5 * minimum:
        return f"screening depth ({minimum}-{5 * minimum - 1} eligible replicates)"
    return f"higher-precision Monte Carlo (>= {5 * minimum} eligible replicates)"


def summarize_binary_operating_characteristics(
    replicate_results: pd.DataFrame,
    *,
    decision_columns: Sequence[str],
    group_columns: Sequence[str] = ("ConditionId", "Engine", "Estimator"),
    truth_column: str = "TruthPositive",
    eligible_column: str = "AnalysisEligible",
    min_study_replicates: int = DEFAULT_MIN_STUDY_REPLICATES,
) -> pd.DataFrame:
    """Summarize decision rate, false-positive rate, and power by condition.

    A decision is eligible only when ``eligible_column`` is true and the
    decision itself is available.  Ineligible and unavailable rows stay in
    the denominator accounting but never become implicit negative decisions.
    """

    if not isinstance(replicate_results, pd.DataFrame) or replicate_results.empty:
        return pd.DataFrame()
    required = {truth_column, eligible_column, *decision_columns}
    missing = required - set(replicate_results.columns)
    if missing:
        raise ValueError(f"replicate results missing required columns: {sorted(missing)}")
    groups = [column for column in group_columns if column in replicate_results.columns]
    if not groups:
        raise ValueError("at least one group column must be present")

    rows: list[dict[str, Any]] = []
    grouped = replicate_results.groupby(groups, dropna=False, sort=False)
    for group_key, part in grouped:
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(groups, keys))
        truth = _coerce_bool(part[truth_column])
        base_eligible = _coerce_bool(part[eligible_column]).fillna(False)
        for decision_column in decision_columns:
            decision = _coerce_bool(part[decision_column])
            eligible = base_eligible & decision.notna() & truth.notna()
            used = part.loc[eligible].copy()
            used_truth = truth.loc[eligible].astype(bool)
            used_decision = decision.loc[eligible].astype(bool)
            tp = int((used_truth & used_decision).sum())
            fn = int((used_truth & ~used_decision).sum())
            fp = int((~used_truth & used_decision).sum())
            tn = int((~used_truth & ~used_decision).sum())
            positive = monte_carlo_rate_summary(tp, tp + fn)
            null = monte_carlo_rate_summary(fp, fp + tn)
            overall = monte_carlo_rate_summary(int(used_decision.sum()), len(used_decision))
            rows.append({
                "SchemaVersion": SCHEMA_VERSION,
                **identity,
                "DecisionRule": decision_column,
                "Attempts": int(len(part)),
                "EligibleDecisions": int(eligible.sum()),
                "UnavailableOrIneligible": int(len(part) - eligible.sum()),
                "TruthPositiveReplicates": int(used_truth.sum()),
                "TruthNullReplicates": int((~used_truth).sum()),
                "TruePositive": tp,
                "FalseNegative": fn,
                "FalsePositive": fp,
                "TrueNegative": tn,
                "DecisionRate": overall["Rate"],
                "DecisionRateMCSE": overall["MonteCarloSE"],
                "DecisionRateWilsonLower95": overall["WilsonLower95"],
                "DecisionRateWilsonUpper95": overall["WilsonUpper95"],
                "Power": positive["Rate"],
                "PowerMCSE": positive["MonteCarloSE"],
                "PowerWilsonLower95": positive["WilsonLower95"],
                "PowerWilsonUpper95": positive["WilsonUpper95"],
                "FalsePositiveRate": null["Rate"],
                "FalsePositiveRateMCSE": null["MonteCarloSE"],
                "FalsePositiveRateWilsonLower95": null["WilsonLower95"],
                "FalsePositiveRateWilsonUpper95": null["WilsonUpper95"],
                "EvidenceTier": evidence_tier(
                    int(eligible.sum()),
                    min_study_replicates=min_study_replicates,
                ),
                "InterpretationBoundary": (
                    "Rates describe only eligible returned decisions in this condition; "
                    "failures and unavailable decisions are reported separately and must not be counted as negatives."
                ),
            })
    return pd.DataFrame(rows)


def summarize_estimation_operating_characteristics(
    estimates: pd.DataFrame,
    *,
    group_columns: Sequence[str] = ("ConditionId", "Engine", "Estimator", "ParameterType"),
    error_column: str = "ErrorAligned",
    se_column: str = "SE",
    included_column: str = "IncludedInSummary",
    nominal: float = DEFAULT_NOMINAL_CI,
) -> pd.DataFrame:
    """Aggregate bias, RMSE, MAE, SE availability, and Wald coverage."""

    if not isinstance(estimates, pd.DataFrame) or estimates.empty:
        return pd.DataFrame()
    missing = {error_column, included_column} - set(estimates.columns)
    if missing:
        raise ValueError(f"estimate rows missing required columns: {sorted(missing)}")
    groups = [column for column in group_columns if column in estimates.columns]
    if not groups:
        raise ValueError("at least one estimation group column must be present")
    if not (0.0 < float(nominal) < 1.0):
        raise ValueError("nominal must be between 0 and 1")

    z_value = 1.959963984540054 if abs(float(nominal) - 0.95) < 1e-12 else float("nan")
    rows: list[dict[str, Any]] = []
    for group_key, part in estimates.groupby(groups, dropna=False, sort=False):
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(groups, keys))
        included = _coerce_bool(part[included_column]).fillna(False)
        error = pd.to_numeric(part[error_column], errors="coerce")
        usable = included & error.notna() & np.isfinite(error)
        error_values = error.loc[usable].to_numpy(dtype=float)
        n = int(error_values.size)
        bias = float(np.mean(error_values)) if n else np.nan
        rmse = float(np.sqrt(np.mean(np.square(error_values)))) if n else np.nan
        mae = float(np.mean(np.abs(error_values))) if n else np.nan
        mcse_bias = float(np.std(error_values, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
        squared = np.square(error_values)
        mcse_rmse = (
            float(np.std(squared, ddof=1) / (2.0 * max(rmse, np.finfo(float).tiny) * np.sqrt(n)))
            if n > 1 and np.isfinite(rmse) and rmse > 0 else np.nan
        )
        se = pd.to_numeric(
            part[se_column] if se_column in part.columns else pd.Series(np.nan, index=part.index),
            errors="coerce",
        )
        se_ok = usable & se.notna() & np.isfinite(se) & (se > 0)
        if np.isfinite(z_value):
            covered = error.loc[se_ok].abs() <= z_value * se.loc[se_ok]
            coverage = float(covered.mean()) if len(covered) else np.nan
            coverage_mcse = (
                math.sqrt(coverage * (1.0 - coverage) / len(covered))
                if len(covered) and np.isfinite(coverage) else np.nan
            )
        else:
            coverage = coverage_mcse = np.nan
        rows.append({
            "SchemaVersion": SCHEMA_VERSION,
            **identity,
            "RowsAttempted": int(len(part)),
            "RowsIncluded": n,
            "RowsUnavailable": int(len(part) - n),
            "Bias": bias,
            "MonteCarloSEBias": mcse_bias,
            "RMSE": rmse,
            "MonteCarloSERMSE": mcse_rmse,
            "MAE": mae,
            "MeanSE": float(se.loc[se_ok].mean()) if int(se_ok.sum()) else np.nan,
            "SEAvailableRate": float(se_ok.sum() / max(int(usable.sum()), 1)),
            "CoverageN": int(se_ok.sum()),
            "NominalCoverage": float(nominal),
            "Coverage": coverage,
            "CoverageMonteCarloSE": coverage_mcse,
            "CoverageError": coverage - float(nominal) if np.isfinite(coverage) else np.nan,
            "InterpretationBoundary": (
                "Coverage is a conditional-Wald diagnostic on the supplied comparison scale "
                "(for example, mean-aligned for an unidentified location or absolute for an "
                "anchor-identified facet). It includes only rows with finite positive SE and "
                "does not include failed fits or unavailable SE rows."
            ),
        })
    return pd.DataFrame(rows)


def summarize_run_accounting(
    runs: pd.DataFrame,
    *,
    group_columns: Sequence[str] = ("ConditionId", "Engine", "Estimator"),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return condition-level run accounting and normalized failure reasons."""

    if not isinstance(runs, pd.DataFrame) or runs.empty:
        return pd.DataFrame(), pd.DataFrame()
    groups = [column for column in group_columns if column in runs.columns]
    if not groups:
        raise ValueError("at least one run-accounting group column must be present")
    work = runs.copy()
    for column in ("FitReturned", "Converged", "InferenceReady", "AnalysisEligible"):
        work[column] = _coerce_bool(
            work[column] if column in work.columns else pd.Series(False, index=work.index)
        ).fillna(False)
    work["FailureStage"] = work.get("FailureStage", pd.Series("", index=work.index)).fillna("").astype(str)
    work["FailureReason"] = work.get("FailureReason", pd.Series("", index=work.index)).fillna("").astype(str)

    summary_rows: list[dict[str, Any]] = []
    reason_rows: list[dict[str, Any]] = []
    for group_key, part in work.groupby(groups, dropna=False, sort=False):
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(groups, keys))
        attempts = int(len(part))
        returned = int(part["FitReturned"].sum())
        converged = int(part["Converged"].sum())
        ready = int(part["InferenceReady"].sum())
        eligible = int(part["AnalysisEligible"].sum())
        summary_rows.append({
            "SchemaVersion": SCHEMA_VERSION,
            **identity,
            "Attempts": attempts,
            "FitsReturned": returned,
            "FitFailures": attempts - returned,
            "Converged": converged,
            "NonConverged": attempts - converged,
            "InferenceReady": ready,
            "AnalysisEligible": eligible,
            "FitReturnRate": returned / attempts if attempts else np.nan,
            "ConvergenceRate": converged / attempts if attempts else np.nan,
            "InferenceReadyRate": ready / attempts if attempts else np.nan,
            "AnalysisEligibleRate": eligible / attempts if attempts else np.nan,
            "EvidenceTier": evidence_tier(eligible),
        })
        failed = part.loc[~part["AnalysisEligible"]].copy()
        if not failed.empty:
            failed["FailureStage"] = failed["FailureStage"].replace("", "unspecified")
            failed["FailureReason"] = failed["FailureReason"].replace("", "unspecified")
            counts = (
                failed.groupby(["FailureStage", "FailureReason"], dropna=False)
                .size()
                .reset_index(name="Count")
            )
            for row in counts.itertuples(index=False):
                reason_rows.append({
                    "SchemaVersion": SCHEMA_VERSION,
                    **identity,
                    "FailureStage": row.FailureStage,
                    "FailureReason": row.FailureReason,
                    "Count": int(row.Count),
                    "ShareOfAttempts": int(row.Count) / attempts if attempts else np.nan,
                })
    reason_columns = [
        "SchemaVersion",
        *groups,
        "FailureStage",
        "FailureReason",
        "Count",
        "ShareOfAttempts",
    ]
    return pd.DataFrame(summary_rows), pd.DataFrame(reason_rows, columns=reason_columns)


def audit_conclusion_sensitivity(
    decisions: pd.DataFrame,
    *,
    raw_column: str,
    comparison_columns: Sequence[str],
    group_columns: Sequence[str] = ("ConditionId", "Engine", "Estimator"),
) -> pd.DataFrame:
    """Count conclusion changes relative to an authoritative raw decision."""

    if not isinstance(decisions, pd.DataFrame) or decisions.empty:
        return pd.DataFrame()
    required = {raw_column, *comparison_columns}
    missing = required - set(decisions.columns)
    if missing:
        raise ValueError(f"decision sensitivity input missing columns: {sorted(missing)}")
    groups = [column for column in group_columns if column in decisions.columns]
    if not groups:
        raise ValueError("at least one sensitivity group column must be present")
    rows: list[dict[str, Any]] = []
    for group_key, part in decisions.groupby(groups, dropna=False, sort=False):
        keys = group_key if isinstance(group_key, tuple) else (group_key,)
        identity = dict(zip(groups, keys))
        raw = _coerce_bool(part[raw_column])
        for comparison_column in comparison_columns:
            comparison = _coerce_bool(part[comparison_column])
            available = raw.notna() & comparison.notna()
            changed = raw.loc[available].astype(bool) != comparison.loc[available].astype(bool)
            rows.append({
                "SchemaVersion": SCHEMA_VERSION,
                **identity,
                "AuthoritativeDecision": raw_column,
                "ComparedDecision": comparison_column,
                "Rows": int(len(part)),
                "ComparableRows": int(available.sum()),
                "UnavailableRows": int(len(part) - available.sum()),
                "ChangedConclusions": int(changed.sum()),
                "ChangedConclusionRate": float(changed.mean()) if len(changed) else np.nan,
                "Boundary": "Descriptive sensitivity only; the authoritative rule remains the unrounded decision contract.",
            })
    return pd.DataFrame(rows)


def build_operating_characteristics_first_read(
    accounting: pd.DataFrame,
    binary_summary: pd.DataFrame,
    estimation_summary: pd.DataFrame,
    *,
    profile: str,
    min_study_replicates: int = DEFAULT_MIN_STUDY_REPLICATES,
) -> pd.DataFrame:
    """Build a compact, UI-ready interpretation surface for a study bundle."""

    columns = [
        "Priority",
        "Check",
        "Status",
        "Evidence",
        "NextAction",
        "DoNotClaim",
    ]
    if not isinstance(accounting, pd.DataFrame) or accounting.empty:
        return pd.DataFrame([{
            "Priority": 1,
            "Check": "Simulation evidence",
            "Status": "Missing",
            "Evidence": "No condition-level run accounting is available.",
            "NextAction": "Run the smoke profile and retain attempted-run accounting.",
            "DoNotClaim": "Do not infer convergence, power, false-positive rate, or coverage.",
        }], columns=columns)

    def numeric_column(frame: pd.DataFrame, column: str) -> pd.Series:
        source = frame[column] if column in frame.columns else pd.Series(0, index=frame.index)
        return pd.to_numeric(source, errors="coerce").fillna(0)

    attempts = int(numeric_column(accounting, "Attempts").sum())
    converged = int(numeric_column(accounting, "Converged").sum())
    eligible = int(numeric_column(accounting, "AnalysisEligible").sum())
    condition_min_eligible = int(
        numeric_column(accounting, "AnalysisEligible").min()
    )
    study_depth = condition_min_eligible >= int(min_study_replicates)
    decision_tiers = sorted(
        set(binary_summary.get("EvidenceTier", pd.Series(dtype=str)).dropna().astype(str))
    ) if isinstance(binary_summary, pd.DataFrame) else []
    recovery_rows = int(
        pd.to_numeric(estimation_summary.get("RowsIncluded", 0), errors="coerce").fillna(0).sum()
    ) if isinstance(estimation_summary, pd.DataFrame) else 0
    coverage_rows = int(
        pd.to_numeric(estimation_summary.get("CoverageN", 0), errors="coerce").fillna(0).sum()
    ) if isinstance(estimation_summary, pd.DataFrame) else 0

    rows = [
        {
            "Priority": 1,
            "Check": "Evidence scope",
            "Status": "Study-depth screen" if study_depth else "Pilot only",
            "Evidence": (
                f"profile={profile}; minimum eligible decisions per condition={condition_min_eligible}; "
                f"configured study-depth floor={int(min_study_replicates)}."
            ),
            "NextAction": (
                "Inspect Monte Carlo SE, failure patterns, and cross-engine agreement before any validation claim."
                if study_depth else
                "Expand only after the smoke matrix, estimands, runtime, and negative controls are accepted."
            ),
            "DoNotClaim": "A replicate-count label is not evidence that the model or thresholds are valid for a new population.",
        },
        {
            "Priority": 2,
            "Check": "Run accounting",
            "Status": "Review" if converged < attempts or eligible < attempts else "Complete smoke path",
            "Evidence": f"{converged}/{attempts} converged; {eligible}/{attempts} eligible decisions.",
            "NextAction": "Open failure_reasons.csv; never recode failed or unavailable decisions as negative results.",
            "DoNotClaim": "Do not report decision rates without the attempted and eligible denominators.",
        },
        {
            "Priority": 3,
            "Check": "Decision operating characteristics",
            "Status": "Pilot only" if not study_depth else "Review Monte Carlo precision",
            "Evidence": f"decision evidence tiers: {', '.join(decision_tiers) or 'unavailable'}.",
            "NextAction": "Report false-positive rate and power separately by truth, design, decision rule, and engine.",
            "DoNotClaim": "Do not pool null and alternative conditions or call an ineligible sparse cell a negative decision.",
        },
        {
            "Priority": 4,
            "Check": "Parameter recovery and coverage",
            "Status": "Pilot only" if recovery_rows else "Missing",
            "Evidence": f"included recovery rows={recovery_rows}; finite positive-SE coverage rows={coverage_rows}.",
            "NextAction": "Keep anchor-identified absolute errors separate from mean-aligned location errors.",
            "DoNotClaim": "Conditional-Wald coverage among returned rows does not include failed fits, alignment uncertainty, or model selection.",
        },
        {
            "Priority": 5,
            "Check": "Public application surface",
            "Status": "Withheld",
            "Evidence": "Repository validation artifact only; no public estimator or automatic design recommendation is enabled.",
            "NextAction": "Add mfrmr/TAM/immer/sirt adapters and study-depth evidence before considering an in-app evidence button.",
            "DoNotClaim": "Do not translate this pilot into a universal sample-size, anchor-share, or pass/fail rule.",
        },
    ]
    return pd.DataFrame(rows, columns=columns)
