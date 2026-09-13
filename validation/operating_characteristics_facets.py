#!/usr/bin/env python3
"""Replay the operating-characteristics bundle through FACETS 4.5.x.

This adapter deliberately does not generate data.  It consumes the same
byte-retained bundle used by the Python, TAM, sirt, and immer lanes, writes one
auditable FACETS specification per RunId, invokes FACETS in documented batch
mode, and normalizes Table 7 score files to the repository's estimator-neutral
parameter-recovery contract.

Direct numerical agreement is interpreted only for FACETS JMLE versus the
application's matched additive-RSM JMLE.  TAM/sirt MML and exact CMLE remain
different-estimand operating-characteristic comparisons.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from validation.facets_visible_launcher import invoke_facets_batch_no


SCHEMA_VERSION = "mfrm-facets-oc-adapter-v1"
FACET_COLUMNS = ((1, "Person"), (2, "Rater"), (3, "Task"), (4, "Criterion"))
RECOVERY_FACETS = ("Rater", "Task", "Criterion")
REQUIRED_BUNDLE_COLUMNS = {
    "manifest.csv": {"RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed", "Categories"},
    "generated_ratings.csv": {"RunId", "Person", "Rater", "Task", "Criterion", "Score"},
    "generated_facet_truth.csv": {"RunId", "Facet", "Level", "Truth"},
    "generated_anchors.csv": {"RunId", "Facet", "Level", "Anchor"},
}
NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)")
FACETS_VERSION_RE = re.compile(r"Facets 64-bit .*? ([0-9]+(?:\.[0-9]+)+)")
FACETS_FIT_DISPLAY_DECIMALS = 2
TABLE8_REPORT_DECIMALS = 2
FACETS_WINDOWS_SAFE_PATH_CHARS = 220
BIAS_ALPHA = 0.05
BIAS_PRACTICAL_LOGIT = 0.50
BIAS_MIN_COUNT = 5


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def validate_bundle(input_dir: Path) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for filename, required in REQUIRED_BUNDLE_COLUMNS.items():
        path = input_dir / filename
        if not path.is_file():
            raise FileNotFoundError(f"Required operating-characteristics file is missing: {path}")
        table = _read_csv(path)
        missing = required.difference(table.columns)
        if missing:
            raise ValueError(f"{filename} is missing columns: {sorted(missing)}")
        tables[filename] = table

    manifest = tables["manifest.csv"]
    if manifest["RunId"].astype(str).duplicated().any():
        raise ValueError("manifest.csv contains duplicate RunId values")
    manifest_ids = set(manifest["RunId"].astype(str))
    for filename in ("generated_ratings.csv", "generated_facet_truth.csv"):
        observed = set(tables[filename]["RunId"].astype(str))
        if observed != manifest_ids:
            raise ValueError(
                f"{filename} RunIds differ from manifest: "
                f"missing={sorted(manifest_ids - observed)[:3]}, "
                f"extra={sorted(observed - manifest_ids)[:3]}"
            )
    anchor_ids = set(tables["generated_anchors.csv"]["RunId"].astype(str))
    if not anchor_ids.issubset(manifest_ids):
        raise ValueError("generated_anchors.csv contains RunIds absent from manifest")
    return tables


def _validate_label(value: object, *, column: str) -> str:
    text = str(value)
    if any(character in text for character in (",", ";", "\n", "\r", "\t")):
        raise ValueError(f"FACETS label in {column} contains a reserved delimiter: {text!r}")
    if not text:
        raise ValueError(f"FACETS label in {column} is empty")
    return text


def _levels_in_order(ratings: pd.DataFrame, column: str) -> list[str]:
    return list(dict.fromkeys(_validate_label(value, column=column) for value in ratings[column]))


def build_facets_spec(
    manifest_row: pd.Series,
    ratings: pd.DataFrame,
    anchors: pd.DataFrame,
    *,
    score_base: Path,
) -> tuple[str, dict[str, dict[str, int]]]:
    """Build a matched additive-RSM FACETS specification and level maps."""

    if ratings.empty:
        raise ValueError(f"RunId {manifest_row['RunId']} has no rating rows")
    scores = pd.to_numeric(ratings["Score"], errors="raise").astype(int)
    categories = sorted(scores.unique().tolist())
    declared_max = int(manifest_row["Categories"]) - 1
    if categories[0] < 0 or categories[-1] > declared_max:
        raise ValueError(f"Observed scores fall outside 0..{declared_max}: {categories}")
    if declared_max < 1:
        raise ValueError("FACETS RSM requires at least two declared categories")

    level_maps: dict[str, dict[str, int]] = {}
    lines = [
        f"Title=FACETS replay {manifest_row['RunId']}",
        "Facets=4",
        "Noncenter=1",
        "Positive=1",
        "Umean=0,1,6",
        "Convergence=0.5,0.01,0,0",
        "Iterations=0",
        "Xtreme=0.3,0.5",
        "CSV=Tabs",
        "Heading lines=Yes",
        f"Scorefile={score_base.resolve()}",
        "Pt-biserial=Measure",
        "Fair=Mean",
        "Models=",
        f"?,?,?,?,R{declared_max}",
        f"?,?B,?B,?,R{declared_max}",
        "*",
        "Labels=",
    ]

    normalized_anchors = anchors.copy()
    if not normalized_anchors.empty:
        normalized_anchors["Facet"] = normalized_anchors["Facet"].astype(str)
        normalized_anchors["Level"] = normalized_anchors["Level"].astype(str)
        normalized_anchors["Anchor"] = pd.to_numeric(normalized_anchors["Anchor"], errors="raise")

    for facet_number, facet_name in FACET_COLUMNS:
        levels = _levels_in_order(ratings, facet_name)
        mapping = {level: index for index, level in enumerate(levels, start=1)}
        level_maps[facet_name] = mapping
        facet_anchors = normalized_anchors.loc[
            normalized_anchors["Facet"].eq(facet_name)
        ] if not normalized_anchors.empty else normalized_anchors
        anchor_lookup = dict(zip(facet_anchors.get("Level", []), facet_anchors.get("Anchor", [])))
        anchor_flag = ",A" if anchor_lookup else ""
        lines.append(f"{facet_number},{facet_name}{anchor_flag}")
        for level, element_number in mapping.items():
            if level in anchor_lookup:
                lines.append(f"{element_number}={level},{float(anchor_lookup[level]):.12g}")
            else:
                lines.append(f"{element_number}={level}")
        lines.append("*")

    lines.append("Data=")
    for row in ratings.itertuples(index=False):
        encoded = [
            level_maps["Person"][_validate_label(row.Person, column="Person")],
            level_maps["Rater"][_validate_label(row.Rater, column="Rater")],
            level_maps["Task"][_validate_label(row.Task, column="Task")],
            level_maps["Criterion"][_validate_label(row.Criterion, column="Criterion")],
            int(row.Score),
        ]
        lines.append(",".join(map(str, encoded)))
    return "\n".join(lines) + "\n", level_maps


def parse_score_file(path: Path, *, facet_number: int, facet_name: str) -> pd.DataFrame:
    """Parse a tab-delimited FACETS score/measure file."""

    if not path.is_file():
        raise FileNotFoundError(f"FACETS score file is missing: {path}")
    first_line = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()[0]
    expected_prefix = f"{facet_number}\t{facet_name}"
    if not first_line.startswith(expected_prefix):
        raise ValueError(f"Unexpected FACETS score header {first_line!r}; expected {expected_prefix!r}")
    table = pd.read_csv(path, sep="\t", skiprows=1, dtype=str)
    element_number_column = str(facet_number)
    if element_number_column not in table or facet_name not in table:
        raise ValueError(
            f"FACETS score file lacks element columns {element_number_column!r}/{facet_name!r}: "
            f"{list(table.columns)}"
        )
    output = pd.DataFrame({
        "Facet": facet_name,
        "ElementNumber": pd.to_numeric(table[element_number_column], errors="coerce").astype("Int64"),
        "Level": table[facet_name].astype(str),
    })
    fields = {
        "T.Score": "TotalScore",
        "T.Count": "TotalCount",
        "Obs.Avge": "ObservedAverage",
        "FairMAvge": "FairAverage",
        "Measure": "Estimate",
        "S.E.": "SE",
        # FACETS 4.5 score files retain these fit fields at two decimal
        # places even when Umean requests more measure decimals.  The names
        # therefore say Displayed explicitly; they must never be treated as
        # raw fit statistics.
        "InfitMS": "InfitMSDisplayed",
        "InfitZ": "InfitZDisplayed",
        "OutfitMS": "OutfitMSDisplayed",
        "OutfitZ": "OutfitZDisplayed",
        "PtBis": "PtMeasure",
        "PtMeExp": "PtMeasureExpected",
        "Displace": "Displacement",
        "Status": "Status",
    }
    for source, target in fields.items():
        output[target] = pd.to_numeric(table[source], errors="coerce") if source in table else np.nan
    return add_fit_display_contract(output)


def displayed_rounding_interval(
    values: pd.Series,
    *,
    decimals: int = FACETS_FIT_DISPLAY_DECIMALS,
) -> tuple[pd.Series, pd.Series]:
    """Return conservative raw-value bounds implied by displayed rounding."""

    numeric = pd.to_numeric(values, errors="coerce")
    half_unit = 0.5 * 10.0 ** (-int(decimals))
    return numeric - half_unit, numeric + half_unit


def classify_displayed_band(
    displayed: pd.Series,
    *,
    lower: float,
    upper: float,
    decimals: int = FACETS_FIT_DISPLAY_DECIMALS,
) -> pd.Series:
    """Classify a threshold decision without inventing unreported precision.

    A displayed value is ``pass`` only when its entire possible raw interval is
    inside the accepted band, ``flag`` only when the entire interval is outside
    on one side, and otherwise ``boundary_uncertain``.
    """

    raw_lower, raw_upper = displayed_rounding_interval(displayed, decimals=decimals)
    result = pd.Series("boundary_uncertain", index=displayed.index, dtype="string")
    result.loc[raw_lower.ge(lower) & raw_upper.le(upper)] = "pass"
    result.loc[raw_upper.lt(lower) | raw_lower.gt(upper)] = "flag"
    result.loc[raw_lower.isna() | raw_upper.isna()] = pd.NA
    return result


def add_fit_display_contract(table: pd.DataFrame) -> pd.DataFrame:
    """Attach interval and decision columns to FACETS' two-decimal fit output."""

    output = table.copy()
    output["FitPrecisionContract"] = "FACETS_display_2dp_not_raw"
    for stem in ("InfitMS", "OutfitMS"):
        displayed = output[f"{stem}Displayed"]
        lower, upper = displayed_rounding_interval(displayed)
        output[f"{stem}RawLowerBound"] = lower
        output[f"{stem}RawUpperBound"] = upper
        output[f"{stem}Decision0p5To1p5"] = classify_displayed_band(
            displayed,
            lower=0.5,
            upper=1.5,
        )
    for stem in ("InfitZ", "OutfitZ"):
        displayed = output[f"{stem}Displayed"]
        lower, upper = displayed_rounding_interval(displayed)
        output[f"{stem}RawLowerBound"] = lower
        output[f"{stem}RawUpperBound"] = upper
        output[f"{stem}DecisionAbs2"] = classify_displayed_band(
            displayed,
            lower=-2.0,
            upper=2.0,
        )
    return output


def parse_display_token(token: object) -> dict[str, Any]:
    """Parse one FACETS display token without inventing hidden precision."""

    raw = str(token if token is not None else "").strip()
    marker_nonmonotonic = "*" in raw
    marker_anchored = raw.upper().endswith("A")
    marker_extreme = ">" if ">" in raw else ("<" if "<" in raw else "")
    cleaned = raw.replace("*", "").replace(">", "").replace("<", "").strip()
    if marker_anchored:
        cleaned = cleaned[:-1].strip()
    parenthesized = cleaned.startswith("(") and cleaned.endswith(")")
    cleaned = cleaned.strip("() ").replace(" ", "")
    if not cleaned or cleaned.lower() in {"low", "high", "nan", "na"}:
        return {
            "token": raw,
            "value": np.nan,
            "decimals": np.nan,
            "lower": np.nan,
            "upper": np.nan,
            "status": "not_numeric",
            "nonmonotonic_marker": marker_nonmonotonic,
            "anchored_marker": marker_anchored,
            "extreme_marker": marker_extreme,
            "parenthesized": parenthesized,
        }
    # A lone decimal point, or an integer ending in a decimal point, is a
    # fixed-width truncation observed when Umean=...,6 is used for Table 8.
    # It is not interpreted as 1.0 or 0.0.
    if cleaned == "." or re.fullmatch(r"[-+]?\d+\.", cleaned):
        return {
            "token": raw,
            "value": np.nan,
            "decimals": np.nan,
            "lower": np.nan,
            "upper": np.nan,
            "status": "ambiguous_truncated_display",
            "nonmonotonic_marker": marker_nonmonotonic,
            "anchored_marker": marker_anchored,
            "extreme_marker": marker_extreme,
            "parenthesized": parenthesized,
        }
    if not re.fullmatch(r"[-+]?(?:\d+(?:\.\d+)?|\.\d+)", cleaned):
        return {
            "token": raw,
            "value": np.nan,
            "decimals": np.nan,
            "lower": np.nan,
            "upper": np.nan,
            "status": "unparsed",
            "nonmonotonic_marker": marker_nonmonotonic,
            "anchored_marker": marker_anchored,
            "extreme_marker": marker_extreme,
            "parenthesized": parenthesized,
        }
    value = float(cleaned)
    decimals = len(cleaned.rsplit(".", 1)[1]) if "." in cleaned else 0
    half_unit = 0.5 * 10.0 ** (-decimals)
    return {
        "token": raw,
        "value": value,
        "decimals": decimals,
        "lower": value - half_unit,
        "upper": value + half_unit,
        "status": "parsed_display",
        "nonmonotonic_marker": marker_nonmonotonic,
        "anchored_marker": marker_anchored,
        "extreme_marker": marker_extreme,
        "parenthesized": parenthesized,
    }


def _display_fields(prefix: str, token: object) -> dict[str, Any]:
    parsed = parse_display_token(token)
    return {
        f"{prefix}Token": parsed["token"],
        f"{prefix}Displayed": parsed["value"],
        f"{prefix}DisplayDecimals": parsed["decimals"],
        f"{prefix}RawLowerBound": parsed["lower"],
        f"{prefix}RawUpperBound": parsed["upper"],
        f"{prefix}ParseStatus": parsed["status"],
        f"{prefix}NonmonotonicMarker": parsed["nonmonotonic_marker"],
        f"{prefix}AnchoredMarker": parsed["anchored_marker"],
        f"{prefix}ExtremeMarker": parsed["extreme_marker"],
        f"{prefix}Parenthesized": parsed["parenthesized"],
    }


def _numeric_tokens(text: str) -> list[str]:
    return re.findall(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[A*<>])?", text)


def _interval_band_decision(
    lower: pd.Series,
    upper: pd.Series,
    *,
    accepted_lower: float,
    accepted_upper: float,
) -> pd.Series:
    result = pd.Series("boundary_uncertain", index=lower.index, dtype="string")
    result.loc[lower.ge(accepted_lower) & upper.le(accepted_upper)] = "pass"
    result.loc[upper.lt(accepted_lower) | lower.gt(accepted_upper)] = "flag"
    result.loc[lower.isna() | upper.isna()] = pd.NA
    return result


def parse_table8_categories(path: Path) -> pd.DataFrame:
    """Parse numeric Table 8 category rows from a FACETS ASCII report."""

    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    headings = [
        index for index, line in enumerate(lines)
        if re.match(r"^Table 8\.\d+\s+Category Statistics\.", line)
    ]
    if not headings:
        raise ValueError("FACETS report does not contain Table 8 Category Statistics")
    rows: list[dict[str, Any]] = []
    for heading_index in headings:
        table_number = lines[heading_index].split()[1]
        model = ""
        for line in lines[heading_index + 1: heading_index + 8]:
            if line.lstrip().startswith("Model ="):
                model = line.split("=", 1)[1].strip()
                break
        for line in lines[heading_index + 1:]:
            if rows and line.startswith("Table "):
                break
            if not line.startswith("|"):
                if rows and line.startswith("+") and "Mean" in line:
                    break
                continue
            segments = line.strip("|").split("|")
            if len(segments) < 7:
                continue
            data_tokens = segments[0].split()
            if len(data_tokens) != 5 or not re.fullmatch(r"[-+]?\d+", data_tokens[0]):
                continue
            quality_tokens = segments[1].split()
            if len(quality_tokens) != 3:
                continue
            threshold_tokens = _numeric_tokens(segments[2])
            expectation_tokens = _numeric_tokens(segments[3])
            category = int(data_tokens[0])
            row: dict[str, Any] = {
                "TableNumber": table_number,
                "Model": model,
                "Category": category,
                "TotalCount": float(data_tokens[1]),
                "UsedCount": float(data_tokens[2]),
                "UsedPercentDisplayed": float(data_tokens[3].rstrip("%")),
                "CumulativePercentDisplayed": float(data_tokens[4].rstrip("%")),
                "MostProbableFromToken": segments[4].strip(),
                "ThurstoneThresholdToken": segments[5].strip(),
                "PeakProbabilityDisplayed": float(segments[6].strip().rstrip("%")),
                "ResponseCategoryName": segments[7].strip() if len(segments) > 7 else "",
                "SourceReport": str(path.resolve()),
                "SourceLine": line,
                "Table8PrecisionContract": "FACETS_auxiliary_Umean_2_display_not_raw",
            }
            row.update(_display_fields("AverageMeasure", quality_tokens[0]))
            row.update(_display_fields("ExpectedMeasure", quality_tokens[1]))
            row.update(_display_fields("CategoryOutfitMS", quality_tokens[2]))
            row.update(_display_fields("ThresholdMeasure", threshold_tokens[0] if threshold_tokens else ""))
            row.update(_display_fields("ThresholdSE", threshold_tokens[1] if len(threshold_tokens) > 1 else ""))
            row.update(_display_fields("ExpectedCategoryMeasure", expectation_tokens[0] if expectation_tokens else ""))
            row.update(_display_fields("ExpectedMinus0p5Measure", expectation_tokens[1] if len(expectation_tokens) > 1 else ""))
            row.update(_display_fields("MostProbableFrom", segments[4].strip()))
            row.update(_display_fields("ThurstoneThreshold", segments[5].strip()))
            rows.append(row)
        if rows:
            # The registered additive RSM has one scale.  Multiple scales are
            # retained by the outer loop, but duplicate row keys fail below.
            continue
    output = pd.DataFrame(rows)
    if output.empty:
        raise ValueError("FACETS Table 8 headings were found but no category rows parsed")
    if output.duplicated(["TableNumber", "Category"]).any():
        raise ValueError("FACETS Table 8 contains duplicate table/category rows")
    output["CategoryOutfitDecision0p5To1p5"] = _interval_band_decision(
        pd.to_numeric(output["CategoryOutfitMSRawLowerBound"], errors="coerce"),
        pd.to_numeric(output["CategoryOutfitMSRawUpperBound"], errors="coerce"),
        accepted_lower=0.5,
        accepted_upper=1.5,
    )
    output["LowCount10"] = pd.to_numeric(output["UsedCount"], errors="coerce").lt(10)
    output["ThresholdOrderDecision"] = "not_applicable_first_threshold"
    output["AverageMeasureOrderDecision"] = "not_applicable_first_category"
    for _, indices in output.groupby("TableNumber", sort=False).groups.items():
        ordered_indices = output.loc[indices].sort_values("Category").index.tolist()
        threshold_indices = [
            index for index in ordered_indices
            if np.isfinite(_number_or_nan(output.at[index, "ThresholdMeasureDisplayed"]))
        ]
        for previous, current in zip(threshold_indices, threshold_indices[1:]):
            prev_lo = _number_or_nan(output.at[previous, "ThresholdMeasureRawLowerBound"])
            prev_hi = _number_or_nan(output.at[previous, "ThresholdMeasureRawUpperBound"])
            curr_lo = _number_or_nan(output.at[current, "ThresholdMeasureRawLowerBound"])
            curr_hi = _number_or_nan(output.at[current, "ThresholdMeasureRawUpperBound"])
            if prev_hi < curr_lo:
                decision = "ordered"
            elif prev_lo >= curr_hi:
                decision = "disordered"
            else:
                decision = "boundary_uncertain"
            output.at[current, "ThresholdOrderDecision"] = decision
        for previous, current in zip(ordered_indices, ordered_indices[1:]):
            prev_lo = _number_or_nan(output.at[previous, "AverageMeasureRawLowerBound"])
            prev_hi = _number_or_nan(output.at[previous, "AverageMeasureRawUpperBound"])
            curr_lo = _number_or_nan(output.at[current, "AverageMeasureRawLowerBound"])
            curr_hi = _number_or_nan(output.at[current, "AverageMeasureRawUpperBound"])
            if not all(np.isfinite(value) for value in (prev_lo, prev_hi, curr_lo, curr_hi)):
                decision = "unavailable"
            elif prev_hi < curr_lo:
                decision = "ordered"
            elif curr_hi <= prev_lo:
                decision = "disordered"
            else:
                decision = "boundary_uncertain"
            output.at[current, "AverageMeasureOrderDecision"] = decision
    return output


def _number_or_nan(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if np.isfinite(number) else float("nan")


def _level_pattern(levels: Iterable[str]) -> str:
    return "(?:" + "|".join(re.escape(str(level)) for level in sorted(levels, key=lambda x: -len(str(x)))) + ")"


def _parse_table13_elements(
    segment: str,
    *,
    level_maps: dict[str, dict[str, int]],
) -> dict[str, Any]:
    rater_pattern = _level_pattern(level_maps["Rater"])
    task_pattern = _level_pattern(level_maps["Task"])
    numeric = r"[<>]?[-+]?(?:\d+(?:\.\d*)?|\.\d+)[<>]?"
    pattern = re.compile(
        rf"^\s*(?P<sequence>\d+)\s+"
        rf"(?P<rater_number>\d+)\s+(?P<rater>{rater_pattern})\s+(?P<rater_measure>{numeric})\s+"
        rf"(?P<task_number>\d+)\s+(?P<task>{task_pattern})\s+(?P<task_measure>{numeric})\s*$"
    )
    match = pattern.match(segment)
    if not match:
        raise ValueError(f"Could not parse FACETS Table 13 element segment: {segment!r}")
    values = match.groupdict()
    if level_maps["Rater"].get(values["rater"]) != int(values["rater_number"]):
        raise ValueError("FACETS Table 13 Rater number/label mismatch")
    if level_maps["Task"].get(values["task"]) != int(values["task_number"]):
        raise ValueError("FACETS Table 13 Task number/label mismatch")
    output: dict[str, Any] = {
        "Sequence": int(values["sequence"]),
        "RaterElementNumber": int(values["rater_number"]),
        "Rater": values["rater"],
        "TaskElementNumber": int(values["task_number"]),
        "Task": values["task"],
    }
    output.update(_display_fields("RaterMeasure", values["rater_measure"]))
    output.update(_display_fields("TaskMeasure", values["task_measure"]))
    return output


def _holm_adjust(p_values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    result = np.full(len(numeric), np.nan, dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(numeric))
    if not len(finite_indices):
        return pd.Series(result, index=p_values.index)
    values = np.clip(numeric[finite_indices], 0.0, 1.0)
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    adjusted_sorted = np.maximum.accumulate(sorted_values * (len(values) - np.arange(len(values))))
    adjusted = np.empty(len(values), dtype=float)
    adjusted[order] = np.clip(adjusted_sorted, 0.0, 1.0)
    result[finite_indices] = adjusted
    return pd.Series(result, index=p_values.index)


def add_table13_decision_contract(table: pd.DataFrame) -> pd.DataFrame:
    output = table.copy()
    output["FamilySize"] = len(output)
    output["HolmPDisplayed"] = _holm_adjust(output["PDisplayed"])
    raw_half_units = (
        0.5 * np.power(10.0, -pd.to_numeric(output["PDisplayDecimals"], errors="coerce"))
    )
    # Holm's sorted/max transform is at most m-Lipschitz in the largest raw-p
    # perturbation.  This deliberately broad bound also covers rank swaps.
    holm_half = len(output) * raw_half_units.max(skipna=True)
    output["HolmPDisplayUncertainty"] = holm_half
    output["HolmPRawLowerBound"] = (output["HolmPDisplayed"] - holm_half).clip(lower=0.0)
    output["HolmPRawUpperBound"] = (output["HolmPDisplayed"] + holm_half).clip(upper=1.0)
    stat = pd.Series("boundary_uncertain", index=output.index, dtype="string")
    stat.loc[output["HolmPRawUpperBound"].lt(BIAS_ALPHA)] = "significant"
    stat.loc[output["HolmPRawLowerBound"].ge(BIAS_ALPHA)] = "not_significant"
    stat.loc[output["HolmPDisplayed"].isna()] = pd.NA
    output["HolmDecision"] = stat

    bias_lo = pd.to_numeric(output["BiasSizeRawLowerBound"], errors="coerce")
    bias_hi = pd.to_numeric(output["BiasSizeRawUpperBound"], errors="coerce")
    minimum_abs = np.where(
        (bias_lo <= 0) & (bias_hi >= 0),
        0.0,
        np.minimum(np.abs(bias_lo), np.abs(bias_hi)),
    )
    maximum_abs = np.maximum(np.abs(bias_lo), np.abs(bias_hi))
    output["AbsBiasRawLowerBound"] = minimum_abs
    output["AbsBiasRawUpperBound"] = maximum_abs
    practical = pd.Series("boundary_uncertain", index=output.index, dtype="string")
    practical.loc[pd.Series(minimum_abs, index=output.index).ge(BIAS_PRACTICAL_LOGIT)] = "practically_large"
    practical.loc[pd.Series(maximum_abs, index=output.index).lt(BIAS_PRACTICAL_LOGIT)] = "below_practical_threshold"
    practical.loc[bias_lo.isna() | bias_hi.isna()] = pd.NA
    output["PracticalDecision"] = practical
    output["SparseCell"] = pd.to_numeric(output["ObservedCountDisplayed"], errors="coerce").lt(BIAS_MIN_COUNT)
    strong = pd.Series("boundary_uncertain", index=output.index, dtype="string")
    strong.loc[output["SparseCell"]] = "sparse_no_claim"
    strong.loc[
        ~output["SparseCell"]
        & output["HolmDecision"].eq("significant")
        & output["PracticalDecision"].eq("practically_large")
    ] = "flag"
    strong.loc[
        ~output["SparseCell"]
        & (
            output["HolmDecision"].eq("not_significant")
            | output["PracticalDecision"].eq("below_practical_threshold")
        )
    ] = "no_flag"
    strong.loc[output["HolmDecision"].isna() | output["PracticalDecision"].isna()] = pd.NA
    output["StrongDecision"] = strong
    output["BiasDecisionContract"] = "Holm_p_lt_0.05_AND_abs_bias_ge_0.50; n_lt_5_sparse"
    return output


def parse_table13_bias(
    path: Path,
    *,
    level_maps: dict[str, dict[str, int]],
) -> pd.DataFrame:
    """Parse the canonical arranged-by-N Table 13 Rater-by-Task report."""

    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    heading_indices = [
        index for index, line in enumerate(lines)
        if re.match(r"^Table 13\.\d+\.\d+\s+Bias/Interaction Report \(arranged by N\)\.", line)
    ]
    if len(heading_indices) != 1:
        raise ValueError(f"Expected one canonical arranged-by-N Table 13, found {len(heading_indices)}")
    heading_index = heading_indices[0]
    interaction = next(
        (line.strip() for line in lines[heading_index + 1: heading_index + 6] if line.startswith("Bias/Interaction:")),
        "",
    )
    if "2. Rater, 3. Task" not in interaction or "higher score = higher bias measure" not in interaction:
        raise ValueError(f"Unexpected FACETS Table 13 interaction/direction: {interaction!r}")
    rows: list[dict[str, Any]] = []
    for line in lines[heading_index + 1:]:
        if line.startswith("Table 14") or line.startswith("Fixed (all = 0)"):
            break
        if not line.startswith("|"):
            continue
        segments = line.strip("|").split("|")
        if len(segments) != 4:
            continue
        observed = segments[0].split()
        bias = segments[1].split()
        fit = segments[2].split()
        if len(observed) != 4 or len(bias) != 5 or len(fit) != 2:
            continue
        if not NUMBER_RE.fullmatch(observed[0]):
            continue
        try:
            elements = _parse_table13_elements(segments[3], level_maps=level_maps)
        except ValueError:
            if "Mean (Count:" in segments[3] or "S.D." in segments[3]:
                continue
            raise
        row: dict[str, Any] = {
            "FacetA": "Rater",
            "FacetA_Level": elements["Rater"],
            "FacetB": "Task",
            "FacetB_Level": elements["Task"],
            "InteractionDirection": "Bias=plus; higher score = positive bias measure",
            "SourceReport": str(path.resolve()),
            "SourceLine": line,
            "Table13Arrangement": "N_canonical",
            "ResidualFitInterpretation": "descriptive_only; Table 13 fit has no usual chi-square/ZSTD properties",
            **elements,
        }
        for prefix, token in (
            ("ObservedScore", observed[0]),
            ("ExpectedScore", observed[1]),
            ("ObservedCount", observed[2]),
            ("ObservedMinusExpectedAverage", observed[3]),
            ("BiasSize", bias[0]),
            ("BiasSE", bias[1]),
            ("TStatistic", bias[2]),
            ("DegreesOfFreedom", bias[3]),
            ("P", bias[4]),
            ("ResidualInfitMS", fit[0]),
            ("ResidualOutfitMS", fit[1]),
        ):
            row.update(_display_fields(prefix, token))
        rows.append(row)
    output = pd.DataFrame(rows)
    if output.empty:
        raise ValueError("FACETS canonical Table 13 was found but no interaction rows parsed")
    if output.duplicated(["FacetA_Level", "FacetB_Level"]).any():
        raise ValueError("FACETS canonical Table 13 contains duplicate Rater-by-Task rows")
    return add_table13_decision_contract(output)


def parse_iteration_report(path: Path) -> dict[str, Any]:
    """Extract the final JMLE convergence row and reported FACETS version."""

    text = path.read_text(encoding="utf-8-sig", errors="replace")
    version_match = FACETS_VERSION_RE.search(text)
    final: dict[str, Any] = {
        "FacetsVersion": version_match.group(1) if version_match else "unknown",
        "Iterations": np.nan,
        "MaxScoreResidual": np.nan,
        "MaxScoreResidualPercent": np.nan,
        "MaxCategoryResidual": np.nan,
        "MaxElementLogitChange": np.nan,
        "MaxStepLogitChange": np.nan,
        "Converged": False,
        "SubsetConnected": "Subset connection O.K." in text,
    }
    for line in text.splitlines():
        if not re.match(r"\s*\|\s*JMLE\s+", line):
            continue
        numbers = [float(value) for value in NUMBER_RE.findall(line)]
        if len(numbers) < 6:
            continue
        final.update({
            "Iterations": int(numbers[0]),
            "MaxScoreResidual": numbers[1],
            "MaxScoreResidualPercent": numbers[2],
            "MaxCategoryResidual": numbers[3],
            "MaxElementLogitChange": numbers[4],
            "MaxStepLogitChange": numbers[5],
        })
    residual = float(final["MaxScoreResidual"])
    change = float(final["MaxElementLogitChange"])
    final["Converged"] = bool(
        np.isfinite(residual)
        and np.isfinite(change)
        and abs(residual) <= 0.5
        and abs(change) <= 0.01
    )
    return final


def _facets_command(
    facets_exe: Path,
    spec_path: Path,
    report_path: Path,
    extra_specs: Iterable[str] | None = None,
    *,
    batch_value: str = "YES",
) -> list[str]:
    return [
        str(facets_exe),
        f"BATCH={batch_value.upper()}",
        str(spec_path),
        str(report_path),
        *[str(specification) for specification in (extra_specs or ())],
    ]


def _facets_specification_value(
    spec_path: Path,
    name: str,
    extra_specs: Iterable[str],
) -> str | None:
    pattern = re.compile(rf"^\s*{re.escape(name)}\s*=\s*(.*?)\s*$", re.IGNORECASE)
    value: str | None = None
    for raw_line in spec_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.split(";", 1)[0]
        match = pattern.match(line)
        if match:
            value = match.group(1).strip().strip('"').strip("'")
    for raw_specification in extra_specs:
        match = pattern.match(str(raw_specification))
        if match:
            value = match.group(1).strip().strip('"').strip("'")
    return value


def _facets_readiness_paths(
    spec_path: Path,
    report_path: Path,
    extra_specs: Iterable[str],
) -> tuple[Path, ...]:
    """Return the report and every expected per-facet Scorefile output."""

    paths = [report_path.resolve()]
    score_value = _facets_specification_value(spec_path, "Scorefile", extra_specs)
    facet_value = _facets_specification_value(spec_path, "Facets", extra_specs)
    if not score_value:
        return tuple(paths)
    if facet_value is None:
        raise ValueError("FACETS Scorefile readiness requires a Facets= specification")
    try:
        facet_count = int(facet_value.split(",", 1)[0].strip())
    except ValueError as exc:
        raise ValueError(f"Invalid FACETS facet count for readiness: {facet_value!r}") from exc
    if facet_count < 1:
        raise ValueError(f"Invalid FACETS facet count for readiness: {facet_count}")
    score_base = Path(score_value)
    if not score_base.is_absolute():
        score_base = spec_path.parent / score_base
    score_base = score_base.resolve()
    for facet_number in range(1, facet_count + 1):
        if score_base.suffix:
            score_path = score_base.with_name(
                f"{score_base.stem}.{facet_number}{score_base.suffix}"
            )
        else:
            score_path = score_base.with_name(f"{score_base.name}.{facet_number}.txt")
        paths.append(score_path)
    return tuple(paths)


def _facets_path_budget_violations(
    paths: Iterable[Path],
    *,
    maximum_characters: int = FACETS_WINDOWS_SAFE_PATH_CHARS,
) -> tuple[tuple[Path, int], ...]:
    """Return paths outside the qualified legacy-Windows FACETS budget."""

    if maximum_characters < 1:
        raise ValueError("FACETS path budget must be positive")
    resolved = tuple(Path(path).resolve() for path in paths)
    return tuple(
        (path, len(str(path)))
        for path in resolved
        if len(str(path)) > maximum_characters
    )


def invoke_facets(
    facets_exe: Path,
    spec_path: Path,
    report_path: Path,
    *,
    timeout_seconds: float,
    extra_specs: Iterable[str] | None = None,
    batch_mode: str | None = None,
) -> subprocess.CompletedProcess[str]:
    """Invoke FACETS through a qualified hidden or visible batch route.

    Windows defaults to visible ``BATCH=NO`` because FACETS 4.5.0 and MINIFAC
    4.5.1 were observed to null-dereference their shared Xojo runtime whenever
    the analysis window was hidden.  Set ``MFRM_FACETS_BATCH_MODE=YES`` or pass
    ``batch_mode="YES"`` only on a separately qualified host.
    """

    specifications = tuple(str(value) for value in (extra_specs or ()))
    selected_mode = (batch_mode or os.environ.get("MFRM_FACETS_BATCH_MODE") or "NO")
    selected_mode = selected_mode.strip().upper()
    if selected_mode not in {"YES", "NO"}:
        raise ValueError(f"Unsupported FACETS batch mode: {selected_mode!r}")
    readiness_paths = _facets_readiness_paths(
        spec_path, report_path, specifications
    )
    if os.name == "nt":
        violations = _facets_path_budget_violations(
            (spec_path, *readiness_paths)
        )
        if violations:
            detail = "; ".join(f"{length}: {path}" for path, length in violations)
            raise OSError(
                "FACETS legacy Windows path budget exceeded "
                f"({FACETS_WINDOWS_SAFE_PATH_CHARS} characters): {detail}"
            )
    command = _facets_command(
        facets_exe,
        spec_path,
        report_path,
        specifications,
        batch_value=selected_mode,
    )
    if selected_mode == "YES":
        return subprocess.run(
            command,
            cwd=spec_path.parent,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
    invocation = invoke_facets_batch_no(
        command,
        cwd=spec_path.parent,
        readiness_paths=readiness_paths,
        timeout_seconds=timeout_seconds,
        stable_seconds=2.0,
        hide_window=False,
    )
    completed = invocation.completed
    launcher_evidence = json.dumps(
        {
            "mode": "BATCH=NO_VISIBLE",
            "stable_seconds": invocation.stable_seconds,
            "windows_closed": invocation.windows_closed,
            "forced_termination": invocation.forced_termination,
            "readiness_paths": [str(path) for path in invocation.readiness_paths],
        },
        sort_keys=True,
    )
    completed.stderr = "\n".join(
        part for part in (completed.stderr or "", f"MFRM_FACETS_LAUNCH={launcher_evidence}")
        if part
    )
    return completed


def align_recovery(
    estimates: pd.DataFrame,
    truth: pd.DataFrame,
    anchors: pd.DataFrame,
    manifest_row: pd.Series,
) -> pd.DataFrame:
    truth_subset = truth.loc[
        truth["Facet"].astype(str).isin(RECOVERY_FACETS),
        ["Facet", "Level", "Truth"],
    ].copy()
    truth_subset["Facet"] = truth_subset["Facet"].astype(str)
    truth_subset["Level"] = truth_subset["Level"].astype(str)
    merged = truth_subset.merge(estimates, on=["Facet", "Level"], how="left")
    for column in ("Truth", "Estimate", "SE", "Status"):
        merged[column] = pd.to_numeric(merged.get(column, np.nan), errors="coerce")
    merged["RawError"] = merged["Estimate"] - merged["Truth"]
    merged["EstimateAligned"] = np.nan
    merged["TruthAligned"] = merged["Truth"]
    merged["ErrorAligned"] = np.nan
    merged["ComparisonScale"] = "mean_aligned_location"

    anchored_facets = set(anchors["Facet"].astype(str)) if not anchors.empty else set()
    for facet, indices in merged.groupby("Facet", sort=False).groups.items():
        subset = merged.loc[indices]
        finite = subset["RawError"].notna()
        if facet in anchored_facets:
            merged.loc[indices, "EstimateAligned"] = subset["Estimate"]
            merged.loc[indices, "ErrorAligned"] = subset["RawError"]
            merged.loc[indices, "ComparisonScale"] = "anchor_identified_absolute"
        else:
            shift = float(subset.loc[finite, "RawError"].mean()) if finite.any() else np.nan
            merged.loc[indices, "EstimateAligned"] = subset["Estimate"] - shift
            merged.loc[indices, "ErrorAligned"] = merged.loc[indices, "EstimateAligned"] - subset["Truth"]

    anchor_keys = set(map(tuple, anchors[["Facet", "Level"]].astype(str).to_numpy())) if not anchors.empty else set()
    merged["Anchored"] = [
        (str(facet), str(level)) in anchor_keys
        for facet, level in zip(merged["Facet"], merged["Level"])
    ]
    merged["CoverageEligible"] = (
        ~merged["Anchored"]
        & merged["SE"].gt(0)
        & merged["SE"].notna()
        & merged["ErrorAligned"].notna()
    )
    merged["Covered95"] = np.where(
        merged["CoverageEligible"],
        merged["ErrorAligned"].abs() <= 1.96 * merged["SE"],
        np.nan,
    )
    merged["IncludedInSummary"] = merged["ErrorAligned"].notna()
    merged["Engine"] = "FACETS"
    merged["Estimator"] = "JMLE"
    merged["Mode"] = "FACETS_4_5_STANDARD_JMLE"
    merged["ParameterType"] = merged["Facet"]
    merged["CoverageMethod"] = "FACETS conditional-Wald diagnostic"
    identity = ["RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed"]
    for column in reversed(identity):
        merged.insert(0, column, manifest_row[column])
    ordered = identity + ["Engine", "Estimator", "Mode"]
    return merged[ordered + [column for column in merged.columns if column not in ordered]]


def fit_one_run(
    manifest_row: pd.Series,
    *,
    ratings: pd.DataFrame,
    truth: pd.DataFrame,
    anchors: pd.DataFrame,
    facets_exe: Path,
    work_root: Path,
    timeout_seconds: float,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    run_id = str(manifest_row["RunId"])
    run_dir = work_root / slugify(run_id)
    run_dir.mkdir(parents=True, exist_ok=False)
    spec_path = run_dir / "analysis.txt"
    report_path = run_dir / "report.out.txt"
    score_base = run_dir / "scores.txt"
    spec, level_maps = build_facets_spec(manifest_row, ratings, anchors, score_base=score_base)
    spec_path.write_text(spec, encoding="utf-8", newline="\n")

    started = time.perf_counter()
    base = {column: manifest_row[column] for column in manifest_row.index}
    run: dict[str, Any] = {
        **base,
        "Engine": "FACETS",
        "Estimator": "JMLE",
        "Mode": "FACETS_4_5_STANDARD_JMLE",
        "FitReturned": False,
        "Converged": False,
        "AnalysisEligible": False,
        "StrictAnalysisEligible": False,
        "FailureStage": "invoke",
        "FailureReason": "",
        "ExitCode": np.nan,
        "ElapsedSeconds": np.nan,
        "Rows": int(len(ratings)),
        "FiniteMainFacetMeasures": 0,
        "ExtremeMainFacetElements": 0,
        "UnmeasurableMainFacetElements": 0,
        "SpecificationSHA256": sha256_file(spec_path),
        "ReportSHA256": "",
        "Table8ReportReturned": False,
        "Table8Parsed": False,
        "Table8Categories": 0,
        "Table8CountAuditPassed": False,
        "Table8FailureReason": "",
        "Table8ReportSHA256": "",
        "Table13Parsed": False,
        "Table13Cells": 0,
        "Table13CompleteFamily": False,
        "Table13FailureReason": "",
    }
    try:
        completed = invoke_facets(
            facets_exe,
            spec_path,
            report_path,
            timeout_seconds=timeout_seconds,
        )
        run["ExitCode"] = int(completed.returncode)
        (run_dir / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
        (run_dir / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"FACETS exited with code {completed.returncode}")
        if not report_path.is_file():
            raise FileNotFoundError("FACETS returned zero but did not write the report")
        run["ReportSHA256"] = sha256_file(report_path)
        run.update(parse_iteration_report(report_path))
        run["FitReturned"] = True
        run["FailureStage"] = "parse"

        score_parts = []
        for facet_number, facet_name in FACET_COLUMNS:
            score_path = run_dir / f"scores.{facet_number}.txt"
            score_parts.append(
                parse_score_file(score_path, facet_number=facet_number, facet_name=facet_name)
            )
        scores = pd.concat(score_parts, ignore_index=True)
        main_scores = scores[scores["Facet"].isin(RECOVERY_FACETS)].copy()
        status = pd.to_numeric(main_scores["Status"], errors="coerce")
        run["FiniteMainFacetMeasures"] = int(pd.to_numeric(main_scores["Estimate"], errors="coerce").notna().sum())
        run["ExtremeMainFacetElements"] = int(status.isin([-2, -3]).sum())
        run["UnmeasurableMainFacetElements"] = int(status.isin([0, -4, -5]).sum() + status.isna().sum())
        expected_main = int(sum(truth["Facet"].astype(str).eq(facet).sum() for facet in RECOVERY_FACETS))
        finite_complete = run["FiniteMainFacetMeasures"] == expected_main
        run["AnalysisEligible"] = bool(run["Converged"] and finite_complete)
        run["StrictAnalysisEligible"] = bool(
            run["AnalysisEligible"]
            and run["ExtremeMainFacetElements"] == 0
            and run["UnmeasurableMainFacetElements"] == 0
        )
        run["FailureStage"] = ""
        recovery = align_recovery(scores, truth, anchors, manifest_row)
        recovery["IncludedInSummary"] &= bool(run["AnalysisEligible"])

        identity = [
            column for column in (
                "RunId", "ConditionId", "Design", "TruthBias", "TruthPositive",
                "BiasCell", "Replicate", "Seed",
            )
            if column in manifest_row.index
        ]
        bias = pd.DataFrame()
        try:
            bias = parse_table13_bias(report_path, level_maps=level_maps)
            for column in reversed(identity):
                bias.insert(0, column, manifest_row[column])
            bias.insert(len(identity), "Engine", "FACETS")
            bias.insert(len(identity) + 1, "Estimator", "JMLE")
            bias.insert(len(identity) + 2, "Mode", "FACETS_4_5_TABLE13_BIAS_PLUS")
            expected_bias_cells = len(level_maps["Rater"]) * len(level_maps["Task"])
            run["Table13Parsed"] = True
            run["Table13Cells"] = int(len(bias))
            run["Table13CompleteFamily"] = bool(
                len(bias) == expected_bias_cells
                and not bias.duplicated(["FacetA_Level", "FacetB_Level"]).any()
            )
        except Exception as bias_exc:  # keep JMLE measurement evidence usable
            run["Table13FailureReason"] = f"{type(bias_exc).__name__}: {bias_exc}"

        categories = pd.DataFrame()
        table8_report_path = run_dir / "table8_u2_report.out.txt"
        table8_score_base = run_dir / "table8_u2_scores.txt"
        try:
            table8_completed = invoke_facets(
                facets_exe,
                spec_path,
                table8_report_path,
                timeout_seconds=timeout_seconds,
                extra_specs=(
                    f"Umean=0,1,{TABLE8_REPORT_DECIMALS}",
                    f"Scorefile={table8_score_base.resolve()}",
                ),
            )
            (run_dir / "table8_stdout.txt").write_text(
                table8_completed.stdout or "", encoding="utf-8"
            )
            (run_dir / "table8_stderr.txt").write_text(
                table8_completed.stderr or "", encoding="utf-8"
            )
            if table8_completed.returncode != 0:
                raise RuntimeError(f"FACETS Table 8 reporting pass exited with code {table8_completed.returncode}")
            if not table8_report_path.is_file():
                raise FileNotFoundError("FACETS Table 8 reporting pass did not write its report")
            run["Table8ReportReturned"] = True
            run["Table8ReportSHA256"] = sha256_file(table8_report_path)
            categories = parse_table8_categories(table8_report_path)
            primary_categories = parse_table8_categories(report_path)
            primary_prefixes = (
                "AverageMeasure", "ExpectedMeasure", "ThresholdMeasure", "ThresholdSE",
                "ExpectedCategoryMeasure", "ExpectedMinus0p5Measure",
                "MostProbableFrom", "ThurstoneThreshold",
            )
            primary_columns = [
                column for column in primary_categories.columns
                if column in {"TableNumber", "Category"}
                or column.startswith(primary_prefixes)
            ]
            primary_categories = primary_categories[primary_columns].rename(columns={
                column: f"PrimaryU6{column}"
                for column in primary_columns
                if column not in {"TableNumber", "Category"}
            })
            categories = categories.merge(
                primary_categories,
                on=["TableNumber", "Category"],
                how="left",
                validate="one_to_one",
            )
            categories["DualPassPrecisionContract"] = (
                "U6 primary for measures/thresholds; U2 auxiliary for category Outfit"
            )
            for column in reversed(identity):
                categories.insert(0, column, manifest_row[column])
            categories.insert(len(identity), "Engine", "FACETS")
            categories.insert(len(identity) + 1, "Estimator", "JMLE")
            categories.insert(len(identity) + 2, "Mode", "FACETS_4_5_TABLE8_U2_DISPLAY")
            raw_counts = ratings.groupby("Score", observed=False).size().to_dict()
            categories["GeneratedRawCount"] = categories["Category"].map(raw_counts).fillna(0).astype(int)
            categories["TotalCountMatchesGenerated"] = (
                pd.to_numeric(categories["TotalCount"], errors="coerce")
                .eq(categories["GeneratedRawCount"])
            )
            categories["UsedCountNotAboveTotal"] = (
                pd.to_numeric(categories["UsedCount"], errors="coerce")
                .le(pd.to_numeric(categories["TotalCount"], errors="coerce"))
            )
            run["Table8Parsed"] = True
            run["Table8Categories"] = int(len(categories))
            run["Table8CountAuditPassed"] = bool(
                categories["TotalCountMatchesGenerated"].all()
                and categories["UsedCountNotAboveTotal"].all()
                and int(pd.to_numeric(categories["TotalCount"], errors="coerce").sum()) == len(ratings)
            )
        except Exception as table8_exc:  # measurement/bias pass remains auditable
            run["Table8FailureReason"] = f"{type(table8_exc).__name__}: {table8_exc}"
        return run, recovery, categories, bias
    except Exception as exc:  # retain a complete attempt ledger
        run["FailureReason"] = f"{type(exc).__name__}: {exc}"
        return run, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    finally:
        run["ElapsedSeconds"] = time.perf_counter() - started


def _normalize_existing_recovery(
    input_dir: Path,
    *,
    run_scope: pd.DataFrame | None = None,
) -> list[pd.DataFrame]:
    sources = [
        ("parameter_recovery.csv", None, "matched_jmle"),
        ("mfrmr_parameter_recovery.csv", "MFRMR_JML_STRICT", "independent_jmle_sensitivity"),
        ("python_mml_parameter_recovery.csv", None, "native_marginal_sensitivity"),
        ("python_cmle_parameter_recovery.csv", None, "conditional_fixed_score"),
        ("tam_facet_recovery.csv", "TAM_MML_Q61_PRIMARY", "marginal_sensitivity"),
        ("sirt_rater_recovery.csv", "SIRT_MML_Q61_PRIMARY", "marginal_sensitivity"),
    ]
    frames: list[pd.DataFrame] = []
    for filename, primary_mode, comparison_class in sources:
        path = input_dir / filename
        if not path.is_file():
            continue
        frame = _read_csv(path)
        if primary_mode is not None and "Mode" in frame:
            frame = frame[frame["Mode"].astype(str).eq(primary_mode)].copy()
        if frame.empty or "ErrorAligned" not in frame:
            continue
        if "Mode" not in frame:
            frame["Mode"] = frame["Estimator"].astype(str)
        if "Facet" not in frame and filename == "sirt_rater_recovery.csv":
            frame["Facet"] = "Rater"
        frame["ComparisonClass"] = comparison_class
        if filename == "parameter_recovery.csv" and run_scope is not None:
            scope = run_scope[["RunId", "ComparisonClass", "DirectComparisonEligible"]].drop_duplicates("RunId")
            frame = frame.drop(columns=["ComparisonClass"], errors="ignore").merge(
                scope,
                on="RunId",
                how="left",
            )
        frames.append(frame)
    return frames


def estimator_bias_summary(recovery: pd.DataFrame) -> pd.DataFrame:
    selected = recovery.loc[
        recovery["IncludedInSummary"].fillna(False).astype(bool)
        & pd.to_numeric(recovery["ErrorAligned"], errors="coerce").notna()
    ].copy()
    selected["ErrorAligned"] = pd.to_numeric(selected["ErrorAligned"], errors="coerce")
    selected["SE"] = pd.to_numeric(selected.get("SE", np.nan), errors="coerce")
    if "CoverageEligible" not in selected:
        selected["CoverageEligible"] = selected["SE"].gt(0)
    if "Covered95" not in selected:
        selected["Covered95"] = np.where(
            selected["CoverageEligible"].fillna(False).astype(bool),
            selected["ErrorAligned"].abs() <= 1.96 * selected["SE"],
            np.nan,
        )
    selected["SquaredError"] = selected["ErrorAligned"].pow(2)
    selected["AbsoluteError"] = selected["ErrorAligned"].abs()

    rows: list[dict[str, Any]] = []
    group_columns = ["Engine", "Estimator", "Mode", "ComparisonClass", "Facet", "Design", "TruthBias"]
    for keys, group in selected.groupby(group_columns, dropna=False, sort=True):
        eligible = group[pd.Series(group["CoverageEligible"], index=group.index).fillna(False).astype(bool)]
        rows.append({
            **dict(zip(group_columns, keys)),
            "Runs": int(group["RunId"].nunique()),
            "Parameters": int(len(group)),
            "MeanError": float(group["ErrorAligned"].mean()),
            "MAE": float(group["AbsoluteError"].mean()),
            "RMSE": float(math.sqrt(group["SquaredError"].mean())),
            "CoverageEligible": int(len(eligible)),
            "Coverage95": float(pd.to_numeric(eligible["Covered95"], errors="coerce").mean()) if len(eligible) else np.nan,
        })
    return pd.DataFrame(rows)


def jmle_agreement(facets: pd.DataFrame, python: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    keys = ["RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Facet", "Level"]
    inclusion_column = (
        "IncludedInDirectComparison"
        if "IncludedInDirectComparison" in facets
        else "IncludedInSummary"
    )
    left = facets[keys + ["EstimateAligned", "SE", inclusion_column]].rename(columns={
        "EstimateAligned": "FACETSEstimate",
        "SE": "FACETSSE",
        inclusion_column: "FACETSIncluded",
    })
    right = python[keys + ["EstimateAligned", "SE", "IncludedInSummary"]].rename(columns={
        "EstimateAligned": "PythonEstimate",
        "SE": "PythonSE",
        "IncludedInSummary": "PythonIncluded",
    })
    pairs = left.merge(right, on=keys, how="inner")
    pairs = pairs[
        pairs["FACETSIncluded"].fillna(False).astype(bool)
        & pairs["PythonIncluded"].fillna(False).astype(bool)
    ].copy()
    pairs["EstimateDifference"] = pd.to_numeric(pairs["FACETSEstimate"], errors="coerce") - pd.to_numeric(pairs["PythonEstimate"], errors="coerce")
    pairs["AbsoluteDifference"] = pairs["EstimateDifference"].abs()
    rows: list[dict[str, Any]] = []
    group_columns = ["Facet", "Design", "TruthBias"]
    for keys_value, group in pairs.groupby(group_columns, dropna=False, sort=True):
        finite = group[["FACETSEstimate", "PythonEstimate"]].apply(pd.to_numeric, errors="coerce").dropna()
        rows.append({
            **dict(zip(group_columns, keys_value)),
            "Runs": int(group["RunId"].nunique()),
            "Parameters": int(len(group)),
            "MeanDifference": float(group["EstimateDifference"].mean()),
            "MAE": float(group["AbsoluteDifference"].mean()),
            "MaxAbsoluteDifference": float(group["AbsoluteDifference"].max()),
            "Pearson": float(finite.corr(method="pearson").iloc[0, 1]) if len(finite) > 1 else np.nan,
            "Spearman": float(finite.corr(method="spearman").iloc[0, 1]) if len(finite) > 1 else np.nan,
        })
    return pairs, pd.DataFrame(rows)


def _python_bias_decisions(input_dir: Path) -> pd.DataFrame:
    path = input_dir / "bias_decision_stability.csv"
    if not path.is_file():
        return pd.DataFrame()
    source = _read_csv(path)
    required = {
        "RunId", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
        "Statistic", "RawValue",
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Python bias decision stability is missing columns: {sorted(missing)}")
    source = source[source["Statistic"].isin(["p_holm", "AbsBias"])].copy()
    keys = ["RunId", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level"]
    if source.duplicated(keys + ["Statistic"]).any():
        raise ValueError("Python bias decision stability has duplicate cell/statistic rows")
    output = source.pivot(index=keys, columns="Statistic", values="RawValue").reset_index()
    output.columns.name = None
    if not {"p_holm", "AbsBias"}.issubset(output.columns):
        raise ValueError("Python bias decision stability does not contain both p_holm and AbsBias")
    output = output.rename(columns={"p_holm": "PythonHolmP", "AbsBias": "PythonAbsBias"})
    output["PythonHolmDecision"] = np.where(
        pd.to_numeric(output["PythonHolmP"], errors="coerce").lt(BIAS_ALPHA),
        "significant",
        "not_significant",
    )
    output["PythonPracticalDecision"] = np.where(
        pd.to_numeric(output["PythonAbsBias"], errors="coerce").ge(BIAS_PRACTICAL_LOGIT),
        "practically_large",
        "below_practical_threshold",
    )
    output["PythonStrongDecision"] = np.where(
        output["PythonHolmDecision"].eq("significant")
        & output["PythonPracticalDecision"].eq("practically_large"),
        "flag",
        "no_flag",
    )
    runs_path = input_dir / "runs.csv"
    if runs_path.is_file():
        runs = _read_csv(runs_path)
        run_columns = [
            column for column in ("RunId", "FitReturned", "Converged", "AnalysisEligible")
            if column in runs
        ]
        if "RunId" in run_columns:
            run_scope = runs[run_columns].drop_duplicates("RunId").rename(columns={
                "FitReturned": "PythonRunFitReturned",
                "Converged": "PythonRunConverged",
                "AnalysisEligible": "PythonRunAnalysisEligible",
            })
            output = output.merge(run_scope, on="RunId", how="left", validate="many_to_one")
    return output


def facets_python_bias_agreement(
    facets_bias: pd.DataFrame,
    input_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    python = _python_bias_decisions(input_dir)
    if facets_bias.empty or python.empty:
        return pd.DataFrame(), pd.DataFrame()
    keys = ["RunId", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level"]
    pairs = facets_bias.merge(python, on=keys, how="inner", validate="one_to_one")
    pairs["FACETSAbsBiasDisplayed"] = pd.to_numeric(
        pairs["BiasSizeDisplayed"], errors="coerce"
    ).abs()
    pairs["AbsBiasDifference"] = pairs["FACETSAbsBiasDisplayed"] - pd.to_numeric(
        pairs["PythonAbsBias"], errors="coerce"
    )
    facets_run_eligible = pairs.get(
        "DirectComparisonEligible", pd.Series(False, index=pairs.index)
    ).fillna(False).astype(bool)
    python_run_eligible = pairs.get(
        "PythonRunAnalysisEligible", pd.Series(False, index=pairs.index)
    ).fillna(False).astype(bool)
    pairs["ComparisonEligible"] = (
        pairs["StrongDecision"].isin(["flag", "no_flag"])
        & pairs["PythonStrongDecision"].isin(["flag", "no_flag"])
        & ~pairs["SparseCell"].fillna(True).astype(bool)
        & facets_run_eligible
        & python_run_eligible
    )
    pairs["StrongDecisionAgreement"] = np.where(
        pairs["ComparisonEligible"],
        pairs["StrongDecision"].eq(pairs["PythonStrongDecision"]),
        np.nan,
    )
    pairs["FocalCell"] = (
        pairs["FacetA_Level"].astype(str)
        + " x "
        + pairs["FacetB_Level"].astype(str)
    ).eq(pairs.get("BiasCell", ""))

    focal_path = input_dir / "bias_decisions.csv"
    if focal_path.is_file():
        focal = _read_csv(focal_path)
        focal_columns = [
            column for column in (
                "RunId", "BiasEstimate", "BiasSE", "p_holm", "AbsBias",
                "DecisionStrongRaw",
            )
            if column in focal
        ]
        if "RunId" in focal_columns:
            focal = focal[focal_columns].drop_duplicates("RunId").rename(columns={
                "BiasEstimate": "PythonFocalBiasEstimate",
                "BiasSE": "PythonFocalBiasSE",
                "p_holm": "PythonFocalHolmP",
                "AbsBias": "PythonFocalAbsBias",
                "DecisionStrongRaw": "PythonFocalStrongFlag",
            })
            pairs = pairs.merge(focal, on="RunId", how="left", validate="many_to_one")

    rows: list[dict[str, Any]] = []
    for keys_value, group in pairs.groupby(["Design", "TruthBias"], dropna=False, sort=True):
        eligible = group[group["ComparisonEligible"]].copy()
        focal = group[group["FocalCell"]].copy()
        focal_eligible = focal[focal["ComparisonEligible"]]
        rows.append({
            "Design": keys_value[0],
            "TruthBias": keys_value[1],
            "RunIds": int(group["RunId"].nunique()),
            "MatchedCells": int(len(group)),
            "ComparisonEligibleCells": int(len(eligible)),
            "FACETSStrongFlags": int(eligible["StrongDecision"].eq("flag").sum()),
            "PythonStrongFlags": int(eligible["PythonStrongDecision"].eq("flag").sum()),
            "StrongDecisionAgreementCells": int(
                pd.to_numeric(eligible["StrongDecisionAgreement"], errors="coerce").sum()
            ),
            "StrongDecisionAgreementRate": float(
                pd.to_numeric(eligible["StrongDecisionAgreement"], errors="coerce").mean()
            ) if len(eligible) else np.nan,
            "BoundaryUncertainCells": int(group["StrongDecision"].eq("boundary_uncertain").sum()),
            "SparseCells": int(group["SparseCell"].fillna(False).astype(bool).sum()),
            "FocalCells": int(len(focal)),
            "FocalEligibleCells": int(len(focal_eligible)),
            "FocalStrongDecisionAgreementRate": float(
                pd.to_numeric(focal_eligible["StrongDecisionAgreement"], errors="coerce").mean()
            ) if len(focal_eligible) else np.nan,
            "ClaimLimit": (
                "parser and decision-mapping qualification; grouped rates are descriptive "
                "and not confirmatory false-positive or power claims"
            ),
        })
    return pairs, pd.DataFrame(rows)


def table13_expected_cell_audit(
    manifest: pd.DataFrame,
    ratings: pd.DataFrame,
    table13: pd.DataFrame,
) -> pd.DataFrame:
    """Account for every expected Rater-by-Task cell, including omissions."""

    rows: list[dict[str, Any]] = []
    identity_columns = [
        column for column in (
            "RunId", "ConditionId", "Design", "TruthBias", "TruthPositive",
            "BiasCell", "Categories", "Replicate", "Seed",
        )
        if column in manifest
    ]
    for manifest_row in manifest.itertuples(index=False):
        run_id = str(manifest_row.RunId)
        run_ratings = ratings[ratings["RunId"].astype(str).eq(run_id)].copy()
        raters = list(dict.fromkeys(run_ratings["Rater"].astype(str)))
        tasks = list(dict.fromkeys(run_ratings["Task"].astype(str)))
        grouped = run_ratings.groupby(["Rater", "Task"], observed=False).agg(
            GeneratedObservationCount=("Score", "size"),
            GeneratedObservedScore=("Score", "sum"),
        )
        identity = {column: getattr(manifest_row, column) for column in identity_columns}
        for rater in raters:
            for task in tasks:
                if (rater, task) in grouped.index:
                    generated = grouped.loc[(rater, task)]
                    count = int(generated["GeneratedObservationCount"])
                    score = float(generated["GeneratedObservedScore"])
                else:
                    count = 0
                    score = 0.0
                rows.append({
                    **identity,
                    "FacetA": "Rater",
                    "FacetA_Level": rater,
                    "FacetB": "Task",
                    "FacetB_Level": task,
                    "GeneratedObservationCount": count,
                    "GeneratedObservedScore": score,
                })
    expected = pd.DataFrame(rows)
    if expected.empty:
        return expected
    reported_columns = [
        "RunId", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level",
        "ObservedCountDisplayed", "ObservedScoreDisplayed", "BiasSizeDisplayed",
        "PDisplayed", "StrongDecision",
    ]
    reported = (
        table13[[column for column in reported_columns if column in table13]].copy()
        if not table13.empty else
        pd.DataFrame(columns=reported_columns)
    )
    keys = ["RunId", "FacetA", "FacetA_Level", "FacetB", "FacetB_Level"]
    output = expected.merge(reported, on=keys, how="left", validate="one_to_one")
    output["FACETSReported"] = pd.to_numeric(
        output.get("ObservedCountDisplayed", np.nan), errors="coerce"
    ).notna()
    output["FACETSCountNotAboveGenerated"] = np.where(
        output["FACETSReported"],
        pd.to_numeric(output["ObservedCountDisplayed"], errors="coerce")
        .le(pd.to_numeric(output["GeneratedObservationCount"], errors="coerce")),
        np.nan,
    )
    output["FACETSObservedScoreNotAbovePossibleMaximum"] = np.where(
        output["FACETSReported"],
        pd.to_numeric(output["ObservedScoreDisplayed"], errors="coerce")
        .le(
            pd.to_numeric(output["GeneratedObservationCount"], errors="coerce")
            * (pd.to_numeric(output["Categories"], errors="coerce") - 1)
        ),
        np.nan,
    )
    output["ReportingStatus"] = np.where(
        output["FACETSReported"],
        "reported",
        "not_reported_possible_extreme_or_unmeasurable_cell",
    )
    output["NegativeFindingAllowed"] = output["FACETSReported"]
    return output


def facets_python_table8_threshold_agreement(
    categories: pd.DataFrame,
    input_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    path = input_dir / "python_jmle_steps.csv"
    if categories.empty or not path.is_file():
        return pd.DataFrame(), pd.DataFrame()
    python = _read_csv(path)
    required = {"RunId", "Category", "ThresholdEstimate", "IncludedInComparison"}
    missing = required.difference(python.columns)
    if missing:
        raise ValueError(f"Python JMLE step output is missing columns: {sorted(missing)}")
    facets_columns = [
        "RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed",
        "Category", "PrimaryU6ThresholdMeasureToken",
        "PrimaryU6ThresholdMeasureDisplayed",
        "PrimaryU6ThresholdMeasureDisplayDecimals",
        "PrimaryU6ThresholdMeasureRawLowerBound",
        "PrimaryU6ThresholdMeasureRawUpperBound",
        "DirectComparisonEligible", "ComparisonClass",
    ]
    left = categories[[column for column in facets_columns if column in categories]].copy()
    right = python[[
        "RunId", "Category", "ThresholdEstimate", "IncludedInComparison",
    ]].copy()
    pairs = left.merge(right, on=["RunId", "Category"], how="inner", validate="one_to_one")
    pairs["IncludedInDirectComparison"] = (
        pairs["DirectComparisonEligible"].fillna(False).astype(bool)
        & pairs["IncludedInComparison"].fillna(False).astype(bool)
        & pd.to_numeric(pairs["PrimaryU6ThresholdMeasureDisplayed"], errors="coerce").notna()
        & pd.to_numeric(pairs["ThresholdEstimate"], errors="coerce").notna()
    )
    pairs["ThresholdDifference"] = (
        pd.to_numeric(pairs["PrimaryU6ThresholdMeasureDisplayed"], errors="coerce")
        - pd.to_numeric(pairs["ThresholdEstimate"], errors="coerce")
    )
    pairs["AbsoluteThresholdDifference"] = pairs["ThresholdDifference"].abs()
    pairs["PythonWithinFACETSDisplayInterval"] = (
        pd.to_numeric(pairs["ThresholdEstimate"], errors="coerce")
        .ge(pd.to_numeric(pairs["PrimaryU6ThresholdMeasureRawLowerBound"], errors="coerce"))
        & pd.to_numeric(pairs["ThresholdEstimate"], errors="coerce")
        .le(pd.to_numeric(pairs["PrimaryU6ThresholdMeasureRawUpperBound"], errors="coerce"))
    )
    rows: list[dict[str, Any]] = []
    selected = pairs[pairs["IncludedInDirectComparison"]]
    for keys, group in selected.groupby(["Design", "TruthBias"], dropna=False, sort=True):
        rows.append({
            "Design": keys[0],
            "TruthBias": keys[1],
            "Runs": int(group["RunId"].nunique()),
            "Thresholds": int(len(group)),
            "MeanDifference": float(group["ThresholdDifference"].mean()),
            "MAE": float(group["AbsoluteThresholdDifference"].mean()),
            "MaxAbsoluteDifference": float(group["AbsoluteThresholdDifference"].max()),
            "WithinFACETSDisplayInterval": int(group["PythonWithinFACETSDisplayInterval"].sum()),
            "WithinFACETSDisplayIntervalRate": float(group["PythonWithinFACETSDisplayInterval"].mean()),
            "InterpretationBoundary": (
                "matched additive-RSM JMLE threshold comparison; FACETS value remains displayed precision"
            ),
        })
    return pairs, pd.DataFrame(rows)


def write_summary_markdown(
    output_dir: Path,
    runs: pd.DataFrame,
    agreement: pd.DataFrame,
    bias_summary: pd.DataFrame,
    categories: pd.DataFrame,
    table13: pd.DataFrame,
    bias_agreement: pd.DataFrame,
    threshold_agreement: pd.DataFrame,
    table13_cell_audit: pd.DataFrame,
) -> None:
    attempts = len(runs)
    returned = int(runs["FitReturned"].fillna(False).astype(bool).sum())
    converged = int(runs["Converged"].fillna(False).astype(bool).sum())
    eligible = int(runs["AnalysisEligible"].fillna(False).astype(bool).sum())
    direct_runs = int(runs["DirectComparisonEligible"].fillna(False).astype(bool).sum())
    rank_blocked = int((runs["ComparisonClass"] == "structural_nonidentification_negative_control").sum())
    anchor_mismatch = int((runs["ComparisonClass"] == "anchor_constraint_mismatch_sensitivity").sum())
    nonconverged = int((runs["ComparisonClass"] == "facets_nonconvergence").sum())
    measurement_ineligible = int((runs["ComparisonClass"] == "facets_measurement_ineligible").sum())
    fit_failed = int((runs["ComparisonClass"] == "facets_fit_failed").sum())
    anchored_direct = int((runs["DirectComparisonEligible"] & runs["HasAnchors"]).sum())
    replicates_per_condition = int(
        runs.groupby("ConditionId", dropna=False)["Replicate"].nunique().max()
    ) if attempts else 0
    direct = agreement[agreement["Parameters"].gt(0)] if not agreement.empty else agreement
    overall_mae = float(
        np.average(direct["MAE"], weights=direct["Parameters"])
    ) if not direct.empty else np.nan
    maximum = float(direct["MaxAbsoluteDifference"].max()) if not direct.empty else np.nan
    table8_runs = int(runs["Table8Parsed"].fillna(False).astype(bool).sum())
    table13_runs = int(runs["Table13Parsed"].fillna(False).astype(bool).sum())
    category_rows = int(len(categories))
    category_count_matches = int(
        categories.get("TotalCountMatchesGenerated", pd.Series(dtype=bool))
        .fillna(False).astype(bool).sum()
    )
    table13_cells = int(len(table13))
    expected_table13_cells = int(len(table13_cell_audit))
    missing_table13_cells = int(
        (~table13_cell_audit.get("FACETSReported", pd.Series(dtype=bool))
         .fillna(False).astype(bool)).sum()
    )
    matched_bias_cells = int(
        pd.to_numeric(
            bias_agreement.get("MatchedCells", pd.Series(dtype=float)), errors="coerce"
        ).sum()
    )
    eligible_bias_cells = int(
        pd.to_numeric(
            bias_agreement.get("ComparisonEligibleCells", pd.Series(dtype=float)),
            errors="coerce",
        ).sum()
    )
    bias_agreement_cells = int(
        pd.to_numeric(
            bias_agreement.get("StrongDecisionAgreementCells", pd.Series(dtype=float)),
            errors="coerce",
        ).sum()
    )
    threshold_count = int(
        pd.to_numeric(
            threshold_agreement.get("Thresholds", pd.Series(dtype=float)), errors="coerce"
        ).sum()
    )
    threshold_mae = float(
        np.average(
            pd.to_numeric(threshold_agreement["MAE"], errors="coerce"),
            weights=pd.to_numeric(threshold_agreement["Thresholds"], errors="coerce"),
        )
    ) if threshold_count else np.nan
    threshold_maximum = float(
        pd.to_numeric(threshold_agreement["MaxAbsoluteDifference"], errors="coerce").max()
    ) if threshold_count else np.nan
    table13_direct = table13[
        table13.get("DirectComparisonEligible", pd.Series(False, index=table13.index))
        .fillna(False).astype(bool)
    ]
    strong_flags = int((table13_direct.get("StrongDecision", pd.Series(dtype=str)) == "flag").sum())
    sparse_no_claim = int(
        (table13.get("StrongDecision", pd.Series(dtype=str)) == "sparse_no_claim").sum()
    )
    bias_mapping_text = (
        f"The eligible mapping includes {strong_flags} non-sparse strong flags; "
        f"{bias_agreement_cells}/{eligible_bias_cells} eligible FACETS/Python decisions agreed."
        if strong_flags else
        f"All {eligible_bias_cells} eligible mappings are no-flag, so agreement is "
        "a degenerate wiring check rather than detection evidence."
    )
    anchor_contract_text = (
        "FACETS fixes supplied anchors and estimates every remaining level freely on "
        "that origin. The replayed Python results use the corrected absolute-origin "
        "anchor semantics, so identified anchored runs enter direct comparison."
        if anchored_direct else
        "FACETS fixes supplied anchors and estimates every remaining level freely on "
        "that origin. The replayed Python results use the legacy constraint that "
        "centers only remaining free levels, so anchored runs are sensitivity evidence."
    )
    text = f"""# FACETS 4.5 operating-characteristics replay

> This is a replay/adapter qualification and pipeline pilot, not a product-
> validation or estimator-superiority claim. The retained bundle has
> {replicates_per_condition} replicates per condition, which is insufficient for
> confirmatory bias, coverage, false-positive, or power estimates.

## Run accounting

- FACETS returned {returned}/{attempts} runs, met its registered convergence
  criteria in {converged}/{attempts}, and yielded eligible main-facet measures
  in {eligible}/{attempts}.
- {direct_runs}/{attempts} runs entered the direct matched-JMLE summary:
  {rank_blocked} structurally unidentified sparse runs were retained as negative
  controls; {nonconverged} nonconverged, {measurement_ineligible} measurement-
  ineligible, and {fit_failed} failed FACETS runs were excluded; and
  {anchor_mismatch} anchored runs were retained as constraint-mismatch
  sensitivity evidence.
- Across direct-eligible matched additive-RSM parameters, the weighted FACETS versus
  Python-app JMLE MAE was {overall_mae:.6g} logits and the maximum absolute
  difference was {maximum:.6g} logits.

## Table 8 and Table 13 qualification

- The auxiliary two-decimal reporting pass parsed Table 8 in {table8_runs}/{attempts}
  runs. {category_count_matches}/{category_rows} displayed category-total rows
  matched the exact generated-data counts.
- Across {threshold_count} structurally eligible additive-RSM thresholds, the
  FACETS versus Python-app JMLE MAE was {threshold_mae:.6g} logits and the
  maximum absolute difference was {threshold_maximum:.6g} logits. No equivalence
  tolerance was preregistered, so this is descriptive numerical alignment.
- The canonical arranged-by-N Table 13 parsed in {table13_runs}/{attempts} runs,
  yielding {table13_cells}/{expected_table13_cells} expected Rater-by-Task cells.
  The {missing_table13_cells} unreported sparse/extreme cells remain unavailable,
  never negative. {matched_bias_cells} cells were matched to retained Python
  Holm/practical decision rows.
- Table 13 produced {strong_flags} non-sparse strong flags and {sparse_no_claim}
  sparse no-claim cells. {bias_mapping_text}
- Table 8 category Outfit is a one-decimal displayed statistic with +/-0.05
  uncertainty. Table 13 bias and p-value tokens retain their individual display
  precision; Holm uncertainty is propagated conservatively.
- The {replicates_per_condition}-replicate bundle qualifies parsing, decision
  mapping, runtime, and failure accounting only. Its descriptive flag counts
  are not confirmatory false-positive rates or power estimates.

## Interpretation contract

- FACETS JMLE versus Python-app JMLE is the direct external implementation
  comparison after explicit scale alignment.
- TAM and sirt integrate over a Person distribution; exact CMLE conditions on
  Person totals.  Their truth-recovery rows are operating-characteristic
  sensitivity evidence, not expected FACETS equality.
- The Python generator is independent of all fitted engines.  FACETS' internal
  simulation facility is reserved for secondary parametric-bootstrap checks so
  that the primary study does not privilege FACETS' own fitted model.
- Anchored facets remain on their absolute supplied scale; unanchored facets are
  mean-aligned.  Re-centering an anchored facet would erase anchor
  contamination and is prohibited.
- {anchor_contract_text}
- FACETS score/residual text exposes fit components to two decimals.  Fit fields
  are stored as displayed values plus x +/- 0.005 bounds; threshold overlap is
  `boundary_uncertain`, never forced to pass or flag.
- FACETS Table 7 contains measurement, reliability, and agreement output.
  Table 8 is the rating-scale/category report; calling agreement "Table 8-style"
  is not compatible with FACETS 4.5.

## Next evidence increment

Freeze the final factor grid, estimands, failure accounting, and decision
thresholds before a minimum 500-eligible-replicate confirmatory run for
false-positive/power claims. Extend Table 8 qualification to PCM/multiple
scales before category-parity claims, and do not pool this pilot into the
confirmatory analysis.
"""
    (output_dir / "FACETS_REPLAY_RESULTS.md").write_text(text, encoding="utf-8", newline="\n")


def run_adapter(args: argparse.Namespace) -> None:
    input_dir = args.input_dir.resolve()
    comparison_dir_arg = getattr(args, "comparison_dir", None)
    comparison_dir = (
        comparison_dir_arg.resolve()
        if comparison_dir_arg is not None
        else input_dir
    )
    output_dir = args.output_dir.resolve()
    facets_exe = args.facets_exe.resolve()
    if not facets_exe.is_file():
        raise FileNotFoundError(f"FACETS executable is missing: {facets_exe}")
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    output_dir.mkdir(parents=True)
    work_root = output_dir / "facets_runs"
    work_root.mkdir()

    tables = validate_bundle(input_dir)
    manifest = tables["manifest.csv"].copy()
    if args.run_id:
        requested = set(args.run_id)
        manifest = manifest[manifest["RunId"].astype(str).isin(requested)]
        missing = requested.difference(manifest["RunId"].astype(str))
        if missing:
            raise ValueError(f"Requested RunIds are absent from manifest: {sorted(missing)}")
    if args.limit is not None:
        manifest = manifest.head(args.limit)
    if manifest.empty:
        raise ValueError("No manifest rows selected")

    ratings_all = tables["generated_ratings.csv"]
    truth_all = tables["generated_facet_truth.csv"]
    anchors_all = tables["generated_anchors.csv"]
    identifiability_path = args.identifiability_runs
    if identifiability_path is None:
        candidate = (
            input_dir.parent
            / "operating_characteristics_identifiability_20260809"
            / "jmle_identifiability_runs.csv"
        )
        identifiability_path = candidate if candidate.is_file() else None
    identifiability = (
        _read_csv(identifiability_path.resolve())
        if identifiability_path is not None
        else pd.DataFrame()
    )
    if not identifiability.empty:
        required_identifiability = {"RunId", "StructurallyIdentified", "StructuralNullity"}
        missing = required_identifiability.difference(identifiability.columns)
        if missing:
            raise ValueError(f"Identifiability ledger is missing columns: {sorted(missing)}")
        identifiability = identifiability.drop_duplicates("RunId").set_index("RunId")
    run_rows: list[dict[str, Any]] = []
    recovery_parts: list[pd.DataFrame] = []
    category_parts: list[pd.DataFrame] = []
    table13_parts: list[pd.DataFrame] = []
    for _, manifest_row in manifest.iterrows():
        run_id = str(manifest_row["RunId"])
        ratings = ratings_all[ratings_all["RunId"].astype(str).eq(run_id)].copy()
        truth = truth_all[truth_all["RunId"].astype(str).eq(run_id)].copy()
        anchors = anchors_all[anchors_all["RunId"].astype(str).eq(run_id)].copy()
        run, recovery, categories, table13 = fit_one_run(
            manifest_row,
            ratings=ratings,
            truth=truth,
            anchors=anchors,
            facets_exe=facets_exe,
            work_root=work_root,
            timeout_seconds=args.timeout_seconds,
        )
        if run_id in identifiability.index:
            ident_row = identifiability.loc[run_id]
            structurally_identified = str(ident_row["StructurallyIdentified"]).lower() == "true"
            structural_nullity = int(ident_row["StructuralNullity"])
            identifiability_source = str(identifiability_path.resolve())
        else:
            structurally_identified = False
            structural_nullity = np.nan
            identifiability_source = "unavailable_fail_closed"
        has_anchors = not anchors.empty
        anchor_constraint_matched = bool(
            not has_anchors
            or args.python_anchor_constraint == "absolute_origin"
        )
        direct_eligible = bool(
            run["AnalysisEligible"]
            and structurally_identified
            and anchor_constraint_matched
        )
        if not structurally_identified:
            comparison_class = "structural_nonidentification_negative_control"
        elif not bool(run["FitReturned"]):
            comparison_class = "facets_fit_failed"
        elif not bool(run["Converged"]):
            comparison_class = "facets_nonconvergence"
        elif not bool(run["AnalysisEligible"]):
            comparison_class = "facets_measurement_ineligible"
        elif not anchor_constraint_matched:
            comparison_class = "anchor_constraint_mismatch_sensitivity"
        else:
            comparison_class = "matched_jmle"
        run.update({
            "StructurallyIdentified": structurally_identified,
            "StructuralNullity": structural_nullity,
            "IdentifiabilitySource": identifiability_source,
            "HasAnchors": has_anchors,
            "AnchorConstraintMatchedToPythonApp": anchor_constraint_matched,
            "DirectComparisonEligible": direct_eligible,
            "ComparisonClass": comparison_class,
        })
        if not recovery.empty:
            recovery["StructurallyIdentified"] = structurally_identified
            recovery["StructuralNullity"] = structural_nullity
            recovery["HasAnchors"] = has_anchors
            recovery["AnchorConstraintMatchedToPythonApp"] = anchor_constraint_matched
            recovery["DirectComparisonEligible"] = direct_eligible
            recovery["IncludedInDirectComparison"] = (
                recovery["IncludedInSummary"].fillna(False).astype(bool)
                & direct_eligible
            )
            recovery["ComparisonClass"] = comparison_class
        if not categories.empty:
            categories["StructurallyIdentified"] = structurally_identified
            categories["DirectComparisonEligible"] = direct_eligible
            categories["ComparisonClass"] = comparison_class
        if not table13.empty:
            table13["StructurallyIdentified"] = structurally_identified
            table13["FACETSRunFitReturned"] = bool(run["FitReturned"])
            table13["FACETSRunConverged"] = bool(run["Converged"])
            table13["FACETSRunAnalysisEligible"] = bool(run["AnalysisEligible"])
            table13["DirectComparisonEligible"] = direct_eligible
            table13["ComparisonClass"] = comparison_class
        run_rows.append(run)
        if not recovery.empty:
            recovery_parts.append(recovery)
        if not categories.empty:
            category_parts.append(categories)
        if not table13.empty:
            table13_parts.append(table13)
        print(
            f"[{len(run_rows):03d}/{len(manifest):03d}] {run_id}: "
            f"returned={run['FitReturned']} converged={run['Converged']} "
            f"eligible={run['AnalysisEligible']} "
            f"table8={run['Table8Parsed']} table13={run['Table13Parsed']}",
            flush=True,
        )

    runs = pd.DataFrame(run_rows)
    recovery = pd.concat(recovery_parts, ignore_index=True) if recovery_parts else pd.DataFrame()
    categories = pd.concat(category_parts, ignore_index=True) if category_parts else pd.DataFrame()
    table13 = pd.concat(table13_parts, ignore_index=True) if table13_parts else pd.DataFrame()
    runs.to_csv(output_dir / "facets_runs.csv", index=False)
    recovery.to_csv(output_dir / "facets_parameter_recovery.csv", index=False)
    categories.to_csv(output_dir / "facets_table8_categories.csv", index=False)
    category_audit_columns = [
        column for column in (
            "RunId", "ConditionId", "Design", "TruthBias", "Replicate", "Seed",
            "Category", "TotalCount", "UsedCount", "GeneratedRawCount",
            "TotalCountMatchesGenerated", "UsedCountNotAboveTotal",
            "CategoryOutfitMSDisplayed", "CategoryOutfitMSDisplayDecimals",
            "CategoryOutfitDecision0p5To1p5", "ThresholdOrderDecision",
            "AverageMeasureOrderDecision", "LowCount10",
        )
        if column in categories
    ]
    categories[category_audit_columns].to_csv(
        output_dir / "facets_category_count_audit.csv", index=False
    ) if not categories.empty else pd.DataFrame().to_csv(
        output_dir / "facets_category_count_audit.csv", index=False
    )
    threshold_pairs, threshold_agreement = facets_python_table8_threshold_agreement(
        categories, comparison_dir
    )
    threshold_pairs.to_csv(
        output_dir / "facets_python_table8_threshold_pairs.csv", index=False
    )
    threshold_agreement.to_csv(
        output_dir / "facets_python_table8_threshold_agreement.csv", index=False
    )
    table13.to_csv(output_dir / "facets_table13_bias.csv", index=False)
    table13_cell_audit = table13_expected_cell_audit(manifest, ratings_all, table13)
    table13_cell_audit.to_csv(output_dir / "facets_table13_cell_audit.csv", index=False)
    bias_pairs, bias_agreement = facets_python_bias_agreement(table13, comparison_dir)
    bias_pairs.to_csv(output_dir / "facets_python_bias_pairs.csv", index=False)
    bias_agreement.to_csv(output_dir / "facets_python_bias_agreement.csv", index=False)

    run_scope = runs[["RunId", "ComparisonClass", "DirectComparisonEligible"]]
    all_recovery = [recovery]
    all_recovery.extend(_normalize_existing_recovery(comparison_dir, run_scope=run_scope))
    combined = pd.concat([frame for frame in all_recovery if not frame.empty], ignore_index=True, sort=False)
    bias_summary = estimator_bias_summary(combined)
    bias_summary.to_csv(output_dir / "estimator_bias_summary.csv", index=False)

    python_path = comparison_dir / "parameter_recovery.csv"
    python_recovery = _read_csv(python_path) if python_path.is_file() else pd.DataFrame()
    if not recovery.empty and not python_recovery.empty:
        all_pairs, all_agreement = jmle_agreement(
            recovery.assign(IncludedInDirectComparison=recovery["IncludedInSummary"]),
            python_recovery,
        )
        pairs, agreement = jmle_agreement(recovery, python_recovery)
    else:
        all_pairs, all_agreement = pd.DataFrame(), pd.DataFrame()
        pairs, agreement = pd.DataFrame(), pd.DataFrame()
    all_pairs.to_csv(output_dir / "facets_python_jmle_pairs_all_conditions.csv", index=False)
    all_agreement.to_csv(output_dir / "facets_python_jmle_agreement_all_conditions.csv", index=False)
    pairs.to_csv(output_dir / "facets_python_jmle_pairs.csv", index=False)
    agreement.to_csv(output_dir / "facets_python_jmle_agreement.csv", index=False)

    first_report = next(work_root.glob("*/report.out.txt"), None)
    reported_version = parse_iteration_report(first_report)["FacetsVersion"] if first_report else "unknown"
    identity = {
        "schema_version": SCHEMA_VERSION,
        "facets_executable": str(facets_exe),
        "facets_executable_sha256": sha256_file(facets_exe),
        "facets_reported_version": reported_version,
        "adapter_sha256": sha256_file(Path(__file__).resolve()),
        "python_version": sys.version,
        "platform": platform.platform(),
        "selected_runs": int(len(manifest)),
        "identifiability_ledger": (
            str(identifiability_path.resolve())
            if identifiability_path is not None
            else None
        ),
        "input_sha256": {
            filename: sha256_file(input_dir / filename)
            for filename in REQUIRED_BUNDLE_COLUMNS
        },
        "comparison_results_dir": str(comparison_dir),
        "comparison_results_sha256": {
            filename: sha256_file(comparison_dir / filename)
            for filename in (
                "parameter_recovery.csv",
                "bias_decisions.csv",
                "bias_decision_stability.csv",
                "python_jmle_steps.csv",
            )
            if (comparison_dir / filename).is_file()
        },
        "facets_contract": {
            "model": "additive RSM with Rater x Task bias declaration",
            "noncenter": ["Person"],
            "positive": ["Person"],
            "convergence": [0.5, 0.01, 0, 0],
            "xtreme": [0.3, 0.5],
            "measure_decimals": 6,
            "table8_auxiliary_measure_decimals": TABLE8_REPORT_DECIMALS,
            "table8_reason": "avoid fixed-width truncation of category Outfit under Umean=...,6",
            "table8_threshold_source": "six-decimal primary pass joined to native Python JMLE steps",
            "table13_arrangement": "N_canonical",
            "table13_bias_direction": "plus; higher score is positive bias measure",
            "table13_holm_alpha": BIAS_ALPHA,
            "table13_practical_logit": BIAS_PRACTICAL_LOGIT,
            "table13_min_count": BIAS_MIN_COUNT,
        },
        "python_anchor_constraint": args.python_anchor_constraint,
    }
    (output_dir / "facets_runtime_identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_summary_markdown(
        output_dir,
        runs,
        agreement,
        bias_summary,
        categories,
        table13,
        bias_agreement,
        threshold_agreement,
        table13_cell_audit,
    )


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument(
        "--comparison-dir",
        type=Path,
        help=(
            "Optional Python-result directory. FACETS always consumes the immutable "
            "bundle from --input-dir; comparison tables may come from a separately "
            "identified refit of content-identical data."
        ),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--facets-exe", type=Path, default=Path(r"C:\Facets\Facets.exe"))
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    parser.add_argument(
        "--identifiability-runs",
        type=Path,
        help="Optional JMLE structural-identifiability RunId ledger; direct comparison fails closed when unavailable.",
    )
    parser.add_argument(
        "--python-anchor-constraint",
        choices=("legacy_free_levels_centered", "absolute_origin"),
        default="legacy_free_levels_centered",
        help=(
            "Declare the Python result bundle's anchor parameterization. "
            "Anchored runs fail closed from direct comparison unless absolute_origin is explicit."
        ),
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--run-id", action="append", help="Repeat to select exact RunIds")
    return parser.parse_args(argv)


def main() -> None:
    run_adapter(parse_args())


if __name__ == "__main__":
    main()
