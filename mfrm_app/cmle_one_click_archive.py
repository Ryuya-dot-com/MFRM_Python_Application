"""Deterministic private archives for repository-only CMLE one-click results.

The archive is intentionally a controlled-access reproduction artifact.  It
contains response rows and identifiers and is neither anonymized nor approved
for public sharing.  Public mode fails closed pending a separate disclosure
contract.
"""

from __future__ import annotations

from collections.abc import Mapping
import hashlib
import io
import json
import math
import zipfile

import numpy as np
import pandas as pd

from mfrm_app import evidence, exports
from mfrm_app.cmle_one_click import run_cmle_one_click_analysis


CMLE_ONE_CLICK_ARCHIVE_SCHEMA_VERSION = "cmle_one_click_private_archive_v1"
CMLE_ONE_CLICK_RESULT_IDENTITY_VERSION = "cmle_one_click_result_identity_v1"
CMLE_ONE_CLICK_INPUT_SCHEMA_VERSION = "cmle_one_click_input_schema_v1"
PRIVATE_REPRODUCTION_MODE = "private_reproduction"
PUBLIC_MODE_ERROR = (
    "Public CMLE one-click export is withheld until a separate de-identification "
    "and disclosure contract is registered."
)
_FIXED_ZIP_TIME = (1980, 1, 1, 0, 0, 0)
_MAX_ARCHIVE_ENTRIES = 64
_MAX_ARCHIVE_ASSET_BYTES = 64 * 1024 * 1024
_MAX_ARCHIVE_TOTAL_BYTES = 256 * 1024 * 1024
_RESULT_ASSET_NAMES = {
    "one_click_summary.csv",
    "result_cards.csv",
    "availability.csv",
    "calibration_summary.csv",
    "calibration_stages.csv",
    "structural_coefficients.csv",
    "facet_parameters.csv",
    "step_parameters.csv",
    "person_results.csv",
    "person_fit_decision_audit.csv",
}


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _json_bytes(payload: object) -> bytes:
    return (evidence.canonical_json(payload) + "\n").encode("utf-8")


def _canonical_frame_bytes(frame: pd.DataFrame) -> bytes:
    if not isinstance(frame, pd.DataFrame):
        raise ValueError("Archive table assets must be pandas DataFrames.")
    if any(not isinstance(column, str) for column in frame.columns):
        raise ValueError("Archive table column names must be strings.")
    if frame.columns.duplicated().any():
        raise ValueError("Archive tables must not contain duplicate columns.")
    return frame.to_csv(
        index=False,
        lineterminator="\n",
        float_format="%.17g",
        na_rep="",
    ).encode("utf-8")


def _input_schema(frame: pd.DataFrame) -> dict[str, object]:
    columns = []
    for name, dtype in frame.dtypes.items():
        kind = getattr(dtype, "kind", "O")
        if kind in {"i", "u"}:
            family = "integer"
        elif kind == "f":
            family = "float"
        elif kind == "b":
            family = "boolean"
        else:
            family = "string"
        columns.append(
            {"name": str(name), "dtype": str(dtype), "family": family}
        )
    return {
        "schema_version": CMLE_ONE_CLICK_INPUT_SCHEMA_VERSION,
        "columns": columns,
    }


def _cell_payload(value: object) -> object:
    if pd.isna(value):
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("CMLE archive input cannot contain NaN or infinity.")
    if isinstance(value, (str, bool, int, float)) or value is None:
        return value
    return str(value)


def _semantic_input_bytes(frame: pd.DataFrame) -> bytes:
    columns = [str(column) for column in frame.columns]
    rows = [
        [_cell_payload(value) for value in row]
        for row in frame.itertuples(index=False, name=None)
    ]
    row_tokens = sorted(evidence.canonical_json(row) for row in rows)
    return _json_bytes({"columns": columns, "sorted_rows": row_tokens})


def _normalize_anchors(value: object) -> pd.DataFrame | None:
    if value is None:
        return None
    if isinstance(value, pd.DataFrame):
        anchors = value.copy()
    else:
        try:
            anchors = pd.DataFrame(list(value))
        except Exception as exc:
            raise ValueError("hard_anchors must be tabular for archive identity.") from exc
    if anchors.empty:
        return None
    required = ["ParameterType", "Facet", "Level", "Value"]
    if not set(required).issubset(anchors.columns):
        raise ValueError("hard_anchors lacks the archive identity columns.")
    out = anchors[required].copy()
    for column in ("ParameterType", "Facet", "Level"):
        out[column] = out[column].astype(str)
    out["Value"] = pd.to_numeric(out["Value"], errors="raise").astype(float)
    if not np.isfinite(out["Value"]).all():
        raise ValueError("hard_anchors contains a non-finite value.")
    return out.sort_values(required, kind="mergesort").reset_index(drop=True)


def _plain_json(value: object) -> object:
    return json.loads(evidence.canonical_json(value))


def _normalized_config(
    calibration_kwargs: Mapping[str, object],
    *,
    engine_version: str,
) -> tuple[dict[str, object], pd.DataFrame | None]:
    if not isinstance(calibration_kwargs, Mapping):
        raise ValueError("calibration_kwargs must be a mapping.")
    kwargs = dict(calibration_kwargs)
    anchors = _normalize_anchors(kwargs.pop("hard_anchors", None))
    if "positive_facets" in kwargs and kwargs["positive_facets"] is not None:
        kwargs["positive_facets"] = sorted(str(value) for value in kwargs["positive_facets"])
    if "facet_cols" in kwargs:
        kwargs["facet_cols"] = [str(value) for value in kwargs["facet_cols"]]
    runner_kwargs = _plain_json(kwargs)
    resolved_settings = {
        "runner_kwargs": runner_kwargs,
        "hard_anchors_asset": "hard_anchors.csv" if anchors is not None else None,
        "hard_anchors_sha256": (
            _sha256(_canonical_frame_bytes(anchors)) if anchors is not None else None
        ),
    }
    return {
        "schema_version": CMLE_ONE_CLICK_ARCHIVE_SCHEMA_VERSION,
        "analysis_type": "cmle.one_click",
        "engine_version": str(engine_version),
        "resolved_settings": resolved_settings,
    }, anchors


def _result_table_assets(result: Mapping[str, object]) -> dict[str, bytes]:
    summary = result.get("summary")
    cards = result.get("cards")
    availability = result.get("availability")
    calibration = result.get("calibration")
    if not all(isinstance(value, pd.DataFrame) for value in (summary, cards, availability)):
        raise ValueError("One-click result lacks summary, cards, or availability tables.")
    if not isinstance(calibration, Mapping):
        raise ValueError("One-click result lacks the retained calibration workflow.")
    calibration_summary = calibration.get("summary")
    calibration_stages = calibration.get("stages")
    if not isinstance(calibration_summary, pd.DataFrame) or not isinstance(
        calibration_stages, pd.DataFrame
    ):
        raise ValueError("Calibration workflow tables are unavailable.")
    stable_calibration_summary = calibration_summary.drop(
        columns=["ElapsedSeconds"], errors="ignore"
    )
    assets = {
        "one_click_summary.csv": _canonical_frame_bytes(summary),
        "result_cards.csv": _canonical_frame_bytes(cards),
        "availability.csv": _canonical_frame_bytes(availability),
        "calibration_summary.csv": _canonical_frame_bytes(
            stable_calibration_summary
        ),
        "calibration_stages.csv": _canonical_frame_bytes(calibration_stages),
    }
    fit = calibration.get("fit")
    if isinstance(fit, Mapping):
        coefficients = fit.get("coefficients")
        if isinstance(coefficients, pd.DataFrame) and not coefficients.empty:
            assets["structural_coefficients.csv"] = _canonical_frame_bytes(
                coefficients
            )
        facets = fit.get("facets")
        facet_parameters = facets.get("others") if isinstance(facets, Mapping) else None
        if isinstance(facet_parameters, pd.DataFrame) and not facet_parameters.empty:
            assets["facet_parameters.csv"] = _canonical_frame_bytes(facet_parameters)
        steps = fit.get("steps")
        if isinstance(steps, pd.DataFrame) and not steps.empty:
            assets["step_parameters.csv"] = _canonical_frame_bytes(steps)
    person_fit = result.get("person_fit")
    if isinstance(person_fit, Mapping):
        persons = person_fit.get("persons")
        decision_audit = person_fit.get("decision_audit")
        if isinstance(persons, pd.DataFrame) and not persons.empty:
            assets["person_results.csv"] = _canonical_frame_bytes(persons)
        if isinstance(decision_audit, pd.DataFrame) and not decision_audit.empty:
            assets["person_fit_decision_audit.csv"] = _canonical_frame_bytes(
                decision_audit
            )
    return assets


def _asset_rows(assets: Mapping[str, bytes]) -> list[dict[str, object]]:
    return [
        {"name": name, "bytes": len(raw), "sha256": _sha256(raw)}
        for name, raw in sorted(assets.items())
    ]


def _deterministic_zip(assets: Mapping[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(
        buffer,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
        strict_timestamps=True,
    ) as archive:
        for name, raw in sorted(assets.items()):
            if "/" in name or "\\" in name or name in {"", ".", ".."}:
                raise ValueError(f"Unsafe deterministic archive name: {name!r}")
            info = zipfile.ZipInfo(name, date_time=_FIXED_ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, raw, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return buffer.getvalue()


def build_cmle_one_click_archive(
    result: Mapping[str, object],
    *,
    input_data: pd.DataFrame,
    calibration_kwargs: Mapping[str, object],
    engine_version: str = "0.2.15-beta",
    privacy_mode: str = PRIVATE_REPRODUCTION_MODE,
) -> dict[str, object]:
    """Build a deterministic private reproduction ZIP and its identities."""

    if privacy_mode != PRIVATE_REPRODUCTION_MODE:
        raise ValueError(PUBLIC_MODE_ERROR)
    if not isinstance(input_data, pd.DataFrame):
        raise ValueError("input_data must be a pandas DataFrame.")
    input_bytes = _canonical_frame_bytes(input_data)
    input_sha = _sha256(input_bytes)
    semantic_sha = _sha256(_semantic_input_bytes(input_data))
    config, anchors = _normalized_config(
        calibration_kwargs, engine_version=engine_version
    )
    resolved_settings = config["resolved_settings"]
    identity = evidence.build_analysis_identity(
        input_data_fingerprint=f"ordered_sha256:{input_sha}",
        resolved_settings=resolved_settings,
        analysis_type="cmle.one_click",
        engine_version=str(engine_version),
    )
    result_assets = _result_table_assets(result)
    result_content_sha = exports.bytes_mapping_fingerprint(result_assets, length=64)
    terminal_status = str(result["summary"].iloc[0]["TerminalStatus"])
    result_identity = {
        "schema_version": CMLE_ONE_CLICK_RESULT_IDENTITY_VERSION,
        "analysis_id": identity.analysis_id,
        "terminal_status": terminal_status,
        "result_content_sha256": result_content_sha,
        "result_assets": _asset_rows(result_assets),
        "decision_input": "finite_unrounded_mnsq",
        "display_decimals": int(result["summary"].iloc[0]["DisplayDecimals"]),
        "public_surface_enabled": False,
    }
    readme = f"""# Read this first / 最初にお読みください

AnalysisID: `{identity.analysis_id}`
Terminal status: `{terminal_status}`

This deterministic archive contains raw response rows and Person/facet identifiers. It is a **private controlled-access reproduction artifact**. It is not anonymized, encrypted, de-identified, or approved for public sharing.

この決定論的アーカイブには、生の応答行とPerson／facet識別子が含まれます。これは**アクセス管理下のprivate再現用成果物**です。匿名化、暗号化、非識別化、公開共有の承認を意味しません。

All fit classifications use finite unrounded MnSq. Display columns are presentation evidence only. A matching hash establishes artifact identity, not model validity, anchor validity, or user comprehension.
"""
    core_assets: dict[str, bytes] = {
        "README_FIRST.md": readme.encode("utf-8"),
        "analysis_config.json": _json_bytes(config),
        "analysis_identity.json": _json_bytes(identity.to_payload()),
        "input_ratings.csv": input_bytes,
        "input_schema.json": _json_bytes(_input_schema(input_data)),
        **result_assets,
        "result_identity.json": _json_bytes(result_identity),
    }
    if anchors is not None:
        core_assets["hard_anchors.csv"] = _canonical_frame_bytes(anchors)
    archive_content_sha = exports.bytes_mapping_fingerprint(core_assets, length=64)
    manifest = {
        "schema_version": CMLE_ONE_CLICK_ARCHIVE_SCHEMA_VERSION,
        "privacy_mode": PRIVATE_REPRODUCTION_MODE,
        "public_surface_enabled": False,
        "analysis_id": identity.analysis_id,
        "ordered_input_sha256": input_sha,
        "semantic_input_sha256": semantic_sha,
        "config_sha256": identity.config_fingerprint,
        "result_content_sha256": result_content_sha,
        "archive_content_sha256": archive_content_sha,
        "terminal_status": terminal_status,
        "asset_count_excluding_manifest": len(core_assets),
        "assets": _asset_rows(core_assets),
        "privacy_warning": (
            "Contains raw response rows and identifiers; private controlled-access "
            "reproduction use only."
        ),
    }
    all_assets = {**core_assets, "archive_manifest.json": _json_bytes(manifest)}
    zip_bytes = _deterministic_zip(all_assets)
    return {
        "zip_bytes": zip_bytes,
        "zip_sha256": _sha256(zip_bytes),
        "archive_content_sha256": archive_content_sha,
        "result_content_sha256": result_content_sha,
        "ordered_input_sha256": input_sha,
        "semantic_input_sha256": semantic_sha,
        "analysis_identity": identity,
        "result_identity": result_identity,
        "manifest": manifest,
        "assets": all_assets,
    }


def _read_zip_assets(zip_bytes: bytes) -> dict[str, bytes]:
    if not isinstance(zip_bytes, (bytes, bytearray)):
        raise ValueError("Archive payload must be bytes.")
    try:
        with zipfile.ZipFile(io.BytesIO(bytes(zip_bytes))) as archive:
            infos = archive.infolist()
            names = [info.filename for info in infos]
            if len(infos) > _MAX_ARCHIVE_ENTRIES:
                raise ValueError("Archive contains too many entries.")
            if len(names) != len(set(names)):
                raise ValueError("Archive contains duplicate entry names.")
            if any(
                "/" in name
                or "\\" in name
                or name in {"", ".", ".."}
                for name in names
            ):
                raise ValueError("Archive contains an unsafe entry name.")
            total = 0
            assets = {}
            for info in infos:
                if info.file_size > _MAX_ARCHIVE_ASSET_BYTES:
                    raise ValueError("Archive entry exceeds the private replay limit.")
                total += info.file_size
                if total > _MAX_ARCHIVE_TOTAL_BYTES:
                    raise ValueError("Archive exceeds the private replay size limit.")
                assets[info.filename] = archive.read(info)
            return assets
    except zipfile.BadZipFile as exc:
        raise ValueError("Archive is not a readable ZIP file.") from exc


def _load_json_asset(assets: Mapping[str, bytes], name: str) -> dict[str, object]:
    try:
        value = json.loads(assets[name].decode("utf-8"))
    except KeyError as exc:
        raise ValueError(f"Archive is missing {name}.") from exc
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Archive {name} is not valid UTF-8 JSON.") from exc
    if not isinstance(value, dict):
        raise ValueError(f"Archive {name} must contain a JSON object.")
    return value


def _read_input_frame(raw: bytes, schema: Mapping[str, object]) -> pd.DataFrame:
    if schema.get("schema_version") != CMLE_ONE_CLICK_INPUT_SCHEMA_VERSION:
        raise ValueError("Unsupported CMLE archive input schema.")
    columns = schema.get("columns")
    if not isinstance(columns, list) or not columns:
        raise ValueError("Archive input schema has no columns.")
    names = [str(item["name"]) for item in columns]
    frame = pd.read_csv(
        io.BytesIO(raw),
        dtype=str,
        keep_default_na=False,
    )
    if frame.columns.tolist() != names:
        raise ValueError("Archive input columns differ from input_schema.json.")
    for item in columns:
        name = str(item["name"])
        family = str(item["family"])
        if family == "integer":
            numeric = pd.to_numeric(frame[name], errors="raise")
            if not np.equal(numeric, np.floor(numeric)).all():
                raise ValueError(f"Archive integer column {name!r} is not integral.")
            frame[name] = numeric.astype(str(item["dtype"]))
        elif family == "float":
            frame[name] = pd.to_numeric(frame[name], errors="raise").astype(float)
        elif family == "boolean":
            mapping = {"True": True, "False": False, "true": True, "false": False}
            if not frame[name].isin(mapping).all():
                raise ValueError(f"Archive boolean column {name!r} is invalid.")
            frame[name] = frame[name].map(mapping).astype(bool)
        elif family != "string":
            raise ValueError(f"Unsupported archive dtype family {family!r}.")
    return frame


def verify_cmle_one_click_archive(zip_bytes: bytes) -> dict[str, object]:
    """Verify every identity and hash before returning archive content."""

    assets = _read_zip_assets(zip_bytes)
    manifest = _load_json_asset(assets, "archive_manifest.json")
    if manifest.get("schema_version") != CMLE_ONE_CLICK_ARCHIVE_SCHEMA_VERSION:
        raise ValueError("Unsupported CMLE one-click archive schema.")
    expected_manifest_keys = {
        "schema_version",
        "privacy_mode",
        "public_surface_enabled",
        "analysis_id",
        "ordered_input_sha256",
        "semantic_input_sha256",
        "config_sha256",
        "result_content_sha256",
        "archive_content_sha256",
        "terminal_status",
        "asset_count_excluding_manifest",
        "assets",
        "privacy_warning",
    }
    if set(manifest) != expected_manifest_keys:
        raise ValueError("Archive manifest shape is not the registered v1 contract.")
    if manifest.get("privacy_mode") != PRIVATE_REPRODUCTION_MODE:
        raise ValueError(PUBLIC_MODE_ERROR)
    if manifest.get("public_surface_enabled") is not False:
        raise ValueError("Archive manifest cannot enable a public surface.")
    rows = manifest.get("assets")
    if not isinstance(rows, list):
        raise ValueError("Archive manifest lacks the asset ledger.")
    declared = {str(row["name"]): row for row in rows}
    if len(declared) != len(rows):
        raise ValueError("Archive manifest contains duplicate asset declarations.")
    core_assets = {
        name: raw for name, raw in assets.items() if name != "archive_manifest.json"
    }
    if set(core_assets) != set(declared):
        missing = sorted(set(declared) - set(core_assets))
        extra = sorted(set(core_assets) - set(declared))
        raise ValueError(f"Archive asset set mismatch; missing={missing}, extra={extra}")
    if int(manifest.get("asset_count_excluding_manifest", -1)) != len(core_assets):
        raise ValueError("Archive manifest asset count does not match the ZIP.")
    for name, raw in core_assets.items():
        row = declared[name]
        if int(row["bytes"]) != len(raw) or str(row["sha256"]) != _sha256(raw):
            raise ValueError(f"Archive asset hash mismatch: {name}")
    archive_content_sha = exports.bytes_mapping_fingerprint(core_assets, length=64)
    if archive_content_sha != manifest.get("archive_content_sha256"):
        raise ValueError("Archive content identity does not match its manifest.")

    identity_payload = _load_json_asset(assets, "analysis_identity.json")
    identity = evidence.AnalysisIdentity.from_payload(identity_payload)
    if manifest.get("analysis_id") != identity.analysis_id:
        raise ValueError("Archive manifest points to a different AnalysisID.")
    config = _load_json_asset(assets, "analysis_config.json")
    if config.get("schema_version") != CMLE_ONE_CLICK_ARCHIVE_SCHEMA_VERSION:
        raise ValueError("Archive analysis configuration schema is unsupported.")
    resolved = config.get("resolved_settings")
    if not isinstance(resolved, dict):
        raise ValueError("Archive resolved settings are unavailable.")
    if evidence.payload_fingerprint(resolved) != identity.config_fingerprint:
        raise ValueError("Archive configuration fingerprint does not match AnalysisID.")
    if manifest.get("config_sha256") != identity.config_fingerprint:
        raise ValueError("Archive manifest configuration identity does not match.")
    input_sha = _sha256(assets["input_ratings.csv"])
    if (
        input_sha != manifest.get("ordered_input_sha256")
        or identity.input_data_fingerprint != f"ordered_sha256:{input_sha}"
    ):
        raise ValueError("Archive ordered input identity does not match AnalysisID.")
    input_schema = _load_json_asset(assets, "input_schema.json")
    input_frame = _read_input_frame(assets["input_ratings.csv"], input_schema)
    semantic_sha = _sha256(_semantic_input_bytes(input_frame))
    if semantic_sha != manifest.get("semantic_input_sha256"):
        raise ValueError("Archive semantic input identity does not match.")
    anchors_name = resolved.get("hard_anchors_asset")
    anchors_sha = resolved.get("hard_anchors_sha256")
    if anchors_name is None:
        if "hard_anchors.csv" in assets or anchors_sha is not None:
            raise ValueError("Archive hard-anchor declaration is inconsistent.")
    elif (
        anchors_name != "hard_anchors.csv"
        or anchors_name not in assets
        or _sha256(assets[anchors_name]) != anchors_sha
    ):
        raise ValueError("Archive hard-anchor identity does not match configuration.")

    result_identity = _load_json_asset(assets, "result_identity.json")
    if result_identity.get("schema_version") != CMLE_ONE_CLICK_RESULT_IDENTITY_VERSION:
        raise ValueError("Unsupported one-click result identity schema.")
    if result_identity.get("analysis_id") != identity.analysis_id:
        raise ValueError("Result identity points to a different AnalysisID.")
    result_rows = result_identity.get("result_assets")
    if not isinstance(result_rows, list):
        raise ValueError("Result identity lacks its asset ledger.")
    result_names = {str(row["name"]) for row in result_rows}
    if len(result_names) != len(result_rows):
        raise ValueError("Result identity contains duplicate asset declarations.")
    if not {
        "one_click_summary.csv",
        "result_cards.csv",
        "availability.csv",
        "calibration_summary.csv",
        "calibration_stages.csv",
    }.issubset(result_names) or not result_names.issubset(_RESULT_ASSET_NAMES):
        raise ValueError("Result identity asset scope is outside the registered contract.")
    result_assets = {name: assets[name] for name in result_names if name in assets}
    if len(result_assets) != len(result_names):
        raise ValueError("A result-identity asset is missing from the archive.")
    for row in result_rows:
        raw = result_assets[str(row["name"])]
        if int(row["bytes"]) != len(raw) or str(row["sha256"]) != _sha256(raw):
            raise ValueError(f"Result asset identity mismatch: {row['name']}")
    result_content_sha = exports.bytes_mapping_fingerprint(result_assets, length=64)
    if (
        result_content_sha != result_identity.get("result_content_sha256")
        or result_content_sha != manifest.get("result_content_sha256")
    ):
        raise ValueError("Result content identity does not match.")
    if manifest.get("terminal_status") != result_identity.get("terminal_status"):
        raise ValueError("Archive terminal status identities do not match.")
    readme = assets["README_FIRST.md"].decode("utf-8")
    if "private controlled-access reproduction artifact" not in readme:
        raise ValueError("Archive privacy warning is missing from README_FIRST.md.")
    return {
        "verified": True,
        "zip_sha256": _sha256(bytes(zip_bytes)),
        "archive_content_sha256": archive_content_sha,
        "result_content_sha256": result_content_sha,
        "analysis_identity": identity,
        "manifest": manifest,
        "result_identity": result_identity,
        "config": config,
        "input_data": input_frame,
        "assets": assets,
    }


def load_cmle_one_click_archive(zip_bytes: bytes) -> dict[str, object]:
    """Verify and load first-read tables without rerunning the estimator."""

    verified = verify_cmle_one_click_archive(zip_bytes)
    assets = verified["assets"]
    tables = {}
    for name in (
        "one_click_summary.csv",
        "result_cards.csv",
        "availability.csv",
        "calibration_summary.csv",
        "calibration_stages.csv",
        "structural_coefficients.csv",
        "facet_parameters.csv",
        "step_parameters.csv",
        "person_results.csv",
        "person_fit_decision_audit.csv",
        "hard_anchors.csv",
    ):
        if name in assets:
            tables[name] = pd.read_csv(
                io.BytesIO(assets[name]), float_precision="round_trip"
            )
    return {**verified, "tables": tables}


def replay_cmle_one_click_archive(zip_bytes: bytes) -> dict[str, object]:
    """Replay a verified private archive and compare every retained identity."""

    original = verify_cmle_one_click_archive(zip_bytes)
    config = original["config"]
    resolved = config["resolved_settings"]
    kwargs = dict(resolved["runner_kwargs"])
    if resolved.get("hard_anchors_asset"):
        kwargs["hard_anchors"] = pd.read_csv(
            io.BytesIO(original["assets"]["hard_anchors.csv"])
        )
    result = run_cmle_one_click_analysis(original["input_data"].copy(), **kwargs)
    rebuilt = build_cmle_one_click_archive(
        result,
        input_data=original["input_data"].copy(),
        calibration_kwargs=kwargs,
        engine_version=str(config["engine_version"]),
        privacy_mode=PRIVATE_REPRODUCTION_MODE,
    )
    original_result_assets = {
        row["name"]: row["sha256"]
        for row in original["result_identity"]["result_assets"]
    }
    rebuilt_result_assets = {
        row["name"]: row["sha256"]
        for row in rebuilt["result_identity"]["result_assets"]
    }
    checks = {
        "AnalysisIDMatch": (
            rebuilt["analysis_identity"].analysis_id
            == original["analysis_identity"].analysis_id
        ),
        "ResultContentSHA256Match": (
            rebuilt["result_content_sha256"] == original["result_content_sha256"]
        ),
        "ResultAssetSHA256Match": rebuilt_result_assets == original_result_assets,
        "TerminalStatusMatch": (
            rebuilt["result_identity"]["terminal_status"]
            == original["result_identity"]["terminal_status"]
        ),
        "CardBytesMatch": (
            rebuilt["assets"]["result_cards.csv"]
            == original["assets"]["result_cards.csv"]
        ),
        "ArchiveContentSHA256Match": (
            rebuilt["archive_content_sha256"]
            == original["archive_content_sha256"]
        ),
        "ZipSHA256Match": rebuilt["zip_sha256"] == original["zip_sha256"],
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "original": original,
        "rebuilt": rebuilt,
        "result": result,
    }


__all__ = [
    "CMLE_ONE_CLICK_ARCHIVE_SCHEMA_VERSION",
    "PRIVATE_REPRODUCTION_MODE",
    "PUBLIC_MODE_ERROR",
    "build_cmle_one_click_archive",
    "load_cmle_one_click_archive",
    "replay_cmle_one_click_archive",
    "verify_cmle_one_click_archive",
]
