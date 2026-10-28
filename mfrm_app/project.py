"""Versioned basic-analysis files. Only explicitly supported data types are read.

This is a result archive, not an optimizer checkpoint or an executable object
format. JSON size and array shapes are checked before constructing arrays.
"""
from __future__ import annotations

import base64
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeResult

from mfrm_app import engine

FORMAT = "mfrm.basic-project"
SCHEMA_VERSION = 1
MAX_PROJECT_BYTES = 128 * 1024 * 1024


class ProjectError(ValueError):
    """The archive is unsupported, incomplete, changed, or malformed."""


def _dtype_encode(dtype):
    if isinstance(dtype, pd.CategoricalDtype):
        return ["category", _encode(dtype.categories), dtype.ordered]
    if isinstance(dtype, pd.StringDtype):
        return ["string", dtype.storage, _encode(dtype.na_value)]
    if isinstance(dtype, np.dtype):
        if dtype.fields or dtype.subdtype:
            raise ProjectError("Structured array dtypes are not supported in schema 1")
        return ["numpy", dtype.str]
    if str(dtype) in {"boolean", "Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64", "Float32", "Float64"}:
        return ["nullable", str(dtype)]
    raise ProjectError(f"Unsupported dtype: {dtype}")


def _dtype_decode(spec):
    kind = spec[0]
    if kind == "category":
        return pd.CategoricalDtype(_decode(spec[1]), ordered=spec[2])
    if kind == "string" and spec[1] in {"python", "pyarrow"}:
        return pd.StringDtype(storage=spec[1], na_value=_decode(spec[2]))
    if kind == "nullable" and spec[1] in {"boolean", "Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64", "Float32", "Float64"}:
        return pd.api.types.pandas_dtype(spec[1])
    if kind == "numpy" and isinstance(spec[1], str):
        dtype = np.dtype(spec[1])
        if not dtype.fields and not dtype.subdtype and dtype.kind in "biufUSOmM" and 0 < dtype.itemsize <= MAX_PROJECT_BYTES:
            return dtype
    raise ProjectError("Unsupported stored dtype")


def _encode(value):
    if value is pd.NA:
        return ["NA"]
    if value is pd.NaT:
        return ["NaT"]
    if value is None or type(value) in (str, int, bool):
        return value
    if type(value) is float:
        return ["float", value.hex()]
    if isinstance(value, np.ndarray):
        data = (_encode(value.ravel().tolist()) if value.dtype.hasobject
                else base64.b64encode(value.tobytes(order="C")).decode("ascii"))
        return ["array", _dtype_encode(value.dtype), list(value.shape), data]
    if isinstance(value, np.generic):
        return ["scalar", _encode(np.asarray(value))]
    if isinstance(value, pd.RangeIndex):
        return ["range", value.start, value.stop, value.step, _encode(value.name)]
    if isinstance(value, pd.MultiIndex):
        return ["multi-index", _encode(list(value.levels)), _encode(list(value.codes)), _encode(list(value.names))]
    if isinstance(value, pd.Index):
        return ["index", _dtype_encode(value.dtype), _encode(value.tolist()), _encode(value.name)]
    if isinstance(value, pd.Series):
        return ["series", _dtype_encode(value.dtype), _encode(value.tolist()), _encode(value.index), _encode(value.name), _encode(value.attrs)]
    if isinstance(value, pd.DataFrame):
        return ["frame", _encode(value.index), _encode(value.columns),
                [_encode(value.iloc[:, i]) for i in range(value.shape[1])], _encode(value.attrs)]
    if isinstance(value, OptimizeResult):
        # The inverse-Hessian operator is optimizer machinery, not a reported
        # uncertainty estimate. It is deliberately absent from result archives.
        return ["optimizer", _encode({k: v for k, v in value.items() if k != "hess_inv"})]
    if type(value) is dict:
        return ["dict", [[_encode(k), _encode(v)] for k, v in value.items()]]
    if type(value) in (list, tuple):
        return [type(value).__name__, [_encode(v) for v in value]]
    if isinstance(value, pd.Timestamp):
        return ["timestamp", value.isoformat()]
    raise ProjectError(f"Unsupported project value: {type(value).__name__}")


def _decode(value, depth=0):
    if depth > 80:
        raise ProjectError("Project nesting exceeds schema 1 limit")
    if value is None or type(value) in (str, int, bool):
        return value
    if not isinstance(value, list) or not value or not isinstance(value[0], str):
        raise ProjectError("Expected a typed project value")
    tag, *parts = value
    decode = lambda item: _decode(item, depth + 1)
    if tag == "NA" and not parts:
        return pd.NA
    if tag == "NaT" and not parts:
        return pd.NaT
    if tag == "float" and len(parts) == 1:
        return float.fromhex(parts[0])
    if tag == "array" and len(parts) == 3:
        dtype, shape, data = _dtype_decode(parts[0]), parts[1], parts[2]
        if (not isinstance(shape, list) or len(shape) > 32
                or any(type(n) is not int or n < 0 or n > MAX_PROJECT_BYTES for n in shape)):
            raise ProjectError("Invalid array shape")
        size = math.prod(shape)
        if size * dtype.itemsize > MAX_PROJECT_BYTES:
            raise ProjectError("Stored array exceeds schema 1 limit")
        if dtype.hasobject:
            items = decode(data)
            if not isinstance(items, list) or len(items) != size:
                raise ProjectError("Array content/shape mismatch")
            array = np.empty(size, dtype=object)
            array[:] = items
        else:
            raw = base64.b64decode(data, validate=True)
            if len(raw) != size * dtype.itemsize:
                raise ProjectError("Array content/shape mismatch")
            array = np.frombuffer(raw, dtype=dtype).copy()
        return array.reshape(shape)
    if tag == "scalar" and len(parts) == 1:
        array = decode(parts[0])
        if not isinstance(array, np.ndarray) or array.shape != ():
            raise ProjectError("Invalid scalar")
        return array[()]
    if tag == "range" and len(parts) == 4:
        index = pd.RangeIndex(*parts[:3], name=decode(parts[3]))
        if len(index) > MAX_PROJECT_BYTES:
            raise ProjectError("Index exceeds schema 1 limit")
        return index
    if tag == "multi-index" and len(parts) == 3:
        return pd.MultiIndex(levels=decode(parts[0]), codes=decode(parts[1]), names=decode(parts[2]), verify_integrity=True)
    if tag == "index" and len(parts) == 3:
        return pd.Index(decode(parts[1]), dtype=_dtype_decode(parts[0]), name=decode(parts[2]), tupleize_cols=False)
    if tag == "series" and len(parts) == 5:
        items, dtype, index = decode(parts[1]), _dtype_decode(parts[0]), decode(parts[2])
        if not isinstance(items, list) or len(items) != len(index):
            raise ProjectError("Series axes/content mismatch")
        if len(items) * getattr(dtype, "itemsize", 8) > MAX_PROJECT_BYTES:
            raise ProjectError("Series exceeds schema 1 limit")
        series = pd.Series(items, dtype=dtype, index=index, name=decode(parts[3]))
        series.attrs = decode(parts[4])
        return series
    if tag == "frame" and len(parts) == 4:
        index, columns = decode(parts[0]), decode(parts[1])
        series = [decode(part) for part in parts[2]]
        if len(series) != len(columns) or any(not s.index.equals(index) for s in series):
            raise ProjectError("Table axes/content mismatch")
        frame = pd.concat(series, axis=1) if series else pd.DataFrame(index=index)
        frame.columns = columns
        frame.attrs = decode(parts[3])
        return frame
    if tag in {"list", "tuple"} and len(parts) == 1:
        items = [decode(item) for item in parts[0]]
        return tuple(items) if tag == "tuple" else items
    if tag == "dict" and len(parts) == 1:
        result = {}
        for key, item in parts[0]:
            key = decode(key)
            if key in result:
                raise ProjectError("Duplicate mapping key")
            result[key] = decode(item)
        return result
    if tag == "optimizer" and len(parts) == 1:
        return OptimizeResult(decode(parts[0]))
    if tag == "timestamp" and len(parts) == 1:
        return pd.Timestamp(parts[0])
    raise ProjectError(f"Unsupported project tag: {tag}")


def _json_bytes(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode("utf-8")


def content_fingerprint(analysis):
    content = {k: v for k, v in analysis.items() if k != "content_fingerprint"}
    return hashlib.sha256(_json_bytes(_encode(content))).hexdigest()


def validate_analysis(analysis):
    if not isinstance(analysis, dict) or analysis.get("kind") != "basic-rsm-jmle-v1":
        raise ProjectError("Only basic RSM/JMLE projects are supported in schema 1")
    if analysis.get("content_fingerprint") != content_fingerprint(analysis):
        raise ProjectError("Analysis content changed: input, settings, results or tables no longer match")
    result = analysis["result"]
    if result["config"].get("model") != "RSM" or result["config"].get("method") != "JMLE":
        raise ProjectError("Project model/method does not match the supported workflow")
    if engine.build_result_analysis_identity(result).analysis_id != analysis.get("analysis_id"):
        raise ProjectError("Analysis identity does not match the stored result")
    if not isinstance(result["prep"].get("data"), pd.DataFrame):
        raise ProjectError("Schema 1 requires the retained input data")
    data, settings = analysis.get("input_data"), analysis.get("fit_settings")
    required = {"person_col", "facet_cols", "score_col", "rating_min", "rating_max"}
    if not isinstance(data, pd.DataFrame) or not isinstance(settings, dict) or not required.issubset(settings):
        raise ProjectError("Schema 1 requires input data and complete column/score settings")
    if not isinstance(settings["facet_cols"], list) or not settings["facet_cols"]:
        raise ProjectError("Project facet settings are missing")
    columns = [settings["person_col"], *settings["facet_cols"], settings["score_col"]]
    if len(set(columns)) != len(columns) or any(column not in data for column in columns):
        raise ProjectError("Project column settings do not match its input data")


def save_project(path, analysis):
    """Atomically save a validated snapshot; never damage an existing file on failure."""
    validate_analysis(analysis)
    document = {"format": FORMAT, "schema_version": SCHEMA_VERSION, "analysis": _encode(analysis)}
    data = _json_bytes(document)
    if len(data) > MAX_PROJECT_BYTES:
        raise ProjectError("Project exceeds the schema 1 size limit (128 MiB)")
    write_bytes_atomic(path, data)


def write_bytes_atomic(path, data):
    """Write a completed project/report, retaining any old file if replacement fails."""
    path = Path(path)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def load_project(path):
    """Restore a completed or failed result without running the estimator."""
    try:
        with Path(path).open("rb") as stream:
            data = stream.read(MAX_PROJECT_BYTES + 1)
        if len(data) > MAX_PROJECT_BYTES:
            raise ProjectError("Project exceeds the schema 1 size limit (128 MiB)")
        document = json.loads(data)
        if not isinstance(document, dict) or document.get("format") != FORMAT:
            raise ProjectError("Not an MFRM basic project")
        if type(document.get("schema_version")) is not int or document["schema_version"] != SCHEMA_VERSION:
            raise ProjectError("Unsupported project schema version")
        if set(document) != {"format", "schema_version", "analysis"}:
            raise ProjectError("Unknown project document fields")
        analysis = _decode(document["analysis"])
        validate_analysis(analysis)
        return analysis
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, OverflowError, RecursionError) as exc:
        if isinstance(exc, ProjectError):
            raise
        raise ProjectError(f"Malformed project: {exc}") from exc
