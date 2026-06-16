"""Table parsing and wide-to-long helpers for MFRM inputs."""

from __future__ import annotations

import io
from pathlib import Path

import pandas as pd


TABLE_DELIMITER_OPTIONS = {
    "Auto": None,
    "Comma (,)": ",",
    "Tab": "\t",
    "Semicolon (;)": ";",
}
TABLE_FILE_UPLOAD_TYPES = ["csv", "tsv", "txt", "xlsx", "xlsm", "parquet", "json"]
TABLE_FILE_UPLOAD_LABEL = "CSV, TSV, TXT, Excel (.xlsx/.xlsm), Parquet, or JSON"

TABLE_FILE_SOFT_WARNING_MB = 50
TABLE_FILE_HARD_LIMIT_MB = 100
TABLE_TEXT_HARD_LIMIT_MB = 10
TABLE_MAX_ROWS = 250_000
TABLE_MAX_COLUMNS = 250
TABLE_MAX_CELLS = 2_500_000
TABLE_FILE_SOFT_WARNING_BYTES = TABLE_FILE_SOFT_WARNING_MB * 1024 * 1024
TABLE_FILE_HARD_LIMIT_BYTES = TABLE_FILE_HARD_LIMIT_MB * 1024 * 1024
TABLE_TEXT_HARD_LIMIT_BYTES = TABLE_TEXT_HARD_LIMIT_MB * 1024 * 1024


def normalize_csv_newlines(value):
    """Normalize CRLF and old-Mac CR line endings before CSV/TSV parsing."""
    if isinstance(value, bytes):
        return value.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return str(value).replace("\r\n", "\n").replace("\r", "\n")


def infer_table_delimiter(text_value: str | bytes | None, file_name: str | None = None) -> str:
    """Infer a simple CSV/TSV delimiter from text content and file name."""
    if text_value is None:
        return "\t" if str(file_name or "").lower().endswith((".tsv", ".txt")) else ","
    if isinstance(text_value, bytes):
        try:
            text = text_value[:8192].decode("utf-8-sig", errors="ignore")
        except Exception:
            text = ""
    else:
        text = str(text_value)[:8192]
    text = normalize_csv_newlines(text)
    lines = [ln for ln in text.split("\n") if ln.strip()][:10]
    if not lines:
        return "\t" if str(file_name or "").lower().endswith((".tsv", ".txt")) else ","
    counts = {
        "\t": sum(ln.count("\t") for ln in lines),
        ";": sum(ln.count(";") for ln in lines),
        ",": sum(ln.count(",") for ln in lines),
    }
    best, n_best = max(counts.items(), key=lambda kv: kv[1])
    if n_best > 0:
        return best
    return "\t" if str(file_name or "").lower().endswith((".tsv", ".txt")) else ","


def resolve_table_delimiter(delimiter: str | None, text_value=None, file_name: str | None = None) -> str:
    """Resolve UI/API delimiter labels to an actual pandas separator."""
    if delimiter in TABLE_DELIMITER_OPTIONS:
        explicit = TABLE_DELIMITER_OPTIONS[delimiter]
        if explicit is not None:
            return explicit
    if delimiter in {",", "\t", ";"}:
        return delimiter
    return infer_table_delimiter(text_value, file_name=file_name)


def normalize_loaded_table(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame with stable string column names after file parsing."""
    if df is None:
        return pd.DataFrame()
    out = df.copy()
    out.columns = [str(c) for c in out.columns]
    return out


def uploaded_file_size_bytes(file_input) -> int | None:
    """Return declared uploaded-file size without reading file contents."""
    if file_input is None:
        return None
    size = getattr(file_input, "size", None)
    try:
        return int(size) if size is not None else None
    except (TypeError, ValueError):
        return None


def enforce_table_file_size(
    size_bytes: int | None,
    *,
    source_label: str,
    hard_limit_bytes: int = TABLE_FILE_HARD_LIMIT_BYTES,
    hard_limit_mb: int = TABLE_FILE_HARD_LIMIT_MB,
) -> None:
    """Reject files above the hard UI parsing limit."""
    if size_bytes is None:
        return
    if int(size_bytes) >= hard_limit_bytes:
        size_mb = int(round(float(size_bytes) / (1024 * 1024)))
        raise ValueError(
            f"{source_label} is {size_mb} MB; files >= {hard_limit_mb} MB "
            "are not parsed in the Streamlit UI. Sample rows locally or run the app on a larger local machine."
        )


def normalize_and_validate_loaded_table(
    df: pd.DataFrame,
    *,
    source_label: str,
    max_rows: int = TABLE_MAX_ROWS,
    max_columns: int = TABLE_MAX_COLUMNS,
    max_cells: int = TABLE_MAX_CELLS,
) -> pd.DataFrame:
    """Normalize parsed data and reject tables that exceed hosted-app budgets."""
    out = normalize_loaded_table(df)
    n_rows, n_cols = out.shape
    n_cells = int(n_rows) * int(n_cols)
    if n_rows > max_rows:
        raise ValueError(
            f"{source_label} has {n_rows:,} rows; the Streamlit UI limit is {max_rows:,} rows."
        )
    if n_cols > max_columns:
        raise ValueError(
            f"{source_label} has {n_cols:,} columns; the Streamlit UI limit is {max_columns:,} columns."
        )
    if n_cells > max_cells:
        raise ValueError(
            f"{source_label} has {n_cells:,} cells; the Streamlit UI limit is {max_cells:,} cells."
        )
    return out


def read_json_table_bytes(raw: bytes) -> pd.DataFrame:
    """Read ordinary JSON records or JSON-lines as a tabular input file."""
    text = normalize_csv_newlines(raw).decode("utf-8-sig", errors="replace").strip()
    if not text:
        return pd.DataFrame()
    try:
        return pd.read_json(io.StringIO(text))
    except ValueError:
        return pd.read_json(io.StringIO(text), lines=True)


def detect_wide_format_columns(df: pd.DataFrame) -> dict:
    """Heuristically detect whether an uploaded table is in wide format."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return {
            "looks_wide": False,
            "probable_id_cols": [],
            "probable_score_cols": [],
            "reason": "empty table",
        }
    numeric_mask: dict[str, bool] = {}
    for c in df.columns:
        coerced = pd.to_numeric(df[c], errors="coerce")
        non_blank = df[c].astype(str).str.strip().replace({"": pd.NA}).dropna()
        if non_blank.empty:
            numeric_mask[c] = False
            continue
        valid = coerced.dropna()
        numeric_mask[c] = len(valid) >= max(1, int(0.7 * len(non_blank)))
    score_candidates = [c for c, is_num in numeric_mask.items() if is_num]
    id_candidates = [c for c, is_num in numeric_mask.items() if not is_num]

    looks_wide = len(score_candidates) >= 3 and len(id_candidates) >= 1
    if looks_wide:
        reason = (
            f"{len(score_candidates)} numeric columns + {len(id_candidates)} "
            "non-numeric columns suggests one row per id with several score "
            "columns (wide-format Excel layout)."
        )
    elif len(score_candidates) <= 1:
        reason = (
            "<=1 numeric column; data is already in long format "
            "(one row per observation)."
        )
    else:
        reason = (
            f"{len(score_candidates)} numeric columns but no non-numeric "
            "id columns; treating as long format."
        )
    return {
        "looks_wide": looks_wide,
        "probable_id_cols": id_candidates,
        "probable_score_cols": score_candidates,
        "reason": reason,
    }


def apply_wide_to_long_pivot(
    df: pd.DataFrame,
    *,
    id_cols: list[str],
    score_cols: list[str],
    new_facet_name: str = "Item",
    score_col_name: str = "Score",
) -> pd.DataFrame:
    """Melt a wide-format rating table into the canonical long format."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return pd.DataFrame()
    id_cols = list(id_cols or [])
    score_cols = list(score_cols or [])
    if not score_cols:
        raise ValueError("apply_wide_to_long_pivot: score_cols is empty.")
    overlap = set(id_cols) & set(score_cols)
    if overlap:
        raise ValueError(
            f"apply_wide_to_long_pivot: id and score columns overlap: {sorted(overlap)}."
        )
    missing = [c for c in (id_cols + score_cols) if c not in df.columns]
    if missing:
        raise ValueError(
            f"apply_wide_to_long_pivot: columns not in input: {missing}."
        )
    if new_facet_name in id_cols:
        raise ValueError(
            f"apply_wide_to_long_pivot: new_facet_name {new_facet_name!r} collides with an id column."
        )
    if score_col_name in id_cols:
        raise ValueError(
            f"apply_wide_to_long_pivot: score_col_name {score_col_name!r} collides with an id column."
        )

    melted = df.melt(
        id_vars=id_cols,
        value_vars=score_cols,
        var_name=new_facet_name,
        value_name=score_col_name,
    )
    melted[score_col_name] = pd.to_numeric(melted[score_col_name], errors="coerce")
    melted = melted.dropna(subset=[score_col_name])
    return melted.reset_index(drop=True)


def read_flexible_table(
    text_value,
    file_input,
    header=True,
    delimiter: str | None = None,
    *,
    file_hard_limit_bytes: int = TABLE_FILE_HARD_LIMIT_BYTES,
    file_hard_limit_mb: int = TABLE_FILE_HARD_LIMIT_MB,
    text_hard_limit_bytes: int = TABLE_TEXT_HARD_LIMIT_BYTES,
    text_hard_limit_mb: int = TABLE_TEXT_HARD_LIMIT_MB,
    max_rows: int = TABLE_MAX_ROWS,
    max_columns: int = TABLE_MAX_COLUMNS,
    max_cells: int = TABLE_MAX_CELLS,
) -> pd.DataFrame:
    """Read a CSV/TSV/TXT/Excel/Parquet/JSON table from text or an upload."""
    def validate(df: pd.DataFrame, *, source_label: str) -> pd.DataFrame:
        return normalize_and_validate_loaded_table(
            df,
            source_label=source_label,
            max_rows=max_rows,
            max_columns=max_columns,
            max_cells=max_cells,
        )

    if file_input is not None:
        name = file_input.name.lower()
        enforce_table_file_size(
            uploaded_file_size_bytes(file_input),
            source_label=name or "uploaded file",
            hard_limit_bytes=file_hard_limit_bytes,
            hard_limit_mb=file_hard_limit_mb,
        )
        try:
            raw = file_input.getvalue()
        except Exception:
            raw = file_input.read()
            try:
                file_input.seek(0)
            except Exception:
                pass
        raw_bytes = raw.encode("utf-8") if isinstance(raw, str) else bytes(raw)
        enforce_table_file_size(
            len(raw_bytes),
            source_label=name or "uploaded file",
            hard_limit_bytes=file_hard_limit_bytes,
            hard_limit_mb=file_hard_limit_mb,
        )
        suffix = Path(name).suffix.lower()
        if suffix in {".xlsx", ".xlsm"}:
            return validate(pd.read_excel(
                io.BytesIO(raw_bytes),
                sheet_name=0,
                header=0 if header else None,
                dtype=str,
                engine="openpyxl",
            ), source_label=name or "uploaded file")
        if suffix == ".parquet":
            return validate(
                pd.read_parquet(io.BytesIO(raw_bytes)),
                source_label=name or "uploaded file",
            )
        if suffix == ".json":
            return validate(
                read_json_table_bytes(raw_bytes),
                source_label=name or "uploaded file",
            )
        raw_norm = normalize_csv_newlines(raw if isinstance(raw, str) else raw_bytes)
        if isinstance(raw_norm, str):
            sep = resolve_table_delimiter(delimiter, raw_norm, file_name=name)
            return validate(
                pd.read_csv(io.StringIO(raw_norm), sep=sep, header=0 if header else None, dtype=str),
                source_label=name or "uploaded file",
            )
        sep = resolve_table_delimiter(delimiter, raw_norm, file_name=name)
        return validate(
            pd.read_csv(io.BytesIO(raw_norm), sep=sep, header=0 if header else None, dtype=str),
            source_label=name or "uploaded file",
        )

    if text_value is None or not str(text_value).strip():
        return pd.DataFrame()
    text_value = normalize_csv_newlines(str(text_value)).strip()
    text_size = len(text_value.encode("utf-8"))
    if text_size >= text_hard_limit_bytes:
        size_mb = int(round(float(text_size) / (1024 * 1024)))
        raise ValueError(
            f"Pasted text is {size_mb} MB; pasted tables are limited to "
            f"{text_hard_limit_mb} MB in the Streamlit UI."
        )
    sep = resolve_table_delimiter(delimiter, text_value)
    return validate(
        pd.read_csv(io.StringIO(text_value), sep=sep, header=0 if header else None, dtype=str),
        source_label="pasted text",
    )
