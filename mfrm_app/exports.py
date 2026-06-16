from __future__ import annotations

import hashlib
import io
import re
import struct
import zipfile

import pandas as pd


def to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


def safe_export_name(name: object, *, default: str = "table", max_len: int = 96) -> str:
    """Sanitize a filename or sheet-name stem for generated export bundles."""
    text = str(name or "").strip()
    text = text.replace("\\", "_").replace("/", "_")
    text = re.sub(r"[^A-Za-z0-9._ -]+", "_", text)
    text = re.sub(r"\s+", "_", text).strip("._- ")
    if not text:
        text = default
    if text.startswith("."):
        text = text.lstrip(".") or default
    return text[:max_len].rstrip("._- ") or default


def safe_zip_entry_name(name: object, *, extension: str | None = None) -> str:
    stem = safe_export_name(name, default="asset")
    ext = ""
    if extension:
        ext = str(extension).strip()
        ext = ext if ext.startswith(".") else f".{ext}"
        ext = re.sub(r"[^A-Za-z0-9.]+", "", ext)
    if ext and not stem.lower().endswith(ext.lower()):
        return f"{stem}{ext}"
    return stem


def unique_excel_sheet_name(name: object, used: set[str]) -> str:
    base = safe_export_name(name, default="Sheet", max_len=31)
    base = re.sub(r"[\[\]:*?/\\]", "_", base).strip("'") or "Sheet"
    candidate = base[:31]
    i = 2
    while candidate in used:
        suffix = f"_{i}"
        candidate = f"{base[:31 - len(suffix)]}{suffix}"
        i += 1
    used.add(candidate)
    return candidate


def to_excel_bytes(frames: dict[str, pd.DataFrame]) -> bytes:
    """Write multiple DataFrames to an in-memory Excel workbook, one sheet per key."""
    buf = io.BytesIO()
    used_sheets: set[str] = set()
    with pd.ExcelWriter(buf, engine="openpyxl") as writer:
        for name, df in frames.items():
            sheet = unique_excel_sheet_name(name, used_sheets)
            df.to_excel(writer, sheet_name=sheet, index=False)
    return buf.getvalue()


def to_html_report(frames: dict[str, pd.DataFrame], title: str = "MFRM Report") -> bytes:
    """Build a self-contained HTML report from multiple DataFrames."""
    parts = [
        "<!DOCTYPE html><html><head>",
        f"<title>{title}</title>",
        "<meta charset='utf-8'>",
        "<style>",
        "body{font-family:system-ui,sans-serif;margin:2em;color:#222}",
        "h1{border-bottom:2px solid #333}",
        "h2{margin-top:1.5em;color:#444}",
        "table{border-collapse:collapse;margin:1em 0;width:100%}",
        "th,td{border:1px solid #ccc;padding:4px 8px;text-align:left;font-size:0.85em}",
        "th{background:#f5f5f5}",
        "tr:nth-child(even){background:#fafafa}",
        "</style></head><body>",
        f"<h1>{title}</h1>",
    ]
    for name, df in frames.items():
        parts.append(f"<h2>{name}</h2>")
        parts.append(df.to_html(index=False, border=0, na_rep=""))
    parts.append("</body></html>")
    return "\n".join(parts).encode("utf-8")


def build_tables_zip(frames: dict[str, pd.DataFrame]) -> bytes:
    """Create a ZIP archive containing one CSV per DataFrame."""
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, df in frames.items():
            zf.writestr(safe_zip_entry_name(name, extension="csv"), df.to_csv(index=False))
    return zip_buf.getvalue()


def build_osf_zip(
    frames: dict[str, pd.DataFrame],
    title: str = "MFRM_Report",
    text_assets: dict[str, str] | None = None,
) -> bytes:
    """Create a ZIP archive with CSV, Excel, HTML, and optional text assets."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        written: set[str] = set()
        for name, df in frames.items():
            entry = safe_zip_entry_name(name, extension="csv")
            if entry not in written:
                zf.writestr(entry, df.to_csv(index=False))
                written.add(entry)
        xlsx_entry = safe_zip_entry_name(title, extension="xlsx")
        if xlsx_entry not in written:
            zf.writestr(xlsx_entry, to_excel_bytes(frames))
            written.add(xlsx_entry)
        html_entry = safe_zip_entry_name(title, extension="html")
        if html_entry not in written:
            zf.writestr(html_entry, to_html_report(frames, title))
            written.add(html_entry)
        for name, text in (text_assets or {}).items():
            entry = safe_zip_entry_name(name)
            if entry not in written:
                zf.writestr(entry, str(text))
                written.add(entry)
    return buf.getvalue()


def bytes_mapping_fingerprint(assets: dict[str, bytes | str], length: int = 16) -> str:
    """Fingerprint a named bytes/string asset bundle for ZIP export caching."""
    digest = hashlib.sha256()
    for name, value in sorted((assets or {}).items(), key=lambda item: str(item[0])):
        raw = value.encode("utf-8") if isinstance(value, str) else bytes(value)
        digest.update(str(name).encode("utf-8"))
        digest.update(struct.pack("<q", len(raw)))
        digest.update(hashlib.sha256(raw).digest())
    return digest.hexdigest()[: int(length)]


def build_named_asset_zip(assets: dict[str, bytes | str], extension: str) -> bytes:
    """Create a ZIP archive for generated assets sharing one extension."""
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, value in assets.items():
            raw = value.encode("utf-8") if isinstance(value, str) else bytes(value)
            zf.writestr(safe_zip_entry_name(name, extension=extension), raw)
    return zip_buf.getvalue()


def build_mixed_asset_zip(assets: dict[str, bytes | str]) -> bytes:
    """Create a ZIP archive for assets whose names already include extensions."""
    zip_buf = io.BytesIO()
    with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, value in assets.items():
            raw = value.encode("utf-8") if isinstance(value, str) else bytes(value)
            zf.writestr(safe_zip_entry_name(name), raw)
    return zip_buf.getvalue()
