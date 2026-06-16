"""Privacy filtering for exported MFRM tables."""

from __future__ import annotations

from collections import OrderedDict

import pandas as pd


PUBLIC_EXPORT_SENSITIVE_FRAME_NAMES: frozenset[str] = frozenset({
    "response_data_row_audit",
    "response_data_excluded_rows",
    "scorefile",
    "fitted_predictions",
    "new_design_predictions",
    "population_person_data",
    "posterior_scores",
    "plausible_values",
    "residuals",
    "person_measures",
    "person_fit_indices",
    "mfrm_stan_id_index_map",
    "mfrm_uto_bayesian_mfrm_id_index_map",
})


def drop_person_level_rows_for_public_export(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Remove person-level measure rows from aggregate tables in public exports."""
    if not isinstance(df, pd.DataFrame) or df.empty or "Facet" not in df.columns:
        return df, 0
    facet_values = df["Facet"].astype(str).str.strip().str.lower()
    mask_person = facet_values.isin({"person", "persons", "student", "students", "learner", "learners"})
    n_removed = int(mask_person.sum())
    if n_removed == 0:
        return df, 0
    return df.loc[~mask_person].reset_index(drop=True), n_removed


def frame_contains_public_identifier_risk(df: pd.DataFrame) -> bool:
    """Heuristic for frames that should not appear in default public exports."""
    if not isinstance(df, pd.DataFrame) or df.empty:
        return False
    identifier_cols = {"person", "student", "learner", "examinee", "candidate", "subject"}
    if any(str(col).strip().lower() in identifier_cols for col in df.columns):
        return True
    facet_cols = {
        "facet1", "facet2", "sourcefacet", "targetfacet",
        "primaryfacet", "secondaryfacet",
    }
    person_tokens = {"person", "persons", "student", "students", "learner", "learners"}
    for col in df.columns:
        if str(col).strip().replace("_", "").lower() not in facet_cols:
            continue
        values = df[col].astype(str).str.strip().str.lower()
        if values.isin(person_tokens).any():
            return True
    return False


def prepare_download_frames_for_privacy(
    frames: dict[str, pd.DataFrame],
    *,
    public_export_mode: bool,
) -> dict[str, pd.DataFrame]:
    """Return export frames with a manifest and public-mode privacy filtering."""
    prepared: OrderedDict[str, pd.DataFrame] = OrderedDict()
    manifest_rows: list[dict] = []
    for name, df in (frames or {}).items():
        frame_name = str(name)
        if not isinstance(df, pd.DataFrame):
            continue
        name_lower = frame_name.lower()
        if public_export_mode and (
            frame_name in PUBLIC_EXPORT_SENSITIVE_FRAME_NAMES
            or "person" in name_lower
            or name_lower.startswith("facets_person")
            or frame_contains_public_identifier_risk(df)
        ):
            manifest_rows.append({
                "Frame": frame_name,
                "Rows": int(len(df)),
                "Columns": int(df.shape[1]),
                "Status": "excluded_public_mode",
                "Reason": "Excluded from public export: row-level or person-level identifiers may be present.",
            })
            continue

        export_df = df.copy()
        removed_person_rows = 0
        if public_export_mode:
            export_df, removed_person_rows = drop_person_level_rows_for_public_export(export_df)
            if export_df.empty and removed_person_rows:
                manifest_rows.append({
                    "Frame": frame_name,
                    "Rows": int(len(df)),
                    "Columns": int(df.shape[1]),
                    "Status": "excluded_public_mode",
                    "Reason": "Excluded from public export: all rows were person-level rows.",
                })
                continue
        prepared[frame_name] = export_df
        manifest_rows.append({
            "Frame": frame_name,
            "Rows": int(len(export_df)),
            "Columns": int(export_df.shape[1]),
            "Status": "included_public_mode" if public_export_mode else "included_private_mode",
            "Reason": (
                f"Included after removing {removed_person_rows} person-level row(s)."
                if removed_person_rows else
                ("Included in public export." if public_export_mode else "Included in complete/private export.")
            ),
        })

    manifest_rows.append({
        "Frame": "export_privacy_manifest",
        "Rows": len(manifest_rows) + 1,
        "Columns": 5,
        "Status": "included",
        "Reason": (
            "Documents public-export filtering decisions. Public mode is not a formal de-identification guarantee."
            if public_export_mode else
            "Complete/private export selected; row-level tables may contain identifiers."
        ),
    })
    prepared["export_privacy_manifest"] = pd.DataFrame(manifest_rows)
    return dict(prepared)

