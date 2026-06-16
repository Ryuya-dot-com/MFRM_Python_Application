from __future__ import annotations

import io
import zipfile
from pathlib import Path

import pandas as pd

from mfrm_app import cli, distribution, exports, frame_bundle, io_tables, preflight, privacy


def test_distribution_parses_requirement_floors_and_ignores_includes():
    parsed = distribution.parse_requirement_floors(
        """
        # runtime dependencies
        -r requirements-base.txt
        numpy>=1.24
        netcdf4>=1.6 ; python_version >= "3.11"
        ignored-package==1.0
        """
    )

    assert parsed["numpy"]["RequirementFloor"] == "1.24"
    assert parsed["netcdf4"]["RequirementPackage"] == "netcdf4"
    assert "ignored-package" not in parsed


def test_distribution_requirements_contract_matches_doctor_runtime_floors():
    requirements_text = Path("requirements.txt").read_text(encoding="utf-8")

    contract = distribution.requirements_floor_contract(
        requirements_text,
        cli.RUNTIME_PACKAGE_FLOORS,
    )

    assert set(contract["DoctorPackage"]) == set(cli.RUNTIME_PACKAGE_FLOORS)
    assert contract["Status"].eq("Ready").all(), contract.to_dict(orient="records")
    netcdf_row = contract.loc[contract["DoctorPackage"].eq("netCDF4")].iloc[0]
    assert netcdf_row["RequirementPackage"] == "netcdf4"


def test_io_tables_read_text_and_detect_wide_layout():
    parsed = io_tables.read_flexible_table(
        "Person,Rater,Score\nP1,R1,4",
        None,
        header=True,
    )

    assert parsed.to_dict(orient="records") == [
        {"Person": "P1", "Rater": "R1", "Score": "4"}
    ]

    wide = pd.DataFrame({
        "Person": ["P1", "P2"],
        "C1": [1, 2],
        "C2": [2, 3],
        "C3": [3, 4],
    })
    detected = io_tables.detect_wide_format_columns(wide)

    assert detected["looks_wide"] is True
    assert detected["probable_id_cols"] == ["Person"]
    assert detected["probable_score_cols"] == ["C1", "C2", "C3"]


def test_io_tables_wide_to_long_pivot_drops_blank_scores():
    wide = pd.DataFrame({
        "Person": ["P1", "P2"],
        "C1": [1, ""],
        "C2": [2, 3],
    })

    long = io_tables.apply_wide_to_long_pivot(
        wide,
        id_cols=["Person"],
        score_cols=["C1", "C2"],
        new_facet_name="Criterion",
        score_col_name="Score",
    )

    assert long.shape == (3, 3)
    assert long["Criterion"].tolist() == ["C1", "C2", "C2"]
    assert long["Score"].tolist() == [1.0, 2.0, 3.0]


def test_privacy_public_export_removes_person_level_outputs():
    frames = {
        "summary": pd.DataFrame({"Metric": ["N"], "Value": [10]}),
        "scorefile": pd.DataFrame({"Person": ["S01"], "Observed": [4]}),
        "measures": pd.DataFrame({
            "Facet": ["Person", "Rater"],
            "Level": ["S01", "R1"],
            "Estimate": [0.2, -0.1],
        }),
    }

    public_frames = privacy.prepare_download_frames_for_privacy(
        frames,
        public_export_mode=True,
    )

    assert "summary" in public_frames
    assert "scorefile" not in public_frames
    assert public_frames["measures"]["Facet"].tolist() == ["Rater"]
    assert "export_privacy_manifest" in public_frames


def test_frame_bundle_adds_only_nonempty_dataframes():
    frames: dict[str, pd.DataFrame] = {}
    nonempty = pd.DataFrame({"A": [1]})
    empty = pd.DataFrame({"A": []})

    assert frame_bundle.safe_frame_key("bias", "Rater x Task") == "bias_Rater_x_Task"
    assert frame_bundle.safe_frame_key("bias", "評価者/課題") == "bias"
    assert frame_bundle.safe_frame_key("bias", "Very Long Label", max_label_length=4) == "bias_Very"

    assert frame_bundle.is_nonempty_frame(nonempty) is True
    assert frame_bundle.is_nonempty_frame(empty) is False
    assert frame_bundle.is_nonempty_frame({"A": [1]}) is False

    assert frame_bundle.add_frame(frames, "empty", empty) is False
    assert "empty" not in frames
    assert frame_bundle.add_frame(frames, "empty", empty, allow_empty=True) is True
    assert frames["empty"].equals(empty)
    assert frame_bundle.add_frame(frames, "table", nonempty) is True
    assert frames["table"].equals(nonempty)
    assert frame_bundle.add_frame(frames, "table", pd.DataFrame({"A": [2]}), overwrite=False) is False
    assert frames["table"].equals(nonempty)

    added = frame_bundle.add_frames(frames, (
        ("second", pd.DataFrame({"B": [3]})),
        ("skip", pd.DataFrame()),
    ))
    assert added == 1
    assert set(frames) == {"empty", "table", "second"}

    assert frame_bundle.add_iterable_frame(frames, "eigen", [1.0, 0.5], column_name="Eigenvalue") is True
    assert frames["eigen"]["Eigenvalue"].tolist() == [1.0, 0.5]
    assert frame_bundle.add_iterable_frame(frames, "letters", "abc", column_name="Value") is False


def test_exports_build_sanitized_bundle_with_expected_assets():
    frames = {
        "../Summary Table": pd.DataFrame({"Metric": ["N"], "Value": [2]}),
        "Long Sheet Name With Invalid/Characters:*?": pd.DataFrame({"A": [1]}),
    }

    assert exports.to_csv_bytes(frames["../Summary Table"]).startswith(b"Metric,Value")
    assert exports.safe_zip_entry_name("../Summary Table", extension="csv") == "Summary_Table.csv"

    html = exports.to_html_report(frames, title="Demo Report").decode("utf-8")
    assert "<title>Demo Report</title>" in html
    assert "Summary Table" in html

    bundle = exports.build_osf_zip(
        frames,
        title="Demo Report",
        text_assets={"README.txt": "hello"},
    )
    with zipfile.ZipFile(io.BytesIO(bundle)) as zf:
        names = set(zf.namelist())

    assert "Summary_Table.csv" in names
    assert "Demo_Report.xlsx" in names
    assert "Demo_Report.html" in names
    assert "README.txt" in names

    table_zip = exports.build_tables_zip(frames)
    with zipfile.ZipFile(io.BytesIO(table_zip)) as zf:
        assert "Summary_Table.csv" in set(zf.namelist())

    named_zip = exports.build_named_asset_zip({"fig 1": b"abc"}, extension="html")
    with zipfile.ZipFile(io.BytesIO(named_zip)) as zf:
        assert zf.read("fig_1.html") == b"abc"

    mixed_zip = exports.build_mixed_asset_zip({"reports/readme.md": "hello"})
    with zipfile.ZipFile(io.BytesIO(mixed_zip)) as zf:
        assert zf.read("reports_readme.md") == b"hello"

    assert exports.bytes_mapping_fingerprint({"a": b"abc"}) == exports.bytes_mapping_fingerprint({"a": "abc"})


def test_preflight_blocks_oversize_hosted_run():
    df = pd.DataFrame({
        "Person": [f"P{i}" for i in range(preflight.ESTIMATION_HOSTED_MAX_OBS + 1)],
        "Rater": ["R1"] * (preflight.ESTIMATION_HOSTED_MAX_OBS + 1),
        "Task": ["T1"] * (preflight.ESTIMATION_HOSTED_MAX_OBS + 1),
        "Score": [1] * (preflight.ESTIMATION_HOSTED_MAX_OBS + 1),
    })

    result = preflight.build_estimation_resource_preflight(
        data=df,
        person_col="Person",
        score_col="Score",
        facet_cols=["Rater", "Task"],
        model_type="RSM",
        est_method="JMLE",
        maxit=400,
        quad_points=15,
        bias_mode="Skip",
        bias_pairs_available=[],
        selected_bias_pair=None,
        compute_strict_marginal=False,
        strict_marginal_pairwise=False,
        strict_marginal_max_pair_cells=400,
        generate_figure_exports=False,
    )

    assert result["block"] is True
    assert preflight.estimation_resource_preflight_status(result) == "Block"
    assert preflight.should_render_estimation_resource_preflight(result) is True


def test_lightweight_doctor_reports_missing_packages_without_streamlit_import(capsys):
    exit_code = cli.run_lightweight_doctor(json_output=True)
    captured = capsys.readouterr()

    assert "package:streamlit" in captured.out
    assert "locale:en.json" in captured.out
    assert "doctor_mode" in captured.out
    assert exit_code in {0, 1}
