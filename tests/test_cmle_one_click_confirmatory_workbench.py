from __future__ import annotations

import io
import zipfile

import pandas as pd

from mfrm_app.cmle_one_click_confirmatory_workbench import (
    DECISION_IDS,
    PREREQUISITES,
    SELECTION_COLUMNS,
    build_decision_dependency_table,
    build_decision_option_catalog,
    build_decision_selection_template,
    build_private_decision_workbench_bundle,
    compute_workbench_content_identity,
    decision_workbench_human_gate_status,
    render_decision_workbench_html,
    validate_decision_selections,
    validate_dependency_table,
    validate_private_decision_workbench_bundle,
    validate_workbench_html,
)


def complete_synthetic_selections() -> pd.DataFrame:
    catalog = build_decision_option_catalog()
    dependencies = build_decision_dependency_table().sort_values("TopologicalOrder")
    selections = build_decision_selection_template()
    for decision_id in dependencies["DecisionId"]:
        option = catalog.loc[catalog["DecisionId"].eq(decision_id)].iloc[0]
        mask = selections["DecisionId"].eq(decision_id)
        selections.loc[mask, "DecisionStatus"] = "resolved_prospectively"
        selections.loc[mask, "SelectedOptionId"] = option["OptionId"]
        selections.loc[mask, "AssumptionSetReference"] = "SYNTHETIC-ASSUMPTIONS"
        selections.loc[mask, "RationaleReference"] = "SYNTHETIC-RATIONALE"
        selections.loc[mask, "OwnerAttestationReference"] = "SYNTHETIC-OWNER"
        selections.loc[mask, "DecisionRecordedBeforeOutcomes"] = True
        if bool(option["RequiresNumericInput"]):
            selections.loc[mask, "NumericValue"] = "100"
            selections.loc[mask, "NumericUnit"] = "synthetic_units"
        if bool(option["ExternalAuthorityRequired"]):
            selections.loc[mask, "ExternalEvidenceReference"] = "SYNTHETIC-EXTERNAL"
            selections.loc[mask, "EvidenceSource"] = "synthetic_external_fixture"
    return selections


def test_catalog_has_three_nonrecommended_options_per_decision() -> None:
    catalog = build_decision_option_catalog()
    assert len(catalog) == 39
    assert set(catalog["DecisionId"]) == set(DECISION_IDS)
    assert catalog.groupby("DecisionId").size().eq(3).all()
    assert catalog["OptionId"].is_unique
    assert not catalog["OptionOrderIsRank"].any()
    assert not catalog["AutoSelectable"].any()
    assert catalog["RecommendationStatus"].eq(
        "not_recommended_by_repository"
    ).all()
    assert catalog[
        [
            "Strength",
            "PrimaryRisk",
            "RequiredEvidence",
            "SampleSizeImpact",
            "EstimandImpact",
            "ParticipantBurdenImpact",
            "AccessibilityImpact",
            "PublicClaimImpact",
        ]
    ].astype(str).apply(lambda column: column.str.strip().ne("").all()).all()


def test_dependency_graph_is_exact_acyclic_and_topological() -> None:
    dependencies = build_decision_dependency_table()
    result = validate_dependency_table(dependencies)
    assert result["valid"]
    assert len(dependencies) == 13
    assert dependencies["TopologicalOrder"].is_unique
    sample = dependencies.loc[dependencies["DecisionId"].eq("SCI-02")].iloc[0]
    assert sample["PrerequisiteCount"] == 7
    assert tuple(sample["PrerequisiteDecisionIds"].split("|")) == PREREQUISITES[
        "SCI-02"
    ]
    assert not dependencies["StageIsPermissionToSkipDependencies"].any()


def test_default_selection_template_is_valid_blank_and_incomplete() -> None:
    selections = build_decision_selection_template()
    result = validate_decision_selections(selections)
    assert tuple(selections.columns) == SELECTION_COLUMNS
    assert tuple(selections["DecisionId"]) == DECISION_IDS
    assert result["valid"]
    assert result["OptionsSelected"] == 0
    assert not result["SelectionsComplete"]
    assert not result["SubstantiveEvidenceVerified"]
    assert not result["RecruitmentReady"]


def test_unresolved_row_cannot_hide_a_selection() -> None:
    selections = build_decision_selection_template()
    selections.loc[0, "SelectedOptionId"] = "SCI-01-A_observed_valid_estimand"
    result = validate_decision_selections(selections)
    assert not result["valid"]
    assert "unresolved_row_contains_selection" in result["failure_codes"]


def test_unknown_and_cross_decision_options_fail() -> None:
    selections = build_decision_selection_template()
    mask = selections["DecisionId"].eq("SCI-01")
    selections.loc[mask, "DecisionStatus"] = "resolved_prospectively"
    selections.loc[mask, "SelectedOptionId"] = "SCI-04-A_primary_all_cells_sensitivity_only"
    selections.loc[mask, ["AssumptionSetReference", "RationaleReference", "OwnerAttestationReference"]] = "SYNTHETIC"
    selections.loc[mask, "DecisionRecordedBeforeOutcomes"] = True
    result = validate_decision_selections(selections)
    assert not result["valid"]
    assert "cross_decision_option" in result["failure_codes"]


def test_child_decision_cannot_resolve_before_prerequisite() -> None:
    selections = build_decision_selection_template()
    mask = selections["DecisionId"].eq("SCI-04")
    selections.loc[mask, "DecisionStatus"] = "resolved_prospectively"
    selections.loc[mask, "SelectedOptionId"] = "SCI-04-A_primary_all_cells_sensitivity_only"
    selections.loc[mask, ["AssumptionSetReference", "RationaleReference", "OwnerAttestationReference"]] = "SYNTHETIC"
    selections.loc[mask, "DecisionRecordedBeforeOutcomes"] = True
    result = validate_decision_selections(selections)
    assert not result["valid"]
    assert "unresolved_prerequisite" in result["failure_codes"]


def test_numeric_option_requires_positive_value_unit_and_assumptions() -> None:
    selections = complete_synthetic_selections()
    mask = selections["DecisionId"].eq("SCI-02")
    selections.loc[mask, "NumericValue"] = ""
    selections.loc[mask, "NumericUnit"] = ""
    result = validate_decision_selections(selections)
    assert not result["valid"]
    assert "missing_or_invalid_numeric_input" in result["failure_codes"]


def test_external_option_rejects_repository_generated_evidence() -> None:
    selections = complete_synthetic_selections()
    mask = selections["DecisionId"].eq("ETH-01")
    selections.loc[mask, "EvidenceSource"] = "repository_generated_fake"
    result = validate_decision_selections(selections)
    assert not result["valid"]
    assert "invalid_external_evidence_source" in result["failure_codes"]


def test_synthetic_complete_worksheet_never_grants_recruitment_readiness() -> None:
    selections = complete_synthetic_selections()
    result = validate_decision_selections(selections)
    assert result["valid"]
    assert result["OptionsSelected"] == 13
    assert result["SelectionsComplete"]
    assert not result["SubstantiveEvidenceVerified"]
    assert not result["RecruitmentReady"]


def test_html_is_self_contained_bilingual_complete_and_blocked() -> None:
    catalog = build_decision_option_catalog()
    dependencies = build_decision_dependency_table()
    selections = build_decision_selection_template()
    identity = compute_workbench_content_identity(catalog, dependencies, selections)
    html = render_decision_workbench_html(catalog, dependencies, identity)
    validation = validate_workbench_html(html, catalog, identity)
    assert validation["valid"]
    assert "BLOCKED" in html
    assert "募集不可" in html
    assert "<form" not in html.lower()
    assert "<script" not in html.lower()
    assert "http://" not in html.lower()
    assert "https://" not in html.lower()


def test_workbench_content_identity_changes_with_catalog_content() -> None:
    catalog = build_decision_option_catalog()
    dependencies = build_decision_dependency_table()
    selections = build_decision_selection_template()
    identity = compute_workbench_content_identity(catalog, dependencies, selections)
    changed = catalog.copy()
    changed.loc[0, "PrimaryRisk"] += " changed"
    assert compute_workbench_content_identity(changed, dependencies, selections) != identity


def test_private_workbench_bundle_is_deterministic_and_valid() -> None:
    first = build_private_decision_workbench_bundle()
    second = build_private_decision_workbench_bundle()
    assert first["bundle_bytes"] == second["bundle_bytes"]
    result = validate_private_decision_workbench_bundle(first["bundle_bytes"])
    assert result["valid"]
    assert len(first["manifest"]) == 6
    assert first["manifest"]["HumanRows"].eq(0).all()
    assert not first["manifest"]["ContainsSelection"].any()


def test_bundle_html_tampering_fails_manifest_validation() -> None:
    bundle = build_private_decision_workbench_bundle()["bundle_bytes"]
    with zipfile.ZipFile(io.BytesIO(bundle), "r") as source:
        payloads = {name: source.read(name) for name in source.namelist()}
    payloads["protocol_decision_workbench.html"] = payloads[
        "protocol_decision_workbench.html"
    ].replace(b"BLOCKED", b"READY  ", 1)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            info.create_system = 3
            archive.writestr(info, payloads[name])
    result = validate_private_decision_workbench_bundle(buffer.getvalue())
    assert not result["valid"]
    assert "manifest_hash_mismatch" in result["failure_codes"]
    assert "invalid_workbench_html" in result["failure_codes"]


def test_human_gate_stays_closed() -> None:
    bundle = build_private_decision_workbench_bundle()
    status = decision_workbench_human_gate_status(bundle["selection_status"]).iloc[0]
    assert status["HumanParticipants"] == 0
    assert status["OptionsSelected"] == 0
    assert not bool(status["SelectionsComplete"])
    assert not bool(status["SubstantiveEvidenceVerified"])
    assert not bool(status["RecruitmentReady"])
    assert not bool(status["MinimumValidPerCellRegistered"])
    assert not bool(status["PlannedRecruitmentNSelected"])
    assert not bool(status["ConfirmatoryResultAvailable"])
    assert not bool(status["PublicSurfaceEnabled"])
