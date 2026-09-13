from __future__ import annotations

import io
import zipfile

import pandas as pd

from mfrm_app.cmle_one_click_confirmatory_estimand import (
    SCI01_OPTION_IDS,
    SCI01_SELECTION_COLUMNS,
    build_private_sci01_bundle,
    build_sci01_downstream_impacts,
    build_sci01_estimand_definitions,
    build_sci01_selection_template,
    compute_sci01_content_identity,
    project_sci01_estimands,
    sci01_human_gate_status,
    validate_private_sci01_bundle,
    validate_sci01_html,
    validate_sci01_selection,
)
from mfrm_app.cmle_one_click_confirmatory_protocol import (
    build_partial_identification_sensitivity,
)


def surface() -> pd.DataFrame:
    return build_partial_identification_sensitivity()


def test_definitions_are_three_unranked_nonrecommended_options() -> None:
    definitions = build_sci01_estimand_definitions()
    assert tuple(definitions["OptionId"]) == SCI01_OPTION_IDS
    assert tuple(definitions["DisplayOrder"]) == (1, 2, 3)
    assert not definitions["DisplayOrderIsRank"].any()
    assert not definitions["AutoSelectable"].any()
    assert definitions["RepositoryRecommendation"].eq("none").all()
    assert not bool(definitions.iloc[1]["IsDirectlyCalculableFromFrozenSurface"])
    assert "without proof" in definitions.iloc[1]["ForbiddenClaim"]


def test_projection_has_exact_same_data_three_way_rows() -> None:
    source = surface()
    projection = project_sci01_estimands(source)
    assert len(source) == 72
    assert len(projection) == 216
    assert projection.groupby("SourceScenarioRow").size().eq(3).all()
    assert set(projection["OptionId"]) == set(SCI01_OPTION_IDS)
    assert not projection["DecisionUsesDisplayedRounding"].any()
    assert not projection["SampleSizeSelected"].any()
    assert projection["HumanParticipants"].eq(0).all()


def test_projection_formulas_use_integer_counts_and_raw_threshold() -> None:
    projection = project_sci01_estimands(surface())
    for _, group in projection.groupby("SourceScenarioRow"):
        d = int(group.iloc[0]["DangerousErrors"])
        v = int(group.iloc[0]["ValidEligibleN"])
        i = int(group.iloc[0]["InvalidEligibleN"])
        a = group.loc[group["OptionId"].eq(SCI01_OPTION_IDS[0])].iloc[0]
        b = group.loc[group["OptionId"].eq(SCI01_OPTION_IDS[1])].iloc[0]
        c = group.loc[group["OptionId"].eq(SCI01_OPTION_IDS[2])].iloc[0]
        assert a["ReportedPointRaw"] == d / v
        assert b["ReportedPointRaw"] == (d + i) / (v + i)
        assert c["ReportedLowerRaw"] == d / (v + i)
        assert c["ReportedUpperRaw"] == (d + i) / (v + i)
        assert bool(a["GoverningWilsonUpper95Raw"] < a["DecisionThresholdRaw"]) == (a["ProjectedDecision"] == "pass")


def test_same_counts_can_change_interpretation_without_rounding() -> None:
    projection = project_sci01_estimands(surface())
    candidates = projection.groupby("SourceScenarioRow").filter(
        lambda group: set(group["ProjectedDecision"]) >= {"pass", "proxy_fail", "observed_pass_not_robust"}
    )
    assert not candidates.empty
    row_id = int(candidates.iloc[0]["SourceScenarioRow"])
    group = projection.loc[projection["SourceScenarioRow"].eq(row_id)]
    assert group[["DangerousErrors", "ValidEligibleN", "InvalidEligibleN"]].nunique().eq(1).all()
    assert set(group["ProjectedDecision"]) == {"pass", "proxy_fail", "observed_pass_not_robust"}


def test_option_b_is_never_labeled_final_scheduled_slot_calculation() -> None:
    projection = project_sci01_estimands(surface())
    option_b = projection.loc[projection["OptionId"].eq(SCI01_OPTION_IDS[1])]
    assert not option_b["FinalScheduledSlotEstimandCalculated"].any()
    assert option_b["ProjectionRole"].eq("all_invalid_dangerous_diagnostic_proxy").all()
    assert option_b["InterpretationLimit"].str.contains("Not the scheduled-slot estimand", regex=False).all()


def test_downstream_map_is_complete_and_nonquantitative() -> None:
    impacts = build_sci01_downstream_impacts()
    assert len(impacts) == 36
    assert impacts.groupby("OptionId").size().eq(12).all()
    assert not impacts["ImpactMagnitudeEstimated"].any()
    assert not impacts["DirectionKnown"].any()
    assert not impacts["DecisionResolved"].any()


def test_blank_selection_is_valid_incomplete_and_fail_closed() -> None:
    selection = build_sci01_selection_template()
    result = validate_sci01_selection(selection)
    assert tuple(selection.columns) == SCI01_SELECTION_COLUMNS
    assert result["valid"]
    assert not result["complete"]
    assert not result["RecruitmentReady"]


def test_resolved_selection_requires_external_owner_evidence_and_timing() -> None:
    selection = build_sci01_selection_template()
    selection.loc[0, "DecisionStatus"] = "resolved_prospectively"
    selection.loc[0, "SelectedOptionId"] = SCI01_OPTION_IDS[2]
    incomplete = validate_sci01_selection(selection)
    assert not incomplete["valid"]
    assert "missing_owner_or_evidence_reference" in incomplete["failure_codes"]
    assert "missing_invalidity_code_map" in incomplete["failure_codes"]
    for column in ("ScientificOwnerRole", "AssumptionSetReference", "RationaleReference", "OwnerAttestationReference", "InvalidityCodeMapReference"):
        selection.loc[0, column] = "SYNTHETIC-EXTERNAL"
    selection.loc[0, "DecisionRecordedBeforeOutcomes"] = True
    complete = validate_sci01_selection(selection)
    assert complete["valid"] and complete["complete"]
    assert not complete["RecruitmentReady"]
    selection.loc[0, "HumanOutcomesInspectedBeforeDecision"] = True
    assert "decision_not_prospectively_locked" in validate_sci01_selection(selection)["failure_codes"]


def test_html_is_read_only_complete_and_nonrecommending() -> None:
    bundle = build_private_sci01_bundle(surface())
    result = validate_sci01_html(bundle["html"], bundle["content_sha256"])
    assert result["valid"]
    lower = bundle["html"].lower()
    assert "<form" not in lower and "<script" not in lower
    assert "http://" not in lower and "https://" not in lower
    assert "scheduled-slot estimand そのものではありません" in bundle["html"]


def test_content_identity_changes_with_any_definition_change() -> None:
    definitions = build_sci01_estimand_definitions(); projection = project_sci01_estimands(surface())
    impacts = build_sci01_downstream_impacts(); selection = build_sci01_selection_template()
    identity = compute_sci01_content_identity(definitions, projection, impacts, selection)
    changed = definitions.copy(); changed.loc[0, "PrimaryRisk"] += " changed"
    assert compute_sci01_content_identity(changed, projection, impacts, selection) != identity


def test_private_bundle_is_deterministic_valid_and_unselected() -> None:
    first = build_private_sci01_bundle(surface()); second = build_private_sci01_bundle(surface())
    assert first["bundle_bytes"] == second["bundle_bytes"]
    assert validate_private_sci01_bundle(first["bundle_bytes"])["valid"]
    assert len(first["manifest"]) == 6
    assert first["manifest"]["HumanRows"].eq(0).all()
    assert not first["manifest"]["ContainsSelection"].any()


def test_bundle_tampering_fails_hash_and_html_contract() -> None:
    bundle = build_private_sci01_bundle(surface())["bundle_bytes"]
    with zipfile.ZipFile(io.BytesIO(bundle), "r") as source:
        payloads = {name: source.read(name) for name in source.namelist()}
    payloads["sci01_estimand_comparison.html"] = payloads["sci01_estimand_comparison.html"].replace(b"BLOCKED", b"READY  ", 1)
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0)); info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16; info.create_system = 3
            archive.writestr(info, payloads[name])
    result = validate_private_sci01_bundle(buffer.getvalue())
    assert not result["valid"]
    assert "manifest_hash_mismatch" in result["failure_codes"]
    assert "invalid_sci01_html" in result["failure_codes"]


def test_human_gate_stays_closed_even_if_sci01_fixture_is_complete() -> None:
    selection = build_sci01_selection_template()
    status = validate_sci01_selection(selection)
    gate = sci01_human_gate_status(status).iloc[0]
    assert gate["HumanParticipants"] == 0
    assert not gate["SCI01Complete"]
    assert not gate["SubstantiveEvidenceVerified"]
    assert not gate["RecruitmentReady"]
    assert not gate["SampleSizeSelected"]
    assert not gate["ConfirmatoryOutcomesAvailable"]
    assert not gate["PublicSurfaceEnabled"]
