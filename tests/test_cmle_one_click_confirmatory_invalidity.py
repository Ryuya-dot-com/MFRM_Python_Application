from __future__ import annotations

import io
import zipfile

from mfrm_app.cmle_one_click_confirmatory_invalidity import (
    BUNDLE_MEMBERS, COMPOSITE_FRACTIONS, SELECTION_COLUMNS,
    build_adjudication_template, build_allowed_state_catalog,
    build_composite_fraction_sensitivity, build_invalidity_code_taxonomy,
    build_owner_review_matrix, build_private_bundle, content_identity,
    validate_adjudication, validate_html, validate_private_bundle,
)
from mfrm_app.cmle_one_click_confirmatory_protocol import (
    REGISTERED_INVALIDITY_CODES, build_invalidity_reason_codebook,
    build_partial_identification_sensitivity,
)


def complete_fixture():
    frame = build_adjudication_template()
    for index, row in frame.iterrows():
        code = row["InvalidityReasonCode"]
        if code == "none": continue
        if code == "other_requires_protocol_amendment":
            frame.loc[index, "AdjudicationStatus"] = "amendment_required"; continue
        frame.loc[index, "AdjudicationStatus"] = "resolved_prospectively"
        frame.loc[index, "TargetPopulationDisposition"] = "in_target"
        frame.loc[index, "CompositeEventDisposition"] = "sensitivity_only"
        frame.loc[index, "PrimaryAnalysisDisposition"] = "exclude_preregistered"
        frame.loc[index, "SensitivitySetDisposition"] = "code_specific_tipping"
        frame.loc[index, "DataUseAuthorityDisposition"] = "not_applicable"
        frame.loc[index, "ScientificOwnerRole"] = "external scientific owner"
        frame.loc[index, ["StatisticalReviewReference", "DomainReviewReference", "RationaleReference"]] = "SYNTHETIC-EXTERNAL"
        frame.loc[index, "RecordedBeforeOutcomes"] = True
        if "withdrawn" in code or code == "declined_before_start":
            frame.loc[index, "DataUseAuthorityDisposition"] = "authorized_use"
            frame.loc[index, "EthicsAuthorityReference"] = "SYNTHETIC-ETHICS"
        if code == "accessibility_barrier":
            frame.loc[index, "ScientificOwnerRole"] = "external accessibility owner"
            frame.loc[index, "SensitivitySetDisposition"] = "separate_stratum_required"
    return frame


def test_taxonomy_preserves_all_codes_without_auto_danger() -> None:
    taxonomy = build_invalidity_code_taxonomy(build_invalidity_reason_codebook())
    assert tuple(taxonomy["InvalidityReasonCode"]) == REGISTERED_INVALIDITY_CODES
    assert len(taxonomy) == 15
    assert not taxonomy["AutomaticallyDangerous"].any()
    assert taxonomy["AutomaticallySafe"].sum() == 1
    assert taxonomy["RepositoryRecommendation"].eq("none").all()


def test_allowed_states_are_unranked_and_nonrecommended() -> None:
    states = build_allowed_state_catalog()
    assert states["Dimension"].nunique() == 5
    assert not states["DisplayOrderIsRank"].any()
    assert not states["RepositoryRecommended"].any()


def test_owner_matrix_requires_ethics_and_accessibility_roles() -> None:
    taxonomy = build_invalidity_code_taxonomy(build_invalidity_reason_codebook())
    owners = build_owner_review_matrix(taxonomy)
    withdrawal = owners.loc[owners["InvalidityReasonCode"].eq("withdrawn_after_response_lock")]
    access = owners.loc[owners["InvalidityReasonCode"].eq("accessibility_barrier")]
    assert "ethics_authority" in set(withdrawal["OwnerRole"])
    assert "accessibility_owner" in set(access["OwnerRole"])
    assert not owners["ReviewCompleted"].any()
    assert not owners["RepositoryMaySelfApprove"].any()


def test_blank_template_is_valid_incomplete_with_only_none_fixed() -> None:
    template = build_adjudication_template(); result = validate_adjudication(template)
    assert tuple(template.columns) == SELECTION_COLUMNS
    assert result["valid"] and not result["complete"]
    assert result["resolved_code_count"] == 0
    assert not result["RecruitmentReady"]
    assert template["AdjudicationStatus"].eq("fixed_parent_contract").sum() == 1


def test_guardrails_reject_auto_duplicate_and_accessibility_pooling() -> None:
    frame = complete_fixture(); assert validate_adjudication(frame)["complete"]
    frame.loc[frame["InvalidityReasonCode"].eq("duplicate_study_record"), "CompositeEventDisposition"] = "dangerous_composite_event"
    assert "duplicate_cannot_be_dangerous_event" in validate_adjudication(frame)["failure_codes"]
    frame = complete_fixture(); mask = frame["InvalidityReasonCode"].eq("accessibility_barrier")
    frame.loc[mask, "SensitivitySetDisposition"] = "code_specific_tipping"
    assert "accessibility_silent_pooling_forbidden" in validate_adjudication(frame)["failure_codes"]
    frame = complete_fixture(); frame.loc[1, "RepositoryAutoClassificationAllowed"] = True
    assert "repository_auto_classification_forbidden" in validate_adjudication(frame)["failure_codes"]


def test_withdrawal_requires_ethics_and_preoutcome_timing() -> None:
    frame = complete_fixture(); mask = frame["InvalidityReasonCode"].eq("withdrawn_after_response_lock")
    frame.loc[mask, "EthicsAuthorityReference"] = ""
    assert "missing_ethics_data_use_authority" in validate_adjudication(frame)["failure_codes"]
    frame = complete_fixture(); frame.loc[mask, "HumanOutcomesInspectedBeforeDecision"] = True
    assert "adjudication_not_prospectively_blinded" in validate_adjudication(frame)["failure_codes"]


def test_fraction_sensitivity_is_exact_72_by_5_and_unrounded() -> None:
    source = build_partial_identification_sensitivity()
    sensitivity = build_composite_fraction_sensitivity(source)
    assert len(source) == 72 and len(sensitivity) == 360
    assert tuple(sorted(sensitivity["CompositeFractionOfInvalidRegistered"].unique())) == COMPOSITE_FRACTIONS
    assert sensitivity.groupby("SourceScenarioRow").size().eq(5).all()
    assert not sensitivity["DecisionUsesDisplayedRounding"].any()
    assert not sensitivity["InfersCodeFrequency"].any()
    assert not sensitivity["SelectsCompositeDefinition"].any()
    assert not sensitivity["FinalScheduledSlotEstimandCalculated"].any()


def test_integer_allocation_and_raw_decision_contract() -> None:
    sensitivity = build_composite_fraction_sensitivity(build_partial_identification_sensitivity())
    for row in sensitivity.itertuples():
        expected = int((row.InvalidEligibleN * row.CompositeFractionOfInvalidRegistered + .5) // 1)
        assert row.InvalidClassifiedCompositeN == expected
        assert (row.CompositeWilsonUpper95Raw < row.DecisionThresholdRaw) == (row.ProjectedDecision == "proxy_pass")


def test_html_and_content_identity_are_fail_closed() -> None:
    bundle = build_private_bundle(build_invalidity_reason_codebook(), build_partial_identification_sensitivity())
    assert validate_html(bundle["html"], bundle["identity"])["valid"]
    assert "無効は危険事象と同義ではありません" in bundle["html"]
    changed = bundle["taxonomy"].copy(); changed.loc[1, "HardGuardrail"] += " changed"
    assert content_identity(changed, bundle["states"], bundle["template"], bundle["owners"], bundle["sensitivity"]) != bundle["identity"]


def test_bundle_is_deterministic_valid_and_tamper_evident() -> None:
    codebook=build_invalidity_reason_codebook(); surface=build_partial_identification_sensitivity()
    first=build_private_bundle(codebook,surface); second=build_private_bundle(codebook,surface)
    assert first["bundle_bytes"] == second["bundle_bytes"]
    assert validate_private_bundle(first["bundle_bytes"])["valid"]
    assert len(first["manifest"]) == len(BUNDLE_MEMBERS)-1
    with zipfile.ZipFile(io.BytesIO(first["bundle_bytes"]),"r") as source:
        payloads={name:source.read(name) for name in source.namelist()}
    payloads["invalidity_adjudication_workbench.html"] = payloads["invalidity_adjudication_workbench.html"].replace(b"BLOCKED",b"READY  ",1)
    buffer=io.BytesIO()
    with zipfile.ZipFile(buffer,"w",compression=zipfile.ZIP_DEFLATED) as archive:
        for name in sorted(payloads):
            info=zipfile.ZipInfo(name,(1980,1,1,0,0,0)); info.compress_type=zipfile.ZIP_DEFLATED; info.external_attr=0o100644<<16; info.create_system=3; archive.writestr(info,payloads[name])
    result=validate_private_bundle(buffer.getvalue())
    assert not result["valid"] and "manifest_hash" in result["failure_codes"] and "invalid_html" in result["failure_codes"]
