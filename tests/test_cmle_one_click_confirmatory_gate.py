"""Contracts for versioned directional CMLE confirmatory gate mechanics."""

from __future__ import annotations

import copy

import pandas as pd
import pytest

from validation.cmle_one_click_comprehension_readiness import run_cases

from mfrm_app.cmle_one_click_cognitive_interview import INSTRUMENT_VERSION
from mfrm_app.cmle_one_click_comprehension import HUMAN_STUDY_STATUS
from mfrm_app.cmle_one_click_confirmatory_gate import (
    ANCHOR_PRIMARY_CASE,
    BLOCKED_PRIMARY_CASES,
    DIRECTIONAL_RULES,
    build_confirmatory_assignment_schedule,
    build_instrument_version_record,
    build_wilson_planning_table,
    compute_instrument_content_identity,
    confirmatory_human_gate_status,
    evaluate_directional_confirmatory_gate,
    maximum_passing_errors,
    validate_instrument_version_ledger,
    wilson_upper_bound,
)


@pytest.fixture(scope="module")
def materials(tmp_path_factory):
    # Build synthetic instrument inputs from tracked code, not local study outputs.
    output = tmp_path_factory.mktemp("comprehension_materials")
    run_cases(output)
    return (
        output / "participant_task_packet.csv",
        output / "researcher_scoring_key.csv",
        output / "previews",
    )


def _key(materials) -> pd.DataFrame:
    return pd.read_csv(materials[1])


def _correct_responses(
    assignment: pd.DataFrame, key: pd.DataFrame
) -> pd.DataFrame:
    rows = assignment.merge(
        key[["CaseId", "Language", "ItemId", "CorrectOption"]],
        on=["CaseId", "Language"],
        validate="many_to_many",
    )
    return rows[
        ["InstrumentVersion", "ConfirmatorySlotId", "Language", "CaseId", "ItemId"]
    ].assign(ResponseCode=rows["CorrectOption"].to_numpy())


def test_instrument_identity_is_reproducible_and_mutation_sensitive(materials) -> None:
    ledger, manifest = build_instrument_version_record(
        *materials
    )
    repeated, repeated_manifest = build_instrument_version_record(
        *materials
    )
    assert ledger.equals(repeated)
    assert manifest.equals(repeated_manifest)
    assert validate_instrument_version_ledger(ledger)["valid"]

    baseline = ledger.iloc[0]["InstrumentContentSHA256"]
    task_mutation = compute_instrument_content_identity(
        materials[0].read_bytes() + b"synthetic-mutation",
        materials[1].read_bytes(),
        manifest,
    )
    assert task_mutation["InstrumentContentSHA256"] != baseline
    rules = copy.deepcopy(list(DIRECTIONAL_RULES))
    rules[0] = {**rules[0], "DangerousResponseCodes": ("ready", "synthetic")}
    rule_mutation = compute_instrument_content_identity(
        materials[0].read_bytes(), materials[1].read_bytes(), manifest, directional_rules=rules
    )
    assert rule_mutation["InstrumentContentSHA256"] != baseline


def test_version_ledger_rejects_alias_pooling_and_post_response_freeze(materials) -> None:
    ledger, _ = build_instrument_version_record(*materials)
    alias = pd.concat(
        [ledger, ledger.assign(InstrumentVersion="alias_version")], ignore_index=True
    )
    assert "duplicate_content_identity" in validate_instrument_version_ledger(alias)[
        "failure_codes"
    ]
    pooling = ledger.assign(PoolingWithDifferentContentAllowed=True)
    assert "cross_version_pooling_allowed" in validate_instrument_version_ledger(
        pooling
    )["failure_codes"]
    post = ledger.assign(HumanResponsesBeforeFreeze=1)
    assert "responses_before_freeze" in validate_instrument_version_ledger(post)[
        "failure_codes"
    ]


def test_confirmatory_assignment_has_two_primary_roles_and_language_parity() -> None:
    schedule = build_confirmatory_assignment_schedule(25)
    assert len(schedule) == 100
    assert schedule["ConfirmatorySlotId"].nunique() == 50
    assert schedule.groupby("ConfirmatorySlotId").size().eq(2).all()
    assert schedule.groupby(["ConfirmatorySlotId", "CaseRole"]).size().eq(1).all()
    assert schedule.loc[schedule["CaseRole"].eq("anchored_primary"), "CaseId"].eq(
        ANCHOR_PRIMARY_CASE
    ).all()
    blocked = schedule.loc[schedule["CaseRole"].eq("blocked_primary")]
    assert set(blocked["CaseId"]) == set(BLOCKED_PRIMARY_CASES)
    assert blocked.groupby(["Language", "CaseId"]).size().groupby("Language").apply(
        lambda values: int(values.max() - values.min()) <= 1
    ).all()
    parity = schedule.pivot_table(
        index=["SlotNumber", "CaseOrder"],
        columns="Language",
        values=["CaseRole", "CaseId"],
        aggfunc="first",
    )
    for field in ("CaseRole", "CaseId"):
        assert parity[(field, "en")].equals(parity[(field, "ja")])


def test_wilson_arithmetic_strict_boundary_and_familywise_sensitivity() -> None:
    assert wilson_upper_bound(0, 24) > 0.10
    assert wilson_upper_bound(0, 25) < 0.10
    assert maximum_passing_errors(24) == -1
    assert maximum_passing_errors(25) == 0
    assert wilson_upper_bound(0, 100) < wilson_upper_bound(1, 100)
    planning = build_wilson_planning_table()
    row24 = planning.set_index("ValidEligibleSlots").loc[24]
    row25 = planning.set_index("ValidEligibleSlots").loc[25]
    assert row24["PrimaryMaximumPassingErrors"] == -1
    assert row25["PrimaryMaximumPassingErrors"] == 0
    assert row25["FamilywiseZeroErrorUpper"] > 0.10
    assert planning["Interpretation"].eq(
        "arithmetic_only_not_sample_size_recommendation"
    ).all()


def test_zero_error_n25_passes_primary_not_familywise_sensitivity(materials) -> None:
    key = _key(materials)
    assignment = build_confirmatory_assignment_schedule(25)
    result = evaluate_directional_confirmatory_gate(
        _correct_responses(assignment, key),
        key,
        assignment,
        registered_instrument_version=INSTRUMENT_VERSION,
        minimum_valid_per_cell=25,
    )
    cells = result["cell_summary"]
    summary = result["study_summary"].iloc[0]
    assert len(cells) == 12
    assert cells["ValidEligibleSlots"].eq(25).all()
    assert cells["DangerousErrors"].eq(0).all()
    assert cells["PrimaryPass"].all()
    assert not cells["FamilywiseSensitivityPass"].any()
    assert bool(summary["OverallPrimaryGatePassed"])
    assert summary["PrimaryPassingCells"] == 12
    assert summary["FamilywiseSensitivityPassingCells"] == 0


def test_directional_error_does_not_pool_language_or_non_dangerous_error(materials) -> None:
    key = _key(materials)
    assignment = build_confirmatory_assignment_schedule(25)
    base = _correct_responses(assignment, key)

    dangerous = base.copy()
    ja_slot = assignment.loc[
        assignment["Language"].eq("ja")
        & assignment["CaseRole"].eq("anchored_primary")
    ].iloc[0]
    mask = (
        dangerous["ConfirmatorySlotId"].eq(ja_slot["ConfirmatorySlotId"])
        & dangerous["CaseId"].eq(ANCHOR_PRIMARY_CASE)
        & dangerous["ItemId"].eq("public_sharing")
    )
    dangerous.loc[mask, "ResponseCode"] = "permitted"
    result = evaluate_directional_confirmatory_gate(
        dangerous,
        key,
        assignment,
        registered_instrument_version=INSTRUMENT_VERSION,
        minimum_valid_per_cell=25,
    )
    cells = result["cell_summary"].set_index(["Language", "Domain"])
    assert cells.loc[("ja", "private_archive_public_sharing"), "DangerousErrors"] == 1
    assert not bool(cells.loc[("ja", "private_archive_public_sharing"), "PrimaryPass"])
    assert cells.loc[("en", "private_archive_public_sharing"), "DangerousErrors"] == 0
    assert bool(cells.loc[("en", "private_archive_public_sharing"), "PrimaryPass"])
    assert not bool(result["study_summary"].iloc[0]["OverallPrimaryGatePassed"])

    conservative = base.copy()
    en_slot = assignment.loc[
        assignment["Language"].eq("en")
        & assignment["CaseRole"].eq("anchored_primary")
    ].iloc[0]
    rounding = (
        conservative["ConfirmatorySlotId"].eq(en_slot["ConfirmatorySlotId"])
        & conservative["CaseId"].eq(ANCHOR_PRIMARY_CASE)
        & conservative["ItemId"].eq("rounding_vignette")
    )
    conservative.loc[rounding, "ResponseCode"] = "cannot_decide"
    conservative_result = evaluate_directional_confirmatory_gate(
        conservative,
        key,
        assignment,
        registered_instrument_version=INSTRUMENT_VERSION,
        minimum_valid_per_cell=25,
    )
    rounding_cell = conservative_result["cell_summary"].set_index(
        ["Language", "Domain"]
    ).loc[("en", "rounded_mnsq_reclassification")]
    assert rounding_cell["AnyComprehensionErrors"] == 1
    assert rounding_cell["DangerousErrors"] == 0
    assert bool(rounding_cell["PrimaryPass"])


def test_invalid_or_mixed_version_slot_is_retained_and_reduces_denominator(materials) -> None:
    key = _key(materials)
    assignment = build_confirmatory_assignment_schedule(25)
    responses = _correct_responses(assignment, key)
    target_slot = responses.loc[responses["Language"].eq("en"), "ConfirmatorySlotId"].iloc[0]
    responses.loc[
        responses["ConfirmatorySlotId"].eq(target_slot), "InstrumentVersion"
    ] = "unregistered_version"
    result = evaluate_directional_confirmatory_gate(
        responses,
        key,
        assignment,
        registered_instrument_version=INSTRUMENT_VERSION,
        minimum_valid_per_cell=25,
    )
    slot = result["slot_audit"].set_index("ConfirmatorySlotId").loc[target_slot]
    assert not bool(slot["ValidSlot"])
    assert "mixed_or_unregistered_instrument_version" in slot["InvalidReasons"]
    en_cells = result["cell_summary"].loc[
        result["cell_summary"]["Language"].eq("en")
    ]
    assert en_cells["ValidEligibleSlots"].eq(24).all()
    assert not en_cells["PrimaryPass"].any()
    assert result["study_summary"].iloc[0]["InvalidAttemptedSlots"] == 1


def test_minimum_n_is_mandatory_and_human_gate_remains_closed(materials) -> None:
    key = _key(materials)
    assignment = build_confirmatory_assignment_schedule(1)
    try:
        evaluate_directional_confirmatory_gate(
            _correct_responses(assignment, key),
            key,
            assignment,
            registered_instrument_version=INSTRUMENT_VERSION,
            minimum_valid_per_cell=0,
        )
    except ValueError as exc:
        assert "prospectively set" in str(exc)
    else:
        raise AssertionError("Confirmatory gate accepted an unregistered minimum n.")
    status = confirmatory_human_gate_status().iloc[0]
    assert status["HumanStudyStatus"] == HUMAN_STUDY_STATUS
    assert status["HumanParticipants"] == 0
    assert not bool(status["PilotResultAvailable"])
    assert not bool(status["ConfirmatoryResultAvailable"])
    assert not bool(status["MinimumValidPerCellRegistered"])
    assert not bool(status["LanguageEquivalenceAvailable"])
    assert not bool(status["PublicSurfaceEnabled"])
