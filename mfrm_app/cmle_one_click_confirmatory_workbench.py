"""Private, non-recommending decision workbench for the CMLE protocol."""

from __future__ import annotations

import hashlib
from html import escape
import io
import math
from typing import Mapping
import zipfile

import pandas as pd

from mfrm_app.cmle_one_click_confirmatory_protocol import (
    build_protocol_decision_register,
)


WORKBENCH_SCHEMA_VERSION = "cmle_confirmatory_decision_workbench_v1"
DECISION_IDS = (
    "SCI-01", "SCI-02", "SCI-03", "SCI-04", "SCI-05", "SCI-06", "SCI-07",
    "SCI-08", "SCI-09", "OPS-01", "ETH-01", "ETH-02", "GOV-01",
)
SELECTION_COLUMNS = (
    "SchemaVersion", "DecisionId", "Decision", "OwnerRole", "DecisionStatus",
    "SelectedOptionId", "NumericValue", "NumericUnit", "AssumptionSetReference",
    "RationaleReference", "OwnerAttestationReference", "ExternalEvidenceReference",
    "EvidenceSource", "DecisionRecordedBeforeOutcomes",
    "HumanOutcomesInspectedBeforeDecision", "RepositoryGeneratedSelectionAllowed",
)
WORKBENCH_BUNDLE_MEMBERS = (
    "BUNDLE_MANIFEST.csv",
    "START_HERE.md",
    "decision_dependency_table.csv",
    "decision_option_catalog.csv",
    "decision_selection_template.csv",
    "decision_workbench_summary.csv",
    "protocol_decision_workbench.html",
)
PREREQUISITES = {
    "SCI-01": (),
    "SCI-04": ("SCI-01",),
    "SCI-05": ("SCI-01",),
    "SCI-08": ("SCI-01", "SCI-04"),
    "SCI-09": ("SCI-01",),
    "SCI-06": ("SCI-01", "SCI-05"),
    "SCI-07": ("SCI-01", "SCI-04"),
    "OPS-01": ("SCI-01", "SCI-04", "SCI-05", "SCI-07"),
    "GOV-01": ("SCI-01", "SCI-05"),
    "SCI-02": ("SCI-01", "SCI-04", "SCI-05", "SCI-06", "SCI-07", "SCI-08", "SCI-09"),
    "SCI-03": ("SCI-02", "SCI-05", "SCI-06", "SCI-09"),
    "ETH-01": ("SCI-02", "SCI-03", "OPS-01", "GOV-01"),
    "ETH-02": ("ETH-01", "GOV-01"),
}
_TOPOLOGICAL_ORDER = {
    "SCI-01": 1, "SCI-04": 2, "SCI-05": 3, "SCI-08": 4, "SCI-09": 5,
    "SCI-06": 6, "SCI-07": 7, "OPS-01": 8, "GOV-01": 9, "SCI-02": 10,
    "SCI-03": 11, "ETH-01": 12, "ETH-02": 13,
}
_STAGES = {
    "SCI-01": (1, "1_estimand"),
    "SCI-04": (2, "2_missingness_multiplicity_language_accessibility"),
    "SCI-05": (2, "2_missingness_multiplicity_language_accessibility"),
    "SCI-08": (2, "2_missingness_multiplicity_language_accessibility"),
    "SCI-09": (2, "2_missingness_multiplicity_language_accessibility"),
    "SCI-06": (3, "3_cluster_mechanism_operations_governance"),
    "SCI-07": (3, "3_cluster_mechanism_operations_governance"),
    "OPS-01": (3, "3_cluster_mechanism_operations_governance"),
    "GOV-01": (3, "3_cluster_mechanism_operations_governance"),
    "SCI-02": (4, "4_sample_size"),
    "SCI-03": (5, "5_recruitment"),
    "ETH-01": (6, "6_ethics_consent"),
    "ETH-02": (6, "6_ethics_consent"),
}


def _truth(value: object) -> bool:
    return value if isinstance(value, bool) else str(value).strip().lower() in {"true", "1", "yes"}


def _csv_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False, float_format="%.17g").encode("utf-8")


def _sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _blocking_register() -> pd.DataFrame:
    register = build_protocol_decision_register()
    result = register.loc[register["Blocking"]].copy()
    result["DecisionId"] = pd.Categorical(
        result["DecisionId"], categories=DECISION_IDS, ordered=True
    )
    return result.sort_values("DecisionId").reset_index(drop=True)


_OPTION_SPECS = (
    ("SCI-01", "SCI-01-A_observed_valid_estimand", "Observed-valid estimand", "Estimate risk among valid eligible records and report denominator limits.", "Simple and directly aligned with the frozen Wilson calculation.", "Can be biased under informative invalidity and may not represent scheduled users.", "Target-population definition and missingness justification", False),
    ("SCI-01", "SCI-01-B_scheduled_slot_composite", "Scheduled-slot composite", "Treat specified invalid outcomes as composite failures in a scheduled-slot estimand.", "Protects against some silent exclusion and denominator drift.", "Composite failures may conflate comprehension, technology, withdrawal, and access barriers.", "Composite-event definition and stakeholder justification", False),
    ("SCI-01", "SCI-01-C_dual_estimand_with_bounds", "Dual estimand with bounds", "Report observed-valid and scheduled-slot bounds as distinct estimands.", "Makes invalid-outcome uncertainty visible without claiming identification.", "More complex interpretation and still requires a governing decision rule.", "Dual-estimand protocol and adjudication rule", False),
    ("SCI-02", "SCI-02-A_per_cell_assurance", "Per-cell assurance target", "Choose numeric n from a prospectively justified per-cell pass-assurance target.", "Transparent link to the frozen primary cell rule.", "May underprotect joint claims and remains assumption-sensitive.", "Assumed rates, assurance target, and sustained-search rule", True),
    ("SCI-02", "SCI-02-B_familywise_assurance", "Familywise assurance target", "Choose numeric n using a simultaneous or conservative familywise target.", "Aligns planning with a joint safety claim.", "Can require substantially larger n and depends on multiplicity choice.", "Familywise method, target, and assumed rates", True),
    ("SCI-02", "SCI-02-C_cluster_simulation_precision", "Cluster simulation and precision", "Choose numeric n from registered cluster-aware simulation and precision criteria.", "Can reflect dependence, imbalance, and discontinuous decisions.", "Model misspecification may give false precision and needs external calibration.", "Frozen simulation model, precision target, and sensitivity envelope", True),
    ("SCI-03", "SCI-03-A_fixed_inflation", "Fixed recruitment inflation", "Inflate selected valid n by one prospectively justified fixed allowance.", "Operationally simple and auditable.", "Fails when retention differs by outcome, site, language, or access need.", "Retention evidence and fixed inflation formula", True),
    ("SCI-03", "SCI-03-B_retention_assurance", "Retention assurance", "Set recruitment to attain valid n with a registered retention-assurance probability.", "Explicitly separates valid target from recruitment uncertainty.", "Homogeneous retention assumptions may be wrong or outcome-dependent.", "Retention model, assurance target, and monitoring plan", True),
    ("SCI-03", "SCI-03-C_blinded_denominator_topup", "Blinded denominator top-up", "Permit prospective top-up using blinded validity denominators only.", "Can protect valid counts without viewing dangerous-error outcomes.", "Operational changes can still induce selection and require a strict cap.", "Blinding proof, top-up cap, trigger, and audit procedure", True),
    ("SCI-04", "SCI-04-A_primary_all_cells_sensitivity_only", "All primary cells; familywise sensitivity", "Keep every per-cell primary gate and report familywise correction as sensitivity.", "Preserves the existing frozen primary interpretation.", "Does not provide a simultaneous 95% familywise coverage claim.", "Claim language and multiplicity justification", False),
    ("SCI-04", "SCI-04-B_bonferroni_governing", "Bonferroni governing rule", "Make the registered Bonferroni sensitivity the governing joint gate.", "Provides transparent conservative multiplicity control.", "May be very conservative and increase n substantially.", "Family definition and governing confidence specification", False),
    ("SCI-04", "SCI-04-C_hierarchical_gatekeeping", "Hierarchical gatekeeping", "Preorder claims and test them through a frozen hierarchy.", "Can align error control with scientifically prioritized claims.", "Ordering can be contestable and post hoc hierarchy changes invalidate control.", "Claim hierarchy, alpha allocation, and lock procedure", False),
    ("SCI-05", "SCI-05-A_complete_case_plus_bounds", "Complete case plus bounds", "Use valid records for the gate and require registered invalid-outcome bounds.", "Retains current calculation while exposing nonidentification.", "A passing observed gate may remain nonrobust and hard to communicate.", "Bound rule, robustness threshold, and denominator audit", False),
    ("SCI-05", "SCI-05-B_all_invalid_dangerous_primary", "All invalid dangerous primary", "Count every invalid eligible outcome as dangerous in the governing analysis.", "Strong protection against outcome-dependent exclusion.", "May be excessively conservative and mix unrelated failure processes.", "Eligibility and invalidity definitions plus stakeholder justification", False),
    ("SCI-05", "SCI-05-C_preregistered_tipping_models", "Preregistered tipping models", "Use multiple frozen missingness models and tipping thresholds.", "Shows when conclusions change across plausible MNAR mechanisms.", "Models are not identified and can invite assumption shopping.", "Model set, parameter ranges, and governing tipping rule", False),
    ("SCI-06", "SCI-06-A_cluster_robust_gee", "Cluster-robust GEE", "Use marginal modeling with prospectively defined cluster-robust inference.", "Targets a population-average effect with familiar robust variance.", "Few or imbalanced clusters can invalidate asymptotics.", "Cluster definition, minimum clusters, and small-sample correction", False),
    ("SCI-06", "SCI-06-B_hierarchical_logistic", "Hierarchical logistic model", "Model participant/site/moderator heterogeneity through random effects.", "Can estimate heterogeneity and partial pooling explicitly.", "Distributional assumptions and sparse clusters can dominate results.", "Model, priors or estimation rules, diagnostics, and sensitivity", False),
    ("SCI-06", "SCI-06-C_cluster_bootstrap", "Cluster bootstrap", "Resample prospectively defined independent clusters for uncertainty.", "Fewer parametric assumptions about within-cluster dependence.", "Needs enough exchangeable clusters and is awkward for discontinuous gates.", "Resampling unit, replicate count, failure handling, and CI rule", False),
    ("SCI-07", "SCI-07-A_pooled_with_heterogeneity_sensitivity", "Pooled with heterogeneity sensitivity", "Keep pooled blocked cells but require mechanism-specific descriptive sensitivity.", "Lower sample burden and continuity with the frozen pooled gate.", "Can conceal an unacceptable mechanism, as synthetic stress demonstrated.", "Pooling rationale and mechanism-specific warning thresholds", False),
    ("SCI-07", "SCI-07-B_mechanism_coprimary", "Mechanism-specific co-primary", "Require each blocked mechanism to satisfy a prospectively defined gate.", "Directly prevents compensation across failure mechanisms.", "Multiplies cells, increases n, and complicates multiplicity.", "Mechanism estimands, minima, and multiplicity rule", False),
    ("SCI-07", "SCI-07-C_hierarchical_partial_pooling", "Hierarchical partial pooling", "Estimate pooled and mechanism-specific risk in one hierarchical model.", "Shares information while retaining mechanism heterogeneity.", "Shrinkage can mask extremes and depends on model assumptions.", "Hierarchical model, tail-risk rule, and sensitivity", False),
    ("SCI-08", "SCI-08-A_separate_pass_no_equivalence", "Separate pass; no equivalence", "Require separate language passes and make no equivalence claim.", "Matches current evidence boundaries and is straightforward.", "Does not answer whether versions perform similarly.", "Language-specific estimands and claim restriction", False),
    ("SCI-08", "SCI-08-B_equivalence_noninferiority", "Equivalence or noninferiority", "Add a prospectively powered cross-language margin-based objective.", "Directly addresses a prespecified similarity claim.", "Margin choice is difficult and substantially increases sample needs.", "Margin justification, effect scale, and powered analysis", True),
    ("SCI-08", "SCI-08-C_invariance_plus_separate_pass", "Invariance plus separate pass", "Combine separate safety gates with measurement-invariance or DIF analyses.", "Examines comparability without allowing compensation across languages.", "Invariance does not prove semantic equivalence and models may be underpowered.", "Translation review, invariance model, DIF thresholds, and sample plan", False),
    ("SCI-09", "SCI-09-A_stratified_minima", "Stratified minima", "Set minimum valid counts for registered access, device, and proficiency strata.", "Makes underrepresented conditions visible and noncompensatory.", "Can sharply increase recruitment and requires defensible strata.", "Population rationale, strata definitions, and numeric minima", True),
    ("SCI-09", "SCI-09-B_universal_qa_monitored_strata", "Universal QA with monitored strata", "Use broad accessibility QA and monitor strata descriptively.", "Lower recruitment burden and broad implementation focus.", "Rare barriers may be missed and descriptive monitoring may not protect claims.", "QA protocol, monitoring thresholds, and escalation rule", False),
    ("SCI-09", "SCI-09-C_accessibility_coprimary", "Accessibility co-primary", "Make selected accessibility conditions co-primary confirmation targets.", "Provides strong evidence for explicitly scoped access claims.", "Needs specialized recruitment, accommodations, and larger n.", "Co-primary conditions, numeric minima, accommodations, and analysis", True),
    ("OPS-01", "OPS-01-A_no_outcome_interim", "No outcome interim", "Allow only integrity checks until the frozen final analysis.", "Minimizes alpha and behavior changes from interim outcomes.", "Cannot respond early to unexpectedly high dangerous-error rates.", "Integrity-monitoring and emergency-stop protocol", False),
    ("OPS-01", "OPS-01-B_blinded_information_monitoring", "Blinded information monitoring", "Monitor denominators, clusters, and validity without dangerous-error outcomes.", "Supports operational correction while limiting outcome influence.", "Blinded metrics can still reveal group patterns and affect selection.", "Blinding boundary, allowed metrics, and change-control log", False),
    ("OPS-01", "OPS-01-C_independent_monitoring", "Independent monitoring", "Use an independent body with frozen access and stopping rules.", "Permits safety oversight separated from investigators.", "Adds governance complexity and requires genuine independence.", "Charter, membership authority, access controls, and stopping rule", False),
    ("ETH-01", "ETH-01-A_single_site_review", "Single-site ethics review", "Seek review for one institution and one frozen protocol scope.", "Clear authority for a single-site implementation.", "Cannot automatically cover other sites or materially changed procedures.", "Independent review decision matching protocol SHA", False),
    ("ETH-01", "ETH-01-B_multisite_reliance", "Multisite reliance", "Use authorized reliance or coordinated review across sites.", "Can align one protocol across multiple institutions.", "Reliance scope and local responsibilities can be ambiguous.", "Reliance agreements and site-specific authorization", False),
    ("ETH-01", "ETH-01-C_authority_exempt_determination", "Authority exemption determination", "Use an exemption/not-human determination only if issued by the proper authority.", "May fit low-risk work when independently determined.", "Investigator or repository self-determination is not acceptable.", "Written authoritative determination matching scope", False),
    ("ETH-02", "ETH-02-A_documented_consent", "Documented consent", "Use a reviewed documented consent process before participation.", "Strong traceability of participant authorization.", "Documentation can add burden and identifying data that need separation.", "Approved consent form, process, storage, and withdrawal rule", False),
    ("ETH-02", "ETH-02-B_econsent", "Electronic consent", "Use an accessible reviewed electronic consent workflow.", "Supports remote participation and version control.", "Identity, accessibility, and platform retention create extra risks.", "Approved e-consent version, accessibility QA, and system controls", False),
    ("ETH-02", "ETH-02-C_authority_waiver_or_alteration", "Authority waiver or alteration", "Use waiver/alteration only when formally authorized.", "May reduce burden for qualifying minimal-risk designs.", "Cannot be inferred from convenience and may narrow permissible scope.", "Written waiver/alteration and approved information process", False),
    ("GOV-01", "GOV-01-A_local_encrypted_minimal", "Local encrypted minimal data", "Keep minimal coded data locally with frozen retention and deletion.", "Reduces transfer and data volume.", "Local security, backup, and deletion controls may be inconsistent.", "Privacy review, encryption, access, backup, retention, incident plan", False),
    ("GOV-01", "GOV-01-B_institutional_managed", "Institutional managed environment", "Use an approved institutional platform and governance controls.", "Provides managed access, audit, and incident processes.", "Platform defaults may retain excess metadata or restrict reproducibility.", "Platform approval, data map, access roles, retention, export policy", False),
    ("GOV-01", "GOV-01-C_multisite_governed_or_federated", "Multisite governed or federated", "Use data agreements or federated summaries across sites.", "Can minimize raw-data transfer and support multisite scope.", "Harmonization, disclosure, and site-level bias remain difficult.", "Data agreements, common schema, disclosure control, incident roles", False),
)


def build_decision_option_catalog() -> pd.DataFrame:
    """Return 39 non-ranked, non-recommended options with qualitative impacts."""
    register = _blocking_register().set_index("DecisionId")
    impact = {
        "SCI-01": ("indirect", "high", "conditional", "conditional", "high"),
        "SCI-02": ("direct", "indirect", "high", "conditional", "high"),
        "SCI-03": ("direct", "indirect", "high", "conditional", "medium"),
        "SCI-04": ("direct", "medium", "indirect", "none", "high"),
        "SCI-05": ("direct", "high", "conditional", "high", "high"),
        "SCI-06": ("direct", "medium", "conditional", "conditional", "high"),
        "SCI-07": ("direct", "high", "high", "conditional", "high"),
        "SCI-08": ("direct", "high", "high", "medium", "high"),
        "SCI-09": ("direct", "medium", "high", "high", "high"),
        "OPS-01": ("conditional", "medium", "medium", "conditional", "medium"),
        "ETH-01": ("conditional", "none", "medium", "medium", "blocking"),
        "ETH-02": ("conditional", "none", "high", "high", "blocking"),
        "GOV-01": ("conditional", "medium", "medium", "high", "blocking"),
    }
    external_decisions = {"ETH-01", "ETH-02", "GOV-01"}
    rows = []
    counts: dict[str, int] = {}
    for decision_id, option_id, label, description, strength, risk, evidence, numeric in _OPTION_SPECS:
        counts[decision_id] = counts.get(decision_id, 0) + 1
        sample, estimand, burden, accessibility, public = impact[decision_id]
        source = register.loc[decision_id]
        rows.append({
            "SchemaVersion": WORKBENCH_SCHEMA_VERSION,
            "DecisionId": decision_id,
            "Decision": source["Decision"],
            "Category": source["Category"],
            "OwnerRole": source["OwnerRole"],
            "OptionOrderWithinDecision": counts[decision_id],
            "OptionId": option_id,
            "OptionLabel": label,
            "Description": description,
            "Strength": strength,
            "PrimaryRisk": risk,
            "RequiredEvidence": evidence,
            "RequiresNumericInput": numeric,
            "ExternalAuthorityRequired": decision_id in external_decisions,
            "SampleSizeImpact": sample,
            "EstimandImpact": estimand,
            "ParticipantBurdenImpact": burden,
            "AccessibilityImpact": accessibility,
            "PublicClaimImpact": public,
            "OptionOrderIsRank": False,
            "AutoSelectable": False,
            "RecommendationStatus": "not_recommended_by_repository",
        })
    return pd.DataFrame(rows)


def build_decision_dependency_table() -> pd.DataFrame:
    """Return the registered acyclic dependency presentation."""
    register = _blocking_register().set_index("DecisionId")
    rows = []
    for decision_id in sorted(DECISION_IDS, key=lambda value: _TOPOLOGICAL_ORDER[value]):
        stage_order, stage = _STAGES[decision_id]
        rows.append({
            "SchemaVersion": WORKBENCH_SCHEMA_VERSION,
            "TopologicalOrder": _TOPOLOGICAL_ORDER[decision_id],
            "BroadStageOrder": stage_order,
            "BroadStageLabel": stage,
            "DecisionId": decision_id,
            "Decision": register.loc[decision_id, "Decision"],
            "PrerequisiteDecisionIds": "|".join(PREREQUISITES[decision_id]),
            "PrerequisiteCount": len(PREREQUISITES[decision_id]),
            "StageIsPermissionToSkipDependencies": False,
            "SelectionResolved": False,
        })
    return pd.DataFrame(rows)


def validate_dependency_table(table: pd.DataFrame) -> dict[str, object]:
    """Validate exact nodes, edges, acyclicity, and strict topological order."""
    failures = []
    required = {"DecisionId", "TopologicalOrder", "PrerequisiteDecisionIds"}
    if not isinstance(table, pd.DataFrame) or not required.issubset(table.columns):
        return {"valid": False, "failure_codes": ("invalid_dependency_schema",)}
    if set(table["DecisionId"].astype(str)) != set(DECISION_IDS) or len(table) != 13:
        failures.append("dependency_node_mismatch")
    order = dict(zip(table["DecisionId"].astype(str), table["TopologicalOrder"].astype(int)))
    for decision_id in DECISION_IDS:
        row = table.loc[table["DecisionId"].astype(str).eq(decision_id)]
        if len(row) != 1:
            continue
        actual = tuple(filter(None, str(row.iloc[0]["PrerequisiteDecisionIds"]).split("|")))
        if actual != PREREQUISITES[decision_id]:
            failures.append("dependency_edge_mismatch")
        if any(order.get(parent, math.inf) >= order.get(decision_id, -math.inf) for parent in actual):
            failures.append("non_topological_dependency")
    if len(set(order.values())) != 13:
        failures.append("duplicate_topological_order")
    unique = tuple(dict.fromkeys(failures))
    return {"valid": not unique, "failure_codes": unique}


def build_decision_selection_template() -> pd.DataFrame:
    """Return one unresolved blank row for each blocking decision."""
    register = _blocking_register().set_index("DecisionId")
    rows = []
    for decision_id in DECISION_IDS:
        source = register.loc[decision_id]
        rows.append({
            "SchemaVersion": WORKBENCH_SCHEMA_VERSION,
            "DecisionId": decision_id,
            "Decision": source["Decision"],
            "OwnerRole": source["OwnerRole"],
            "DecisionStatus": "unresolved",
            "SelectedOptionId": "",
            "NumericValue": "",
            "NumericUnit": "",
            "AssumptionSetReference": "",
            "RationaleReference": "",
            "OwnerAttestationReference": "",
            "ExternalEvidenceReference": "",
            "EvidenceSource": "",
            "DecisionRecordedBeforeOutcomes": False,
            "HumanOutcomesInspectedBeforeDecision": False,
            "RepositoryGeneratedSelectionAllowed": False,
        })
    return pd.DataFrame(rows, columns=SELECTION_COLUMNS)


def validate_decision_selections(
    selections: pd.DataFrame,
    catalog: pd.DataFrame | None = None,
) -> dict[str, object]:
    """Validate prospective selections while never granting recruitment readiness."""
    catalog = build_decision_option_catalog() if catalog is None else catalog
    failures = []
    if not isinstance(selections, pd.DataFrame) or tuple(selections.columns) != SELECTION_COLUMNS:
        return {"valid": False, "failure_codes": ("invalid_selection_schema",), "SelectionsComplete": False, "RecruitmentReady": False}
    if tuple(selections["DecisionId"].astype(str)) != DECISION_IDS:
        failures.append("decision_coverage_or_order_mismatch")
    if selections["DecisionId"].astype(str).duplicated().any():
        failures.append("duplicate_decision_selection")
    status = dict(zip(selections["DecisionId"].astype(str), selections["DecisionStatus"].astype(str)))
    allowed_status = {"unresolved", "resolved_prospectively"}
    if not selections["DecisionStatus"].astype(str).isin(allowed_status).all():
        failures.append("unknown_decision_status")
    selected_count = 0
    for _, row in selections.iterrows():
        decision_id = str(row["DecisionId"])
        resolved = str(row["DecisionStatus"]) == "resolved_prospectively"
        option_id = str(row["SelectedOptionId"]).strip()
        if not resolved:
            fields = ("SelectedOptionId", "NumericValue", "NumericUnit", "AssumptionSetReference", "RationaleReference", "OwnerAttestationReference", "ExternalEvidenceReference", "EvidenceSource")
            if any(str(row[field]).strip() for field in fields):
                failures.append("unresolved_row_contains_selection")
            continue
        selected_count += 1
        option = catalog.loc[catalog["OptionId"].astype(str).eq(option_id)]
        if len(option) != 1:
            failures.append("unknown_selected_option")
            continue
        option_row = option.iloc[0]
        if str(option_row["DecisionId"]) != decision_id:
            failures.append("cross_decision_option")
        if any(status.get(parent) != "resolved_prospectively" for parent in PREREQUISITES.get(decision_id, ())):
            failures.append("unresolved_prerequisite")
        required_fields = ("AssumptionSetReference", "RationaleReference", "OwnerAttestationReference")
        if any(not str(row[field]).strip() for field in required_fields):
            failures.append("missing_resolution_evidence")
        if not _truth(row["DecisionRecordedBeforeOutcomes"]) or _truth(row["HumanOutcomesInspectedBeforeDecision"]):
            failures.append("decision_not_prospectively_locked")
        if _truth(row["RepositoryGeneratedSelectionAllowed"]):
            failures.append("repository_selection_forbidden")
        if _truth(option_row["RequiresNumericInput"]):
            try:
                numeric = float(row["NumericValue"])
            except (TypeError, ValueError):
                numeric = math.nan
            if not math.isfinite(numeric) or numeric <= 0 or not str(row["NumericUnit"]).strip() or not str(row["AssumptionSetReference"]).strip():
                failures.append("missing_or_invalid_numeric_input")
        if _truth(option_row["ExternalAuthorityRequired"]):
            if not str(row["ExternalEvidenceReference"]).strip():
                failures.append("missing_external_authority_evidence")
            if not str(row["EvidenceSource"]).strip() or str(row["EvidenceSource"]).startswith("repository_generated"):
                failures.append("invalid_external_evidence_source")
    unique = tuple(dict.fromkeys(failures))
    complete = selected_count == 13 and not unique
    return {
        "valid": not unique,
        "failure_codes": unique,
        "OptionsSelected": selected_count,
        "SelectionsComplete": complete,
        "SubstantiveEvidenceVerified": False,
        "RecruitmentReady": False,
    }


def build_workbench_summary(selection_status: Mapping[str, object]) -> pd.DataFrame:
    """Return ordered private first-read cards."""
    cards = (
        (1, "overall", "blocked", "No option is selected; recruitment remains prohibited."),
        (2, "options", "catalog_ready", "39 non-ranked options expose strengths, risks, evidence, and impacts."),
        (3, "dependencies", "guarded", "13 decisions follow a strict acyclic prerequisite order."),
        (4, "sample_size", "blocked", "Numeric n is disabled until seven scientific prerequisites resolve."),
        (5, "external_authority", "blocked", "Ethics, consent, and governance require genuine external evidence."),
        (6, "evidence", "unverified", "Worksheet references are not substantively verified by this tool."),
        (7, "public_release", "withheld", "No public route, writable form, or evidence button is enabled."),
    )
    return pd.DataFrame([{
        "SchemaVersion": WORKBENCH_SCHEMA_VERSION,
        "CardOrder": order, "CardId": card_id, "Status": status, "FirstRead": read,
        "OptionsSelected": int(selection_status["OptionsSelected"]),
        "SelectionsComplete": bool(selection_status["SelectionsComplete"]),
        "RecruitmentReady": False, "PublicSurfaceEnabled": False,
    } for order, card_id, status, read in cards])


def compute_workbench_content_identity(catalog: pd.DataFrame, dependencies: pd.DataFrame, selections: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    for label, frame in (("catalog", catalog), ("dependencies", dependencies), ("selection_template", selections)):
        label_bytes, payload = label.encode(), _csv_bytes(frame)
        digest.update(len(label_bytes).to_bytes(8, "big")); digest.update(label_bytes)
        digest.update(len(payload).to_bytes(8, "big")); digest.update(payload)
    return digest.hexdigest()


def render_decision_workbench_html(catalog: pd.DataFrame, dependencies: pd.DataFrame, content_sha256: str) -> str:
    """Render a self-contained read-only bilingual workbench."""
    sections = []
    for _, dependency in dependencies.sort_values("TopologicalOrder").iterrows():
        decision_id = str(dependency["DecisionId"])
        options = catalog.loc[catalog["DecisionId"].eq(decision_id)].sort_values("OptionOrderWithinDecision")
        option_html = []
        for _, option in options.iterrows():
            option_html.append(
                f"<article class='option'><h4>{escape(str(option['OptionId']))}: {escape(str(option['OptionLabel']))}</h4>"
                f"<p>{escape(str(option['Description']))}</p><dl>"
                f"<dt>Strength</dt><dd>{escape(str(option['Strength']))}</dd>"
                f"<dt>Primary risk</dt><dd>{escape(str(option['PrimaryRisk']))}</dd>"
                f"<dt>Required evidence</dt><dd>{escape(str(option['RequiredEvidence']))}</dd>"
                f"<dt>Impacts</dt><dd>sample={escape(str(option['SampleSizeImpact']))}; estimand={escape(str(option['EstimandImpact']))}; burden={escape(str(option['ParticipantBurdenImpact']))}; accessibility={escape(str(option['AccessibilityImpact']))}; public claim={escape(str(option['PublicClaimImpact']))}</dd>"
                "</dl><p class='warning'>Not recommended, ranked, defaulted, or auto-selectable by the repository.</p></article>"
            )
        prerequisites = str(dependency["PrerequisiteDecisionIds"]) or "none"
        sections.append(
            f"<details><summary>{escape(decision_id)} — {escape(str(dependency['Decision']))}</summary>"
            f"<p><b>Broad stage:</b> {escape(str(dependency['BroadStageLabel']))}; <b>prerequisites:</b> {escape(prerequisites)}</p>"
            + "".join(option_html) + "</details>"
        )
    return """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>BLOCKED — CMLE protocol decision workbench</title><style>
body{font-family:system-ui,sans-serif;max-width:1100px;margin:auto;padding:1.2rem;color:#172033;background:#f7f8fb}header{background:#8b1e2d;color:white;padding:1rem 1.3rem;border-radius:.5rem}.notice{border-left:.4rem solid #8b1e2d;background:white;padding:1rem;margin:1rem 0}details{background:white;margin:.7rem 0;padding:.8rem;border:1px solid #ccd2dc;border-radius:.4rem}summary{font-weight:700;cursor:pointer}.option{border-top:1px solid #dde2ea;padding:.5rem 0}dt{font-weight:700}dd{margin-bottom:.4rem}.warning{color:#7a2530;font-weight:650}code{overflow-wrap:anywhere}@media(max-width:650px){body{padding:.6rem}}
</style></head><body><header><h1>BLOCKED — recruitment prohibited / 募集不可</h1><p>Private read-only decision aid. No option is selected or recommended.</p></header>
<section class='notice'><p>This offline workbench is not ethics approval, scientific approval, recruitment authority, or a selected sample size.</p><p>この非公開資料は選択肢を比較するためのもので、倫理承認、科学的承認、募集許可、標本数決定ではありません。</p><p>WorkbenchContentSHA256: <code>""" + escape(content_sha256) + "</code></p></section><main>" + "".join(sections) + "</main></body></html>"


def validate_workbench_html(html: str, catalog: pd.DataFrame, content_sha256: str) -> dict[str, object]:
    failures = []
    for required in (
        "<title>BLOCKED — CMLE protocol decision workbench</title>",
        "<h1>BLOCKED — recruitment prohibited / 募集不可</h1>",
        "No option is selected or recommended",
        content_sha256,
    ):
        if required not in html: failures.append("missing_blocked_or_identity_content")
    lower = html.lower()
    if any(token in lower for token in ("<form", "<input", "<script", "http://", "https://")):
        failures.append("writable_or_external_html_content")
    if any(str(option) not in html for option in catalog["OptionId"]):
        failures.append("missing_option_in_html")
    unique = tuple(dict.fromkeys(failures))
    return {"valid": not unique, "failure_codes": unique}


def _zip(payloads: Mapping[str, bytes]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name in sorted(payloads):
            info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0)); info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16; info.create_system = 3
            archive.writestr(info, payloads[name], compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return buffer.getvalue()


def build_private_decision_workbench_bundle() -> dict[str, object]:
    catalog = build_decision_option_catalog(); dependencies = build_decision_dependency_table()
    selections = build_decision_selection_template(); status = validate_decision_selections(selections, catalog)
    summary = build_workbench_summary(status); content_sha = compute_workbench_content_identity(catalog, dependencies, selections)
    html = render_decision_workbench_html(catalog, dependencies, content_sha)
    start = f"""# BLOCKED — private confirmatory decision workbench

OptionsSelected=0. SelectionsComplete=false. RecruitmentReady=false. PublicSurfaceEnabled=false.

This seven-member offline bundle presents alternatives and dependencies. It does not recommend or select an option and is not ethics approval or recruitment authority.

募集不可。選択肢は未選択です。リポジトリは推奨、倫理承認、募集許可、標本数決定を行いません。

WorkbenchContentSHA256: `{content_sha}`
"""
    payloads = {
        "START_HERE.md": start.encode(),
        "decision_dependency_table.csv": _csv_bytes(dependencies),
        "decision_option_catalog.csv": _csv_bytes(catalog),
        "decision_selection_template.csv": _csv_bytes(selections),
        "decision_workbench_summary.csv": _csv_bytes(summary),
        "protocol_decision_workbench.html": html.encode(),
    }
    manifest = pd.DataFrame([{"Member": name, "ByteCount": len(payload), "SHA256": _sha(payload), "HumanRows": 0, "ContainsSelection": False, "PublicArtifact": False} for name, payload in sorted(payloads.items())])
    payloads["BUNDLE_MANIFEST.csv"] = _csv_bytes(manifest)
    return {"catalog": catalog, "dependencies": dependencies, "selections": selections, "selection_status": status, "summary": summary, "content_sha256": content_sha, "html": html, "manifest": manifest, "bundle_bytes": _zip(payloads)}


def validate_private_decision_workbench_bundle(bundle_bytes: bytes) -> dict[str, object]:
    failures = []
    try:
        with zipfile.ZipFile(io.BytesIO(bundle_bytes), "r") as archive:
            infos = archive.infolist(); names = tuple(info.filename for info in infos)
            if names != tuple(sorted(WORKBENCH_BUNDLE_MEMBERS)): failures.append("member_order_or_membership_mismatch")
            for info in infos:
                if info.date_time != (1980,1,1,0,0,0): failures.append("non_deterministic_timestamp")
                if (info.external_attr >> 16) & 0o777 != 0o644: failures.append("non_deterministic_permissions")
            manifest = pd.read_csv(io.BytesIO(archive.read("BUNDLE_MANIFEST.csv")), dtype=str)
            expected = set(names) - {"BUNDLE_MANIFEST.csv"}
            if set(manifest["Member"]) != expected: failures.append("manifest_member_mismatch")
            for _, row in manifest.iterrows():
                payload = archive.read(str(row["Member"]))
                if len(payload) != int(row["ByteCount"]): failures.append("manifest_byte_count_mismatch")
                if _sha(payload) != str(row["SHA256"]): failures.append("manifest_hash_mismatch")
                if int(row["HumanRows"]) != 0 or _truth(row["ContainsSelection"]): failures.append("human_rows_or_selection_present")
                if _truth(row["PublicArtifact"]): failures.append("public_artifact_present")
            start = archive.read("START_HERE.md").decode(); html = archive.read("protocol_decision_workbench.html").decode()
            if any(value not in start for value in ("BLOCKED", "募集不可", "OptionsSelected=0", "RecruitmentReady=false")): failures.append("first_read_not_blocked")
            catalog = pd.read_csv(io.BytesIO(archive.read("decision_option_catalog.csv")))
            identity_marker = "WorkbenchContentSHA256: `"; content_sha = start.split(identity_marker,1)[1].split("`",1)[0] if identity_marker in start else ""
            if not validate_workbench_html(html, catalog, content_sha)["valid"]: failures.append("invalid_workbench_html")
    except (zipfile.BadZipFile, KeyError, ValueError, UnicodeDecodeError, IndexError) as error:
        failures.append(f"invalid_zip:{type(error).__name__}")
    unique = tuple(dict.fromkeys(failures))
    return {"valid": not unique, "failure_codes": unique, "bundle_sha256": _sha(bundle_bytes)}


def decision_workbench_human_gate_status(selection_status: Mapping[str, object]) -> pd.DataFrame:
    return pd.DataFrame([{
        "SchemaVersion": WORKBENCH_SCHEMA_VERSION, "HumanParticipants": 0,
        "HumanStudyStatus": "not_started_no_human_data", "OptionsSelected": int(selection_status["OptionsSelected"]),
        "SelectionsComplete": bool(selection_status["SelectionsComplete"]), "SubstantiveEvidenceVerified": False,
        "RecruitmentReady": False, "MinimumValidPerCellRegistered": False,
        "PlannedRecruitmentNSelected": False, "ConfirmatoryResultAvailable": False,
        "PublicSurfaceEnabled": False,
    }])
