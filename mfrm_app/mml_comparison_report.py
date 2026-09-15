"""Four separate numerical comparison answers; never scientific qualification.

Callers supply measurements and their provenance. This presentation layer does
not authenticate external runs or infer evidence from a native convergence flag.
"""
from __future__ import annotations

from html import escape
import math
from numbers import Real

import pandas as pd

from .evidence import canonical_json
from .exports import to_html_report


QUESTIONS = {
    "calibration": ("校正値は近いか", "Are calibrations close?"),
    "stationarity": ("返却点で推定は停留しているか", "Is the returned point stationary?"),
    "integration": ("積分精度は目安を満たすか", "Does integration meet the stated targets?"),
    "native_scores": ("標準得点の安定性を確認できたか", "Has native-score stability been assessed?"),
}
STATUS = {
    "within_targets": ("記載した目安内", "Within stated targets"),
    "exceeds_targets": ("目安未達", "Targets not met"),
    "not_assessed": ("判定保留", "Not assessed"),
}


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be nonempty text")
    return value


def _number(value, name, *, positive=False):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number or None")
    if value < 0 or (positive and value == 0):
        raise ValueError(f"{name} must be {'positive' if positive else 'nonnegative'}")
    return float(value)


def build_comparison_report(*, title: str, checks: dict, details: dict) -> dict:
    """Use strict upper targets on nonnegative error metrics, independently by question.

    Missing values, absent targets and missing prerequisites never produce a pass.
    A measured target failure remains visible even when other evidence is missing.
    Targets are supplied explicitly; this function invents no equivalence margin.
    """
    _text(title, "title")
    if not isinstance(checks, dict) or set(checks) != set(QUESTIONS):
        raise ValueError("Exactly four comparison questions are required")
    if not isinstance(details, dict) or not details:
        raise ValueError("Evaluation definitions and provenance details are required")
    answer = {}
    for key in QUESTIONS:
        check = checks[key]
        if not isinstance(check, dict) or set(check) != {"scope", "target_basis", "missing_evidence", "metrics"}:
            raise ValueError(f"Malformed comparison check: {key}")
        scope = _text(check["scope"], "scope")
        target_basis = _text(check["target_basis"], "target_basis")
        missing = check["missing_evidence"]
        if not isinstance(missing, list):
            raise ValueError("missing_evidence must be a list")
        missing = [_text(v, "missing_evidence") for v in missing]
        if not isinstance(check["metrics"], list):
            raise ValueError("metrics must be a list")
        metrics = []
        for metric in check["metrics"]:
            if not isinstance(metric, dict) or set(metric) != {"name", "value", "upper_target", "unit"}:
                raise ValueError("Malformed comparison metric")
            value = _number(metric["value"], "value")
            limit = _number(metric["upper_target"], "upper_target", positive=True)
            met = None if value is None or limit is None else value < limit
            metrics.append(dict(name=_text(metric["name"], "metric name"), value=value,
                upper_target=limit, unit=_text(metric["unit"], "unit"), target_met=met))
        names = [m["name"] for m in metrics]
        if len(set(names)) != len(names):
            raise ValueError("Metric names must be unique within each question")
        status = ("exceeds_targets" if any(m["target_met"] is False for m in metrics) else
                  "not_assessed" if missing or not metrics or any(m["target_met"] is None for m in metrics) else
                  "within_targets")
        answer[key] = dict(scope=scope, target_basis=target_basis, missing_evidence=missing,
                           metrics=metrics, status=status)
    result = dict(schema_version="mfrm_numerical_comparison_v1", title=title,
        scientific_inference_ready=False, qualification_eligible=False, checks=answer, details=details)
    # Defensive JSON copy; rejects nonfinite values anywhere in the provenance.
    import json
    return json.loads(canonical_json(result))


def _checked(report):
    """Rebuild decisions when rendering/exporting; never trust a stored pass label."""
    if report.get("schema_version") != "mfrm_numerical_comparison_v1":
        raise ValueError("Unsupported comparison report version")
    if report.get("scientific_inference_ready") is not False or report.get("qualification_eligible") is not False:
        raise ValueError("Numerical comparison reports cannot grant inference qualification")
    checks = {key:{name:value for name,value in check.items() if name != "status"}
              for key,check in report["checks"].items()}
    for check in checks.values():
        check["metrics"] = [{k:v for k,v in m.items() if k != "target_met"} for m in check["metrics"]]
    return build_comparison_report(title=report["title"], checks=checks, details=report["details"])


def comparison_frames(report, *, language="ja"):
    report = _checked(report)
    if language not in ("ja", "en"):
        raise ValueError("language must be ja or en")
    lang = int(language == "en")
    summary, metrics = [], []
    for key in QUESTIONS:
        check = report["checks"][key]
        summary.append(dict(Question=QUESTIONS[key][lang], Status=STATUS[check["status"]][lang],
            Scope=check["scope"], MissingEvidence="; ".join(check["missing_evidence"])))
        for metric in check["metrics"]:
            metrics.append(dict(Question=QUESTIONS[key][lang], Metric=metric["name"],
                Value=metric["value"], UpperTarget=metric["upper_target"], Operator="<",
                Unit=metric["unit"], TargetMet=metric["target_met"], TargetBasis=check["target_basis"]))
    return {"comparison_summary": pd.DataFrame(summary), "comparison_metrics": pd.DataFrame(metrics,
        columns=["Question", "Metric", "Value", "UpperTarget", "Operator", "Unit", "TargetMet", "TargetBasis"])}


def comparison_html(report, *, language="ja") -> bytes:
    report = _checked(report)
    frames = comparison_frames(report, language=language)
    ja = language == "ja"
    summary = frames["comparison_summary"].rename(columns={
        "Question":"確認すること", "Status":"今回の結果", "Scope":"対象と意味", "MissingEvidence":"未確認のこと"
    } if ja else {})
    heading = "比較の要点" if ja else "Comparison answers"
    html = to_html_report({heading:summary}, title=escape(report["title"])).decode()
    boundary = ("数値比較の記録です。未評価の項目を含むことがあります。推論資格・標準誤差・信頼区間・被覆率の認定ではありません。"
                if ja else "Numerical comparison record; checks may be unassessed. No qualification of inference, standard errors, intervals or coverage.")
    label = "数値と目安を見る" if ja else "Measurements and targets"
    provenance = "評価点・計算条件・保存記録を見る" if ja else "Evaluation points, settings and provenance"
    detail = (f'<details><summary>{label}</summary>{frames["comparison_metrics"].to_html(index=False, na_rep="—")}</details>'
              f'<details><summary>{provenance}</summary><pre>{escape(canonical_json(report["details"]))}</pre></details>')
    html = html.replace("<html>", f'<html lang="{language}">').replace("</head>",
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        '<style>body{max-width:1100px;margin:2rem auto;padding:0 1rem;line-height:1.5}'
        'details{margin:1.2rem 0;padding:.8rem;border:1px solid #bbb;border-radius:.4rem;overflow:auto}'
        'summary{cursor:pointer;font-weight:600}pre{white-space:pre-wrap;overflow-wrap:anywhere}'
        'summary:focus-visible{outline:3px solid #17617d}th{white-space:normal}</style></head>')
    html = html.replace("<body>", f'<body><p>{boundary}</p>').replace("</body>", detail+"</body>")
    return html.encode("utf-8")
