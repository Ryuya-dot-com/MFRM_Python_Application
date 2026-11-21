"""Public Python interface to the shared MFRM application calculations.

Use :func:`check_design` and :func:`check_rater_contrast` for optimization-free
preflight checks, :func:`analyze`
for RSM/PCM or observed-score G/D studies, :func:`summary`, :func:`tables` and
:func:`export_csv` for qualified results, and :func:`save_project` /
:func:`load_project` for portable .mfrm archives. Raw retained results do not
override output limits.
"""
from contextlib import contextmanager
import hashlib
from importlib.resources import files
import os
import tempfile
import warnings
from pathlib import Path

import pandas as pd

from mfrm_app import basic_analysis, desktop_analysis, exports, project
from mfrm_app.desktop_analysis import check_design, run_analysis as analyze
from mfrm_app.project import ProjectError, load_project, save_project
from mfrm_app.rater_contrast import check_rater_contrast

__version__ = "0.2.2"

__all__ = [
    "__version__", "analyze", "check_design", "check_rater_contrast", "summary", "tables", "export_csv",
    "load_example", "save_project", "load_project", "ProjectError", "plot", "report_html",
]

_TABLE_ALIASES = {
    "estimates": "measures",
    "fit": "measures",
    "category": "categories",
    "step": "steps",
    "thresholds": "steps",
    "linking": "connectivity",
    "anchors": "connectivity",
    "separation": "reliability",
    "interval": "intervals",
    "joint_interaction": "joint",
}

_VALID_TABLE_KINDS = frozenset({
    "design", "measures", "categories", "steps", "connectivity", "reliability",
    "agreement", "agreement_summary", "pairwise_agreement", "rescaled",
    "bias", "intervals", "joint", "gd", "gd_observed", "gd_selected",
    "gd_variance", "gd_d_study", "gd_metadata", "pca",
    *_TABLE_ALIASES,
})


def _canonical_table_kind(title: str) -> str:
    if title == "Analysis Design and Retained Ratings":
        return "design"
    if title == "Measures and Fit":
        return "measures"
    if title == "Observed Category Use":
        return "categories"
    if title == "Step Thresholds":
        return "steps"
    if title == "Connectivity and Scale Linking":
        return "connectivity"
    if title == "Separation and Reliability by Facet":
        return "reliability"
    if title == "Inter-rater Agreement and Krippendorff's Alpha":
        return "agreement_summary"
    if title.startswith("Pairwise Agreement:"):
        return "pairwise_agreement"
    if title == "Rescaled Measures (Umean / Uscale)":
        return "rescaled"
    if title.startswith("Bias:"):
        return "bias"
    if title.startswith("Specified intervals:"):
        return "intervals"
    if title.startswith("Joint interaction:") or title == "Joint model: main effects and steps":
        return "joint"
    if title == "G/D: observed_coefficients":
        return "gd_observed"
    if title == "G/D: selected_d_study":
        return "gd_selected"
    if title == "G/D: variance_components":
        return "gd_variance"
    if title == "G/D: d_study":
        return "gd_d_study"
    if title == "G/D: metadata":
        return "gd_metadata"
    if title.startswith("Residual PCA:"):
        return "pca"
    return "other"


def _available_plot_kinds(analysis) -> list[str]:
    kinds: list[str] = []
    if "result" in analysis:
        kinds.extend([
            "wright", "pathway", "category", "facet", "rater_severity", "step",
            "bubble", "expected", "information", "fair_average", "unexpected", "residual_qq",
        ])
    if analysis.get("bias_analyses"):
        kinds.extend(["bias", "bias_ranked"])
    if analysis.get("interval_analyses"):
        kinds.append("interval")
    if "gd_analysis" in analysis:
        kinds.append("gd")
    if "joint_analysis" in analysis:
        kinds.append("joint")
    if (analysis.get("diagnostics") or {}).get("pca_enabled"):
        kinds.extend(["pca", "pca_loadings"])
    return kinds


def tables(analysis, kind=None):
    """Return (title, DataFrame, interpretation note) tuples without refitting.

    Pass ``kind`` (e.g. ``"measures"``, ``"steps"``, ``"reliability"``,
    ``"agreement"``, ``"categories"``, ``"design"``, ``"bias"``, ``"intervals"``,
    ``"joint"``, ``"gd"``, ``"pca"``) or an iterable of kinds to filter specific tables.
    """
    all_tables = desktop_analysis.current_tables(analysis)
    if kind is None:
        return all_tables
    requested = [kind] if isinstance(kind, str) else list(kind)
    if not requested or any(not isinstance(k, str) or not k.strip() for k in requested):
        raise ValueError("kind must be a nonblank table kind string or iterable of strings")
    title_lookup = {t.lower(): t for t, _, _ in all_tables}
    target_kinds: set[str] = set()
    target_titles: set[str] = set()
    for raw_k in requested:
        norm_k = raw_k.strip().lower()
        mapped = _TABLE_ALIASES.get(norm_k, norm_k)
        if mapped == "agreement":
            target_kinds.update({"agreement_summary", "pairwise_agreement"})
        elif mapped == "gd":
            target_kinds.update({"gd_observed", "gd_selected", "gd_variance", "gd_d_study", "gd_metadata"})
        elif mapped in _VALID_TABLE_KINDS:
            target_kinds.add(mapped)
        elif norm_k in title_lookup:
            target_titles.add(title_lookup[norm_k])
        else:
            raise ValueError(f"Unknown table kind: {raw_k!r}")
    return [
        (t, f, n)
        for t, f, n in all_tables
        if _canonical_table_kind(t) in target_kinds or t in target_titles
    ]


def summary(analysis) -> pd.DataFrame:
    """Return a concise one-screen DataFrame of design, qualification, and available outputs without refitting."""
    if reason := project.input_refit_reason(analysis):
        project.validate_analysis(analysis)
        return pd.DataFrame([("AnalysisID", analysis["analysis_id"]), ("RunQualificationState", "refit_required"),
                             ("Reason", reason), ("SupportsFormalInference", "False")], columns=["Field", "Value"])
    all_tables = desktop_analysis.current_tables(analysis)
    avail_tables = ", ".join(dict.fromkeys(
        ("gd" if k.startswith("gd_") else ("agreement" if k in {"agreement_summary", "pairwise_agreement"} else k))
        for k in (_canonical_table_kind(t) for t, _, _ in all_tables)
    ))
    avail_plots = ", ".join(_available_plot_kinds(analysis))

    if analysis.get("kind") == "gd-analysis-v2":
        gd_entry = analysis["gd_analysis"]
        gd_cfg = gd_entry["settings"]
        gd_res = gd_entry["result"]
        avail = bool(gd_res.get("available"))
        obs_c = gd_res.get("observed_coefficients") or {}
        rows = [
            ("Workflow", "gd"),
            ("AnalysisID", str(analysis.get("analysis_id", ""))),
            ("Available", str(avail)),
            ("ReasonCode", str(gd_res.get("reason_code") or "none")),
            ("ObjectFacet", str(gd_cfg.get("object_facet", ""))),
            ("RandomFacets", ", ".join(map(str, gd_cfg.get("random_facets") or ()))),
            ("Observations", str((gd_res.get("design") or {}).get("n_observations", 0))),
            ("ObservedG", f"{float(obs_c['G']):.3f}" if avail and pd.notna(obs_c.get("G")) else ""),
            ("ObservedPhi", f"{float(obs_c['Phi']):.3f}" if avail and pd.notna(obs_c.get("Phi")) else ""),
            ("SupportsFormalInference", "False"),
            ("AvailableTables", avail_tables),
            ("AvailablePlots", avail_plots),
        ]
        return pd.DataFrame(rows, columns=["Field", "Value"])

    result = analysis["result"]
    cfg = result.get("config") or {}
    prep = result.get("prep") or {}
    snap = basic_analysis.qualification_snapshot(result)
    fnames = list(cfg.get("facet_names") or ())
    flevels = cfg.get("facet_levels") or {}
    fsigns = cfg.get("facet_signs") or {}
    pos_facets = [f for f in fnames if int(fsigns.get(f, -1)) == 1]
    unused = prep.get("unused_score_categories") or []
    rows = [
        ("Workflow", "mfrm"),
        ("AnalysisID", str(analysis.get("analysis_id", ""))),
        ("Model", str(cfg.get("model", ""))),
        ("Method", str(cfg.get("method", ""))),
        ("RetainedRatings", str(len(prep.get("data", ())))),
        ("Persons", str(cfg.get("n_person", ""))),
        ("Facets", ", ".join(f"{f} ({len(flevels.get(f, ()))})" for f in fnames)),
        ("NoncenterFacet", str(cfg.get("noncenter_facet") or "Person")),
        ("PositiveFacets", ", ".join(pos_facets) if pos_facets else "none"),
        ("ScoreRange", f"{prep.get('rating_min')}–{prep.get('rating_max')} (unused: {', '.join(map(str, unused)) if unused else 'none'})"),
        ("OptimizerTerminated", str(bool(snap.get("optimizer_terminated")))),
        ("RunQualificationState", str(snap.get("run_state", ""))),
        ("IdentifiabilityState", str(snap.get("identifiability_state") or "n/a")),
        ("FormalHoldReason", str(snap.get("calibration_hold_reason") or "none")),
        ("ScreeningBlockReason", str(snap.get("screening_block_reason") or "none")),
        ("SupportsFormalInference", "False"),
        ("AvailableTables", avail_tables),
        ("AvailablePlots", avail_plots),
    ]
    return pd.DataFrame(rows, columns=["Field", "Value"])


@contextmanager
def _manifest_lock(directory):
    # OS locks are released if a writer exits, including on Windows.
    with (directory / ".mfrmpy-export.lock").open("a+b") as lock:
        if os.name == "nt":
            import msvcrt
            lock.seek(0)
            msvcrt.locking(lock.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl
            fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def _previous_export(manifest):
    if not manifest.exists():
        return []
    previous = pd.read_csv(manifest, dtype=str, keep_default_na=False)
    required = {"AnalysisID", "TableIndex", "File", "Title", "Rows", "RunQualificationState"}
    if not required.issubset(previous):
        raise ValueError(f"{manifest} already exists and is not an mfrmpy export manifest")
    managed = []
    for record in previous.to_dict("records"):
        relative = Path(record["File"])
        index = int(record["TableIndex"])
        slug = exports.safe_export_name(record["Title"], default=f"table_{index}", max_len=64)
        if relative.name != f"{index:02d}_{slug}.csv" or relative.is_absolute():
            continue
        if len(relative.parts) > 1 and (len(relative.parts) != 2 or not relative.parts[0].startswith("mfrmpy-tables-")):
            continue
        path = manifest.parent / relative
        if path.is_symlink() or path.parent.is_symlink() or not path.is_file():
            continue
        if "SHA256" in record:
            if hashlib.sha256(path.read_bytes()).hexdigest() != record["SHA256"]:
                continue
        else:
            # 0.2.0 wrote flat tables without digests. Verify their provenance.
            try:
                table = pd.read_csv(path, dtype=str, usecols=["AnalysisID", "Table"], keep_default_na=False)
            except (ValueError, UnicodeError):
                continue
            if (len(table) != int(record["Rows"]) or not table.AnalysisID.eq(record["AnalysisID"]).all()
                    or not table.Table.eq(record["Title"]).all()):
                continue
        managed.append(path)
    return managed


def export_csv(analysis, directory) -> list[Path]:
    """Write qualified CSV tables and a provenance manifest without refitting.

    Applies the same shared qualification rules as :func:`tables`, :func:`plot`
    and :func:`report_html`, while recording both current and save-time
    qualification reasons in ``00_manifest.csv``.
    """
    qualified_tables = tables(analysis)
    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)

    analysis_id = str(analysis.get("analysis_id", ""))
    has_mfrm = "result" in analysis
    current_snap = basic_analysis.qualification_snapshot(analysis["result"]) if has_mfrm else {}
    saved_snap = analysis.get("qualification_at_save") if isinstance(analysis.get("qualification_at_save"), dict) else {}

    generation = Path(tempfile.mkdtemp(prefix="mfrmpy-tables-", dir=out_dir))
    try:
        manifest_rows: list[dict[str, object]] = []
        written_paths: list[Path] = [out_dir / "00_manifest.csv"]
        for idx, (title, frame, note) in enumerate(qualified_tables, start=1):
            slug = exports.safe_export_name(title, default=f"table_{idx}", max_len=64)
            filename = f"{idx:02d}_{slug}.csv"
            target = generation / filename
            out_frame = frame.copy()
            out_frame["AnalysisID"] = analysis_id
            out_frame["ContentFingerprint"] = analysis["content_fingerprint"]
            out_frame["Table"] = title
            out_frame["Note"] = note
            payload = out_frame.to_csv(index=False).encode("utf-8")
            project.write_bytes_atomic(target, payload)
            written_paths.append(target)
            manifest_rows.append({
                "AnalysisID": analysis_id,
                "TableIndex": idx,
                "File": target.relative_to(out_dir).as_posix(),
                "SHA256": hashlib.sha256(payload).hexdigest(),
                "ContentFingerprint": analysis["content_fingerprint"],
                "PackageVersion": __version__,
                "Title": title,
                "Rows": int(len(frame)),
                "SupportsFormalInference": False,
                "RunQualificationState": str(current_snap.get("run_state", "gd_observed_scores" if not has_mfrm else "")),
                "FormalHoldReason": str(current_snap.get("calibration_hold_reason", "")),
                "ScreeningBlockReason": str(current_snap.get("screening_block_reason", "")),
                "SavedHoldReason": str(saved_snap.get("calibration_hold_reason", "")),
                "SavedHoldMessage": str(saved_snap.get("calibration_hold_message", "")),
                "Note": note,
            })

        manifest_df = pd.DataFrame(manifest_rows)
        with _manifest_lock(out_dir):
            previous = _previous_export(written_paths[0])
            project.write_bytes_atomic(written_paths[0], manifest_df.to_csv(index=False).encode("utf-8"))
    except Exception:
        for path in generation.iterdir():
            path.unlink()
        generation.rmdir()
        raise

    for path in previous:
        try:
            path.unlink(missing_ok=True)
        except OSError as exc:
            warnings.warn(f"Previous export table retained: {exc}", stacklevel=2)
    for directory in {p.parent for p in previous if p.parent != out_dir}:
        try:
            directory.rmdir()
        except OSError:
            pass  # Retain any user files added to that directory.
    return written_paths


def load_example():
    """Return a fresh copy of the synthetic 0–5 writing ratings as a DataFrame."""
    resource = files(__package__).joinpath("data", "writing.csv")
    if not resource.is_file():
        resource = Path(__file__).resolve().parents[1] / "examples" / "desktop_writing.csv"
    with resource.open("r", encoding="utf-8") as handle:
        return pd.read_csv(handle)


def plot(analysis, kind="wright", *, level=None, index=0, facet=None,
         statistic="Infit", presentation="explained", gd_swap=False, top_n=40):
    """Return a Matplotlib Figure; install the ``figures`` extra first.

    Supported ``kind`` values include ``"wright"``, ``"pathway"`` (alias
    ``"fit_pathway"``; ``statistic="Infit"|"Outfit"|"InfitZSTD"|"OutfitZSTD"``),
    ``"category"`` (alias ``"ccc"``), ``"facet"``, ``"rater_severity"``,
    ``"step"`` (alias ``"thresholds"``), ``"bubble"``, ``"expected"`` (alias
    ``"tcc"``), ``"information"``, ``"fair_average"``, ``"unexpected"``,
    ``"residual_qq"``, ``"bias"``, ``"bias_ranked"``, ``"interval"``,
    ``"joint"``, ``"gd"``, ``"pca"`` (alias ``"residual_pca"``), and
    ``"pca_loadings"``. PCM Wright maps and category curves require ``level``.
    No estimator is rerun to draw results.
    """
    from mfrm_app.desktop_outputs import analysis_figure
    return analysis_figure(analysis, kind, level=level, index=index, facet=facet,
                           statistic=statistic, presentation=presentation,
                           gd_swap=gd_swap, top_n=top_n)


def report_html(analysis, *, pathway_facet=None, pathway_statistic="Infit",
                facet=None, statistic=None, presentation=None, language=None, gd_swap=None):
    """Return an offline HTML report; requires the ``figures`` extra.

    ``facet`` and ``statistic`` are accepted as aliases for ``pathway_facet``
    and ``pathway_statistic`` to match :func:`plot`. Reports and .mfrm archives
    can contain the original ratings and identifiers; handle them as research data.
    """
    from mfrm_app.desktop_outputs import report_html as render
    return render(analysis, pathway_facet=pathway_facet,
                  pathway_statistic=pathway_statistic, facet=facet, statistic=statistic,
                  presentation=presentation, language=language, gd_swap=gd_swap)
