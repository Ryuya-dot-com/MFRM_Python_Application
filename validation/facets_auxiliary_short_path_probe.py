"""Diagnose FACETS 4.5 auxiliary-report behavior at a short work path.

This is an operational probe only.  It does not read measurement estimates or
fit statistics, and its output is not eligible for a scientific endpoint.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from validation.operating_characteristics_facets import invoke_facets, sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_spec", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--facets", type=Path, default=Path(r"C:\Facets\Facets.exe"))
    parser.add_argument("--timeout", type=float, default=60.0)
    args = parser.parse_args()

    source_spec = args.source_spec.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite probe output: {output_dir}")
    output_dir.mkdir(parents=True)

    spec_path = output_dir / "a.txt"
    report_path = output_dir / "r.txt"
    primary_score_base = output_dir / "s.txt"
    auxiliary_score_base = output_dir / "u.txt"
    source = source_spec.read_text(encoding="utf-8", errors="strict")
    source, replacements = re.subn(
        r"(?im)^\s*Scorefile\s*=.*$",
        lambda _match: f"Scorefile={primary_score_base}",
        source,
    )
    if replacements != 1:
        raise ValueError(f"Expected one Scorefile specification, found {replacements}")
    spec_path.write_text(source, encoding="utf-8", newline="\n")

    started = time.perf_counter()
    completed = invoke_facets(
        args.facets.resolve(),
        spec_path,
        report_path,
        timeout_seconds=args.timeout,
        extra_specs=("Umean=0,1,2", f"Scorefile={auxiliary_score_base}"),
        batch_mode="NO",
    )
    elapsed = time.perf_counter() - started
    (output_dir / "stdout.txt").write_text(completed.stdout or "", encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr or "", encoding="utf-8")

    expected_scores = [output_dir / f"u.{number}.txt" for number in range(1, 5)]
    result = {
        "schema_version": "facets-auxiliary-short-path-probe-v1",
        "operational_only": True,
        "scientific_endpoint_read": False,
        "source_spec": str(source_spec),
        "source_spec_sha256": sha256_file(source_spec),
        "probe_script_sha256": sha256_file(Path(__file__).resolve()),
        "facets_executable": str(args.facets.resolve()),
        "facets_executable_sha256": sha256_file(args.facets.resolve()),
        "returncode": completed.returncode,
        "elapsed_seconds": elapsed,
        "path_lengths": {
            "spec": len(str(spec_path)),
            "report": len(str(report_path)),
            "score_base": len(str(auxiliary_score_base)),
            "longest_expected_score": max(len(str(path)) for path in expected_scores),
        },
        "report": {
            "exists": report_path.is_file(),
            "size_bytes": report_path.stat().st_size if report_path.is_file() else 0,
            "sha256": sha256_file(report_path) if report_path.is_file() else None,
        },
        "scores": {
            path.name: {
                "exists": path.is_file(),
                "size_bytes": path.stat().st_size if path.is_file() else 0,
                "sha256": sha256_file(path) if path.is_file() else None,
            }
            for path in expected_scores
        },
        "pass": bool(
            completed.returncode == 0
            and report_path.is_file()
            and report_path.stat().st_size > 0
            and all(path.is_file() and path.stat().st_size > 0 for path in expected_scores)
        ),
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, ensure_ascii=False))
    if not result["pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
