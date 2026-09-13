"""Static guardrails for the standalone, Streamlit-light Python core."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from pathlib import Path
import re

import pytest


CORE_PACKAGE = Path(__file__).resolve().parents[1] / "mfrm_app"
REPOSITORY_ROOT = CORE_PACKAGE.parent
STREAMLIT_ENTRYPOINT = REPOSITORY_ROOT / "streamlit_app.py"
LEGACY_COMPATIBILITY_PACKAGE = CORE_PACKAGE / "legacy_compat"

FORBIDDEN_IMPORT_ROOTS = {
    "subprocess",
    "streamlit",
    "rpy2",
    "mfrmr",
    "tam",
    "sirt",
    "mirt",
    "facets",
    "conquest",
    "cmdstanpy",
    "pystan",
    "stan",
    "juliacall",
    "julia",
    "requests",
    "httpx",
    "httpcore",
    "aiohttp",
    "urllib3",
    "socket",
    "ftplib",
    "websockets",
    "grpc",
    "ctypes",
    "cffi",
}
FORBIDDEN_IMPORT_NAMES = {"urllib.request", "http.client"}
FORBIDDEN_CALLS = {
    "subprocess.run",
    "subprocess.Popen",
    "subprocess.call",
    "subprocess.check_call",
    "subprocess.check_output",
    "os.system",
    "os.popen",
    "asyncio.create_subprocess_exec",
    "asyncio.create_subprocess_shell",
    "importlib.import_module",
    "importlib.util.spec_from_file_location",
    "importlib.util.module_from_spec",
    "runpy.run_path",
    "runpy.run_module",
    "importlib.machinery.SourceFileLoader",
    "__import__",
    "builtins.__import__",
    "exec",
    "eval",
    "compile",
    "os.posix_spawn",
    "os.posix_spawnp",
    "os.startfile",
    "pty.spawn",
    "ctypes.CDLL",
    "ctypes.PyDLL",
}
FORBIDDEN_STRING_PATTERNS = {
    "Rscript command": re.compile(r"(?:^|\n)\s*(?:#![^\n]*\bRscript\b|Rscript\s+\S+)"),
    "R package namespace call": re.compile(r"\b(?:TAM|sirt|mirt|mfrmr)::[A-Za-z0-9_.]+\s*\("),
    "R package loader": re.compile(
        r"\b(?:library|requireNamespace)\s*\(\s*['\"]?(?:TAM|sirt|mirt|mfrmr)\b"
    ),
    "CmdStan runner": re.compile(
        r"(?:^|\n)\s*(?:from|import)\s+cmdstanpy\b|\bcmdstanr::|\blibrary\s*\(\s*['\"]?cmdstanr\b"
    ),
}
FORBIDDEN_RUNTIME_DISTRIBUTIONS = {
    "rpy2",
    "mfrmr",
    "cmdstanpy",
    "pystan",
    "stan",
    "juliacall",
    "julia",
    "requests",
    "httpx",
    "httpcore",
    "aiohttp",
    "urllib3",
    "websockets",
    "grpcio",
    "tam",
    "sirt",
    "mirt",
    "facets",
    "conquest",
    "ctypes",
    "cffi",
}
FORBIDDEN_CORE_ASSET_SUFFIXES = {
    ".r",
    ".rmd",
    ".qmd",
    ".jl",
    ".stan",
    ".sh",
    ".bash",
    ".zsh",
    ".fish",
    ".ksh",
    ".csh",
    ".tcsh",
    ".command",
    ".bat",
    ".cmd",
    ".ps1",
}


@dataclass(frozen=True)
class BoundaryViolation:
    path: str
    line: int
    kind: str
    detail: str


def _qualified_name(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _qualified_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return ""


def _import_aliases(tree: ast.AST) -> dict[str, str]:
    aliases: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                bound = alias.asname or alias.name.split(".", 1)[0]
                aliases[bound] = alias.name if alias.asname else bound
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            for alias in node.names:
                bound = alias.asname or alias.name
                aliases[bound] = ".".join(part for part in (module, alias.name) if part)
    return aliases


def _resolve_alias(name: str, aliases: dict[str, str]) -> str:
    head, separator, tail = name.partition(".")
    resolved_head = aliases.get(head, head)
    return f"{resolved_head}.{tail}" if separator else resolved_head


def _forbidden_import(name: str) -> bool:
    root = name.split(".", 1)[0].lower()
    lowered = name.lower()
    return (
        root in FORBIDDEN_IMPORT_ROOTS
        or lowered in FORBIDDEN_IMPORT_NAMES
        or "legacy_compat" in lowered.split(".")
    )


def scan_core_source(source: str, *, path: str = "<memory>") -> list[BoundaryViolation]:
    """Return external-runtime and UI-dependency violations in core source."""
    tree = ast.parse(source, filename=path)
    aliases = _import_aliases(tree)
    violations: list[BoundaryViolation] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if _forbidden_import(alias.name):
                    violations.append(
                        BoundaryViolation(path, node.lineno, "import", alias.name)
                    )
        elif isinstance(node, ast.ImportFrom):
            module = node.module or ""
            imported_names = [
                ".".join(part for part in (module, alias.name) if part)
                for alias in node.names
            ]
            if node.level:
                forbidden = [
                    name for name in imported_names
                    if "legacy_compat" in name.lower().split(".")
                ]
                module_forbidden = "legacy_compat" in module.lower().split(".")
            else:
                forbidden = [name for name in imported_names if _forbidden_import(name)]
                module_forbidden = _forbidden_import(module)
            if module_forbidden or forbidden:
                violations.append(
                    BoundaryViolation(
                        path,
                        node.lineno,
                        "import",
                        forbidden[0] if forbidden else ("." * node.level + module),
                    )
                )
        elif isinstance(node, ast.Call):
            name = _resolve_alias(_qualified_name(node.func), aliases)
            if name in FORBIDDEN_CALLS or name.startswith("os.spawn") or name.startswith("os.exec"):
                violations.append(
                    BoundaryViolation(path, node.lineno, "runtime call", name)
                )
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            for label, pattern in FORBIDDEN_STRING_PATTERNS.items():
                if pattern.search(node.value):
                    violations.append(
                        BoundaryViolation(path, node.lineno, "external handoff marker", label)
                    )
    return violations


def _normalize_distribution_name(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name.strip().lower())


def _name_from_direct_reference(reference: str) -> str:
    egg = re.search(r"[#&]egg=([A-Za-z0-9_.-]+)", reference)
    if egg:
        return _normalize_distribution_name(egg.group(1))
    filename = Path(reference.split("?", 1)[0].rstrip("/")).name
    wheel = re.match(r"([A-Za-z0-9_.-]+?)-\d[^/]*\.whl$", filename, flags=re.IGNORECASE)
    if wheel:
        return _normalize_distribution_name(wheel.group(1))
    return ""


def scan_requirement_files(path: Path) -> tuple[set[str], list[str]]:
    """Resolve requirement includes and fail closed on unidentified references."""
    names: set[str] = set()
    unresolved: list[str] = []
    seen: set[Path] = set()

    def scan(current: Path) -> None:
        resolved = current.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        if not resolved.is_file():
            unresolved.append(f"missing requirement file: {current}")
            return
        for raw_line in resolved.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            include = re.match(r"^(?:-r|--requirement|-c|--constraint)\s+(.+)$", line)
            if include:
                scan((resolved.parent / include.group(1).strip()).resolve())
                continue
            editable = re.match(r"^(?:-e|--editable)\s+(.+)$", line)
            if editable:
                name = _name_from_direct_reference(editable.group(1).strip())
                if name:
                    names.add(name)
                else:
                    unresolved.append(line)
                continue
            requirement = line.split(";", 1)[0].strip()
            named_reference = re.match(r"^([A-Za-z0-9_.-]+)(?:\[[^]]+\])?\s*@\s*.+$", requirement)
            if named_reference:
                names.add(_normalize_distribution_name(named_reference.group(1)))
                continue
            if requirement.startswith(("git+", "hg+", "svn+", "bzr+", "http://", "https://", ".", "/")):
                name = _name_from_direct_reference(requirement)
                if name:
                    names.add(name)
                else:
                    unresolved.append(line)
                continue
            name_match = re.match(r"^([A-Za-z0-9_.-]+)(?:\[[^]]+\])?(?:\s|[<>=!~]|$)", requirement)
            if name_match:
                names.add(_normalize_distribution_name(name_match.group(1)))
            else:
                unresolved.append(line)

    scan(path)
    return names, unresolved


def _forbidden_asset_reason(path: Path, text: str = "") -> str:
    if path.suffix.lower() in FORBIDDEN_CORE_ASSET_SUFFIXES:
        return f"forbidden executable/handoff suffix {path.suffix}"
    if not path.suffix and path.exists() and path.stat().st_mode & 0o111:
        return "extensionless executable"
    for label, pattern in FORBIDDEN_STRING_PATTERNS.items():
        if pattern.search(text):
            return f"external handoff marker: {label}"
    return ""


def _core_python_files() -> list[Path]:
    return [
        path
        for path in sorted(CORE_PACKAGE.rglob("*.py"))
        if LEGACY_COMPATIBILITY_PACKAGE not in path.parents
    ]


def test_core_package_has_no_external_engine_or_streamlit_runtime_surface():
    violations: list[BoundaryViolation] = []
    for path in _core_python_files():
        violations.extend(
            scan_core_source(
                path.read_text(encoding="utf-8"),
                path=str(path.relative_to(CORE_PACKAGE.parent)),
            )
        )

    assert not violations, "\n".join(
        f"{item.path}:{item.line}: {item.kind}: {item.detail}"
        for item in violations
    )


def test_streamlit_entrypoint_has_no_repository_external_python_engine_loader():
    """Keep the supported estimator on the repository-local Python path."""
    tree = ast.parse(
        STREAMLIT_ENTRYPOINT.read_text(encoding="utf-8"),
        filename=str(STREAMLIT_ENTRYPOINT),
    )
    aliases = _import_aliases(tree)
    forbidden_calls: list[tuple[int, str]] = []
    forbidden_names: list[tuple[int, str]] = []
    legacy_symbols = {
        "_load_shared_simulation_engine",
        "SHARED_SIM_ENGINE",
        "SHARED_SIM_ENGINE_LOAD_ERROR",
        "USE_EMBEDDED_ENGINE_ONLY",
    }
    dynamic_loader_calls = {
        "importlib.util.spec_from_file_location",
        "importlib.util.module_from_spec",
    }
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = _resolve_alias(_qualified_name(node.func), aliases)
            if name in dynamic_loader_calls or name.endswith(".exec_module"):
                forbidden_calls.append((node.lineno, name))
        elif isinstance(node, ast.Name) and node.id in legacy_symbols:
            forbidden_names.append((node.lineno, node.id))

    assert forbidden_calls == []
    assert forbidden_names == []


def test_core_runtime_dependencies_exclude_engines_and_remote_workers():
    requirement_names, unresolved = scan_requirement_files(
        REPOSITORY_ROOT / "requirements.txt"
    )
    forbidden = sorted(
        name
        for name in requirement_names
        if name in FORBIDDEN_RUNTIME_DISTRIBUTIONS or name.startswith("rpy2-")
    )
    assert unresolved == []
    assert forbidden == []


def test_core_package_contains_no_external_runner_assets():
    forbidden_assets: list[str] = []
    for path in CORE_PACKAGE.rglob("*"):
        if not path.is_file() or LEGACY_COMPATIBILITY_PACKAGE in path.parents or path.suffix == ".py":
            continue
        try:
            text = path.read_text(encoding="utf-8") if path.stat().st_size <= 2_000_000 else ""
        except UnicodeDecodeError:
            text = ""
        reason = _forbidden_asset_reason(path, text)
        if reason:
            forbidden_assets.append(f"{path.relative_to(REPOSITORY_ROOT)}: {reason}")
    assert forbidden_assets == []


@pytest.mark.parametrize(
    "source, expected_detail",
    [
        ("import rpy2.robjects\n", "rpy2.robjects"),
        ("import subprocess\nsubprocess.run(['Rscript', 'fit.R'])\n", "subprocess"),
        ("import os as operating_system\noperating_system.system('solver')\n", "os.system"),
        ("from os import system as launch\nlaunch('solver')\n", "os.system"),
        ("import importlib\nimportlib.import_module('solver')\n", "importlib.import_module"),
        (
            "import importlib.util\nimportlib.util.spec_from_file_location('engine', '/tmp/engine.py')\n",
            "importlib.util.spec_from_file_location",
        ),
        ("SCRIPT = 'fit <- TAM::tam.mml.mfr(resp)'\n", "R package namespace call"),
        ("from mfrm_app.legacy_compat import tam_handoff\n", "mfrm_app.legacy_compat"),
        ("from .legacy_compat import tam_handoff\n", "legacy_compat"),
        ("from mfrm_app import legacy_compat\n", "mfrm_app.legacy_compat"),
        ("engine = __import__('rpy' + '2')\n", "__import__"),
        ("import requests\nrequests.post('https://solver.example')\n", "requests"),
        ("import urllib3\nurllib3.PoolManager().request('POST', 'https://solver.example')\n", "urllib3"),
        ("from urllib import request\nrequest.urlopen('https://solver.example')\n", "urllib.request"),
        ("import os\nos.posix_spawn('Rscript', ['Rscript', 'fit.R'], {})\n", "os.posix_spawn"),
        ("import ctypes\nctypes.CDLL('libexternal_engine.so')\n", "ctypes"),
        ("exec('import rpy2')\n", "exec"),
        ("import streamlit as st\n", "streamlit"),
    ],
)
def test_boundary_scanner_detects_forbidden_surfaces(source: str, expected_detail: str):
    violations = scan_core_source(source)
    assert any(expected_detail in item.detail for item in violations), violations


def test_boundary_scanner_allows_native_statistics_facets_and_report_archives():
    source = """
import scipy.optimize
import zipfile

FACET_LABEL = "FACETS-style display convention"
METHOD_REFERENCE = "sirt::pcm.fit"

def fit_facet(values):
    return scipy.optimize.minimize(lambda x: sum((x - value) ** 2 for value in values), 0.0)
"""
    assert scan_core_source(source) == []


def test_boundary_scanner_allows_relative_native_facets_module():
    assert scan_core_source("from .facets import normalize\n") == []


def test_requirement_scanner_resolves_includes_and_direct_references(tmp_path: Path):
    included = tmp_path / "external.txt"
    included.write_text(
        "-e git+https://example.invalid/rpy2.git#egg=rpy2\n"
        "https://example.invalid/cmdstanpy-1.2.3-py3-none-any.whl\n",
        encoding="utf-8",
    )
    root = tmp_path / "requirements.txt"
    root.write_text("-r external.txt\nnumpy>=1.24\n", encoding="utf-8")

    names, unresolved = scan_requirement_files(root)

    assert unresolved == []
    assert {"numpy", "rpy2", "cmdstanpy"}.issubset(names)


def test_asset_guard_rejects_notebook_handoff_and_external_code_text():
    assert _forbidden_asset_reason(Path("runner.Rmd"), "")
    assert _forbidden_asset_reason(Path("runner.txt"), "fit <- TAM::tam.mml.mfr(resp)")


def test_standalone_make_and_workflows_do_not_expose_parity_exporter():
    makefile = (REPOSITORY_ROOT / "Makefile").read_text(encoding="utf-8")
    lines = makefile.splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith("verify:"))
    verify_block: list[str] = []
    for line in lines[start:]:
        if verify_block and line and not line[0].isspace() and ":" in line:
            break
        verify_block.append(line)
    workflows = "\n".join(
        path.read_text(encoding="utf-8")
        for suffix in ("*.yml", "*.yaml")
        for path in sorted((REPOSITORY_ROOT / ".github" / "workflows").glob(suffix))
    )

    assert "--export-parity-fixture" not in makefile
    assert not re.search(r"(?m)^(?:parity|compatibility-parity)\s*:", makefile)
    assert "compatibility-parity" not in "\n".join(verify_block)
    assert "--export-parity-fixture" not in workflows
    assert "compatibility-parity" not in workflows
    assert not re.search(r"\b(?:make|gmake)\s+parity\b", workflows)
    assert 'NATIVE_PYTEST_ARGS ?= -m "not legacy_compat and not retained_evidence"' in makefile
    assert "make apptest" in workflows
    assert "make ux-contracts" in workflows
    assert "Verify tracked checkout stayed clean" in workflows
    assert "git diff --exit-code" in workflows
    assert "git status --porcelain=v1 --untracked-files=all" in workflows
    assert not re.search(r"python\s+-m\s+pytest\s+tests", workflows)


def test_standalone_cli_and_self_test_registry_do_not_call_legacy_generators():
    tree = ast.parse(
        STREAMLIT_ENTRYPOINT.read_text(encoding="utf-8"),
        filename=str(STREAMLIT_ENTRYPOINT),
    )
    functions = {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }

    cli_flags_node = next(
        node
        for node in tree.body
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "CLI_CHECK_FLAGS" for target in node.targets)
    )
    assert "--export-parity-fixture" not in ast.literal_eval(cli_flags_node.value)

    main_guard = next(
        node
        for node in tree.body
        if isinstance(node, ast.If)
        and any(isinstance(part, ast.Name) and part.id == "__name__" for part in ast.walk(node.test))
    )
    dispatch_names = {
        part.id for part in ast.walk(main_guard) if isinstance(part, ast.Name)
    }
    assert "export_reference_parity_fixture" not in dispatch_names

    legacy_self_tests = {
        "_self_test_cross_package_validation_plan",
        "_self_test_posterior_viewer_loaders",
        "_self_test_advanced_model_generators",
        "_self_test_posterior_load_netcdf",
        "_self_test_posterior_load_cmdstan_csvs",
        "_self_test_cross_engine_bundle",
    }
    run_self_tests = functions["run_self_tests"]
    registered_names = {
        part.id for part in ast.walk(run_self_tests) if isinstance(part, ast.Name)
    }
    assert registered_names.isdisjoint(legacy_self_tests)

    legacy_generators = {
        "_generate_repro_r_script",
        "external_simulation_reference_inventory",
        "external_simulation_template_inventory",
        "external_simulation_template_scripts",
        "external_validation_artifact_checklist",
        "external_validation_report_template",
        "reproducibility_script_export_matrix",
        "bayesian_mfrm_stan_refinement_plan",
        "bayesian_stan_runner_templates",
        "generate_advanced_model_stan_code",
        "posterior_viewer_example_package_assets",
        "build_cross_engine_validation_bundle",
    }
    registered_calls: set[str] = set()
    for function_name in registered_names:
        function = functions.get(function_name)
        if function is None:
            continue
        registered_calls.update(
            part.func.id
            for part in ast.walk(function)
            if isinstance(part, ast.Call) and isinstance(part.func, ast.Name)
        )
    assert registered_calls.isdisjoint(legacy_generators)


def test_public_streamlit_routes_do_not_expose_legacy_compatibility_tools():
    source = STREAMLIT_ENTRYPOINT.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(STREAMLIT_ENTRYPOINT))
    functions = {
        node.name: ast.get_source_segment(source, node) or ""
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }

    forbidden_by_route = {
        "main": {"Posterior Viewer (upload)", "render_posterior_viewer_mode"},
        "run_facets_mode": {
            "generate_advanced_model_stan_code",
            "facets_mode_advanced_generate",
        },
        "show_report_section": {"_render_stan_code", '"Stan Code"'},
        "show_help_section": {
            "guided_stan_posterior_reproducibility_help_table",
            "external_simulation_reference_inventory",
            "mfrmr_020_migration_coverage_table",
            'selected_help_label == "Model Capability"',
            'selected_help_label == "Public Beta"',
        },
        "show_classical_dif_section": {"build_dif_validation_bundle", "mfrm_difR_crosscheck_bundle.zip"},
        "show_tutorial": {"install.packages", "TAM", "sirt", "mirt", "eRm"},
        "generate_method_appendix_text": {"functional parity target", "external-validation claims"},
    }
    for function_name, forbidden_tokens in forbidden_by_route.items():
        route_source = functions[function_name]
        assert forbidden_tokens.isdisjoint(
            token for token in forbidden_tokens if token in route_source
        ), function_name

    assert "standalone_release_limitations_table()" in functions["render_app_scope_badges"]
    assert "python_yardstick_reproducibility_assets()" in functions["_draw_yardstick"]
    assert "yardstick_reproducibility_scripts(" not in functions["_draw_yardstick"]
    assert "python_rating_scale_recode_assets(" in functions["show_categories_section"]
    assert "rating_scale_recode_script_assets(" not in functions["show_categories_section"]

    assert "python_only=True" in functions["python_yardstick_reproducibility_assets"]
    assert "python_only=True" in functions["python_rating_scale_recode_assets"]


def test_native_claim_and_visual_builders_exclude_legacy_product_rows():
    source = STREAMLIT_ENTRYPOINT.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(STREAMLIT_ENTRYPOINT))
    function_sources = {
        node.name: ast.get_source_segment(source, node) or ""
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
    }

    claim_source = function_sources["build_manuscript_claim_guide"]
    assert '"ManuscriptArea": "External package comparison"' not in claim_source
    assert '"ManuscriptArea": "Bayesian / Uto-family extension"' not in claim_source

    visual_source = function_sources["visual_method_evidence_table"]
    assert "Posterior Viewer" not in visual_source
    assert "TAM/mirt-style" not in visual_source
