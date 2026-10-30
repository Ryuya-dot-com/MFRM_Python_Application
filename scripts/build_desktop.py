"""Build a local prototype; does not sign, notarize, publish, or install it."""
import argparse
import ast
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import sysconfig


def collect_notices(output):
    """Retain the installed distributions' notices, including bundled-library text."""
    from packaging.requirements import Requirement
    directory = output / 'third_party_notices'
    directory.mkdir(parents=True, exist_ok=True)
    packages, files = {}, {}

    def copy(source, relative):
        payload = source.read_bytes()
        target = directory / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        files[str(relative)] = hashlib.sha256(payload).hexdigest()

    pending = ['numpy', 'pandas', 'scipy', 'matplotlib', 'pyarrow', 'pyinstaller']
    while pending:
        distribution = importlib.metadata.distribution(pending.pop())
        name = distribution.metadata['Name'].lower().replace('_', '-')
        if name in packages:
            continue
        packages[name] = distribution.version
        notices = [f for f in distribution.files or [] if
                   any(word in str(f).lower() for word in ('license', 'notice', 'copying', 'copyright'))
                   and f.suffix.lower() not in ('.py', '.pyc', '.so', '.dll', '.dylib')]
        if not notices:
            raise RuntimeError(f'No installed license/notice files found for {name}')
        for path in notices:
            # Keep relative subdirectories to preserve identically named notices.
            if '..' in path.parts or path.is_absolute():
                raise RuntimeError(f'Unexpected notice path: {path}')
            copy(Path(distribution.locate_file(path)), Path(name) / path)
        for value in distribution.requires or []:
            requirement = Requirement(value)
            if requirement.marker is None or requirement.marker.evaluate({'extra': ''}):
                pending.append(requirement.name)

    base = Path(sys.base_prefix)
    python_license_dir = base if sys.platform == 'win32' else Path(sysconfig.get_path('stdlib'))
    copy(python_license_dir / 'LICENSE.txt', Path('Python/LICENSE.txt'))
    if sys.platform == 'darwin':
        copy(base / 'Resources/English.lproj/Documentation/_sources/license.rst.txt',
             Path('Python/incorporated-software-notices.rst'))
        for name in ('Tcl', 'Tk'):
            copy(base / f'Frameworks/{name}.framework/Versions/8.6/Resources/license.terms',
                 Path(name) / 'license.terms')
    elif sys.platform == 'win32':
        # CPython's Windows installer omits the separate Tcl/Tk license files.
        for name in ('tcl', 'tk'):
            copy(Path(__file__).parent / 'runtime_notices' / f'{name}-8.6.15.terms',
                 Path(name) / 'license.terms')
    else:
        for name in ('tcl', 'tk'):
            copy(base / f'tcl/{name}8.6/license.terms', Path(name) / 'license.terms')
    index = ('Third-party notices for this build\n\n'
             'Package notices include their embedded-library acknowledgements.\n'
             'PyInstaller and its build dependencies are listed for build/bootloader attribution;\n'
             'this list does not mean every build tool is part of the runtime.\n'
             'System fonts are selected from the host OS and are not redistributed.\n\n'
             f'Python: {sys.version}\n' + '\n'.join(f'{name} {version}' for name, version in sorted(packages.items()))
             + '\n\nFiles (SHA-256):\n' + '\n'.join(f'{digest}  {name}' for name, digest in sorted(files.items())) + '\n')
    (directory / 'INDEX.txt').write_text(index, encoding='utf-8')
    return dict(packages=packages, files=files)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if (sys.platform, platform.machine().lower()) not in {("darwin", "arm64"), ("win32", "amd64"), ("win32", "x86_64")}:
        parser.error("Initial targets are macOS arm64 and Windows x64; build on the target OS/CPU")
    root = Path(__file__).resolve().parents[1]
    app_version = next(ast.literal_eval(node.value)
                       for node in ast.parse((root / "mfrm_app" / "engine.py").read_text(encoding="utf-8")).body
                       if isinstance(node, ast.Assign)
                       and any(isinstance(target, ast.Name) and target.id == "APP_VERSION" for target in node.targets))
    if json.loads((root / "CITATION.cff").read_text(encoding="utf-8"))["version"] != app_version:
        parser.error("APP_VERSION and CITATION.cff must agree before packaging")
    version_info_path = root / "scripts" / "windows_version_info.txt"
    if f"'{app_version}'" not in version_info_path.read_text(encoding="utf-8"):
        parser.error("APP_VERSION and scripts/windows_version_info.txt must agree before packaging")
    bundle_version = app_version.split("-", 1)[0]
    if not re.fullmatch(r"\d+\.\d+\.\d+", bundle_version):
        parser.error("The bundle version must have three numeric components")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    notices = collect_notices(output)
    source_paths = [root / "desktop_app.py", root / "examples" / "desktop_writing.csv",
                    root / "LICENSE", root / "LICENSE_NOTICE.md", root / "CITATION.cff", root / "locales" / "en.json",
                    root / "scripts" / "desktop.spec", version_info_path,
                    Path(__file__).resolve(), root / "requirements-desktop-build.txt",
                    *sorted((root / "scripts" / "icons").glob("*")),
                    *sorted((root / "scripts" / "runtime_notices").glob("*")),
                    *sorted((root / "mfrm_app").glob("*.py"))]
    record = {
        "scope": "local prototype; developer signing/notarization and platform acceptance not completed",
        "app_version": app_version, "bundle_version": bundle_version,
        "platform": platform.platform(), "machine": platform.machine(), "python": sys.version,
        "source_sha256": {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_paths},
        "packages": {name: importlib.metadata.version(name) for name in
                     ("pyinstaller", "numpy", "pandas", "scipy", "matplotlib", "pyarrow")},
        "third_party_notices": notices,
    }
    manifest = output / "desktop_build.json"
    manifest.write_text(json.dumps(record, indent=2), encoding="utf-8")
    command = [sys.executable, "-m", "PyInstaller", "--noconfirm",
               "--distpath", str(output / "dist"), "--workpath", str(output / "work"),
               str(root / "scripts" / "desktop.spec")]
    environment = dict(os.environ, PYINSTALLER_CONFIG_DIR=str(output / "cache"),
                       MFRM_DESKTOP_BUILD_MANIFEST=str(manifest))
    environment.setdefault("MPLCONFIGDIR", str(output / "matplotlib-cache"))
    subprocess.run(command, cwd=root, env=environment, check=True)


if __name__ == "__main__":
    main()
