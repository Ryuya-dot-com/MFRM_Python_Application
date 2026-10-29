"""Check the built Windows executable, then package the same verified directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile


def check_build(build, output):
    build, output = build.resolve(), output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    bundle = build / 'dist/MFRM-Desktop'
    exe = bundle / 'MFRM-Desktop.exe'
    if not exe.is_file():
        raise FileNotFoundError(exe)
    manifest = json.loads((build / 'desktop_build.json').read_text(encoding='utf-8'))
    root = Path(__file__).resolve().parents[1]
    for relative, expected in manifest['source_sha256'].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Source changed since build: {relative}')
    smoke = output / 'analysis'
    # Do not let imports or assets resolve from the checkout / current directory.
    with tempfile.TemporaryDirectory(prefix='mfrm-frozen-check-') as working:
        environment = dict(os.environ)
        environment.pop('PYTHONPATH', None)
        environment.pop('PYTHONHOME', None)
        def run(*args):
            with (output / 'process.log').open('ab') as log:
                subprocess.run([str(exe), *map(str, args)], cwd=working, env=environment,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        run('--smoke-test', smoke)
        record = json.loads((smoke / 'smoke.json').read_text(encoding='utf-8'))
        expected = {('RSM', 'JMLE'), ('PCM', 'JMLE'), ('RSM', 'MML'), ('PCM', 'MML'), ('G/D', 'ANOVA')}
        if {(r['model'], r['method']) for r in record['extended']} != expected:
            raise ValueError('Frozen checks did not cover all shipped workflows')
        reopened = []
        for project in sorted(smoke.glob('*/analysis.mfrm')):
            before = hashlib.sha256(project.read_bytes()).hexdigest()
            destination = output / 'reopened' / project.parent.name
            run('--check-project', project, '--output-dir', destination)
            result = json.loads((destination / 'reopen.json').read_text(encoding='utf-8'))
            if result['estimation_calls'] != 0:
                raise ValueError('Restore re-estimated a model')
            if (destination / 'report.html').read_bytes() != (project.parent / 'report.html').read_bytes():
                raise ValueError(f'Restored report differs: {project}')
            if hashlib.sha256(project.read_bytes()).hexdigest() != before:
                raise ValueError(f'Original project was changed: {project}')
            reopened.append(dict(workflow=project.parent.name, **result))
        if len(reopened) != 5:
            raise ValueError('Expected five independently restored projects')
    # ZIP contains the complete onedir app; never distribute the exe on its own.
    shutil.copyfile(root / 'docs/desktop_windows_installation.md', bundle / 'START_HERE.txt')
    archive = Path(shutil.make_archive(str(build / 'MFRM-Desktop-Windows-x64'), 'zip', bundle.parent, bundle.name))
    with zipfile.ZipFile(archive) as zipped:
        if zipped.testzip() is not None:
            raise ValueError('ZIP integrity check failed')
        for path in bundle.rglob('*'):
            if path.is_file() and zipped.read(path.relative_to(bundle.parent).as_posix()) != path.read_bytes():
                raise ValueError(f'ZIP differs from tested app: {path}')
    sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
    (build / 'SHA256SUMS.txt').write_text(f'{sha256}  {archive.name}\n', encoding='utf-8')
    result = dict(workflows=record['extended'], reopened=reopened, archive=archive.name,
                  bytes=archive.stat().st_size, sha256=sha256, windows_11_gui_accepted=False,
                  limits=['Build-host process checks; no clean Windows 11 GUI or browser-download acceptance',
                          'Unsigned validation candidate; no public-release approval implied'])
    (output / 'checks.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(check_build(args.build_dir, args.output), indent=2))
