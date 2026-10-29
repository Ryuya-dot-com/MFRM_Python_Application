"""Check main's frozen desktop and archive the exact tested application folder."""
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
    manifest = json.loads((build / 'desktop_build.json').read_text(encoding='utf-8'))
    root = Path(__file__).resolve().parents[1]
    for relative, expected in manifest['source_sha256'].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Source changed since build: {relative}')
    if json.loads((bundle / '_internal/mfrm_app/desktop_build.json').read_text(encoding='utf-8')) != manifest:
        raise ValueError('Bundled manifest differs from the build record')
    for relative, expected in manifest['third_party_notices']['files'].items():
        if hashlib.sha256((bundle / '_internal/third_party_notices' / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Bundled notice differs: {relative}')
    with tempfile.TemporaryDirectory(prefix='mfrm-frozen-check-') as working:
        environment = dict(os.environ)
        environment.pop('PYTHONPATH', None)
        environment.pop('PYTHONHOME', None)
        def run(*args):
            with (output / 'process.log').open('ab') as log:
                subprocess.run([str(exe), *map(str, args)], cwd=working, env=environment,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
        smoke = output / 'analysis'
        run('--smoke-test', smoke)
        record = json.loads((smoke / 'smoke.json').read_text(encoding='utf-8'))
        if (record['worker_result'] != 'complete' or record['retained_rows'] != 47
                or not record['missing_row_excluded'] or record['report_svg_count'] != 3):
            raise ValueError('Frozen missing-score workflow did not pass')
        project = smoke / '作文.mfrm'
        before = hashlib.sha256(project.read_bytes()).hexdigest()
        destination = output / 'reopened'
        run('--check-project', project, '--output-dir', destination)
        reopened = json.loads((destination / 'reopen.json').read_text(encoding='utf-8'))
        if reopened['estimation_calls'] != 0:
            raise ValueError('Restore re-estimated a model')
        if (destination / 'report.html').read_bytes() != (smoke / 'report.html').read_bytes():
            raise ValueError('Restored report differs')
        if hashlib.sha256(project.read_bytes()).hexdigest() != before:
            raise ValueError('Original project changed during restore')
    shutil.copyfile(root / 'docs/desktop_windows_installation.md', bundle / 'START_HERE.txt')
    archive = Path(shutil.make_archive(str(build / 'MFRM-Desktop-Windows-x64'), 'zip', bundle.parent, bundle.name))
    with zipfile.ZipFile(archive) as zipped:
        if zipped.testzip() is not None:
            raise ValueError('ZIP integrity check failed')
        for path in bundle.rglob('*'):
            if path.is_file() and zipped.read(path.relative_to(bundle.parent).as_posix()) != path.read_bytes():
                raise ValueError(f'ZIP differs from the tested bundle: {path}')
    sha256 = hashlib.sha256(archive.read_bytes()).hexdigest()
    (build / 'SHA256SUMS.txt').write_text(f'{sha256}  {archive.name}\n', encoding='utf-8')
    result = dict(source_base='d78fa706a7e2dc5f01574b3101009f97c4daeb69', workflow=record,
                  reopened=reopened, archive=archive.name, bytes=archive.stat().st_size,
                  sha256=sha256, exe_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                  windows_11_gui_accepted=False,
                  limits=['Windows runner process checks; Windows 11 GUI/download acceptance remains open',
                          'Current main desktop: RSM/JMLE; unsigned validation candidate'])
    (output / 'checks.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--build-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(check_build(args.build_dir, args.output), indent=2))
