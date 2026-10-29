"""Notice inclusion and refusal to distribute incomplete frozen checks."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import zipfile

import pytest

import desktop_app
from scripts import build_desktop, check_windows_desktop


@pytest.mark.parametrize('system', ['win32', 'darwin'])
def test_runtime_notices_use_the_installer_layout(tmp_path, monkeypatch, system):
    base, stdlib = tmp_path / 'python', tmp_path / 'stdlib'
    files = [base / 'LICENSE.txt'] if system == 'win32' else [stdlib / 'LICENSE.txt',
        base / 'Resources/English.lproj/Documentation/_sources/license.rst.txt']
    if system == 'darwin':
        for name in ('Tcl', 'Tk'):
            files.append(base / f'Frameworks/{name}.framework/Versions/8.6/Resources/license.terms')
    for path in files:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f'License from {path.name}')
    package_license = tmp_path / 'package-LICENSE.txt'
    package_license.write_text('Package license')
    monkeypatch.setattr(build_desktop.importlib.metadata, 'distribution', lambda name:
        SimpleNamespace(metadata={'Name': name}, version='1', files=[Path('LICENSE.txt')],
                        requires=[], locate_file=lambda path: package_license))
    monkeypatch.setattr(build_desktop.sys, 'base_prefix', str(base))
    monkeypatch.setattr(build_desktop.sys, 'platform', system)
    monkeypatch.setattr(build_desktop.sysconfig, 'get_path', lambda name: str(stdlib))
    result = build_desktop.collect_notices(tmp_path / 'build')
    assert {'Python/LICENSE.txt', 'Tcl/license.terms' if system == 'darwin' else 'tcl/license.terms',
            'Tk/license.terms' if system == 'darwin' else 'tk/license.terms'} <= result['files'].keys()
    if system == 'win32':
        for name in ('tcl', 'tk'):
            source = Path(build_desktop.__file__).parent / 'runtime_notices' / f'{name}-8.6.15.terms'
            assert (tmp_path / 'build/third_party_notices' / name / 'license.terms').read_bytes() == source.read_bytes()

@pytest.mark.parametrize('mode', ['--smoke-test', '--check-project'])
def test_windowed_validation_failure_exits_with_a_durable_traceback(tmp_path, monkeypatch, mode):
    def broken(*args):
        raise ValueError('deliberate validation failure')
    monkeypatch.setattr(desktop_app, 'smoke_check', broken)
    monkeypatch.setattr(desktop_app, 'smoke_reopen', broken)
    args = [mode, str(tmp_path)] if mode == '--smoke-test' else [mode, 'sample.mfrm', '--output-dir', str(tmp_path)]
    monkeypatch.setattr(sys, 'argv', ['MFRM-Desktop.exe', *args])
    with monkeypatch.context() as patch:
        patch.setattr(sys, 'stderr', None)
        with pytest.raises(SystemExit) as error:
            desktop_app.main()
    assert error.value.code == 1
    assert 'deliberate validation failure' in (tmp_path / 'failure.txt').read_text()


@pytest.mark.parametrize('broken', [None, 'worker', 'report', 'original'])
def test_packaging_requires_successful_restore_and_preserves_the_bundle(tmp_path, monkeypatch, broken):
    monkeypatch.setattr(check_windows_desktop, '__file__', str(tmp_path / 'scripts/check_windows_desktop.py'))
    (tmp_path / 'docs').mkdir()
    (tmp_path / 'docs/desktop_windows_installation.md').write_text('Install notes')
    (tmp_path / 'source.py').write_text('source')
    build = tmp_path / 'build'
    bundle = build / 'dist/MFRM-Desktop'
    (bundle / '_internal/mfrm_app').mkdir(parents=True)
    (bundle / 'MFRM-Desktop.exe').write_bytes(b'candidate')
    (bundle / '_internal/runtime.dll').write_bytes(b'runtime')
    manifest = {'source_sha256': {'source.py': hashlib.sha256(b'source').hexdigest()},
                'third_party_notices': {'files': {}}}
    for path in (build / 'desktop_build.json', bundle / '_internal/mfrm_app/desktop_build.json'):
        path.write_text(json.dumps(manifest))
    original = None

    def process(command, **kwargs):
        nonlocal original
        assert Path(kwargs['cwd']) != tmp_path
        if broken == 'worker':
            raise subprocess.CalledProcessError(1, command)
        if command[1] == '--smoke-test':
            out = Path(command[2]); out.mkdir()
            (out / 'smoke.json').write_text(json.dumps(dict(worker_result='complete', retained_rows=47,
                                                         missing_row_excluded=True, report_svg_count=3)))
            (out / 'report.html').write_text('three figures')
            original = out / '作文.mfrm'
            original.write_bytes(b'project')
        else:
            out = Path(command[4]); out.mkdir()
            (out / 'reopen.json').write_text(json.dumps({'analysis_id': 'saved', 'estimation_calls': 0}))
            (out / 'report.html').write_text('changed' if broken == 'report' else 'three figures')
            if broken == 'original':
                original.write_bytes(b'changed')

    monkeypatch.setattr(check_windows_desktop.subprocess, 'run', process)
    if broken:
        with pytest.raises((ValueError, subprocess.CalledProcessError)):
            check_windows_desktop.check_build(build, tmp_path / 'checks')
        assert not (build / 'MFRM-Desktop-Windows-x64.zip').exists()
    else:
        result = check_windows_desktop.check_build(build, tmp_path / 'checks')
        assert result['reopened']['estimation_calls'] == 0
        with zipfile.ZipFile(build / result['archive']) as zipped:
            assert zipped.read('MFRM-Desktop/_internal/runtime.dll') == b'runtime'
            assert zipped.read('MFRM-Desktop/START_HERE.txt') == b'Install notes'
        assert original.read_bytes() == b'project'


