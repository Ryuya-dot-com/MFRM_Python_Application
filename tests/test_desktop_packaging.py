"""Build/launcher contracts that can be checked before running the Windows job."""
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
            'Tk/license.terms' if system == 'darwin' else 'tk/license.terms'} <= {
                Path(name).as_posix() for name in result['files']}
    if system == 'win32':
        for name in ('tcl', 'tk'):
            source = Path(build_desktop.__file__).parent / 'runtime_notices' / f'{name}-8.6.15.terms'
            assert (tmp_path / 'build/third_party_notices' / name / 'license.terms').read_bytes() == source.read_bytes()


def test_project_argument_reaches_the_existing_reader_after_ui_initialization(tmp_path, monkeypatch):
    import tkinter
    events = []
    root = SimpleNamespace(after_idle=lambda callback, path: events.append((callback, path)),
                           mainloop=lambda: events.append('loop'))
    app = SimpleNamespace(open_project=lambda path: None)
    monkeypatch.setattr(tkinter, 'Tk', lambda: root)
    monkeypatch.setattr(desktop_app, 'DesktopApp', lambda window: app)
    path = tmp_path / '作文 sample.mfrm'
    monkeypatch.setattr(sys, 'argv', ['MFRM-Desktop.exe', str(path)])
    desktop_app.main()
    assert events == [(app.open_project, path), 'loop']


@pytest.mark.parametrize('closed_at', ['poll', 'recv'])
@pytest.mark.parametrize('terminal,exitcode', [(('complete', ''), 0), (('complete', ''), 1),
                                             (('failed', 'invalid input'), 0), (None, 1)])
def test_closed_worker_pipe_preserves_terminal_status_and_checks_process_exit(tmp_path, monkeypatch, closed_at, terminal, exitcode):
    from mfrm_app import project
    messages = [terminal] if terminal else []
    def poll():
        if not messages and closed_at == 'poll':
            raise BrokenPipeError('The pipe has been ended')
        return True
    def recv():
        if messages:
            return messages.pop()
        raise EOFError()
    alive = True
    worker = desktop_app.AnalysisWorker()
    worker.receiver = SimpleNamespace(poll=poll, recv=recv, close=lambda: None)
    worker.process = SimpleNamespace(is_alive=lambda: alive, join=lambda: None,
                                     close=lambda: None, exitcode=exitcode)
    worker.output_path = tmp_path / 'result.mfrm'
    loaded = []
    monkeypatch.setattr(project, 'load_project', lambda path: loaded.append(path) or 'saved result')
    assert worker.poll() == []  # A pipe closing alone is not successful completion.
    alive = False
    events = worker.poll()
    if terminal == ('complete', '') and exitcode == 0:
        assert events == [('complete', 'saved result')] and loaded == [worker.output_path]
    else:
        assert events[0][0] == 'failed' and loaded == []
        if terminal and terminal[0] == 'failed':
            assert events == [terminal]
    assert worker.process is worker.receiver is None


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
def test_windows_packaging_requires_complete_restore_and_preserves_the_bundle(tmp_path, monkeypatch, broken):
    monkeypatch.setattr(check_windows_desktop, '__file__', str(tmp_path / 'scripts/check_windows_desktop.py'))
    (tmp_path / 'docs').mkdir()
    (tmp_path / 'docs/desktop_windows_installation.md').write_text('Install notes')
    (tmp_path / 'source.py').write_text('source')
    build = tmp_path / 'build'
    bundle = build / 'dist/MFRM-Desktop'
    bundle.mkdir(parents=True)
    (bundle / 'MFRM-Desktop.exe').write_bytes(b'candidate')
    (bundle / '_internal').mkdir()
    (bundle / '_internal/runtime.dll').write_bytes(b'runtime')
    (build / 'desktop_build.json').write_text(json.dumps({'source_sha256': {
        'source.py': hashlib.sha256(b'source').hexdigest()}}))
    workflows = [('RSM', 'JMLE'), ('PCM', 'JMLE'), ('RSM', 'MML'), ('PCM', 'MML'), ('G/D', 'ANOVA')]

    def process(command, **kwargs):
        assert Path(kwargs['cwd']) != tmp_path
        if broken == 'worker':
            raise subprocess.CalledProcessError(1, command)
        if command[1] == '--smoke-test':
            out = Path(command[2]); out.mkdir()
            (out / 'smoke.json').write_text(json.dumps({'extended': [dict(model=m, method=e) for m,e in workflows]}))
            (out / '作文.mfrm').write_bytes(b'original')
            (out / 'report.html').write_bytes(b'<svg/>')
            for model, method in workflows:
                directory = out / (model.replace('/', '') + '_' + method); directory.mkdir()
                (directory / 'analysis.mfrm').write_bytes(b'original')
                (directory / 'report.html').write_bytes(b'<svg/>')
        else:
            source, out = Path(command[2]), Path(command[4]); out.mkdir(parents=True)
            (out / 'reopen.json').write_text(json.dumps({'analysis_id': 'id', 'estimation_calls': 0}))
            (out / 'report.html').write_bytes(b'changed' if broken == 'report' else b'<svg/>')
            if broken == 'original': source.write_bytes(b'changed')

    monkeypatch.setattr(check_windows_desktop.subprocess, 'run', process)
    if broken:
        with pytest.raises((ValueError, subprocess.CalledProcessError)):
            check_windows_desktop.check_build(build, tmp_path / 'checks')
        assert not (build / 'MFRM-Desktop-Windows-x64.zip').exists()
    else:
        result = check_windows_desktop.check_build(build, tmp_path / 'checks')
        assert len(result['reopened']) == 6 and not result['windows_11_gui_accepted']
        with zipfile.ZipFile(build / result['archive']) as archive:
            assert archive.read('MFRM-Desktop/_internal/runtime.dll') == b'runtime'
            assert archive.read('MFRM-Desktop/START_HERE.txt') == b'Install notes'
