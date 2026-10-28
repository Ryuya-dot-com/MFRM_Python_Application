"""An implementation switch must bind a verified receipt and exact sources."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

from validation import mml_execution_amendment as amendment


@pytest.fixture
def registered(tmp_path, monkeypatch):
    # Preserve the sealed macOS root naming/checks inside pytest's isolated directory.
    monkeypatch.setattr(amendment, 'Path', lambda value: tmp_path if value == '/private/tmp' else Path(value))
    def write(path,value): path.write_text(json.dumps(value,sort_keys=True))
    study=tmp_path/'study';study.mkdir();write(study/'T0.json',dict(test=str(tmp_path)))
    write(study/'source_manifest.json',dict(original='retained'))
    identity=amendment.digest(study/'T0.json')
    folder=tmp_path/'amendment';folder.mkdir()
    binary=tmp_path/'openssl';binary.write_bytes(b'unit-test executable identity')
    (folder/'cacert.pem').write_bytes(b'unit-test CA identity')
    source=b'checked source';source_sha=__import__('hashlib').sha256(source).hexdigest()
    write(folder/'sources.json',{'code.py':source_sha})
    write(folder/'amendment.json',dict(change_scope='implementation_only_same_estimand_design_seeds_solver_options_and_criteria',
        studies={identity:dict(original_source_manifest_sha256=amendment.digest(study/'source_manifest.json'),candidate_source_manifest='sources.json')}))
    seal=dict(schema='mfrm-numerical-implementation-amendment-v1',openssl_sha256=amendment.digest(binary),
        files={p.name:amendment.digest(p) for p in folder.iterdir()})
    write(folder/'T0.json',seal);sha=amendment.digest(folder/'T0.json')
    write(study/'implementation_activation.json',dict(schema='mfrm-implementation-activation-v1',study_T0_sha256=identity,
        amendment_T0_sha256=sha,amendment_directory=str(folder)))
    root=amendment.expected_root(study,identity,tmp_path/'old');root.mkdir(exist_ok=False)
    (root/'code.py').write_bytes(source)
    calls=[]
    monkeypatch.setattr(amendment.subprocess,'run',lambda command,**kwargs:calls.append(command))
    try:yield dict(study=study,root=root,folder=folder,binary=binary,calls=calls,identity=identity,
                   ca_sha=amendment.digest(folder/'cacert.pem'))
    finally:shutil.rmtree(root)


def verify(case):
    return amendment.verified_sources(case['study'],case['root'],case['binary'],case['ca_sha'])


def test_no_activation_preserves_original_route(tmp_path):
    (tmp_path/'T0.json').write_text('original')
    original=tmp_path/'original'
    assert amendment.expected_root(tmp_path,amendment.digest(tmp_path/'T0.json'),original)==original
    assert amendment.verified_sources(tmp_path,original,tmp_path/'absent','unused') is None


def test_amended_route_verifies_query_and_exact_t0_and_returns_binding(registered):
    hashes,binding=verify(registered)
    assert set(hashes)=={'code.py'} and binding['study_T0_sha256']==registered['identity']
    assert binding['activation_sha256']==amendment.digest(registered['study']/'implementation_activation.json')
    assert len(registered['calls'])==2
    assert '-queryfile' in registered['calls'][0] and '-data' in registered['calls'][1]


@pytest.mark.parametrize('target',['code','extra','manifest','plan','seal','ca','binary','activation'])
def test_tampered_activation_sources_or_signed_artifacts_fail(registered,target):
    paths=dict(code=registered['root']/'code.py',extra=registered['root']/'extra.py',
        manifest=registered['folder']/'sources.json',plan=registered['folder']/'amendment.json',
        seal=registered['folder']/'T0.json',ca=registered['folder']/'cacert.pem',binary=registered['binary'],
        activation=registered['study']/'implementation_activation.json')
    paths[target].write_text('{}')
    with pytest.raises(ValueError):verify(registered)


def test_receipt_failure_cannot_return_a_source_binding(registered,monkeypatch):
    def failed(command,**kwargs):raise subprocess.CalledProcessError(1,command,stderr='receipt mismatch')
    monkeypatch.setattr(amendment.subprocess,'run',failed)
    with pytest.raises(subprocess.CalledProcessError):verify(registered)
