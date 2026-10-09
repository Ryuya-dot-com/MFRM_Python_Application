"""Verify a signed implementation amendment without replacing a study seal.

Original protocols, seeds and append-only ledgers retain their original T0.
An explicit activation selects an independently sealed source bundle; new
attempts record both identities. No activation means the original execution.
"""
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    path=Path(path);raw=path.read_bytes()
    if len(raw)!=path.stat().st_size:raise OSError('Incomplete file read: '+str(path))
    return hashlib.sha256(raw).hexdigest()


def activation(study,identity):
    path=Path(study)/'implementation_activation.json'
    if not path.exists():return None
    value=json.loads(path.read_bytes())
    if value.get('schema')!='mfrm-implementation-activation-v1' or value.get('study_T0_sha256')!=identity:
        raise ValueError('Execution activation belongs to another study or schema')
    sha=value.get('amendment_T0_sha256','')
    if len(sha)!=64 or any(c not in '0123456789abcdef' for c in sha):
        raise ValueError('Invalid amendment identity')
    return value


def expected_root(study,identity,original_root):
    value=activation(study,identity)
    return (Path('/private/tmp')/('mfrm-efficient-'+value['amendment_T0_sha256'][:16]+'-'+identity[:16])
            if value else Path(original_root))


def verified_sources(study,root,openssl,trusted_ca_sha256):
    """Verify both binding and receipt before the caller creates a verifier.

The caller already verifies the original study seal and its protocol. Source
and runtime checks still run before every numerical phase in that caller.
"""
    study=Path(study);identity=digest(study/'T0.json');value=activation(study,identity)
    if value is None:return None
    directory=Path(value['amendment_directory'])
    if digest(directory/'T0.json')!=value['amendment_T0_sha256']:
        raise ValueError('Amendment seal changed')
    seal=json.loads((directory/'T0.json').read_bytes())
    if seal.get('schema')!='mfrm-numerical-implementation-amendment-v1':
        raise ValueError('Unknown implementation amendment schema')
    if not {'amendment.json','cacert.pem'}.issubset(seal['files']):
        raise ValueError('Incomplete amendment inventory')
    if digest(openssl)!=seal['openssl_sha256'] or digest(directory/'cacert.pem')!=trusted_ca_sha256:
        raise ValueError('Timestamp verifier or trusted root changed')
    for name,sha in seal['files'].items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts or digest(directory/relative)!=sha:
            raise ValueError('Amendment artifact changed: '+name)
    for option,target in [('-queryfile',directory/'query.tsq'),('-data',directory/'T0.json')]:
        subprocess.run([str(openssl),'ts','-verify',option,str(target),'-in',str(directory/'response.tsr'),
                        '-CAfile',str(directory/'cacert.pem')],capture_output=True,text=True,check=True)
    plan=json.loads((directory/'amendment.json').read_bytes())
    if plan.get('change_scope')!='implementation_only_same_estimand_design_seeds_solver_options_and_criteria':
        raise ValueError('Amendment is outside the numerical implementation scope')
    registered=plan['studies'].get(identity)
    if registered is None or registered['original_source_manifest_sha256']!=digest(study/'source_manifest.json'):
        raise ValueError('Amendment does not bind the original study sources')
    source_name=registered['candidate_source_manifest']
    if source_name not in seal['files']:raise ValueError('Candidate sources are not sealed')
    hashes=json.loads((directory/source_name).read_bytes())
    if Path(root).resolve()!=expected_root(study,identity,root):
        raise ValueError('Amended execution requires its isolated source root')
    if {str(p.relative_to(root)) for p in Path(root).rglob('*') if p.is_file()}!=set(hashes):
        raise ValueError('Unexpected file inventory in amended source root')
    for name,sha in hashes.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts or digest(Path(root)/relative)!=sha:
            raise ValueError('Amended source changed: '+name)
    binding=dict(schema='mfrm-run-implementation-v1',study_T0_sha256=identity,
        amendment_T0_sha256=value['amendment_T0_sha256'],
        activation_sha256=digest(study/'implementation_activation.json'),
        source_manifest_sha256=seal['files'][source_name])
    return hashes,binding
