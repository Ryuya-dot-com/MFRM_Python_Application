"""Observed development fixtures, never prospective qualification data."""
import hashlib

import numpy as np
import pandas as pd
import pytest

import streamlit_app as app
from tests.test_mml_free_sd_qualification_registrar import _plan, _inputs
from validation import mml_free_sd_prefit_transport as transport
from validation import mml_free_sd_qualification_registrar as registrar
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize('change', ['version', 'added', 'removed'])
def test_package_inventory_preserves_identity_and_detects_changes(tmp_path, monkeypatch, change):
    distributions = []
    for index, (name, version) in enumerate([('Mixed_Package', '2.0'), ('Alpha', '1.0'), ('Alpha', '1.1')]):
        directory = tmp_path / f'package{index}.dist-info'
        directory.mkdir()
        (directory / 'METADATA').write_text(f'Name: {name}\nVersion: {version}\n\n')
        distributions.append(transport.importlib.metadata.PathDistribution(directory))
    monkeypatch.setattr(transport.importlib.metadata, 'distributions', lambda: iter(distributions))
    expected = sorted((d.metadata['Name'].lower(), d.version) for d in distributions)
    environment = transport.live_python_environment()
    assert environment['runtime']['installed_packages'] == expected
    assert len(expected) == 3  # Duplicate names remain part of the inventory.
    payload = registrar._json_bytes(environment)
    path = tmp_path / 'environment.json'
    path.write_bytes(payload)
    transport._verified_environment(path, digest(payload))
    if change == 'version':
        (tmp_path / 'package0.dist-info/METADATA').write_text('Name: Mixed_Package\nVersion: 2.1\n\n')
    elif change == 'added':
        distributions.append(distributions[0])
    else:
        distributions.pop()
    with pytest.raises(ValueError, match='Live Python source/runtime differs'):
        transport._verified_environment(path, digest(payload))


def fixture(tmp_path, *, environment_fault=None, canonical_fault=False, same_model=False, primary_q=31):
    tmp_path = tmp_path.resolve()
    plan = _plan()
    for i, run in enumerate(plan['planned_runs']):
        run['seed'] = 2**100+i  # Preserve integer identity beyond float64 and int64.
    plan['planned_runs_sha256'] = registrar._records_digest('MFRM_QUALIFICATION_PLANNED_RUNS', 'v1', plan['planned_runs'])
    environment = transport.live_python_environment()
    if environment_fault == 'source': environment['sources']['streamlit_app.py'] = '0'*64
    if environment_fault == 'runtime': environment['runtime']['numpy'] = '99.0.0'
    env_path = tmp_path/'environment.json'; env_bytes = registrar._json_bytes(environment); env_path.write_bytes(env_bytes)
    plan['source_manifest'] = {p: environment['sources'][p] for p in registrar.REQUIRED_SOURCE_FILES}
    for key in ['python', 'python_implementation', 'os', 'architecture', 'byteorder', 'numpy', 'scipy', 'blas']:
        plan['runtime_manifest'][key] = environment['runtime'][key]
    plan['runtime_manifest']['package_lock_sha256'] = digest(registrar._json_bytes(environment['runtime']['installed_packages']))
    plan['numerical_contracts']['objective_evaluator_implementation_sha256'] = environment['evaluator_implementation_sha256']
    for name, domain in [('source_manifest', 'SOURCE_MANIFEST'), ('runtime_manifest', 'RUNTIME_MANIFEST'),
                         ('numerical_contracts', 'NUMERICAL_CONTRACTS')]:
        plan[name+'_sha256'] = registrar._domain_sha256('MFRM_QUALIFICATION_'+domain, 'v1', plan[name])
    config = dict(person_col='Person', facet_cols=['Rater', 'Criterion'], score_col='Score',
        rating_min=0, rating_max=2, model='RSM', noncenter_facet='Criterion', step_facet=None,
        positive_facets=[], weight_col=None, quadrature_points=primary_q, sigma_start=1.)
    bindings = []
    for i, run in enumerate(plan['planned_runs']):
        cfg = dict(config)
        if i == 1 and not same_model: cfg.update(model='PCM', step_facet='Criterion')
        bindings.append(dict(run_id=run['run_id'], seed=run['seed'], input_file=run['run_id']+'.csv', estimator_config=cfg))
    rules = dict(schema_version='prefit-reconstruction-rules-v1', environment_sha256=digest(env_bytes), runs=bindings)
    rules_path = tmp_path/'rules.json'; rules_bytes = registrar._json_bytes(rules); rules_path.write_bytes(rules_bytes)
    plan['input_validity_rules_sha256'] = digest(rules_bytes)
    plan['input_canonicalization_sha256'] = transport.CANONICALIZATION_SHA256
    identity = registrar.pre_generation_plan_identity_sha256(plan)
    registrar.publish_sealed_artifact(output_dir=tmp_path/'T0', document=plan,
        validator=registrar.validate_pre_generation_plan, semantic_identity_function=registrar.pre_generation_plan_identity_sha256)
    inputs = _inputs(plan)
    raw_root = tmp_path/'raw'; raw_root.mkdir()
    data = pd.DataFrame([(person, f'R{r}', f'C{c}', (p+r+c) % 3)
        for p, person in enumerate(['001', '010', 'NA']) for r in range(2) for c in range(2)],
        columns=['Person', 'Rater', 'Criterion', 'Score'])
    expected = {}
    for record, binding in zip(inputs['records'], bindings, strict=True):
        if record['status'] != 'GENERATED': continue
        payload = data.to_csv(index=False, lineterminator='\n').encode(); (raw_root/binding['input_file']).write_bytes(payload)
        prep, problem = prepare_app_free_sd_problem_from_data(data, **binding['estimator_config'])
        canonical = prep['data'][['Person', 'Rater', 'Criterion', 'Score', 'Weight']].to_csv(index=False, lineterminator='\n').encode()
        record.update(raw_input_artifact_sha256=digest(payload), canonical_response_sha256='0'*64 if canonical_fault else digest(canonical), rows=len(data))
        record['input_instance_sha256'] = registrar.input_instance_sha256(plan_identity_sha256=identity, run_id=record['run_id'],
            seed=record['seed'], raw_input_artifact_sha256=record['raw_input_artifact_sha256'], canonical_response_sha256=record['canonical_response_sha256'])
        expected[record['run_id']] = problem
    inputs['records_sha256'] = registrar._records_digest('MFRM_QUALIFICATION_REALIZED_INPUT_RECORDS', 'v1', inputs['records'])
    registrar.publish_sealed_artifact(output_dir=tmp_path/'T1', document=inputs,
        validator=lambda obj: registrar.validate_realized_inputs(plan, obj),
        semantic_identity_function=lambda obj: registrar.realized_inputs_identity_sha256(plan, obj))
    kwargs = dict(plan_dir=tmp_path/'T0', inputs_dir=tmp_path/'T1', rules_path=rules_path,
        environment_path=env_path, input_root=raw_root, expected_plan_identity_sha256=identity)
    return kwargs, expected


def test_reconstructs_both_models_from_retained_bytes_without_fitting(tmp_path, monkeypatch):
    kwargs, expected = fixture(tmp_path)
    def forbidden(*args, **kw): pytest.fail('Reconstruction called an optimizer')
    monkeypatch.setattr(app, 'mfrm_estimate', forbidden); monkeypatch.setattr(app, 'minimize', forbidden)
    result = transport.reconstruct_fit_authorization(**kwargs)
    assert not result['scientific_inference_ready']
    assert not result['authorization']['pure_prefit_problem_factory_qualified']
    assert len(result['authorization']['records']) == 4
    assert len(result['failure_records_requiring_independent_audit']) == 2
    assert [r['seed'] for r in result['authorization']['records']] == [2**100+i for i in range(4)]
    for run, problem in result['problems'].items():
        point = np.r_[np.linspace(-.3, .4, len(problem.structural_start)), np.log(1.2)]
        nll, gradient = problem.joint_value_gradient(point)
        expected_nll, expected_gradient = expected[run].joint_value_gradient(point)
        assert nll == expected_nll
        np.testing.assert_array_equal(gradient, expected_gradient)
        assert problem.config['n_person'] == 3


def test_coincident_likelihoods_keep_distinct_registered_input_instances(tmp_path):
    kwargs, _ = fixture(tmp_path, same_model=True)
    rows = transport.reconstruct_fit_authorization(**kwargs)['authorization']['records'][:2]
    assert rows[0]['likelihood_problem_digest'] == rows[1]['likelihood_problem_digest']
    assert rows[0]['input_instance_sha256'] != rows[1]['input_instance_sha256']


def test_factory_cannot_authorize_a_different_quadrature_from_the_registered_lane(tmp_path):
    kwargs, _ = fixture(tmp_path, primary_q=9)
    with pytest.raises(ValueError, match='primary quadrature'):
        transport.reconstruct_fit_authorization(**kwargs)


@pytest.mark.parametrize('fault', ['raw_bytes', 'rules', 'environment_bytes', 'source', 'runtime', 'canonical', 'parent_identity', 'link'])
def test_reconstruction_rejects_changed_or_forged_evidence(tmp_path, fault):
    kwargs, _ = fixture(tmp_path, environment_fault=fault, canonical_fault=fault == 'canonical')
    if fault == 'raw_bytes': (kwargs['input_root']/'run-001.csv').write_text('changed')
    elif fault == 'rules': kwargs['rules_path'].write_bytes(kwargs['rules_path'].read_bytes()+b' ')
    elif fault == 'environment_bytes': kwargs['environment_path'].write_bytes(kwargs['environment_path'].read_bytes()+b' ')
    elif fault == 'parent_identity': kwargs['expected_plan_identity_sha256'] = '0'*64
    elif fault == 'link':
        original = kwargs['input_root']/'run-001.csv'; moved = kwargs['input_root']/'moved.csv'
        original.rename(moved); original.symlink_to(moved)
    with pytest.raises(ValueError): transport.reconstruct_fit_authorization(**kwargs)


def test_changed_input_during_factory_is_rejected(tmp_path, monkeypatch):
    kwargs, _ = fixture(tmp_path)
    original = transport.prepare_app_free_sd_problem_from_data
    def mutating(data, **cfg):
        result = original(data, **cfg)
        path = kwargs['input_root']/'run-001.csv'; path.write_bytes(path.read_bytes()+b'\n')
        return result
    monkeypatch.setattr(transport, 'prepare_app_free_sd_problem_from_data', mutating)
    with pytest.raises(ValueError, match='byte digest'):
        transport.reconstruct_fit_authorization(**kwargs)


@pytest.mark.parametrize('payload', [b'Person,Rater,Criterion,Score,Score\np,r,c,1,1\n',
    b'Person,Rater,Criterion,Score\np,r,c,1,2\n', b'Person,Rater,Criterion,Score\n"p,r,c,1\n'])
def test_strict_csv_does_not_silently_repair_rows_or_columns(payload):
    config = dict(person_col='Person', facet_cols=['Rater', 'Criterion'], score_col='Score', rating_min=0, rating_max=2,
        model='RSM', noncenter_facet='Criterion', step_facet=None, positive_facets=[], weight_col=None, quadrature_points=9, sigma_start=1.)
    with pytest.raises(ValueError): transport._configured_data(payload, config)
