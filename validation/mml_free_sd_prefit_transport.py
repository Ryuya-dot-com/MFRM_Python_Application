"""Registrar-owned reconstruction from retained, pre-fit input bytes.

Development only. The T0 input-validity artifact binds every run's estimator
configuration and the live Python environment before input generation. This
does not supply an external timestamp, a single-attempt executor, or R/FACETS
runtime attestation; the existing qualification flags remain false.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import io
import os
from pathlib import Path
import platform
import stat
import sys

import numpy as np
import pandas as pd
import scipy

from validation import mml_free_sd_qualification_registrar as registrar
from validation.mml_free_sd_quadrature_adapter import (
    app_free_sd_evaluator_implementation_sha256, app_free_sd_problem_digest,
)
from validation.mml_free_sd_stationarity_adapter import prepare_app_free_sd_problem_from_data


ROOT = Path(__file__).resolve().parents[1]
CONFIG_KEYS = {'person_col', 'facet_cols', 'score_col', 'rating_min', 'rating_max', 'model',
    'noncenter_facet', 'step_facet', 'positive_facets', 'weight_col', 'quadrature_points', 'sigma_start'}
CANONICALIZATION = dict(version='prefit-response-csv-v1', encoding='UTF-8',
    columns='Person, declared facets in order, Score, Weight', line_ending='LF',
    row_order='retained input order', identifiers='strings; normalized by frozen app factory',
    response='integer supported category', weights='unit float', missing='absent rows only')
CANONICALIZATION_SHA256 = hashlib.sha256(registrar._json_bytes(CANONICALIZATION)).hexdigest()


def _read_snapshot(path: Path, expected_sha256: str) -> bytes:
    """Read a regular file once; parsing and hashing use the same byte snapshot."""
    registrar._sha256(expected_sha256, 'expected retained-file digest')
    path = Path(path).absolute()
    if path != path.resolve():
        raise ValueError('Retained artifact path contains a link or noncanonical component')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(fd, 'rb') as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError('Retained artifact is not a regular file')
        payload = handle.read()
        after = os.fstat(handle.fileno())
    identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
    if identity(before) != identity(after) or identity(after) != identity(path.stat(follow_symlinks=False)):
        raise ValueError('Retained artifact changed while reading')
    if hashlib.sha256(payload).hexdigest() != expected_sha256:
        raise ValueError('Retained artifact byte digest differs')
    return payload


def live_python_environment() -> dict:
    """Observe sources and Python packages; never accept a worker's environment."""
    evaluator = app_free_sd_evaluator_implementation_sha256()  # also detects loaded-source drift
    paths = sorted({ROOT/'streamlit_app.py', Path(__file__).resolve(),
        *[ROOT/p for p in registrar.REQUIRED_SOURCE_FILES], *list((ROOT/'mfrm_app').glob('*.py'))})
    blas = registrar._json_bytes(np.show_config(mode='dicts')).decode('utf-8')
    # Distribution.version reparses METADATA; use one fresh parse per package.
    packages = sorted((metadata['Name'].lower(), metadata['Version'])
        for metadata in (d.metadata for d in importlib.metadata.distributions()))
    return dict(schema_version='prefit-python-environment-v1',
        sources={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        evaluator_implementation_sha256=evaluator,
        runtime=dict(python=platform.python_version(), python_implementation=platform.python_implementation(),
            os=platform.platform(), architecture=platform.machine(), byteorder=sys.byteorder,
            numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__, blas=blas,
            python_executable_sha256=hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest(),
            installed_packages=packages))


def _verified_environment(path, digest):
    payload = _read_snapshot(path, digest)
    expected = registrar._load_strict_json_object(payload, 'retained Python environment')
    actual = live_python_environment()
    if registrar._json_bytes(expected) != registrar._json_bytes(actual):
        raise ValueError('Live Python source/runtime differs from the pre-generation environment')
    return actual


def _configured_data(payload, config):
    registrar._exact_keys(config, CONFIG_KEYS, 'frozen estimator config')
    for key in ['person_col', 'score_col', 'noncenter_facet', 'model']:
        registrar._text(config[key], key)
    for key in ['facet_cols', 'positive_facets']:
        if not isinstance(config[key], list) or not all(isinstance(v, str) and v for v in config[key]):
            raise ValueError(key+' must be an explicit string list')
    if config['weight_col'] is not None:
        registrar._text(config['weight_col'], 'weight_col')
    columns = [config['person_col'], *config['facet_cols'], config['score_col']]
    if config['weight_col'] is not None:
        columns.append(config['weight_col'])
    if len(set(columns)) != len(columns):
        raise ValueError('Input roles must name distinct columns')
    try:
        rows = list(csv.reader(io.StringIO(payload.decode('utf-8'), newline=''), strict=True))
    except (UnicodeDecodeError, csv.Error) as exc:
        raise ValueError('Input must be strict UTF-8 CSV') from exc
    if len(rows) < 2 or len(set(rows[0])) != len(rows[0]) or set(rows[0]) != set(columns):
        raise ValueError('Input header or observed rows differ from the declared columns')
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError('Input CSV has ragged rows')
    data = pd.DataFrame(rows[1:], columns=rows[0])
    for name in [config['score_col'], config['weight_col']]:
        if name is not None:
            data[name] = pd.to_numeric(data[name], errors='raise')
    return data


def reconstruct_fit_authorization(*, plan_dir, inputs_dir, rules_path, environment_path,
                                  input_root, expected_plan_identity_sha256):
    """Reconstruct ALL generated runs without fitting or worker callbacks.

    The caller must obtain the expected T0 identity from its registration
    record. Local identity consistency alone is not proof of prospective time.
    Invalid/generated-failure statuses are retained as unauthorized records;
    validating their failure evidence remains a separate registrar step.
    """
    plan = registrar.validate_published_artifact(Path(plan_dir),
        validator=registrar.validate_pre_generation_plan,
        semantic_identity_function=registrar.pre_generation_plan_identity_sha256)['artifact']
    if registrar.pre_generation_plan_identity_sha256(plan) != expected_plan_identity_sha256:
        raise ValueError('Published T0 differs from the expected registration identity')
    inputs = registrar.validate_published_artifact(Path(inputs_dir),
        validator=lambda obj: registrar.validate_realized_inputs(plan, obj),
        semantic_identity_function=lambda obj: registrar.realized_inputs_identity_sha256(plan, obj))['artifact']
    if plan['input_canonicalization_sha256'] != CANONICALIZATION_SHA256:
        raise ValueError('Registered input canonicalization differs')
    rules_bytes = _read_snapshot(rules_path, plan['input_validity_rules_sha256'])
    rules = registrar._load_strict_json_object(rules_bytes, 'pre-generation reconstruction rules')
    registrar._exact_keys(rules, {'schema_version', 'environment_sha256', 'runs'}, 'reconstruction rules')
    if rules['schema_version'] != 'prefit-reconstruction-rules-v1':
        raise ValueError('Reconstruction rules schema differs')
    environment = _verified_environment(environment_path, rules['environment_sha256'])
    for name, digest in plan['source_manifest'].items():
        if environment['sources'].get(name) != digest:
            raise ValueError('T0 source manifest differs from the live reconstruction sources')
    if plan['numerical_contracts']['objective_evaluator_implementation_sha256'] != environment['evaluator_implementation_sha256']:
        raise ValueError('T0 evaluator identity differs from the live evaluator')
    for name in ['python', 'python_implementation', 'os', 'architecture', 'byteorder', 'numpy', 'scipy', 'blas']:
        if plan['runtime_manifest'][name] != environment['runtime'][name]:
            raise ValueError('T0 Python runtime differs: '+name)
    package_digest = hashlib.sha256(registrar._json_bytes(environment['runtime']['installed_packages'])).hexdigest()
    if plan['runtime_manifest']['package_lock_sha256'] != package_digest:
        raise ValueError('T0 package inventory digest differs')
    planned = plan['planned_runs']
    sensitivity = plan['numerical_contracts']['sensitivity']
    if (sensitivity['primary_quadrature_points'], sensitivity['sensitivity_quadrature_points']) != (31, 61):
        raise ValueError('The v1 registrar lanes require the declared Q31/Q61 contract')
    if not isinstance(rules['runs'], list) or len(rules['runs']) != len(planned):
        raise ValueError('Reconstruction rules do not cover every planned run')
    root = Path(input_root).resolve(strict=True)
    records, problems, evidence = [], {}, []
    for run, realized, binding in zip(planned, inputs['records'], rules['runs'], strict=True):
        registrar._exact_keys(binding, {'run_id', 'seed', 'input_file', 'estimator_config'}, 'run binding')
        if binding['run_id'] != run['run_id'] or type(binding['seed']) is not int or binding['seed'] != run['seed']:
            raise ValueError('Reconstruction rules RunId/seed differ from the plan')
        relative = Path(registrar._text(binding['input_file'], 'input_file'))
        if relative.is_absolute() or '..' in relative.parts or not relative.parts:
            raise ValueError('Input file must be a relative path inside the retained input root')
        registrar._exact_keys(binding['estimator_config'], CONFIG_KEYS, 'frozen estimator config')
        if (binding['estimator_config']['model'] not in {'RSM', 'PCM'}
            or type(binding['estimator_config']['quadrature_points']) is not int
            or binding['estimator_config']['quadrature_points'] != sensitivity['primary_quadrature_points']):
            raise ValueError('Pre-fit scope or primary quadrature differs from the registered contract')
        record = dict(run_id=run['run_id'], seed=run['seed'], input_instance_sha256=None,
                      likelihood_problem_digest=None, estimator_input_config_sha256=None)
        if realized['status'] != 'GENERATED':
            record['status'] = ('NOT_AUTHORIZED_GENERATION_FAILED' if realized['status'] == 'GENERATION_FAILED'
                                else 'NOT_AUTHORIZED_INPUT_INVALID')
        else:
            payload = _read_snapshot(root/relative, realized['raw_input_artifact_sha256'])
            config = binding['estimator_config']
            data = _configured_data(payload, config)
            prep, problem = prepare_app_free_sd_problem_from_data(data, **config)
            columns = ['Person', *problem.config['facet_names'], 'Score', 'Weight']
            canonical = prep['data'][columns].to_csv(index=False, lineterminator='\n').encode('utf-8')
            canonical_digest = hashlib.sha256(canonical).hexdigest()
            if len(data) != realized['rows'] or canonical_digest != realized['canonical_response_sha256']:
                raise ValueError('Retained response rows or canonical identity differ from T1')
            record.update(status='AUTHORIZED_FOR_SINGLE_ATTEMPT_DEVELOPMENT_ONLY',
                input_instance_sha256=realized['input_instance_sha256'],
                likelihood_problem_digest=app_free_sd_problem_digest({'prep': prep}, problem),
                estimator_input_config_sha256=hashlib.sha256(registrar._json_bytes(config)).hexdigest())
            problems[run['run_id']] = problem
            evidence.append(dict(run_id=run['run_id'], seed=run['seed'], raw_input_sha256=hashlib.sha256(payload).hexdigest(),
                raw_input_bytes=len(payload), canonical_response_sha256=canonical_digest, rows=len(data),
                structural_start=list(problem.structural_start), sigma_start=problem.sigma_start,
                quadrature_points=problem.quadrature_points, scientific_inference_ready=False))
        records.append(record)
    # Detect mutations during parsing/reconstruction; use no caller evaluator.
    _verified_environment(environment_path, rules['environment_sha256'])
    _read_snapshot(rules_path, plan['input_validity_rules_sha256'])
    for realized, binding in zip(inputs['records'], rules['runs'], strict=True):
        if realized['status'] == 'GENERATED':
            _read_snapshot(root/binding['input_file'], realized['raw_input_artifact_sha256'])
    final_plan = registrar.validate_published_artifact(Path(plan_dir),
        validator=registrar.validate_pre_generation_plan,
        semantic_identity_function=registrar.pre_generation_plan_identity_sha256)['artifact']
    final_inputs = registrar.validate_published_artifact(Path(inputs_dir),
        validator=lambda obj: registrar.validate_realized_inputs(plan, obj),
        semantic_identity_function=lambda obj: registrar.realized_inputs_identity_sha256(plan, obj))['artifact']
    if final_plan != plan or final_inputs != inputs:
        raise ValueError('Published T0/T1 changed during reconstruction')
    authorization = dict(schema_version=registrar.FIT_SCHEMA, status='DEVELOPMENT_ONLY_PREFIT_FACTORY_NOT_QUALIFIED',
        plan_identity_sha256=expected_plan_identity_sha256,
        realized_inputs_identity_sha256=registrar.realized_inputs_identity_sha256(plan, inputs),
        planned_runs_sha256=plan['planned_runs_sha256'],
        problem_digest_factory_implementation_sha256=hashlib.sha256(registrar._json_bytes(environment)).hexdigest(),
        records=records, records_sha256=registrar._records_digest('MFRM_QUALIFICATION_FIT_AUTHORIZATION_RECORDS', 'v1', records),
        estimator_attempts_started_before_authorization=False, pure_prefit_problem_factory_qualified=False,
        scientific_inference_enabled=False)
    registrar.validate_fit_authorization(plan, inputs, authorization)
    return dict(authorization=authorization, problems=problems, reconstruction_evidence=evidence,
                failure_records_requiring_independent_audit=[r for r in inputs['records'] if r['status'] != 'GENERATED'],
                environment_sha256=rules['environment_sha256'], scientific_inference_ready=False)
