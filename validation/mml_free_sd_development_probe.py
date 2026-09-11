#!/usr/bin/env python3
"""Replay observed RSM/PCM development fixtures at Q31/Q61; never qualification.

Reuse the retained test fixture and its numerical tolerances, with Q changed
from 9/15 to 31/61. Requires the development dependencies (including pytest).
The output is a new, exclusive JSON file; frozen studies are not changed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import runpy
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import streamlit_app as app
from mfrm_app.mml_stationarity import make_joint_free_sd_functions, JointPolishOptions
from mfrm_app.mml_quadrature_sensitivity import QuadratureSensitivityContract
from validation.mml_free_sd_quadrature_adapter import run_app_free_sd_quadrature_sensitivity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output already exists; use a new path")
    source_paths = [
        "streamlit_app.py", "mfrm_app/mml_stationarity.py", "mfrm_app/mml_engine_v2.py",
        "mfrm_app/mml_quadrature_sensitivity.py", "mfrm_app/mml_qualification_batch.py",
        "validation/mml_free_sd_stationarity_adapter.py",
        "validation/mml_free_sd_quadrature_adapter.py",
        "tests/test_mml_free_sd_quadrature_adapter.py",
        "validation/mml_free_sd_development_probe.py",
    ]
    hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in source_paths}
    # Reuse an already-observed development fixture and its numerical contract.
    test = runpy.run_path(str(ROOT / 'tests/test_mml_free_sd_quadrature_adapter.py'))
    source = test['native_result'].__wrapped__()
    stationarity_contract = test['development_stationarity_contract']()
    sensitivity_contract = QuadratureSensitivityContract(
        primary_quadrature_points=31, sensitivity_quadrature_points=61,
        max_structural_parameter_difference=0.1, max_log_sigma_difference=0.1,
        max_sensitivity_optimization_gain_per_observation=1e-4,
        max_negative_sensitivity_gain_per_observation=1e-8,
        max_abs_primary_back_evaluation_change_per_observation=1e-4,
        max_objective_reconstruction_disagreement=1e-8,
    )
    records = []
    for model in ('RSM', 'PCM'):
        result = app.mfrm_estimate(
            source['prep']['data'], person_col='Person', facet_cols=['Rater', 'Task', 'Criterion'],
            score_col='Score', rating_min=0, rating_max=3, model=model, method='MML',
            noncenter_facet='Criterion', step_facet='Criterion', mml_engine='EM',
            quad_points=31, estimate_population_sd=True, population_prior_sd=1.0,
            maxit=300, reltol=1e-6, min_obs_per_element=1, min_obs_per_category=1,
        )
        bundle = run_app_free_sd_quadrature_sensitivity(
            result, primary_quadrature_points=31, sensitivity_quadrature_points=61,
            stationarity_contract=stationarity_contract, sensitivity_contract=sensitivity_contract,
            options=JointPolishOptions(maxiter=250, gtol=1e-8, ftol=1e-15),
        )
        problem = bundle.primary_problem
        _, score = make_joint_free_sd_functions(problem.value, problem.value_gradient)
        _, gradient = score(np.r_[problem.structural_start, np.log(problem.sigma_start)])
        assessment = bundle.assessment
        payload = {
            'classification': 'OBSERVED_DEVELOPMENT_ONLY_EXCLUDED_FROM_QUALIFICATION',
            'model': model, 'n_person': result['config']['n_person'], 'n_ratings': result['prep']['n_obs'],
            'legacy_success': bool(result['opt'].success),
            'legacy_inference_ready': bool(result['summary'].iloc[0]['InferenceReady']), 'legacy_sigma': problem.sigma_start,
            'legacy_score_supnorm': float(np.max(np.abs(gradient))),
            'q31_sigma': bundle.primary_run.final_sigma, 'q61_sigma': bundle.sensitivity_run.final_sigma,
            'q31_polish_gain': bundle.primary_run.primary_polish.objective_improvement,
            'q31_score_supnorm': bundle.primary_run.final_projected_gradient_supnorm,
            'q61_score_supnorm': bundle.sensitivity_run.final_projected_gradient_supnorm,
            'parameter_difference': assessment.maximum_structural_parameter_difference,
            'log_sigma_difference': assessment.absolute_log_sigma_difference,
            'q61_gain_per_observation': assessment.sensitivity_optimization_gain_per_observation,
            'q31_stationarity': assessment.primary_stationarity_pass,
            'q61_stationarity': assessment.sensitivity_stationarity_pass,
            'numerical_sensitivity_pass': assessment.numerical_sensitivity_pass,
            'scientific_inference_ready': assessment.scientific_inference_ready,
            'raw': assessment.to_dict(),
        }
        assert not payload['legacy_inference_ready']
        assert not payload['scientific_inference_ready']
        records.append(payload)
        print(json.dumps({key: value for key, value in payload.items() if key != 'raw'}), flush=True)

    if hashes != {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in source_paths}:
        raise RuntimeError("source changed during the probe")
    payload = {
        "classification": "OBSERVED_DEVELOPMENT_ONLY_EXCLUDED_FROM_QUALIFICATION",
        "qualification_eligible": False,
        "scientific_inference_ready": False,
        "source_sha256": hashes,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "fixture": {"seed": 20260811, "persons": 24, "ratings": 192, "generation_model": "RSM"},
        "records": records,
    }
    encoded = json.dumps(payload, indent=2, allow_nan=False) + "\n"
    with args.output.open("x", encoding="utf-8") as handle:
        handle.write(encoded)


if __name__ == "__main__":
    main()
