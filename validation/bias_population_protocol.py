"""Draft G3 operating-characteristic specification; no data generation.

Not yet sealed or authorized for qualification execution. This complements,
and cannot retroactively amend, the separately frozen population-recovery T0.
"""
from dataclasses import asdict
import numpy as np
from validation.mml_population_protocol import design,scenario_id,make_plan
from validation.bias_population_design import contrast_matrix,GENERATOR_VERSION


def make_bias_plan():
    scenarios=[]
    for model in ['RSM','PCM','GPCM']:
        for assignment in ['regular','sparse','unbalanced']:
            for persons,categories in [(48,4),(192,4),(192,2),(192,6)]:
                spec=design(model,assignment,persons=persons,categories=categories)
                for magnitude in [0.,.8]:
                    scenarios.append(dict(id=scenario_id(spec)+f'_B{round(100*magnitude):03d}',
                        design=asdict(spec),eta_shift=contrast_matrix(3,2,magnitude).tolist(),
                        null=magnitude==0.,replicates=2000,independent_reference_replicates=[1,1001]))
    assert len(scenarios)==72 and len({s['id'] for s in scenarios})==72
    return dict(schema='bias-population-operating-draft-v1',study_id='bias_population_operating_20261004_v1',
        classification='DRAFT_NOT_REGISTERED',scientific_inference_ready=False,
        qualification_execution_allowed=False,master_seed=2**119+307922016491,
        generator=GENERATOR_VERSION,scenarios=scenarios,planned_population_datasets=144000,
        seed_rule='replicate_seed(master_seed, study_id+":"+scenario_id, replicate); fixed order-independent 128-bit seed',
        generation='Fresh iid normal Person effects; independent assignment and response streams; no sample centering or replacement',
        scope=dict(pair=['Rater','Task'],calibration='unit-weight unpenalized joint free-SD MML in the additive baseline model',
            facets='native sum-zero Rater/Task and noncentered Criterion; Criterion-owned steps/slopes',
            shifts='row/column sum-zero 2x2 contrast embedded in the declared 3x2 pair; no fitted model interaction added',
            lanes=['known_calibration_known_person','known_calibration_EAP','fitted_calibration_known_person','fitted_calibration_EAP'],
            primary_lane='fitted_calibration_EAP',
            EAP_model='Both known/fitted EAP lanes score under the additive baseline, excluding the planted interaction. Known baseline EAP is not the posterior under the entire alternative DGP.',
            oracle='Only known additive baseline parameters AND true Person values identify each scalar target with its planted eta shift',
            other_lanes='Residual-screen operating behavior; do not assign DGP-effect coverage or sampling-SE validity to plugin intervals',
            remaining_scopes=['JMLE','independent calibration samples','Person-containing pairs','other facet pairs','anchors','penalized or nonunit-weight inference']),
        estimator=make_plan()['estimator'],
        independent_reference=dict(replicates=[1,1001],
            limits=dict(estimate=1e-8,se=1e-8,lr=1e-8,p=1e-8,ci=1e-7,adjustment=1e-14),
            family_input='Complete retained per-row eta/score/slope/steps for one calibration/Person basis; R derives observed and extreme-level cells independently',
            cells='Every declared Rater x Task cell in all four lanes; include missing screens and unobserved cells in required-slot accounting',
            adjustments='Holm and descriptive BH on the observed-cell family including failed/excluded cells',
            derived_adjustments='Also use independently computed R cell p-values; tolerance = observed family size * p tolerance + adjustment tolerance, with every threshold-crossing mismatch retained',
            failures='One reference attempt per selected run/lane; no retry. Missing input, failed R, missing/invalid output and numerical mismatch all prevent reference qualification.',
            nominal_required_families=576,nominal_required_scalar_cells=3456,nominal_required_adjustments=1152,
            nominal_required_derived_adjustments=1152),
        deterministic_controls=dict(scalar_cases=17,family_cases=['RSM','PCM','GPCM'],
            source='validation/bias_negative_controls.py',
            role='Check native/independent-R formulas and guards. Fixed inputs do not contribute independent Monte Carlo repetitions.',
            distinctions=['No finite scalar MLE vs finite root outside search range vs range-limited interval',
                'Nonunit-weight curvature without a chi-square p/interval',
                'Unobserved cell vs extreme entire facet level vs all-extreme scores in only one cell',
                'Known fixed thresholds permit a scalar MLE with missing internal observed categories; this differs from estimating those thresholds']),
        screen=dict(alpha=.05,practical_logit=.5,min_positive_weight_observations=5,max_abs=10.,omit_extreme=True,
            primary_rule='Strong review: finite unit-weight LR candidate, Holm p<.05, abs(bias)>=.5, support>=5',
            secondary_rules=['raw LR p<.05','Holm p<.05','Strong review or Statistical review'],
            review_flag='Practical, sparse and unavailable review labels are not statistical rejections',
            multiplicity='One observed Rater x Task family per replicate before excluded/failed-cell handling; BH descriptive only',
            missing='Declared but unobserved cells remain unavailable in planned-cell evidence, outside the actual observed-cell Holm family'),
        endpoints=dict(null='Global-null false-positive probability per cell and any cell in the family, for every declared rule',
            effect='Correct-direction detection per planted cell and any planted cell; report no-shift-cell flags as off-target detections',
            family='Partially available families with an observed rejection have a known union event. Report fully available-family rates separately.',
            missing='Conditional rates plus all-planned partial-identification bounds and confidence envelope; never unavailable=nonrejection',
            recovery='Only the oracle lane: cell bias/RMSE, RMS SE/empirical SD, LR-cutoff interval coverage, interval availability',
            pairwise='Record native target sign, (b1-b2)/sign and conditional disjoint-cell SE; no pairwise p-value',
            failure='Generation, calibration, quadrature, unobserved cells, extreme-level exclusion, no finite cell MLE, range limitation, reference failure',
            paired_lanes='Shared input across lanes; any numerical between-lane contrast uses replicate-level paired differences'),
        draft_interpretation_limits=dict(global_null_family_error_upper=.065,unavailability_upper=.02,
            global_null_cell_error_upper=.065,
            family_error_rules=['holm','statistical_review','strong_review'],
            raw_family_error='Descriptive only: a per-cell alpha does not set the unadjusted union-event error to alpha',
            alternative_any_correct_strong_detection_lower=.8,oracle_absolute_bias=.1,
            oracle_SE_ratio=[.9,1.1],oracle_coverage=[.925,.975],
            decisions='Pointwise 95% MC intervals/envelopes must lie within the declared limits; overlap is inconclusive. No automatic promotion or simultaneous condition-grid guarantee.',
            rationale='Family error allows 1.5 percentage points above alpha; positive usefulness also requires power, avoiding qualification of a never-rejecting screen. Any-cell power is not every-cell sensitivity.'),
        precision=dict(independent_unit='replicate',replicates=2000,
            nominal_null_mcse=float(np.sqrt(.05*.95/2000)),nominal_power_mcse=float(np.sqrt(.8*.2/2000)),
            binary='Exact binomial intervals; unknown outcomes give identification bounds and a 95% confidence envelope',
            continuous='Replicate-level MCSE; joint-moment delta method for RMS model-SE/empirical-SD ratio',
            no_interim_endpoint_peeking=True,no_replacement=True,no_optional_extension=True),
        prerequisites=['Separate source/runtime/seed/analysis seal before qualification data',
            'Single-attempt retained-input execution and complete-scenario aggregation',
            'Independent-reference failure denominator and negative controls',
            'Resource and completed-endpoint verification before qualification decisions'],
        exclusions='All earlier development samples, exact-binomial enumeration and the 18-case transport control are excluded; no post-hoc addition of bias endpoints to the frozen MML recovery cohort',
        literature=[dict(zotero_key='U9PYTVPL',title='Myford & Wolfe, Detecting and measuring rater effects, Part II',
                year=2004,read_scope='Zotero metadata and public abstract/bibliographic record, not attachment full text',
                url='https://pubmed.ncbi.nlm.nih.gov/15064538/'),
            dict(zotero_key='PKQMUBH7',title='Morris, White & Crowther (2019)',
                url='https://arxiv.org/html/1712.03198v3',role='Simulation targets, failure denominators and Monte Carlo precision')])


def make_registered_bias_plan():
    """Execution specification; a local dictionary is not a signed receipt."""
    plan=make_bias_plan()
    plan.update(schema='bias-population-operating-plan-v1',classification='PROSPECTIVE_PLAN_REQUIRES_VERIFIED_RECEIPT',
        qualification_execution_allowed=True,
        execution=dict(maximum_workers=2,disk_reserve_bytes=4*1024**3,
            planning_bytes_per_remaining_replicate=128*1024,reference_storage_allowance_bytes=1024**3,
            capacity_policy='Before a scenario starts, require headroom for all remaining runs plus references and the reserve. Before each phase enforce the reserve. Interrupted attempts stay unavailable; no replacement.',
            source_policy='Isolated read-only source bundle, Python/R runtime identity, signed exact T0 bytes and matching RFC3161 query required before any generation',
            retry_policy='Started attempts are never numerically replayed; only unstarted runs may continue',
            controls='All fixed formula/guard controls must pass before any population run',
            completion='Every planned replicate is terminal before scenario endpoints; all required reference slots remain in the evidence packet; no automatic public qualification'))
    return plan


def planned_bias_runs(plan):
    from validation.mfrm_recovery_design import replicate_seed
    for scenario in plan['scenarios']:
        for replicate in range(1,scenario['replicates']+1):
            yield dict(run_id=f"{scenario['id']}__{replicate:04d}",scenario_id=scenario['id'],replicate=replicate,
                seed=replicate_seed(plan['master_seed'],plan['study_id']+':'+scenario['id'],replicate))
