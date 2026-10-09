"""Pre-generation specification for the high-order population recovery study.

This is a new study of joint free-SD MML. It does not rewrite the historical
Q31/Q61 v2 registration or qualify JMLE, fixed-SD, penalties, anchors or screens.
No random data are generated while constructing the plan.
"""
from dataclasses import asdict
import numpy as np

from validation.mfrm_recovery_design import RecoveryDesign, GENERATOR_VERSION, replicate_seed
from validation.mml_grouped_fit import OPTIONS, LIMITS


STUDY_ID = 'joint_mml_population_20261004_v1'
MASTER_SEED = 2**120+202610040173
REPLICATES = 2000
QUADRATURE = (641, 1281)


def design(model, assignment, persons=192, categories=4, sigma=1.2, spread=.4, levels=(3, 2, 3)):
    r, t, c = levels
    step = np.linspace(-.9, .9, categories-1) if categories > 2 else np.zeros(1)
    shape = np.linspace(-1., 1., categories-1)**2
    shape -= shape.mean()
    steps = [step] if model == 'RSM' else [step+.15*(j-(c-1)/2)*shape for j in range(c)]
    spec = RecoveryDesign(model=model, assignment=assignment, persons=persons, population_sd=sigma,
        rater=tuple(np.linspace(-.45, .45, r)), task=tuple(np.linspace(-.25, .25, t)),
        criterion=tuple(np.linspace(-.5, .7, c)), steps=tuple(tuple(row) for row in steps),
        log_slopes=tuple(np.linspace(-spread, spread, c)) if model == 'GPCM' else ())
    spec.validate()
    return spec


def scenario_id(spec):
    spread = max(abs(x) for x in spec.log_slopes) if spec.log_slopes else 0.
    return (f'{spec.model}_{spec.assignment}_N{spec.persons}_K{len(spec.steps[0])+1}'
        f'_S{round(spec.population_sd*100):03d}_A{round(spread*100):03d}'
        f'_R{len(spec.rater)}T{len(spec.task)}C{len(spec.criterion)}')


def make_plan():
    scenarios = []
    def add(spec, purpose):
        paired = spec.model in {'RSM', 'PCM'} and len(spec.steps[0]) > 1
        scenarios.append(dict(id=scenario_id(spec), cohort='population_recovery', purpose=purpose,
            design=asdict(spec), replicates=REPLICATES,
            fit_models=['RSM', 'PCM'] if paired else [spec.model],
            primary_model=spec.model, paired_rsm_pcm_comparison=paired,
            independent_R_replicates=[1, 1001], TAM_replicates=[1] if spec.model != 'GPCM' else []))
    # Core factorial: models x assignment x population sample size.
    for model in ['RSM', 'PCM', 'GPCM']:
        for assignment in ['regular', 'sparse', 'unbalanced']:
            for n in [48, 192]: add(design(model, assignment, persons=n), 'core')
    # Category support, slope spread and scale are explicit additional cells;
    # this is not a claim to have crossed every factor with every other factor.
    for model in ['RSM', 'PCM', 'GPCM']:
        for assignment in ['regular', 'sparse', 'unbalanced']:
            for k in [2, 6]: add(design(model, assignment, categories=k), 'category_grid')
    for assignment in ['regular', 'sparse', 'unbalanced']:
        for spread in [0., .8]: add(design('GPCM', assignment, spread=spread), 'slope_grid')
    for model in ['RSM', 'PCM', 'GPCM']:
        for sigma in [.6, 2.]: add(design(model, 'regular', sigma=sigma), 'population_scale')
        add(design(model, 'sparse', levels=(6, 4, 4)), 'larger_facet_design')
    controls = []
    for model in ['RSM', 'PCM', 'GPCM']:
        for mutation, expected in [('all_minimum', 'input_single_observed_category'),
            ('all_maximum', 'input_single_observed_category'), ('unused_internal_category', 'prefit_no_finite_mle'),
            ('missing_declared_rater', 'input_missing_declared_level'),
            ('fractional_weight', 'input_nonunit_weight'), ('aliased_rater_task', 'prefit_structural_rank_deficient')]:
            controls.append(dict(id=f'NEG_{model}_{mutation}', cohort='negative_control', mutation=mutation,
                design=asdict(design(model, 'regular', persons=48)), replicates=1, fit_models=[model],
                primary_model=model, expected_status=expected,
                independent_R_replicates=[], TAM_replicates=[], paired_rsm_pcm_comparison=False))
    assert len(scenarios) == 51 and len({s['id'] for s in scenarios}) == 51
    return dict(schema='joint-mml-population-study-v1', study_id=STUDY_ID,
        classification='PROSPECTIVE_PROTOCOL_PENDING_EXTERNAL_RECEIPT', scientific_inference_ready=False,
        master_seed=MASTER_SEED, generator=GENERATOR_VERSION, scenarios=scenarios+controls,
        population_replicates_per_cell=REPLICATES, population_datasets=51*REPLICATES, negative_controls=18,
        seed_rule='mfrm_recovery_design.replicate_seed(master_seed, study_id+":"+scenario_id, replicate); full 128-bit integer',
        ordered_runs='scenario order, then replicate 1..replicates; execution order cannot change inputs',
        scope=dict(estimator='native identified unpenalized unit-weight normal-mean-zero joint free-SD MML',
            slopes='GPCM: Criterion-owned, geometric mean one, all expanded log-slopes bounded [-3,3]',
            sigma_bounds=[.05, 10.], start='native input-only cold start; sigma=1; no DGP truth used in initialization',
            not_qualified_by_this_study=['JMLE','fixed-SD MML','anchors','penalties','nonunit weights',
                'latent regression','GPCM LRT','bias screens','causal/fairness claims','public app integration']),
        estimator=dict(estimator_id='PYTHON_MML_JOINT_GH641_GROUPED_V1', quadrature=list(QUADRATURE), quadrature_implementation='scipy.special.roots_hermitenorm; zero weights contribute exactly zero',
            primary_q=QUADRATURE[0], check_q=QUADRATURE[1], independent_cold_starts=True,
            stages=list(OPTIONS), numerical_limits=LIMITS,
            quadrature_limits=dict(coordinate=1e-4, cross_score_supnorm=1e-4, nll_per_person=1e-6,
                posterior_moment=1e-4, covariance_scaled=.005, expanded_se_relative=.005),
            all_stages_preserved=True, choose_best_run=False, retry=False),
        numerical_reference=dict(independent_R='mml_joint_information_reference.R with independently constructed R cell counts and a rowwise same-point NLL check, at both primary/check returned points for planned subsets',
            finite_R_limits=dict(nll_per_person=1e-9, gradient_per_person=5e-7, information_per_person=2e-6,
                hessian_step_per_person=2e-6, covariance_scaled=.005, posterior_moment=1e-8,
                expanded_se_relative=.005, expansion=1e-8),
            continuous='Independent Python full information and R continuous moments at each subset primary point',
            true_calibration='At the planned independent-R subset, also check continuous moments at DGP truth for the primary model; every run checks true-calibration moments under both finite rules',
            continuous_settings=[dict(bound=12., rel_tol=1e-9, abs_tol=1e-11),dict(bound=16., rel_tol=1e-11, abs_tol=1e-13)],
            continuous_limits=dict(nll_per_person=1e-6, score_supnorm=1e-4, posterior_moment=1e-4, expanded_se_relative=.005,
                reference_nll=1e-8, reference_gradient=1e-8, reference_information_per_person=2e-8,
                reference_moments=1e-9, relative_mass_error=1e-7, tail_mass=1e-10, tail_gradient=1e-8,
                tail_information=1e-7, numeric_information_allowance_per_person=1e-6,
                r_nll=1e-8, r_moments=1e-9, r_sigma_score=1e-8),
            native_same_point='All independent-R subset points also use native rowwise app NLL/gradient on identical positive-weight nodes',
            native_limits=dict(nll=1e-8, gradient=1e-8),
            TAM=dict(version='4.3.25', models=['RSM','PCM'], grids=[dict(bound=12.,spacing=.1),dict(bound=16.,spacing=.05)],
                maxiter=2000, conv=1e-10, convD=1e-12, convM=1e-12, Msteps=20,
                limits=dict(coordinate=1e-5, nll_per_person=1e-8, joint_score=1e-4, expanded_se_relative=.005,
                    design=1e-14, fresh_moments=1e-8),
                se_definition='Native TAM SE retained separately; qualification uses the matching full joint inverse information')),
        endpoints=dict(
            calibration='Native expanded facets, steps, log-slopes, slopes, log-population-SD and population-SD; fixed coordinates excluded',
            alignment='None: retain Criterion common location; sum-zero block mean bias is a constraint identity, use level-specific bias',
            intervals='95% full joint observed-information Wald; positive slope/SD via exponentiated log-Wald interval',
            primary='Level-specific bias, RMSE, empirical SD, sqrt(mean(model SE squared))/empirical SD, interval coverage and availability',
            covariance='Never floor eigenvalues or substitute an SD-fixed block; ineligible fits have unavailable SE/CI',
            person='EAP bias/RMSE and existing normal posterior interval, with fitted calibration and known-true calibration reported separately',
            person_quantiles='Fixed N(0,sigma^2) quartile cutpoints; no empirical ranks or sample mean alignment',
            person_interpretation='No unbiasedness requirement for EAP within true-theta bins, and no interpretation of posterior SD as sampling SE',
            optimizer_truth_diagnostic='Report fitted NLL minus NLL at the true DGP point. A worse value is evidence against maximization; never use truth to restart/select fits or to remove operationally eligible estimates from recovery summaries.',
            comparison='RSM versus PCM on the same data, both numerically eligible; LR=2*(NLL_RSM-NLL_PCM); df=(Criterion_levels-1)*(categories-2)',
            comparison_status='No LRT for binary equal models; failure to fit either model means test unavailable, not nonrejection',
            comparison_metrics='Null rejection at alpha=.05 under RSM; power under PCM; AIC/BIC and K descriptive; both available-test and all-planned denominators'),
        interpretation_limits=dict(calibration_absolute_bias=.1, positive_parameter_relative_bias=.1,
            rms_model_se_over_empirical_sd=[.9,1.1], coverage_95=[.925,.975], eligible_failure_rate_max=.02,
            null_LRT_rejection=[.035,.065],
            rationale='0.1 on native logit/log scale (or relative positive scale), 10% SE error, 2.5 percentage points from nominal coverage and 1.5 points from nominal type-I rate; engineering acceptance targets, not universal psychometric standards',
            decision='Pass only when a pointwise 95% Monte Carlo interval lies inside the prespecified acceptable region; outside is fail, overlap is inconclusive. Pointwise evidence cannot be advertised as simultaneous qualification over selected cells.',
            publication='No automatic promotion. Keep output-specific uncertainty and unsupported conditions; all cells/endpoints reported, including inconclusive and failed ones.'),
        monte_carlo=dict(independent_unit='replicate, not response/Person/parameter row', replicates=REPLICATES,
            expected_mcse_at_95=float(np.sqrt(.95*.05/REPLICATES)), worst_binary_mcse=float(np.sqrt(.25/REPLICATES)),
            rationale='At 95% coverage, 2000 independent replicates yield MCSE about 0.00487 before unavailable intervals; report realized conditional and failure-aware precision.',
            interval_policy='Exact binomial 95% intervals for one binary endpoint per replicate; replicate-cluster delta-method intervals for population Person ratios; robust joint-moment delta method for SE ratios',
            failure_policy='No replacement, no numerical retry or optional extension. Every planned RunId remains; known targets with missing results and unknown truth from generation failure are distinguished.',
            interim_analysis='Execution progress only within a scenario; endpoint summaries only after all its planned replicates are terminal',
            source='Morris, White & Crowther (2019), Zotero PKQMUBH7; https://arxiv.org/html/1712.03198v3'),
        execution=dict(input_order='Plan+source+runtime receipt before any fresh data; retain and hash each raw input and truth before its fit authorization',
            attempts='Append-only SQLite events; unique RunId/phase; started execution without completion is a retained interrupted attempt, never replayed',
            restart='Only unstarted registered runs may be dispatched after an operational interruption; sealed inputs/finished fits are never replaced',
            storage='Compressed retained CSV, truth and numerical arrays; posterior grids reconstructable from input, final coordinates and frozen quadrature',
            disk_reserve_bytes=4*1024**3, maximum_workers=2,
            source_policy='Run an isolated frozen source copy; verify source/runtime before generation and before fitting, retain hashes and executable identities',
            timestamp='RFC3161 SHA256+nonce signed under pinned FreeTSA root, both query and exact plan-manifest bytes verified before execution'))


def planned_runs(plan):
    for scenario in plan['scenarios']:
        for replicate in range(1, scenario['replicates']+1):
            yield dict(run_id=f"{scenario['id']}__{replicate:04d}", scenario_id=scenario['id'],
                replicate=replicate, seed=replicate_seed(plan['master_seed'], plan['study_id']+':'+scenario['id'], replicate),
                cohort=scenario['cohort'])
