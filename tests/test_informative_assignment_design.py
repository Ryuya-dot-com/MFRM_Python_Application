from __future__ import annotations

import numpy as np

from mfrm_app.operating_characteristics import deterministic_replicate_seed
from validation.estimand_distribution_study import (
    BASE_SEED,
    _latent_rows,
    generate_person_shapes,
)
from validation.facets_pcm_known_truth_smoke import (
    CRITERION_TRUTH,
    N_PERSONS,
    RATER_TRUTH,
    TASK_TRUTH,
    THRESHOLD_CONDITIONS,
    apply_threshold_condition,
)
from validation.informative_assignment_design import (
    ALIGNED_PAIRS,
    DESIGNS,
    ability_quartiles,
    apply_informative_assignment_design,
    informative_assignment_audit,
)


def _complete_replicate(replicate: int = 221):
    seed = deterministic_replicate_seed(BASE_SEED, "person-shape-study", replicate)
    persons = generate_person_shapes(seed)["normal"]
    uniform_seed = deterministic_replicate_seed(BASE_SEED, "row-uniforms", replicate)
    rows = N_PERSONS * len(RATER_TRUTH) * len(TASK_TRUTH) * len(CRITERION_TRUTH)
    latent = _latent_rows(persons, np.random.default_rng(uniform_seed).random(rows))
    complete = apply_threshold_condition(latent, THRESHOLD_CONDITIONS["heterogeneous"])
    theta = latent[["Person", "Theta"]].drop_duplicates("Person")
    return complete.merge(theta, on="Person", how="left", validate="many_to_one")


def test_ability_quartiles_are_stable_and_exactly_balanced() -> None:
    quartiles = ability_quartiles(_complete_replicate())
    assert len(quartiles) == 80
    assert quartiles.groupby("AbilityQuartile").size().to_dict() == {
        1: 20,
        2: 20,
        3: 20,
        4: 20,
    }
    assert quartiles["Theta"].is_monotonic_increasing
    assert quartiles["AbilityRank"].tolist() == list(range(1, 81))


def test_aligned_assignment_has_registered_pair_per_quartile() -> None:
    complete = _complete_replicate()
    aligned = apply_informative_assignment_design(
        complete, "ability_severity_aligned_connected"
    )
    quartiles = ability_quartiles(complete)
    pairs = (
        aligned.groupby("Person")["Rater"]
        .agg(lambda values: tuple(sorted(set(values))))
        .rename("Pair")
        .reset_index()
        .merge(quartiles[["Person", "AbilityQuartile"]], on="Person")
    )
    for quartile, expected in enumerate(ALIGNED_PAIRS, start=1):
        assert set(pairs.loc[pairs["AbilityQuartile"].eq(quartile), "Pair"]) == {expected}


def test_sparse_designs_match_density_exposure_and_identifiability() -> None:
    complete = _complete_replicate()
    audits = {design: informative_assignment_audit(complete, design) for design in DESIGNS}
    assert audits["complete"]["Rows"] == 1920
    for design in ("planned_connected", "ability_severity_aligned_connected"):
        audit = audits[design]
        assert audit["Rows"] == 960
        assert audit["RaterPersonCounts"] == {level: 40 for level in RATER_TRUTH}
        assert audit["RaterRowCounts"] == {level: 240 for level in RATER_TRUTH}
        assert audit["PersonRaterComponents"] == 1
        assert audit["ConstrainedPCMNullity"] == 0
        assert audit["AllInvariantsPass"]
    assert audits["ability_severity_aligned_connected"][
        "ThetaAssignedSeveritySpearman"
    ] > 0.9


def test_all_designs_are_exact_subsets_of_same_complete_scores() -> None:
    complete = _complete_replicate()
    keys = ["Person", "Rater", "Task", "Criterion"]
    complete_scores = complete.set_index(keys)["Score"]
    for design in DESIGNS:
        observed = apply_informative_assignment_design(complete, design).set_index(keys)["Score"]
        assert observed.equals(complete_scores.loc[observed.index])
