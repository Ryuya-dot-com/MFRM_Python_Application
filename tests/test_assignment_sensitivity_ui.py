from __future__ import annotations

from streamlit.testing.v1 import AppTest


HARNESS = r'''
import pandas as pd
import streamlit_app as app

assignment = {
    "p0": ("r2", "r3"), "p1": ("r2", "r3"),
    "p2": ("r1", "r3"), "p3": ("r1", "r2"),
    "p4": ("r0", "r3"), "p5": ("r0", "r2"),
    "p6": ("r0", "r1"), "p7": ("r0", "r1"),
}
rows = []
for person, raters in assignment.items():
    for rater in raters:
        for task in ("t1", "t2"):
            rows.append({
                "Person": person,
                "Rater": rater,
                "Task": task,
                "Score": (int(person[1:]) + int(rater[1:])) % 4,
                "Weight": 1.0,
            })
data = pd.DataFrame(rows)
result = app.mfrm_estimate(
    data,
    person_col="Person",
    facet_cols=["Rater", "Task"],
    score_col="Score",
    rating_min=0,
    rating_max=3,
    weight_col="Weight",
    keep_original=True,
    model="RSM",
    method="JMLE",
    maxit=120,
    reltol=1e-5,
)
app.render_fixed_density_assignment_sensitivity(result)
'''


MILP_HARNESS = r'''
import numpy as np
import pandas as pd
import streamlit_app as app

base = [
    (("r0", "A"), ("r1", "B")), (("r0", "A"), ("r2", "B")),
    (("r1", "A"), ("r0", "B")), (("r1", "A"), ("r2", "B")),
    (("r2", "A"), ("r0", "B")), (("r2", "A"), ("r1", "B")),
]
rows = []
for person_index in range(12):
    for rater, signature in base[person_index % 6]:
        tasks = ["t1"] if signature == "A" else ["t1", "t2"]
        for task in tasks:
            for criterion_index in range(5):
                latent = (
                    1.5 + 0.22 * (person_index - 5.5) - 0.25 * int(rater[1])
                    + 0.35 * (criterion_index - 2) + (0.15 if task == "t2" else -0.15)
                )
                rows.append({
                    "Person": f"p{person_index:02d}", "Rater": rater,
                    "Task": task, "Criterion": f"c{criterion_index}",
                    "Score": int(np.clip(np.rint(latent), 0, 3)), "Weight": 1.0,
                })
result = app.mfrm_estimate(
    pd.DataFrame(rows), person_col="Person",
    facet_cols=["Rater", "Task", "Criterion"], score_col="Score",
    rating_min=0, rating_max=3, weight_col="Weight", keep_original=True,
    model="RSM", method="JMLE", maxit=250, reltol=1e-5,
)
app.render_fixed_density_assignment_sensitivity(result)
'''


MML_HARNESS = HARNESS.replace('method="JMLE"', 'method="MML"').replace(
    "maxit=120", "maxit=300"
)


def test_fixed_density_runner_renders_and_executes_without_ui_exception():
    at = AppTest.from_string(HARNESS, default_timeout=45).run()

    assert not at.exception
    at.button(key="run_assignment_sensitivity").click()
    at.run(timeout=45)

    assert not at.exception
    bundle = at.session_state["mfrm_assignment_sensitivity_bundle"]
    assert bundle["available"] is True
    assert bundle["completion"].iloc[0]["PairedContrastsComplete"] == 2
    assert at.download_button(key="dl_assignment_sensitivity_bundle")


def test_unequal_context_ui_selects_endpoint_only_milp_fallback():
    at = AppTest.from_string(MILP_HARNESS, default_timeout=60).run()

    assert not at.exception
    assert not [item for item in at.multiselect if item.key == "assignment_sensitivity_doses"]
    at.number_input(key="assignment_sensitivity_replicates").set_value(1)
    at.button(key="run_assignment_sensitivity").click()
    at.run(timeout=60)

    assert not at.exception
    bundle = at.session_state["mfrm_assignment_sensitivity_bundle"]
    assert bundle["available"] is True
    assert bundle["design_engine"] == "context_margin_milp"
    assert bundle["dose_table"]["AchievedAlignmentDose"].tolist() == [0.0, 1.0]
    assert bundle["solver_audit"].iloc[0]["GlobalOptimumCertified"]
    assert at.download_button(key="dl_assignment_sensitivity_bundle")


def test_mml_ui_runs_rank_preserving_population_generator():
    at = AppTest.from_string(MML_HARNESS, default_timeout=60).run()

    assert not at.exception
    at.selectbox(key="assignment_sensitivity_person_generator").select(
        "mml_population_rank_preserving"
    )
    at.number_input(key="assignment_sensitivity_replicates").set_value(1)
    at.run(timeout=60)
    at.button(key="run_assignment_sensitivity").click()
    at.run(timeout=60)

    assert not at.exception
    bundle = at.session_state["mfrm_assignment_sensitivity_bundle"]
    assert bundle["available"] is True
    assert bundle["person_generation_mode"] == "mml_population_rank_preserving"
    assert bundle["person_generation_summary"].iloc[0]["RankPreservedExactly"]
