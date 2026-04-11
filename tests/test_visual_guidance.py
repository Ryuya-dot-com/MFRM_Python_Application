import streamlit_app as app


def test_visual_interpretation_checklist_contract():
    checklist = app.visual_interpretation_checklist()
    assert not checklist.empty
    assert checklist["Priority"].is_monotonic_increasing

    required = {
        "Visualization",
        "Where",
        "PrimaryQuestion",
        "ReadFirst",
        "ReviewTrigger",
        "BeginnerAction",
        "Caveat",
    }
    assert required.issubset(checklist.columns)

    joined = " ".join(checklist["Visualization"].astype(str))
    for expected in [
        "Wright map",
        "Category probability curves",
        "Fit scatter",
        "Residual PCA",
        "Bias heatmap",
        "Strict marginal",
    ]:
        assert expected in joined
