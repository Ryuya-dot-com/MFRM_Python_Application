"""Conservative score-boundary screen for free JMLE facet effects."""

import numpy as np
import pandas as pd


HOLD_REASON = "stat.jmle.facet_score_boundary_review"
HOLD_MESSAGE = (
    "A free facet level has only minimum or only maximum scores. "
    "Optimizer termination does not establish a finite unpenalized estimate; "
    "review score coverage and constraints before inference."
)


def facet_score_boundary_audit(prep, config):
    """Flag extreme free levels, respecting fixed and group-mean constraints.

    This is a screening hold, not a separation certificate: coupled constraints
    or penalties can yield finite estimates. No flag does not prove existence
    of a finite MLE or rule out separation involving several facets together.
    """
    columns = ["Facet", "Level", "ObservedRows", "MinimumScore", "MaximumScore",
               "FreeEffect", "ExtremeScorePattern", "RequiresReview"]
    if str(config.get("method", "")).upper() != "JMLE":
        return pd.DataFrame(columns=columns)
    data = prep["data"]
    if "Weight" in data:
        data = data.loc[data["Weight"].gt(0)]
    rows = []
    for facet in config["facet_names"]:
        spec = config["facet_specs"][facet]
        free = np.isnan(spec["anchors"])
        groups = [None if pd.isna(g) or g == "" else g for g in spec["groups"]]
        group_counts = {}
        for is_free, group in zip(free, groups, strict=True):
            if is_free:
                group_counts[group] = group_counts.get(group, 0) + 1
        scores = data.groupby(facet, observed=True)["Score"].agg(["size", "min", "max"])
        scores.index = scores.index.astype(str)
        for index, level in enumerate(spec["levels"]):
            group = groups[index]
            varies = bool(free[index]) and (
                group_counts.get(group, 0) > 1
                if group is not None or spec["centered"] else True
            )
            observed = scores.loc[str(level)] if str(level) in scores.index else None
            minimum = observed["min"] if observed is not None else np.nan
            maximum = observed["max"] if observed is not None else np.nan
            extreme = bool(minimum == maximum and minimum in (prep["rating_min"], prep["rating_max"]))
            rows.append(dict(Facet=facet, Level=str(level),
                ObservedRows=int(observed["size"]) if observed is not None else 0,
                MinimumScore=minimum, MaximumScore=maximum, FreeEffect=varies,
                ExtremeScorePattern=extreme, RequiresReview=varies and extreme))
    return pd.DataFrame(rows, columns=columns)
