"""Observation-level design support for a prespecified rater severity contrast."""
from collections import defaultdict, deque

import numpy as np
import pandas as pd

from mfrm_app import basic_analysis, engine


def _pair(values, levels, name):
    if not isinstance(values, (tuple, list)) or len(values) != 2:
        raise ValueError(f"{name} requires an ordered pair of two levels")
    if any(not isinstance(x, str) or x.strip() not in levels for x in values):
        raise ValueError(f"{name} must refer to two observed level labels")
    left, right = (x.strip() for x in values)
    if left == right:
        raise ValueError(f"{name} requires two distinct levels")
    return left, right


def _observed_assignments(data, person, facets, score, rating_min, rating_max, *, allow_empty=False):
    cols = [person, *facets, score]
    if len(set(cols)) != len(cols) or any(f in {"Person", "Score", "Weight"} for f in facets):
        raise ValueError("Select distinct roles and use facet names other than Person, Score and Weight")
    if not data.columns.is_unique or any(col not in data for col in cols):
        raise ValueError("Selected columns must exist and column names must be unique")
    audit = engine.build_response_data_audit(data, person_col=person, facet_cols=facets, score_col=score)
    engine.validate_rating_scores(audit["rows"])
    keep = audit["rows"]["IncludedInLikelihood"].to_numpy(dtype=bool)
    if rating_min is not None or rating_max is not None:
        if not (basic_analysis._is_int(rating_min) and basic_analysis._is_int(rating_max) and rating_min < rating_max):
            raise ValueError("Specify both integer rubric endpoints with rating_min < rating_max")
        if not pd.to_numeric(data.loc[keep, score]).round().between(rating_min, rating_max).all():
            raise ValueError("Observed scores fall outside the supplied rating range")
    rows = data.loc[keep, [person, *facets]].rename(columns={person: "Person"}).copy()
    if rows.empty and not allow_empty:
        raise ValueError("No observed ratings remain after input exclusions")
    for col in rows:
        rows[col] = pd.Categorical(engine.normalized_rating_ids(rows[col]))
    return rows.drop_duplicates(), len(rows)


def _overlap(rows, rater, pair, group_level):
    # Traverse the bipartite graph so one commonly rated work costs O(raters) edges.
    by_rater, by_person = defaultdict(set), defaultdict(set)
    for person_id, rater_id in rows.loc[:, ["Person", rater]].drop_duplicates().itertuples(index=False, name=None):
        by_person[person_id].add(rater_id)
        by_rater[rater_id].add(person_id)
    left, right = pair
    queue = deque([left])
    parent = {left: None}
    seen_people = set()
    while queue and right not in parent:
        current = queue.popleft()
        for person_id in sorted(by_rater[current] - seen_people):
            seen_people.add(person_id)
            for other in sorted(by_person[person_id]):
                if other not in parent:
                    parent[other] = current
                    queue.append(other)
    path = []
    if right in parent:
        node = right
        while node != left:
            previous = parent[node]
            path.append(dict(from_rater=previous, to_rater=node,
                             persons=sorted(by_rater[previous] & by_rater[node])))
            node = previous
        path.reverse()
    kind = "none" if not path else ("direct" if len(path) == 1 else "indirect")
    return dict(group_level=group_level, kind=kind, path=path)


def _in_row_space(matrix, target):
    # Economy SVD also covers wide designs whose omitted right singular vectors
    # span part of the null space. Project onto the retained row space.
    basis, tol = _row_space(matrix)
    residual = target - (target @ basis.T) @ basis
    return bool(np.linalg.norm(residual) <= tol * np.linalg.norm(target))


def _row_space(matrix):
    _, singular, vh = np.linalg.svd(matrix, full_matrices=False)
    tol = max(matrix.shape) * np.finfo(float).eps
    return (vh[singular > tol * singular[0]] if len(singular) else vh), tol


def _constraint_rows(levels, offsets, origin, fixed, means):
    indices = {f: {name: i + offsets[f] for i, name in enumerate(names)} for f, names in levels.items()}
    rows = [[indices[row["Facet"]][row["Level"]]] for row in fixed]
    by_mean = defaultdict(list)
    for row in means:
        by_mean[row["Facet"], row["Group"]].append(indices[row["Facet"]][row["Level"]])
    rows.extend(by_mean.values())
    anchored = {row["Facet"] for row in fixed + means}
    for f, names in levels.items():
        if f != origin and f not in anchored:
            rows.append(range(offsets[f], offsets[f] + len(names)))
    return rows


def _design(rows, levels, signs, rater, group, weights, origin=None, fixed=(), means=()):
    offsets = {}
    n_main = 0
    for facet, names in levels.items():
        offsets[facet] = n_main
        n_main += len(names)
    n_raters = len(levels[rater])
    n_groups = len(levels[group]) if group is not None else 0
    dimension = n_main + n_raters * n_groups
    constraints = []
    if fixed or means:
        constraints = _constraint_rows(levels, offsets, origin, fixed, means)
        if group is not None:
            constraints.extend(range(n_main+i, dimension, n_raters) for i in range(n_raters))
            for j in range(n_groups):
                first = n_main + j*n_raters
                constraints.append(range(first, first+n_raters))
    if (len(rows) + len(constraints)) * dimension > 5_000_000:
        return None, None
    matrix = np.zeros((len(rows) + len(constraints), dimension))
    positions = np.arange(len(rows))
    maps = {f: {level: i for i, level in enumerate(names)} for f, names in levels.items()}
    codes = {f: rows[f].map(index).to_numpy(dtype=int) for f, index in maps.items()}
    for f, start in offsets.items():
        matrix[positions, start+codes[f]] = signs[f]
    if group is not None:
        matrix[positions, n_main+codes[rater]+n_raters*codes[group]] = 1
    target = np.zeros(dimension)
    for term in weights:
        i = maps[rater][term["rater_level"]]
        j = maps[group][term["group_level"]] if group is not None else None
        col = n_main+j*n_raters+i if group is not None else offsets[rater]+i
        target[col] = term["weight"]
    for i, cols in enumerate(constraints, len(rows)):
        matrix[i, cols] = 1
    return matrix, target


def _planned_rows(data, person, facets):
    cols = [person, *facets]
    if (not isinstance(data, pd.DataFrame) or data.empty or not data.columns.is_unique
            or any(col not in data for col in cols)):
        raise ValueError("planned_assignments requires a nonempty table with the work and all selected facet columns")
    if any(engine._blank_or_na_mask(data[col]).any() for col in cols):
        raise ValueError("Planned assignments require complete work and facet IDs")
    rows = data[cols].rename(columns={person: "Person"})
    return rows.apply(engine.normalized_rating_ids).drop_duplicates().reset_index(drop=True)


def _restoring_rows(observed, missing, target):
    basis, tol = _row_space(observed)
    residual = target - (target @ basis.T) @ basis
    outside = missing - (missing @ basis.T) @ basis
    # One additional row restores support exactly when its residual spans the missing direction.
    lengths = np.linalg.norm(outside, axis=1)
    unit = np.divide(outside, lengths[:, None], out=np.zeros_like(outside),
                     where=lengths[:, None] > tol*np.linalg.norm(missing, axis=1)[:, None])
    remainder = residual - (unit @ residual)[:, None]*unit
    return (np.linalg.norm(remainder, axis=1) <= tol*np.linalg.norm(target)).tolist()


def _compare_plan(observed, planned, levels, signs, rater, group, weights, person, observed_support):
    keys = list(planned.columns)
    roster = pd.MultiIndex.from_frame(planned)
    scored = pd.MultiIndex.from_frame(observed)
    missing = planned.loc[~roster.isin(scored)].sort_values(keys)
    extra = observed.loc[~scored.isin(roster)].sort_values(keys)
    planned_support = observed_support
    if not missing.empty or not extra.empty:
        matrix, target = _design(planned, levels, signs, rater, group, weights)
        planned_support = None if matrix is None else _in_row_space(matrix, target)
        del matrix
    completed_support = observed_support
    repair = [None]*len(missing)
    if observed_support is not True and not missing.empty:
        completed = pd.concat([observed, missing], ignore_index=True)
        matrix, target = _design(completed, levels, signs, rater, group, weights)
        completed_support = None if matrix is None else _in_row_space(matrix, target)
        if observed_support is False and completed_support is not None:
            repair = (_restoring_rows(matrix[:len(observed)], matrix[len(observed):], target)
                      if completed_support else [False]*len(missing))
    by_group = []
    groups = list(dict.fromkeys(term["group_level"] for term in weights))
    left, right = weights[0]["rater_level"], weights[1]["rater_level"]
    for code in groups:
        rows = [part if group is None else part.loc[part[group].eq(code)]
                for part in (planned, observed, missing, extra)]
        plan_works, scored_works = (set(part.loc[part[rater].eq(left), "Person"])
                                 & set(part.loc[part[rater].eq(right), "Person"]) for part in rows[:2])
        by_group.append(dict(group_level=code, planned_cells=len(rows[0]), observed_cells=len(rows[1]),
                             missing_cells=len(rows[2]), unplanned_cells=len(rows[3]),
                             planned_shared_works=sorted(plan_works), observed_shared_works=sorted(scored_works),
                             missing_shared_works=sorted(plan_works - scored_works)))
    return dict(planned_cells=len(planned), observed_planned_cells=len(planned)-len(missing),
                planned_observation_estimable=planned_support,
                completed_observation_estimable=completed_support,
                missing=[dict(assignment=row, restores_contrast=restores) for row, restores in
                         zip(missing.rename(columns={"Person": person}).to_dict("records"), repair, strict=True)],
                unplanned=extra.rename(columns={"Person": person}).to_dict("records"), by_group=by_group)


def check_rater_contrast(data, *, person, facets, score, rater, difference,
                         group=None, group_difference=None, rating_min=None, rating_max=None, positive_facets=None,
                         noncenter_facet="Person", anchors=None, group_anchors=None, planned_assignments=None):
    """Trace shared persons/works and check structural support without fitting.

    `difference=(a, b)` always means severity of a minus severity of b.
    With `group` and `group_difference=(g, h)`, compare that difference in
    g minus h using a rater-by-group interaction. Include every modeled facet
    in `facets`. Use a work ID as `person` to trace repeated scoring of works.
    Supply `planned_assignments` with work and facet columns to compare the
    declared roster with the scored cells, without adding any responses.
    See the package guide for the meanings of the returned support states.
    """
    if not isinstance(data, pd.DataFrame) or (data.empty and planned_assignments is None):
        raise ValueError("A nonempty long-form ratings table is required")
    if not isinstance(facets, (list, tuple)) or rater not in facets:
        raise ValueError("rater must be one of the selected facets")
    if (group is None) != (group_difference is None):
        raise ValueError("Supply group and group_difference together")
    if group is not None and (group not in facets or group == rater):
        raise ValueError("group must be a selected facet distinct from rater")
    facets = list(facets)
    rows, included = _observed_assignments(data, person, facets, score, rating_min, rating_max,
                                          allow_empty=planned_assignments is not None)
    levels = {f: list(rows[f].cat.categories) for f in rows}
    planned = None if planned_assignments is None else _planned_rows(planned_assignments, person, facets)
    known = levels if planned is None else {f: sorted(set(levels[f]).union(planned[f])) for f in levels}
    pair = _pair(difference, known[rater], "difference")
    groups = _pair(group_difference, known[group], "group_difference") if group is not None else ()
    roles = dict(person_col=person, facet_cols=facets)
    origin, positive = basic_analysis.normalize_scale_orientation(noncenter_facet, positive_facets, **roles)
    observed = {f: set(names) for f, names in levels.items()}
    fixed = basic_analysis.normalize_anchor_records(anchors, observed_levels=observed, **roles)
    means = basic_analysis.normalize_group_anchor_records(group_anchors, observed_levels=observed,
                                                          fixed_records=fixed, **roles)
    signs = {f: 1 if f == "Person" or f in positive else -1 for f in levels}
    if groups:
        overlap = [_overlap(rows.loc[rows[group].eq(g)], rater, pair, g) for g in groups]
    else:
        overlap = [_overlap(rows, rater, pair, None)]
    # The joint interaction enters eta positively. Greater severity lowers eta.
    multiplier = -1 if groups else -signs[rater]
    weights = []
    for g, direction in zip(groups or (None,), (1, -1)):
        for level, side in zip(pair, (1, -1)):
            weights.append(dict(rater_level=level, group_level=g, weight=direction * side * multiplier))

    result = dict(status="audit_size_limit", person=person, rater=rater, group=group,
                  observation_estimable=None, fixed_by_anchors=None,
                  estimable_with_anchors=None, overlap=overlap, contrast_weights=weights,
                  facet_signs=signs, noncenter_facet=origin, anchors=fixed, group_anchors=means,
                  included_rows=included, excluded_rows=len(data) - included)

    # Planned levels only enter the observed design if a requested level has no scores.
    levels[rater] = sorted(set(levels[rater]).union(pair))
    if groups:
        levels[group] = sorted(set(levels[group]).union(groups))
    matrix, target = _design(rows, levels, signs, rater, group, weights, origin, fixed, means)
    if matrix is not None:
        obs = matrix[:len(rows)]
        data_support = _in_row_space(obs, target)
        constrained = len(matrix) > len(rows)
        fixed_target = constrained and _in_row_space(matrix[len(rows):], target)
        conditional = (data_support or fixed_target or _in_row_space(matrix, target)) if constrained else None
        if fixed_target:
            status = "fixed_by_anchors"
        elif data_support:
            status = "observed"
        elif conditional:
            status = "anchor_dependent"
        else:
            start = len(levels["Person"])
            status = "common_location_only" if _in_row_space(obs[:, start:], target[start:]) else "unsupported"
        result.update(status=status, observation_estimable=data_support, fixed_by_anchors=fixed_target,
                      estimable_with_anchors=conditional)
    if planned is not None:
        result["plan"] = _compare_plan(rows.astype(str), planned, known, signs, rater, group, weights, person,
                                       result["observation_estimable"])
    return result
