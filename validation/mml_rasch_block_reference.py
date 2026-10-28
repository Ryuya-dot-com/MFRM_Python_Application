"""Additive continuous reference, with ten persons per adaptive integral call.

The independent RaschIntegral kernel, global parameterization and tolerances
are unchanged. This helper neither fits parameters nor computes covariances.
"""
import numpy as np

from validation.mml_rasch_integral import RaschIntegral


class BlockRaschIntegral(RaschIntegral):
    """Preserve global coordinates and person order while limiting vector width."""

    def __init__(self, data, config):
        super().__init__(data, config)
        self.data, self.config = data, config

    def continuous(self, point, *, bound=12., rel_tol=1e-10, abs_tol=1e-12, information=False):
        values, row_work = [], 0
        for start in range(0, self.n, 10):
            rows = self.data.iloc[np.flatnonzero((self.person >= start) & (self.person < start + 10))]
            value = RaschIntegral(rows, self.config).continuous(point, bound=bound,
                rel_tol=rel_tol, abs_tol=abs_tol, information=information)
            values.append(value)
            row_work += value["quadrature_neval"] * len(rows)
        combined = {}
        for key in values[0]:
            if key in {"eap", "sd"}:
                combined[key] = np.concatenate([value[key] for value in values])
            elif key == "sigma":
                combined[key] = values[0][key]
            elif key == "quadrature_status":
                combined[key] = max(value[key] for value in values)
            else:
                # Sum full information before any downstream fixed-SD inverse.
                # Summing blockwise maximum-norm allowances is conservative.
                combined[key] = sum(value[key] for value in values)
        combined.update(person_labels=self.person_labels, integration_blocks=len(values),
                        integration_person_block_limit=10, integrand_response_row_evaluations=row_work)
        return combined
