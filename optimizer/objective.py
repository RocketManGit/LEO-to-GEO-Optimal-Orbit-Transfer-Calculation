# optimizer/objective.py
# ─────────────────────────────────────────────────────────────────────────────
# Scalarized objective function for the Hohmann transfer optimizer.
#
# The optimization problem:
#   minimize:  J(r_apo) = w * (Δv / Δv_ref) + (1-w) * (TOF / TOF_ref)
#   over:      r_apo ∈ [r2, r_apo_max]
#
# Normalizing by reference values (Hohmann Δv and TOF) ensures both terms
# are dimensionless and on the same scale — so the weight w is meaningful.
# ─────────────────────────────────────────────────────────────────────────────

import numpy as np
from physics.maneuvers import fast_transfer_physics


def build_objective(r1_km: float, r2_km: float,
                    dv_ref: float, tof_ref: float,
                    weight_dv: float):
    """
    Build and return the scalarized objective function for a given weight.

    This is a closure — it captures the mission parameters and weight, and
    returns a callable f(x) that scipy.optimize.minimize can call directly.

    Args:
        r1_km     : Departure orbit radius [km]
        r2_km     : Target   orbit radius  [km]
        dv_ref    : Normalisation reference for Δv  (use Hohmann Δv)  [km/s]
        tof_ref   : Normalisation reference for TOF (use Hohmann TOF) [s]
        weight_dv : w ∈ [0, 1] — weight placed on Δv minimization
                    (1 - weight_dv goes to TOF minimization)

    Returns:
        callable: f(x) where x = [r_apo_km]
                  Returns a large penalty (1e9) for physically invalid inputs.
    """
    if not (0.0 <= weight_dv <= 1.0):
        raise ValueError(f"weight_dv must be in [0, 1], got {weight_dv}.")

    def objective(x: list) -> float:
        r_apo = x[0]
        dv, tof = fast_transfer_physics(r_apo, r1_km, r2_km)

        # Return a large penalty if the geometry is invalid
        if not np.isfinite(dv) or not np.isfinite(tof):
            return 1e9

        # Normalised scalarized cost — both terms are dimensionless
        J = weight_dv * (dv / dv_ref) + (1 - weight_dv) * (tof / tof_ref)
        return J

    return objective
