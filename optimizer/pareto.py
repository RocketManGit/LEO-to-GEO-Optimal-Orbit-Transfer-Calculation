# optimizer/pareto.py
# ─────────────────────────────────────────────────────────────────────────────
# Pareto front construction via scalarization.
#
# Strategy: sweep weight w from 0 → 1, solve a fresh scipy.optimize problem
# at each weight. The collection of solutions traces the Pareto front.
#
#   w = 1.0  →  pure Δv minimizer  →  recovers the Hohmann solution
#   w = 0.0  →  pure TOF minimizer →  fastest physically achievable transfer
#   w = 0.5  →  balanced trade-off
# ─────────────────────────────────────────────────────────────────────────────

import numpy as np
from scipy.optimize import minimize, Bounds

from physics.maneuvers       import fast_transfer_physics
from physics.orbital_mechanics import altitude_to_radius
from optimizer.objective     import build_objective


def solve_pareto_front(alt1_km: float, alt2_km: float,
                       n_points: int   = 60,
                       r_apo_max_factor: float = 8.0) -> dict:
    """
    Trace the Pareto front by solving a constrained L-BFGS-B optimization
    problem for each weight w in [0, 1].

    Each individual solve finds the r_apo (transfer ellipse apoapsis) that
    minimises the scalarized cost for that particular weight. Running this
    across many weights maps out the full trade space between Δv and TOF.

    Args:
        alt1_km          : Departure orbit altitude above surface [km]
        alt2_km          : Target   orbit altitude above surface  [km]
        n_points         : Number of Pareto points (= number of optimizer runs)
        r_apo_max_factor : Upper bound on r_apo = r2 * this factor

    Returns:
        dict with keys:
            weights       — array of w values used
            dv_km_s       — optimal Δv at each weight [km/s]
            tof_hours     — optimal TOF at each weight [hours]
            r_apo_km      — optimal apoapsis at each weight [km]
            hohmann_dv    — Hohmann reference Δv  [km/s]
            hohmann_tof   — Hohmann reference TOF [hours]
            r1_km         — departure orbit radius [km]
            r2_km         — target   orbit radius  [km]

    Raises:
        ValueError: if alt1_km >= alt2_km.
    """
    r1_km = altitude_to_radius(alt1_km)
    r2_km = altitude_to_radius(alt2_km)

    if r1_km >= r2_km:
        raise ValueError("alt1_km must be less than alt2_km.")

    # ── Reference point: Hohmann (r_apo == r2) ───────────────────────────────
    dv_ref, tof_ref   = fast_transfer_physics(r2_km, r1_km, r2_km)
    hohmann_tof_hours = tof_ref / 3600.0

    print(f"\n  Hohmann reference  →  Δv = {dv_ref:.4f} km/s  |  "
          f"TOF = {hohmann_tof_hours:.3f} hrs")
    print(f"  Solving {n_points} optimizer problems across weight range...\n")

    # ── Bounds: r_apo ∈ [r2, r2 * factor] ───────────────────────────────────
    r_apo_max = r2_km * r_apo_max_factor
    bounds    = Bounds(lb=r2_km, ub=r_apo_max)

    # ── Sweep weights and solve ───────────────────────────────────────────────
    # Sweep from w=1 (Hohmann) down to w=0 (time-optimal) so each solution
    # can warm-start from the previous one — much more stable than always
    # starting from r2.
    weights      = np.linspace(1.0, 0.0, n_points)   # high→low weight
    dv_results   = []
    tof_results  = []
    rapo_results = []

    # Initial guess at w=1: Hohmann is optimal, start there
    x0_current = r2_km

    for i, w in enumerate(weights):

        obj = build_objective(r1_km, r2_km, dv_ref, tof_ref, weight_dv=w)

        # Warm-start: use the previous solution as the initial guess.
        # As w decreases, the optimizer nudges r_apo progressively higher
        # rather than having to jump from r2 all the way to the optimum.
        result = minimize(
            obj,
            x0      = [x0_current],
            method  = 'L-BFGS-B',
            bounds  = bounds,
            options = {'ftol': 1e-14, 'gtol': 1e-12, 'maxiter': 2000}
        )

        # Update warm-start for next iteration
        x0_current = result.x[0]

        r_apo_opt       = result.x[0]
        dv_opt, tof_opt = fast_transfer_physics(r_apo_opt, r1_km, r2_km)

        dv_results.append(dv_opt)
        tof_results.append(tof_opt / 3600.0)
        rapo_results.append(r_apo_opt)

        if (i + 1) % 10 == 0:
            print(f"  [{i+1:3d}/{n_points}]  w={w:.2f}  →  "
                  f"Δv={dv_opt:.4f} km/s  |  TOF={tof_opt/3600:.3f} hrs  |  "
                  f"r_apo={r_apo_opt:.0f} km")

    print(f"\n  Done. {n_points} Pareto-optimal solutions found.")

    # Reverse so arrays run w=0→1 (time-optimal → Δv-optimal) for clean plotting
    weights      = weights[::-1]
    dv_results   = list(reversed(dv_results))
    tof_results  = list(reversed(tof_results))
    rapo_results = list(reversed(rapo_results))

    return {
        'weights'     : weights,
        'dv_km_s'     : np.array(dv_results),
        'tof_hours'   : np.array(tof_results),
        'r_apo_km'    : np.array(rapo_results),
        'hohmann_dv'  : dv_ref,
        'hohmann_tof' : hohmann_tof_hours,
        'r1_km'       : r1_km,
        'r2_km'       : r2_km,
    }
