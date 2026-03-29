# physics/maneuvers.py
# ─────────────────────────────────────────────────────────────────────────────
# Orbital maneuver calculations built on top of orbital_mechanics.py.
# Covers: Hohmann delta-v, generalised fast transfers, and time of flight.
# ─────────────────────────────────────────────────────────────────────────────

import math
import numpy as np

from .orbital_mechanics import (
    circular_velocity,
    ellipse_velocity,
    transfer_semimajor_axis,
    altitude_to_radius,
)


def hohmann_delta_v(alt1_km: float, alt2_km: float, mu: float = 398600.4418) -> dict:
    """
    Compute the two impulsive burns for a standard Hohmann transfer.

    The Hohmann transfer is the minimum delta-v solution for moving between
    two circular, coplanar orbits using exactly two burns:
        Burn 1 (dv1): at r1, accelerate into the transfer ellipse
        Burn 2 (dv2): at r2, circularise into the target orbit

    Args:
        alt1_km : Initial orbit altitude above Earth's surface [km]
        alt2_km : Target  orbit altitude above Earth's surface [km]
        mu      : Gravitational parameter [km3/s2]

    Returns:
        dict with keys:
            dv1_km_s    - first burn  [km/s]
            dv2_km_s    - second burn [km/s]
            total_km_s  - total dv   [km/s]
            r1_km       - inner orbit radius [km]
            r2_km       - outer orbit radius [km]
            a_transfer  - transfer ellipse semi-major axis [km]
            v_circ_1    - circular speed at r1 [km/s]
            v_circ_2    - circular speed at r2 [km/s]
    """
    if math.isclose(alt1_km, alt2_km, rel_tol=1e-6):
        raise ValueError("Initial and target orbits are identical.")

    r1_km = altitude_to_radius(alt1_km)
    r2_km = altitude_to_radius(alt2_km)

    if r1_km > r2_km:
        r1_km, r2_km = r2_km, r1_km

    v_circ_1     = circular_velocity(r1_km, mu)
    v_circ_2     = circular_velocity(r2_km, mu)
    a_transfer   = transfer_semimajor_axis(r1_km, r2_km)
    v_trans_peri = ellipse_velocity(r1_km, a_transfer, mu)
    v_trans_apo  = ellipse_velocity(r2_km, a_transfer, mu)

    dv1 = v_trans_peri - v_circ_1
    dv2 = v_circ_2     - v_trans_apo

    return {
        'dv1_km_s'   : dv1,
        'dv2_km_s'   : dv2,
        'total_km_s' : dv1 + dv2,
        'r1_km'      : r1_km,
        'r2_km'      : r2_km,
        'a_transfer' : a_transfer,
        'v_circ_1'   : v_circ_1,
        'v_circ_2'   : v_circ_2,
    }


def fast_transfer_physics(r_apo_km: float, r1_km: float, r2_km: float,
                           mu: float = 398600.4418) -> tuple:
    """
    Compute delta-v and TOF for a generalised two-burn transfer.

    When r_apo == r2: standard Hohmann — spacecraft coasts to apoapsis exactly
    at r2. Both burns are purely tangential (no radial velocity at apse points).

    When r_apo > r2: fast transfer — spacecraft crosses r2 altitude on the
    ASCENDING leg before reaching apoapsis. At this crossing point the velocity
    vector has BOTH tangential and radial components. The delta-v for Burn 2
    must be computed as a VECTOR difference, not a scalar difference.

    Physics of Burn 2 when r_apo > r2:
    ─────────────────────────────────
    At the crossing point on the ellipse at radius r2:

      Angular momentum:       h = sqrt(mu * a * (1 - e^2))
      Tangential velocity:    v_tan = h / r2
      Total speed (vis-viva): v_total = sqrt(mu * (2/r2 - 1/a))
      Radial velocity:        v_rad = sqrt(v_total^2 - v_tan^2)

    The target circular orbit has:
      v_circ = sqrt(mu / r2)   (purely tangential, zero radial component)

    The vector delta-v required:
      dv2 = sqrt((v_tan - v_circ)^2 + v_rad^2)

    This is always >= the Hohmann dv2, as physically required.

    Args:
        r_apo_km : Apoapsis of the transfer ellipse [km]
        r1_km    : Departure orbit radius [km]
        r2_km    : Target orbit radius [km]
        mu       : Gravitational parameter [km3/s2]

    Returns:
        (delta_v_total [km/s], tof_seconds [s])
        Returns (inf, inf) if geometry is invalid.
    """
    if r_apo_km < r2_km:
        return (np.inf, np.inf)

    # Transfer ellipse parameters
    a   = (r1_km + r_apo_km) / 2.0
    ecc = (r_apo_km - r1_km) / (r_apo_km + r1_km)

    # ── Burn 1: purely tangential at periapsis (r1) ───────────────────────────
    # At periapsis the radial velocity is zero by definition, so this is exact.
    v_circ_1 = circular_velocity(r1_km, mu)
    v_peri   = ellipse_velocity(r1_km, a, mu)
    dv1      = v_peri - v_circ_1          # always positive (prograde)

    # ── Burn 2: vector burn at r2 ─────────────────────────────────────────────
    # Angular momentum magnitude (constant along the orbit)
    h = math.sqrt(mu * a * (1.0 - ecc**2))

    # Tangential velocity component at r2
    v_tan = h / r2_km

    # Total speed at r2 from Vis-Viva
    v_total = ellipse_velocity(r2_km, a, mu)

    # Radial velocity component (Pythagorean decomposition)
    # v_rad = 0 when r_apo == r2 (apoapsis), grows as r_apo increases
    v_rad_sq = max(0.0, v_total**2 - v_tan**2)   # clamp for float safety
    v_rad    = math.sqrt(v_rad_sq)

    # Target: circular orbit at r2 is purely tangential
    v_circ_2 = circular_velocity(r2_km, mu)

    # Vector delta-v: difference between ellipse velocity vector and
    # circular velocity vector (which has zero radial component)
    dv2 = math.sqrt((v_tan - v_circ_2)**2 + v_rad**2)

    delta_v_total = dv1 + dv2

    # ── Time of flight via Kepler's equation ──────────────────────────────────
    def true_anomaly(r: float) -> float:
        """True anomaly at radius r on this ellipse."""
        cos_theta = (a * (1.0 - ecc**2) / r - 1.0) / ecc
        cos_theta = max(-1.0, min(1.0, cos_theta))
        return math.acos(cos_theta)

    def eccentric_anomaly(theta: float) -> float:
        """Eccentric anomaly from true anomaly."""
        return 2.0 * math.atan(
            math.sqrt((1.0 - ecc) / (1.0 + ecc)) * math.tan(theta / 2.0)
        )

    theta1 = true_anomaly(r1_km)    # = 0 at periapsis
    theta2 = true_anomaly(r2_km)    # angle at r2 crossing
    E1     = eccentric_anomaly(theta1)
    E2     = eccentric_anomaly(theta2)
    M1     = E1 - ecc * math.sin(E1)
    M2     = E2 - ecc * math.sin(E2)

    T_period    = 2.0 * math.pi * math.sqrt(a**3 / mu)
    tof_seconds = (M2 - M1) / (2.0 * math.pi) * T_period

    return (delta_v_total, tof_seconds)


def transfer_time(a_transfer_km: float, mu: float = 398600.4418) -> dict:
    """
    Time of flight for a standard Hohmann transfer (half an ellipse).

    Kepler's Third Law:
        T = 2*pi * sqrt(a^3 / mu)
        t_transfer = T / 2

    Args:
        a_transfer_km : Semi-major axis of the transfer ellipse [km]
        mu            : Gravitational parameter [km3/s2]

    Returns:
        dict with tof_seconds, tof_minutes, tof_hours
    """
    period_s = 2.0 * math.pi * math.sqrt(a_transfer_km**3 / mu)
    tof_s    = period_s / 2.0

    return {
        'tof_seconds' : tof_s,
        'tof_minutes' : tof_s / 60.0,
        'tof_hours'   : tof_s / 3600.0,
    }
