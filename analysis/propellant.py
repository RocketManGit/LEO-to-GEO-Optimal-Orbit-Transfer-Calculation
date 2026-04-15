# analysis/propellant.py
# ─────────────────────────────────────────────────────────────────────────────
# Propellant mass estimation via the Tsiolkovsky rocket equation.
# Sits in analysis/ because it consumes optimizer outputs, not raw physics.
# ─────────────────────────────────────────────────────────────────────────────

import math


def propellant_mass(total_dv_km_s: float, m_initial_kg: float,
                    isp_s: float, g0: float = 9.80665e-3) -> dict:
    """
    Estimate propellant consumed for a given delta-v budget.

    The Tsiolkovsky rocket equation:
        Δv = Isp * g₀ * ln(m_initial / m_final)

    Rearranged:
        m_final = m_initial * exp(−Δv / (Isp * g₀))
        m_prop  = m_initial − m_final

    Typical Isp values for reference:
        ~220 s  — cold gas thrusters
        ~311 s  — hypergolic engines (e.g. R-4D, used on many GEO sats)
        ~450 s  — cryogenic LH2/LOX (e.g. RL-10)

    Args:
        total_dv_km_s : Total delta-v required [km/s]
        m_initial_kg  : Spacecraft wet mass (including propellant) [kg]
        isp_s         : Engine specific impulse [seconds]
        g0            : Standard gravity in km/s² (default: 9.80665e-3 km/s²)

    Returns:
        dict with keys:
            m_initial_kg  — wet mass [kg]
            m_final_kg    — dry mass after burns [kg]
            m_prop_kg     — propellant consumed [kg]
            mass_fraction — propellant / wet mass  [dimensionless]

    Raises:
        ValueError: if Isp or initial mass are non-positive.
    """
    if isp_s <= 0:
        raise ValueError(f"Isp must be positive, got {isp_s} s.")
    if m_initial_kg <= 0:
        raise ValueError(f"Initial mass must be positive, got {m_initial_kg} kg.")

    v_exhaust  = isp_s * g0                                      # effective exhaust velocity [km/s]
    m_final_kg = m_initial_kg * math.exp(-total_dv_km_s / v_exhaust)
    m_prop_kg  = m_initial_kg - m_final_kg

    return {
        'm_initial_kg'  : m_initial_kg,
        'm_final_kg'    : m_final_kg,
        'm_prop_kg'     : m_prop_kg,
        'mass_fraction' : m_prop_kg / m_initial_kg,
    }
