# physics/orbital_mechanics.py
# ─────────────────────────────────────────────────────────────────────────────
# Fundamental orbital mechanics equations.
# These are stateless, unit-pure functions — no mission logic here.
#
# Unit convention throughout: distances in km, time in seconds, speed in km/s.
# ─────────────────────────────────────────────────────────────────────────────

import math

# ── Constants ─────────────────────────────────────────────────────────────────

MU_EARTH        = 398600.4418   # Standard gravitational parameter  [km³/s²]
EARTH_RADIUS_KM = 6378.137      # Equatorial Earth radius [km]
                                 # GEO altitude (35,786 km) is defined relative
                                 # to this value, so it must match here.


def circular_velocity(radius_km: float, mu: float = MU_EARTH) -> float:
    """
    Speed of a spacecraft in a circular orbit at a given radius.

    Derived by equating gravitational and centripetal acceleration:
        v = sqrt(μ / r)

    Args:
        radius_km : Distance from Earth's centre [km]
        mu        : Gravitational parameter [km³/s²]

    Returns:
        Circular orbital speed [km/s]
    """
    if radius_km <= 0:
        raise ValueError(f"Orbital radius must be positive, got {radius_km} km.")
    return math.sqrt(mu / radius_km)


def ellipse_velocity(radius_km: float, semimajor_km: float, mu: float = MU_EARTH) -> float:
    """
    Speed at any point on an elliptical orbit via the Vis-Viva equation:
        v² = μ * (2/r - 1/a)

    This is the central equation of orbital mechanics — it relates speed,
    current position, and orbit shape without needing to know the time.

    Args:
        radius_km    : Current distance from Earth's centre [km]
        semimajor_km : Semi-major axis of the ellipse [km]
        mu           : Gravitational parameter [km³/s²]

    Returns:
        Speed at that point [km/s]
    """
    if radius_km <= 0 or semimajor_km <= 0:
        raise ValueError("Radius and semi-major axis must be positive.")
    return math.sqrt(mu * (2.0 / radius_km - 1.0 / semimajor_km))


def transfer_semimajor_axis(r1_km: float, r2_km: float) -> float:
    """
    Semi-major axis of the Hohmann transfer ellipse.

    The Hohmann ellipse has periapsis at r1 and apoapsis at r2, so its
    semi-major axis is simply their average:
        a = (r1 + r2) / 2

    Args:
        r1_km : Inner orbit radius [km]
        r2_km : Outer orbit radius [km]

    Returns:
        Semi-major axis [km]
    """
    return (r1_km + r2_km) / 2.0


def altitude_to_radius(altitude_km: float) -> float:
    """
    Convert an altitude above Earth's surface to a radius from Earth's centre.

    Args:
        altitude_km : Altitude above mean surface [km]

    Returns:
        Radius from Earth's centre [km]
    """
    if altitude_km < 0:
        raise ValueError(f"Altitude must be non-negative, got {altitude_km} km.")
    return EARTH_RADIUS_KM + altitude_km
