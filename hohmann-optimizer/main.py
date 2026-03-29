# main.py
# ─────────────────────────────────────────────────────────────────────────────
# Hohmann Transfer Optimizer — Entry Point
#
# This is the only file you need to run. Edit the MISSION PARAMETERS block
# below, then execute:
#
#     python main.py
#
# What it does:
#   1. Computes the reference Hohmann transfer (Δv, TOF, propellant)
#   2. Runs the scalarized Pareto front optimizer (scipy.optimize)
#   3. Prints a mission report
#   4. Generates three plots: orbital diagram, Pareto front, Pareto orbit evolution
# ─────────────────────────────────────────────────────────────────────────────

import os

# Ensure outputs directory exists
os.makedirs("outputs", exist_ok=True)

from physics    import hohmann_delta_v, transfer_time
from analysis   import propellant_mass
from optimizer  import solve_pareto_front
from visualization import plot_transfer, plot_pareto_front, plot_pareto_orbits


# ══════════════════════════════════════════════════════════════════════════════
# MISSION PARAMETERS — edit these
# ══════════════════════════════════════════════════════════════════════════════

ALT1_KM            = 400       # Departure orbit altitude [km]  (e.g. ISS-like LEO)
ALT2_KM            = 35786     # Target   orbit altitude [km]  (GEO)

SPACECRAFT_MASS_KG = 5000      # Wet mass — spacecraft + full propellant load [kg]
ENGINE_ISP_S       = 311       # Specific impulse [s]  (hypergolic, e.g. R-4D)

N_PARETO_POINTS    = 60        # Pareto front resolution (more = smoother curve)
R_APO_MAX_FACTOR   = 8.0       # Optimizer search ceiling: r_apo ≤ r2 * this


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Hohmann reference transfer
# ══════════════════════════════════════════════════════════════════════════════

dv   = hohmann_delta_v(ALT1_KM, ALT2_KM)
tof  = transfer_time(dv['a_transfer'])
prop = propellant_mass(dv['total_km_s'], SPACECRAFT_MASS_KG, ENGINE_ISP_S)

_W = 56
print("\n" + "="*_W)
print("  HOHMANN REFERENCE — MISSION REPORT")
print("="*_W)
print(f"  Initial altitude        : {ALT1_KM:>10.1f} km")
print(f"  Target  altitude        : {ALT2_KM:>10.1f} km")
print(f"  Transfer ellipse SMA    : {dv['a_transfer']:>10.1f} km")
print("-"*_W)
print(f"  Δv₁  (departure burn)   : {dv['dv1_km_s']:>10.4f} km/s")
print(f"  Δv₂  (arrival burn)     : {dv['dv2_km_s']:>10.4f} km/s")
print(f"  Total Δv                : {dv['total_km_s']:>10.4f} km/s")
print("-"*_W)
print(f"  Time of flight          : {tof['tof_hours']:>10.2f} hrs")
print("-"*_W)
print(f"  Initial (wet) mass      : {prop['m_initial_kg']:>10.1f} kg")
print(f"  Propellant consumed     : {prop['m_prop_kg']:>10.1f} kg")
print(f"  Final   (dry) mass      : {prop['m_final_kg']:>10.1f} kg")
print(f"  Propellant mass frac.   : {prop['mass_fraction']*100:>10.1f} %")
print("="*_W)


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Pareto front optimizer
# ══════════════════════════════════════════════════════════════════════════════

pareto_data = solve_pareto_front(
    ALT1_KM, ALT2_KM,
    n_points         = N_PARETO_POINTS,
    r_apo_max_factor = R_APO_MAX_FACTOR,
)

# Print a summary table of selected Pareto-optimal solutions
h_dv  = pareto_data['hohmann_dv']
h_tof = pareto_data['hohmann_tof']

print("\n" + "="*72)
print("  PARETO-OPTIMAL SOLUTIONS  (selected weights)")
print(f"  Hohmann reference → Δv = {h_dv:.4f} km/s  |  TOF = {h_tof:.3f} hrs")
print("="*72)
print(f"  {'w':>5}  {'Δv [km/s]':>12}  {'TOF [hrs]':>10}"
      f"  {'Δv penalty':>12}  {'TOF saving':>12}  {'r_apo [km]':>12}")
print("-"*72)

import numpy as np
indices = np.linspace(0, N_PARETO_POINTS - 1, 10, dtype=int)
for i in indices:
    w       = pareto_data['weights'][i]
    dv_i    = pareto_data['dv_km_s'][i]
    tof_i   = pareto_data['tof_hours'][i]
    rapo_i  = pareto_data['r_apo_km'][i]
    dv_pen  = (dv_i  - h_dv)  / h_dv  * 100
    tof_sav = (h_tof - tof_i) / h_tof * 100
    print(f"  {w:>5.2f}  {dv_i:>12.4f}  {tof_i:>10.3f}"
          f"  {dv_pen:>+11.1f}%  {tof_sav:>+11.1f}%  {rapo_i:>12.1f}")

print("="*72)
print(f"\n  → w=1.0 recovers Hohmann exactly")
print(f"  → w=0.0 minimizes time at a Δv cost of "
      f"{(pareto_data['dv_km_s'][0]-h_dv)/h_dv*100:.1f}%")
print(f"  → the 'knee' of the curve is the natural engineering sweet spot\n")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Plots
# ══════════════════════════════════════════════════════════════════════════════

plot_transfer(ALT1_KM, ALT2_KM, dv,
              save_path="outputs/orbit.png")

plot_pareto_front(pareto_data, ALT1_KM, ALT2_KM,
                  save_path="outputs/pareto_front.png")

plot_pareto_orbits(pareto_data, ALT1_KM, ALT2_KM,
                   n_samples=7,
                   save_path="outputs/pareto_orbits.png")
