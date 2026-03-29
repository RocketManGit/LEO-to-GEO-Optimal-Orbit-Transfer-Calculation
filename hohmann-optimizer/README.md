# Hohmann Transfer Optimizer

A Python tool for computing, optimizing, and visualizing Hohmann orbital transfers — built step-by-step as a learning project for aerospace engineering students.

---

## What it does

Given a departure orbit and a target orbit, this tool:

1. **Computes the reference Hohmann transfer** — minimum delta-v solution, time of flight, and propellant mass via the Tsiolkovsky rocket equation
2. **Solves a scalarized Pareto front** — uses `scipy.optimize` to find the optimal transfer ellipse apoapsis for every trade-off between fuel and time
3. **Plots the results** — orbital diagram and a three-panel Pareto front visualization

---

## The physics

A Hohmann transfer moves a spacecraft between two circular, coplanar orbits using exactly two impulsive burns:

- **Burn 1** at the inner orbit — accelerates into a transfer ellipse
- **Burn 2** at the outer orbit — circularizes into the target orbit

This is the fuel-optimal solution for the two-burn case. But it is not time-optimal. By raising the transfer ellipse apoapsis above the target orbit, the spacecraft crosses the target altitude sooner — at the cost of additional delta-v.

The optimizer finds the Pareto front of this trade-off.

---

## Project structure

```
hohmann-optimizer/
│
├── main.py                        ← entry point — run this
├── requirements.txt
│
├── physics/
│   ├── orbital_mechanics.py       ← vis-viva, circular velocity, geometry
│   └── maneuvers.py               ← hohmann Δv, fast transfer, Kepler TOF
│
├── optimizer/
│   ├── objective.py               ← scalarized cost function
│   └── pareto.py                  ← Pareto front solver (scipy L-BFGS-B)
│
├── analysis/
│   └── propellant.py              ← Tsiolkovsky rocket equation
│
└── visualization/
    ├── orbit_plot.py              ← 2D orbital diagram
    └── pareto_plot.py             ← 3-panel Pareto front plot
```

---

## Quickstart

```bash
# 1. Clone the repository
git clone https://github.com/your-username/hohmann-optimizer.git
cd hohmann-optimizer

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run
python main.py
```

Edit the **MISSION PARAMETERS** block at the top of `main.py` to change orbits, spacecraft mass, or engine Isp.

---

## Example output (LEO → GEO)

```
══════════════════════════════════════════════════════════
  HOHMANN REFERENCE — MISSION REPORT
══════════════════════════════════════════════════════════
  Initial altitude        :      400.0 km
  Target  altitude        :    35786.0 km
  Total Δv                :     3.9042 km/s
  Time of flight          :       5.26 hrs
  Propellant consumed     :    3963.3 kg  (79.3% of wet mass)
══════════════════════════════════════════════════════════
```

| Weight w | Δv [km/s] | TOF [hrs] | Δv penalty | TOF saving |
|---|---|---|---|---|
| 1.00 | 3.9042 | 5.256 | +0.0% | +0.0% |
| 0.50 | ~4.10  | ~3.8  | +5.1% | +27.7% |
| 0.00 | ~4.60  | ~2.1  | +17.8% | +60.2% |

---

## Key equations

| Quantity | Equation |
|---|---|
| Circular velocity | `v = sqrt(μ / r)` |
| Vis-Viva | `v² = μ * (2/r - 1/a)` |
| Kepler's 3rd Law | `T = 2π * sqrt(a³ / μ)` |
| Rocket equation | `Δv = Isp * g₀ * ln(m₀ / mf)` |
| Scalarized cost | `J = w * (Δv/Δv_ref) + (1-w) * (TOF/TOF_ref)` |

---

## Dependencies

- `numpy` — array operations
- `scipy` — L-BFGS-B optimizer
- `matplotlib` — plotting

---

## References

- Bate, Mueller & White — *Fundamentals of Astrodynamics* (Dover, 1971)
- Curtis, H. — *Orbital Mechanics for Engineering Students* (Elsevier, 2013)
- Vallado, D. — *Fundamentals of Astrodynamics and Applications* (Microcosm, 2013)
