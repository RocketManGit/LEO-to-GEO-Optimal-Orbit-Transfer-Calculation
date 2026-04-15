# Hohmann Transfer Optimizer 🛰️
### Finding the Smartest Path from Low Earth Orbit to Geostationary Orbit

---

## The Story Behind This Project

Every satellite launched into space faces the same fundamental challenge: getting from where the rocket drops it off to where it actually needs to be.

Most commercial communications satellites — the ones powering your television, internet, and weather forecasts — need to reach **Geostationary Orbit (GEO)**, a precise ring 35,786 km above the equator where a satellite's orbital speed perfectly matches the Earth's rotation, making it appear to "hover" over the same spot on the ground indefinitely.

But rockets don't deliver satellites there directly. They drop them off in a much lower orbit — **Low Earth Orbit (LEO)**, just a few hundred kilometres up — because reaching GEO in a single burn would require enormous amounts of fuel. Instead, satellites perform a carefully choreographed series of engine burns to climb from LEO to GEO in the most fuel-efficient way possible.

The classical solution to this problem is called a **Hohmann Transfer** — an elegantly simple two-burn manoeuvre first described mathematically by German engineer Walter Hohmann in 1925. It remains, a century later, the backbone of how we move satellites through space.

**This tool computes, optimises, and visualises that manoeuvre** — and goes one step further by exploring the trade-off between fuel efficiency and transfer time.

---

## The Core Problem: Fuel vs. Speed

The Hohmann Transfer is the most *fuel-efficient* way to move between two circular orbits using two engine burns. But fuel efficiency comes at a cost: time.

A standard LEO-to-GEO Hohmann Transfer takes approximately **5.26 hours**. That's acceptable for most missions. But what if a satellite operator needs their asset in position faster? Perhaps due to an emergency coverage gap, a launch window constraint, or a commercial deadline?

By pushing the transfer ellipse higher than the target orbit, the spacecraft can arrive at GEO altitude sooner — shaving hours off the journey. The catch: it burns more fuel to do so, which means either carrying more propellant at launch (expensive) or having less fuel left for the satellite's operational lifetime (costly in a different way).

**This is the trade-off this tool is designed to explore and optimise.**

---

## What This Tool Does

Given a starting orbit and a target orbit, the Hohmann Transfer Optimizer does three things:

**1. Computes the reference Hohmann transfer**
The classical, fuel-optimal solution — including total velocity change required (Δv), flight time, and how much propellant the spacecraft will consume.

**2. Maps the full Pareto frontier**
Using mathematical optimisation, it systematically finds every meaningful combination of fuel cost and transfer time — from the most fuel-efficient option all the way to the fastest possible transfer. This gives mission planners a complete picture of their options, not just a single answer.

**3. Visualises the results**
Produces a clear orbital diagram showing the transfer geometry, alongside a three-panel chart mapping the Pareto frontier so decision-makers can immediately see the consequences of any trade-off.

---

## What the Numbers Look Like

For a standard LEO → GEO transfer (400 km to 35,786 km altitude):

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

And here is what the trade-off looks like across different mission priorities:

| Mission Priority | Fuel Cost (Δv) | Transfer Time | Fuel Penalty | Time Saved |
|---|---|---|---|---|
| Maximum fuel efficiency (classical Hohmann) | 3.90 km/s | 5.26 hrs | Baseline | Baseline |
| Balanced (equal weight on fuel and time) | ~4.10 km/s | ~3.8 hrs | +5.1% | ~28% faster |
| Maximum speed (minimum time) | ~4.60 km/s | ~2.1 hrs | +17.8% | ~60% faster |

A mission planner can read this table and make an informed decision: *"We can get there in half the time if we're willing to burn 18% more fuel. Is that worth it given our remaining propellant budget?"*

---

## Why It Matters

| Without This Tool | With This Tool |
|---|---|
| Engineers compute a single transfer option | A full spectrum of options is generated automatically |
| Trade-offs are estimated manually or ignored | The Pareto frontier is mapped precisely |
| Fuel vs. time decisions are made on intuition | Decisions are backed by quantitative analysis |
| Visualisation requires separate tools | Orbital diagrams and trade-off charts are produced in one run |

---

## The Physics (For the Curious)

A Hohmann Transfer uses exactly two engine burns:

- **Burn 1** fires at the low orbit, accelerating the spacecraft onto an elliptical path that climbs toward the target altitude
- **Burn 2** fires at the high point of that ellipse, circularising the spacecraft into the target orbit

The key insight is that burns performed at high speed (close to Earth) are disproportionately effective — a concept known as the **Oberth Effect**. Hohmann's manoeuvre exploits this naturally, which is why it remains the fuel-optimal two-burn solution after a century.

The optimiser in this tool extends the classical approach by asking: *what if we raised the apex of the ellipse above the target?* The spacecraft would then intersect the target altitude on the way up, arriving sooner. The optimiser finds the mathematically exact cost of every possible version of this idea.

---

## Project Structure (For the Technical Team)

```
hohmann-optimizer/
│
├── main.py                        ← Entry point — run this
├── requirements.txt
│
├── physics/
│   ├── orbital_mechanics.py       ← Core orbital velocity and geometry calculations
│   └── maneuvers.py               ← Hohmann Δv, fast transfer, and flight time equations
│
├── optimizer/
│   ├── objective.py               ← The trade-off cost function
│   └── pareto.py                  ← Pareto frontier solver
│
├── analysis/
│   └── propellant.py              ← Propellant mass calculations (Tsiolkovsky rocket equation)
│
└── visualization/
    ├── orbit_plot.py              ← 2D orbital diagram
    └── pareto_plot.py             ← Three-panel Pareto trade-off chart
```

---

## Getting Started (For Developers)

```bash
# 1. Clone the repository
git clone https://github.com/RocketManGit/LEO-to-GEO-Optimal-Orbit-Transfer-Calculation.git
cd LEO-to-GEO-Optimal-Orbit-Transfer-Calculation

# 2. Install dependencies
pip install -r requirements.txt

# 3. Run
python main.py
```

To change the mission parameters (orbit altitudes, spacecraft mass, engine efficiency), edit the **MISSION PARAMETERS** block at the top of `main.py`.

---

## Dependencies

```
numpy       — numerical calculations
scipy       — optimisation engine (L-BFGS-B algorithm)
matplotlib  — plotting and visualisation
```

---

## References & Further Reading

- Bate, Mueller & White — *Fundamentals of Astrodynamics* (Dover, 1971)
- Curtis, H. — *Orbital Mechanics for Engineering Students* (Elsevier, 2013)
- Vallado, D. — *Fundamentals of Astrodynamics and Applications* (Microcosm, 2013)

---

*Built as a learning project for aerospace engineering students — bridging orbital mechanics theory and practical mission design.*
