# visualization/orbit_plot.py
# ─────────────────────────────────────────────────────────────────────────────
# 2D orbital diagrams for the Hohmann transfer project.
#
# Two plots:
#   plot_transfer()      — standard Hohmann transfer diagram
#   plot_pareto_orbits() — overlaid transfer ellipses across the Pareto front
# ─────────────────────────────────────────────────────────────────────────────

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from physics.orbital_mechanics import EARTH_RADIUS_KM


# ── Shared helpers ────────────────────────────────────────────────────────────

def _ellipse_upper_half(r1_km: float, r_apo_km: float) -> tuple:
    """
    Return (x, y) arrays for the upper half of a transfer ellipse.

    Earth sits at the left focus. Periapsis is at (+r1, 0) on the right,
    apoapsis at (-r_apo, 0) on the left. Sweep theta 0 → pi.
    """
    a  = (r1_km + r_apo_km) / 2.0
    c  = a - r1_km
    b  = math.sqrt(a**2 - c**2)
    cx = -c
    theta = np.linspace(0, np.pi, 600)
    return cx + a * np.cos(theta), b * np.sin(theta)


def _draw_earth_and_orbits(ax, r1_km: float, r2_km: float) -> None:
    """Draw Earth circle, LEO ring, and GEO ring."""
    theta = np.linspace(0, 2 * np.pi, 600)
    ax.add_patch(plt.Circle((0, 0), EARTH_RADIUS_KM,
                             color='#1a6fa8', zorder=3))
    ax.text(0, 0, 'Earth', fontsize=9, ha='center', va='center',
            color='white', fontweight='bold', zorder=4)
    ax.plot(0, 0, '+', color='white', markersize=7, zorder=5, linewidth=1.2)
    ax.plot(r1_km * np.cos(theta), r1_km * np.sin(theta),
            color='#4ade80', linewidth=1.4, label=f'LEO  ({r1_km - EARTH_RADIUS_KM:.0f} km alt)')
    ax.plot(r2_km * np.cos(theta), r2_km * np.sin(theta),
            color='#60a5fa', linewidth=1.4, label=f'GEO  ({r2_km - EARTH_RADIUS_KM:.0f} km alt)')


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 1 — Standard Hohmann Transfer
# ═══════════════════════════════════════════════════════════════════════════════

def plot_transfer(alt1_km: float, alt2_km: float, dv_result: dict,
                  save_path: str = "outputs/orbit.png") -> None:
    """
    Draw the standard Hohmann transfer: Earth, two circular orbits,
    the transfer ellipse (upper half), and the two burn markers.
    """
    r1 = dv_result['r1_km']
    r2 = dv_result['r2_km']

    x_ellipse, y_ellipse = _ellipse_upper_half(r1, r2)

    fig, ax = plt.subplots(figsize=(10, 9), facecolor='#0d0d1a')
    ax.set_facecolor('#0d0d1a')
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(
        f"Hohmann Transfer  |  {alt1_km:.0f} km \u2192 {alt2_km:.0f} km altitude",
        color='white', fontsize=13, pad=14, fontweight='bold'
    )

    _draw_earth_and_orbits(ax, r1, r2)

    # Transfer ellipse
    ax.plot(x_ellipse, y_ellipse,
            color='#fb923c', linewidth=2.2, linestyle='--',
            label='Transfer ellipse')

    # Burn 1 — periapsis, right side
    ax.plot(r1, 0, 'o', color='#facc15', markersize=11, zorder=6)
    ax.annotate(
        f"Burn 1  +{dv_result['dv1_km_s']:.3f} km/s",
        xy=(r1, 0), xytext=(r1 + r2 * 0.18, r2 * 0.18),
        color='#facc15', fontsize=9, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color='#facc15', lw=1.3)
    )

    # Burn 2 — apoapsis, left side
    ax.plot(-r2, 0, 'o', color='#f472b6', markersize=11, zorder=6)
    ax.annotate(
        f"Burn 2  +{dv_result['dv2_km_s']:.3f} km/s",
        xy=(-r2, 0), xytext=(-r2 + r2 * 0.18, r2 * 0.18),
        color='#f472b6', fontsize=9, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color='#f472b6', lw=1.3)
    )

    margin = r2 * 1.22
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin * 0.6, margin)
    ax.legend(loc='lower right', facecolor='#1a1a2e',
              edgecolor='#444', labelcolor='white', fontsize=9)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.show()
    print(f"  [Orbit plot saved → {save_path}]")


# ═══════════════════════════════════════════════════════════════════════════════
# PLOT 2 — Pareto Front Orbital Evolution
# ═══════════════════════════════════════════════════════════════════════════════

def plot_pareto_orbits(pareto_data: dict, alt1_km: float, alt2_km: float,
                       n_samples: int = 6,
                       save_path: str = "outputs/pareto_orbits.png") -> None:
    """
    Overlay transfer ellipses from across the Pareto front on one diagram.

    Samples n_samples evenly spaced points from the Pareto front (always
    including the Hohmann solution). Each ellipse is colored by its Δv
    penalty: dark purple = Hohmann (0%) → bright yellow = most expensive.

    Labels are placed at the top of each ellipse arc (the midpoint, where
    y is maximum) to avoid overlap with each other and with the apoapsis
    annotation cluster on the left.
    """
    r1     = pareto_data['r1_km']
    r2     = pareto_data['r2_km']
    dv_arr = pareto_data['dv_km_s']
    tof_arr= pareto_data['tof_hours']
    r_apo  = pareto_data['r_apo_km']
    h_dv   = pareto_data['hohmann_dv']

    # ── Sample Pareto front evenly, always include Hohmann (last index) ───────
    indices     = np.linspace(0, len(dv_arr) - 1, n_samples, dtype=int)
    # Deduplicate while preserving order (can happen with small n_points)
    seen = set()
    indices = [i for i in indices if not (i in seen or seen.add(i))]

    penalties   = [(dv_arr[i] - h_dv) / h_dv * 100 for i in indices]
    max_penalty = max(p for p in penalties if p > 0) if any(p > 0 for p in penalties) else 1.0

    cmap      = cm.plasma
    norm      = mcolors.Normalize(vmin=0, vmax=max_penalty)
    color_map = cm.ScalarMappable(norm=norm, cmap=cmap)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 9), facecolor='#0d0d1a')
    ax.set_facecolor('#0d0d1a')
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(
        f"Transfer Ellipse Evolution Across Pareto Front\n"
        f"{alt1_km:.0f} km \u2192 {alt2_km:.0f} km  |  "
        f"Colored by \u0394v penalty vs Hohmann",
        color='white', fontsize=12, pad=14, fontweight='bold'
    )

    _draw_earth_and_orbits(ax, r1, r2)

    # ── Draw ellipses ─────────────────────────────────────────────────────────
    for idx, i in enumerate(indices):
        r_apo_i   = r_apo[i]
        penalty_i = penalties[idx]
        tof_i     = tof_arr[i]
        color_i   = color_map.to_rgba(penalty_i)
        is_hohmann = (i == indices[-1])

        lw    = 2.4 if is_hohmann else 1.6
        ls    = '-'  if is_hohmann else '--'
        alpha = 1.0  if is_hohmann else 0.9

        x_ell, y_ell = _ellipse_upper_half(r1, r_apo_i)
        ax.plot(x_ell, y_ell, color=color_i, linewidth=lw,
                linestyle=ls, alpha=alpha)

        # ── Label at the arc top (maximum y point) ────────────────────────────
        # This keeps all labels in the upper half of the figure,
        # spread naturally by their differing arc heights.
        top_idx   = np.argmax(y_ell)
        x_top     = x_ell[top_idx]
        y_top     = y_ell[top_idx]

        if is_hohmann:
            label = f"Hohmann\n\u0394v=+0.0%  {tof_i:.2f} hrs"
        else:
            label = f"\u0394v=+{penalty_i:.1f}%\n{tof_i:.2f} hrs"

        # Offset label slightly above the arc top
        ax.text(x_top, y_top + r2 * 0.04, label,
                color=color_i, fontsize=7.5, ha='center', va='bottom',
                fontweight='bold' if is_hohmann else 'normal')

    # ── Burn 1 marker ─────────────────────────────────────────────────────────
    ax.plot(r1, 0, 'o', color='#facc15', markersize=10, zorder=6)
    ax.annotate(
        'Burn 1\n(all transfers)',
        xy=(r1, 0), xytext=(r1 + r2 * 0.25, r2 * 0.14),
        color='#facc15', fontsize=8, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color='#facc15', lw=1.2)
    )

    # ── Orbit labels ──────────────────────────────────────────────────────────
    ax.text(r2 * 0.06, r1 * 1.18, 'LEO', color='#4ade80',
            fontsize=8, ha='center', va='bottom')
    ax.text(r2 * 0.06, r2 * 1.03, 'GEO', color='#60a5fa',
            fontsize=8, ha='center', va='bottom')

    # ── Colorbar ──────────────────────────────────────────────────────────────
    cbar = fig.colorbar(color_map, ax=ax, fraction=0.022, pad=0.01)
    cbar.set_label('\u0394v penalty vs Hohmann  [%]', color='#ccc', fontsize=9)
    cbar.ax.yaxis.set_tick_params(color='#888')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='#aaa', fontsize=8)

    # ── Axis limits ───────────────────────────────────────────────────────────
    r_apo_max = r_apo[indices[0]]          # widest ellipse = time-optimal
    x_margin  = max(r_apo_max, r2) * 1.15
    y_top_max = max(                        # tallest arc height across samples
        math.sqrt(((r1 + r_apo[i]) / 2)**2 -
                  (((r1 + r_apo[i]) / 2) - r1)**2)
        for i in indices
    )
    ax.set_xlim(-x_margin, r2 * 1.15)
    ax.set_ylim(-r2 * 0.15, y_top_max * 1.22)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.show()
    print(f"  [Pareto orbit plot saved → {save_path}]")
