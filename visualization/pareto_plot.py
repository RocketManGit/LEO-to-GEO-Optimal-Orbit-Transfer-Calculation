# visualization/pareto_plot.py
# ─────────────────────────────────────────────────────────────────────────────
# Three-panel Pareto front visualization.
#   Panel 1: Δv vs TOF — the Pareto front, colored by weight w
#   Panel 2: Optimal r_apo vs weight w
#   Panel 3: % Δv penalty vs % TOF saving relative to Hohmann
# ─────────────────────────────────────────────────────────────────────────────

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def plot_pareto_front(data: dict, alt1_km: float, alt2_km: float,
                      save_path: str = "outputs/pareto_front.png") -> None:
    """
    Render and save the three-panel Pareto front plot.

    Args:
        data      : Output dict from solve_pareto_front()
        alt1_km   : Departure orbit altitude [km]  — for the title
        alt2_km   : Target   orbit altitude [km]   — for the title
        save_path : File path to save the figure
    """
    dv      = data['dv_km_s']
    tof     = data['tof_hours']
    weights = data['weights']
    r_apo   = data['r_apo_km']
    h_dv    = data['hohmann_dv']
    h_tof   = data['hohmann_tof']
    r2      = data['r2_km']

    fig = plt.figure(figsize=(15, 6), facecolor='#0d0d1a')
    gs  = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38)

    # ── Panel 1: Pareto Front ─────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0])
    ax1.set_facecolor('#0d0d1a')

    sc = ax1.scatter(tof, dv, c=weights, cmap='plasma',
                     s=28, zorder=3, edgecolors='none')
    ax1.scatter(h_tof, h_dv, color='#4ade80', s=120, zorder=5,
                marker='*', label='Hohmann (Δv-optimal)')
    ax1.scatter(tof[0], dv[0], color='#f472b6', s=120, zorder=5,
                marker='D', label='Time-optimal (w=0)')

    cbar = fig.colorbar(sc, ax=ax1, pad=0.02)
    cbar.set_label('Weight w  (1 = minimize Δv)', color='#aaa', fontsize=8)
    cbar.ax.yaxis.set_tick_params(color='#aaa')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='#aaa')

    _style_axes(ax1,
                xlabel='Time of Flight  [hours]',
                ylabel='Total Δv  [km/s]',
                title='Pareto Front\n(each point = one scipy.optimize solve)')
    ax1.legend(fontsize=7, facecolor='#1a1a2e', edgecolor='#444',
               labelcolor='white', loc='upper right')

    # ── Panel 2: Optimal r_apo vs weight ─────────────────────────────────────
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor('#0d0d1a')

    ax2.plot(weights, r_apo / 1000, color='#fb923c', linewidth=2)
    ax2.axhline(r2 / 1000, color='#60a5fa', linewidth=1.2, linestyle='--',
                label=f'Target orbit r₂ = {r2/1000:.1f} Mm')

    _style_axes(ax2,
                xlabel='Weight w',
                ylabel='Optimal apoapsis  [×10³ km]',
                title='Optimal Apoapsis Radius\nvs. Scalarization Weight')
    ax2.legend(fontsize=7, facecolor='#1a1a2e', edgecolor='#444',
               labelcolor='white')

    # ── Panel 3: % trade-off ──────────────────────────────────────────────────
    ax3 = fig.add_subplot(gs[2])
    ax3.set_facecolor('#0d0d1a')

    dv_penalty_pct = (dv   - h_dv)  / h_dv  * 100
    tof_saving_pct = (h_tof - tof)  / h_tof * 100

    ax3.scatter(tof_saving_pct, dv_penalty_pct,
                c=weights, cmap='plasma', s=28, zorder=3, edgecolors='none')
    ax3.axhline(0, color='#4ade80', linewidth=0.8, linestyle='--')
    ax3.axvline(0, color='#4ade80', linewidth=0.8, linestyle='--')

    _style_axes(ax3,
                xlabel='TOF saving vs. Hohmann  [%]',
                ylabel='Δv penalty vs. Hohmann  [%]',
                title='Trade-off: Time Saved\nvs. Fuel Cost Penalty')

    fig.suptitle(
        f"Hohmann Transfer Optimizer  |  {alt1_km:.0f} km → {alt2_km:.0f} km altitude",
        color='white', fontsize=12, fontweight='bold', y=1.02
    )

    plt.savefig(save_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.show()
    print(f"  [Pareto plot saved → {save_path}]")


def _style_axes(ax, xlabel: str, ylabel: str, title: str) -> None:
    """Apply consistent dark-theme styling to an axes object."""
    ax.set_xlabel(xlabel, color='#ccc', fontsize=9)
    ax.set_ylabel(ylabel, color='#ccc', fontsize=9)
    ax.set_title(title,   color='white', fontsize=9, fontweight='bold')
    ax.tick_params(colors='#888')
    ax.spines[:].set_color('#333')
    ax.grid(True, color='#222', linewidth=0.6)
