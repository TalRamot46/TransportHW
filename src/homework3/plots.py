"""The order scan shared by Questions 3 and 4: the same three panels in either geometry."""

import numpy as np
from typing import NamedTuple
from homework1.criticality import critical_dimensions
from homework3.figures import subplots, panel, finish, savefig
from homework3.sn import multiplying_medium

C_VALUES = (1.2, 1.5, 1.8)
ORDERS = (2, 4, 6, 10)
REFERENCE_C = 1.5              # the c whose k curve and flux shapes the figure shows
CURVE_POINTS = 21

COLORS = {2: '#2c3e50', 4: '#e67e22', 6: '#27ae60', 10: '#c0392b'}

class OrderScan(NamedTuple):
    """Critical sizes for every (c, N), and the k curve and flux shapes at REFERENCE_C."""
    sizes: dict          # {c: {N: critical size}}
    grid: np.ndarray     # sizes at which k was evaluated, at REFERENCE_C
    curves: dict         # {N: k on that grid}
    fluxes: dict         # {N: (mesh, scalar flux)} of the critical system
    index: int           # which entry of critical_dimensions is this geometry's size

def reference_size(c, index):
    """Assignment 1's exact-transport critical half-thickness (index 2) or radius (3)."""
    return critical_dimensions(c, 'transport-ref')[index]

def scan_orders(critical_size, k_result, index):
    """Runs the whole (c, N) scan for one geometry; `index` picks its reference size."""
    sizes = {c: {N: critical_size(multiplying_medium(c), N, reference_size(c, index))
                 for N in ORDERS}
             for c in C_VALUES}

    medium = multiplying_medium(REFERENCE_C)
    critical = sizes[REFERENCE_C]
    grid = np.linspace(0.5 * min(critical.values()), 1.6 * max(critical.values()),
                       CURVE_POINTS)

    curves, fluxes = {}, {}
    for N in ORDERS:
        curves[N] = np.array([k_result(size, medium, N).k for size in grid])
        result = k_result(critical[N], medium, N)
        fluxes[N] = (result.x / critical[N], result.phi / result.phi[0])

    return OrderScan(sizes, grid, curves, fluxes, index)

def _panel_orders(ax, scan, symbol):
    """Critical size against N, with Assignment 1's exact-transport value as a dashed line."""
    for i, c in enumerate(C_VALUES):
        color = f'C{i}'
        ax.plot(ORDERS, [scan.sizes[c][N] for N in ORDERS], 'o-', color=color,
                linewidth=2.0, label=f'$c = {c}$')
        ax.axhline(reference_size(c, scan.index), color=color, linestyle='--',
                   linewidth=1.2, alpha=0.7)
    panel(ax, f'Convergence of {symbol} with $N$', 'Order $N$',
          f'{symbol} [mean free paths]', legend='upper right', fontsize=9)

def _panel_curves(ax, scan, symbol):
    """k against system size at REFERENCE_C, one curve per order, with the k = 1 crossings."""
    for N in ORDERS:
        ax.plot(scan.grid, scan.curves[N], color=COLORS[N], linewidth=2.0, label=f'$S_{{{N}}}$')
        ax.plot([scan.sizes[REFERENCE_C][N]], [1.0], color=COLORS[N], marker='o',
                markersize=7, markerfacecolor='none', markeredgewidth=1.8)
    ax.axhline(1.0, color='grey', linestyle=':', linewidth=1.4)
    panel(ax, f'$k$ against size, $c = {REFERENCE_C}$', f'{symbol} [mean free paths]',
          'Multiplication factor $k$', legend='lower right', fontsize=9)

def _panel_fluxes(ax, scan):
    """Scalar flux of the critical system, normalised, one curve per order."""
    for N in ORDERS:
        x, phi = scan.fluxes[N]
        ax.plot(x, phi, color=COLORS[N], linewidth=2.0, label=f'$S_{{{N}}}$')
    panel(ax, f'Critical flux shape, $c = {REFERENCE_C}$', 'Position / critical size',
          r'$\phi / \phi(0)$', legend='lower left', fontsize=9)

def plot_orders(scan, symbol, suptitle, save_path):
    """Writes the three-panel order-scan figure of one geometry."""
    fig, axes = subplots(1, 3, width=5.2, height=4.2)
    _panel_orders(axes[0][0], scan, symbol)
    _panel_curves(axes[0][1], scan, symbol)
    _panel_fluxes(axes[0][2], scan)
    finish(fig, suptitle)
    savefig(fig, save_path)
