"""Question 4: the critical radius of a bare sphere from the k = 1 search."""

import os
import logging
import numpy as np
from homework1.spherical import (build_medium, buckling, k_eigenvalue, critical_radius,
                                 analytic_critical_radius, mesh_convergence)
from homework1.figures import subplots, panel, finish, savefig
from homework1.tables import log_section, log_table

logger = logging.getLogger(__name__)

C_VALUES = (1.02, 1.05, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0)
PROFILE_C = (1.05, 1.5, 2.0)
APPROXIMATIONS = ('classical', 'asymptotic')
LABELS = {'classical': 'Classical', 'asymptotic': 'Asymptotic'}
COLORS = ('#2c3e50', '#e74c3c', '#3498db')
N_CELLS = 400

def _radii(c_values, approximation, n_cells=N_CELLS):
    """(numerical, analytic) critical radii over a list of c, in mean free paths."""
    media = [build_medium(1.0, 1.0, c, approximation) for c in c_values]
    return (np.array([critical_radius(m, n_cells=n_cells) for m in media]),
            np.array([analytic_critical_radius(m) for m in media]))

def plot_critical_radius(save_path):
    """Numerical radius against the analytic pi/B - z0, with the difference below."""
    c = np.asarray(C_VALUES, dtype=float)
    fig, axes = subplots(2, 2)

    for j, approximation in enumerate(APPROXIMATIONS):
        numerical, analytic = _radii(c, approximation)

        axes[0][j].plot(c, analytic, label=r'analytic, $\pi / B - z_0$',
                        color='#2c3e50', linewidth=2.6)
        axes[0][j].plot(c, numerical, label=f'numerical, $k = 1$ ({N_CELLS} cells)',
                        color='#e74c3c', linestyle='none', marker='o', markersize=5,
                        markerfacecolor='none', markeredgewidth=1.4)
        panel(axes[0][j], f'{LABELS[approximation]} diffusion', '$c$',
              r'$\Sigma_t R_c$ [mean free paths]', log=True, legend='upper right')

        axes[1][j].plot(c, np.abs(numerical - analytic) / analytic * 100.0,
                        color='#8e44ad', linewidth=2.0, marker='o', markersize=4)
        panel(axes[1][j], None, '$c$', 'Relative difference (%)', log=True)

    finish(fig, 'Question 4: critical radius from the $k = 1$ search vs. analytic')
    savefig(fig, save_path)

def plot_mesh_convergence(save_path):
    """Relative error of the critical radius against mesh refinement, expected 1/N^2."""
    fig, axes = subplots(1, 2, height=5.0)

    for j, approximation in enumerate(APPROXIMATIONS):
        for c, color in zip(PROFILE_C, COLORS):
            cells, _, errors = mesh_convergence(build_medium(1.0, 1.0, c, approximation))
            axes[0][j].loglog(cells, errors, marker='o', color=color, label=f'$c = {c}$')
        axes[0][j].loglog(cells, errors[0] * (cells[0] / cells)**2, 'k--',
                          alpha=0.6, label=r'$\propto N^{-2}$')
        panel(axes[0][j], f'{LABELS[approximation]} diffusion', 'Number of cells $N$',
              'Relative error of $R_c$', legend='lower left')

    finish(fig, 'Question 4: mesh convergence of the critical radius')
    savefig(fig, save_path)

def plot_flux_profiles(save_path):
    """Converged critical flux against the fundamental mode sin(Br)/Br, normalised at r = 0."""
    fig, axes = subplots(1, 2, height=5.0)

    for j, approximation in enumerate(APPROXIMATIONS):
        for c, color in zip(PROFILE_C, COLORS):
            medium = build_medium(1.0, 1.0, c, approximation)
            R = critical_radius(medium, n_cells=N_CELLS)
            _, r, phi = k_eigenvalue(R, medium, n_cells=N_CELLS)

            # np.sinc(y) = sin(pi y) / (pi y), which supplies the r = 0 limit.
            shape = np.sinc(buckling(medium) * r / np.pi)
            axes[0][j].plot(r / R, phi / phi[0], color=color, linewidth=2.2,
                            label=f'$c = {c}$, numerical')
            axes[0][j].plot(r / R, shape / shape[0], color=color, linestyle='--',
                            linewidth=1.4, alpha=0.8, label=f'$c = {c}$, $\\sin(Br)/Br$')

        axes[0][j].axvline(1.0, color='grey', linestyle=':', linewidth=1.2)
        panel(axes[0][j], f'{LABELS[approximation]} diffusion',
              '$r / R_c$   (the mesh runs out to $R_c + z_0$)',
              r'$\phi(r) / \phi(0)$', legend='lower left', fontsize=7)

    finish(fig, 'Question 4: critical flux shape vs. the fundamental mode')
    savefig(fig, save_path)

def report(figs):
    """Prints the radius table and the mesh-refinement check, and writes the figures."""
    log_section('Homework 1 Question 4',
                f'k by Bell & Glasstone source iteration on {N_CELLS} cells, critical',
                'radius by bisection on k(R) - 1. Lengths in mean free paths.')

    rows = []
    for approximation in APPROXIMATIONS:
        numerical, analytic = _radii(C_VALUES, approximation)
        for c, R_num, R_ana in zip(C_VALUES, numerical, analytic):
            rows.append([f'{c:.2f}', approximation, f'{R_num:.8f}', f'{R_ana:.8f}',
                         f'{(R_num - R_ana) / R_ana * 100.0:+.2e}'])
    log_table(['c', 'approximation', 'R_c numerical', 'R_c analytic', 'diff %'], rows)

    cells, _, errors = mesh_convergence(build_medium(1.0, 1.0, 1.5, 'classical'))
    ratios = errors[:-1] / errors[1:]
    logger.info(f"Mesh refinement at c = 1.5: error {errors[0]:.2e} -> {errors[-1]:.2e} "
                f"over N = {cells[0]} -> {cells[-1]}, ratio per doubling "
                f"{ratios.min():.2f}-{ratios.max():.2f}")

    plot_critical_radius(os.path.join(figs, 'q4_critical_radius.pdf'))
    plot_mesh_convergence(os.path.join(figs, 'q4_mesh_convergence.pdf'))
    plot_flux_profiles(os.path.join(figs, 'q4_flux_profiles.pdf'))
