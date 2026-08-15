"""Question 5: bare critical radius and mass of the fissile benchmark materials."""

import os
import logging
import numpy as np
from homework1.materials import BENCHMARK, FISSILE, PROMPT_U235, critical_mass, solve_material
from homework1.spherical import build_medium, critical_radius, k_eigenvalue
from homework1.figures import subplots, panel, finish, savefig
from homework1.tables import log_section, log_table

logger = logging.getLogger(__name__)

APPROXIMATIONS = ('classical', 'asymptotic')
LABELS = {'classical': 'Classical', 'asymptotic': 'Asymptotic'}
COLORS = {'classical': '#2c3e50', 'asymptotic': '#e74c3c'}
N_CELLS = 400

def plot_criticality(save_path):
    """k against sphere radius, with the k = 1 crossing that fixes the critical mass."""
    fig, axes = subplots(1, len(FISSILE), height=5.0)

    for j, name in enumerate(FISSILE):
        material = BENCHMARK[name]
        ax = axes[0][j]

        for i, approximation in enumerate(APPROXIMATIONS):
            medium = build_medium(material.sigma_t, material.sigma_a,
                                  material.nu_sigma_f, approximation)
            R_c = critical_radius(medium, n_cells=N_CELLS)
            radii = np.linspace(0.4 * R_c, 1.8 * R_c, 40)

            ax.plot(radii, [k_eigenvalue(R, medium, n_cells=N_CELLS).k for R in radii],
                    color=COLORS[approximation], linewidth=2.2,
                    label=LABELS[approximation])
            ax.plot([R_c], [1.0], color=COLORS[approximation], marker='o', markersize=7,
                    markerfacecolor='none', markeredgewidth=1.8)
            ax.annotate(f'$R_c = {R_c:.3f}$ cm,  '
                        f'$M_c = {critical_mass(R_c, material.density):.2f}$ kg',
                        xy=(0.03, 0.90 - 0.08 * i), xycoords='axes fraction',
                        fontsize=9, color=COLORS[approximation])

        ax.axhline(1.0, color='grey', linestyle=':', linewidth=1.4)
        panel(ax, f'{name}   ($c = {material.c:.2f}$, '
                  rf'$\rho = {material.density}$ g/cm$^3$)',
              'Sphere radius $R$ [cm]', 'Multiplication factor $k$',
              legend='lower right', fontsize=9)

    finish(fig, 'Question 5: criticality of the bare benchmark spheres')
    savefig(fig, save_path)

def _mass_table():
    """Critical radius and mass of each fissile row, both approximations."""
    rows = []
    for material in [BENCHMARK[name] for name in FISSILE] + [PROMPT_U235]:
        # Sigma_t must equal the sum of its parts, or the row has a typo.
        if abs(material.sigma_t_sum - material.sigma_t) > 1e-12:
            logger.warning(f"  {material.name}: Sigma_t is not the sum of its parts.")
        for approximation in APPROXIMATIONS:
            r = solve_material(material, approximation, n_cells=N_CELLS)
            rows.append([material.name, approximation, f'{material.c:.4f}',
                         f'{r.R_numerical:.4f}', f'{r.R_analytic:.4f}',
                         f'{r.mass_numerical:.3f}'])
    log_table(['material', 'approximation', 'c', 'R_c [cm]', 'R analytic [cm]',
               'M_c [kg]'], rows)

def report(figs):
    """Prints the critical radii and masses, and writes the Question 5 figure."""
    log_section('Homework 1 Question 5',
                'Bare critical sphere from the Sood et al. one-group benchmark data.')
    _mass_table()
    logger.info("The prompt variant is the U-235 row given in the task prompt rather than "
                "in the assignment PDF; its own cross sections give c = 1.365, not 1.50.")
    subcritical = ', '.join(name for name in BENCHMARK if name not in FISSILE)
    logger.info(f"No bare critical sphere exists for {subcritical}: c <= 1, non-multiplying.")

    plot_criticality(os.path.join(figs, 'q5_criticality.pdf'))
