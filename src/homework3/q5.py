"""Question 5: critical radius and mass of U-235 and Pu-239 by S_N."""

import os
import logging
import numpy as np
from homework1.criticality import critical_dimensions
from homework1.materials import BENCHMARK, FISSILE, critical_mass, solve_material
from homework1.tables import log_section, log_table
from homework3.figures import subplots, panel, finish, savefig
from homework3.sn import Medium
from homework3 import sphere

logger = logging.getLogger(__name__)

ORDERS = (2, 4, 8, 16)
CURVE_POINTS = 13
COLORS = {'Pu-239': '#c0392b', 'U-235': '#2c3e50'}

def sn_medium(material):
    """The benchmark cross-section row as an S_N medium."""
    return Medium(material.sigma_t, material.sigma_s, material.nu_sigma_f)

def radius_guess(material):
    """Assignment 1's exact-transport critical radius, converted from mfp to cm."""
    return critical_dimensions(material.c, 'transport-ref')[3] / material.sigma_t

def critical_radii(material):
    """{N: critical radius in cm} of one material, over the orders of the scan."""
    return {N: sphere.critical_radius(sn_medium(material), N, radius_guess(material))
            for N in ORDERS}

def _table(radii):
    """Critical radius and mass of every (material, N), against Assignment 1's diffusion."""
    rows = []
    for name in FISSILE:
        material = BENCHMARK[name]
        diffusion = solve_material(material, 'asymptotic')
        for N in ORDERS:
            R = radii[name][N]
            rows.append([name, f'S{N}', f'{material.c:.3f}', f'{R:.4f}',
                         f'{critical_mass(R, material.density):.3f}',
                         f'{diffusion.R_numerical:.4f}',
                         f'{R / diffusion.R_numerical - 1.0:+.2%}'])
    log_table(['material', 'order', 'c', 'R_c [cm]', 'M_c [kg]',
               'asympt. diffusion R_c [cm]', 'departure'], rows)

def plot_criticality(radii, save_path):
    """k against sphere radius at the highest order, with the k = 1 crossing."""
    fig, axes = subplots(1, len(FISSILE), height=5.0)
    order = ORDERS[-1]

    for j, name in enumerate(FISSILE):
        material, ax = BENCHMARK[name], axes[0][j]
        R_c = radii[name][order]
        grid = np.linspace(0.6 * R_c, 1.5 * R_c, CURVE_POINTS)

        ax.plot(grid, [sphere.k_eigenvalue(R, sn_medium(material), order).k for R in grid],
                color=COLORS[name], linewidth=2.2, label=f'$S_{{{order}}}$')
        ax.plot([R_c], [1.0], color=COLORS[name], marker='o', markersize=7,
                markerfacecolor='none', markeredgewidth=1.8)
        ax.axhline(1.0, color='grey', linestyle=':', linewidth=1.4)
        ax.annotate(f'$R_c = {R_c:.3f}$ cm,  '
                    f'$M_c = {critical_mass(R_c, material.density):.2f}$ kg',
                    xy=(0.03, 0.90), xycoords='axes fraction', fontsize=9,
                    color=COLORS[name])
        panel(ax, f'{name}   ($c = {material.c:.2f}$, '
                  rf'$\rho = {material.density}$ g/cm$^3$)',
              'Sphere radius $R$ [cm]', 'Multiplication factor $k$',
              legend='lower right', fontsize=9)

    finish(fig, 'Question 5: criticality of the bare benchmark spheres by $S_N$')
    savefig(fig, save_path)

def report(figs):
    """Prints the Question 5 table and writes its figure."""
    log_section('Assignment 3 Question 5',
                'Critical radius and mass of the fissile benchmark rows by S_N, with the '
                'scattering (inner) iteration now doing real work: Sigma_s > 0.')

    radii = {name: critical_radii(BENCHMARK[name]) for name in FISSILE}
    _table(radii)
    logger.info("Diffusion overestimates both radii: it cannot reproduce the forward-peaked "
                "angular flux at a surface only two mean free paths from the centre.")

    plot_criticality(radii, os.path.join(figs, 'q5_critical_mass.pdf'))
