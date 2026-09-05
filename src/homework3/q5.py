"""Question 5: critical radius and mass of U-235 and Pu-239 by S_N."""

import os
import logging
import numpy as np
from homework1.criticality import critical_dimensions
from homework1.materials import BENCHMARK, FISSILE, critical_mass
from homework1.spherical import build_medium, analytic_critical_radius
from homework1.tables import log_section, log_table
from homework3.figures import subplots, panel, finish, savefig
from homework3.sn.core import Medium
from homework3.sn import sphere

logger = logging.getLogger(__name__)

# The ladder is doubled to S_64 so that the last two rungs settle the radius to better
# than 1e-3 cm; the table prints the change so the reader can see that, not take it.
ORDERS = (2, 4, 8, 16, 32, 64)
CONVERGED = ORDERS[-1]
PLOT_ORDER = 16                # the k(R) curve costs one solve per point; S_16 is enough
CURVE_POINTS = 13
COLORS = {'Pu-239': '#c0392b', 'U-235': '#2c3e50'}

def sn_medium(material):
    """The benchmark cross-section row as an S_N medium."""
    return Medium(material.sigma_t, material.sigma_s, material.nu_sigma_f)

def radius_guess(material):
    """The asymptotic + Milne critical radius in cm, used to start every search."""
    return critical_dimensions(material.c, 'transport-ref')[3] / material.sigma_t

def diffusion_radius(material, approximation):
    """Critical radius in cm from R_c = pi/B - z0 under 'classical' or 'asymptotic'."""
    return analytic_critical_radius(build_medium(
        material.sigma_t, material.sigma_a, material.nu_sigma_f, approximation))

def references(material):
    """{label: critical radius in cm} of the three closed forms, in decreasing crudeness."""
    return {'classical diffusion': diffusion_radius(material, 'classical'),
            'asymptotic diffusion': diffusion_radius(material, 'asymptotic'),
            'asymptotic + Milne': radius_guess(material)}

def critical_radii(material):
    """{N: critical radius in cm} of one material, over the orders of the ladder."""
    return {N: sphere.critical_radius(sn_medium(material), N, radius_guess(material))
            for N in ORDERS}

def _convergence_table(radii):
    """The S_N ladder, with the change per rung: this is what fixes the converged order."""
    rows = []
    for name in FISSILE:
        material, previous = BENCHMARK[name], None
        for N in ORDERS:
            R = radii[name][N]
            rows.append([name, f'S{N}', f'{R:.5f}',
                         '' if previous is None else f'{R - previous:+.5f}',
                         f'{critical_mass(R, material.density):.4f}'])
            previous = R
    log_table(['material', 'order', 'R_c [cm]', 'change [cm]', 'M_c [kg]'], rows)

def _comparison_table(radii):
    """Converged S_N against the three closed forms, radius and mass."""
    rows = []
    for name in FISSILE:
        material = BENCHMARK[name]
        converged = radii[name][CONVERGED]
        entries = [(f'converged S_{CONVERGED}', converged)] + list(references(material).items())
        for label, R in entries:
            # Four decimals, because the asymptotic + Milne row is inside 1e-5 of S_N and
            # two would print it as a flat zero, hiding the whole point of the table.
            rows.append([name, f'{material.c:.3f}', label, f'{R:.4f}',
                         f'{critical_mass(R, material.density):.3f}',
                         f'{R / converged - 1.0:+.4%}'])
    log_table(['material', 'c', 'method', 'R_c [cm]', 'M_c [kg]', 'vs converged S_N'], rows)

def plot_criticality(radii, save_path):
    """k against sphere radius at PLOT_ORDER, with the k = 1 crossing."""
    fig, axes = subplots(1, len(FISSILE), width=3.4, height=4.0)

    for j, name in enumerate(FISSILE):
        material, ax = BENCHMARK[name], axes[0][j]
        R_c = radii[name][PLOT_ORDER]
        grid = np.linspace(0.6 * R_c, 1.5 * R_c, CURVE_POINTS)

        ax.plot(grid, [sphere.k_eigenvalue(R, sn_medium(material), PLOT_ORDER).k for R in grid],
                color=COLORS[name], linewidth=2.2, label=f'$S_{{{PLOT_ORDER}}}$')
        ax.plot([R_c], [1.0], color=COLORS[name], marker='o', markersize=7,
                markerfacecolor='none', markeredgewidth=1.8)
        ax.axhline(1.0, color='grey', linestyle=':', linewidth=1.4)
        # The material and its parameters ride in the annotation, not the axis label:
        # at half the text width a label carrying them overruns the panel.
        ax.annotate(rf'{name},  $c = {material.c:.2f}$' '\n'
                    rf'$R_c = {R_c:.3f}$ cm' '\n'
                    rf'$M_c = {critical_mass(R_c, material.density):.2f}$ kg',
                    xy=(0.04, 0.96), xycoords='axes fraction', fontsize=11,
                    va='top', color=COLORS[name])
        panel(ax, 'Sphere radius $R$ [cm]', 'Multiplication factor $k$',
              legend='lower right')

    finish(fig)
    savefig(fig, save_path)

def report(figs):
    """Prints the Question 5 tables and writes its figure."""
    log_section('Assignment 3 Question 5',
                'Critical radius and mass of the fissile benchmark rows by S_N, with the '
                'scattering (inner) iteration now doing real work: Sigma_s > 0.',
                'Converged S_N against classical diffusion, asymptotic diffusion and the '
                'asymptotic + Milne relation R_c = (pi nu0 - z0)/Sigma_t.')

    radii = {name: critical_radii(BENCHMARK[name]) for name in FISSILE}
    _convergence_table(radii)
    _comparison_table(radii)

    logger.info("Classical diffusion is 9-12% high and asymptotic diffusion 0.1-0.2%. The "
                "asymptotic + Milne relation agrees with converged S_N to better than 0.01% "
                "on both materials -- it is not a transport solve, but for a bare sphere it "
                "is as good as one. Note that it and asymptotic diffusion are the same "
                "relation, R = pi nu0 - z0; they differ only in where z0 comes from. "
                "See report Question 5.")

    plot_criticality(radii, os.path.join(figs, 'q5_critical_mass.pdf'))
