"""Question 1: critical radius of a reflected sphere, in three diffusion approximations."""

import os
import logging

from homework1.materials import BENCHMARK, FISSILE, critical_mass
from homework1.spherical import build_medium, analytic_critical_radius
from homework1.tables import log_section, log_table
from homework3.figures import subplots, panel, finish, savefig
from homework3 import reflected

logger = logging.getLogger(__name__)

REFLECTORS = ('H2O', 'Fe', 'Na')
REFLECTOR_LABELS = {'H2O': 'Water', 'Fe': 'Iron', 'Na': 'Sodium'}
DEPTHS = (1, 2, 3, 10)
FLUX_DEPTH = 3                 # the reflector thickness whose flux profiles are drawn

COLORS = {'classic': '#2c3e50', 'asymptotic': '#e67e22', 'zimmerman': '#c0392b'}

# The bare limit each theory is measured against. Both asymptotic theories share
# Assignment 1's asymptotic sphere: with no reflector there is no interface, so no jump.
BARE = {'classic': 'classical', 'asymptotic': 'asymptotic', 'zimmerman': 'asymptotic'}

def bare_radius(material, theory):
    """Assignment 1's bare critical radius in cm, under the matching approximation."""
    return analytic_critical_radius(
        build_medium(material.sigma_t, material.sigma_a, material.nu_sigma_f, BARE[theory]))

def radii(core_name):
    """{(reflector, d, theory): critical core radius in cm} of one core."""
    core = BENCHMARK[core_name]
    return {(name, d, theory): reflected.critical_radius(core, BENCHMARK[name], d, theory)
            for name in REFLECTORS for d in DEPTHS for theory in reflected.THEORIES}

def _parameter_table():
    """The asymptotic transport parameters of every material in play."""
    rows = []
    for name in FISSILE + REFLECTORS:
        region = reflected.region(BENCHMARK[name], 'asymptotic')
        rows.append([name, f'{region.c:.4f}', f'{1.0 / region.sigma_t:.3f}',
                     f'{region.rate:.5f}', f'{region.D0:.5f}',
                     f'{region.D0 / region.sigma_t:.4f}',
                     f'{region.mu0:.5f}', f'{region.z0:.5f}'])
    log_table(['material', 'c', 'mfp [cm]', 'k0 or 1/nu0', 'D0', 'D [cm]', 'mu0', 'z0'], rows)

def _bare_table():
    """Assignment 1's bare spheres, the reference for every reflected radius below."""
    rows = []
    for name in FISSILE:
        material = BENCHMARK[name]
        R = [bare_radius(material, theory) for theory in ('classic', 'asymptotic')]
        rows.append([name] + [f'{v:.4f}' for v in R]
                    + [f'{critical_mass(v, material.density):.3f}' for v in R])
    log_table(['core', 'classic R_c [cm]', 'asympt R_c [cm]',
               'classic M_c [kg]', 'asympt M_c [kg]'], rows)

def _core_table(core_name, table):
    """Critical radius and mass of one core behind every reflector, thickness and theory."""
    core = BENCHMARK[core_name]
    rows = []
    for name in REFLECTORS:
        for d in DEPTHS:
            R = [table[(name, d, theory)] for theory in reflected.THEORIES]
            rows.append([REFLECTOR_LABELS[name], f'{d}']
                        + [f'{v:.4f}' for v in R]
                        + [f'{critical_mass(v, core.density):.3f}' for v in R]
                        + [f'{R[2] / bare_radius(core, "zimmerman") - 1.0:+.1%}'])
    log_table(['reflector', 'd [mfp]', 'R_c (a) [cm]', 'R_c (b) [cm]', 'R_c (c) [cm]',
               'M_c (a) [kg]', 'M_c (b) [kg]', 'M_c (c) [kg]', '(c) vs bare'], rows)

def plot_radii(tables, save_path):
    """Critical radius against reflector thickness, one panel per core and reflector."""
    fig, axes = subplots(len(FISSILE), len(REFLECTORS), width=5.0, height=4.0)

    for i, core_name in enumerate(FISSILE):
        core = BENCHMARK[core_name]
        for j, name in enumerate(REFLECTORS):
            ax = axes[i][j]
            for theory in reflected.THEORIES:
                ax.plot(DEPTHS, [tables[core_name][(name, d, theory)] for d in DEPTHS],
                        'o-', color=COLORS[theory], linewidth=2.0,
                        label=reflected.THEORY_LABELS[theory])
            for theory, style in (('classic', '--'), ('asymptotic', ':')):
                ax.axhline(bare_radius(core, theory), color=COLORS[theory], linestyle=style,
                           linewidth=1.4, alpha=0.8,
                           label=f'bare, {reflected.THEORY_LABELS[theory].lower()}')
            panel(ax, f'{core_name} + {REFLECTOR_LABELS[name]}',
                  'Reflector thickness $d$ [mfp]', 'Critical core radius $R_c$ [cm]')

    finish(fig, 'Question 1: critical core radius of the reflected spheres')
    # One legend for all six panels: the curves leave no free corner inside them.
    fig.legend(*axes[0][0].get_legend_handles_labels(), loc='lower center',
               ncol=5, fontsize=9, frameon=True, bbox_to_anchor=(0.5, -0.03))
    savefig(fig, save_path)

def plot_fluxes(save_path):
    """Critical flux across core and reflector; the Zimmerman curve jumps at the interface."""
    fig, axes = subplots(1, len(REFLECTORS), width=5.0, height=4.2)
    core = BENCHMARK[FISSILE[0]]

    for j, name in enumerate(REFLECTORS):
        ax = axes[0][j]
        for theory in reflected.THEORIES:
            r, phi = reflected.flux_profile(core, BENCHMARK[name], FLUX_DEPTH, theory)
            R = reflected.critical_radius(core, BENCHMARK[name], FLUX_DEPTH, theory)
            ax.plot(r, phi, color=COLORS[theory], linewidth=2.0,
                    label=reflected.THEORY_LABELS[theory])
            ax.axvline(R, color=COLORS[theory], linestyle=':', linewidth=1.0, alpha=0.6)
        panel(ax, f'{FISSILE[0]} + {REFLECTOR_LABELS[name]},  $d = {FLUX_DEPTH}$ mfp',
              'Radius $r$ [cm]', r'$\phi(r) / \phi(0)$', legend='upper right', fontsize=8)

    finish(fig, f'Question 1: critical flux profiles behind {FLUX_DEPTH} mfp of reflector')
    savefig(fig, save_path)

def report(figs):
    """Prints the Question 1 tables and writes its two figures."""
    log_section('Assignment 3 Question 1',
                'Critical radius of a reflected sphere, from the two-region diffusion '
                'criticality condition in (a) classic, (b) asymptotic and (c) discontinuous '
                'asymptotic (Zimmerman) form.',
                'Flux continuity at the interface in (a) and (b); mu0 phi continuous in (c). '
                'The net current -D dphi/dr is continuous in all three.')

    _parameter_table()
    _bare_table()

    tables = {name: radii(name) for name in FISSILE}
    for name in FISSILE:
        logger.info(f'\n{name} core:')
        _core_table(name, tables[name])

    logger.info("\nWater and iron behave: every radius falls below the bare one and keeps "
                "falling with d. Sodium does not -- its radii sit ABOVE the bare sphere. "
                "That is diffusion theory failing, not the reflector: sodium's mean free "
                "path is 11.58 cm and its D is 3.86 cm, against a core radius near 6 cm, "
                "so the spherical term D_R/R alone leaks more than a vacuum. "
                "See report §1.")

    plot_radii(tables, os.path.join(figs, 'q1_critical_radius.pdf'))
    plot_fluxes(os.path.join(figs, 'q1_flux_profiles.pdf'))
