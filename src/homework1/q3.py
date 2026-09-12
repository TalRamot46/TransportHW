"""Question 3: critical dimensions from the transport relations and from diffusion."""

import os
import logging
import numpy as np
from homework1.criticality import (critical_dimensions, critical_dimensions_applied_bc,
                                   extrapolation_distance, METHOD_LABELS,
                                   MARSHAK_EXTRAPOLATION, MARK_EXTRAPOLATION,
                                   CASE_TABLE_8_C, CASE_TABLE_8_K0,
                                   CASE_TABLE_23_C, CASE_TABLE_23_CZ0)
from homework1.figures import (subplots, panel, finish, savefig,
                               NAVY, ORANGE, GREEN, RED, BLUE, PURPLE)
from homework1.tables import log_section, log_table

logger = logging.getLogger(__name__)

C_VALUES = (1.02, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0)
C_GRID = np.linspace(1.02, 2.0, 400)

HALF, RADIUS = 2, 3      # indices of a/2 and Sigma_t R_c in the critical_dimensions tuple
# Index, y-label of the dimension, y-label of its departure panel, and where that
# panel's legend clears its curves.
DIMENSIONS = ((HALF, r'$a/2$ [mfp]', r'Departure of $a/2$ (%)', 'lower left'),
              (RADIUS, r'$\Sigma_t R_c$ [mfp]',
               r'Departure of $\Sigma_t R_c$ (%)', 'upper left'))

METHOD_STYLES = {'transport': (NAVY, '-'),
                 'marshak': (ORANGE, '-'),
                 'mark': (GREEN, '-')}
APPLIED_STYLES = {'marshak': (ORANGE, ':'), 'mark': (GREEN, ':')}

def _draw(ax, c, curves):
    """Plots a list of (values, label, color, linestyle) on one axis."""
    for values, label, color, style in curves:
        ax.plot(c, values, label=label, color=color, linestyle=style, linewidth=2.0)

def _error(values, reference):
    """Relative error in per cent."""
    return np.abs(values - reference) / reference * 100.0

def plot_critical_dimensions(save_path):
    """The two dimensions from the approximate inputs, against a tabulated reference."""
    fit = critical_dimensions(C_GRID, 'transport')
    ref = critical_dimensions(C_GRID, 'transport-ref')
    alt = critical_dimensions(C_GRID, 'transport-q+')

    fig, axes = subplots(2, 2)

    for column, (index, ylabel, _, _legend) in enumerate(DIMENSIONS):
        _draw(axes[0][column], C_GRID, [
            (ref[index], r'exact $|\nu_0|$, tabulated $z_0$', NAVY, '-'),
            (fit[index], r'fit, $q = -0.0199$', RED, '--'),
            (alt[index], r'fit, $q = +0.0199$', BLUE, ':')])
        panel(axes[0][column], '$c$', ylabel, log=True, legend='upper right')

    _draw(axes[1][0], C_GRID, [
        (_error(fit[0], ref[0]), r'$|\nu_0|$: fit vs. root', PURPLE, '-'),
        (_error(fit[1], ref[1]), r'$z_0$: $q = -0.0199$ vs. Table 23', RED, '--'),
        (_error(extrapolation_distance(C_GRID, 0.0199), ref[1]),
         r'$z_0$: $q = +0.0199$ vs. Table 23', BLUE, ':')])
    panel(axes[1][0], '$c$', 'Relative error (%)', log=True, legend='lower right')

    _draw(axes[1][1], C_GRID, [
        (_error(fit[HALF], ref[HALF]), r'$a/2$, $q = -0.0199$', RED, '--'),
        (_error(fit[RADIUS], ref[RADIUS]), r'$\Sigma_t R_c$, $q = -0.0199$', ORANGE, '-'),
        (_error(alt[HALF], ref[HALF]), r'$a/2$, $q = +0.0199$', BLUE, ':'),
        (_error(alt[RADIUS], ref[RADIUS]), r'$\Sigma_t R_c$, $q = +0.0199$', GREEN, '-.')])
    panel(axes[1][1], '$c$', 'Relative error (%)', log=True, legend='upper left')

    finish(fig)
    savefig(fig, save_path)

def plot_method_comparison(save_path):
    """Parts 3(a)-3(c) against each other, with the applied-BC variants dotted."""
    dimensions = {m: critical_dimensions(C_GRID, m) for m in METHOD_STYLES}
    applied = {'marshak': critical_dimensions_applied_bc(C_GRID, MARSHAK_EXTRAPOLATION),
               'mark': critical_dimensions_applied_bc(C_GRID, MARK_EXTRAPOLATION)}

    fig, axes = subplots(2, 2)

    for column, (index, ylabel, departure, legend) in enumerate(DIMENSIONS):
        _draw(axes[0][column], C_GRID,
              [(dimensions[m][index], METHOD_LABELS[m], *METHOD_STYLES[m])
               for m in METHOD_STYLES]
              + [(applied[m][column], f'{METHOD_LABELS[m]}, BC applied', *APPLIED_STYLES[m])
                 for m in APPLIED_STYLES])
        panel(axes[0][column], '$c$', ylabel, log=True, legend='upper right')

        # Signed departure from transport, so over- and under-estimation are distinct.
        reference = dimensions['transport'][index]
        _draw(axes[1][column], C_GRID,
              [((dimensions[m][index] - reference) / reference * 100.0,
                METHOD_LABELS[m], *METHOD_STYLES[m]) for m in APPLIED_STYLES]
              + [((applied[m][column] - reference) / reference * 100.0,
                  f'{METHOD_LABELS[m]}, BC applied', *APPLIED_STYLES[m])
                 for m in APPLIED_STYLES])
        axes[1][column].axhline(0.0, color=NAVY, linewidth=1.2)
        panel(axes[1][column], '$c$', departure, legend=legend)

    finish(fig)
    savefig(fig, save_path)

def plot_extrapolation_distance(save_path):
    """The tabulated product c z0(c) against both signs of the quadratic term."""
    # The expansion is about c = 1 and is not meant to hold as c -> 0, where the
    # tabulated product climbs to 1; this window keeps the near-critical structure legible.
    c = np.linspace(0.5, CASE_TABLE_23_C[-1], 600)

    fig, axes = subplots(1, 1, width=5.3, height=3.5)
    ax = axes[0][0]
    ax.plot(CASE_TABLE_23_C, CASE_TABLE_23_CZ0, 'o', color=NAVY,
            markersize=5, label='Case et al., Table 23')
    _draw(ax, c, [(c * extrapolation_distance(c, -0.0199), 'fit, $q = -0.0199$', RED, '--'),
                  (c * extrapolation_distance(c, +0.0199), 'fit, $q = +0.0199$', BLUE, ':')])
    ax.axhline(0.710446, color='#7f8c8d', linewidth=1.0, alpha=0.8)
    ax.set_xlim(0.5, 3.0)
    ax.set_ylim(0.700, 0.726)
    panel(ax, '$c$', r'$c\,z_0(c)$', legend='upper right')

    finish(fig)
    savefig(fig, save_path)

def _transport_table():
    """The part 3(a) table: approximate inputs and dimensions against the reference."""
    rows = []
    for c in C_VALUES:
        nu_fit, z0_fit, half_fit, radius_fit = critical_dimensions(c, 'transport')
        nu_ref, z0_ref, half_ref, _ = critical_dimensions(c, 'transport-ref')
        rows.append([f'{c:.2f}', nu_fit, nu_ref, z0_fit, z0_ref, half_fit, radius_fit,
                     f'{_error(half_fit, half_ref):.4f}'])
    log_table(['c', '|nu0| fit', '|nu0| ref', 'z0 fit', 'z0 ref', 'a/2', 'R_c',
               'a/2 err %'], rows)

def _method_table():
    """The parts 3(a)-3(c) comparison table."""
    rows = []
    for c in C_VALUES:
        values = [critical_dimensions(c, m) for m in ('transport', 'marshak', 'mark')]
        rows.append([f'{c:.2f}'] + [v[HALF] for v in values] + [v[RADIUS] for v in values])
    log_table(['c', 'a/2 transport', 'a/2 Marshak', 'a/2 Mark',
               'R transport', 'R Marshak', 'R Mark'], rows)

def report(figs):
    """Prints both Question 3 tables with their checks, and writes the figures."""
    log_section('Homework 1 Question 3',
                'a/2 = (pi/2)|nu0| - z0,  Sigma_t R_c = pi |nu0| - z0, in mean free paths.',
                "'fit' uses the approximate |nu0| and z0; 'ref' uses the transcendental",
                "root and Case's Table 23.")
    _transport_table()

    k0_error = max(abs(1.0 / critical_dimensions(c, 'transport-ref')[0] - k0)
                   for c, k0 in zip(CASE_TABLE_8_C[1:], CASE_TABLE_8_K0[1:]))
    logger.info(f"Eigenvalue check: max |k0 - Case Table 8| = {k0_error:.2e}")

    # The printed z0 fit carries a minus sign on the quadratic term, which has the
    # wrong sign against Table 23; both are reported. See report §3.
    for q in (-0.0199, 0.0199):
        error = max(_error(extrapolation_distance(c, q),
                           critical_dimensions(c, 'transport-ref')[1]) for c in C_VALUES)
        logger.info(f"z0 fit with q = {q:+.4f}: max error vs. Table 23 = {error:.4f} %")

    log_section('Homework 1 Question 3(b), 3(c)',
                'Diffusion with Marshak (z0 = 2/3) and Mark (z0 = 1/sqrt(3)).')
    _method_table()

    for name, l0 in (('Marshak', MARSHAK_EXTRAPOLATION), ('Mark', MARK_EXTRAPOLATION)):
        for c in (1.02, 2.0):
            half_bc, radius_bc = critical_dimensions_applied_bc(c, l0)
            _, _, half, radius = critical_dimensions(c, name.lower())
            logger.info(f"{name} at c = {c}: a/2 = {half:.6f} extrapolated vs. "
                        f"{half_bc:.6f} with the BC applied; R_c = {radius:.6f} vs. {radius_bc:.6f}")

    plot_critical_dimensions(os.path.join(figs, 'q3_critical_dimensions.pdf'))
    plot_method_comparison(os.path.join(figs, 'q3_method_comparison.pdf'))
    plot_extrapolation_distance(os.path.join(figs, 'q3_extrapolation_distance.pdf'))
