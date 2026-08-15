"""Question 4: critical sphere radius by S_N, for c = 1.2, 1.5, 1.8."""

import os
import logging
from homework1.tables import log_section, log_table
from homework3 import plots, sphere

logger = logging.getLogger(__name__)

RADIUS = 3          # index of Sigma_t R_c in critical_dimensions

def _table(scan):
    """Critical radius of every (c, N), against Assignment 1's exact-transport value."""
    rows = []
    for c in plots.C_VALUES:
        reference = plots.reference_size(c, RADIUS)
        for N in plots.ORDERS:
            radius = scan.sizes[c][N]
            rows.append([f'{c}', f'S{N}', f'{radius:.5f}', f'{reference:.5f}',
                         f'{radius / reference - 1.0:+.2%}'])
    log_table(['c', 'order', 'R_c [mfp]', 'exact transport [mfp]', 'departure'], rows)

def report(figs):
    """Prints the Question 4 table and writes its figure."""
    log_section('Assignment 3 Question 4',
                'Critical sphere radius by S_N, with the angular redistribution term '
                'differenced by the alpha recursion and a mu = -1 starting direction.',
                f'Reflective at r = 0, vacuum at r = R, {sphere.N_CELLS} cells.')

    scan = plots.scan_orders(sphere.critical_radius, sphere.k_eigenvalue, RADIUS)
    _table(scan)
    logger.info("The slab and the sphere share an extrapolation distance, so R_c and a/2 "
                "should satisfy R_c + z0 = 2 (a/2 + z0); see report §4.")

    plots.plot_orders(scan, '$R_c$',
                      'Question 4: critical sphere radius by $S_N$',
                      os.path.join(figs, 'q4_sphere_orders.pdf'))
