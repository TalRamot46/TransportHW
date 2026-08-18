"""Question 2: critical slab half-thickness by P_N, for c = 1.2, 1.5, 1.8."""

import os
import logging
from homework1.tables import log_section, log_table
from homework3 import plots, pn_box, pn_modal
from homework3.sn import multiplying_medium

logger = logging.getLogger(__name__)

HALF_THICKNESS = 2          # index of a/2 in critical_dimensions
ORDERS = (1, 3, 5, 9)

def _table(scan):
    """Both methods at every (c, N), against Assignment 1's exact-transport value."""
    rows = []
    for c in plots.C_VALUES:
        reference = plots.reference_size(c, HALF_THICKNESS)
        medium = multiplying_medium(c)
        for N in ORDERS:
            box = scan.sizes[c][N]
            modal = pn_modal.critical_half_thickness(medium, N)
            rows.append([f'{c}', f'P{N}', f'{box:.6f}', f'{modal:.6f}',
                         f'{abs(box / modal - 1.0):.1e}', f'{reference:.5f}',
                         f'{modal / reference - 1.0:+.2%}'])
    log_table(['c', 'order', 'method 1 [mfp]', 'method 2 [mfp]', 'gap',
               'exact transport [mfp]', 'departure'], rows)

def report(figs):
    """Prints the Question 2 table and writes its figure."""
    log_section('Assignment 3 Question 2',
                'Critical slab half-thickness by P_N: Marshak conditions at the vacuum '
                'face, two independent solutions.',
                f'Method 1, box scheme on {pn_box.N_CELLS} cells with the k power '
                'iteration; Method 2, modal elimination and det H(a) = 0.',
                'Reflective at x = 0, vacuum at x = a.')

    scan = plots.scan_orders(pn_box.critical_half_thickness, pn_box.pn_k_eigenvalue,
                             HALF_THICKNESS, orders=ORDERS, family='P')
    _table(scan)

    plots.plot_orders(scan, '$a$',
                      'Question 2: critical slab half-thickness by $P_N$',
                      os.path.join(figs, 'q2_pn_slab_orders.pdf'))
