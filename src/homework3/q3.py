"""Question 3: critical slab half-thickness by S_N, for c = 1.2, 1.5, 1.8."""

import os
import logging
from homework1.tables import log_section, log_table
from homework3 import plots, slab

logger = logging.getLogger(__name__)

HALF_THICKNESS = 2          # index of a/2 in critical_dimensions

# Exact one-speed transport value at c = 1.5, from the Sood et al. slab benchmark
# whose Pu-239 half-thickness is 1.853722 cm at Sigma_t = 0.32640 cm^-1.
BENCHMARK_C15 = 0.605055

def _table(scan):
    """Critical half-thickness of every (c, N), against Assignment 1's exact-transport value."""
    rows = []
    for c in plots.C_VALUES:
        reference = plots.reference_size(c, HALF_THICKNESS)
        for N in plots.ORDERS:
            half = scan.sizes[c][N]
            rows.append([f'{c}', f'S{N}', f'{half:.5f}', f'{reference:.5f}',
                         f'{half / reference - 1.0:+.2%}'])
    log_table(['c', 'order', 'a/2 [mfp]', 'exact transport [mfp]', 'departure'], rows)

def report(figs):
    """Prints the Question 3 table and writes its figure."""
    log_section('Assignment 3 Question 3',
                'Critical slab half-thickness by S_N: Gauss-Legendre ordinates, diamond '
                'difference, negative-flux fixup.',
                f'Reflective at x = 0, vacuum at x = a/2, {slab.N_CELLS} cells.')

    scan = plots.scan_orders(slab.critical_half_thickness, slab.k_eigenvalue, HALF_THICKNESS)
    _table(scan)
    logger.info(f"At c = 1.5 the S_N sequence runs down towards the exact one-speed value "
                f"{BENCHMARK_C15:.6f} mfp; S10 is off by "
                f"{scan.sizes[1.5][10] / BENCHMARK_C15 - 1.0:+.2%}.")

    plots.plot_orders(scan, '$a/2$',
                      'Question 3: critical slab half-thickness by $S_N$',
                      os.path.join(figs, 'q3_slab_orders.pdf'))
