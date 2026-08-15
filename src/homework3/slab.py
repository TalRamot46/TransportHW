"""Question 3: S_N in a slab, reflective at x = 0 and vacuum at x = a."""

import numpy as np
from homework3 import sn

N_CELLS = 200

class SlabSolver:
    """One transport sweep over the half-slab [0, a], on N_CELLS diamond-difference cells."""

    def __init__(self, half_thickness, medium, n_ordinates, n_cells=N_CELLS):
        self.medium = medium
        self.mu, self.weights = sn.ordinates(n_ordinates)
        self.n_cells = n_cells
        self.dx = half_thickness / n_cells
        self.centres = (np.arange(n_cells) + 0.5) * self.dx
        self.volumes = np.full(n_cells, self.dx)

    def _sweep(self, m, source, incoming):
        """One sweep: marches ordinate m across the mesh, returning its contribution
        to the scalar flux."""
        mu = self.mu[m]
        inward = mu < 0.0
        cells = range(self.n_cells - 1, -1, -1) if inward else range(self.n_cells)

        # Vacuum at x = a for the inward half, the reflected flux at x = 0 for the outward one.
        psi_in = 0.0 if inward else incoming[m]
        psi = np.empty(self.n_cells)
        for i in cells:
            # One face, the spatial one, weighted by |mu| on both sides: report eq. (3).
            psi[i], (psi_in,) = sn.cell_flux(self.medium.sigma_t * self.dx, source[i] * self.dx,
                                             [sn.Face(abs(mu), abs(mu), psi_in)])
        if inward:
            incoming[len(self.mu) - 1 - m] = psi_in   # reflective boundary at x = 0
        return self.weights[m] * psi

    def sn_iteration(self, source):
        """One S_N iteration: one sweep per ordinate, summed into the scalar flux.
        The source S is isotropic, so each ordinate carries S/2."""
        q = 0.5 * source
        incoming = np.zeros(len(self.mu))
        # Ascending mu, so every inward ordinate is swept before the outward one it feeds.
        return sum(self._sweep(m, q, incoming) for m in range(len(self.mu)))

def slab_k_eigenvalue(half_thickness, medium, n_ordinates, n_cells=N_CELLS):
    """KResult of a slab of the given half-thickness, in mean free paths."""
    return sn.k_eigenvalue(SlabSolver(half_thickness, medium, n_ordinates, n_cells))

def critical_half_thickness(medium, n_ordinates, guess, n_cells=N_CELLS):
    """Half-thickness at which k = 1."""
    return sn.critical_size(
        lambda a: slab_k_eigenvalue(a, medium, n_ordinates, n_cells).k, guess)
