"""Question 3: S_N in a slab, reflective at x = 0 and vacuum at x = a/2."""

import numpy as np
from homework3.sn import core

N_CELLS = 200

def cell_flux(removal, source, mu_abs, psi_in):
    """Cell-centre and outgoing flux of one slab cell, with the negative-flux fixup.
    Report eqs. (22) forward and (23) backward."""
    psi = (source + 2.0 * mu_abs * psi_in) / (removal + 2.0 * mu_abs)
    psi_out = 2.0 * psi - psi_in

    if psi_out < 0.0:
        # Hold the face at zero and re-solve: the outgoing term leaves the balance
        # and only the inflow remains, which lowers psi -- see report Fixup.
        psi_out = 0.0
        psi = (source + mu_abs * psi_in) / removal
    return psi, psi_out

class SlabSolver:
    """One transport sweep over the half-slab [0, a/2], on N_CELLS diamond cells."""

    def __init__(self, half_thickness, medium, n_ordinates, n_cells=N_CELLS):
        self.medium = medium
        self.mu, self.weights = core.ordinates(n_ordinates)
        self.n_cells = n_cells
        self.dx = half_thickness / n_cells
        self.centres = (np.arange(n_cells) + 0.5) * self.dx
        self.volumes = np.full(n_cells, self.dx)

    def _cell(self, i, m, q, psi_in):
        """One cell of a march, in the units the balance is written in."""
        return cell_flux(self.medium.sigma_t * self.dx, q[i] * self.dx,
                         abs(self.mu[m]), psi_in)

    def sweep_backward(self, m, q):
        """mu_m < 0: marches inwards from the vacuum face at x = a/2 to x = 0. Returns
        the ordinate's flux and what it leaves at the reflective boundary."""
        psi, psi_in = np.empty(self.n_cells), 0.0          # vacuum: nothing enters
        for i in range(self.n_cells - 1, -1, -1):
            psi[i], psi_in = self._cell(i, m, q, psi_in)
        return psi, psi_in

    def sweep_forward(self, m, q, psi_in):
        """mu_m > 0: marches outwards from x = 0 to the vacuum face, starting from the
        flux its mirror ordinate -mu_m left at the reflective boundary."""
        psi = np.empty(self.n_cells)
        for i in range(self.n_cells):
            psi[i], psi_in = self._cell(i, m, q, psi_in)
        return psi

    def sweep_all_angles(self, source):
        """One S_N iteration: every backward ordinate, then every forward one. The
        source S is isotropic, so each ordinate carries S/2."""
        q = 0.5 * source
        phi, reflected = np.zeros(self.n_cells), np.zeros(len(self.mu))
        mirror = len(self.mu) - 1

        # Backward first: each -mu_m fills the reflective inflow of its mirror +mu_m.
        for m in np.flatnonzero(self.mu < 0.0):
            psi, reflected[mirror - m] = self.sweep_backward(m, q)
            phi += self.weights[m] * psi
        for m in np.flatnonzero(self.mu > 0.0):
            phi += self.weights[m] * self.sweep_forward(m, q, reflected[m])
        return phi

def k_eigenvalue(half_thickness, medium, n_ordinates, n_cells=N_CELLS):
    """KResult of a slab of the given half-thickness, in mean free paths."""
    return core.k_eigenvalue(SlabSolver(half_thickness, medium, n_ordinates, n_cells))

def critical_half_thickness(medium, n_ordinates, guess, n_cells=N_CELLS):
    """Half-thickness at which k = 1."""
    return core.critical_size(
        lambda a: k_eigenvalue(a, medium, n_ordinates, n_cells).k, guess)
