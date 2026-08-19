"""Question 2, Method 1: the midpoint box discretisation of the P_N system, driven by
the Bell & Glasstone k iteration."""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from homework3 import pn, sn

N_CELLS = 200
HOPF_CONSTANT = 0.7104        # extrapolation distance of the starting flux guess

class BoxSystem:
    """The k-independent box matrix of report eq. (13) on [0, a/2], factorised once."""

    def __init__(self, half_thickness, medium, N, n_cells=N_CELLS):
        self.medium, self.width, self.n_cells = medium, N + 1, n_cells
        self.n_conditions = (N + 1) // 2
        self.dx = half_thickness / n_cells
        self.nodes = np.arange(n_cells + 1) * self.dx
        self.lu = splu(self._matrix(N).tocsc())

    def _matrix(self, N):
        """Symmetry rows, then one block per cell, then the Marshak rows."""
        width, first = self.width, self.n_conditions
        matrix = sparse.lil_matrix(((self.n_cells + 1) * width,) * 2)

        for row, n in enumerate(pn.parity_indices(N)[1]):
            matrix[row, n] = 1.0                                 # phi_odd(0) = 0

        collision = 0.5 * self.medium.sigma_t * np.eye(width)
        streaming = pn.streaming_matrix(N) / self.dx
        block = np.hstack([collision - streaming, collision + streaming])
        for j in range(self.n_cells):
            row = first + j * width
            matrix[row:row + width, j * width:(j + 2) * width] = block

        matrix[-first:, -width:] = pn.marshak_matrix(N)
        return matrix

    def solve(self, source):
        """Every moment at every node, given the cell-midpoint source of report eq. (13)."""
        rhs = np.zeros((self.n_cells + 1) * self.width)
        rhs[self.n_conditions::self.width][:self.n_cells] = source
        return self.lu.solve(rhs).reshape(self.n_cells + 1, self.width)

def _midpoints(phi):
    """Cell-midpoint values of a nodal profile: the average of its two faces."""
    return 0.5 * (phi[1:] + phi[:-1])

def pn_k_eigenvalue(half_thickness, medium, N, n_cells=N_CELLS, tol=1e-9, max_iter=2000):
    """KResult of the slab by Method 1, the power iteration of report eqs. (14)-(15)."""
    system = BoxSystem(half_thickness, medium, N, n_cells)
    production = medium.sigma_t * medium.c
    phi = np.cos(0.5 * np.pi * system.nodes / (half_thickness + HOPF_CONSTANT))
    k = 1.0

    for outer in range(1, max_iter + 1):
        source = _midpoints(phi)
        phi_new = system.solve(production / k * source)[:, 0]

        # The same midpoint average the source uses, so this is the ratio of the two
        # discrete fission integrals and not a second quadrature of them.
        k_new = k * _midpoints(phi_new).sum() / source.sum()
        converged = abs(k_new - k) < tol * abs(k_new)

        phi, k = phi_new / phi_new.max(), k_new
        if converged:
            return sn.KResult(k, system.nodes, phi, outer)

    raise RuntimeError(f"The P_N k iteration did not converge in {max_iter} outers.")

def critical_half_thickness(medium, N, guess, n_cells=N_CELLS):
    """Half-thickness at which k = 1, by brentq on Method 1's k."""
    return sn.critical_size(
        lambda a: pn_k_eigenvalue(a, medium, N, n_cells).k, guess)
