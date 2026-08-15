"""Questions 4 and 5: S_N in a sphere, reflective at r = 0 and vacuum at r = R."""

import numpy as np
from homework3 import sn

# 100 cells is already mesh converged: at c = 1.5 the S10 radius moves by 7e-6 mfp
# between 50 and 800 cells. See explanations/06.
N_CELLS = 100

def angular_coefficients(mu, weights):
    """Carlson's alpha_{m+1/2} = alpha_{m-1/2} - w_m mu_m, zero at both ends; the
    recursion is what makes a flat flux an exact solution. See explanations/02."""
    alpha = np.concatenate(([0.0], -np.cumsum(weights * mu)))
    alpha[-1] = 0.0   # exactly zero by sum(w mu) = 0; set so no current leaks at mu = +1
    return alpha

class SphereSolver:
    """One transport sweep over [0, R], with the angular redistribution term of
    spherical geometry differenced in angle as well as in space."""

    def __init__(self, radius, medium, n_ordinates, n_cells=N_CELLS):
        self.medium = medium
        self.mu, self.weights = sn.ordinates(n_ordinates)
        self.alpha = angular_coefficients(self.mu, self.weights)
        self.n_cells = n_cells
        self.dr = radius / n_cells

        faces = np.linspace(0.0, radius, n_cells + 1)
        self.centres = 0.5 * (faces[:-1] + faces[1:])
        self.areas = 4.0 * np.pi * faces**2
        self.volumes = (4.0 * np.pi / 3.0) * (faces[1:]**3 - faces[:-1]**3)

    def _starting_direction(self, source):
        """Half-angle flux at mu = -1, where the angular derivative drops out and the
        sweep is the plain slab one; see explanations/02."""
        psi, psi_in = np.empty(self.n_cells), 0.0
        for i in range(self.n_cells - 1, -1, -1):
            psi[i], (psi_in,) = sn.cell_flux(self.medium.sigma_t * self.dr,
                                             source[i] * self.dr, [(1.0, 1.0, psi_in)])
        return psi

    def _direction(self, m, source, psi_low, incoming):
        """Sweeps one ordinate and returns its flux contribution and its upper half-angle flux."""
        mu, weight = self.mu[m], self.weights[m]
        inward = mu < 0.0
        cells = range(self.n_cells - 1, -1, -1) if inward else range(self.n_cells)

        # Vacuum at r = R for the inward half, the reflected flux at r = 0 for the outward one.
        psi_in = 0.0 if inward else incoming[m]
        psi, psi_high = np.empty(self.n_cells), np.empty(self.n_cells)
        for i in cells:
            inner, outer = self.areas[i], self.areas[i + 1]
            a_out, a_in = (inner, outer) if inward else (outer, inner)
            angular = (outer - inner) / weight
            psi[i], (psi_in, psi_high[i]) = sn.cell_flux(
                self.medium.sigma_t * self.volumes[i], source[i] * self.volumes[i],
                [(abs(mu) * a_out, abs(mu) * a_in, psi_in),
                 (angular * self.alpha[m + 1], angular * self.alpha[m], psi_low[i])])

        if inward:
            incoming[len(self.mu) - 1 - m] = psi_in   # reflective boundary at r = 0
        return weight * psi, psi_high

    def sweep(self, source):
        """Scalar flux produced by one isotropic source S; the ordinates carry S/2."""
        q = 0.5 * source
        psi_low = self._starting_direction(q)
        incoming = np.zeros(len(self.mu))

        phi = np.zeros(self.n_cells)
        # Ascending mu: the angular sweep runs from mu = -1 upwards, and every inward
        # ordinate is swept before the outward one it feeds through the r = 0 reflection.
        for m in range(len(self.mu)):
            contribution, psi_low = self._direction(m, q, psi_low, incoming)
            phi += contribution
        return phi

def k_eigenvalue(radius, medium, n_ordinates, n_cells=N_CELLS):
    """KResult of a bare sphere of the given radius."""
    return sn.k_eigenvalue(SphereSolver(radius, medium, n_ordinates, n_cells))

def critical_radius(medium, n_ordinates, guess, n_cells=N_CELLS):
    """Radius at which k = 1."""
    return sn.critical_size(
        lambda R: k_eigenvalue(R, medium, n_ordinates, n_cells).k, guess)
