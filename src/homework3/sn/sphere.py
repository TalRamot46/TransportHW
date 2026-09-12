"""Questions 4 and 5: S_N in a sphere, reflective at r = 0 and vacuum at r = R."""

import numpy as np
from typing import NamedTuple
from homework3.sn import core
from homework3.sn.slab import cell_flux as plane_cell_flux

# 100 cells is already mesh converged: at c = 1.5 the S10 radius moves by 7e-6 mfp
# between 50 and 800 cells. See explanations/06.
N_CELLS = 100

def angular_coefficients(mu, weights):
    """Carlson's alpha_{m+1/2} = alpha_{m-1/2} - w_m mu_m, zero at both ends; the
    recursion is what makes a flat flux an exact solution. See explanations/02."""
    alpha = np.concatenate(([0.0], -np.cumsum(weights * mu)))
    alpha[-1] = 0.0   # exactly zero by sum(w mu) = 0; set so no current leaks at mu = +1
    return alpha

class Face(NamedTuple):
    """One outgoing face of a spherical cell: the weight its outgoing flux carries in
    the balance, the weight its inflow carries, and that inflow. The two weights differ
    because a shell has two areas and an angular bin two alphas."""
    a_out: float
    a_in: float
    psi_in: float

def _balance(removal, source, radial_face, angular_face, clamp_radial, clamp_angular):
    """Report eq. (30) solved for psi_{i,m}. A clamped face keeps its inflow but loses
    its outgoing term from both sides."""
    denominator = removal
    numerator = (source + radial_face.a_in * radial_face.psi_in
                 + angular_face.a_in * angular_face.psi_in)

    if not clamp_radial:
        denominator += 2.0 * radial_face.a_out
        numerator += radial_face.a_out * radial_face.psi_in
    if not clamp_angular:
        denominator += 2.0 * angular_face.a_out
        numerator += angular_face.a_out * angular_face.psi_in
    return numerator / denominator

def cell_flux(removal, source, radial_face, angular_face):
    """Cell-centre flux and both outgoing fluxes of one spherical cell, with the
    negative-flux fixup. Report eq. (30), closed by (27) and (28)."""
    clamp_radial = clamp_angular = False

    # Clamping one face lowers psi and can drive the other negative, so the test is
    # repeated; with two faces it settles in at most two extra passes.
    for _ in range(3):
        psi = _balance(removal, source, radial_face, angular_face,
                       clamp_radial, clamp_angular)
        psi_radial = 0.0 if clamp_radial else 2.0 * psi - radial_face.psi_in
        psi_angular = 0.0 if clamp_angular else 2.0 * psi - angular_face.psi_in
        if psi_radial >= 0.0 and psi_angular >= 0.0:
            return psi, psi_radial, psi_angular
        clamp_radial = clamp_radial or psi_radial < 0.0
        clamp_angular = clamp_angular or psi_angular < 0.0
    return psi, max(psi_radial, 0.0), max(psi_angular, 0.0)

class SphereSolver:
    """One transport sweep over [0, R], with the angular redistribution term of
    spherical geometry differenced in angle as well as in space."""

    def __init__(self, radius, medium, n_ordinates, n_cells=N_CELLS):
        self.medium = medium
        self.mu, self.weights = core.ordinates(n_ordinates)
        self.alpha = angular_coefficients(self.mu, self.weights)
        self.n_cells = n_cells
        self.dr = radius / n_cells

        faces = np.linspace(0.0, radius, n_cells + 1)
        self.centres = 0.5 * (faces[:-1] + faces[1:])
        self.areas = 4.0 * np.pi * faces**2
        self.volumes = (4.0 * np.pi / 3.0) * (faces[1:]**3 - faces[:-1]**3)

    def _angular_face(self, i, m, psi_in):
        """The mu-face of cell (i, m): the two alphas of the bin, area-weighted."""
        beta = (self.areas[i + 1] - self.areas[i]) / self.weights[m]
        return Face(beta * self.alpha[m + 1], beta * self.alpha[m], psi_in)

    def _cell(self, i, m, q, radial_face, psi_low_i):
        """One cell of a march: its radial face is the direction's, its angular face
        always runs from m-1/2 to m+1/2."""
        return cell_flux(self.medium.sigma_t * self.volumes[i], q[i] * self.volumes[i],
                         radial_face, self._angular_face(i, m, psi_low_i))

    def _starting_direction(self, q):
        """Half-angle flux at mu = -1, where 1 - mu^2 = 0 removes the angular term and
        the balance is the plane one at |mu| = 1; see explanations/02."""
        psi, psi_in = np.empty(self.n_cells), 0.0
        for i in range(self.n_cells - 1, -1, -1):
            psi[i], psi_in = plane_cell_flux(self.medium.sigma_t * self.dr,
                                             q[i] * self.dr, 1.0, psi_in)
        return psi

    def sweep_backward(self, m, q, psi_low):
        """mu_m < 0: marches inwards from the vacuum face at r = R to r = 0. The radial
        outflow is through the inner face, the inflow through the outer one."""
        mu_abs = abs(self.mu[m])
        psi, psi_high, psi_in = np.empty(self.n_cells), np.empty(self.n_cells), 0.0
        for i in range(self.n_cells - 1, -1, -1):
            radial_face = Face(mu_abs * self.areas[i], mu_abs * self.areas[i + 1], psi_in)
            psi[i], psi_in, psi_high[i] = self._cell(i, m, q, radial_face, psi_low[i])
        return psi, psi_high, psi_in

    def sweep_forward(self, m, q, psi_low, psi_in):
        """mu_m > 0: marches outwards from r = 0 to the vacuum face, starting from the
        flux its mirror ordinate -mu_m left at the centre. The radial outflow is through
        the outer face, the inflow through the inner one."""
        mu_abs = self.mu[m]
        psi, psi_high = np.empty(self.n_cells), np.empty(self.n_cells)
        for i in range(self.n_cells):
            radial_face = Face(mu_abs * self.areas[i + 1], mu_abs * self.areas[i], psi_in)
            psi[i], psi_in, psi_high[i] = self._cell(i, m, q, radial_face, psi_low[i])
        return psi, psi_high

    def sweep_all_angles(self, source):
        """One S_N iteration: every backward ordinate, then every forward one. That is
        ascending mu, which the angular recursion needs -- it runs from mu = -1 upwards,
        and it also puts every backward ordinate before the forward one it feeds through
        the reflection at r = 0."""
        q = 0.5 * source
        psi_low = self._starting_direction(q)
        phi, reflected = np.zeros(self.n_cells), np.zeros(len(self.mu))
        mirror = len(self.mu) - 1

        for m in np.flatnonzero(self.mu < 0.0):
            psi, psi_low, reflected[mirror - m] = self.sweep_backward(m, q, psi_low)
            phi += self.weights[m] * psi
        for m in np.flatnonzero(self.mu > 0.0):
            psi, psi_low = self.sweep_forward(m, q, psi_low, reflected[m])
            phi += self.weights[m] * psi
        return phi

def k_eigenvalue(radius, medium, n_ordinates, n_cells=N_CELLS):
    """KResult of a bare sphere of the given radius."""
    return core.k_eigenvalue(SphereSolver(radius, medium, n_ordinates, n_cells))

def critical_radius(medium, n_ordinates, guess, n_cells=N_CELLS):
    """Radius at which k = 1."""
    return core.critical_size(
        lambda R: k_eigenvalue(R, medium, n_ordinates, n_cells).k, guess)
