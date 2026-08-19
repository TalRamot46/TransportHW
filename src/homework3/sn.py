"""Geometry-independent S_N machinery: the quadrature, the diamond-difference cell
solve with its negative-flux fixup, and the Bell & Glasstone k iteration."""

import numpy as np
from dataclasses import dataclass
from typing import NamedTuple
from scipy.optimize import brentq

@dataclass(frozen=True)
class Medium:
    """Total, scattering and production cross sections of one homogeneous medium, in cm^-1."""
    sigma_t: float
    sigma_s: float
    nu_sigma_f: float

    @property
    def c(self):
        """Secondaries per collision, (Sigma_s + nu Sigma_f) / Sigma_t."""
        return (self.sigma_s + self.nu_sigma_f) / self.sigma_t

def multiplying_medium(c, sigma_t=1.0):
    """Medium of c secondaries per collision, all counted as fission; the critical
    size depends on c alone, so the split is free. See explanations/02."""
    return Medium(sigma_t, 0.0, c * sigma_t)

def ordinates(n):
    """Gauss-Legendre ordinates and weights on [-1, 1], ascending, summing to 2."""
    return np.polynomial.legendre.leggauss(n)

class Face(NamedTuple):
    """
    One outgoing face of a cell, as it appears in the cell balance.

    Its contribution is `a_out psi_out - a_in psi_in`. The two coefficients differ only
    where the face is weighted differently on the way in and on the way out -- the two
    areas of a spherical shell, the two alphas of an angular bin -- so in the slab they
    are both |mu|. See report eq. (31).
    """
    a_out: float
    a_in: float
    psi_in: float

def cell_flux(removal, source, faces):
    """
    Diamond-difference cell-centre flux, with the set-to-zero negative-flux fixup.

    Solves the cell balance over any number of outgoing faces,

        sum_f (a_out psi_out - a_in psi_in) + removal psi = source

    by closing each face with the diamond relation psi_out = 2 psi - psi_in. Collecting
    psi gives report eq. (32), which is what the loop below evaluates:

        psi = [ source + sum_f (a_out + a_in) psi_in ] / [ removal + 2 sum_f a_out ]

    A face whose outgoing flux has been clamped to zero contributes no a_out to either
    sum, leaving only its -a_in psi_in inflow. Returns the cell-centre flux and the
    outgoing fluxes, in the order of `faces`.
    """
    clamped = [False] * len(faces)

    # One pass per face at most: clamping every outgoing flux to zero terminates it.
    for _ in range(len(faces) + 1):
        denominator, numerator = removal, source
        for face, off in zip(faces, clamped):
            if off:
                numerator += face.a_in * face.psi_in
            else:
                denominator += 2.0 * face.a_out
                numerator += (face.a_out + face.a_in) * face.psi_in

        psi = numerator / denominator
        outgoing = [0.0 if off else 2.0 * psi - face.psi_in
                    for face, off in zip(faces, clamped)]
        if all(value >= 0.0 for value in outgoing):
            return psi, outgoing
        clamped = [off or value < 0.0 for off, value in zip(clamped, outgoing)]

    return psi, [max(value, 0.0) for value in outgoing]

def run_sn(solver, medium, fission, phi, tol=1e-8, max_iter=2000):
    """The S_N solve itself: repeats the angular iteration at a fixed fission source
    until the scalar flux stops moving. `k_eigenvalue` only wraps this."""
    for _ in range(max_iter):
        phi_new = solver.sn_iteration(medium.sigma_s * phi + fission)
        settled = np.max(np.abs(phi_new - phi)) <= tol * np.max(phi_new)
        phi = phi_new
        # Without scattering one iteration already inverts the transport operator exactly.
        if settled or medium.sigma_s == 0.0:
            return phi
    raise RuntimeError("The S_N iteration did not converge.")

class KResult(NamedTuple):
    """One k calculation: the eigenvalue, the flux on its mesh, and its cost."""
    k: float
    x: np.ndarray
    phi: np.ndarray
    outers: int

def k_eigenvalue(solver, tol=1e-9, max_iter=2000):
    """KResult of one geometry, by the Bell & Glasstone outer iteration
    k <- k P_new / P_old; see explanations/03."""
    medium = solver.medium
    phi = np.ones(solver.n_cells)
    production = (medium.nu_sigma_f * phi * solver.volumes).sum()
    k = 1.0

    for outer in range(1, max_iter + 1):
        # The S_N solve takes a source density, the eigenvalue a source integral.
        phi_new = run_sn(solver, medium, medium.nu_sigma_f * phi / k, phi)
        production_new = (medium.nu_sigma_f * phi_new * solver.volumes).sum()

        k_new = k * production_new / production
        converged = abs(k_new - k) < tol * abs(k_new)

        # Renormalise, or the amplitude drifts by a factor k per outer and the
        # convergence test loses its significant digits.
        scale = 1.0 / phi_new.max()
        phi, production, k = phi_new * scale, production_new * scale, k_new

        if converged:
            return KResult(k, solver.centres, phi, outer)

    raise RuntimeError(f"The k iteration did not converge in {max_iter} outers.")

def _bracket(f, guess, growth=1.03, max_steps=40):
    """Widens an interval about `guess` until f changes sign; starts narrow, since
    each evaluation of f is a whole power iteration."""
    lo, hi = guess / growth, guess * growth
    for _ in range(max_steps):
        if f(lo) * f(hi) < 0.0:
            return lo, hi
        lo, hi = lo / growth, hi * growth
    raise RuntimeError("Could not bracket a size with k = 1.")

def critical_size(k_of_size, guess, rtol=1e-7):
    """Size at which k = 1, by brentq on k(size) - 1 bracketed around `guess`. The
    tolerance is relative, since the sizes are mean free paths in one question and
    centimetres in the next."""
    def residual(size):
        return k_of_size(size) - 1.0

    lo, hi = _bracket(residual, guess)
    return brentq(residual, lo, hi, xtol=rtol * guess)
