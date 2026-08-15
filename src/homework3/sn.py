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

def cell_flux(removal, source, links):
    """
    Diamond-difference cell-centre flux, with the set-to-zero negative-flux fixup.

    Each link is an (out coefficient, in coefficient, incoming flux) triple of one
    outgoing face -- the spatial one, plus the angular one in the sphere -- and the
    balance solved is sum(c_out psi_out - c_in psi_in) + removal psi = source.
    Returns the cell-centre flux and the outgoing fluxes, in the order of `links`.
    """
    clamped = [False] * len(links)

    # One pass per link at most: clamping every outgoing flux to zero terminates it.
    for _ in range(len(links) + 1):
        lhs, rhs = removal, source
        for (c_out, c_in, psi_in), off in zip(links, clamped):
            rhs += c_in * psi_in
            if not off:
                lhs += 2.0 * c_out
                rhs += c_out * psi_in

        psi = rhs / lhs
        outgoing = [0.0 if off else 2.0 * psi - psi_in
                    for (_, _, psi_in), off in zip(links, clamped)]
        if all(value >= 0.0 for value in outgoing):
            return psi, outgoing
        clamped = [off or value < 0.0 for off, value in zip(clamped, outgoing)]

    return psi, [max(value, 0.0) for value in outgoing]

def inner_iteration(solver, medium, fission, phi, tol=1e-8, max_iter=2000):
    """Scattering iteration: sweeps until the scalar flux stops moving."""
    for _ in range(max_iter):
        phi_new = solver.sweep(medium.sigma_s * phi + fission)
        settled = np.max(np.abs(phi_new - phi)) <= tol * np.max(phi_new)
        phi = phi_new
        # Without scattering the sweep already inverts the transport operator exactly.
        if settled or medium.sigma_s == 0.0:
            return phi

    raise RuntimeError("The scattering iteration did not converge.")

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
        # The sweep takes a source density, the eigenvalue a source integral.
        phi_new = inner_iteration(solver, medium, medium.nu_sigma_f * phi / k, phi)
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
