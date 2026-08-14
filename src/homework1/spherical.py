"""Question 4: the k-eigenvalue of a bare sphere, by finite volume and source iteration."""

import numpy as np
from dataclasses import dataclass
from typing import NamedTuple
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from homework1.criticality import extrapolation_distance
from homework1.exact_solution import compute_nu0_magnitude

@dataclass(frozen=True)
class SphericalMedium:
    """Diffusion constants of a bare multiplying medium, in cm^-1 and cm."""
    D: float
    sigma_a: float
    nu_sigma_f: float
    # Surface length: Marshak's 2D for the classical branch, Milne's z0(c) for the
    # asymptotic one. Used as an extrapolation distance by 'extrapolated' and as a
    # linear extrapolation length by 'robin'; the two coincide only as B z0 -> 0,
    # see explanations/07.
    z0: float
    c: float

class KResult(NamedTuple):
    """One k calculation: the eigenvalue, the flux on its mesh, and its cost."""
    k: float
    r: np.ndarray
    phi: np.ndarray
    sweeps: int

def build_medium(sigma_t, sigma_a, nu_sigma_f, approximation='classical'):
    """Diffusion constants of one approximation from three cross sections."""
    c = (sigma_t - sigma_a + nu_sigma_f) / sigma_t

    if approximation == 'classical':
        D = 1.0 / (3.0 * sigma_t)
        return SphericalMedium(D, sigma_a, nu_sigma_f, z0=2.0 * D, c=c)
    if approximation == 'asymptotic':
        if c <= 1.0:
            raise ValueError("The asymptotic approximation needs c > 1 here.")
        # D = (c-1)|nu0|^2 is the continuation above c = 1 of Question 2's
        # D = (1-c) nu0^2; both factors flip sign. See explanations/02.
        D = (c - 1.0) * compute_nu0_magnitude(c, method='numerical')**2 / sigma_t
        z0 = float(extrapolation_distance(c)) / sigma_t
        return SphericalMedium(D, sigma_a, nu_sigma_f, z0=z0, c=c)
    raise ValueError("approximation must be 'classical' or 'asymptotic'")

def buckling(medium):
    """Critical buckling B = sqrt((nu Sigma_f - Sigma_a) / D)."""
    return np.sqrt((medium.nu_sigma_f - medium.sigma_a) / medium.D)

def analytic_critical_radius(medium):
    """Critical radius from the extrapolated-zero condition, R_c = pi/B - z0."""
    return np.pi / buckling(medium) - medium.z0

def _mesh(r_outer, n_cells):
    """Cell-centred spherical mesh: centres, face areas, cell volumes, spacing."""
    h = r_outer / n_cells
    faces = np.linspace(0.0, r_outer, n_cells + 1)
    centres = 0.5 * (faces[:-1] + faces[1:])
    volumes = (4.0 * np.pi / 3.0) * (faces[1:]**3 - faces[:-1]**3)
    return centres, 4.0 * np.pi * faces**2, volumes, h

def _outer_boundary(medium, R, boundary):
    """(outer mesh radius, extrapolation length) for the two treatments: an extrapolated
    zero at R + z0, or a Marshak-type condition applied at R. See explanations/07."""
    if boundary == 'extrapolated':
        return R + medium.z0, 0.0
    if boundary == 'robin':
        return R, medium.z0
    raise ValueError("boundary must be 'extrapolated' or 'robin'")

def _banded_operator(medium, areas, volumes, h, l0):
    """Removal operator -D grad^2 + Sigma_a in scipy's banded layout, row = cell balance."""
    conductance = medium.D * areas / h
    # The innermost face has zero area, which is the symmetry condition at r = 0.
    # The outer face carries the boundary condition; see explanations/06.
    conductance[-1] = medium.D * areas[-1] / (l0 + 0.5 * h)

    ab = np.zeros((3, len(volumes)))
    ab[0, 1:] = -conductance[1:-1]
    ab[1, :] = conductance[:-1] + conductance[1:] + medium.sigma_a * volumes
    ab[2, :-1] = -conductance[1:-1]
    return ab

def k_eigenvalue(R, medium, n_cells=400, boundary='extrapolated', tol=1e-10, max_iter=20000):
    """KResult of a sphere of radius R, by Bell & Glasstone source iteration;
    see explanations/08."""
    r_outer, l0 = _outer_boundary(medium, R, boundary)
    centres, areas, volumes, h = _mesh(r_outer, n_cells)
    ab = _banded_operator(medium, areas, volumes, h, l0)

    phi = np.ones(n_cells)
    fission = medium.nu_sigma_f * phi * volumes
    k = 1.0

    for sweep in range(1, max_iter + 1):
        phi_new = solve_banded((1, 1), ab, fission / k)
        fission_new = medium.nu_sigma_f * phi_new * volumes

        k_new = k * fission_new.sum() / fission.sum()
        converged = abs(k_new - k) < tol * abs(k_new)

        # Renormalise, or the amplitude drifts by a factor k per sweep and the
        # convergence test loses its significant digits.
        scale = 1.0 / phi_new.max()
        phi, fission, k = phi_new * scale, fission_new * scale, k_new

        if converged:
            return KResult(k, centres, phi, sweep)

    raise RuntimeError(f"k iteration did not converge in {max_iter} sweeps.")

def dominance_ratio(medium, R, boundary='extrapolated'):
    """Convergence rate of the source iteration, k_1/k_0, from the first two
    spherical modes; at a critical radius it reduces to c/(4c-3)."""
    r_outer, _ = _outer_boundary(medium, R, boundary)

    def k_mode(n):
        return medium.nu_sigma_f / (medium.sigma_a
                                    + medium.D * (n * np.pi / r_outer)**2)

    return k_mode(2) / k_mode(1)

def neutron_balance(R, medium, n_cells=400, boundary='extrapolated'):
    """(production, absorption, leakage, relative residual) of the converged flux,
    which must balance at a critical radius; see explanations/09."""
    result = k_eigenvalue(R, medium, n_cells, boundary)
    r_outer, l0 = _outer_boundary(medium, R, boundary)
    _, areas, volumes, h = _mesh(r_outer, n_cells)

    production = (medium.nu_sigma_f * result.phi * volumes).sum() / result.k
    absorption = (medium.sigma_a * result.phi * volumes).sum()
    leakage = areas[-1] * medium.D * result.phi[-1] / (l0 + 0.5 * h)
    return production, absorption, leakage, abs(production - absorption - leakage) / production

def _bracket_radius(f, guess, growth=1.05, max_steps=40):
    """Widens an interval about `guess` until f changes sign; starts narrow, since k
    away from criticality costs thousands of sweeps."""
    lo, hi = guess / growth, guess * growth
    for _ in range(max_steps):
        if f(lo) * f(hi) < 0.0:
            return lo, hi
        lo, hi = lo / growth, hi * growth
    raise RuntimeError("Could not bracket a radius with k = 1.")

def critical_radius(medium, n_cells=400, boundary='extrapolated', xtol=1e-12):
    """Radius at which k = 1, by brentq on k(R) - 1, bracketed around the analytic one."""
    def residual(R):
        return k_eigenvalue(R, medium, n_cells, boundary).k - 1.0

    lo, hi = _bracket_radius(residual, analytic_critical_radius(medium))
    return brentq(residual, lo, hi, xtol=xtol)

def mesh_convergence(medium, cell_counts=(25, 50, 100, 200, 400, 800)):
    """(cells, radii, relative errors) against the analytic radius, which is the
    extrapolated-boundary reference; the scheme is second order."""
    exact = analytic_critical_radius(medium)
    radii = np.array([critical_radius(medium, n_cells=n) for n in cell_counts])
    return np.array(cell_counts), radii, np.abs(radii - exact) / exact
