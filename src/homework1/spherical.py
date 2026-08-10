import numpy as np
from dataclasses import dataclass
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from homework1.criticality import extrapolation_distance
from homework1.exact_solution import compute_nu0_magnitude

# ---------------------------------------------------------------------------
# Question 4: the k-eigenvalue problem for a bare sphere
#
#     -D (1/r^2) d/dr (r^2 dphi/dr) + Sigma_a phi = (1/k) nu Sigma_f phi
#
# The fission source is divided by k, so a solution exists for every radius;
# criticality is the radius at which k = 1. The two approximations differ only
# in D and in the extrapolation distance -- see build_medium.
#
# Derivations: explanations/01-spherical-finite-volume.md (discretisation),
# 02-bell-glasstone-k-iteration.md (the iteration), 03-boundary-and-analytic.md
# (the boundary conditions and the analytic radius compared against).
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SphericalMedium:
    """Diffusion constants of a bare multiplying medium, in cm^-1 and cm."""
    D: float
    sigma_a: float
    nu_sigma_f: float
    z0: float
    c: float

def build_medium(sigma_t, sigma_a, nu_sigma_f, approximation='classical'):
    """
    Diffusion constants for one of the two approximations, from three cross
    sections. Lengths come out in whatever unit sigma_t is reciprocal to.
    """
    c = (sigma_t - sigma_a + nu_sigma_f) / sigma_t

    if approximation == 'classical':
        D = 1.0 / (3.0 * sigma_t)
        z0 = 2.0 * D                      # Marshak, z0 = l0 = 2D
    elif approximation == 'asymptotic':
        if c <= 1.0:
            raise ValueError("The asymptotic approximation needs c > 1 here.")
        # D is fixed by making the diffusion buckling equal the transport one,
        # B = 1/|nu0| in mean free paths; see homework1.diffusion for the
        # c < 1 form D = (1 - c) nu0^2, of which this is the continuation.
        nu0 = compute_nu0_magnitude(c, method='numerical')
        D = (c - 1.0) * nu0**2 / sigma_t
        z0 = float(extrapolation_distance(c)) / sigma_t
    else:
        raise ValueError("approximation must be 'classical' or 'asymptotic'")

    return SphericalMedium(D=D, sigma_a=sigma_a, nu_sigma_f=nu_sigma_f, z0=z0, c=c)

def buckling(medium):
    """Critical buckling B = sqrt((nu Sigma_f - Sigma_a) / D)."""
    return np.sqrt((medium.nu_sigma_f - medium.sigma_a) / medium.D)

def analytic_critical_radius(medium):
    """Critical radius from the extrapolated-zero condition, R_c = pi/B - z0."""
    return np.pi / buckling(medium) - medium.z0

# ---------------------------------------------------------------------------
# Discretisation
# ---------------------------------------------------------------------------

def _mesh(r_outer, n_cells):
    """Cell-centred spherical mesh: centres, face areas, cell volumes, spacing."""
    h = r_outer / n_cells
    faces = np.linspace(0.0, r_outer, n_cells + 1)
    centres = 0.5 * (faces[:-1] + faces[1:])
    areas = 4.0 * np.pi * faces**2
    volumes = (4.0 * np.pi / 3.0) * (faces[1:]**3 - faces[:-1]**3)
    return centres, areas, volumes, h

def _outer_boundary(medium, R, boundary):
    """
    (outer mesh radius, linear extrapolation length) for the two treatments:
    a zero of the flux at the extrapolated radius R + z0, or the condition
    phi + z0 phi' = 0 applied at the physical surface R.
    """
    if boundary == 'extrapolated':
        return R + medium.z0, 0.0
    if boundary == 'robin':
        return R, medium.z0
    raise ValueError("boundary must be 'extrapolated' or 'robin'")

def _banded_operator(medium, areas, volumes, h, l0):
    """
    The removal operator -D grad^2 + Sigma_a on the mesh, in scipy's banded
    (3, N) layout. Row i is the neutron balance over cell i.
    """
    n = len(volumes)
    D = medium.D

    # Conductance of each interior face, D A / h. Face 0 has zero area, which
    # is exactly the symmetry condition dphi/dr = 0 at the origin -- it needs
    # no special case. The outer face carries the boundary condition instead:
    # eliminating the surface flux gives the leakage D A / (l0 + h/2).
    conductance = D * areas / h
    conductance[-1] = D * areas[-1] / (l0 + 0.5 * h)

    ab = np.zeros((3, n))
    ab[0, 1:] = -conductance[1:-1]                       # upper diagonal
    ab[1, :] = conductance[:-1] + conductance[1:] + medium.sigma_a * volumes
    ab[2, :-1] = -conductance[1:-1]                      # lower diagonal
    return ab

# ---------------------------------------------------------------------------
# Bell & Glasstone source iteration (pp. 189-192)
# ---------------------------------------------------------------------------

def k_eigenvalue(R, medium, n_cells=400, boundary='extrapolated',
                 tol=1e-10, max_iter=20000):
    """
    Multiplication factor of a bare sphere of radius R, by source iteration:
    solve the removal operator against the current fission source, then update
    k by the ratio of successive fission-source integrals.

    Returns (k, cell centres, normalised flux).
    """
    r_outer, l0 = _outer_boundary(medium, R, boundary)
    centres, areas, volumes, h = _mesh(r_outer, n_cells)
    ab = _banded_operator(medium, areas, volumes, h, l0)

    phi = np.ones(n_cells)
    fission = medium.nu_sigma_f * phi * volumes
    k = 1.0

    for _ in range(max_iter):
        phi_new = solve_banded((1, 1), ab, fission / k)
        fission_new = medium.nu_sigma_f * phi_new * volumes

        k_new = k * fission_new.sum() / fission.sum()
        converged = abs(k_new - k) < tol * abs(k_new)

        # Renormalise, otherwise the flux drifts by a factor of k each sweep
        # and the iteration loses precision long before it converges.
        scale = 1.0 / phi_new.max()
        phi, fission, k = phi_new * scale, fission_new * scale, k_new

        if converged:
            break
    else:
        raise RuntimeError(f"k iteration did not converge in {max_iter} sweeps.")

    return k, centres, phi

def _bracket_radius(f, guess, growth=1.05, max_steps=40):
    """
    Widens an interval about `guess` until f changes sign across it.

    The initial bracket is deliberately narrow. Source iteration converges at
    the dominance ratio k_1/k_0, which tends to 1 for a large weakly
    multiplying sphere, so evaluating k far from the critical radius costs
    thousands of sweeps; the numerical radius is within a per cent of the
    analytic guess, so +/- 5 % almost always brackets on the first try.
    """
    lo, hi = guess / growth, guess * growth
    for _ in range(max_steps):
        if f(lo) * f(hi) < 0.0:
            return lo, hi
        lo, hi = lo / growth, hi * growth
    raise RuntimeError("Could not bracket a radius with k = 1.")

def critical_radius(medium, n_cells=400, boundary='extrapolated', xtol=1e-10):
    """
    Radius at which k = 1, found by bisection (brentq) on k(R) - 1, which is
    monotonically increasing in R. Bracketed around the analytic radius.
    """
    def residual(R):
        return k_eigenvalue(R, medium, n_cells, boundary)[0] - 1.0

    lo, hi = _bracket_radius(residual, analytic_critical_radius(medium))
    return brentq(residual, lo, hi, xtol=xtol)

def mesh_convergence(medium, cell_counts=(25, 50, 100, 200, 400, 800)):
    """
    Critical radius against mesh refinement, with the relative error of each
    against the analytic radius. The finite-volume scheme is second order, so
    the error should fall by four per doubling.
    """
    exact = analytic_critical_radius(medium)
    radii = np.array([critical_radius(medium, n_cells=n) for n in cell_counts])
    return np.array(cell_counts), radii, np.abs(radii - exact) / exact
