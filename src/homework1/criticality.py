"""Question 3: critical slab half-thickness and sphere radius of a bare multiplying medium."""

import numpy as np
from scipy.optimize import brentq
from homework1.exact_solution import compute_nu0_magnitude

HOPF_CONSTANT = 0.710446           # extrapolation distance of the c = 1 Milne problem
Z0_FIT_CORRECTION = -0.0199        # quadratic coefficient q as printed; report §3
MARSHAK_EXTRAPOLATION = 2.0 / 3.0
MARK_EXTRAPOLATION = 1.0 / np.sqrt(3.0)

# Case, de Hoffmann & Placzek (1953), Table 23: the product c z0(c), tabulated
# because it is nearly constant. Reference value for the z0 fit.
CASE_TABLE_23_C = np.arange(0.0, 3.05, 0.1)
CASE_TABLE_23_CZ0 = np.array([
    1.0000, 0.8539, 0.7851, 0.7491, 0.7305, 0.7207, 0.7155, 0.7127, 0.7113,
    0.7106, 0.7104, 0.7106, 0.7109, 0.7113, 0.7118, 0.7123, 0.7129, 0.7134,
    0.7140, 0.7145, 0.7151, 0.7156, 0.7162, 0.7167, 0.7172, 0.7177, 0.7182,
    0.7186, 0.7191, 0.7195, 0.7199,
])

# Same reference, Table 8 Part II: the root k0 = 1/|nu0| of c arctan(k0) = k0.
CASE_TABLE_8_C = np.arange(1.0, 2.05, 0.1)
CASE_TABLE_8_K0 = np.array([
    0.00000, 0.56926, 0.83454, 1.05708, 1.25981, 1.45110, 1.63500, 1.81378,
    1.98883, 2.16107, 2.33112,
])

def extrapolation_distance(c, correction=Z0_FIT_CORRECTION):
    """z0(c) = 0.710446 [1 + q (1-c)^2] / c, the two-term expansion about c = 1."""
    c = np.asarray(c, dtype=float)
    return HOPF_CONSTANT * (1.0 + correction * (1.0 - c)**2) / c

def extrapolation_distance_table(c):
    """z0(c) interpolated from Case's Table 23, interpolating the flat product c z0."""
    c = np.asarray(c, dtype=float)
    if np.any(c < CASE_TABLE_23_C[0]) or np.any(c > CASE_TABLE_23_C[-1]):
        raise ValueError(f"Table 23 covers c in [{CASE_TABLE_23_C[0]}, {CASE_TABLE_23_C[-1]}].")
    return np.interp(c, CASE_TABLE_23_C, CASE_TABLE_23_CZ0) / c

def diffusion_relaxation_length(c):
    """The diffusion counterpart of |nu0(c)|, 1/B = 1/sqrt(3(c-1)) in mean free paths."""
    c = np.asarray(c, dtype=float)
    return 1.0 / np.sqrt(3.0 * (c - 1.0))

def _relaxation_length(c, source):
    """|nu0(c)| from the fit or the transcendental root, or the diffusion 1/B."""
    if source == 'diffusion':
        return diffusion_relaxation_length(c)
    method = {'fit': 'approx', 'exact': 'numerical'}[source]
    return np.array([compute_nu0_magnitude(float(ci), method=method)
                     for ci in np.atleast_1d(c)])

def _z0(c, source):
    """z0 from the printed fit, the sign-flipped fit, Case's table, or a constant."""
    if source == 'fit':
        return extrapolation_distance(c)
    if source == 'fit+':
        return extrapolation_distance(c, correction=-Z0_FIT_CORRECTION)
    if source == 'table':
        return extrapolation_distance_table(c)
    constant = {'marshak': MARSHAK_EXTRAPOLATION, 'mark': MARK_EXTRAPOLATION}[source]
    return np.full_like(np.asarray(c, dtype=float), constant)

# Each method pairs a relaxation length with an extrapolation distance. 'transport'
# is part 3(a), 'marshak' 3(b) and 'mark' 3(c); the rest re-evaluate 3(a).
METHODS = {
    'transport': ('fit', 'fit'),
    'transport-q+': ('fit', 'fit+'),
    'transport-ref': ('exact', 'table'),
    'marshak': ('diffusion', 'marshak'),
    'mark': ('diffusion', 'mark'),
}

METHOD_LABELS = {
    'transport': 'Exact transport',
    'transport-q+': r'Exact transport ($q = +0.0199$)',
    'transport-ref': 'Exact transport (reference)',
    'marshak': 'Diffusion, Marshak',
    'mark': 'Diffusion, Mark',
}

def critical_dimensions(c, method='transport'):
    """
    (nu, z0, a/2, Sigma_t R_c) in mean free paths, from a/2 = (pi/2) nu - z0 and
    Sigma_t R_c = pi nu - z0. Every method shares those relations; see METHODS.
    """
    nu_source, z0_source = METHODS[method]

    c_arr = np.atleast_1d(np.asarray(c, dtype=float))
    if np.any(c_arr <= 1.0):
        raise ValueError("A bare system is critical only for c > 1.")

    nu, z0 = _relaxation_length(c_arr, nu_source), _z0(c_arr, z0_source)
    result = (nu, z0, 0.5 * np.pi * nu - z0, np.pi * nu - z0)

    return tuple(float(v[0]) for v in result) if np.ndim(c) == 0 else result

def critical_dimensions_applied_bc(c, extrapolation):
    """
    (a/2, Sigma_t R_c) with phi + l0 phi' = 0 imposed on the flux shape itself, giving
    B a/2 = arctan(1/(B l0)) and u cot u = 1 - u/(B l0), instead of an extrapolated zero.
    """
    c_arr = np.atleast_1d(np.asarray(c, dtype=float))
    if np.any(c_arr <= 1.0):
        raise ValueError("A bare system is critical only for c > 1.")

    B = 1.0 / diffusion_relaxation_length(c_arr)
    half = np.arctan(1.0 / (B * extrapolation)) / B

    def sphere_root(Bi):
        # g rises from 0+ at u -> 0 and falls to -infinity at u -> pi, so (0, pi) brackets.
        def g(u):
            return u / np.tan(u) - 1.0 + u / (Bi * extrapolation)
        return brentq(g, 1e-6, np.pi - 1e-12, xtol=1e-14, rtol=8.9e-16)

    radius = np.array([sphere_root(Bi) for Bi in B]) / B

    return (float(half[0]), float(radius[0])) if np.ndim(c) == 0 else (half, radius)
