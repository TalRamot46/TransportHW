import numpy as np
from scipy.optimize import brentq
from homework1.exact_solution import compute_nu0_magnitude

# ---------------------------------------------------------------------------
# Question 3: critical dimensions from the asymptotic (exact-transport) relations
#
# For c > 1 the discrete eigenvalue is imaginary, nu0 = i |nu0|, so the
# asymptotic scalar flux is trigonometric rather than exponential:
#
#     planar     phi_as(x) = A cos(Sigma_t x / |nu0|)      (symmetry kills the sine)
#     spherical  phi_as(r) = A sin(Sigma_t r / |nu0|) / (Sigma_t r)
#
# A bare system is critical when the asymptotic flux extrapolates to zero at the
# extrapolation distance z0(c) outside the physical surface -- phi_as vanishes at
# -z0 measured from the boundary. Requiring the first zero of each shape to fall
# there gives the two criticality relations of this question, in mean free paths
# (Sigma_t = 1):
#
#     a / 2        = (pi / 2) |nu0(c)| - z0(c)          (planar half-thickness)
#     Sigma_t R_c  =  pi      |nu0(c)| - z0(c)          (critical sphere radius)
#
# The planar relation places the quarter-wave zero of the cosine at the
# extrapolated surface; the spherical one places the first zero of the sine
# there, at argument pi. Both are "exact transport" in the sense that they use
# the exact transport eigenvalue and the exact Milne extrapolation distance -- no
# diffusion coefficient enters -- but they still neglect the transient part of
# the flux near the boundary, so they are asymptotic rather than exact.
# ---------------------------------------------------------------------------

# Hopf's constant: the extrapolation distance of the c = 1 Milne problem.
HOPF_CONSTANT = 0.710446

# Coefficient of the quadratic term of the z0 fit, as printed in the course
# notes. See extrapolation_distance for why the sign is a keyword.
Z0_FIT_CORRECTION = -0.0199

# Case, de Hoffmann & Placzek (1953), Table 23: the product c z0(c), which is
# the natural tabulation because z0 ~ 0.7104 / c. Used as the reference value of
# the extrapolation distance, against which the fit below is measured.
CASE_TABLE_23_C = np.arange(0.0, 3.05, 0.1)
CASE_TABLE_23_CZ0 = np.array([
    1.0000, 0.8539, 0.7851, 0.7491, 0.7305, 0.7207, 0.7155, 0.7127, 0.7113,
    0.7106, 0.7104, 0.7106, 0.7109, 0.7113, 0.7118, 0.7123, 0.7129, 0.7134,
    0.7140, 0.7145, 0.7151, 0.7156, 0.7162, 0.7167, 0.7172, 0.7177, 0.7182,
    0.7186, 0.7191, 0.7195, 0.7199,
])

# Case, de Hoffmann & Placzek (1953), Table 8 Part II: the root k0 = 1 / |nu0|
# of c arctan(k0) = k0 for c >= 1. Used to verify the eigenvalue solver and its
# analytic fit on the multiplying branch.
CASE_TABLE_8_C = np.arange(1.0, 2.05, 0.1)
CASE_TABLE_8_K0 = np.array([
    0.00000, 0.56926, 0.83454, 1.05708, 1.25981, 1.45110, 1.63500, 1.81378,
    1.98883, 2.16107, 2.33112,
])

def extrapolation_distance(c, correction=Z0_FIT_CORRECTION):
    """
    Extrapolation distance z0(c) in mean free paths, from the expansion about
    c = 1

        z0(c) = 0.710446 [ 1/c + q (1 - c)^2 / c + O((1 - c)^3 / c) ]

    The leading term is exact in the sense that 0.710446 / c = 0.710446 nu0
    arctanh(1 / nu0), so all of the c-dependence at first order comes from the
    dispersion relation itself.

    Parameters:
    correction : float, optional
        The coefficient q. The course notes print q = -0.0199, which is the
        default; Case's own Table 23 is reproduced by q = +0.0199. See
        `extrapolation_distance_table` and the discussion in the report: the
        tabulated c z0 rises on *both* sides of c = 1, so a negative q has the
        wrong sign, while |q| = 0.0199 matches the table to its four printed
        digits at c = 0.9 and c = 1.1. The choice moves a/2 by less than 0.05 %
        below c = 1.2 but by 3 % at c = 2.
    """
    c = np.asarray(c, dtype=float)
    return HOPF_CONSTANT * (1.0 + correction * (1.0 - c)**2) / c

def extrapolation_distance_table(c):
    """
    Extrapolation distance z0(c) interpolated from Case's Table 23, for use as
    the reference value.

    The interpolation is done on c z0 rather than on z0, because c z0 varies by
    less than 1 % over the whole table while z0 itself varies by a factor of
    three; linear interpolation of the flat quantity is accurate to about 1e-6,
    below the resolution of the four tabulated digits.
    """
    c = np.asarray(c, dtype=float)
    if np.any(c < CASE_TABLE_23_C[0]) or np.any(c > CASE_TABLE_23_C[-1]):
        raise ValueError(f"Table 23 covers c in "
                         f"[{CASE_TABLE_23_C[0]}, {CASE_TABLE_23_C[-1]}].")
    return np.interp(c, CASE_TABLE_23_C, CASE_TABLE_23_CZ0) / c

# ---------------------------------------------------------------------------
# The diffusion approximations, parts 3(b) and 3(c)
#
# In mean free paths (Sigma_t = 1) classical diffusion has D = 1/3, and a
# multiplying medium has a negative absorption cross section, Sigma_a = 1 - c
# < 0. The critical equation
#
#     -D phi'' + (1 - c) phi = 0    =>    phi'' + B^2 phi = 0,   B^2 = 3 (c - 1)
#
# has exactly the flux shapes of the transport problem, cos(Bx) and
# sin(Br)/r, with the transport relaxation length |nu0(c)| replaced by its
# diffusion counterpart 1/B. The criticality relations are therefore unchanged
# in form; only the two inputs move:
#
#     |nu0(c)| -> 1/B = 1/sqrt(3(c-1))
#     z0(c)    -> 2/3        (Marshak)   or   1/sqrt(3)   (Mark)
#
# The two constants are the extrapolation distances of the two standard
# diffusion boundary conditions, and for each of them the extrapolation distance
# and the linear extrapolation length coincide, z0 = l0. Marshak sets the
# incoming partial current to zero, J^-(s) = phi/4 + D phi'/2 = 0, giving
# z0 = 2D = 2/3; Mark sets the incoming angular flux to zero along the P1
# quadrature direction mu = 1/sqrt(3), giving z0 = sqrt(3) D = 1/sqrt(3).
#
# Both are properties of the diffusion operator alone: neither depends on c,
# which is the essential difference from the transport z0(c) above and the
# reason the two approximations diverge as the system shrinks.
# ---------------------------------------------------------------------------

MARSHAK_EXTRAPOLATION = 2.0 / 3.0
MARK_EXTRAPOLATION = 1.0 / np.sqrt(3.0)

def diffusion_relaxation_length(c):
    """
    The diffusion counterpart of |nu0(c)|, 1/B = 1/sqrt(3(c-1)), in mean free
    paths.

    Expanding the transport dispersion relation about c = 1 gives
    |nu0| = 1/sqrt(3(c-1)) + O(1), so the two agree in the near-critical limit
    and separate as c grows: 1/B is 0.8 % high at c = 1.02 and 34 % high at
    c = 2.
    """
    c = np.asarray(c, dtype=float)
    return 1.0 / np.sqrt(3.0 * (c - 1.0))

def _relaxation_length(c, source):
    """
    |nu0(c)| on the multiplying branch, either from the analytic fit ('fit') or
    from the root of the transcendental equation ('exact'), or the diffusion
    counterpart 1/B ('diffusion').
    """
    if source == 'diffusion':
        return diffusion_relaxation_length(c)

    method = {'fit': 'approx', 'exact': 'numerical'}.get(source)
    if method is None:
        raise ValueError("nu_source must be 'fit', 'exact' or 'diffusion'")
    return np.array([compute_nu0_magnitude(float(ci), method=method)
                     for ci in np.atleast_1d(c)])

def _z0(c, source):
    """
    z0(c) from the fit as printed ('fit'), from the fit with the sign of the
    quadratic term flipped ('fit+'), interpolated from Case's table ('table'),
    or one of the two constant diffusion values ('marshak', 'mark').
    """
    if source == 'fit':
        return extrapolation_distance(c)
    if source == 'fit+':
        return extrapolation_distance(c, correction=-Z0_FIT_CORRECTION)
    if source == 'table':
        return extrapolation_distance_table(c)
    if source == 'marshak':
        return np.full_like(np.asarray(c, dtype=float), MARSHAK_EXTRAPOLATION)
    if source == 'mark':
        return np.full_like(np.asarray(c, dtype=float), MARK_EXTRAPOLATION)
    raise ValueError("z0_source must be 'fit', 'fit+', 'table', 'marshak' or 'mark'")

# Each method pairs a source for the relaxation length with a source for the
# extrapolation distance. 'transport' is part 3(a) of the assignment, 'marshak'
# is 3(b) and 'mark' is 3(c); the other two are reference evaluations of 3(a).
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
    Critical planar half-thickness and critical sphere radius for a bare
    multiplying medium, in mean free paths.

    Every method evaluates the same pair of relations,

        a / 2 = (pi / 2) nu - z0,      Sigma_t R_c = pi nu - z0

    and differs only in where the relaxation length nu and the extrapolation
    distance z0 come from; see METHODS.

    Parameters:
    c : float or array_like
        Scattering ratio, which must exceed 1 for a critical system to exist.
    method : str, optional
        A key of METHODS.

    Returns:
    nu : relaxation length, |nu0(c)| or its diffusion counterpart 1/B
    z0 : extrapolation distance
    half_thickness : a/2
    radius : Sigma_t R_c

    Each is a float for scalar `c` and an array otherwise. Values are returned
    even where they are negative, which happens once the medium is so strongly
    multiplying that the extrapolation distance exceeds the quarter wavelength;
    the caller decides what to do with them.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {sorted(METHODS)}")
    nu_source, z0_source = METHODS[method]

    c_arr = np.atleast_1d(np.asarray(c, dtype=float))
    if np.any(c_arr <= 1.0):
        raise ValueError("A bare system is critical only for c > 1.")

    nu = _relaxation_length(c_arr, nu_source)
    z0 = _z0(c_arr, z0_source)

    result = (nu, z0, 0.5 * np.pi * nu - z0, np.pi * nu - z0)

    if np.ndim(c) == 0:
        return tuple(float(value[0]) for value in result)
    return result

def critical_dimensions_applied_bc(c, extrapolation):
    """
    The same diffusion problem as 'marshak' / 'mark', but with the boundary
    condition imposed on the flux shape itself instead of being replaced by an
    extrapolated zero.

    Writing the condition as phi(s) + l0 phi'(s) = 0 at the surface s, with
    l0 = 2/3 for Marshak and 1/sqrt(3) for Mark, and applying it to the
    diffusion shapes cos(Bx) and sin(Br)/r gives

        slab     B a/2 = arctan(1 / (B l0))
        sphere   u cot u = 1 - u / (B l0),        u = B R in (0, pi)

    rather than B a/2 = pi/2 - B l0 and u = pi - B l0. The two agree to first
    order in B l0 -- arctan(1/y) = pi/2 - y + O(y^3) -- so they coincide as
    c -> 1, where the system is many mean free paths across and the boundary
    layer is a small correction. They part company as the system shrinks: at
    c = 2 the slab half-thickness differs by a factor of 1.7.

    Returns (half_thickness, radius) in mean free paths.
    """
    c_arr = np.atleast_1d(np.asarray(c, dtype=float))
    if np.any(c_arr <= 1.0):
        raise ValueError("A bare system is critical only for c > 1.")

    B = 1.0 / diffusion_relaxation_length(c_arr)
    half = np.arctan(1.0 / (B * extrapolation)) / B

    def sphere_root(Bi):
        # g(0+) = 0 from above (u cot u -> 1 - u^2/3) and g -> -infinity as
        # u -> pi, so the single physical root is bracketed by (0, pi).
        def g(u):
            return u / np.tan(u) - 1.0 + u / (Bi * extrapolation)
        return brentq(g, 1e-6, np.pi - 1e-12, xtol=1e-14, rtol=8.9e-16)

    radius = np.array([sphere_root(Bi) for Bi in B]) / B

    if np.ndim(c) == 0:
        return float(half[0]), float(radius[0])
    return half, radius
