import numpy as np
from scipy.optimize import brentq
import scipy.integrate as integrate

def compute_nu0_numerical(c):
    """
    Computes the discrete eigenvalue nu0 numerically.
    Solves: 1/c = arctanh(k0)/k0 for k0 = 1/nu0 in (0, 1).
    """
    if c <= 0.0:
        return None
    if c >= 1.0:
        raise ValueError("Scattering ratio c must be less than 1.0 for subcritical infinite medium.")

    # Objective function: f(k0) = arctanh(k0)/k0 - 1/c
    def f(k0):
        if k0 == 0.0:
            return 1.0 - 1.0 / c
        return np.arctanh(k0) / k0 - 1.0 / c

    # Bracket the root. Since c < 1, 1/c > 1.
    # At k0 = 0, f(k0) = 1 - 1/c < 0.
    # As k0 -> 1, f(k0) -> +inf.
    # We use 1.0 - 1e-15 as the upper bound for brentq.
    try:
        k0 = brentq(f, 0.0, 1.0 - 1e-15)
        return 1.0 / k0
    except ValueError:
        # Fallback with slightly larger safety margin if needed
        k0 = brentq(f, 0.0, 1.0 - 1e-12)
        return 1.0 / k0

def _nu0_fit_exponent(c):
    """
    Exponent p(c) of the analytic fit for the discrete eigenvalue,
    nu0 = 1 / sqrt(1 - c^p(c)).
    """
    return 2.47412 + 0.00363081 / (c**2) - 0.0352458 * c + 0.557498 / c

def compute_nu0_approx(c):
    """
    Computes the discrete eigenvalue nu0 using the high-accuracy analytic approximation:
    nu0 approx 1 / sqrt(1 - c^p(c))
    """
    if c <= 0.0:
        return None
    if c >= 1.0:
        raise ValueError("Scattering ratio c must be less than 1.0.")

    return 1.0 / np.sqrt(1.0 - c**_nu0_fit_exponent(c))

def compute_nu0(c, method='numerical'):
    """
    Computes nu0 using the specified method: 'numerical' or 'approx'.
    """
    if method == 'numerical':
        return compute_nu0_numerical(c)
    elif method == 'approx':
        return compute_nu0_approx(c)
    else:
        raise ValueError("Method must be 'numerical' or 'approx'")

# ---------------------------------------------------------------------------
# The multiplying branch, c > 1
#
# The transcendental equation for the discrete eigenvalue,
#
#     c nu0 arctanh(1 / nu0) = 1
#
# has no real root once c > 1: the eigenvalue moves onto the imaginary axis,
# nu0 = i |nu0|. Writing k0 = 1 / |nu0| and using arctanh(i k) = i arctan(k),
# the equation becomes real again,
#
#     c arctan(k0) = k0
#
# which is the form tabulated by Case, de Hoffmann & Placzek (Table 8, Part II).
# An imaginary nu0 turns the infinite-medium mode exp(-x / nu0) into cos(x/|nu0|),
# i.e. the flux shape of a critical system, which is what Question 3 needs.
# ---------------------------------------------------------------------------

def compute_nu0_magnitude_numerical(c):
    """
    Computes |nu0| for a multiplying medium (c > 1) by solving c arctan(k0) = k0
    for k0 = 1 / |nu0|, and returns 1 / k0.
    """
    if c <= 1.0:
        raise ValueError("The imaginary eigenvalue branch requires c > 1.")

    def f(k0):
        return c * np.arctan(k0) - k0

    # f(0) = 0 with f'(0) = c - 1 > 0, so f rises away from the trivial root, and
    # f(c pi/2) = c (arctan(c pi/2) - pi/2) < 0 since arctan is bounded by pi/2.
    # The interval below therefore brackets the single positive root for any c > 1.
    return 1.0 / brentq(f, 1e-12, c * np.pi / 2.0, xtol=1e-15, rtol=8.9e-16)

def compute_nu0_magnitude_approx(c):
    """
    Computes |nu0| for a multiplying medium (c > 1) from the same analytic fit
    used below c = 1.

    For c > 1 the fitted radicand 1 - c^p(c) turns negative, which is exactly the
    statement that nu0 has become imaginary; taking the magnitude gives
    |nu0| = 1 / sqrt(c^p(c) - 1). Measured against the transcendental root this
    is accurate to 0.0001 % at c = 1.02 and 0.10 % at c = 2.
    """
    if c <= 1.0:
        raise ValueError("The imaginary eigenvalue branch requires c > 1.")

    return 1.0 / np.sqrt(c**_nu0_fit_exponent(c) - 1.0)

def compute_nu0_magnitude(c, method='approx'):
    """
    Computes |nu0| for c > 1 using the specified method: 'numerical' or 'approx'.
    """
    if method == 'numerical':
        return compute_nu0_magnitude_numerical(c)
    elif method == 'approx':
        return compute_nu0_magnitude_approx(c)
    else:
        raise ValueError("Method must be 'numerical' or 'approx'")

def compute_N0_plus(c, nu0):
    """
    Computes the discrete normalization factor N0^+.
    """
    if c <= 0.0 or nu0 is None:
        return None
    return 0.5 * c * (nu0**3) * (c / (nu0**2 - 1.0) - 1.0 / (nu0**2))

def compute_lambda(nu, c):
    """
    Computes the dispersion function lambda(nu) for nu in [0, 1].
    """
    if nu == 0.0:
        return 1.0
    if nu >= 1.0:
        return -np.inf
    return 1.0 - c * nu * np.arctanh(nu)

def compute_N_nu(nu, c):
    """
    Computes the continuous normalization factor N_nu.
    """
    if nu == 0.0:
        return 0.0
    lam = compute_lambda(nu, c)
    return nu * (lam**2 + 0.25 * (np.pi * c * nu)**2)

def phi_asymptotic(x, c, method='numerical'):
    """
    Computes the asymptotic component of the scalar flux:
    phi_as(x) = exp(-|x| / nu0) / (2 * N0^+)
    """
    if c == 0.0:
        if isinstance(x, np.ndarray):
            return np.zeros_like(x)
        return 0.0
    
    nu0 = compute_nu0(c, method=method)
    N0_plus = compute_N0_plus(c, nu0)
    
    return np.exp(-np.abs(x) / nu0) / (2.0 * N0_plus)

def transient_integrand(nu, x, c):
    """
    Integrand for the transient scalar flux:
    g(nu) = exp(-|x| / nu) / N_nu
    """
    # Safeguard for nu near 0 to avoid division by zero or underflow NaN
    if nu < 1e-12:
        return 0.0
    # Safeguard for nu near 1
    if nu >= 1.0:
        return 0.0
        
    N_nu = compute_N_nu(nu, c)
    if N_nu == 0.0:
        return 0.0
        
    return np.exp(-np.abs(x) / nu) / N_nu

def phi_transient(x, c):
    """
    Computes the transient component of the scalar flux:
    phi_tr(x) = 0.5 * integral_0^1 (exp(-|x| / nu) / N_nu) dnu
    """
    if isinstance(x, (list, np.ndarray)):
        res = []
        for xi in x:
            # We integrate over (0, 1) using scipy.integrate.quad.
            # epsabs and epsrel can be configured for accuracy.
            val, _ = integrate.quad(transient_integrand, 0.0, 1.0, args=(xi, c), epsabs=1e-12, epsrel=1e-10)
            res.append(0.5 * val)
        return np.array(res)
    else:
        val, _ = integrate.quad(transient_integrand, 0.0, 1.0, args=(x, c), epsabs=1e-12, epsrel=1e-10)
        return 0.5 * val

def phi_exact(x, c, method='numerical'):
    """
    Computes the exact scalar flux:
    phi(x) = phi_as(x) + phi_tr(x)
    """
    phi_as = phi_asymptotic(x, c, method=method)
    phi_tr = phi_transient(x, c)
    return phi_as + phi_tr
