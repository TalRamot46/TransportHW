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

def compute_nu0_approx(c):
    """
    Computes the discrete eigenvalue nu0 using the high-accuracy analytic approximation:
    nu0 approx 1 / sqrt(1 - c^p(c))
    """
    if c <= 0.0:
        return None
    if c >= 1.0:
        raise ValueError("Scattering ratio c must be less than 1.0.")
    
    p = 2.47412 + 0.00363081 / (c**2) - 0.0352458 * c + 0.557498 / c
    return 1.0 / np.sqrt(1.0 - c**p)

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
