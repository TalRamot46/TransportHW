"""Case's exact scalar flux for a plane source in a slab, and its discrete eigenvalue."""

import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad

def _fit_exponent(c):
    """Exponent p(c) of the analytic fit nu0 = 1 / sqrt(1 - c^p(c))."""
    return 2.47412 + 0.00363081 / c**2 - 0.0352458 * c + 0.557498 / c

def compute_nu0_numerical(c):
    """Discrete eigenvalue nu0 (c < 1) from the root of arctanh(k0)/k0 = 1/c."""
    if c <= 0.0:
        return None
    if c >= 1.0:
        raise ValueError("nu0 is real only for c < 1.")

    def f(k0):
        return 1.0 - 1.0 / c if k0 == 0.0 else np.arctanh(k0) / k0 - 1.0 / c

    return 1.0 / brentq(f, 0.0, 1.0 - 1e-15)

def compute_nu0_approx(c):
    """Discrete eigenvalue nu0 (c < 1) from the analytic fit."""
    if c <= 0.0:
        return None
    if c >= 1.0:
        raise ValueError("nu0 is real only for c < 1.")
    return 1.0 / np.sqrt(1.0 - c**_fit_exponent(c))

def compute_nu0(c, method='numerical'):
    """nu0 by the requested route, 'numerical' or 'approx'."""
    return compute_nu0_numerical(c) if method == 'numerical' else compute_nu0_approx(c)

def compute_nu0_magnitude_numerical(c):
    """|nu0| for c > 1 from the root of c arctan(k0) = k0; see report §3."""
    if c <= 1.0:
        raise ValueError("The imaginary eigenvalue branch requires c > 1.")

    # f(0) = 0 with f'(0) = c - 1 > 0 and f(c pi/2) < 0, so this brackets the
    # single positive root without catching the trivial one at the origin.
    def f(k0):
        return c * np.arctan(k0) - k0

    return 1.0 / brentq(f, 1e-12, c * np.pi / 2.0, xtol=1e-15, rtol=8.9e-16)

def compute_nu0_magnitude_approx(c):
    """|nu0| for c > 1 from the same fit, whose radicand turns negative there."""
    if c <= 1.0:
        raise ValueError("The imaginary eigenvalue branch requires c > 1.")
    return 1.0 / np.sqrt(c**_fit_exponent(c) - 1.0)

def compute_nu0_magnitude(c, method='approx'):
    """|nu0| for c > 1 by the requested route, 'numerical' or 'approx'."""
    if method == 'numerical':
        return compute_nu0_magnitude_numerical(c)
    return compute_nu0_magnitude_approx(c)

def compute_N0_plus(c, nu0):
    """Normalisation of the discrete mode, N0+."""
    return 0.5 * c * nu0**3 * (c / (nu0**2 - 1.0) - 1.0 / nu0**2)

def compute_lambda(nu, c):
    """Dispersion function lambda(nu) = 1 - c nu arctanh(nu) on [0, 1]."""
    if nu == 0.0:
        return 1.0
    if nu >= 1.0:
        return -np.inf
    return 1.0 - c * nu * np.arctanh(nu)

def compute_N_nu(nu, c):
    """Normalisation of the continuum mode, N_nu."""
    if nu == 0.0:
        return 0.0
    return nu * (compute_lambda(nu, c)**2 + 0.25 * (np.pi * c * nu)**2)

def phi_asymptotic(x, c, method='numerical'):
    """Asymptotic component, exp(-|x|/nu0) / (2 N0+); zero at c = 0."""
    if c == 0.0:
        return np.zeros_like(x, dtype=float) if np.ndim(x) else 0.0

    nu0 = compute_nu0(c, method=method)
    return np.exp(-np.abs(x) / nu0) / (2.0 * compute_N0_plus(c, nu0))

def _transient_integrand(nu, x, c):
    """exp(-|x|/nu) / N_nu, cut off at the endpoints where N_nu is 0 or infinite."""
    if nu < 1e-12 or nu >= 1.0:
        return 0.0
    return np.exp(-np.abs(x) / nu) / compute_N_nu(nu, c)

def phi_transient(x, c):
    """Transient component, (1/2) integral over nu in (0, 1) of the continuum modes."""
    values = [0.5 * quad(_transient_integrand, 0.0, 1.0, args=(xi, c),
                         epsabs=1e-12, epsrel=1e-10)[0]
              for xi in np.atleast_1d(x)]
    return np.array(values) if np.ndim(x) else values[0]

def phi_exact(x, c, method='numerical'):
    """Exact scalar flux, the sum of the asymptotic and transient components."""
    return phi_asymptotic(x, c, method=method) + phi_transient(x, c)
