"""
Classical and asymptotic diffusion solutions for the same planar pulse source.

Both solve  dn/dt - D n'' + (1-c) n = delta(x) delta(t)  with Sigma_t = v = 1, and
differ only through D. See explanations/03-diffusion-coefficients.md.
"""

import numpy as np
from homework1.exact_solution import compute_nu0_numerical, compute_nu0_magnitude_numerical

# Slab-geometry classical value, D = 1/(3 Sigma_t).
D_CLASSICAL = 1.0 / 3.0

def nu0_squared(c):
    """
    Signed square of the discrete transport eigenvalue: nu0^2 for c < 1, -|nu0|^2 for c > 1.

    The sign carries the fact that nu0 turns imaginary above c = 1, which is what
    keeps D0 = (1-c) nu0^2 positive on both sides.
    """
    if c < 1.0:
        return compute_nu0_numerical(c) ** 2
    if c > 1.0:
        return -compute_nu0_magnitude_numerical(c) ** 2
    return np.inf

def diffusion_coefficient(c, approximation="classical"):
    """D of the requested approximation: 1/3 for classical, D0(c) = (1-c) nu0^2 for asymptotic."""
    if approximation == "classical":
        return D_CLASSICAL
    if approximation != "asymptotic":
        raise ValueError("approximation must be 'classical' or 'asymptotic'")

    # At c = 1 there is no discrete eigenvalue; the limit of (1-c) nu0^2 from either side is 1/3.
    return D_CLASSICAL if c == 1.0 else (1.0 - c) * nu0_squared(c)

def _phi_diffusion(x, t, c, D):
    """Green's function of the diffusion equation: a Gaussian decaying at the absorption rate."""
    x = np.asarray(x, dtype=float)
    return np.exp(-(1.0 - c) * t) / np.sqrt(4.0 * np.pi * D * t) * np.exp(-x**2 / (4.0 * D * t))

def phi_classical_diffusion(x, t, c):
    """Time-dependent classical diffusion flux, D = 1/3."""
    return _phi_diffusion(x, t, c, D_CLASSICAL)

def phi_asymptotic_diffusion(x, t, c):
    """Time-dependent asymptotic diffusion flux, D = D0(c) = (1-c) nu0^2."""
    return _phi_diffusion(x, t, c, diffusion_coefficient(c, "asymptotic"))

def _phi_steady(x, c, D):
    """Steady planar Green's function e^{-kappa|x|}/(2 D kappa), the time-integral of _phi_diffusion."""
    if c >= 1.0:
        raise ValueError("The steady solution exists only for c < 1; above it the flux grows without bound.")
    kappa = np.sqrt((1.0 - c) / D)
    return np.exp(-kappa * np.abs(np.asarray(x, dtype=float))) / (2.0 * D * kappa)

def phi_steady_classical(x, c):
    """Steady classical diffusion flux, sqrt(3)/(2 sqrt(1-c)) e^{-sqrt(3(1-c))|x|}."""
    return _phi_steady(x, c, D_CLASSICAL)

def phi_steady_asymptotic(x, c):
    """Steady asymptotic diffusion flux, e^{-|x|/nu0} / (2(1-c) nu0)."""
    return _phi_steady(x, c, diffusion_coefficient(c, "asymptotic"))
