"""Question 1: the reflected sphere, in classic, asymptotic and Zimmerman diffusion."""

import numpy as np
from dataclasses import dataclass
from scipy.optimize import brentq

from homework1.criticality import extrapolation_distance, MARSHAK_EXTRAPOLATION
from homework1.exact_solution import compute_nu0_numerical, compute_nu0_magnitude_numerical

THEORIES = ('classic', 'asymptotic', 'zimmerman')

THEORY_LABELS = {
    'classic': 'Continuous classic',
    'asymptotic': 'Continuous asymptotic',
    'zimmerman': 'Discontinuous asymptotic',
}

@dataclass(frozen=True)
class Region:
    """Diffusion parameters of one medium: D0, mu0 and z0 dimensionless, rate in mfp^-1."""
    sigma_t: float
    c: float
    D0: float
    rate: float
    mu0: float
    z0: float

def relaxation_rate(c):
    """1/|nu0(c)|: the Helmholtz wavenumber above c = 1, the decay rate below it."""
    if c > 1.0:
        return 1.0 / compute_nu0_magnitude_numerical(c)
    if c < 1.0:
        return 1.0 / compute_nu0_numerical(c)
    # A pure scatterer neither grows nor decays; r phi is linear there. See explanations/04.
    return 0.0

def partial_current_factor(c, rate):
    """Zimmerman's mu0, the partial currents per unit flux; see explanations/04."""
    if rate == 0.0:
        return 0.5
    if c > 1.0:
        return c * np.log1p(rate**2) / (2.0 * rate**2)
    return -c * np.log1p(-rate**2) / (2.0 * rate**2)

def region(material, theory):
    """The Region of one benchmark material under one of THEORIES."""
    c = material.c
    if theory == 'classic':
        return Region(material.sigma_t, c, 1.0 / 3.0, np.sqrt(3.0 * abs(c - 1.0)),
                      0.5, MARSHAK_EXTRAPOLATION)

    rate = relaxation_rate(c)
    D0 = abs(c - 1.0) / rate**2 if rate else 1.0 / 3.0
    return Region(material.sigma_t, c, D0, rate, partial_current_factor(c, rate),
                  float(extrapolation_distance(c)))

def jump_ratio(core, reflector, theory):
    """phi_R/phi_C at the interface: mu0_C/mu0_R for Zimmerman, 1 for the continuous pair."""
    return core.mu0 / reflector.mu0 if theory == 'zimmerman' else 1.0

def _coth_over_length(rate, thickness):
    """kappa coth(kappa L) in reflector mfp, and its 1/L limit at c = 1 where kappa = 0."""
    return rate / np.tanh(rate * thickness) if rate else 1.0 / thickness

def _decay(rate, depth):
    """sinh(kappa s)/kappa in reflector mfp, which is s itself at c = 1."""
    return np.sinh(rate * depth) / rate if rate else depth

def _residual(R, core, reflector, thickness, g):
    """Interface balance of the two regions; its zero in (0, pi/B) is the critical radius."""
    D_C, D_R = core.D0 / core.sigma_t, reflector.D0 / reflector.sigma_t
    B = core.rate * core.sigma_t
    return (D_C * B / np.tan(B * R)
            + g * D_R * reflector.sigma_t * _coth_over_length(reflector.rate, thickness)
            - (D_C - g * D_R) / R)

def _setup(core_material, reflector_material, theory):
    """(core, reflector, jump ratio) of one pair under one theory."""
    core = region(core_material, theory)
    reflector = region(reflector_material, theory)
    return core, reflector, jump_ratio(core, reflector, theory)

def critical_radius(core_material, reflector_material, d, theory):
    """Critical core radius in cm behind d mean free paths of reflector."""
    core, reflector, g = _setup(core_material, reflector_material, theory)

    # The residual runs from +infinity at R -> 0 to -infinity at R = pi/B, the bare
    # unreflected limit, so the fundamental mode is always bracketed by that interval.
    span = np.pi / (core.rate * core.sigma_t)
    return brentq(lambda R: _residual(R, core, reflector, d + reflector.z0, g),
                  1e-6 * span, span * (1.0 - 1e-12), xtol=1e-13, rtol=8.9e-16)

def flux_profile(core_material, reflector_material, d, theory, n_points=400):
    """(r, phi) in cm across core and reflector at criticality, normalised to phi(0) = 1."""
    core, reflector, g = _setup(core_material, reflector_material, theory)
    R = critical_radius(core_material, reflector_material, d, theory)
    thickness = d + reflector.z0

    r_core = np.linspace(0.0, R, n_points)
    # sinc(x) = sin(pi x)/(pi x), which supplies the 1 at r = 0 that sin(Br)/(Br) cannot.
    phi_core = np.sinc(core.rate * core.sigma_t * r_core / np.pi)

    r_ref = np.linspace(R, R + d / reflector.sigma_t, n_points)
    depth = thickness - reflector.sigma_t * (r_ref - R)
    phi_ref = (g * phi_core[-1] * (R / r_ref)
               * _decay(reflector.rate, depth) / _decay(reflector.rate, thickness))

    return np.concatenate([r_core, r_ref]), np.concatenate([phi_core, phi_ref])
