"""Question 1: the reflected sphere, in classic, asymptotic and Zimmerman diffusion."""

import numpy as np
from dataclasses import dataclass
from scipy.optimize import brentq

from homework1.criticality import extrapolation_distance, MARSHAK_EXTRAPOLATION
from homework1.exact_solution import compute_nu0_numerical, compute_k0_numerical

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
    rate: float # nu0 or k0 depending on c < 1 or c > 1 
    mu0: float
    z0: float

def relaxation_rate(c):
    """1/nu0(c) or 1/k0(c)."""
    if c > 1.0:
        return 1.0 / compute_k0_numerical(c)
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
    """The mu0 ratio in front of the reflector in report eq. (13); 1 for the continuous pair."""
    return core.mu0 / reflector.mu0 if theory == 'zimmerman' else 1.0

# Report eqs. (8) and (13), each written as cot(k0 a) minus its own right-hand side, so that
# the critical a is a zero of whichever of the two applies.

def _cot(x):
    """cot(x), which numpy does not provide."""
    return 1.0 / np.tan(x)

def _far_face_return(reflector, thickness):
    """coth([d + z0]/nu0)/nu0 in reflector mfp, and its 1/[d + z0] limit at c = 1."""
    rate = reflector.rate
    return rate / np.tanh(rate * thickness) if rate else 1.0 / thickness

def _leakage(core, reflector):
    """D_R/(D_C k0), the factor both equations put in front of the reflector."""
    return reflector.D0 / (core.D0 * core.rate)

def _interface_in_reflector_mfp(a, core, reflector):
    """b - d: the interface radius a, counted in reflector mean free paths instead of core."""
    return a * reflector.sigma_t / core.sigma_t

def _continuous(a, core, reflector, thickness):
    """Report eq. (8): phi and the net current continuous, theories (a) and (b)."""
    k0a = core.rate * a
    curvature = 1.0 / _interface_in_reflector_mfp(a, core, reflector)
    return (_cot(k0a) - 1.0 / k0a
            + _leakage(core, reflector) * (_far_face_return(reflector, thickness) + curvature))

def _fixed(a, core, reflector, thickness, mu_ratio):
    """Report eq. (13): the jump imposed on r phi with j_2/1 = 1, theory (c)."""
    return (_cot(core.rate * a)
            + mu_ratio * _leakage(core, reflector) * _far_face_return(reflector, thickness))

def _residual(a, core, reflector, thickness, mu_ratio, theory):
    """Whichever of the two criticality equations this theory is solved from."""
    if theory == 'zimmerman':
        return _fixed(a, core, reflector, thickness, mu_ratio)
    return _continuous(a, core, reflector, thickness)

def _decay(rate, depth):
    """nu0 sinh(depth/nu0) in reflector mfp, which is depth itself at c = 1."""
    return np.sinh(rate * depth) / rate if rate else depth

def _setup(core_material, reflector_material, theory):
    """(core, reflector, mu0_C/mu0_R) of one pair under one theory."""
    core = region(core_material, theory)
    reflector = region(reflector_material, theory)
    return core, reflector, jump_ratio(core, reflector, theory)

def critical_radius(core_material, reflector_material, d, theory):
    """Critical core radius in cm behind d mean free paths of reflector."""
    core, reflector, mu_ratio = _setup(core_material, reflector_material, theory)

    # Both residuals run from +infinity at a -> 0 to -infinity at a = pi/k0, the bare
    # unreflected limit, so the fundamental mode is always bracketed by that interval.
    span = np.pi / core.rate
    a = brentq(lambda a: _residual(a, core, reflector, d + reflector.z0, mu_ratio, theory),
               1e-6 * span, span * (1.0 - 1e-12), xtol=1e-13, rtol=8.9e-16)
    return a / core.sigma_t

def flux_profile(core_material, reflector_material, d, theory, n_points=400):
    """(r, phi) in cm across core and reflector at criticality, normalised to phi(0) = 1."""
    core, reflector, mu_ratio = _setup(core_material, reflector_material, theory)
    R = critical_radius(core_material, reflector_material, d, theory)
    thickness = d + reflector.z0

    r_core = np.linspace(0.0, R, n_points)
    # sinc(x) = sin(pi x)/(pi x), which supplies the 1 at r = 0 that sin(k0 a)/(k0 a) cannot.
    phi_core = np.sinc(core.rate * core.sigma_t * r_core / np.pi)

    r_ref = np.linspace(R, R + d / reflector.sigma_t, n_points)
    depth = thickness - reflector.sigma_t * (r_ref - R)
    phi_ref = (mu_ratio * phi_core[-1] * (R / r_ref)
               * _decay(reflector.rate, depth) / _decay(reflector.rate, thickness))

    return np.concatenate([r_core, r_ref]), np.concatenate([phi_core, phi_ref])
