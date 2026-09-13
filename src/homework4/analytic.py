"""Question 5: closed-form critical radii of a bare sphere, for comparison with Monte Carlo."""

import numpy as np
from homework1.criticality import critical_dimensions

RADIUS = 3          # index of Sigma_t R_c in the critical_dimensions tuple

def classical_radius(c):
    """Classical diffusion, D = 1/3 with z0 = 2D: R = pi/sqrt(3(c-1)) - 2/3, in mfp."""
    return critical_dimensions(np.atleast_1d(c), 'marshak')[RADIUS]

def milne_radius(c):
    """Asymptotic diffusion: D0 = (c-1)|nu0|^2 with Case's exact Milne z0(c), in mfp."""
    return critical_dimensions(np.atleast_1d(c), 'transport-ref')[RADIUS]

# The two bare-sphere theories of Question 5, in the order they are plotted. Both are
# R = pi/B - z0; they differ in the pair (B, z0). See report §5.
THEORIES = (('Classical diffusion', classical_radius),
            ('Asymptotic Milne (exact)', milne_radius))

def critical_mass(radius_cm, rho):
    """Mass in grams of a sphere of that radius at density rho."""
    return rho * (4.0 / 3.0) * np.pi * radius_cm**3