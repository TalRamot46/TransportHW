"""
The summary numbers of part 3(d): one for the whole profile, one for the front, one for the
late-time limit the profile settles onto.

All three compare a diffusion solution to the exact transport flux without having to name a
position, which the pointwise error of the 3(d) figure cannot avoid. Report §3, part 3(d)
derives them; the code side is explanations/07.
"""

import numpy as np
from scipy.special import erf, erfc
from homework2.exact import phi_exact
from homework2.diffusion import D_CLASSICAL, diffusion_coefficient
from homework2.solver import V

# The exact flux is discontinuous at the front, so the L1 integral is taken on a grid far
# finer than the solver's, with the numeric profile interpolated onto it. Otherwise the one
# cell straddling the jump carries an O(h) error of its own.
N_FINE = 8001

def population(t, c):
    """Total particles at t: the same e^{-(1-c)t} for the exact flux and for both diffusions."""
    return np.exp(-(1.0 - c) * V * t)

def total_variation(x, phi, t, c):
    """Fraction of the population that diffusion puts in the wrong place: half the L1 distance."""
    fine = np.linspace(0.0, x[-1], N_FINE)
    difference = np.abs(np.interp(fine, x, phi) - phi_exact(fine, t, c, "series"))
    # Half of the full-line L1 distance is the whole half-line one, both being even in x.
    return np.trapezoid(difference, fine) / population(t, c)

def front_leakage(t, c, approximation="classical"):
    """Fraction of the diffusion population already past the causal front, erfc(sqrt(vt/4D))."""
    return erfc(np.sqrt(V * t / (4.0 * diffusion_coefficient(c, approximation))))

def measured_front_leakage(x, phi, t, c):
    """The same fraction read off the solver's grid, as the check on `front_leakage`."""
    beyond = x >= V * t
    return 2.0 * np.trapezoid(phi[beyond], x[beyond]) / population(t, c)

def front_in_sigmas(t, c, approximation="classical"):
    """Distance to the front in standard deviations of the diffusion packet, sqrt(vt/2D)."""
    return np.sqrt(V * t / (2.0 * diffusion_coefficient(c, approximation)))

def late_time_coefficient(c):
    """The coefficient the exact flux actually spreads with at late times, 1/(3c Sigma_t)."""
    return D_CLASSICAL / c

def gaussian_gap(first, second):
    """Half the L1 distance between two centred Gaussians of variances 2Dt. Free of t."""
    narrow, wide = sorted((first, second))
    if wide == narrow:
        return 0.0

    # The two cross at x*, and the narrower exceeds the wider inside it; both erf arguments
    # below are x*/(sigma sqrt 2), in which the t of sigma^2 = 2Dt cancels against the t in x*^2.
    ratio = wide / narrow
    log_a = 0.5 * np.log(ratio)
    return erf(np.sqrt(ratio * log_a / (ratio - 1.0))) - erf(np.sqrt(log_a / (ratio - 1.0)))

def late_time_floor(c, approximation="classical"):
    """The part of `total_variation` that never decays: wrong D against the true 1/(3c)."""
    return gaussian_gap(diffusion_coefficient(c, approximation), late_time_coefficient(c))
