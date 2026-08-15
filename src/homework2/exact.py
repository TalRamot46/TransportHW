"""
Exact planar-source scalar flux of Assignment 2, Question 2.

Units are those of the assignment throughout: Sigma_t = v = 1, so the Paasschens
argument a = v Sigma_t t reduces to t, and the flux equals the number density.
Derivation: report §2. The numerically stable form: explanations/02.
"""

import numpy as np
from scipy.special import dawsn, gammaln

# Paasschens' interpolation constant, G(w) ~ e^w sqrt(1 + b/w).
B_INTERP = 2.026

# Terms of the exact series for G; the ratio of factorials decays fast enough
# that this many terms is converged to machine precision for w0 <= 30.
SERIES_TERMS = 200

def _series_log_coefficients(n):
    """Log of Gamma(3N/4 + 3/2) / Gamma(3N/4) / N! for N = 1..n, for the exact G series."""
    N = np.arange(1, n + 1, dtype=float)
    return gammaln(0.75 * N + 1.5) - gammaln(0.75 * N) - gammaln(N + 1.0)

def _safe_log(w):
    """log(w) with the zeros replaced by 1, so that w^N is formed without a divide warning."""
    return np.log(np.where(w > 0.0, w, 1.0))

def G(w, form="interp"):
    """Paasschens' G: the interpolation e^w sqrt(1 + b/w), or the exact series."""
    w = np.asarray(w, dtype=float)
    if form == "interp":
        return np.exp(w) * np.sqrt(1.0 + B_INTERP / w)
    if form != "series":
        raise ValueError("form must be 'interp' or 'series'")

    N = np.arange(1, SERIES_TERMS + 1, dtype=float)
    terms = np.exp(_series_log_coefficients(SERIES_TERMS) + _safe_log(w)[..., None] * N)
    return 8.0 * (3.0 * w) ** -1.5 * terms.sum(axis=-1)

def _dawson_antiderivative(s):
    """The bounded part of int e^w sqrt(w+b) dw, i.e. sqrt(s) - D(sqrt(s)) at s = w + b."""
    root = np.sqrt(s)
    return root - dawsn(root)

def collided_integral(w0, form="interp"):
    """
    int_0^w0 sqrt(w) G(w) dw, scaled by e^{-w0} to keep it bounded.

    For the interpolation this is elementary, since sqrt(w) G(w) = e^w sqrt(w+b);
    the Dawson form is used instead of erfi because erfi overflows past w0 ~ 20.
    """
    w0 = np.asarray(w0, dtype=float)
    if form == "interp":
        return _dawson_antiderivative(w0 + B_INTERP) - np.exp(-w0) * _dawson_antiderivative(B_INTERP)

    # Exact series, integrated term by term: sum_N coeff_N w0^N / N.
    N = np.arange(1, SERIES_TERMS + 1, dtype=float)
    log_terms = _series_log_coefficients(SERIES_TERMS) - np.log(N) + _safe_log(w0)[..., None] * N
    total = 8.0 / 3.0 ** 1.5 * np.exp(log_terms - w0[..., None]).sum(axis=-1)
    return np.where(w0 > 0.0, total, 0.0)

def front_distance(x, t):
    """w0 = t (1 - x^2/t^2)^{3/4}, the Paasschens variable at the field point, 0 at the front."""
    return t * np.maximum(1.0 - (x / t) ** 2, 0.0) ** 0.75

def phi_c1(x, t, form="interp"):
    """
    Planar-source scalar flux for c = 1, phi(x,t;1), from the Q2 closed form.

    The leading e^{-t} is the uncollided plateau; the bracket adds the collided
    part, measured in units of that plateau.
    """
    x = np.asarray(x, dtype=float)
    inside = np.abs(x) < t

    w0 = front_distance(np.where(inside, x, 0.0), t)
    collided = np.sqrt(3.0 / np.pi) * np.exp(-(t - w0)) * collided_integral(w0, form)

    return np.where(inside, (np.exp(-t) + collided) / (2.0 * t), 0.0)

def phi_exact(x, t, c, form="interp"):
    """Scalar flux for any c, from the Q1 scaling identity phi(x,t;c) = c e^{-(1-c)t} phi(cx,ct;1)."""
    return c * np.exp(-(1.0 - c) * t) * phi_c1(c * np.asarray(x, dtype=float), c * t, form)
