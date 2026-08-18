"""Question 2, Method 2: the same P_N system solved in closed form, criticality from
det H(a) = 0. Carries no spatial error, so it is the benchmark for Method 1."""

import numpy as np
from scipy.optimize import brentq
from homework3 import pn

SCAN_POINTS = 200

def modes(c, N):
    """Eigenpairs of K^2 = -(AB)^-1 Sigma_0 at k = 1, and the odd-to-even block B;
    report eq. (23)."""
    A, B = pn.parity_blocks(N)
    sigma_0 = np.eye(len(A))
    sigma_0[0, 0] = 1.0 - c

    values, vectors = np.linalg.eig(-np.linalg.solve(A @ B, sigma_0))

    # K^2 is similar to a symmetric matrix, so any imaginary part is round-off.
    values, vectors = values.real, vectors.real
    if values.max() <= 0.0:
        raise ValueError("No oscillatory mode: a bare slab is critical only for c > 1.")
    return values, vectors, B

def _boundary_matrix(a, values, vectors, B, m_even, m_odd):
    """H(a) of report eq. (26), one column per mode."""
    columns = []
    for value, v in zip(values, vectors.T):
        flux, current = m_even @ v, m_odd @ (B @ v)
        if value > 0.0:
            root = np.sqrt(value)
            columns.append(flux * np.cos(root * a) + current * root * np.sin(root * a))
        else:
            # Divided through by cosh(root a), which cannot move a zero of the
            # determinant and keeps the boundary-layer columns from overflowing.
            root = np.sqrt(-value)
            columns.append(flux - current * root * np.tanh(root * a))
    return np.column_stack(columns)

def _first_root(determinant, upper):
    """First sign change of det H on (0, upper), refined by brentq."""
    grid = np.linspace(0.0, upper, SCAN_POINTS)
    values = np.array([determinant(a) for a in grid])

    crossings = np.flatnonzero(values[:-1] * values[1:] < 0.0)
    if len(crossings) == 0:
        raise RuntimeError("det H(a) does not change sign below the fundamental cutoff.")

    first = crossings[0]
    return brentq(determinant, grid[first], grid[first + 1], xtol=1e-13)

def critical_half_thickness(medium, N):
    """Half-thickness at which det H(a) = 0, the smallest positive root."""
    values, vectors, B = modes(medium.c, N)
    marshak = pn.marshak_matrix(N)
    even, odd = pn.parity_indices(N)
    m_even, m_odd = marshak[:, even], marshak[:, odd]

    def determinant(a):
        return np.linalg.det(_boundary_matrix(a, values, vectors, B, m_even, m_odd))

    # The critical slab is always thinner than the first zero of its fundamental cosine.
    cutoff = 0.5 * np.pi / np.sqrt(values.max())
    return _first_root(determinant, cutoff) / medium.sigma_t
