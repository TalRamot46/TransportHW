"""The P_N algebra shared by both solutions of Question 2: the truncated Legendre
recurrence and the Marshak conditions at the vacuum face."""

import numpy as np

def streaming_matrix(N):
    """The matrix of report eq. (14): A[n, n+1] = (n+1)/(2n+1), A[n, n-1] = n/(2n+1)."""
    n = np.arange(N + 1)
    matrix = np.zeros((N + 1, N + 1))
    matrix[n[:-1], n[:-1] + 1] = (n[:-1] + 1) / (2 * n[:-1] + 1)
    matrix[n[1:], n[1:] - 1] = n[1:] / (2 * n[1:] + 1)
    return matrix

def parity_indices(N):
    """The even and odd moment indices, ascending."""
    return np.arange(0, N + 1, 2), np.arange(1, N + 1, 2)

def parity_blocks(N):
    """The even-to-odd and odd-to-even blocks A and B of report eq. (22)."""
    matrix = streaming_matrix(N)
    even, odd = parity_indices(N)
    return matrix[np.ix_(even, odd)], matrix[np.ix_(odd, even)]

def marshak_matrix(N):
    """The (N+1)/2 Marshak rows c[m, n] of report eq. (17), one row per condition."""
    nodes, weights = np.polynomial.legendre.leggauss(N + 1)
    mu, weights = 0.5 * (nodes + 1.0), 0.5 * weights        # the half range [0, 1]
    powers = np.array([mu ** (2 * m - 1) for m in range(1, (N + 3) // 2)])

    # Exact, not approximate: the integrand has degree at most 2N and the rule has N+1 nodes.
    half_range = (powers * weights) @ np.polynomial.legendre.legvander(mu, N)

    n = np.arange(N + 1)
    return (-1.0) ** (n + 1) * (2 * n + 1) / 2 * half_range
