"""
Explicit finite-difference solution of the time-dependent diffusion equation, part 3(c).

The absorption is removed analytically and the pure heat equation is marched forward with
FTCS on the half domain, the symmetry at x = 0 entering as a ghost value. Report §3, part
3(c) derives the scheme, its stability limit and its conservation identity; the code side,
the defaults and the traps are in explanations/06.
"""

import numpy as np
from homework2.diffusion import diffusion_coefficient

# Sigma_t = v = 1 in this assignment. Kept explicit in r, in the source and in the
# prefactor, so that the expressions stay right if the units ever change.
V = 1.0

SAFETY = 0.9        # fraction of the r = 1/2 stability limit that the default step takes
TAIL_SIGMAS = 8.0   # truncation distance beyond x_max, in standard deviations
X_REACH = 1.25      # the grid covers this multiple of the causal front, as the plots do
N_NODES = 2000
WARM_T0 = 0.05      # start time of the diagnostic warm start

def _grid(D, t_max, x_max, n_nodes):
    """Node-centred half domain [0, L], truncated where the Gaussian tail is negligible."""
    length = x_max + TAIL_SIGMAS * np.sqrt(2.0 * D * V * t_max)
    x = np.linspace(0.0, length, n_nodes + 1)
    return x, x[1] - x[0]

def _time_step(D, h, safety):
    """Largest stable step from r = D v dt / h^2 <= 1/2, reduced by `safety`."""
    return safety * h**2 / (2.0 * D * V)

def _pulse(x, h):
    """The delta as v/h in the first cell; node 0 carries weight h once symmetry is counted."""
    u = np.zeros_like(x)
    u[0] = V / h
    return u

def _warm(x, t0, D):
    """Analytic Gaussian at t0 -- a diagnostic start, not the answer to 3(c)."""
    return V * np.exp(-x**2 / (4.0 * D * V * t0)) / np.sqrt(4.0 * np.pi * D * V * t0)

def step(u, r):
    """One FTCS sweep. The j = 0 row carries the factor of two from the ghost value u_-1 = u_1."""
    new = u.copy()
    new[1:-1] += r * (u[:-2] - 2.0 * u[1:-1] + u[2:])
    new[0] += 2.0 * r * (u[1] - u[0])
    new[-1] = 0.0
    return new

def _march(u, r, n_steps):
    """Applies `n_steps` sweeps at fixed r."""
    for _ in range(n_steps):
        u = step(u, r)
    return u

def solve(c, approximation='classical', times=(1.0,), start='pulse',
          n_nodes=N_NODES, safety=SAFETY):
    """Marches to each requested time; returns (x, {t: phi}) with the absorption restored."""
    D = diffusion_coefficient(c, approximation)
    x, h = _grid(D, max(times), X_REACH * max(times), n_nodes)
    dt = _time_step(D, h, safety)

    t = 0.0 if start == 'pulse' else WARM_T0
    u = _pulse(x, h) if start == 'pulse' else _warm(x, WARM_T0, D)

    fluxes = {}
    for target in sorted(times):
        # The step is shrunk to land exactly on the target, which only lowers r further;
        # overshooting by up to one step would otherwise show up as a time error.
        n_steps = int(np.ceil((target - t) / dt))
        h_t = (target - t) / n_steps
        u = _march(u, D * V * h_t / h**2, n_steps)
        t = target
        fluxes[target] = np.exp(-(1.0 - c) * V * t) * u

    return x, fluxes

def bulk(x, t, D, sigmas=3.0):
    """Mask of the region within `sigmas` standard deviations, where a relative error means something."""
    return x <= sigmas * np.sqrt(2.0 * D * V * t)

def mass(x, u):
    """int u dx over the whole line, from the half-domain values; report eq. (42)."""
    return 2.0 * np.trapezoid(u, x)
