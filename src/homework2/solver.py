"""
Numerical time-dependent diffusion solver for the planar pulse source (Question 3(c)).

Centred second differences in space on the half domain [0, L], Crank-Nicolson in
time with a backward-Euler startup. The absorption is factored out analytically, so
what is marched is the pure heat equation. See explanations/07 and 08.
"""

import numpy as np
from scipy.sparse import diags, identity
from scipy.sparse.linalg import splu
from homework2.diffusion import (
    diffusion_coefficient,
    phi_classical_diffusion,
    phi_asymptotic_diffusion,
)

ANALYTIC = {'classical': phi_classical_diffusion, 'asymptotic': phi_asymptotic_diffusion}

RANNACHER_STEPS = 4    # backward-Euler startup steps, of dt/RANNACHER_STEPS each
DOMAIN_MARGIN = 8.0    # far boundary, in Gaussian widths beyond the plotted range
N_NODES = 4000

def _grid(D, t_max, x_max, n_nodes):
    """Node-centred half domain, long enough that the pulse never reaches the far boundary."""
    length = x_max + DOMAIN_MARGIN * np.sqrt(2.0 * D * t_max)
    x = np.linspace(0.0, length, n_nodes)
    return x, x[1] - x[0]

def _time_step(target):
    """Step size near `target`, snapped so that every integer output time is an exact multiple."""
    return 1.0 / np.ceil(1.0 / target)

def _laplacian(n_nodes, h, D):
    """D d2/dx2 with reflecting ghost-node rows at both ends: symmetry at x=0, zero flux at L."""
    upper, lower = np.ones(n_nodes - 1), np.ones(n_nodes - 1)
    upper[0] = 2.0     # ghost u_{-1} = u_1
    lower[-1] = 2.0    # ghost u_{N+1} = u_{N-1}
    return diags([lower, np.full(n_nodes, -2.0), upper], [-1, 0, 1]) * (D / h**2)

def _pulse(n_nodes, h):
    """The delta as initial data, all of it in the node at the origin."""
    u = np.zeros(n_nodes)
    # Half the source on the half domain, spread over the half cell the node owns:
    # the two factors of two cancel, and 2*trapezoid(u, x) comes out 1.
    u[0] = 1.0 / h
    return u

def _gaussian(x, t, D):
    """Analytic solution of the heat equation at time t, the absorption factored out."""
    return np.exp(-x**2 / (4.0 * D * t)) / np.sqrt(4.0 * np.pi * D * t)

def _rannacher(A, u, dt):
    """Backward-Euler startup, which damps the rough modes of the delta that CN would ring on."""
    solve = splu((identity(A.shape[0]) - (dt / RANNACHER_STEPS) * A).tocsc()).solve
    for _ in range(RANNACHER_STEPS):
        u = solve(u)
    return u

def _crank_nicolson(A, u, dt, n_steps):
    """Yields u after each Crank-Nicolson step; both operators are built once."""
    unit = identity(A.shape[0])
    solve = splu((unit - 0.5 * dt * A).tocsc()).solve
    advance = (unit + 0.5 * dt * A).tocsc()
    for _ in range(n_steps):
        u = solve(advance @ u)
        yield u

def _initial_state(A, x, h, dt, D, start):
    """
    The state at t = dt, from which the Crank-Nicolson march is identical for both starts.

    'pulse' seeds the delta at t = 0 and lets the startup carry it to dt, which is what
    Question 3(c) asks for; 'warm' seeds the analytic Gaussian there instead and skips the
    startup, which removes the source-smearing error and so isolates the other two.
    """
    if start == "pulse":
        return _rannacher(A, _pulse(len(x), h), dt)
    if start == "warm":
        return _gaussian(x, dt, D)
    raise ValueError("start must be 'pulse' or 'warm'")

def solve_diffusion(c, approximation, times, x_max, n_nodes=N_NODES, dt=None, start="pulse"):
    """Numerical scalar flux at each requested time; returns (x, {t: phi}) on the half domain."""
    D = diffusion_coefficient(c, approximation)
    x, h = _grid(D, max(times), x_max, n_nodes)
    dt = _time_step(h if dt is None else dt)
    A = _laplacian(n_nodes, h, D)

    u = _initial_state(A, x, h, dt, D, start)
    wanted = {int(round(t / dt)): t for t in times}

    states = {1: u} if 1 in wanted else {}
    for step, state in enumerate(_crank_nicolson(A, u, dt, max(wanted) - 1), start=2):
        if step in wanted:
            states[step] = state

    return x, {wanted[step]: np.exp(-(1.0 - c) * wanted[step]) * state
               for step, state in states.items()}

def max_relative_error(c, approximation, times, x_max, **kwargs):
    """Largest departure of the numerical flux from the analytic one, scaled by the peak."""
    x, solution = solve_diffusion(c, approximation, times, x_max, **kwargs)
    analytic = ANALYTIC[approximation]
    return max(np.abs(phi - analytic(x, t, c)).max() / analytic(0.0, t, c)
               for t, phi in solution.items())

def particle_balance(c, approximation, times, x_max, **kwargs):
    """2 * integral(phi dx) over the half domain, divided by the expected exp(-(1-c)t)."""
    x, solution = solve_diffusion(c, approximation, times, x_max, **kwargs)
    return {t: 2.0 * np.trapezoid(phi, x) / np.exp(-(1.0 - c) * t)
            for t, phi in solution.items()}
