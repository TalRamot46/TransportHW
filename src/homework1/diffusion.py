"""Planar diffusion for a plane source: coefficients, Green's function and solver."""

import numpy as np
from scipy.integrate import solve_ivp
from homework1.exact_solution import compute_nu0

CLASSICAL_D = 1.0 / 3.0
SOURCE_CURRENT = 0.5      # J(0+), half the source streaming each way; report §2

def diffusion_coefficients(c, approximation='classical', D=None, method='numerical'):
    """(D, Sigma_a) of one approximation; see report §1 for D_asymptotic."""
    if c < 0.0 or c >= 1.0:
        raise ValueError("Scattering ratio c must be in [0, 1).")

    sigma_a = 1.0 - c
    if approximation == 'classical':
        return (CLASSICAL_D if D is None else D), sigma_a
    if approximation == 'asymptotic':
        if c == 0.0:
            raise ValueError("The asymptotic approximation requires c > 0.")
        return sigma_a * compute_nu0(c, method=method)**2, sigma_a
    raise ValueError("approximation must be 'classical' or 'asymptotic'")

def phi_diffusion_analytic(x, c, approximation='classical', D=None, method='numerical'):
    """Diffusion Green's function exp(-kappa|x|) / (2 D kappa), kappa = sqrt(Sigma_a/D)."""
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D, method)
    kappa = np.sqrt(sigma_a / D_eff)
    return np.exp(-kappa * np.abs(x)) / (2.0 * D_eff * kappa)

def solve_diffusion_shooting(c, approximation='classical', D=None, num_points=500,
                             n_diffusion_lengths=10.0, rtol=1e-10, atol=1e-14,
                             method='numerical', x_eval=None):
    """
    Solves the half-domain [0, a] by one integration from a to 0, rescaled to J(0+) = 1/2.

    The equation is linear and homogeneous away from the source, so the solution is
    proportional to its starting amplitude and no root-finding is needed. Integrating
    inwards runs in the direction the solution grows, which is the stable one.
    """
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D, method)
    kappa = np.sqrt(sigma_a / D_eff)
    a = n_diffusion_lengths / kappa

    x_grid = np.linspace(0.0, a, num_points) if x_eval is None else np.asarray(x_eval, float)
    if x_grid[0] < 0.0 or x_grid[-1] > a:
        raise ValueError(f"x_eval must lie within [0, {a:.4f}] for c = {c}.")

    # The rescaling reads phi'(0), so x = 0 is integrated even if the caller omits it.
    at_origin = x_grid[0] == 0.0
    points = x_grid if at_origin else np.concatenate(([0.0], x_grid))

    def system(x, Y):
        return [Y[1], (sigma_a / D_eff) * Y[0]]

    # Start from the radiation condition phi'(a) = -kappa phi(a) with unit amplitude.
    sol = solve_ivp(system, [a, 0.0], [1.0, -kappa], t_eval=points[::-1],
                    method='RK45', rtol=rtol, atol=atol)

    phi = sol.y[0][::-1] * (-SOURCE_CURRENT / D_eff) / sol.y[1][-1]
    return x_grid, phi if at_origin else phi[1:]

def absorption_balance(x, phi, c, approximation='classical', D=None):
    """Sigma_a * integral(phi dx) over the line, which must equal the unit source."""
    _, sigma_a = diffusion_coefficients(c, approximation, D)
    integral = np.trapezoid(phi, x)
    return sigma_a * (2.0 * integral if x[0] >= 0.0 else integral)

def convergence_study(c, approximation='classical',
                      rtols=(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11)):
    """Relative L2 error against the closed form as the integrator tolerance tightens."""
    errors = []
    for rtol in rtols:
        x, phi = solve_diffusion_shooting(c, approximation, rtol=rtol, atol=1e-16)
        exact = phi_diffusion_analytic(x, c, approximation)
        errors.append(np.linalg.norm(phi - exact) / np.linalg.norm(exact))

    return np.array(rtols), np.array(errors)
