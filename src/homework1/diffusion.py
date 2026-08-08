import numpy as np
from scipy.integrate import solve_ivp
from homework1.exact_solution import compute_nu0

# ---------------------------------------------------------------------------
# Diffusion parameters
#
# Both diffusion approximations solve the same equation for a unit isotropic
# plane source at the origin, in units of the mean free path (Sigma_t = 1):
#
#     -D phi''(x) + Sigma_a phi(x) = delta(x)
#
# and differ only in the diffusion coefficient:
#
#   classical    D = 1/3,             Sigma_a = 1 - c
#   asymptotic   D = (1 - c) nu0^2,   Sigma_a = 1 - c
#
# The asymptotic coefficient is fixed by requiring the closed-form solution to
# reproduce Case's asymptotic mode, phi = exp(-|x|/nu0) / (2 (1-c) nu0):
# matching the decay rate gives kappa = sqrt(Sigma_a / D) = 1/nu0, and matching
# the amplitude 1/(2 D kappa) = 1/(2 (1-c) nu0) gives D = (1-c) nu0^2.
# ---------------------------------------------------------------------------

CLASSICAL_D = 1.0 / 3.0

def diffusion_coefficients(c, approximation='classical', D=None, method='numerical'):
    """
    Returns (D, Sigma_a) for the requested diffusion approximation.

    Parameters:
    c : float
        Scattering ratio, in (0, 1) for the asymptotic approximation and
        [0, 1) for the classical one.
    approximation : str, optional
        'classical' or 'asymptotic'.
    D : float, optional
        Overrides the diffusion coefficient of the classical approximation.
        Ignored for the asymptotic approximation, whose coefficient is fixed
        by nu0.
    method : str, optional
        How nu0 is obtained for the asymptotic approximation, 'numerical' or
        'approx'. The analytic solution and the solver must be given the same
        value, otherwise the ~5e-3% difference between the two nu0 evaluations
        appears as a spurious floor when they are compared against each other.
    """
    if c < 0.0 or c >= 1.0:
        raise ValueError("Scattering ratio c must be in [0, 1).")

    sigma_a = 1.0 - c

    if approximation == 'classical':
        return (CLASSICAL_D if D is None else D), sigma_a

    if approximation == 'asymptotic':
        if c == 0.0:
            raise ValueError("The asymptotic approximation requires c > 0.")
        nu0 = compute_nu0(c, method=method)
        return sigma_a * nu0**2, sigma_a

    raise ValueError("approximation must be 'classical' or 'asymptotic'")

def phi_diffusion_analytic(x, c, approximation='classical', D=None, method='numerical'):
    """
    Closed-form Green's function of the diffusion equation for an isotropic
    plane delta-source:

        phi(x) = exp(-kappa |x|) / (2 D kappa),    kappa = sqrt(Sigma_a / D)

    which is 1 / (2 sqrt(D (1 - c))) * exp(-|x| sqrt((1 - c) / D)).
    """
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D, method)
    kappa = np.sqrt(sigma_a / D_eff)
    return np.exp(-kappa * np.abs(x)) / (2.0 * D_eff * kappa)

def phi_classical_diffusion(x, c, D=1.0/3.0):
    """
    Classical diffusion solution for an isotropic plane delta-source.
    """
    if c >= 1.0:
        raise ValueError("Scattering ratio c must be less than 1.0 for classical diffusion.")
    if c < 0.0:
        raise ValueError("Scattering ratio c must be non-negative.")

    return phi_diffusion_analytic(x, c, 'classical', D=D)

def phi_asymptotic_diffusion(x, c, method='numerical'):
    """
    Asymptotic diffusion solution for an isotropic plane delta-source:
    phi(x) = 1 / (2 * (1 - c) * nu0) * exp(-|x| / nu0)

    For c = 0, returns 0.0 (or zeros of the same shape as x), since without
    scattering there is no discrete eigenvalue and hence no asymptotic mode.
    """
    if c < 0.0 or c >= 1.0:
        raise ValueError("Scattering ratio c must be in [0, 1).")

    if c == 0.0:
        if isinstance(x, np.ndarray):
            return np.zeros_like(x, dtype=float)
        return 0.0

    nu0 = compute_nu0(c, method=method)
    if nu0 is None:
        if isinstance(x, np.ndarray):
            return np.zeros_like(x, dtype=float)
        return 0.0

    coeff = 1.0 / (2.0 * (1.0 - c) * nu0)
    return coeff * np.exp(-np.abs(x) / nu0)

# ---------------------------------------------------------------------------
# Modelling the delta source
#
# Integrating -D phi'' + Sigma_a phi = delta(x) across (-eps, +eps) gives the
# jump condition
#
#     phi'(0-) - phi'(0+) = 1 / D
#
# The Green's function of an infinite homogeneous medium is even, so
# phi'(0+) = -phi'(0-) and therefore
#
#     phi'(0+) = -1 / (2 D)     equivalently     J(0+) = -D phi'(0+) = 1/2
#
# i.e. half the source neutrons stream in each direction. This is the answer to
# the assignment's "how should a delta-function source be modelled in a
# numerical scheme?": it is not discretised at all. Symmetry turns it into a
# boundary condition, the problem is solved on the half-line x >= 0 with a
# prescribed current at x = 0, and the singularity never enters the numerics.
#
# The outer boundary is a truncation of an infinite medium, not a physical
# surface, so it carries the radiation (Robin) condition
#
#     phi'(a) = -kappa phi(a)
#
# which the decaying solution exp(-kappa x) satisfies identically. This
# annihilates the growing mode and makes the truncation exact for any a. Using
# a zero-flux condition phi(a) = 0 instead forces the flux to vanish where the
# true Green's function is small but non-zero, producing a 100% relative error
# at the outer boundary that does not shrink as the solver tolerance is reduced.
# ---------------------------------------------------------------------------

SOURCE_CURRENT = 0.5  # J(0+), from the symmetry of the Green's function

def _default_half_width(kappa, n_diffusion_lengths=10.0):
    """
    Half-width of the truncated domain, in diffusion lengths 1/kappa.
    """
    return n_diffusion_lengths / kappa

def solve_diffusion_shooting(c, approximation='classical', D=None, num_points=500,
                             n_diffusion_lengths=10.0, rtol=1e-10, atol=1e-14,
                             method='numerical', x_eval=None):
    """
    Solves the diffusion equation on the half-domain [0, a] by integrating the
    equivalent first-order system, using the symmetry-derived current condition
    at x = 0 and the radiation condition at x = a.

    The equation is linear and homogeneous away from the source, so the solution
    is exactly proportional to its starting amplitude. A single integration
    therefore suffices: integrate from a to 0 with an arbitrary amplitude, then
    rescale so that -D phi'(0) = 1/2. No root-finding is required.

    Integration runs from a towards 0, i.e. in the direction in which the
    solution grows, which is the numerically stable direction.

    Parameters:
    rtol, atol : float, optional
        Tolerances passed to the ODE integrator. These, not `num_points`,
        control the accuracy of the result; `num_points` only sets the
        resolution of the output grid.
    method : str, optional
        How nu0 is obtained for the asymptotic approximation. Must match the
        value used for any analytic solution it is compared against.
    x_eval : array_like, optional
        Output points, which must lie within [0, a]. Defaults to a uniform grid
        of `num_points` points spanning the domain. Supplying the comparison
        grid directly avoids interpolating the result afterwards.

    Returns:
    x_grid : numpy.ndarray
        The output grid, either `x_eval` or a uniform grid on [0, a].
    phi : numpy.ndarray
        Scalar flux on that grid.
    """
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D, method)
    kappa = np.sqrt(sigma_a / D_eff)
    a = _default_half_width(kappa, n_diffusion_lengths)

    # 1st-order system: dy/dx = [y2, (Sigma_a / D) y1]
    def system(x, Y):
        return [Y[1], (sigma_a / D_eff) * Y[0]]

    if x_eval is None:
        x_grid = np.linspace(0.0, a, num_points)
    else:
        x_grid = np.asarray(x_eval, dtype=float)
        if x_grid[0] < 0.0 or x_grid[-1] > a:
            raise ValueError(f"x_eval must lie within [0, {a:.4f}] for c={c}.")

    # The normalisation below reads the derivative at x = 0, so that point must
    # be evaluated even when the caller's output grid starts elsewhere.
    starts_at_origin = x_grid[0] == 0.0
    eval_points = x_grid if starts_at_origin else np.concatenate(([0.0], x_grid))

    # Start from the radiation condition at x = a with unit amplitude:
    # phi(a) = 1, phi'(a) = -kappa phi(a).
    sol = solve_ivp(system, [a, 0.0], [1.0, -kappa], t_eval=eval_points[::-1],
                    method='RK45', rtol=rtol, atol=atol)

    phi_unscaled = sol.y[0][::-1]
    dphi0_unscaled = sol.y[1][::-1][0]

    # Rescale to satisfy -D phi'(0) = J(0+) = 1/2.
    scale = (-SOURCE_CURRENT / D_eff) / dphi0_unscaled
    phi = scale * phi_unscaled

    return x_grid, phi if starts_at_origin else phi[1:]

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def absorption_balance(x, phi, c, approximation='classical', D=None):
    """
    Returns Sigma_a * integral(phi dx) over the whole line, which must equal the
    unit source strength. `x` may cover either the half-domain or the full
    domain; a half-domain profile is doubled by symmetry.

    The residual deficit at the default domain size is the physical tail of the
    Green's function beyond the truncation, not solver error.
    """
    _, sigma_a = diffusion_coefficients(c, approximation, D)
    integral = np.trapezoid(phi, x) if hasattr(np, 'trapezoid') else np.trapz(phi, x)
    if x[0] >= 0.0:
        integral *= 2.0
    return sigma_a * integral

def convergence_study(c, approximation='classical',
                      rtols=(1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11),
                      n_diffusion_lengths=10.0):
    """
    Tightens the integrator tolerance and measures the relative L2 error against
    the closed-form solution.

    For a shooting method the accuracy is set by the ODE tolerance rather than
    by a mesh spacing, so this is the meaningful analogue of a grid-refinement
    study and is what "converged numerical results" means here.

    Returns:
    rtols : numpy.ndarray
        Requested relative tolerance for each run.
    errors : numpy.ndarray
        Relative L2 error against the analytic solution for each run.
    """
    errors = []
    for rtol in rtols:
        x, phi_num = solve_diffusion_shooting(
            c, approximation, rtol=rtol, atol=1e-16,
            n_diffusion_lengths=n_diffusion_lengths)
        phi_ana = phi_diffusion_analytic(x, c, approximation)
        errors.append(np.linalg.norm(phi_num - phi_ana) / np.linalg.norm(phi_ana))

    return np.array(rtols), np.array(errors)
