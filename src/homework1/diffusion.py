import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from scipy.optimize import brentq
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

def diffusion_coefficients(c, approximation='classical', D=None):
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
    """
    if c < 0.0 or c >= 1.0:
        raise ValueError("Scattering ratio c must be in [0, 1).")

    sigma_a = 1.0 - c

    if approximation == 'classical':
        return (CLASSICAL_D if D is None else D), sigma_a

    if approximation == 'asymptotic':
        if c == 0.0:
            raise ValueError("The asymptotic approximation requires c > 0.")
        nu0 = compute_nu0(c)
        return sigma_a * nu0**2, sigma_a

    raise ValueError("approximation must be 'classical' or 'asymptotic'")

def phi_diffusion_analytic(x, c, approximation='classical', D=None):
    """
    Closed-form Green's function of the diffusion equation for an isotropic
    plane delta-source:

        phi(x) = exp(-kappa |x|) / (2 D kappa),    kappa = sqrt(Sigma_a / D)

    which is 1 / (2 sqrt(D (1 - c))) * exp(-|x| sqrt((1 - c) / D)).
    """
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D)
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
# i.e. half the source neutrons stream in each direction. This turns the delta
# into a *boundary condition*: the problem is solved on the half-line x >= 0
# with a prescribed current at x = 0, and the singularity never enters the
# discretisation. `solve_diffusion_shooting` and `solve_diffusion_fv` below use
# this. `solve_diffusion_fv_full` instead keeps the whole line and smears the
# delta over the single cell containing the origin, for comparison.
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
# at the outer boundary that does not shrink under mesh refinement.
# ---------------------------------------------------------------------------

def _source_current():
    """
    Current at the source-side boundary, J(0+) = 1/2, from the symmetry of the
    Green's function. Independent of D and c.
    """
    return 0.5

def _default_half_width(kappa, n_diffusion_lengths=10.0):
    """
    Half-width of the truncated domain, in diffusion lengths 1/kappa.
    """
    return n_diffusion_lengths / kappa

def solve_diffusion_shooting(c, approximation='classical', D=None, num_points=500,
                             n_diffusion_lengths=10.0):
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

    Returns:
    x_grid : numpy.ndarray
        Uniform grid on [0, a].
    phi : numpy.ndarray
        Scalar flux on that grid.
    """
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D)
    kappa = np.sqrt(sigma_a / D_eff)
    a = _default_half_width(kappa, n_diffusion_lengths)

    # 1st-order system: dy/dx = [y2, (Sigma_a / D) y1]
    def system(x, Y):
        return [Y[1], (sigma_a / D_eff) * Y[0]]

    x_grid = np.linspace(0.0, a, num_points)

    # Start from the radiation condition at x = a with unit amplitude:
    # phi(a) = 1, phi'(a) = -kappa phi(a).
    sol = solve_ivp(system, [a, 0.0], [1.0, -kappa], t_eval=x_grid[::-1],
                    method='RK45', rtol=1e-10, atol=1e-14)

    phi_unscaled = sol.y[0][::-1]
    dphi0_unscaled = sol.y[1][-1]

    # Rescale to satisfy -D phi'(0) = J(0+) = 1/2.
    dphi0_target = -_source_current() / D_eff
    scale = dphi0_target / dphi0_unscaled

    return x_grid, scale * phi_unscaled

def _solve_tridiagonal(lower, diag, upper, rhs):
    """
    Solves a tridiagonal system given its three diagonals.
    `lower[i]` multiplies unknown i-1 in row i, `upper[i]` multiplies unknown
    i+1 in row i; `lower[0]` and `upper[-1]` are unused.
    """
    n = len(diag)
    ab = np.zeros((3, n))
    ab[0, 1:] = upper[:-1]
    ab[1, :] = diag
    ab[2, :-1] = lower[1:]
    return solve_banded((1, 1), ab, rhs)

def _robin_face_coefficient(D_eff, kappa, dx):
    """
    Coefficient k such that the outward current at a radiation boundary is
    k * phi_edge_cell.

    With a ghost cell the Robin condition phi' = -kappa phi, discretised with
    a centred difference and a centred average at the face, gives
    phi_ghost = phi_edge (2 - kappa dx) / (2 + kappa dx), and hence an outward
    current of 2 D kappa / (2 + kappa dx) times the edge-cell flux.
    """
    return 2.0 * D_eff * kappa / (2.0 + kappa * dx)

def solve_diffusion_fv(c, approximation='classical', D=None, n_cells=500,
                       n_diffusion_lengths=10.0):
    """
    Cell-centred finite-volume solution on the half-domain [0, a].

    The delta source enters only through the face current J(0) = 1/2 at the
    left boundary, so the kink of the Green's function sits exactly on a face
    and is represented exactly. The right boundary uses the radiation
    condition.

    Balance over cell i, with faces i and i+1:
        J_{i+1} - J_i + Sigma_a dx phi_i = 0,   J = -D dphi/dx

    Returns:
    x_centers : numpy.ndarray
        Cell-centre coordinates on [0, a].
    phi : numpy.ndarray
        Cell-averaged scalar flux.
    """
    D_eff, sigma_a = diffusion_coefficients(c, approximation, D)
    kappa = np.sqrt(sigma_a / D_eff)
    a = _default_half_width(kappa, n_diffusion_lengths)

    dx = a / n_cells
    x_centers = (np.arange(n_cells) + 0.5) * dx

    t = D_eff / dx  # face transmission coefficient
    diag = np.full(n_cells, 2.0 * t + sigma_a * dx)
    lower = np.full(n_cells, -t)
    upper = np.full(n_cells, -t)
    rhs = np.zeros(n_cells)

    # Left boundary (x = 0): prescribed current from the source symmetry.
    diag[0] = t + sigma_a * dx
    rhs[0] = _source_current()

    # Right boundary (x = a): radiation condition.
    diag[-1] = t + sigma_a * dx + _robin_face_coefficient(D_eff, kappa, dx)

    return x_centers, _solve_tridiagonal(lower, diag, upper, rhs)

def solve_diffusion_fv_full(c, approximation='classical', D=None, n_cells=501,
                            n_diffusion_lengths=10.0):
    """
    Cell-centred finite-volume solution on the full domain [-a, a], with the
    delta source smeared over the single cell containing the origin.

    This is the alternative answer to "how should a delta-function source be
    modelled in a numerical scheme?": rather than converting the source into a
    boundary condition, keep the whole line and set the cell-integrated source
    to 1 in the origin cell. It generalises to sources that are not delta
    functions, at the cost of smearing the kink at x = 0 over one cell.

    `n_cells` is forced odd so that the origin lies at a cell centre and the
    discretisation stays symmetric.

    Returns:
    x_centers : numpy.ndarray
        Cell-centre coordinates on [-a, a].
    phi : numpy.ndarray
        Cell-averaged scalar flux.
    """
    if n_cells % 2 == 0:
        n_cells += 1

    D_eff, sigma_a = diffusion_coefficients(c, approximation, D)
    kappa = np.sqrt(sigma_a / D_eff)
    a = _default_half_width(kappa, n_diffusion_lengths)

    dx = 2.0 * a / n_cells
    x_centers = -a + (np.arange(n_cells) + 0.5) * dx

    t = D_eff / dx
    diag = np.full(n_cells, 2.0 * t + sigma_a * dx)
    lower = np.full(n_cells, -t)
    upper = np.full(n_cells, -t)
    rhs = np.zeros(n_cells)

    # Radiation condition at both outer boundaries.
    robin = _robin_face_coefficient(D_eff, kappa, dx)
    diag[0] = t + sigma_a * dx + robin
    diag[-1] = t + sigma_a * dx + robin

    # Cell-integrated delta source: the whole unit source in the origin cell.
    rhs[n_cells // 2] = 1.0

    return x_centers, _solve_tridiagonal(lower, diag, upper, rhs)

def solve_diffusion_numerical(c, D=1.0/3.0, num_points=500):
    """
    Original shooting solution on [-a, 0] with a zero-flux outer boundary,
    retained for reference and for the existing Question 2 figure.

    Prefer `solve_diffusion_shooting` or `solve_diffusion_fv`: the zero-flux
    condition phi(-a) = 0 used here is a truncation artifact, and forces a 100%
    relative error at x = -a where the true Green's function is small but
    non-zero. Its effect on the interior is only about -2e-9 relative for the
    default domain size.
    """
    if c >= 1.0 or c <= 0.0:
        raise ValueError("c must be in (0, 1) for the shooting method.")

    kappa = np.sqrt((1.0 - c) / D)
    a = 10.0 / kappa  # 10 diffusion lengths

    # 1st-order system: dy/dx = [y2, (1-c)/D * y1]
    def system(x, Y):
        y1, y2 = Y
        return [y2, ((1.0 - c) / D) * y1]

    # Objective function for shooting: f(s) = y2(0; s) - 1/(2D)
    target = 1.0 / (2.0 * D)

    def shoot(s):
        # Integrate from -a to 0
        sol = solve_ivp(system, [-a, 0.0], [0.0, s], method='RK45', rtol=1e-10, atol=1e-12)
        # Return value of y2 at the end (x=0) minus target
        return sol.y[1, -1] - target

    # Bracket the root. s_min = 0.0 gives shoot(s_min) = -target < 0.
    # Since exact s = target * exp(-kappa * a), target * 1.5 is a safe upper bound.
    s_min = 0.0
    s_max = target * 1.5

    # Solve for shooting parameter s
    s_opt = brentq(shoot, s_min, s_max)

    # Re-integrate with optimal shooting parameter to get dense output on a uniform grid
    x_grid = np.linspace(-a, 0.0, num_points)
    sol = solve_ivp(system, [-a, 0.0], [0.0, s_opt], t_eval=x_grid, method='RK45', rtol=1e-10, atol=1e-12)

    return x_grid, sol.y[0]

# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def absorption_balance(x, phi, c, approximation='classical', D=None,
                       quadrature='trapezoid'):
    """
    Returns Sigma_a * integral(phi dx) over the whole line, which must equal the
    unit source strength. `x` may cover either the half-domain or the full
    domain; a half-domain profile is doubled by symmetry.

    Parameters:
    quadrature : str
        'trapezoid' for point values sampled on a grid (the shooting solvers),
        'midpoint' for cell averages on cell centres (the finite-volume
        solvers). The distinction matters: a trapezoid rule over cell centres
        omits the half-cell slivers at each end of the domain, which for the
        half-domain solver discards a fraction kappa dx / 2 of the integral --
        1% at the default resolution, purely as a quadrature artifact.
    """
    _, sigma_a = diffusion_coefficients(c, approximation, D)

    if quadrature == 'midpoint':
        integral = np.sum(phi) * (x[1] - x[0])
    elif quadrature == 'trapezoid':
        integral = np.trapezoid(phi, x) if hasattr(np, 'trapezoid') else np.trapz(phi, x)
    else:
        raise ValueError("quadrature must be 'trapezoid' or 'midpoint'")

    if x[0] >= 0.0:
        integral *= 2.0
    return sigma_a * integral

def convergence_study(c, approximation='classical', solver='fv',
                      n_cells_list=(25, 50, 100, 200, 400, 800, 1600),
                      n_diffusion_lengths=10.0):
    """
    Refines the mesh and measures the relative L2 error against the closed-form
    solution, returning the observed order of accuracy between successive
    refinements.

    Parameters:
    solver : str
        'fv' for the half-domain finite-volume solver, 'fv_full' for the
        full-domain solver with a smeared source.

    Returns:
    dx : numpy.ndarray
        Mesh spacing for each refinement.
    errors : numpy.ndarray
        Relative L2 error for each refinement.
    orders : numpy.ndarray
        Observed convergence order between successive refinements; one shorter
        than `dx`.
    """
    solvers = {'fv': solve_diffusion_fv, 'fv_full': solve_diffusion_fv_full}
    if solver not in solvers:
        raise ValueError("solver must be 'fv' or 'fv_full'")
    solve = solvers[solver]

    dx_values = []
    errors = []
    for n_cells in n_cells_list:
        x, phi_num = solve(c, approximation, n_cells=n_cells,
                           n_diffusion_lengths=n_diffusion_lengths)
        phi_ana = phi_diffusion_analytic(x, c, approximation)
        dx_values.append(x[1] - x[0])
        errors.append(np.linalg.norm(phi_num - phi_ana) / np.linalg.norm(phi_ana))

    dx_values = np.array(dx_values)
    errors = np.array(errors)
    orders = np.log(errors[:-1] / errors[1:]) / np.log(dx_values[:-1] / dx_values[1:])

    return dx_values, errors, orders
