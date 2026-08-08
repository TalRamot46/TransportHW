import os
import logging
import numpy as np
import matplotlib

# These figures are written to file and never displayed, so select the
# non-interactive backend before pyplot is imported. This keeps the run
# headless and avoids starting a Tk event loop for nothing.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from homework1.exact_solution import phi_asymptotic, phi_transient, phi_exact
from homework1.diffusion import (
    phi_diffusion_analytic,
    solve_diffusion_shooting,
    convergence_study,
)
from homework1.criticality import (
    critical_dimensions,
    critical_dimensions_applied_bc,
    extrapolation_distance,
    METHOD_LABELS,
    MARSHAK_EXTRAPOLATION,
    MARK_EXTRAPOLATION,
    CASE_TABLE_23_C,
    CASE_TABLE_23_CZ0,
)

logger = logging.getLogger(__name__)

N_COLS = 2

def create_figs_dir():
    figs_dir = os.path.join("docs", "homework1", "figs")
    os.makedirs(figs_dir, exist_ok=True)
    return figs_dir

def _savefig(fig, save_path):
    """
    Writes `fig` to `save_path`, doing nothing if no path was given.

    The target is deleted first. Overwriting an existing PDF in place fails
    intermittently on this machine with "OSError: [Errno 22] Invalid argument",
    striking a different figure on each run, while a run into an empty directory
    always succeeds -- the signature of a scanner or indexer still holding a
    handle on the previous file. Writing to a name that does not exist yet
    sidesteps it.
    """
    if save_path is None:
        return

    if os.path.exists(save_path):
        os.remove(save_path)

    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    logger.info(f"Saved figure to: {save_path}")

def _make_grid(n_panels):
    """
    Creates a two-column grid of subplots with one panel per value of c.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    n_rows = -(-n_panels // N_COLS)  # ceiling division
    fig, axes = plt.subplots(n_rows, N_COLS, figsize=(12, 4 * n_rows), sharex=True)
    return fig, axes.flatten()

def _finalize_grid(fig, axes, n_panels, ylabel):
    """
    Removes the unused panels of the grid and applies the axis labels.

    The trailing panels are deleted first: with an odd number of c values the last
    row is only half filled. Because the grid is created with sharex=True, matplotlib
    hides the x tick labels on every panel that is not at the bottom of its column,
    so deleting a panel exposes a new bottom panel whose tick labels must be turned
    back on explicitly.
    """
    for j in range(n_panels, len(axes)):
        fig.delaxes(axes[j])

    for i in range(n_panels):
        if i % N_COLS == 0:
            axes[i].set_ylabel(ylabel, fontsize=10)
        if i >= n_panels - N_COLS:
            axes[i].tick_params(labelbottom=True)
            axes[i].set_xlabel('$x$ [mean free paths]', fontsize=10)

def plot_flux_components(x, c_values, method='numerical', save_path=None):
    """
    Plots the exact scalar flux, its asymptotic component, and its transient component
    for each c in c_values.
    """
    fig, axes = _make_grid(len(c_values))

    for i, c in enumerate(c_values):
        ax = axes[i]
        
        # Calculate components
        phi_as = phi_asymptotic(x, c, method=method)
        phi_tr = phi_transient(x, c)
        phi_ex = phi_as + phi_tr
        
        # Plot
        ax.plot(x, phi_ex, label=r'$\phi(x)$ (Exact)', color='#2c3e50', linewidth=2.0)
        
        if c > 0.0:
            ax.plot(x, phi_as, label=r'$\phi_{as}(x)$ (Asymptotic)', color='#e74c3c', linestyle='--', linewidth=1.8)
        
        ax.plot(x, phi_tr, label=r'$\phi_{tr}(x)$ (Transient)', color='#3498db', linestyle=':', linewidth=1.8)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', shadow=True)

    _finalize_grid(fig, axes, len(c_values), 'Scalar Flux (log scale)')

    plt.suptitle(f'Scalar Flux Components ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    _savefig(fig, save_path)
    plt.close(fig)

def plot_relative_contributions(x, c_values, method='numerical', save_path=None):
    """
    Plots the relative contribution of the transient and asymptotic components:
    phi_as/phi and phi_tr/phi as a function of x.
    """
    fig, axes = _make_grid(len(c_values))

    for i, c in enumerate(c_values):
        ax = axes[i]

        # Calculate components
        phi_as = phi_asymptotic(x, c, method=method)
        phi_tr = phi_transient(x, c)
        phi_ex = phi_as + phi_tr

        # Relative contributions
        # For c = 0, phi_as = 0, so contribution of phi_tr is 100%
        if c == 0.0:
            rel_as = np.zeros_like(x)
            rel_tr = np.ones_like(x)
        else:
            rel_as = phi_as / phi_ex
            rel_tr = phi_tr / phi_ex
            
        # Plot
        ax.plot(x, rel_as, label=r'$\phi_{as}(x) / \phi(x)$', color='#e74c3c', linewidth=2.0)
        ax.plot(x, rel_tr, label=r'$\phi_{tr}(x) / \phi(x)$', color='#3498db', linewidth=2.0)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        ax.legend(loc='center right', frameon=True, facecolor='white', edgecolor='none', shadow=True)

    _finalize_grid(fig, axes, len(c_values), 'Relative Contribution')

    plt.suptitle(f'Relative Contributions of Components ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    _savefig(fig, save_path)
    plt.close(fig)

DIFFUSION_COLORS = {'classical': '#e67e22', 'asymptotic': '#2ecc71'}
DIFFUSION_LABELS = {'classical': 'Classical Diff.', 'asymptotic': 'Asymptotic Diff.'}

def _diffusion_solutions(x, c, method='numerical'):
    """
    Returns {approximation: (phi_analytic, phi_numerical)} on the grid `x`, for
    every diffusion approximation applicable at this c.

    The numerical solver is evaluated directly on `x` rather than on its own
    grid, so no interpolation enters the comparison, and it is given the same
    `method` as the analytic form so that both use identical nu0.

    The asymptotic approximation is omitted at c = 0, where there is no discrete
    eigenvalue and hence no asymptotic mode.
    """
    approximations = ('classical',) if c == 0.0 else ('classical', 'asymptotic')

    solutions = {}
    for approximation in approximations:
        phi_ana = phi_diffusion_analytic(x, c, approximation, method=method)
        _, phi_num = solve_diffusion_shooting(c, approximation, method=method, x_eval=x)
        solutions[approximation] = (phi_ana, phi_num)
    return solutions

def plot_diffusion_comparison(x, c_values, method='numerical', save_path=None):
    """
    Plots the exact transport solution against the classical and asymptotic
    diffusion approximations, showing both the analytic closed form and the
    numerical solution of each, for every c in c_values.

    The analytic and numerical curves of a given approximation coincide to
    within about 2e-8%, so they overlie each other; that they do is the point of
    the figure. Their separation is resolved in the companion error figure.
    """
    fig, axes = _make_grid(len(c_values))

    for i, c in enumerate(c_values):
        ax = axes[i]

        # Exact transport solution
        phi_ex_val = phi_asymptotic(x, c, method=method) + phi_transient(x, c)
        ax.plot(x, phi_ex_val, label=r'$\phi(x)$ (Exact Transport)',
                color='#2c3e50', linewidth=2.4)

        for approximation, (phi_ana, phi_num) in _diffusion_solutions(x, c, method).items():
            color = DIFFUSION_COLORS[approximation]
            label = DIFFUSION_LABELS[approximation]
            ax.plot(x, phi_ana, label=f'{label} (analytic)',
                    color=color, linestyle='-', linewidth=1.6)
            ax.plot(x, phi_num, label=f'{label} (numerical)',
                    color=color, linestyle='--', linewidth=2.4, alpha=0.6)

        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)

        ax.legend(loc='upper right', fontsize=8, frameon=True, facecolor='white', edgecolor='none', shadow=True)

    _finalize_grid(fig, axes, len(c_values), 'Scalar Flux (log scale)')

    plt.suptitle(f'Exact Transport vs. Diffusion Approximations ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    _savefig(fig, save_path)
    plt.close(fig)

def plot_diffusion_errors(x, c_values, method='numerical', save_path=None):
    """
    Plots two relative errors per diffusion approximation:

    - numerical solution vs. the exact transport solution, which measures the
      error of the diffusion *approximation itself* and is the quantity
      Question 1d asks for;
    - numerical solution vs. the analytic diffusion solution, which measures the
      error of the *solver* and is the verification Question 2 asks for.

    The two differ by many orders of magnitude, which is the point: the solver
    error is negligible compared with the modelling error, so the upper curves
    can be read as properties of the approximation rather than of the code.
    """
    fig, axes = _make_grid(len(c_values))

    for i, c in enumerate(c_values):
        ax = axes[i]

        phi_ex_val = phi_asymptotic(x, c, method=method) + phi_transient(x, c)

        for approximation, (phi_ana, phi_num) in _diffusion_solutions(x, c, method).items():
            color = DIFFUSION_COLORS[approximation]
            label = DIFFUSION_LABELS[approximation]

            err_vs_exact = np.abs(phi_num - phi_ex_val) / phi_ex_val * 100.0
            err_vs_analytic = np.abs(phi_num - phi_ana) / phi_ana * 100.0

            ax.plot(x, err_vs_exact, label=f'{label}: numeric vs. exact',
                    color=color, linestyle='-', linewidth=2.0)
            ax.plot(x, err_vs_analytic, label=f'{label}: numeric vs. analytic',
                    color=color, linestyle=':', linewidth=1.8)

        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)

        ax.legend(loc='center right', fontsize=7.5, frameon=True, facecolor='white', edgecolor='none', shadow=True)

    _finalize_grid(fig, axes, len(c_values), 'Relative Error (%) (log scale)')

    plt.suptitle(f'Relative Error of Diffusion Approximations ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    _savefig(fig, save_path)
    plt.close(fig)

# ---------------------------------------------------------------------------
# Question 2: numerical diffusion solvers
# ---------------------------------------------------------------------------

APPROXIMATION_LABELS = {'classical': 'Classical', 'asymptotic': 'Asymptotic'}

SHOOTING_LABEL = 'Shooting (symmetry + radiation BC)'
SHOOTING_STYLE = dict(color='#e74c3c', linestyle='--', linewidth=1.8)

def plot_q2_solution_comparison(c_values, approximations=('classical', 'asymptotic'),
                                save_path=None):
    """
    Compares the numerical solution against the closed-form diffusion Green's
    function, one row per value of c and one column per approximation.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    n_rows, n_cols = len(c_values), len(approximations)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 3.5 * n_rows), squeeze=False)

    for i, c in enumerate(c_values):
        for j, approximation in enumerate(approximations):
            ax = axes[i][j]
            x, phi = solve_diffusion_shooting(c, approximation)

            ax.plot(x, phi_diffusion_analytic(x, c, approximation),
                    label='Analytic', color='#2c3e50', linewidth=2.5)
            ax.plot(x, phi, label=SHOOTING_LABEL, **SHOOTING_STYLE)

            ax.set_yscale('log')
            ax.set_title(f'{APPROXIMATION_LABELS[approximation]} diffusion, $c = {c}$',
                         fontsize=11, fontweight='bold')
            ax.grid(True, which='both', ls='--', alpha=0.5)
            if j == 0:
                ax.set_ylabel('Scalar Flux (log scale)', fontsize=10)
            if i == n_rows - 1:
                ax.set_xlabel('$x$ [mean free paths]', fontsize=10)
            if i == 0 and j == 0:
                ax.legend(loc='upper right', fontsize=8, frameon=True)

    fig.suptitle('Question 2: Numerical Diffusion Solvers vs. Analytic Green\'s Function',
                 fontsize=15, fontweight='bold', y=0.995)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

def plot_q2_error_profiles(c_values, approximations=('classical', 'asymptotic'),
                           save_path=None):
    """
    Relative error of the numerical solution against the closed form, as a
    function of x. With the radiation boundary condition the error stays flat
    across the domain instead of diverging at the outer edge.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    n_rows, n_cols = len(c_values), len(approximations)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 3.5 * n_rows), squeeze=False)

    for i, c in enumerate(c_values):
        for j, approximation in enumerate(approximations):
            ax = axes[i][j]

            x, phi = solve_diffusion_shooting(c, approximation)
            phi_ana = phi_diffusion_analytic(x, c, approximation)
            ax.plot(x, np.abs(phi - phi_ana) / phi_ana * 100.0,
                    label=SHOOTING_LABEL, **SHOOTING_STYLE)

            ax.set_yscale('log')
            ax.set_title(f'{APPROXIMATION_LABELS[approximation]} diffusion, $c = {c}$',
                         fontsize=11, fontweight='bold')
            ax.grid(True, which='both', ls='--', alpha=0.5)
            if j == 0:
                ax.set_ylabel('Relative Error (%)', fontsize=10)
            if i == n_rows - 1:
                ax.set_xlabel('$x$ [mean free paths]', fontsize=10)
            if i == 0 and j == 0:
                ax.legend(loc='lower right', fontsize=8, frameon=True)

    fig.suptitle('Question 2: Relative Error of the Numerical Solvers',
                 fontsize=15, fontweight='bold', y=0.995)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

# ---------------------------------------------------------------------------
# Question 3: critical dimensions of a multiplying medium
# ---------------------------------------------------------------------------

def plot_q3_critical_dimensions(c, save_path=None):
    """
    Critical planar half-thickness and critical sphere radius against c, from the
    asymptotic criticality relations, together with the error incurred by using
    the approximate formulas for the two inputs |nu0(c)| and z0(c).

    Both dimensions are plotted from the approximate inputs (the quantity the
    question asks for) and from the reference inputs -- the root of the
    transcendental equation for |nu0| and Case's Table 23 for z0 -- so that the
    cost of the approximation is visible rather than asserted.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    c = np.asarray(c, dtype=float)
    nu0_fit, z0_fit, half_fit, radius_fit = critical_dimensions(c, 'transport')
    nu0_ref, z0_ref, half_ref, radius_ref = critical_dimensions(c, 'transport-ref')
    _, _, half_alt, radius_alt = critical_dimensions(c, 'transport-q+')

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    panels = (
        (axes[0][0], half_fit, half_ref, half_alt, r'$a/2$',
         'Critical Half-Thickness (Planar)'),
        (axes[0][1], radius_fit, radius_ref, radius_alt, r'$\Sigma_t R_c$',
         'Critical Radius (Spherical)'),
    )
    for ax, fit, ref, alt, symbol, title in panels:
        ax.plot(c, ref, label=r'exact $|\nu_0|$, tabulated $z_0$',
                color='#2c3e50', linewidth=2.6)
        ax.plot(c, fit, label=r'fitted $|\nu_0|$, fitted $z_0$ ($q = -0.0199$)',
                color='#e74c3c', linestyle='--', linewidth=2.0)
        ax.plot(c, alt, label=r'fitted $|\nu_0|$, fitted $z_0$ ($q = +0.0199$)',
                color='#3498db', linestyle=':', linewidth=2.0)
        ax.set_yscale('log')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('$c$', fontsize=10)
        ax.set_ylabel(f'{symbol} [mean free paths]', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=8, frameon=True)

    # Error of the two inputs, against the transcendental root and Case's table.
    ax = axes[1][0]
    ax.plot(c, np.abs(nu0_fit - nu0_ref) / nu0_ref * 100.0,
            label=r'$|\nu_0|$: fit vs. transcendental root',
            color='#8e44ad', linewidth=2.0)
    ax.plot(c, np.abs(z0_fit - z0_ref) / z0_ref * 100.0,
            label=r'$z_0$: fit ($q = -0.0199$) vs. Table 23',
            color='#e74c3c', linestyle='--', linewidth=2.0)
    ax.plot(c, np.abs(extrapolation_distance(c, correction=0.0199) - z0_ref) / z0_ref * 100.0,
            label=r'$z_0$: fit ($q = +0.0199$) vs. Table 23',
            color='#3498db', linestyle=':', linewidth=2.0)
    ax.set_yscale('log')
    ax.set_title('Error of the Approximate Inputs', fontsize=12, fontweight='bold')
    ax.set_xlabel('$c$', fontsize=10)
    ax.set_ylabel('Relative Error (%)', fontsize=10)
    ax.grid(True, which='both', ls='--', alpha=0.5)
    ax.legend(loc='upper left', fontsize=8, frameon=True)

    # Error the approximate inputs induce in the two critical dimensions.
    ax = axes[1][1]
    for values, reference, label, color, style in (
            (half_fit, half_ref, r'$a/2$, $q = -0.0199$', '#e74c3c', '--'),
            (radius_fit, radius_ref, r'$\Sigma_t R_c$, $q = -0.0199$', '#c0392b', '-'),
            (half_alt, half_ref, r'$a/2$, $q = +0.0199$', '#3498db', ':'),
            (radius_alt, radius_ref, r'$\Sigma_t R_c$, $q = +0.0199$', '#2980b9', '-.')):
        ax.plot(c, np.abs(values - reference) / reference * 100.0,
                label=label, color=color, linestyle=style, linewidth=2.0)
    ax.set_yscale('log')
    ax.set_title('Error Induced in the Critical Dimensions', fontsize=12, fontweight='bold')
    ax.set_xlabel('$c$', fontsize=10)
    ax.set_ylabel('Relative Error (%)', fontsize=10)
    ax.grid(True, which='both', ls='--', alpha=0.5)
    ax.legend(loc='upper left', fontsize=8, frameon=True)

    fig.suptitle('Question 3: Critical Dimensions from the Asymptotic Transport Relations',
                 fontsize=15, fontweight='bold', y=0.995)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

Q3_METHOD_STYLES = {
    'transport': dict(color='#2c3e50', linestyle='-', linewidth=2.6),
    'marshak': dict(color='#e67e22', linestyle='-', linewidth=2.0),
    'mark': dict(color='#27ae60', linestyle='-', linewidth=2.0),
}

Q3_APPLIED_BC_STYLES = {
    'marshak': dict(color='#e67e22', linestyle=':', linewidth=1.6),
    'mark': dict(color='#27ae60', linestyle=':', linewidth=1.6),
}

def plot_q3_method_comparison(c, save_path=None):
    """
    The three criticality treatments of parts 3(a)-3(c) against each other:
    exact transport, diffusion with a Marshak boundary condition, and diffusion
    with a Mark boundary condition.

    The dotted curves are the same two diffusion problems with the boundary
    condition applied to the flux shape itself rather than replaced by an
    extrapolated zero, which is where the shortcut z0 = l0 stops being
    equivalent to solving the boundary condition.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    c = np.asarray(c, dtype=float)
    dimensions = {method: critical_dimensions(c, method)
                  for method in Q3_METHOD_STYLES}
    applied = {'marshak': critical_dimensions_applied_bc(c, MARSHAK_EXTRAPOLATION),
               'mark': critical_dimensions_applied_bc(c, MARK_EXTRAPOLATION)}

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # index 2 is a/2 and index 3 is Sigma_t R_c in the critical_dimensions tuple;
    # critical_dimensions_applied_bc returns the same pair on its own.
    for column, (index, applied_index, symbol, title, legend_loc) in enumerate((
            (2, 0, r'$a/2$', 'Critical Half-Thickness (Planar)', 'lower left'),
            (3, 1, r'$\Sigma_t R_c$', 'Critical Radius (Spherical)', 'upper left'))):

        ax = axes[0][column]
        for method, style in Q3_METHOD_STYLES.items():
            ax.plot(c, dimensions[method][index], label=METHOD_LABELS[method], **style)
        for method, style in Q3_APPLIED_BC_STYLES.items():
            ax.plot(c, applied[method][applied_index],
                    label=f'{METHOD_LABELS[method]}, BC applied', **style)
        ax.set_yscale('log')
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.set_xlabel('$c$', fontsize=10)
        ax.set_ylabel(f'{symbol} [mean free paths]', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=8, frameon=True)

        # Departure from the transport result, signed, so that over- and
        # under-estimation are distinguishable.
        reference = dimensions['transport'][index]
        ax = axes[1][column]
        for method in ('marshak', 'mark'):
            ax.plot(c, (dimensions[method][index] - reference) / reference * 100.0,
                    label=METHOD_LABELS[method], **Q3_METHOD_STYLES[method])
            ax.plot(c, (applied[method][applied_index] - reference) / reference * 100.0,
                    label=f'{METHOD_LABELS[method]}, BC applied',
                    **Q3_APPLIED_BC_STYLES[method])
        ax.axhline(0.0, color='#2c3e50', linewidth=1.2)
        ax.set_title(f'Departure from Exact Transport, {symbol}',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('$c$', fontsize=10)
        ax.set_ylabel('Relative Difference (\\%)', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc=legend_loc, fontsize=8, frameon=True)

    fig.suptitle('Question 3: Exact Transport vs. Diffusion with Marshak and Mark Conditions',
                 fontsize=15, fontweight='bold', y=0.995)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

def plot_q3_extrapolation_distance(save_path=None):
    """
    The extrapolation distance across the tabulated range, as c z0(c).

    This is the figure that settles the sign of the quadratic term: the
    tabulated product rises away from 0.710446 on *both* sides of c = 1, so the
    correction must be positive, and |q| = 0.0199 then reproduces the table to
    its printed precision near c = 1.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    # The expansion is about c = 1 and is not meant to hold as c -> 0, where the
    # tabulated c z0 climbs to 1; the window below keeps the near-critical
    # structure, which is a few parts in 10^4, legible.
    c = np.linspace(0.5, CASE_TABLE_23_C[-1], 600)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.plot(CASE_TABLE_23_C, CASE_TABLE_23_CZ0, 'o', color='#2c3e50',
            markersize=5, label='Case et al., Table 23')
    ax.plot(c, c * extrapolation_distance(c, correction=-0.0199),
            color='#e74c3c', linestyle='--', linewidth=2.0, label='fit, $q = -0.0199$')
    ax.plot(c, c * extrapolation_distance(c, correction=+0.0199),
            color='#3498db', linestyle=':', linewidth=2.0, label='fit, $q = +0.0199$')
    ax.axhline(0.710446, color='#7f8c8d', linewidth=1.0, alpha=0.8)

    ax.set_xlim(0.5, 3.0)
    ax.set_ylim(0.700, 0.726)
    ax.set_xlabel('$c$', fontsize=10)
    ax.set_ylabel('$c\\,z_0(c)$', fontsize=10)
    ax.set_title('Extrapolation Distance: Fit vs. Tabulated Values',
                 fontsize=13, fontweight='bold')
    ax.grid(True, which='both', ls='--', alpha=0.5)
    ax.legend(loc='upper right', fontsize=9, frameon=True)

    fig.tight_layout()
    _savefig(fig, save_path)
    plt.close(fig)

def plot_q2_convergence(c_values, approximations=('classical', 'asymptotic'),
                        save_path=None):
    """
    Tolerance-refinement study. A shooting method has no mesh spacing to refine:
    its accuracy is set by the integrator tolerance, so convergence is
    demonstrated by tightening rtol and showing the error follow it down to the
    floor set by round-off.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    fig, axes = plt.subplots(1, len(approximations),
                             figsize=(6 * len(approximations), 5), squeeze=False)

    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(c_values)))

    for j, approximation in enumerate(approximations):
        ax = axes[0][j]

        for c, color in zip(c_values, colors):
            rtols, errors = convergence_study(c, approximation)
            ax.loglog(rtols, errors, marker='o', color=color, alpha=0.9,
                      label=f'$c = {c}$')
            logger.info(f"  {approximation:11s} c={c}  "
                        f"error {errors[0]:.2e} (rtol {rtols[0]:.0e}) -> "
                        f"{errors[-1]:.2e} (rtol {rtols[-1]:.0e})")

        # Reference line: error tracking the requested tolerance one-for-one.
        ref = np.array([rtols[0], rtols[-1]])
        ax.loglog(ref, ref, 'k--', alpha=0.6, linewidth=1.2,
                  label='error = rtol')

        ax.set_title(f'{APPROXIMATION_LABELS[approximation]} diffusion',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('Integrator relative tolerance', fontsize=10)
        if j == 0:
            ax.set_ylabel('Relative $L_2$ error', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='upper left', fontsize=9, frameon=True)

    fig.suptitle('Question 2: Tolerance Convergence of the Shooting Solver',
                 fontsize=15, fontweight='bold', y=1.0)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)
