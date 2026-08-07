import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from homework1.exact_solution import phi_asymptotic, phi_transient, phi_exact
from homework1.diffusion import (
    phi_classical_diffusion,
    phi_asymptotic_diffusion,
    phi_diffusion_analytic,
    solve_diffusion_numerical,
    solve_diffusion_shooting,
    convergence_study,
)

logger = logging.getLogger(__name__)

N_COLS = 2

def create_figs_dir():
    figs_dir = os.path.join("docs", "homework1", "figs")
    os.makedirs(figs_dir, exist_ok=True)
    return figs_dir

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
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

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
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_diffusion_comparison(x, c_values, method='numerical', save_path=None):
    """
    Plots a comparison of exact transport, asymptotic transport, classical diffusion,
    and asymptotic diffusion for each c in c_values.
    """
    fig, axes = _make_grid(len(c_values))

    for i, c in enumerate(c_values):
        ax = axes[i]

        # Calculate solutions
        phi_as_val = phi_asymptotic(x, c, method=method)
        phi_tr_val = phi_transient(x, c)
        phi_ex_val = phi_as_val + phi_tr_val

        phi_class_val = phi_classical_diffusion(x, c)
        phi_as_diff_val = phi_asymptotic_diffusion(x, c, method=method)

        # Plot
        ax.plot(x, phi_ex_val, label=r'$\phi(x)$ (Exact Transport)', color='#2c3e50', linewidth=2.0)
        
        if c > 0.0:
            ax.plot(x, phi_as_val, label=r'$\phi_{as}(x)$ (Asymptotic)', color='#e74c3c', linestyle='--', linewidth=1.5)
            ax.plot(x, phi_as_diff_val, label=r'$\phi_{as,diff}(x)$ (Asymptotic Diff)', color='#2ecc71', linestyle='-.', linewidth=1.8)
        
        ax.plot(x, phi_class_val, label=r'$\phi_{class}(x)$ (Classical Diff)', color='#e67e22', linestyle=':', linewidth=1.8)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', shadow=True)

    _finalize_grid(fig, axes, len(c_values), 'Scalar Flux (log scale)')

    plt.suptitle(f'Comparison of Transport and Diffusion Solutions ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_diffusion_errors(x, c_values, method='numerical', save_path=None):
    """
    Plots the relative error (%) of classical and asymptotic diffusion approximations
    relative to the exact transport solution.
    """
    fig, axes = _make_grid(len(c_values))

    for i, c in enumerate(c_values):
        ax = axes[i]

        # Calculate solutions
        phi_as_val = phi_asymptotic(x, c, method=method)
        phi_tr_val = phi_transient(x, c)
        phi_ex_val = phi_as_val + phi_tr_val
        
        phi_class_val = phi_classical_diffusion(x, c)
        phi_as_diff_val = phi_asymptotic_diffusion(x, c, method=method)
        
        # Compute relative error (%)
        err_class = np.abs(phi_class_val - phi_ex_val) / phi_ex_val * 100.0
        
        ax.plot(x, err_class, label='Classical Diffusion', color='#e67e22', linewidth=2.0)
        
        if c > 0.0:
            err_as_diff = np.abs(phi_as_diff_val - phi_ex_val) / phi_ex_val * 100.0
            ax.plot(x, err_as_diff, label='Asymptotic Diffusion', color='#2ecc71', linewidth=2.0)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', shadow=True)

    _finalize_grid(fig, axes, len(c_values), 'Relative Error (%) (log scale)')

    plt.suptitle(f'Relative Error of Diffusion Approximations ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_numerical_vs_analytic_diffusion(c_values=[0.5, 0.7, 0.9], D=1.0/3.0, save_path_comp=None, save_path_err=None):
    """
    Compares the numerical diffusion solution against the analytical classical diffusion solution
    for a subset of c values on [-a, 0]. Generates:
    1. A plot of the analytic vs numerical solutions.
    2. A plot of the relative error between them.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # 1. Comparison Plot
    fig1, axes1 = plt.subplots(len(c_values), 1, figsize=(8, 12), sharex=False)
    if len(c_values) == 1:
        axes1 = [axes1]
        
    # 2. Error Plot
    fig2, axes2 = plt.subplots(len(c_values), 1, figsize=(8, 12), sharex=False)
    if len(c_values) == 1:
        axes2 = [axes2]
        
    for i, c in enumerate(c_values):
        # Solve numerically
        x_num, phi_num = solve_diffusion_numerical(c, D=D)
        
        # Compute analytical solution on the same grid
        phi_ana = phi_classical_diffusion(x_num, c, D=D)
        
        # Calculate relative error (%)
        rel_err = np.abs(phi_num - phi_ana) / phi_ana * 100.0
        
        # Plot Comparison
        ax1 = axes1[i]
        ax1.plot(x_num, phi_ana, label='Analytical Classical Diffusion', color='#3498db', linewidth=2.0)
        ax1.plot(x_num, phi_num, label='Numerical (Shooting Method)', color='#e74c3c', linestyle='--', linewidth=2.0)
        ax1.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Scalar Flux', fontsize=10)
        ax1.set_xlabel('$x$ [mean free paths]', fontsize=10)
        ax1.legend(loc='lower left', frameon=True)
        ax1.grid(True, which="both", ls="--", alpha=0.5)
        
        # Plot Error
        ax2 = axes2[i]
        ax2.plot(x_num, rel_err, color='#9b59b6', linewidth=2.0)
        ax2.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Relative Error (%)', fontsize=10)
        ax2.set_yscale('log')
        ax2.set_xlabel('$x$ [mean free paths]', fontsize=10)
        ax2.grid(True, which="both", ls="--", alpha=0.5)
        
    fig1.suptitle('Numerical vs. Analytical Diffusion Solutions', fontsize=14, fontweight='bold', y=0.98)
    fig1.tight_layout()
    if save_path_comp:
        fig1.savefig(save_path_comp, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path_comp}")
    plt.close(fig1)
    
    fig2.suptitle('Relative Error of Numerical Diffusion Solution', fontsize=14, fontweight='bold', y=0.98)
    fig2.tight_layout()
    if save_path_err:
        fig2.savefig(save_path_err, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path_err}")
    plt.close(fig2)

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

    if save_path:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
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

    if save_path:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
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

    if save_path:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close(fig)
