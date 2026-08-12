"""Question 1: the exact flux, its two components, and the diffusion approximations."""

import os
import numpy as np
from homework1.exact_solution import (compute_nu0_numerical, compute_nu0_approx,
                                      phi_asymptotic, phi_transient)
from homework1.diffusion import phi_diffusion_analytic, solve_diffusion_shooting
from homework1.figures import make_grid, panel, finish, savefig
from homework1.tables import log_section, log_table

C_VALUES = (0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
X = np.linspace(1e-3, 5.0, 500)
XLABEL = '$x$ [mean free paths]'

DIFFUSION_COLORS = {'classical': '#e67e22', 'asymptotic': '#2ecc71'}
DIFFUSION_LABELS = {'classical': 'Classical diff.', 'asymptotic': 'Asymptotic diff.'}

def _components(c, method):
    """(asymptotic, transient, exact) scalar flux on the X grid."""
    phi_as = phi_asymptotic(X, c, method=method)
    phi_tr = phi_transient(X, c)
    return phi_as, phi_tr, phi_as + phi_tr

def _diffusion_solutions(c, method):
    """{approximation: (analytic, numerical)} on X, omitting the asymptotic one at c = 0.

    The solver is evaluated on X directly and with the same nu0 as the closed form,
    so neither interpolation nor a mismatched eigenvalue enters the comparison.
    """
    approximations = ('classical',) if c == 0.0 else ('classical', 'asymptotic')
    return {a: (phi_diffusion_analytic(X, c, a, method=method),
                solve_diffusion_shooting(c, a, method=method, x_eval=X)[1])
            for a in approximations}

def plot_flux_components(method, save_path):
    """Exact flux with its asymptotic and transient components, one panel per c."""
    fig, axes = make_grid(len(C_VALUES))

    for ax, c in zip(axes, C_VALUES):
        phi_as, phi_tr, phi = _components(c, method)
        ax.plot(X, phi, label=r'$\phi$ (exact)', color='#2c3e50', linewidth=2.0)
        if c > 0.0:
            ax.plot(X, phi_as, label=r'$\phi_{as}$', color='#e74c3c',
                    linestyle='--', linewidth=1.8)
        ax.plot(X, phi_tr, label=r'$\phi_{tr}$', color='#3498db',
                linestyle=':', linewidth=1.8)
        panel(ax, f'$c = {c}$', XLABEL, 'Scalar flux', log=True, legend='upper right')

    finish(fig, f'Scalar flux components ({method})')
    savefig(fig, save_path)

def plot_relative_contributions(method, save_path):
    """Fractions phi_as/phi and phi_tr/phi; at c = 0 the transient carries everything."""
    fig, axes = make_grid(len(C_VALUES))

    for ax, c in zip(axes, C_VALUES):
        phi_as, phi_tr, phi = _components(c, method)
        ax.plot(X, phi_as / phi, label=r'$\phi_{as} / \phi$', color='#e74c3c', linewidth=2.0)
        ax.plot(X, phi_tr / phi, label=r'$\phi_{tr} / \phi$', color='#3498db', linewidth=2.0)
        ax.set_ylim(-0.05, 1.05)
        panel(ax, f'$c = {c}$', XLABEL, 'Relative contribution', legend='center right')

    finish(fig, f'Relative contributions of the components ({method})')
    savefig(fig, save_path)

def plot_diffusion_comparison(method, save_path):
    """Exact transport against both diffusion approximations, analytic and numerical."""
    fig, axes = make_grid(len(C_VALUES))

    for ax, c in zip(axes, C_VALUES):
        ax.plot(X, _components(c, method)[2], label=r'$\phi$ (exact transport)',
                color='#2c3e50', linewidth=2.4)
        for approximation, (analytic, numerical) in _diffusion_solutions(c, method).items():
            color, label = DIFFUSION_COLORS[approximation], DIFFUSION_LABELS[approximation]
            ax.plot(X, analytic, label=f'{label} (analytic)', color=color, linewidth=1.6)
            ax.plot(X, numerical, label=f'{label} (numerical)', color=color,
                    linestyle='--', linewidth=2.4, alpha=0.6)
        panel(ax, f'$c = {c}$', XLABEL, 'Scalar flux', log=True, legend='upper right')

    finish(fig, f'Exact transport vs. diffusion ({method})')
    savefig(fig, save_path)

def plot_diffusion_errors(method, save_path):
    """Error of each approximation against exact transport (1d) and against its own
    closed form (the solver check of Question 2). The two differ by orders of magnitude."""
    fig, axes = make_grid(len(C_VALUES))

    for ax, c in zip(axes, C_VALUES):
        exact = _components(c, method)[2]
        for approximation, (analytic, numerical) in _diffusion_solutions(c, method).items():
            color, label = DIFFUSION_COLORS[approximation], DIFFUSION_LABELS[approximation]
            ax.plot(X, np.abs(numerical - exact) / exact * 100.0,
                    label=f'{label}: vs. exact', color=color, linewidth=2.0)
            ax.plot(X, np.abs(numerical - analytic) / analytic * 100.0,
                    label=f'{label}: vs. analytic', color=color,
                    linestyle=':', linewidth=1.8)
        panel(ax, f'$c = {c}$', XLABEL, 'Relative error (%)', log=True,
              legend='center right', fontsize=7.5)

    finish(fig, f'Relative error of the diffusion approximations ({method})')
    savefig(fig, save_path)

def report(figs):
    """Prints the nu0 comparison table and writes the Question 1 figures, both methods."""
    log_section('Homework 1 Question 1', 'Discrete eigenvalue nu0 of Case (1953).')
    rows = []
    for c in C_VALUES[1:]:
        numerical, approximate = compute_nu0_numerical(c), compute_nu0_approx(c)
        rows.append([f'{c}', f'{numerical:.10f}', f'{approximate:.10f}',
                     f'{abs(numerical - approximate) / numerical * 100.0:.8f}'])
    log_table(['c', 'nu0 numerical', 'nu0 approx', 'rel. error %'], rows)

    for method in ('numerical', 'approx'):
        for plot, name in ((plot_flux_components, 'flux_components'),
                           (plot_relative_contributions, 'relative_contributions'),
                           (plot_diffusion_comparison, 'diffusion_comparison'),
                           (plot_diffusion_errors, 'diffusion_errors')):
            plot(method, os.path.join(figs, f'{name}_{method}.pdf'))
