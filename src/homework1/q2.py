"""Question 2: the planar diffusion solver against its closed form."""

import os
import numpy as np
from homework1.diffusion import (phi_diffusion_analytic, solve_diffusion_shooting,
                                 absorption_balance, convergence_study)
from homework1.figures import subplots, panel, finish, savefig
from homework1.tables import log_section, log_table

C_VALUES = (0.5, 0.7, 0.9)
APPROXIMATIONS = ('classical', 'asymptotic')
LABELS = {'classical': 'Classical', 'asymptotic': 'Asymptotic'}
XLABEL = '$x$ [mean free paths]'
NUMERICAL = dict(color='#e74c3c', linestyle='--', linewidth=1.8)

def plot_solution_comparison(save_path):
    """Numerical solution against the analytic Green's function, one panel per case."""
    fig, axes = subplots(len(C_VALUES), len(APPROXIMATIONS), height=3.5)

    for i, c in enumerate(C_VALUES):
        for j, approximation in enumerate(APPROXIMATIONS):
            x, phi = solve_diffusion_shooting(c, approximation)
            axes[i][j].plot(x, phi_diffusion_analytic(x, c, approximation),
                            label='Analytic', color='#2c3e50', linewidth=2.5)
            axes[i][j].plot(x, phi, label='Shooting solver', **NUMERICAL)
            panel(axes[i][j], f'{LABELS[approximation]} diffusion, $c = {c}$',
                  XLABEL, 'Scalar flux', log=True, legend='upper right')

    finish(fig, "Question 2: numerical solver vs. analytic Green's function")
    savefig(fig, save_path)

def plot_error_profiles(save_path):
    """Relative error across the domain; the radiation condition keeps it flat."""
    fig, axes = subplots(len(C_VALUES), len(APPROXIMATIONS), height=3.5)

    for i, c in enumerate(C_VALUES):
        for j, approximation in enumerate(APPROXIMATIONS):
            x, phi = solve_diffusion_shooting(c, approximation)
            analytic = phi_diffusion_analytic(x, c, approximation)
            axes[i][j].plot(x, np.abs(phi - analytic) / analytic * 100.0, **NUMERICAL)
            panel(axes[i][j], f'{LABELS[approximation]} diffusion, $c = {c}$',
                  XLABEL, 'Relative error (%)', log=True)

    finish(fig, 'Question 2: relative error of the numerical solver')
    savefig(fig, save_path)

def plot_convergence(save_path):
    """Error against integrator tolerance: a shooting method has no mesh to refine."""
    fig, axes = subplots(1, len(APPROXIMATIONS), height=5.0)
    colors = ['#2c3e50', '#e74c3c', '#3498db']

    for j, approximation in enumerate(APPROXIMATIONS):
        for c, color in zip(C_VALUES, colors):
            rtols, errors = convergence_study(c, approximation)
            axes[0][j].loglog(rtols, errors, marker='o', color=color, label=f'$c = {c}$')
        axes[0][j].loglog(rtols, rtols, 'k--', alpha=0.6, label='error = rtol')
        panel(axes[0][j], f'{LABELS[approximation]} diffusion',
              'Integrator relative tolerance', 'Relative $L_2$ error', legend='upper left')

    finish(fig, 'Question 2: tolerance convergence of the shooting solver')
    savefig(fig, save_path)

def report(figs):
    """Prints the neutron balance and the tolerance study, and writes the figures."""
    log_section('Homework 1 Question 2',
                'Neutron balance, Sigma_a * integral(phi dx), which must be 1.')
    rows = []
    for approximation in APPROXIMATIONS:
        for c in C_VALUES:
            x, phi = solve_diffusion_shooting(c, approximation)
            rows.append([f'{c}', approximation,
                         f'{absorption_balance(x, phi, c, approximation):.10f}'])
    log_table(['c', 'approximation', 'balance'], rows)

    rows = []
    for approximation in APPROXIMATIONS:
        for c in C_VALUES:
            rtols, errors = convergence_study(c, approximation)
            rows.append([f'{c}', approximation, f'{errors[0]:.2e}', f'{errors[-1]:.2e}'])
    log_table(['c', 'approximation', 'L2 error at rtol 1e-4', 'at rtol 1e-11'], rows)

    plot_solution_comparison(os.path.join(figs, 'q2_solution_comparison.pdf'))
    plot_error_profiles(os.path.join(figs, 'q2_error_profiles.pdf'))
    plot_convergence(os.path.join(figs, 'q2_convergence.pdf'))
