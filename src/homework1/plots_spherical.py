import numpy as np
import matplotlib

# Same reasoning as in plots.py: these figures are written to file, never shown.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from homework1.plots import _savefig, APPROXIMATION_LABELS
from homework1.spherical import (
    build_medium,
    buckling,
    critical_radius,
    analytic_critical_radius,
    k_eigenvalue,
    mesh_convergence,
)
from homework1.materials import BENCHMARK, FISSILE, critical_mass

# Figures for Questions 4 and 5. Kept out of plots.py, which is already long.

CURVE_COLOURS = ('#2c3e50', '#e74c3c', '#3498db')

def _use_style():
    """Applies the same grid style as the Question 1-3 figures."""
    plt.style.use('seaborn-v0_8-whitegrid'
                  if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

def _critical_radii(c_values, approximation, n_cells):
    """Numerical and analytic critical radii over a list of c, in mean free paths."""
    media = [build_medium(1.0, 1.0, c, approximation) for c in c_values]
    numerical = np.array([critical_radius(m, n_cells=n_cells) for m in media])
    analytic = np.array([analytic_critical_radius(m) for m in media])
    return numerical, analytic

def plot_q4_critical_radius(c_values, n_cells=400, save_path=None):
    """
    Critical radius from the k = 1 search against the analytic radius pi/B - z0,
    for both diffusion approximations, with the relative difference below.
    """
    _use_style()
    c_values = np.asarray(c_values, dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    for j, approximation in enumerate(('classical', 'asymptotic')):
        numerical, analytic = _critical_radii(c_values, approximation, n_cells)

        ax = axes[0][j]
        ax.plot(c_values, analytic, label=r'analytic, $\pi / B - z_0$',
                color='#2c3e50', linewidth=2.6)
        ax.plot(c_values, numerical, label=f'numerical, $k = 1$ ({n_cells} cells)',
                color='#e74c3c', linestyle='none', marker='o', markersize=5,
                markerfacecolor='none', markeredgewidth=1.4)
        ax.set_yscale('log')
        ax.set_title(f'{APPROXIMATION_LABELS[approximation]} diffusion',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('$c$', fontsize=10)
        ax.set_ylabel(r'$\Sigma_t R_c$ [mean free paths]', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='upper right', fontsize=9, frameon=True)

        ax = axes[1][j]
        ax.plot(c_values, np.abs(numerical - analytic) / analytic * 100.0,
                color='#8e44ad', linewidth=2.0, marker='o', markersize=4)
        ax.set_yscale('log')
        ax.set_xlabel('$c$', fontsize=10)
        ax.set_ylabel('Relative difference (%)', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)

    fig.suptitle('Question 4: Critical Radius from the $k = 1$ Search vs. Analytic',
                 fontsize=15, fontweight='bold', y=0.995)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

def plot_q4_mesh_convergence(c_values=(1.05, 1.5, 2.0), save_path=None):
    """
    Relative error of the numerical critical radius against mesh refinement.
    The finite-volume scheme is second order, so the error falls as 1/N^2.
    """
    _use_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), squeeze=False)

    for j, approximation in enumerate(('classical', 'asymptotic')):
        ax = axes[0][j]
        for c, colour in zip(c_values, CURVE_COLOURS):
            cells, _, errors = mesh_convergence(build_medium(1.0, 1.0, c, approximation))
            ax.loglog(cells, errors, marker='o', color=colour, linewidth=1.8,
                      label=f'$c = {c}$')

        # Second-order guide, anchored on the coarsest point of the last curve.
        ax.loglog(cells, errors[0] * (cells[0] / cells)**2, 'k--', alpha=0.6,
                  linewidth=1.2, label=r'$\propto N^{-2}$')

        ax.set_title(f'{APPROXIMATION_LABELS[approximation]} diffusion',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('Number of cells $N$', fontsize=10)
        if j == 0:
            ax.set_ylabel('Relative error of $R_c$', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='lower left', fontsize=9, frameon=True)

    fig.suptitle('Question 4: Mesh Convergence of the Critical Radius',
                 fontsize=15, fontweight='bold', y=1.0)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

def plot_q4_flux_profiles(c_values=(1.05, 1.5, 2.0), n_cells=400, save_path=None):
    """
    Converged critical flux shape against the analytic fundamental mode
    sin(Br)/(Br), both normalised to unity at the centre.
    """
    _use_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), squeeze=False)

    for j, approximation in enumerate(('classical', 'asymptotic')):
        ax = axes[0][j]
        for c, colour in zip(c_values, CURVE_COLOURS):
            medium = build_medium(1.0, 1.0, c, approximation)
            R = critical_radius(medium, n_cells=n_cells)
            _, r, phi = k_eigenvalue(R, medium, n_cells=n_cells)

            # np.sinc(y) = sin(pi y) / (pi y), which supplies the r = 0 limit.
            shape = np.sinc(buckling(medium) * r / np.pi)
            ax.plot(r / R, phi / phi[0], color=colour, linewidth=2.2,
                    label=f'$c = {c}$, numerical')
            ax.plot(r / R, shape / shape[0], color=colour, linestyle='--',
                    linewidth=1.4, alpha=0.8, label=f'$c = {c}$, $\\sin(Br)/Br$')

        ax.axvline(1.0, color='grey', linestyle=':', linewidth=1.2)
        ax.set_title(f'{APPROXIMATION_LABELS[approximation]} diffusion',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('$r / R_c$   (the mesh runs out to $R_c + z_0$)', fontsize=10)
        if j == 0:
            ax.set_ylabel(r'$\phi(r) / \phi(0)$', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='lower left', fontsize=7, frameon=True, ncol=2)

    fig.suptitle('Question 4: Critical Flux Shape vs. the Fundamental Mode',
                 fontsize=15, fontweight='bold', y=1.0)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)

def plot_q5_criticality(n_cells=400, save_path=None):
    """
    k against sphere radius for the two fissile benchmark materials, with the
    k = 1 crossing that fixes the critical radius and hence the critical mass.
    """
    _use_style()
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), squeeze=False)

    for j, name in enumerate(FISSILE):
        material = BENCHMARK[name]
        ax = axes[0][j]

        for approximation, colour in (('classical', '#2c3e50'),
                                      ('asymptotic', '#e74c3c')):
            medium = build_medium(material.sigma_t, material.sigma_a,
                                  material.nu_sigma_f, approximation)
            R_c = critical_radius(medium, n_cells=n_cells)
            radii = np.linspace(0.4 * R_c, 1.8 * R_c, 40)
            k = [k_eigenvalue(R, medium, n_cells=n_cells)[0] for R in radii]

            ax.plot(radii, k, color=colour, linewidth=2.2,
                    label=APPROXIMATION_LABELS[approximation])
            ax.plot([R_c], [1.0], color=colour, marker='o', markersize=7,
                    markerfacecolor='none', markeredgewidth=1.8)
            ax.annotate(f'$R_c = {R_c:.3f}$ cm,  '
                        f'$M_c = {critical_mass(R_c, material.density):.2f}$ kg',
                        xy=(0.03, 0.90 if approximation == 'classical' else 0.82),
                        xycoords='axes fraction', fontsize=9, color=colour)

        ax.axhline(1.0, color='grey', linestyle=':', linewidth=1.4)
        ax.set_title(f'{name}   ($c = {material.c:.2f}$, '
                     rf'$\rho = {material.density}$ g/cm$^3$)',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('Sphere radius $R$ [cm]', fontsize=10)
        if j == 0:
            ax.set_ylabel('Multiplication factor $k$', fontsize=10)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        ax.legend(loc='lower right', fontsize=9, frameon=True)

    fig.suptitle('Question 5: Criticality of the Bare Benchmark Spheres',
                 fontsize=15, fontweight='bold', y=1.0)
    fig.tight_layout()

    _savefig(fig, save_path)
    plt.close(fig)
