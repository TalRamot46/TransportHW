"""
Question 3(a)(b) figures: the exact planar flux against the two diffusion solutions.

One figure per scattering ratio c, one panel per time t.
"""

import numpy as np
from homework2.exact import phi_exact
from homework2.diffusion import phi_classical_diffusion, phi_asymptotic_diffusion
from homework2.solver import solve_diffusion
from homework2.figures import make_grid, label_grid, savefig, close

C_VALUES = (0.6, 0.8, 1.0, 1.2, 1.5)
TIMES = (1.0, 2.0, 3.0, 4.0, 7.0, 15.0)

# Fraction of the causal front covered by the x grid, so the diffusion tails that
# leak past |x| = vt stay visible.
X_REACH = 1.25
X_POINTS = 1500

# Vertical padding, as a factor below the smallest and above the largest value the
# exact curve reaches inside the front, plus a floor on the total range so that an
# early-time panel -- where the exact solution barely varies -- is not blown up.
Y_PAD_BELOW, Y_PAD_ABOVE, Y_MIN_DECADES = 0.25, 3.0, 2.0

CURVES = (
    ('Exact (Paasschens)', '#2c3e50', '-', 2.2),
    ('Classical diffusion', '#e67e22', '--', 1.8),
    ('Asymptotic diffusion', '#27ae60', ':', 2.0),
)

# Numerical solutions are drawn as markers over the analytic curve of the same colour,
# so that agreement reads as markers sitting on a line.
NUMERICAL_COLORS = {'classical': '#e67e22', 'asymptotic': '#27ae60'}
MARKERS_PER_PANEL = 14

def _panel_curves(x, t, c):
    """The three fluxes on the grid `x`; the exact one is masked beyond the causal front."""
    exact = phi_exact(x, t, c)
    return (np.where(x < t, exact, np.nan),
            phi_classical_diffusion(x, t, c),
            phi_asymptotic_diffusion(x, t, c))

def _y_limits(exact):
    """Log y-limits framing the exact curve, widened downwards to at least Y_MIN_DECADES."""
    top = np.nanmax(exact) * Y_PAD_ABOVE
    bottom = np.nanmin(exact) * Y_PAD_BELOW
    return min(bottom, top * 10.0 ** -Y_MIN_DECADES), top

def _draw_panel(ax, t, c):
    """Draws one (c, t) panel: the three curves, the causal front, and the axis scaling."""
    x = np.linspace(0.0, X_REACH * t, X_POINTS)
    curves = _panel_curves(x, t, c)

    for values, (label, color, style, width) in zip(curves, CURVES):
        ax.plot(x, values, label=label, color=color, linestyle=style, linewidth=width)

    ax.axvline(t, color='#7f8c8d', linewidth=1.0, alpha=0.8)
    ax.set_yscale('log')
    ax.set_ylim(*_y_limits(curves[0]))
    ax.set_xlim(0.0, X_REACH * t)
    ax.set_title(f'$t = {t:g}$', fontsize=11, fontweight='bold')
    ax.grid(True, which='both', ls='--', alpha=0.5)

def _markers(ax, x, y, color, label):
    """Draws a curve as sparse open markers, thinned to MARKERS_PER_PANEL points."""
    every = max(1, len(x) // MARKERS_PER_PANEL)
    ax.plot(x[::every], y[::every], linestyle='none', marker='o', markersize=4.0,
            markerfacecolor='none', markeredgewidth=1.1, color=color, label=label)

def _draw_numerical(ax, x, phi, t, color, label):
    """Overlays one numerical solution as markers, over the panel's x range."""
    inside = x <= X_REACH * t
    _markers(ax, x[inside], phi[inside], color, label)

def plot_numerical_for_c(c, times=TIMES, save_path=None):
    """
    Question 3(c): the numerical diffusion solutions over the analytic ones and the exact flux.

    The solver runs once per approximation, out to the largest time, and every panel reads
    off the same run.
    """
    fig, axes = make_grid(len(times))
    runs = {approximation: solve_diffusion(c, approximation, times, X_REACH * max(times))
            for approximation in NUMERICAL_COLORS}

    for i, (ax, t) in enumerate(zip(axes, times)):
        _draw_panel(ax, t, c)
        for approximation, (x, solution) in runs.items():
            label = f'{approximation.capitalize()} (numerical)' if i == 0 else None
            _draw_numerical(ax, x, solution[t], t, NUMERICAL_COLORS[approximation], label)

    axes[0].legend(loc='lower left', fontsize=7.5, frameon=True, facecolor='white')
    label_grid(axes, len(times), '$x$ [mean free paths]', r'$\phi(x,t)$ (log scale)')

    fig.suptitle(f'Numerical vs. analytic diffusion and exact transport, $c = {c:g}$',
                 fontsize=14, fontweight='bold', y=0.995)
    fig.tight_layout()

    savefig(fig, save_path)
    close(fig)

ANALYTIC_DIFFUSION = {'classical': phi_classical_diffusion, 'asymptotic': phi_asymptotic_diffusion}

# The error is evaluated short of the causal front, where the exact solution drops to zero
# and the relative error diverges for reasons that say nothing about the approximation.
ERROR_REACH = 0.99
ERROR_FLOOR = 1e-2      # percent; below this the approximation is exact for any practical purpose

def _relative_error(phi, exact):
    """Relative error in percent."""
    return np.abs(phi - exact) / exact * 100.0

def plot_errors_for_c(c, times=TIMES, save_path=None):
    """Question 3(d): relative error of each diffusion solution against the exact transport flux."""
    fig, axes = make_grid(len(times))
    runs = {approximation: solve_diffusion(c, approximation, times, X_REACH * max(times))
            for approximation in NUMERICAL_COLORS}

    for i, (ax, t) in enumerate(zip(axes, times)):
        x = np.linspace(0.0, ERROR_REACH * t, X_POINTS)
        for approximation, color in NUMERICAL_COLORS.items():
            analytic = ANALYTIC_DIFFUSION[approximation]
            ax.plot(x, _relative_error(analytic(x, t, c), phi_exact(x, t, c)), color=color,
                    linewidth=2.0, label=f'{approximation.capitalize()}' if i == 0 else None)

            x_num, solution = runs[approximation]
            inside = x_num <= ERROR_REACH * t
            _markers(ax, x_num[inside], _relative_error(solution[t][inside],
                                                        phi_exact(x_num[inside], t, c)),
                     color, f'{approximation.capitalize()} (numerical)' if i == 0 else None)

        ax.set_yscale('log')
        # A common floor across panels; the downward spikes are sign changes, where a
        # diffusion curve crosses the exact one, and carry no information below this.
        ax.set_ylim(bottom=ERROR_FLOOR)
        ax.set_xlim(0.0, ERROR_REACH * t)
        ax.set_title(f'$t = {t:g}$', fontsize=11, fontweight='bold')
        ax.grid(True, which='both', ls='--', alpha=0.5)

    axes[0].legend(loc='lower left', fontsize=7.5, frameon=True, facecolor='white')
    label_grid(axes, len(times), '$x$ [mean free paths]', 'Relative error (\\%)')

    fig.suptitle(f'Relative error of the diffusion approximations, $c = {c:g}$',
                 fontsize=14, fontweight='bold', y=0.995)
    fig.tight_layout()

    savefig(fig, save_path)
    close(fig)

def plot_comparison_for_c(c, times=TIMES, save_path=None):
    """Builds the exact-vs-diffusion comparison figure for one c, one panel per time."""
    fig, axes = make_grid(len(times))

    for ax, t in zip(axes, times):
        _draw_panel(ax, t, c)

    axes[0].legend(loc='lower left', fontsize=8, frameon=True, facecolor='white')
    label_grid(axes, len(times), '$x$ [mean free paths]', r'$\phi(x,t)$ (log scale)')

    fig.suptitle(f'Exact transport vs. diffusion, $c = {c:g}$ '
                 r'(vertical line: causal front $|x| = vt$)',
                 fontsize=14, fontweight='bold', y=0.995)
    fig.tight_layout()

    savefig(fig, save_path)
    close(fig)
