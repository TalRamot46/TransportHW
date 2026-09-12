"""
Question 3 figures: the exact planar flux against diffusion (3a, 3b), the explicit solver
against the closed form it discretises (3c), and the solver's error against exact transport (3d).

One panel per case throughout; the case each panel shows is annotated inside it.
"""

import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, ScalarFormatter
from homework2.exact import phi_exact
from homework2.diffusion import (phi_classical_diffusion, phi_asymptotic_diffusion,
                                 diffusion_coefficient)
from homework2.solver import solve, bulk, X_REACH
from homework2.figures import (make_grid, subplots, panel, case_label, label_grid,
                               finish, savefig, close, LEGEND_SIZE,
                               NAVY, ORANGE, GREEN, RED, GREY)

C_VALUES = (0.6, 0.8, 1.0, 1.2, 1.5)
TIMES = (1.0, 2.0, 3.0, 4.0, 7.0, 15.0)

X_POINTS = 1500
SOLVER_C = 1.0                       # the 3(c) figures fix c: the scheme does not depend on it
NODE_COUNTS = (250, 500, 1000, 2000)
MARKERS_PER_PANEL = 12               # thinning, so the solver reads as markers over a line

ERROR_TIMES = (1.0, 4.0, 15.0)
ERROR_STYLES = (':', '--', '-')      # earliest time dotted, latest solid
FRONT_FRACTION = 0.98                # the exact flux vanishes at |x| = vt; stop just inside it
# The classical line is drawn wider so that at c = 1, where the two D coincide, it still
# shows from under the asymptotic one instead of vanishing.
APPROXIMATIONS = (('Classical', 'classical', ORANGE, 2.8), ('Asymptotic', 'asymptotic', GREEN, 1.6))
ERROR_FLOOR = 1e-3                   # below this the panel is all sign-change dip, not error

# Vertical padding, as a factor below the smallest and above the largest value the exact
# curve reaches inside the front, plus a floor on the total range so that an early-time
# panel -- where the exact solution barely varies -- is not blown up.
Y_PAD_BELOW, Y_PAD_ABOVE, Y_MIN_DECADES = 0.25, 3.0, 2.0

CURVES = (
    ('Exact (Paasschens)', NAVY, '-', 2.2),
    ('Classical diffusion', ORANGE, '--', 1.8),
    ('Asymptotic diffusion', GREEN, ':', 2.0),
)

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

    ax.axvline(t, color=GREY, linewidth=1.0, alpha=0.8)
    ax.set_ylim(*_y_limits(curves[0]))
    ax.set_xlim(0.0, X_REACH * t)
    case_label(ax, f'$t = {t:g}$')
    panel(ax, log=True)

def plot_comparison_for_c(c, times=TIMES, save_path=None):
    """Parts 3(a) and 3(b): exact transport against both diffusion curves, one panel per time."""
    fig, axes = make_grid(len(times))

    for ax, t in zip(axes, times):
        _draw_panel(ax, t, c)

    axes[0].legend(loc='lower left', fontsize=9, frameon=True)
    label_grid(axes, len(times), '$x$ [mean free paths]', r'$\phi(x,t)$')

    finish(fig)
    savefig(fig, save_path)
    close(fig)

def plot_solver(c=SOLVER_C, times=TIMES, save_path=None):
    """Part 3(c): the explicit solver over the closed form it discretises, one panel per time."""
    x, fluxes = solve(c, 'classical', times=times)
    fig, axes = make_grid(len(times))

    for ax, t in zip(axes, times):
        window = x <= X_REACH * t
        # The window holds a different number of nodes in every panel, so the thinning is
        # set per panel rather than as one stride, or t = 1 would carry two markers.
        every = max(1, int(window.sum()) // MARKERS_PER_PANEL)
        ax.plot(x[window], phi_classical_diffusion(x[window], t, c),
                label='Closed form', color=NAVY, linewidth=2.2)
        ax.plot(x[window][::every], fluxes[t][window][::every],
                label='Explicit solver', color=RED, linestyle='none', marker='o',
                markersize=4.5, markerfacecolor='none', markeredgewidth=1.4)
        ax.set_xlim(0.0, X_REACH * t)
        panel(ax, log=True)
        # After the scale is log, not before: on the linear axis the limits mean
        # something else entirely and the headroom swallows the whole panel.
        ax.set_ylim(top=ax.get_ylim()[1] * 4.0)
        case_label(ax, f'$t = {t:g}$')

    axes[0].legend(loc='lower left', fontsize=9, frameon=True)
    label_grid(axes, len(times), '$x$ [mean free paths]', r'$\phi(x,t)$')

    finish(fig)
    savefig(fig, save_path)
    close(fig)

def _error_profile(ax, c, times):
    """Relative error against the closed form, over x, one curve per time."""
    x, fluxes = solve(c, 'classical', times=times)
    D = diffusion_coefficient(c, 'classical')

    for color, t in zip((NAVY, ORANGE, GREEN, RED), times):
        inside = bulk(x, t, D)
        exact = phi_classical_diffusion(x[inside], t, c)
        ax.plot(x[inside], np.abs(fluxes[t][inside] / exact - 1.0),
                color=color, linewidth=2.0, label=f'$t = {t:g}$')

    panel(ax, '$x$ [mean free paths]', 'Relative error', log=True, legend='lower right')

def _order_panel(ax, c, t):
    """Peak relative error against the node count, with the second-order guide."""
    errors = []
    for n_nodes in NODE_COUNTS:
        x, fluxes = solve(c, 'classical', times=(t,), n_nodes=n_nodes)
        inside = bulk(x, t, diffusion_coefficient(c, 'classical'))
        errors.append(np.abs(fluxes[t][inside] / phi_classical_diffusion(x[inside], t, c)
                             - 1.0).max())

    nodes = np.array(NODE_COUNTS, dtype=float)
    ax.loglog(nodes, errors, color=RED, linewidth=2.0, marker='o', markersize=5,
              label=f'measured, $t = {t:g}$')
    ax.loglog(nodes, errors[0] * (nodes[0] / nodes) ** 2, color=GREY, linestyle='--',
              linewidth=1.4, label=r'$\propto N^{-2}$')

    # Under two decades matplotlib labels the minor ticks too, and they collide; the node
    # counts themselves are the informative ticks here.
    ax.set_xticks(nodes)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    panel(ax, 'Number of nodes $N$', 'Peak relative error', legend='upper right')

def plot_solver_error(c=SOLVER_C, times=(1.0, 4.0, 15.0), save_path=None):
    """Part 3(c): where the solver's error lives, and how it falls under refinement."""
    fig, axes = subplots(1, 2, width=4.6, height=3.4)
    _error_profile(axes[0][0], c, times)
    _order_panel(axes[0][1], c, 4.0)

    finish(fig)
    savefig(fig, save_path)
    close(fig)

def _error_curves(c, approximation, times):
    """Relative error of the numeric diffusion flux against exact transport, one pair per time."""
    x, fluxes = solve(c, approximation, times=times)
    for t in times:
        inside = x <= FRONT_FRACTION * t
        # The exact series G, not the default interpolation: the interpolation's own ~1.5%
        # would otherwise sit inside the very error being measured.
        exact = phi_exact(x[inside], t, c, form='series')
        yield x[inside] / t, np.abs(fluxes[t][inside] / exact - 1.0)

def _draw_error_panel(ax, c, times):
    """One panel of the 3(d) figure: both approximations, one line style per time."""
    for _, approximation, color, width in APPROXIMATIONS:
        for (scaled, error), style in zip(_error_curves(c, approximation, times), ERROR_STYLES):
            ax.plot(scaled, error, color=color, linestyle=style, linewidth=width)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=ERROR_FLOOR)
    case_label(ax, f'$c = {c:g}$')
    panel(ax, log=True)

def _error_legend(fig, times):
    """Colour names the approximation, line style the time; both go in the grid's empty slot."""
    handles = [Line2D([], [], color=color, linewidth=width, label=f'{name} diffusion')
               for name, _, color, width in APPROXIMATIONS]
    handles += [Line2D([], [], color=GREY, linestyle=style, linewidth=1.8, label=f'$t = {t:g}$')
                for t, style in zip(times, ERROR_STYLES)]
    # Figure coordinates of the sixth panel of a 2x3 grid, which carries no data.
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.84, 0.26),
               fontsize=LEGEND_SIZE, frameon=False)

def plot_diffusion_error(c_values=C_VALUES, times=ERROR_TIMES, save_path=None):
    """Part 3(d): relative error of both numeric diffusion solutions against exact transport."""
    fig, axes = make_grid(len(c_values))

    for ax, c in zip(axes, c_values):
        _draw_error_panel(ax, c, times)

    label_grid(axes, len(c_values), '$x / vt$',
               r'$\left| \phi_{\mathrm{diff}} / \phi_{\mathrm{exact}} - 1 \right|$')

    finish(fig)
    _error_legend(fig, times)
    savefig(fig, save_path)
    close(fig)
