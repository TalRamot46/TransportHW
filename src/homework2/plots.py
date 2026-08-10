"""
Question 3(a)(b) figures: the exact planar flux against the two diffusion solutions.

One figure per scattering ratio c, one panel per time t.
"""

import numpy as np
from homework2.exact import phi_exact
from homework2.diffusion import phi_classical_diffusion, phi_asymptotic_diffusion
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
