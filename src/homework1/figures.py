"""Shared matplotlib helpers: Assignment 3's styling, a headless backend, and a safe save.

Assignment 3's look, adopted here: serif text, ticks turned inwards on all four sides, a
light solid grid behind the data, one shared palette and a transparent save. No titles at
all, panel or figure: the report caption carries the description, and whatever a title used
to identify -- the case `c`, the approximation, the material -- is now in the axis label,
the legend, or an annotation inside the panel.
"""

import os
import logging
import matplotlib

# The figures are written to file and never displayed, so the non-interactive
# backend is selected before pyplot is imported.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

logger = logging.getLogger(__name__)

N_COLS = 2

# Panels are smaller than the page they land on, so a figure scaled to \textwidth is
# reduced and a 13 pt label arrives near the report's 11 pt body. Raising these numbers
# makes the *text* smaller, not larger. The ratios are those of the earlier sizes, so the
# page layout of every figure is unchanged.
PANEL_WIDTH = 4.6
PANEL_HEIGHT = 3.2
GRID_WIDTH = 4.2
GRID_HEIGHT = 2.8

LABEL_SIZE = 13
TICK_SIZE = 11
LEGEND_SIZE = 11

# Assignment 3's palette, so both reports read as one set of figures.
NAVY, ORANGE, GREEN, RED, BLUE, PURPLE = ('#2c3e50', '#e67e22', '#27ae60', '#c0392b',
                                          '#2980b9', '#8e44ad')
COLORS = (NAVY, ORANGE, GREEN, RED)

def use_style():
    """The shared rc block: family, sizes and axis weight, applied per figure."""
    plt.rcParams.update({
        'font.family': 'serif',
        'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
        'mathtext.fontset': 'dejavusans',
        'axes.labelsize': LABEL_SIZE,
        'axes.linewidth': 0.8,
        'xtick.labelsize': TICK_SIZE,
        'ytick.labelsize': TICK_SIZE,
        'legend.fontsize': LEGEND_SIZE,
        'figure.facecolor': 'none',
    })

def subplots(n_rows, n_cols, width=PANEL_WIDTH, height=PANEL_HEIGHT):
    """Styled subplot array of the given shape, always two-dimensional."""
    use_style()
    return plt.subplots(n_rows, n_cols, figsize=(width * n_cols, height * n_rows),
                        squeeze=False)

def make_grid(n_panels, width=GRID_WIDTH, height=GRID_HEIGHT):
    """Two-column grid with one panel per case; any spare panel is removed."""
    use_style()
    n_rows = -(-n_panels // N_COLS)  # ceiling division
    fig, axes = plt.subplots(n_rows, N_COLS, figsize=(width * N_COLS, height * n_rows))
    axes = axes.flatten()
    for ax in axes[n_panels:]:
        fig.delaxes(ax)
    return fig, axes[:n_panels]

def case_label(ax, text, xy=(0.03, 0.94)):
    """Names the case a panel shows -- its `c`, its material -- inside the axes.

    A grid with one panel per case needs that label somewhere, and the y-label is the wrong
    place: prefixed with the case it grows long enough to collide with the row above.
    """
    ax.annotate(text, xy=xy, xycoords='axes fraction', fontsize=LABEL_SIZE,
                va='top', ha='left',
                bbox=dict(boxstyle='square,pad=0.2', facecolor='white',
                          edgecolor='none', alpha=0.75))

def panel(ax, xlabel=None, ylabel=None, log=False, legend=None, fontsize=LEGEND_SIZE):
    """Applies the shared panel styling. No title: the report caption carries it."""
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale('log')

    ax.set_axisbelow(True)
    ax.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    if legend:
        ax.legend(loc=legend, fontsize=fontsize, frameon=True)
    #

def finish(fig):
    """Tightens the layout. There is no figure title, by design."""
    fig.tight_layout()

def savefig(fig, save_path):
    """
    Writes `fig` to `save_path`, deleting the target first.

    Overwriting an existing PDF in place fails intermittently on this machine with
    "OSError: [Errno 22] Invalid argument"; writing to a name that does not exist
    yet sidesteps it. Same workaround as docs/build.ps1 uses for the LaTeX output.
    """
    if save_path is None:
        return

    if os.path.exists(save_path):
        os.remove(save_path)

    fig.savefig(save_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close(fig)
    logger.info(f"Saved figure to: {save_path}")

def figs_dir():
    """Creates and returns the figure directory of the report."""
    path = os.path.join("docs", "homework1", "figs")
    os.makedirs(path, exist_ok=True)
    return path
