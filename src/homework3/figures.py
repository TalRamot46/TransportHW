"""Figure helpers for this report: Assignment 4's plot styling.

Assignment 4's look, adopted here: serif text, no panel titles (the caption carries them,
and `column_title` is the one exception, for grid headers), ticks turned inwards on all
four sides, a light solid grid behind the data, and a transparent save. Mathtext stays on the sans set, so symbols read a little lighter than
the words around them; switch mathtext.fontset to 'cm' to match Assignment 4 exactly.
"""

import os
import logging
import matplotlib

# The figures are written to file and never displayed, so the non-interactive
# backend is selected before pyplot is imported.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

logger = logging.getLogger(__name__)

# Panels are sized so that a figure scaled to \textwidth is barely reduced: at three
# columns the scale is about 0.7, so a 13 pt label lands near 9 pt on the page, against
# the report's 11 pt body. Raising these numbers makes the *text* smaller, not larger.
PANEL_WIDTH = 3.6
PANEL_HEIGHT = 3.4

LABEL_SIZE = 13
TICK_SIZE = 11
LEGEND_SIZE = 11

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

def two_over_one(width=PANEL_WIDTH, height=PANEL_HEIGHT):
    """Two panels side by side with a third spanning the width beneath them.

    Three panels in a single row would each be a third of the text width, which is what
    made the old figures unreadable; this gives every panel half of it or more.
    """
    use_style()
    fig = plt.figure(figsize=(2.0 * width, 2.0 * height))
    grid = fig.add_gridspec(2, 2)
    return fig, (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
                 fig.add_subplot(grid[1, :]))

def column_title(ax, text):
    """A column header on the top panel of a grid: the one title the style allows.

    A grid whose columns differ by one variable would otherwise have to repeat that
    variable in every y-label, which is what it is here to avoid.
    """
    ax.set_title(text, fontsize=LABEL_SIZE, pad=8)

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
    path = os.path.join("docs", "homework3", "figs")
    os.makedirs(path, exist_ok=True)
    return path
