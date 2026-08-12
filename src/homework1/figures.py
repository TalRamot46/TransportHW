"""Shared matplotlib helpers: a headless backend, panel styling, and a safe save."""

import os
import logging
import matplotlib

# The figures are written to file and never displayed, so the non-interactive
# backend is selected before pyplot is imported.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

logger = logging.getLogger(__name__)

N_COLS = 2

def use_style():
    """Applies the same grid style as the other homework figures."""
    style = 'seaborn-v0_8-whitegrid'
    plt.style.use(style if style in plt.style.available else 'default')

def subplots(n_rows, n_cols, width=6.5, height=4.5):
    """Styled subplot array of the given shape, always two-dimensional."""
    use_style()
    return plt.subplots(n_rows, n_cols, figsize=(width * n_cols, height * n_rows),
                        squeeze=False)

def make_grid(n_panels, width=6.0, height=4.0):
    """Two-column grid with one panel per case; any spare panel is removed."""
    use_style()
    n_rows = -(-n_panels // N_COLS)  # ceiling division
    fig, axes = plt.subplots(n_rows, N_COLS, figsize=(width * N_COLS, height * n_rows))
    axes = axes.flatten()
    for ax in axes[n_panels:]:
        fig.delaxes(ax)
    return fig, axes[:n_panels]

def panel(ax, title=None, xlabel=None, ylabel=None, log=False, legend=None, fontsize=8):
    """Applies the shared panel styling: title, labels, log scale, grid and legend."""
    if title:
        ax.set_title(title, fontsize=12, fontweight='bold')
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=10)
    if log:
        ax.set_yscale('log')
    ax.grid(True, which='both', ls='--', alpha=0.5)
    if legend:
        ax.legend(loc=legend, fontsize=fontsize, frameon=True)

def finish(fig, suptitle=None):
    """Adds the figure title and tightens the layout."""
    if suptitle:
        fig.suptitle(suptitle, fontsize=15, fontweight='bold', y=0.995)
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

    fig.savefig(save_path, bbox_inches='tight', dpi=300)
    plt.close(fig)
    logger.info(f"Saved figure to: {save_path}")

def figs_dir():
    """Creates and returns the figure directory of the report."""
    path = os.path.join("docs", "homework1", "figs")
    os.makedirs(path, exist_ok=True)
    return path
