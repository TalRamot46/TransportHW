"""Shared matplotlib helpers: a headless backend, a panel grid, and a safe save."""

import os
import logging
import matplotlib

# The figures are written to file and never displayed, so the non-interactive
# backend is selected before pyplot is imported.
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

logger = logging.getLogger(__name__)

N_COLS = 3

def use_style():
    """Applies the same grid style as the other homework figures."""
    style = 'seaborn-v0_8-whitegrid'
    plt.style.use(style if style in plt.style.available else 'default')

def make_grid(n_panels):
    """Creates a three-column grid of panels and returns (fig, flattened axes)."""
    use_style()
    n_rows = -(-n_panels // N_COLS)  # ceiling division
    fig, axes = plt.subplots(n_rows, N_COLS, figsize=(5.0 * N_COLS, 3.6 * n_rows))
    return fig, axes.flatten()

def label_grid(axes, n_panels, xlabel, ylabel):
    """Removes the unused panels and labels the outer edges of the grid."""
    for j in range(n_panels, len(axes)):
        axes[j].remove()

    for i in range(n_panels):
        if i % N_COLS == 0:
            axes[i].set_ylabel(ylabel, fontsize=10)
        if i >= n_panels - N_COLS:
            axes[i].set_xlabel(xlabel, fontsize=10)

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
    logger.info(f"Saved figure to: {save_path}")

def close(fig):
    """Releases the figure, so that a batch of them does not accumulate in memory."""
    plt.close(fig)

def figs_dir():
    """Creates and returns the figure directory of the Question 3 report."""
    path = os.path.join("docs", "homework2", "figs")
    os.makedirs(path, exist_ok=True)
    return path
