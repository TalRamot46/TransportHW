"""Assignment 1's matplotlib helpers, pointed at this report's figure directory."""

import os
from homework1.figures import subplots, panel, finish, savefig  # noqa: F401

def figs_dir():
    """Creates and returns the figure directory of the report."""
    path = os.path.join("docs", "homework3", "figs")
    os.makedirs(path, exist_ok=True)
    return path
