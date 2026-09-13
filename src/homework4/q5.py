"""Question 5: the Monte Carlo critical radii against the bare-sphere analytic theories."""

import os
import logging
import numpy as np
from collections import namedtuple

from homework1.figures import subplots, panel, finish, savefig, NAVY, ORANGE, GREEN, RED, BLUE
from homework1.tables import log_section, log_table
from homework4.analytic import THEORIES, critical_mass
from homework4.config import SimulationConfig

logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIGS_DIR = os.path.join(BASE_DIR, "docs", "homework4", "figs")
DATA_DIR = os.path.join(BASE_DIR, "docs", "homework4", "data")

# The cached Monte Carlo sweeps of Questions 3 and 4, plotted as two point sets.
SOURCES = (('Question 3', 'probability_study.csv', 'x', RED),
           ('Question 4', 'deep_probability_study.csv', 'x', BLUE))
CURVE_STYLES = ((NAVY, '-'), (GREEN, '-'))
C_GRID = np.linspace(1.02, 2.0, 400)

Study = namedtuple('Study', 'label marker colour c radius mass')
Theory = namedtuple('Theory', 'label colour style radius mass')

def load_study(filename, marker, colour, label):
    """One cached sweep as a Study, with the sub-critical (c <= 1) rows dropped."""
    data = np.loadtxt(os.path.join(DATA_DIR, filename), delimiter=',', skiprows=1)
    data = data[~np.isnan(data[:, 4])]
    # The stored c carries the round-off of P1 + 2 P2, which would split the groups
    # of the table below; the tabulated values are exact multiples of 0.01.
    return Study(label, marker, colour, np.round(data[:, 3], 6), data[:, 4], data[:, 5])

def build_theories(c, config):
    """The three theory curves on the grid `c`, converted from mfp to cm and grams."""
    curves = []
    for (label, radius_mfp), (colour, style) in zip(THEORIES, CURVE_STYLES):
        radius = radius_mfp(c) * config.mfp
        curves.append(Theory(label, colour, style, radius, critical_mass(radius, config.rho)))
    return curves

def _draw(ax, theories, studies, attr, ylabel, log, legend):
    """One panel: both theory curves, then both Monte Carlo point sets over them."""
    for t in theories:
        ax.plot(C_GRID, getattr(t, attr), label=t.label, color=t.colour,
                linestyle=t.style, linewidth=2.0)
    for s in studies:
        # 'x' is a stroked marker, so it takes a colour and a line width, not a face.
        ax.scatter(s.c, getattr(s, attr), marker=s.marker, s=52, color=s.colour,
                   linewidths=1.6, zorder=3, label=f'Monte Carlo, {s.label}')
    panel(ax, '$c = P_1 + 2P_2$', ylabel, log=log, legend=legend)

def _draw_departure(ax, studies, config):
    """Signed departure of every Monte Carlo radius from each theory, in per cent."""
    c = np.concatenate([s.c for s in studies])
    radius = np.concatenate([s.radius for s in studies])
    for (label, radius_mfp), (colour, _) in zip(THEORIES, CURVE_STYLES):
        reference = radius_mfp(c) * config.mfp
        ax.scatter(c, (radius - reference) / reference * 100.0, s=34, color=colour,
                   alpha=0.75, zorder=3, label=label)
    ax.axhline(0.0, color=NAVY, linewidth=1.2)
    panel(ax, '$c = P_1 + 2P_2$', 'Departure of Monte Carlo (%)', legend='upper right')

def plot_comparison(save_path, studies, config):
    """The Question 5 figure: critical radius, critical mass, and the departure panel."""
    theories = build_theories(C_GRID, config)

    fig, axes = subplots(1, 3, width=4.7, height=3.4)
    _draw(axes[0][0], theories, studies, 'radius', 'Critical radius [cm]', False, 'upper right')
    _draw(axes[0][1], theories, studies, 'mass', 'Critical mass [g]', True, 'upper right')
    _draw_departure(axes[0][2], studies, config)

    finish(fig)
    savefig(fig, save_path)

def _comparison_table(studies, config):
    """One row per distinct c: the mean Monte Carlo radius against the three theories."""
    c_all = np.concatenate([s.c for s in studies])
    radius_all = np.concatenate([s.radius for s in studies])

    rows = []
    for c in np.unique(c_all):
        mc = radius_all[c_all == c].mean()
        theory = [f(c)[0] * config.mfp for _, f in THEORIES]
        rows.append([f'{c:.2f}', mc] + theory
                    + [f'{(mc - t) / t * 100.0:+.2f}' for t in theory])
    log_table(['c', 'R_MC [cm]', 'R_class [cm]', 'R_Milne [cm]',
               'dev class %', 'dev Milne %'], rows)

def run_analytic_comparison():
    """Question 5 end to end: load both cached sweeps, print the table, write the figure."""
    config = SimulationConfig()
    log_section('Homework 4 Question 5',
                'Monte Carlo critical radii against the bare-sphere analytic theories.',
                f'Sigma_t = {config.Sigma_t:.4f} cm^-1, mfp = {config.mfp:.4f} cm, '
                f'rho = {config.rho:.1f} g/cm^3.',
                'R = pi/B - z0 for both. They differ in the pair (B, z0).')

    studies = [load_study(filename, marker, colour, label)
               for label, filename, marker, colour in SOURCES]
    _comparison_table(studies, config)
    plot_comparison(os.path.join(FIGS_DIR, 'analytic_comparison.pdf'), studies, config)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    run_analytic_comparison()
