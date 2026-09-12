"""Entry point: runs the verification checks, then writes the Question 3 figures."""

import os
import logging
import numpy as np
from scipy.integrate import quad
from homework2.exact import phi_c1, phi_exact
from homework2.diffusion import (
    diffusion_coefficient,
    phi_classical_diffusion,
    phi_asymptotic_diffusion,
    phi_steady_classical,
    phi_steady_asymptotic,
)
from homework2.solver import solve, step, mass, bulk
from homework2.figures import figs_dir
from homework2.plots import (plot_comparison_for_c, plot_solver, plot_solver_error,
                             plot_diffusion_error, C_VALUES, TIMES)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

SOLVER_TIMES = (1.0, 4.0, 15.0)   # a subset of TIMES: each solve marches the whole way

def _bulk_error(x, phi, t, c):
    """Largest relative departure of the solver from the closed form, inside the bulk."""
    inside = bulk(x, t, diffusion_coefficient(c, 'classical'))
    return np.abs(phi[inside] / phi_classical_diffusion(x, t, c)[inside] - 1.0).max()

def check_normalisation():
    """int phi(x,t;c) dx must equal e^{-(1-c)t}: one particle emitted, decaying by absorption."""
    logger.info("Normalisation, int phi dx vs. exp(-(1-c)t)  [exact series G]")
    logger.info(f"{'c':<6} | {'t':<6} | {'integral':<14} | {'expected':<14} | {'rel. err':<10}")
    for c in C_VALUES:
        for t in (1.0, 4.0, 15.0):
            value, _ = quad(lambda x: float(phi_exact(x, t, c, "series")), -t, t, limit=400)
            expected = np.exp(-(1.0 - c) * t)
            logger.info(f"{c:<6} | {t:<6.0f} | {value:<14.8f} | {expected:<14.8f} | "
                        f"{abs(value / expected - 1.0):<10.2e}")

def check_interpolation_error():
    """The interpolated G costs a few tenths of a percent in normalisation; Paasschens quotes ~2%."""
    logger.info("\nCost of the interpolated G, as the departure of int phi dx from the exact G")
    for t in (0.3, 1.0, 3.0, 10.0, 30.0):
        value, _ = quad(lambda x: float(phi_c1(x, t)), -t, t, limit=400)
        logger.info(f"  t = {t:<5.1f} integral = {value:.8f}  ({(value - 1.0) * 100:+.4f} %)")

def check_front():
    """At the causal front w0 -> 0, the collided bracket cancels and only the uncollided plateau is left."""
    logger.info("\nValue at the front, phi(|x| -> vt) vs. the uncollided plateau e^{-t}/(2t)")
    for t in (1.0, 4.0, 15.0):
        logger.info(f"  t = {t:<5.0f} {float(phi_c1(t * (1.0 - 1e-9), t)):.10f} vs. "
                    f"{np.exp(-t) / (2.0 * t):.10f}")

def check_steady_identity():
    """The steady solution quoted in the assignment is the time-integral of the diffusion Green's function."""
    logger.info("\nSteady identity, int_0^inf phi_diff(x,t) dt vs. the closed-form steady solution")
    logger.info(f"{'c':<6} | {'x':<6} | {'approximation':<13} | {'integral':<14} | {'closed form':<14} | {'rel. err':<10}")
    for c in (0.6, 0.8):
        for x in (0.5, 2.0):
            for name, transient, steady in (('classical', phi_classical_diffusion, phi_steady_classical),
                                            ('asymptotic', phi_asymptotic_diffusion, phi_steady_asymptotic)):
                value, _ = quad(lambda t: float(transient(x, t, c)), 0.0, np.inf, limit=400)
                closed = float(steady(x, c))
                logger.info(f"{c:<6} | {x:<6.1f} | {name:<13} | {value:<14.10f} | {closed:<14.10f} | "
                            f"{abs(value / closed - 1.0):<10.2e}")

def check_diffusion_coefficients():
    """D0(c) = (1-c) nu0^2 must stay positive on both sides of c = 1 and tend to 1/3 there."""
    logger.info("\nAsymptotic diffusion coefficient D0(c) = (1-c) nu0^2  (classical is 1/3)")
    for c in C_VALUES:
        logger.info(f"  c = {c:<5} D0 = {diffusion_coefficient(c, 'asymptotic'):.8f}")
    for c in (1.0 - 1e-6, 1.0 + 1e-6):
        logger.info(f"  c = {c:<10.6f} D0 = {diffusion_coefficient(c, 'asymptotic'):.8f} "
                    f"(limit 1/3 = {1/3:.8f})")

def check_diffusion_limit():
    """At late times the exact solution must relax onto the classical diffusion peak."""
    logger.info("\nDiffusion limit at c = 1, phi(0,t) vs. (4 pi t/3)^{-1/2}")
    for t in (10.0, 100.0, 300.0):
        exact = float(phi_c1(0.0, t))
        diffusive = float(phi_classical_diffusion(0.0, t, 1.0))
        logger.info(f"  t = {t:<6.0f} {exact:.6f} vs. {diffusive:.6f}  "
                    f"({(exact / diffusive - 1.0) * 100:+.2f} %)")

def check_solver_conservation():
    """The explicit scheme leaks only through the far boundary, so the mass is exact; report eq. (42)."""
    logger.info("\nSolver conservation, 2 int u dx vs. exp(-(1-c)t)")
    logger.info(f"{'c':<6} | {'t':<6} | {'mass':<16} | {'expected':<16} | {'rel. err':<10}")
    for c in (0.6, 1.0, 1.5):
        x, fluxes = solve(c, 'classical', times=SOLVER_TIMES)
        for t, phi in fluxes.items():
            expected = np.exp(-(1.0 - c) * t)
            logger.info(f"{c:<6} | {t:<6.0f} | {mass(x, phi):<16.10f} | {expected:<16.10f} | "
                        f"{abs(mass(x, phi) / expected - 1.0):<10.2e}")

def check_solver_against_closed_form():
    """The solver must reproduce the Green's function it discretises, inside the bulk."""
    logger.info("\nSolver vs. the closed form, max relative error within three sigma")
    for c in (0.6, 1.0, 1.5):
        x, fluxes = solve(c, 'classical', times=SOLVER_TIMES)
        errors = '  '.join(f"t = {t:<4.0f} {_bulk_error(x, phi, t, c):.2e}"
                           for t, phi in fluxes.items())
        logger.info(f"  c = {c:<5} {errors}")

def check_solver_order():
    """Halving h at fixed r must quarter the error: the scheme is second order in space."""
    logger.info("\nSolver order, bulk error at t = 4 against the node count")
    previous = None
    for n_nodes in (500, 1000, 2000):
        x, fluxes = solve(1.0, 'classical', times=(4.0,), n_nodes=n_nodes)
        error = _bulk_error(x, fluxes[4.0], 4.0, 1.0)
        ratio = '' if previous is None else f"  (ratio {previous / error:.2f}, expected 4)"
        logger.info(f"  N = {n_nodes:<6} error = {error:.3e}{ratio}")
        previous = error

def check_solver_stability_edge():
    """Past r = 1/2 the scheme does not degrade, it explodes; report eq. (39)."""
    logger.info("\nStability edge, peak amplitude after 400 sweeps of an isolated spike")
    for r in (0.49, 0.51):
        u = np.zeros(101)
        u[50] = 1.0
        for _ in range(400):
            u = step(u, r)
        logger.info(f"  r = {r}  max|u| = {np.abs(u).max():.3e}")

def check_solver_source_treatment():
    """Starting from the smeared delta or from the analytic Gaussian must agree; report eq. (35)."""
    logger.info("\nSource treatment, pulse start vs. warm start at t = 4")
    x, pulse = solve(1.0, 'classical', times=(4.0,), start='pulse')
    _, warm = solve(1.0, 'classical', times=(4.0,), start='warm')
    inside = bulk(x, 4.0, diffusion_coefficient(1.0, 'classical'))
    difference = np.abs(warm[4.0][inside] / pulse[4.0][inside] - 1.0).max()
    logger.info(f"  max relative difference = {difference:.2e}, against a discretisation error of "
                f"{_bulk_error(x, pulse[4.0], 4.0, 1.0):.2e}")

def _error_at(x, phi, t, c, fraction):
    """Relative error of a numeric flux against exact transport, at the node nearest x = fraction vt."""
    j = int(np.argmin(np.abs(x - fraction * t)))
    return phi[j] / float(phi_exact(x[j], t, c, "series")) - 1.0

def check_diffusion_error():
    """Part 3(d): how far each numeric diffusion solution sits from the exact transport flux."""
    logger.info("\nDiffusion against exact transport, relative error at x = 0 and at x = 0.9 vt")
    logger.info(f"{'c':<6} | {'t':<6} | {'classical':<22} | {'asymptotic':<22}")
    for c in C_VALUES:
        runs = {a: solve(c, a, times=SOLVER_TIMES) for a in ('classical', 'asymptotic')}
        for t in SOLVER_TIMES:
            cells = [f"{_error_at(x, fluxes[t], t, c, 0.0):+8.2%} {_error_at(x, fluxes[t], t, c, 0.9):+11.2%}"
                     for x, fluxes in runs.values()]
            logger.info(f"{c:<6} | {t:<6.0f} | {cells[0]:<22} | {cells[1]:<22}")

def generate_figures():
    """Writes one comparison figure per value of c into the report's figure directory."""
    directory = figs_dir()
    logger.info(f"\nGenerating Question 3 figures into {directory} ...")
    for c in C_VALUES:
        plot_comparison_for_c(c, TIMES, os.path.join(directory, f"q3_comparison_c{c:g}.pdf"))
    plot_solver(save_path=os.path.join(directory, "q3c_solver.pdf"))
    plot_solver_error(save_path=os.path.join(directory, "q3c_solver_error.pdf"))
    plot_diffusion_error(save_path=os.path.join(directory, "q3d_diffusion_error.pdf"))

def main():
    logger.info("=== Assignment 2, Question 3: exact transport vs. diffusion ===")
    logger.info(f"Sigma_t = v = 1;  c = {C_VALUES};  t = {tuple(int(t) for t in TIMES)}\n")

    check_normalisation()
    check_interpolation_error()
    check_front()
    check_diffusion_coefficients()
    check_steady_identity()
    check_diffusion_limit()

    logger.info("\n--- Part 3(c): the explicit solver ---")
    check_solver_conservation()
    check_solver_against_closed_form()
    check_solver_order()
    check_solver_stability_edge()
    check_solver_source_treatment()

    logger.info("\n--- Part 3(d): the numeric diffusion solutions against exact transport ---")
    check_diffusion_error()

    generate_figures()

    logger.info("Done.")

if __name__ == "__main__":
    main()
