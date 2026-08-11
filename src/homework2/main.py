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
from homework2.solver import max_relative_error, particle_balance
from homework2.figures import figs_dir
from homework2.plots import (
    plot_comparison_for_c,
    plot_numerical_for_c,
    plot_errors_for_c,
    ANALYTIC_DIFFUSION,
    C_VALUES,
    TIMES,
    X_REACH,
)

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

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

APPROXIMATIONS = ('classical', 'asymptotic')
X_MAX = X_REACH * max(TIMES)

def check_solver():
    """The numerical solver must reproduce the analytic diffusion solution and conserve particles."""
    logger.info("\n=== Question 3(c): Crank-Nicolson solver vs. the analytic diffusion solution ===")
    logger.info(f"{'c':<6} | {'approximation':<13} | {'max rel err':<12} | {'particle balance':<16}")
    for c in C_VALUES:
        for approximation in APPROXIMATIONS:
            error = max_relative_error(c, approximation, TIMES, X_MAX)
            balance = particle_balance(c, approximation, TIMES, X_MAX)
            logger.info(f"{c:<6} | {approximation:<13} | {error:<12.3e} | "
                        f"{min(balance.values()):.8f}")

def _refinement_table(label, values, **run):
    """Prints an error-versus-refinement table with the ratio between successive rows."""
    logger.info(f"\n{label}")
    previous = None
    for name, override in values:
        error = max_relative_error(1.0, 'classical', (1.0, 4.0), X_REACH * 4.0, **run, **override)
        ratio = f"ratio {previous / error:.2f}" if previous else ""
        logger.info(f"  {name:<14} error {error:.3e}   {ratio}")
        previous = error

def check_solver_order():
    """Both discretisations must be second order; each is refined with the other held fixed."""
    _refinement_table("Spatial order (warm start, dt = 1/640 fixed):",
                      [(f"n = {n}", {'n_nodes': n}) for n in (500, 1000, 2000, 4000)],
                      dt=1 / 640, start='warm')
    _refinement_table("Temporal order (warm start, n_nodes = 8000 fixed):",
                      [(f"dt = 1/{int(1/dt)}", {'dt': dt}) for dt in (1/20, 1/40, 1/80, 1/160)],
                      n_nodes=8000, start='warm')

def check_delta_cost():
    """
    Seeding the delta costs a larger error constant than starting from a resolved Gaussian,
    but not a worse order: both fall as h^2.
    """
    logger.info("\nCost of seeding the delta, against a warm start on the same mesh:")
    for n in (200, 400, 800, 1600):
        errors = {start: max_relative_error(1.0, 'classical', (1.0, 4.0), X_REACH * 4.0,
                                            n_nodes=n, start=start)
                  for start in ('warm', 'pulse')}
        logger.info(f"  n = {n:<5} warm {errors['warm']:.3e}   pulse {errors['pulse']:.3e}   "
                    f"({errors['pulse'] / errors['warm']:.1f}x)")

def check_diffusion_errors():
    """Part 3(d): the relative error at the source, where the diffusion solutions peak."""
    logger.info("\n=== Question 3(d): relative error at x = 0, in percent ===")
    logger.info("             " + "".join(f"t = {t:<8.0f}" for t in TIMES))
    for c in C_VALUES:
        for approximation, analytic in ANALYTIC_DIFFUSION.items():
            errors = "".join(
                f"{abs(analytic(0.0, t, c) / phi_exact(0.0, t, c) - 1.0) * 100:<12.2f}"
                for t in TIMES)
            logger.info(f"c={c:<4} {approximation[:5]:<6} {errors}")

def generate_figures():
    """Writes the 3(a)(b), 3(c) and 3(d) figures, one of each per value of c."""
    directory = figs_dir()
    logger.info(f"\nGenerating Question 3 figures into {directory} ...")
    for c in C_VALUES:
        plot_comparison_for_c(c, TIMES, os.path.join(directory, f"q3_comparison_c{c:g}.pdf"))
        plot_numerical_for_c(c, TIMES, os.path.join(directory, f"q3c_numerical_c{c:g}.pdf"))
        plot_errors_for_c(c, TIMES, os.path.join(directory, f"q3d_errors_c{c:g}.pdf"))

def main():
    logger.info("=== Assignment 2, Question 3: exact transport vs. diffusion ===")
    logger.info(f"Sigma_t = v = 1;  c = {C_VALUES};  t = {tuple(int(t) for t in TIMES)}\n")

    check_normalisation()
    check_interpolation_error()
    check_front()
    check_diffusion_coefficients()
    check_steady_identity()
    check_diffusion_limit()
    check_solver()
    check_solver_order()
    check_delta_cost()
    check_diffusion_errors()
    generate_figures()

    logger.info("Done.")

if __name__ == "__main__":
    main()
