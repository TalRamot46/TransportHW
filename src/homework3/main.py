"""Entry point: runs the verification checks, then writes the Question 3 figures."""

import os
import logging
import numpy as np
from scipy.integrate import quad
from homework3.exact import phi_c1, phi_exact
from homework3.diffusion import (
    diffusion_coefficient,
    phi_classical_diffusion,
    phi_asymptotic_diffusion,
    phi_steady_classical,
    phi_steady_asymptotic,
)
from homework3.figures import figs_dir
from homework3.plots import plot_comparison_for_c, C_VALUES, TIMES

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

def generate_figures():
    """Writes one comparison figure per value of c into the report's figure directory."""
    directory = figs_dir()
    logger.info(f"\nGenerating Question 3 figures into {directory} ...")
    for c in C_VALUES:
        plot_comparison_for_c(c, TIMES, os.path.join(directory, f"q3_comparison_c{c:g}.pdf"))

def main():
    logger.info("=== Assignment 2, Question 3: exact transport vs. diffusion ===")
    logger.info(f"Sigma_t = v = 1;  c = {C_VALUES};  t = {tuple(int(t) for t in TIMES)}\n")

    check_normalisation()
    check_interpolation_error()
    check_front()
    check_diffusion_coefficients()
    check_steady_identity()
    check_diffusion_limit()
    generate_figures()

    logger.info("Done.")

if __name__ == "__main__":
    main()
