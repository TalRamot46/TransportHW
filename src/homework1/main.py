import os
import logging
import numpy as np
from homework1.exact_solution import compute_nu0_numerical, compute_nu0_approx
from homework1.diffusion import absorption_balance, solve_diffusion_shooting
from homework1.criticality import (
    critical_dimensions,
    critical_dimensions_applied_bc,
    extrapolation_distance,
    MARSHAK_EXTRAPOLATION,
    MARK_EXTRAPOLATION,
    CASE_TABLE_8_C,
    CASE_TABLE_8_K0,
)
from homework1.plots import (
    create_figs_dir,
    plot_flux_components,
    plot_relative_contributions,
    plot_diffusion_comparison,
    plot_diffusion_errors,
    plot_q2_solution_comparison,
    plot_q2_error_profiles,
    plot_q2_convergence,
    plot_q3_critical_dimensions,
    plot_q3_method_comparison,
    plot_q3_extrapolation_distance
)

# Setup logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)

def main():
    # Local configuration variables
    method = "both"
    x_min = 1e-3
    x_max = 5.0
    x_points = 500

    c_values = [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    x = np.linspace(x_min, x_max, x_points)

    logger.info("=== Homework 1 Question 1a ===")
    logger.info("Slab-geometry Boltzmann equation analyzed by Case (1953)")
    logger.info("-" * 65)
    
    # Print comparison table for nu0
    logger.info(f"{'c':<6} | {'nu0 (numerical)':<18} | {'nu0 (approx)':<18} | {'Relative Error (%)':<20}")
    logger.info("-" * 65)
    for c in c_values:
        if c == 0.0:
            logger.info(f"{c:<6} | {'None':<18} | {'None':<18} | {'N/A':<20}")
            continue
        nu0_num = compute_nu0_numerical(c)
        nu0_app = compute_nu0_approx(c)
        rel_err = np.abs(nu0_num - nu0_app) / nu0_num * 100.0
        logger.info(f"{c:<6} | {nu0_num:<18.10f} | {nu0_app:<18.10f} | {rel_err:<20.8f}%")
    logger.info("-" * 65)

    # Generate plots
    figs_dir = create_figs_dir()
    
    methods_to_run = ["numerical", "approx"] if method == "both" else [method]
    
    for method_item in methods_to_run:
        logger.info(f"\nGenerating plots using method: {method_item}...")

        flux_path = os.path.join(figs_dir, f"flux_components_{method_item}.pdf")
        plot_flux_components(x, c_values, method=method_item, save_path=flux_path)

        contrib_path = os.path.join(figs_dir, f"relative_contributions_{method_item}.pdf")
        plot_relative_contributions(x, c_values, method=method_item, save_path=contrib_path)

        diff_path = os.path.join(figs_dir, f"diffusion_comparison_{method_item}.pdf")
        plot_diffusion_comparison(x, c_values, method=method_item, save_path=diff_path)

        error_path = os.path.join(figs_dir, f"diffusion_errors_{method_item}.pdf")
        plot_diffusion_errors(x, c_values, method=method_item, save_path=error_path)

    # Question 2: solver using the symmetry-derived current condition at the
    # source and a radiation condition at the outer boundary.
    q2_c_values = [0.5, 0.7, 0.9]
    q2_approximations = ("classical", "asymptotic")

    logger.info("\n=== Homework 1 Question 2 ===")
    logger.info("Neutron balance check, Sigma_a * integral(phi dx) (should be 1):")
    logger.info(f"{'c':<6} | {'approximation':<13} | {'balance':<14}")
    logger.info("-" * 40)
    for approximation in q2_approximations:
        for c in q2_c_values:
            x_num, phi_num = solve_diffusion_shooting(c, approximation)
            balance = absorption_balance(x_num, phi_num, c, approximation)
            logger.info(f"{c:<6} | {approximation:<13} | {balance:<14.10f}")
    logger.info("-" * 40)

    logger.info("\nGenerating Question 2 solution comparison, error profiles and convergence...")
    plot_q2_solution_comparison(
        q2_c_values, q2_approximations,
        save_path=os.path.join(figs_dir, "q2_solution_comparison.pdf"))
    plot_q2_error_profiles(
        q2_c_values, q2_approximations,
        save_path=os.path.join(figs_dir, "q2_error_profiles.pdf"))
    plot_q2_convergence(
        q2_c_values, q2_approximations,
        save_path=os.path.join(figs_dir, "q2_convergence.pdf"))


    # Question 3: critical dimensions of a bare multiplying medium, from the
    # asymptotic criticality relations with the approximate |nu0(c)| and z0(c).
    q3_c_values = [1.02, 1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

    logger.info("\n=== Homework 1 Question 3 ===")
    logger.info("a/2 = (pi/2)|nu0(c)| - z0(c),  Sigma_t R_c = pi |nu0(c)| - z0(c)")
    logger.info("All lengths in mean free paths. 'fit' uses the approximate formulas for")
    logger.info("both inputs; 'ref' uses the transcendental root for |nu0| and Case's")
    logger.info("Table 23 for z0.")
    logger.info("-" * 88)
    logger.info(f"{'c':<6} | {'|nu0| fit':<11} | {'|nu0| ref':<11} | {'z0 fit':<9} | "
                f"{'z0 ref':<9} | {'a/2 fit':<9} | {'R_c fit':<9} | {'a/2 err %':<9}")
    logger.info("-" * 88)
    reference = {c: critical_dimensions(c, 'transport-ref') for c in q3_c_values}
    for c in q3_c_values:
        nu0_fit, z0_fit, half_fit, radius_fit = critical_dimensions(c, 'transport')
        nu0_ref, z0_ref, half_ref, _ = reference[c]
        err = abs(half_fit - half_ref) / half_ref * 100.0
        logger.info(f"{c:<6.2f} | {nu0_fit:<11.6f} | {nu0_ref:<11.6f} | {z0_fit:<9.6f} | "
                    f"{z0_ref:<9.6f} | {half_fit:<9.6f} | {radius_fit:<9.6f} | {err:<9.4f}")
    logger.info("-" * 88)

    # Verification of the eigenvalue branch against Case's Table 8 Part II.
    k0_errors = [abs(1.0 / critical_dimensions(c, 'transport-ref')[0] - k0)
                 for c, k0 in zip(CASE_TABLE_8_C[1:], CASE_TABLE_8_K0[1:])]
    logger.info(f"Eigenvalue check: max |k0 - Case Table 8| = {max(k0_errors):.2e} "
                f"over c = 1.1 ... 2.0")

    # The printed z0 fit carries a minus sign on the quadratic term, which has
    # the wrong sign against Case's Table 23; both are reported.
    for correction in (-0.0199, 0.0199):
        errors = [abs(extrapolation_distance(c, correction) - reference[c][1])
                  / reference[c][1] * 100.0 for c in q3_c_values]
        logger.info(f"z0 fit with q = {correction:+.4f}: max error vs. Table 23 = "
                    f"{max(errors):.4f} %")

    # Parts 3(b) and 3(c): the same relations under classical diffusion, with
    # the constant extrapolation distances of the Marshak and Mark conditions.
    logger.info("\nDiffusion approximations, Marshak (z0 = 2/3) and Mark (z0 = 1/sqrt(3)),")
    logger.info("against the exact-transport result of part 3(a). 'BC applied' solves the")
    logger.info("boundary condition on the flux shape instead of using an extrapolated zero.")
    logger.info("-" * 100)
    logger.info(f"{'c':<6} | {'a/2 transp':<11} | {'a/2 Marshak':<12} | {'a/2 Mark':<11} | "
                f"{'R transp':<10} | {'R Marshak':<11} | {'R Mark':<10} | {'Marshak d%':<10}")
    logger.info("-" * 100)
    for c in q3_c_values:
        _, _, half_t, radius_t = critical_dimensions(c, 'transport')
        _, _, half_ma, radius_ma = critical_dimensions(c, 'marshak')
        _, _, half_mk, radius_mk = critical_dimensions(c, 'mark')
        delta = (half_ma - half_t) / half_t * 100.0
        logger.info(f"{c:<6.2f} | {half_t:<11.6f} | {half_ma:<12.6f} | {half_mk:<11.6f} | "
                    f"{radius_t:<10.6f} | {radius_ma:<11.6f} | {radius_mk:<10.6f} | "
                    f"{delta:<+10.2f}")
    logger.info("-" * 100)

    for name, extrapolation in (('Marshak', MARSHAK_EXTRAPOLATION),
                                ('Mark', MARK_EXTRAPOLATION)):
        for c in (1.02, 2.0):
            half_bc, radius_bc = critical_dimensions_applied_bc(c, extrapolation)
            _, _, half, radius = critical_dimensions(
                c, 'marshak' if name == 'Marshak' else 'mark')
            logger.info(f"{name} at c = {c}: extrapolated a/2 = {half:.6f} vs. "
                        f"BC applied {half_bc:.6f}; R_c = {radius:.6f} vs. {radius_bc:.6f}")

    logger.info("\nGenerating Question 3 critical-dimension figures...")
    q3_c_grid = np.linspace(1.02, 2.0, 400)
    plot_q3_critical_dimensions(
        q3_c_grid, save_path=os.path.join(figs_dir, "q3_critical_dimensions.pdf"))
    plot_q3_method_comparison(
        q3_c_grid, save_path=os.path.join(figs_dir, "q3_method_comparison.pdf"))
    plot_q3_extrapolation_distance(
        save_path=os.path.join(figs_dir, "q3_extrapolation_distance.pdf"))

    logger.info(f"\nAll plots have been successfully generated and saved to:")
    logger.info(f"  {figs_dir}")
    logger.info("==============================")

if __name__ == "__main__":
    main()
