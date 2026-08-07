import os
import logging
import numpy as np
from homework1.exact_solution import compute_nu0_numerical, compute_nu0_approx
from homework1.diffusion import absorption_balance, solve_diffusion_shooting
from homework1.plots import (
    create_figs_dir,
    plot_flux_components,
    plot_relative_contributions,
    plot_diffusion_comparison,
    plot_diffusion_errors,
    plot_numerical_vs_analytic_diffusion,
    plot_q2_solution_comparison,
    plot_q2_error_profiles,
    plot_q2_convergence
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

    # Question 2 reference figures (original shooting solver, zero-flux outer BC)
    logger.info("\nGenerating Question 2 numerical vs. analytical diffusion comparison...")
    num_comp_path = os.path.join(figs_dir, "diffusion_numerical_vs_analytic.pdf")
    num_err_path = os.path.join(figs_dir, "diffusion_numerical_errors.pdf")
    plot_numerical_vs_analytic_diffusion(c_values=[0.5, 0.7, 0.9], D=1.0/3.0, save_path_comp=num_comp_path, save_path_err=num_err_path)

    # Question 2: solvers using the symmetry-derived current condition at the
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


    logger.info(f"\nAll plots have been successfully generated and saved to:")
    logger.info(f"  {figs_dir}")
    logger.info("==============================")

if __name__ == "__main__":
    main()
