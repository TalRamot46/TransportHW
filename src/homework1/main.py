<<<<<<< Updated upstream
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
)  # MARSHAK_EXTRAPOLATION and critical_dimensions_applied_bc are reused by Q4
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
from homework1.plots_spherical import (
    plot_q4_critical_radius,
    plot_q4_mesh_convergence,
    plot_q4_flux_profiles,
    plot_q5_criticality,
)
from homework1.spherical import (
    build_medium,
    critical_radius,
    analytic_critical_radius,
    mesh_convergence,
    k_eigenvalue,
    dominance_ratio,
    neutron_balance,
)
from homework1.materials import (
    BENCHMARK,
    FISSILE,
    PROMPT_U235,
    solve_material,
)
=======
"""Entry point: runs every question's checks and writes its figures."""
>>>>>>> Stashed changes

import logging
from homework1 import q1, q2, q3, q4, q5
from homework1.figures import figs_dir

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

<<<<<<< Updated upstream
def _report_q4_iteration(c_values, n_cells):
    """Sweeps taken by the source iteration, against the predicted dominance ratio."""
    logger.info("\nSource-iteration cost at the critical radius. The error falls by the")
    logger.info("dominance ratio per sweep, so a ratio near 1 means slow convergence.")
    logger.info("-" * 62)
    logger.info(f"{'c':<6} | {'approximation':<13} | {'dominance ratio':<16} | {'sweeps':<8}")
    logger.info("-" * 62)
    for approximation in ('classical', 'asymptotic'):
        for c in c_values:
            medium = build_medium(1.0, 1.0, c, approximation)
            R = critical_radius(medium, n_cells=n_cells)
            result = k_eigenvalue(R, medium, n_cells=n_cells)
            logger.info(f"{c:<6.2f} | {approximation:<13} | "
                        f"{dominance_ratio(medium, R):<16.6f} | {result.sweeps:<8}")
    logger.info("-" * 62)

def _report_q4_balance(c_values, n_cells):
    """Production against absorption + leakage for the converged critical flux."""
    logger.info("\nNeutron balance at the critical radius: production must equal")
    logger.info("absorption + leakage.")
    logger.info("-" * 88)
    logger.info(f"{'c':<6} | {'approximation':<13} | {'production':<12} | "
                f"{'absorption':<12} | {'leakage':<12} | {'residual':<10}")
    logger.info("-" * 88)
    for approximation in ('classical', 'asymptotic'):
        for c in c_values:
            medium = build_medium(1.0, 1.0, c, approximation)
            R = critical_radius(medium, n_cells=n_cells)
            production, absorption, leakage, residual = neutron_balance(
                R, medium, n_cells=n_cells)
            logger.info(f"{c:<6.2f} | {approximation:<13} | {production:<12.6f} | "
                        f"{absorption:<12.6f} | {leakage:<12.6f} | {residual:<10.2e}")
    logger.info("-" * 88)

def _report_q4_boundary(c_values, n_cells):
    """
    The extrapolated zero against the Robin condition applied at the physical
    surface, with the analytic counterpart of each taken from Question 3.
    """
    logger.info("\nBoundary treatment: extrapolated zero at R + z0 against the condition")
    logger.info("phi + z0 phi' = 0 applied at R. Analytic references are Question 3's")
    logger.info("critical_dimensions ('marshak') and critical_dimensions_applied_bc.")
    logger.info("-" * 78)
    logger.info(f"{'c':<6} | {'extrapolated':<13} | {'Robin':<13} | "
                f"{'Robin analytic':<15} | {'difference %':<12}")
    logger.info("-" * 78)
    for c in c_values:
        medium = build_medium(1.0, 1.0, c, 'classical')
        R_ext = critical_radius(medium, n_cells=n_cells)
        R_rob = critical_radius(medium, n_cells=n_cells, boundary='robin')
        _, R_rob_ana = critical_dimensions_applied_bc(c, MARSHAK_EXTRAPOLATION)
        delta = (R_rob - R_ext) / R_ext * 100.0
        logger.info(f"{c:<6.2f} | {R_ext:<13.8f} | {R_rob:<13.8f} | "
                    f"{R_rob_ana:<15.8f} | {delta:<+12.3f}")
    logger.info("-" * 78)

def report_q4(figs_dir, n_cells=400):
    """
    Question 4: critical radius from the spherical k = 1 search, against the
    analytic radius of the same approximation.
    """
    c_values = [1.02, 1.05, 1.1, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0]

    logger.info("\n=== Homework 1 Question 4 ===")
    logger.info("Spherical diffusion, k by Bell & Glasstone source iteration, critical")
    logger.info(f"radius by bisection on k(R) - 1 with {n_cells} cells. Lengths in mfp.")
    logger.info("-" * 82)
    logger.info(f"{'c':<6} | {'approximation':<13} | {'R_c numerical':<14} | "
                f"{'R_c analytic':<14} | {'diff %':<10}")
    logger.info("-" * 82)
    for approximation in ('classical', 'asymptotic'):
        for c in c_values:
            medium = build_medium(1.0, 1.0, c, approximation)
            R_num = critical_radius(medium, n_cells=n_cells)
            R_ana = analytic_critical_radius(medium)
            diff = (R_num - R_ana) / R_ana * 100.0
            logger.info(f"{c:<6.2f} | {approximation:<13} | {R_num:<14.8f} | "
                        f"{R_ana:<14.8f} | {diff:<+10.2e}")
    logger.info("-" * 82)

    # Second-order convergence: the error should fall by four per doubling.
    cells, _, errors = mesh_convergence(build_medium(1.0, 1.0, 1.5, 'classical'))
    ratios = errors[:-1] / errors[1:]
    logger.info(f"Mesh refinement at c = 1.5: error {errors[0]:.2e} -> {errors[-1]:.2e} "
                f"over N = {cells[0]} -> {cells[-1]}, ratio per doubling "
                f"{ratios.min():.2f}-{ratios.max():.2f}")

    _report_q4_iteration(c_values, n_cells)
    _report_q4_balance(c_values, n_cells)
    _report_q4_boundary(c_values, n_cells)

    logger.info("\nGenerating Question 4 figures...")
    plot_q4_critical_radius(
        c_values, n_cells=n_cells,
        save_path=os.path.join(figs_dir, "q4_critical_radius.pdf"))
    plot_q4_mesh_convergence(
        save_path=os.path.join(figs_dir, "q4_mesh_convergence.pdf"))
    plot_q4_flux_profiles(
        n_cells=n_cells, save_path=os.path.join(figs_dir, "q4_flux_profiles.pdf"))

def report_q5(figs_dir, n_cells=400):
    """
    Question 5: bare critical radius and mass of the two fissile benchmark
    materials, under both diffusion approximations.
    """
    logger.info("\n=== Homework 1 Question 5 ===")
    logger.info("Bare critical sphere from the Sood et al. one-group benchmark data.")
    logger.info("-" * 94)
    logger.info(f"{'material':<22} | {'approx':<11} | {'c':<7} | {'R_c [cm]':<10} | "
                f"{'R_ana [cm]':<11} | {'M_c [kg]':<9}")
    logger.info("-" * 94)

    materials = [BENCHMARK[name] for name in FISSILE] + [PROMPT_U235]
    for material in materials:
        # Sigma_t must equal the sum of its parts, or the row has a typo.
        residual = abs(material.sigma_t_sum - material.sigma_t)
        if residual > 1e-12:
            logger.warning(f"  {material.name}: Sigma_t off by {residual:.2e}")
        for approximation in ('classical', 'asymptotic'):
            r = solve_material(material, approximation, n_cells=n_cells)
            logger.info(f"{material.name:<22} | {approximation:<11} | {material.c:<7.4f} | "
                        f"{r.R_numerical:<10.4f} | {r.R_analytic:<11.4f} | "
                        f"{r.mass_numerical:<9.3f}")
    logger.info("-" * 94)
    logger.info("'U-235 (prompt variant)' is the cross-section row given in the task")
    logger.info("prompt rather than in the assignment PDF; its own cross sections give")
    logger.info("c = 1.365, not the c = 1.50 quoted alongside them.")

    # The mass goes as R^3, so the boundary treatment is worth quantifying.
    logger.info("\nSensitivity of the mass to the boundary treatment:")
    logger.info("-" * 76)
    logger.info(f"{'material':<10} | {'approx':<11} | {'M extrapolated':<15} | "
                f"{'M Robin':<10} | {'difference %':<12}")
    logger.info("-" * 76)
    for name in FISSILE:
        material = BENCHMARK[name]
        for approximation in ('classical', 'asymptotic'):
            ext = solve_material(material, approximation, n_cells=n_cells)
            rob = solve_material(material, approximation, n_cells=n_cells,
                                 boundary='robin')
            delta = (rob.mass_numerical - ext.mass_numerical) / ext.mass_numerical * 100.0
            logger.info(f"{name:<10} | {approximation:<11} | "
                        f"{ext.mass_numerical:<15.3f} | {rob.mass_numerical:<10.3f} | "
                        f"{delta:<+12.2f}")
    logger.info("-" * 76)

    # The three non-fissile rows of the benchmark table have c <= 1 and so have
    # no bare critical size at all; they are listed for completeness.
    subcritical = [name for name in BENCHMARK if name not in FISSILE]
    logger.info(f"No bare critical sphere exists for {', '.join(subcritical)} "
                f"(c <= 1, non-multiplying).")

    logger.info("\nGenerating Question 5 figure...")
    plot_q5_criticality(n_cells=n_cells,
                        save_path=os.path.join(figs_dir, "q5_criticality.pdf"))

=======
>>>>>>> Stashed changes
def main():
    """Runs Questions 1 to 5 in order."""
    figs = figs_dir()
    for question in (q1, q2, q3, q4, q5):
        question.report(figs)
    logger.info(f"\nAll figures written to {figs}")

if __name__ == "__main__":
    main()
