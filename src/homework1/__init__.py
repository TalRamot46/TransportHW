from homework1.exact_solution import (
    compute_nu0_numerical,
    compute_nu0_approx,
    compute_nu0,
    compute_nu0_magnitude_numerical,
    compute_nu0_magnitude_approx,
    compute_nu0_magnitude,
    compute_N0_plus,
    compute_lambda,
    compute_N_nu,
    phi_asymptotic,
    phi_transient,
    phi_exact
)

from homework1.diffusion import (
    diffusion_coefficients,
    phi_diffusion_analytic,
    phi_classical_diffusion,
    phi_asymptotic_diffusion,
    solve_diffusion_shooting,
    absorption_balance,
    convergence_study
)

from homework1.criticality import (
    critical_dimensions,
    critical_dimensions_applied_bc,
    diffusion_relaxation_length,
    extrapolation_distance,
    extrapolation_distance_table
)

from homework1.spherical import (
    SphericalMedium,
    KResult,
    build_medium,
    buckling,
    analytic_critical_radius,
    k_eigenvalue,
    critical_radius,
    dominance_ratio,
    neutron_balance,
    mesh_convergence
)

from homework1.materials import (
    Material,
    MaterialResult,
    BENCHMARK,
    FISSILE,
    PROMPT_U235,
    critical_mass,
    solve_material,
    solve_all
)

