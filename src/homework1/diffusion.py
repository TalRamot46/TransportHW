import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from homework1.exact_solution import compute_nu0

def phi_classical_diffusion(x, c, D=1.0/3.0):
    """
    Computes the classical diffusion solution for an isotropic plane delta-source:
    phi(x) = 1 / (2 * sqrt(D * (1 - c))) * exp(-|x| * sqrt((1 - c) / D))
    
    Parameters:
    x : float or numpy.ndarray
        Distance from the source (in units of mean free path).
    c : float
        Scattering ratio (c < 1).
    D : float, optional
        Diffusion coefficient. Default is 1.0/3.0 (standard P1 approximation).
    """
    if c >= 1.0:
        raise ValueError("Scattering ratio c must be less than 1.0 for classical diffusion.")
    if c < 0.0:
        raise ValueError("Scattering ratio c must be non-negative.")
    
    kappa = np.sqrt((1.0 - c) / D)
    coeff = 1.0 / (2.0 * np.sqrt(D * (1.0 - c)))
    return coeff * np.exp(-kappa * np.abs(x))

def phi_asymptotic_diffusion(x, c, method='numerical'):
    """
    Computes the asymptotic diffusion solution for an isotropic plane delta-source:
    phi(x) = 1 / (2 * (1 - c) * nu0) * exp(-|x| / nu0)
    
    For c = 0, returns 0.0 (or zeros of the same shape as x).
    """
    if c < 0.0 or c >= 1.0:
        raise ValueError("Scattering ratio c must be in [0, 1).")
    
    if c == 0.0:
        if isinstance(x, np.ndarray):
            return np.zeros_like(x, dtype=float)
        return 0.0
    
    nu0 = compute_nu0(c, method=method)
    if nu0 is None:
        if isinstance(x, np.ndarray):
            return np.zeros_like(x, dtype=float)
        return 0.0
    
    coeff = 1.0 / (2.0 * (1.0 - c) * nu0)
    return coeff * np.exp(-np.abs(x) / nu0)

def solve_diffusion_numerical(c, D=1.0/3.0, num_points=500):
    """
    Solves the time-independent diffusion equation numerically on [-a, 0]
    using a shooting method.
    
    Parameters:
    c : float
        Scattering ratio.
    D : float, optional
        Diffusion coefficient. Default is 1.0/3.0.
    num_points : int, optional
        Number of points in the output grid. Default is 500.
        
    Returns:
    x_grid : numpy.ndarray
        Coordinate grid from -a to 0.
    phi_num : numpy.ndarray
        Numerical solution for the scalar flux.
    """
    if c >= 1.0 or c <= 0.0:
        raise ValueError("c must be in (0, 1) for the shooting method.")
    
    kappa = np.sqrt((1.0 - c) / D)
    a = 10.0 / kappa  # 3 decay lengths
    
    # 1st-order system: dy/dx = [y2, (1-c)/D * y1]
    def system(x, Y):
        y1, y2 = Y
        return [y2, ((1.0 - c) / D) * y1]
    
    # Objective function for shooting: f(s) = y2(0; s) - 1/(2D)
    target = 1.0 / (2.0 * D)
    
    def shoot(s):
        # Integrate from -a to 0
        sol = solve_ivp(system, [-a, 0.0], [0.0, s], method='RK45', rtol=1e-10, atol=1e-12)
        # Return value of y2 at the end (x=0) minus target
        return sol.y[1, -1] - target
    
    # Bracket the root. s_min = 0.0 gives shoot(s_min) = -target < 0.
    # Since exact s = target * exp(-kappa * a), target * 1.5 is a safe upper bound.
    s_min = 0.0
    s_max = target * 1.5
    
    # Solve for shooting parameter s
    s_opt = brentq(shoot, s_min, s_max)
    
    # Re-integrate with optimal shooting parameter to get dense output on a uniform grid
    x_grid = np.linspace(-a, 0.0, num_points)
    sol = solve_ivp(system, [-a, 0.0], [0.0, s_opt], t_eval=x_grid, method='RK45', rtol=1e-10, atol=1e-12)
    
    return x_grid, sol.y[0]

