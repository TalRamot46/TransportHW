import numpy as np
from dataclasses import dataclass

from homework1.spherical import build_medium, critical_radius, analytic_critical_radius

# ---------------------------------------------------------------------------
# Question 5: the Sood, Forster & Parsons one-group benchmark cross sections
# (Progress in Nuclear Energy 42(1), 55-106, 2003), as tabulated in the
# assignment. Cross sections are in cm^-1, densities in g/cm^3.
#
# Only the two fissile materials carry a density, because only they have a
# bare critical mass; the moderators are listed because the assignment table
# lists them, and they exercise the c < 1 guard in build_medium.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Material:
    """One row of the benchmark cross-section table."""
    name: str
    nu: float
    sigma_f: float
    sigma_c: float
    sigma_s: float
    sigma_t: float
    density: float = None

    @property
    def sigma_a(self):
        """Absorption, Sigma_c + Sigma_f."""
        return self.sigma_c + self.sigma_f

    @property
    def nu_sigma_f(self):
        """Fission production, nu Sigma_f."""
        return self.nu * self.sigma_f

    @property
    def c(self):
        """Scattering ratio including fission multiplication."""
        return (self.nu_sigma_f + self.sigma_s) / self.sigma_t

    @property
    def sigma_t_sum(self):
        """Sigma_f + Sigma_c + Sigma_s, for checking the tabulated Sigma_t."""
        return self.sigma_f + self.sigma_c + self.sigma_s

BENCHMARK = {
    'Pu-239': Material('Pu-239', 3.24, 0.081600, 0.019584, 0.225216, 0.32640, 15.7),
    'U-235': Material('U-235', 2.70, 0.065280, 0.013056, 0.248064, 0.32640, 19.0),
    'H2O': Material('H2O', 0.0, 0.0, 0.032640, 0.293760, 0.32640),
    'Fe': Material('Fe', 0.0, 0.0, 0.00046512, 0.23209488, 0.23256),
    'Na': Material('Na', 0.0, 0.0, 0.0, 0.086368032, 0.086368032),
}

# The U-235 row as given in the task prompt rather than in the assignment PDF.
# It is self-consistent in Sigma_t but not in c: its own cross sections give
# c = 1.365, not the 1.50 quoted with it. Kept so both can be reported.
PROMPT_U235 = Material('U-235 (prompt variant)', 2.70, 0.065280, 0.015672,
                       0.180448, 0.26140, 19.0)

FISSILE = ('Pu-239', 'U-235')

def critical_mass(radius, density):
    """Mass of a sphere of the given radius, in kg for radius in cm."""
    return density * (4.0 / 3.0) * np.pi * radius**3 / 1000.0

def solve_material(material, approximation, n_cells=400):
    """
    Numerical and analytic critical radius (cm) and mass (kg) for one material
    under one diffusion approximation.
    """
    medium = build_medium(material.sigma_t, material.sigma_a,
                          material.nu_sigma_f, approximation)
    R_num = critical_radius(medium, n_cells=n_cells)
    R_ana = analytic_critical_radius(medium)
    return {
        'material': material.name,
        'approximation': approximation,
        'c': material.c,
        'D': medium.D,
        'z0': medium.z0,
        'R_numerical': R_num,
        'R_analytic': R_ana,
        'mass_numerical': critical_mass(R_num, material.density),
        'mass_analytic': critical_mass(R_ana, material.density),
    }

def solve_all(materials=None, approximations=('classical', 'asymptotic')):
    """Runs solve_material over every material and approximation requested."""
    if materials is None:
        materials = [BENCHMARK[name] for name in FISSILE]
    return [solve_material(m, a) for m in materials for a in approximations]
