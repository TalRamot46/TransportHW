import numpy as np

class SimulationConfig:
    def __init__(self, A=300.0, rho=30.0, NA=0.6e24, Sigma_t_barn=10.0, 
                 p_absorb=0.2, p_scatter=0.5, p_fission=0.3):
        self.A = A  # g/mol
        self.rho = rho  # g/cm^3
        self.NA = NA  # atoms/mol
        self.Sigma_t_barn = Sigma_t_barn  # barn
        
        # Probabilities
        self.p_absorb = p_absorb
        self.p_scatter = p_scatter
        self.p_fission = p_fission
        assert np.isclose(p_absorb + p_scatter + p_fission, 1.0)
        
        # Calculate derived properties
        # rho_A: atoms/cm^3
        self.rho_A = (self.NA * self.rho) / self.A
        # Sigma_t: cm^-1
        # 1 barn = 10^-24 cm^2
        self.Sigma_t = self.rho_A * (self.Sigma_t_barn * 1e-24)
        # mean free path: cm
        self.mfp = 1.0 / self.Sigma_t
