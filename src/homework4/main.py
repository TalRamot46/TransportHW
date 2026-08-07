import logging
import numpy as np
from homework4 import SimulationConfig
from homework4 import find_critical_radius

# Setup logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)

def main():
    config = SimulationConfig()
    
    logger.info("=== Studention Transport Monte Carlo Simulation ===")
    logger.info("Physical Parameters:")
    logger.info(f"  A = {config.A} g/mol")
    logger.info(f"  density (rho) = {config.rho} g/cm^3")
    logger.info(f"  NA = {config.NA:.2e} atoms/mol")
    logger.info(f"  Sigma_t (micro) = {config.Sigma_t_barn} barn")
    logger.info(f"  Atomic Density (rho_A) = {config.rho_A:.4e} atoms/cm^3")
    logger.info(f"  Macroscopic Cross Section (Sigma_t) = {config.Sigma_t:.4f} cm^-1")
    logger.info(f"  Mean Free Path (mfp) = {config.mfp:.4f} cm")
    logger.info(f"  Probabilities: Absorb={config.p_absorb:.2f}, Scatter={config.p_scatter:.2f}, Fission={config.p_fission:.2f}")
    logger.info("-" * 50)
    
    r_crit = find_critical_radius(config, num_histories=100)
    
    v_crit = (4.0 / 3.0) * np.pi * (r_crit ** 3)
    m_crit = config.rho * v_crit
    
    logger.info("-" * 50)
    logger.info("Results:")
    logger.info(f"  Critical Radius = {r_crit:.4f} cm")
    logger.info(f"  Critical Volume = {v_crit:.4f} cm^3")
    logger.info(f"  Critical Mass   = {m_crit:.4f} g")
    logger.info("==================================================")

if __name__ == "__main__":
    main()
