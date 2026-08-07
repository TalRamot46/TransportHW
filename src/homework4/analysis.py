import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from homework4 import SimulationConfig, find_critical_radius

# Setup logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)

# Directory setup
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIGS_DIR = os.path.join(BASE_DIR, "docs", "homework4", "figs")
DATA_DIR = os.path.join(BASE_DIR, "docs", "homework4", "data")
os.makedirs(FIGS_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Apply professional LaTeX-style plotting settings matching the reference plot
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'

def run_density_study():
    logger.info("\n=== TASK 2: Density Study (10 to 100 g/cm^3) ===")
    data_path = os.path.join(DATA_DIR, "density_study.csv")
    
    if os.path.exists(data_path):
        logger.info(f"Loading cached density study data from {data_path}...")
        data = np.loadtxt(data_path, delimiter=',', skiprows=1)
        densities = data[:, 0]
        radii = data[:, 1]
        masses = data[:, 2]
    else:
        densities = np.linspace(10, 100, 11)
        radii = []
        masses = []
        
        for rho in densities:
            logger.info(f"\nEvaluating Density: {rho:.1f} g/cm^3")
            config = SimulationConfig(rho=rho)
            r_crit = find_critical_radius(config, num_histories=100)
            v_crit = (4.0 / 3.0) * np.pi * (r_crit ** 3)
            m_crit = rho * v_crit
            
            radii.append(r_crit)
            masses.append(m_crit)
            logger.info(f"Density = {rho:.1f} g/cm^3 -> R_crit = {r_crit:.4f} cm, M_crit = {m_crit:.4f} g")
        
        # Save cache
        header = 'density,radius,mass'
        np.savetxt(data_path, np.column_stack((densities, radii, masses)), delimiter=',', header=header, comments='')
        logger.info(f"Saved density study data to {data_path}")

    # Plot results in two separate subplots side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    
    # Left subplot: Critical Radius
    ax1.scatter(densities, radii, c=densities, cmap='rainbow', marker='o', s=50, 
                edgecolors='black', linewidths=0.7, zorder=3)
    ax1.set_xlabel('Density (g/cm$^3$)', fontsize=12)
    ax1.set_ylabel('Critical Radius (cm)', fontsize=12)
    ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    ax1.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    
    # Right subplot: Critical Mass
    ax2.scatter(densities, masses, c=densities, cmap='rainbow', marker='s', s=50, 
                edgecolors='black', linewidths=0.7, zorder=3)
    ax2.set_yscale('log')
    ax2.set_xlabel('Density (g/cm$^3$)', fontsize=12)
    ax2.set_ylabel('Critical Mass (g)', fontsize=12)
    ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    ax2.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    
    fig.tight_layout()
    plot_path = os.path.join(FIGS_DIR, "density_analysis.pdf")
    plt.savefig(plot_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close()
    logger.info(f"Density study plot saved as PDF to: {plot_path}")


def run_probability_study():
    logger.info("\n=== TASK 3: Probability Sweep ===")
    data_path = os.path.join(DATA_DIR, "probability_study.csv")
    
    if os.path.exists(data_path):
        logger.info(f"Loading cached probability study data from {data_path}...")
        data = np.loadtxt(data_path, delimiter=',', skiprows=1)
        p0 = data[:, 0]
        p1 = data[:, 1]
        p2 = data[:, 2]
        c_vals_all = data[:, 3]
        radii_all = data[:, 4]
        masses_all = data[:, 5]
    else:
        prob_configs = [
            (0.1, 0.75, 0.15),
            (0.1, 0.70, 0.20),
            (0.1, 0.65, 0.25),
            (0.2, 0.55, 0.25),
            (0.2, 0.50, 0.30),
            (0.2, 0.45, 0.35),
            (0.3, 0.35, 0.30),
            (0.3, 0.30, 0.40),
            (0.3, 0.25, 0.45),
            (0.2, 0.47, 0.33),
            (0.2, 0.49, 0.31),
            (0.2, 0.51, 0.29),
            (0.2, 0.53, 0.27)
        ]
        
        results = []
        for p_abs, p_scat, p_fiss in prob_configs:
            c = p_scat + 2.0 * p_fiss
            logger.info(f"\nEvaluating: P0={p_abs:.2f}, P1={p_scat:.2f}, P2={p_fiss:.2f} (c = {c:.2f})")
            
            if c <= 1.0:
                logger.info("  c <= 1.0: System is noncritical in infinite medium. Infinite critical radius.")
                results.append((p_abs, p_scat, p_fiss, c, np.nan, np.nan))
                continue
                
            config = SimulationConfig(p_absorb=p_abs, p_scatter=p_scat, p_fission=p_fiss)
            r_crit = find_critical_radius(config, num_histories=100)
            v_crit = (4.0 / 3.0) * np.pi * (r_crit ** 3)
            m_crit = config.rho * v_crit
            
            results.append((p_abs, p_scat, p_fiss, c, r_crit, m_crit))
            logger.info(f"  R_crit = {r_crit:.4f} cm, M_crit = {m_crit:.4f} g")
        
        results = np.array(results)
        p0, p1, p2, c_vals_all, radii_all, masses_all = (
            results[:, 0], results[:, 1], results[:, 2], results[:, 3], results[:, 4], results[:, 5]
        )
        
        # Save cache
        header = 'p0,p1,p2,c,radius,mass'
        np.savetxt(data_path, results, delimiter=',', header=header, comments='')
        logger.info(f"Saved probability study data to {data_path}")

    # Filter out non-critical/infinite cases for plotting
    valid_mask = ~np.isnan(radii_all)
    c_vals = c_vals_all[valid_mask]
    radii = radii_all[valid_mask]
    masses = masses_all[valid_mask]
    
    # Plot results in two separate subplots side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    
    # Left subplot: Critical Radius
    ax1.scatter(c_vals, radii, c=c_vals, cmap='rainbow', marker='o', s=50, 
                edgecolors='black', linewidths=0.7, zorder=3)
    ax1.set_xlabel('Mean Secondary Neutrons per Collision ($c = P_1 + 2P_2$)', fontsize=12)
    ax1.set_ylabel('Critical Radius (cm)', fontsize=12)
    ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    ax1.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    
    # Right subplot: Critical Mass
    ax2.scatter(c_vals, masses, c=c_vals, cmap='rainbow', marker='s', s=50, 
                edgecolors='black', linewidths=0.7, zorder=3)
    ax2.set_yscale('log')
    ax2.set_xlabel('Mean Secondary Neutrons per Collision ($c = P_1 + 2P_2$)', fontsize=12)
    ax2.set_ylabel('Critical Mass (g)', fontsize=12)
    ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    ax2.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    
    fig.tight_layout()
    plot_path = os.path.join(FIGS_DIR, "probability_analysis.pdf")
    plt.savefig(plot_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close()
    logger.info(f"Probability sweep plot saved as PDF to: {plot_path}")


def run_deep_probability_study():
    logger.info("\n=== TASK 4: Deep Probability Study ===")
    data_path = os.path.join(DATA_DIR, "deep_probability_study.csv")
    
    if os.path.exists(data_path):
        logger.info(f"Loading cached deep probability study data from {data_path}...")
        data = np.loadtxt(data_path, delimiter=',', skiprows=1)
        p0 = data[:, 0]
        p1 = data[:, 1]
        p2 = data[:, 2]
        c_vals_all = data[:, 3]
        radii_all = data[:, 4]
        masses_all = data[:, 5]
    else:
        deep_prob_configs = [
            # Left side of PDF table
            (0.0, 0.95, 0.05),
            (0.1, 0.75, 0.15),
            (0.2, 0.55, 0.25),
            (0.3, 0.35, 0.35),
            (0.0, 0.90, 0.10),
            (0.1, 0.70, 0.20),
            (0.2, 0.50, 0.30),
            (0.3, 0.30, 0.40),
            (0.0, 0.85, 0.15),
            (0.1, 0.65, 0.25),
            (0.2, 0.45, 0.35),
            (0.3, 0.25, 0.45),
            (0.0, 0.80, 0.20),
            (0.1, 0.60, 0.30),
            (0.2, 0.40, 0.40),
            (0.3, 0.20, 0.50),
            (0.0, 0.70, 0.30),
            # Right side of PDF table
            (0.1, 0.50, 0.40),
            (0.2, 0.30, 0.50),
            (0.3, 0.10, 0.60),
            (0.0, 0.60, 0.40),
            (0.1, 0.40, 0.50),
            (0.2, 0.20, 0.60),
            (0.3, 0.00, 0.70),
            (0.0, 0.50, 0.50),
            (0.1, 0.30, 0.60),
            (0.2, 0.10, 0.70),
            (0.0, 0.40, 0.60),
            (0.1, 0.20, 0.70),
            (0.2, 0.00, 0.80),
            (0.0, 0.20, 0.80),
            (0.1, 0.00, 0.90),
            (0.0, 0.00, 1.00)
        ]
        
        results = []
        for p0_val, p1_val, p2_val in deep_prob_configs:
            c = p1_val + 2.0 * p2_val
            logger.info(f"\nEvaluating: P0={p0_val:.2f}, P1={p1_val:.2f}, P2={p2_val:.2f} (c = {c:.2f})")
            
            if c <= 1.0:
                logger.info("  c <= 1.0: System is noncritical in infinite medium. Infinite critical radius.")
                results.append((p0_val, p1_val, p2_val, c, np.nan, np.nan))
                continue
                
            config = SimulationConfig(p_absorb=p0_val, p_scatter=p1_val, p_fission=p2_val)
            r_crit = find_critical_radius(config, num_histories=100)
            v_crit = (4.0 / 3.0) * np.pi * (r_crit ** 3)
            m_crit = config.rho * v_crit
            
            results.append((p0_val, p1_val, p2_val, c, r_crit, m_crit))
            logger.info(f"  R_crit = {r_crit:.4f} cm, M_crit = {m_crit:.4f} g")
        
        results = np.array(results)
        p0, p1, p2, c_vals_all, radii_all, masses_all = (
            results[:, 0], results[:, 1], results[:, 2], results[:, 3], results[:, 4], results[:, 5]
        )
        
        # Save cache
        header = 'p0,p1,p2,c,radius,mass'
        np.savetxt(data_path, results, delimiter=',', header=header, comments='')
        logger.info(f"Saved deep probability study data to {data_path}")

    # Filter out non-critical/infinite cases for plotting
    valid_mask = ~np.isnan(radii_all)
    c_vals = c_vals_all[valid_mask]
    radii = radii_all[valid_mask]
    masses = masses_all[valid_mask]
    
    # Plot results in two separate subplots side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.set_axisbelow(True)
    ax2.set_axisbelow(True)
    
    # Left subplot: Critical Radius
    ax1.scatter(c_vals, radii, c=c_vals, cmap='rainbow', marker='o', s=50, 
                edgecolors='black', linewidths=0.7, zorder=3)
    ax1.set_xlabel('Mean Secondary Neutrons per Collision ($c = P_1 + 2P_2$)', fontsize=12)
    ax1.set_ylabel('Critical Radius (cm)', fontsize=12)
    ax1.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    ax1.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    
    # Right subplot: Critical Mass
    ax2.scatter(c_vals, masses, c=c_vals, cmap='rainbow', marker='s', s=50, 
                edgecolors='black', linewidths=0.7, zorder=3)
    ax2.set_yscale('log')
    ax2.set_xlabel('Mean Secondary Neutrons per Collision ($c = P_1 + 2P_2$)', fontsize=12)
    ax2.set_ylabel('Critical Mass (g)', fontsize=12)
    ax2.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=10)
    ax2.grid(True, linestyle='-', linewidth=0.5, color='#cccccc')
    
    fig.tight_layout()
    plot_path = os.path.join(FIGS_DIR, "deep_probability_analysis.pdf")
    plt.savefig(plot_path, format='pdf', bbox_inches='tight', transparent=True)
    plt.close()
    logger.info(f"Deep probability sweep plot saved as PDF to: {plot_path}")

def main():
    run_density_study()
    run_probability_study()
    run_deep_probability_study()

if __name__ == "__main__":
    main()
