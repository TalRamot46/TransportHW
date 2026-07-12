import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from homework1.exact_solution import phi_asymptotic, phi_transient, phi_exact
from homework1.diffusion import phi_classical_diffusion, phi_asymptotic_diffusion, solve_diffusion_numerical

logger = logging.getLogger(__name__)

def create_figs_dir():
    figs_dir = os.path.join("docs", "homework1", "figs")
    os.makedirs(figs_dir, exist_ok=True)
    return figs_dir

def plot_flux_components(x, c_values, method='numerical', save_path=None):
    """
    Plots the exact scalar flux, its asymptotic component, and its transient component
    for each c in c_values.
    """
    # Use a professional, clean style
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
    axes = axes.flatten()
    
    for i, c in enumerate(c_values):
        ax = axes[i]
        
        # Calculate components
        phi_as = phi_asymptotic(x, c, method=method)
        phi_tr = phi_transient(x, c)
        phi_ex = phi_as + phi_tr
        
        # Plot
        ax.plot(x, phi_ex, label=r'$\phi(x)$ (Exact)', color='#2c3e50', linewidth=2.0)
        
        if c > 0.0:
            ax.plot(x, phi_as, label=r'$\phi_{as}(x)$ (Asymptotic)', color='#e74c3c', linestyle='--', linewidth=1.8)
        
        ax.plot(x, phi_tr, label=r'$\phi_{tr}(x)$ (Transient)', color='#3498db', linestyle=':', linewidth=1.8)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        if i % 2 == 0:
            ax.set_ylabel('Scalar Flux (log scale)', fontsize=10)
        if i >= 6:
            ax.set_xlabel('$x$ [mean free paths]', fontsize=10)
            
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', shadow=True)
    
    # Hide the 8th subplot (since we have 7 c values)
    fig.delaxes(axes[7])
    
    plt.suptitle(f'Scalar Flux Components ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_relative_contributions(x, c_values, method='numerical', save_path=None):
    """
    Plots the relative contribution of the transient and asymptotic components:
    phi_as/phi and phi_tr/phi as a function of x.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
    axes = axes.flatten()
    
    for i, c in enumerate(c_values):
        ax = axes[i]
        
        # Calculate components
        phi_as = phi_asymptotic(x, c, method=method)
        phi_tr = phi_transient(x, c)
        phi_ex = phi_as + phi_tr
        
        # Relative contributions
        # For c = 0, phi_as = 0, so contribution of phi_tr is 100%
        if c == 0.0:
            rel_as = np.zeros_like(x)
            rel_tr = np.ones_like(x)
        else:
            rel_as = phi_as / phi_ex
            rel_tr = phi_tr / phi_ex
            
        # Plot
        ax.plot(x, rel_as, label=r'$\phi_{as}(x) / \phi(x)$', color='#e74c3c', linewidth=2.0)
        ax.plot(x, rel_tr, label=r'$\phi_{tr}(x) / \phi(x)$', color='#3498db', linewidth=2.0)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        if i % 2 == 0:
            ax.set_ylabel('Relative Contribution', fontsize=10)
        if i >= 6:
            ax.set_xlabel('$x$ [mean free paths]', fontsize=10)
            
        ax.legend(loc='center right', frameon=True, facecolor='white', edgecolor='none', shadow=True)
        
    # Hide the 8th subplot
    fig.delaxes(axes[7])
    
    plt.suptitle(f'Relative Contributions of Components ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_diffusion_comparison(x, c_values, method='numerical', save_path=None):
    """
    Plots a comparison of exact transport, asymptotic transport, classical diffusion,
    and asymptotic diffusion for each c in c_values.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
    axes = axes.flatten()
    
    for i, c in enumerate(c_values):
        ax = axes[i]
        
        # Calculate solutions
        phi_as_val = phi_asymptotic(x, c, method=method)
        phi_tr_val = phi_transient(x, c)
        phi_ex_val = phi_as_val + phi_tr_val
        
        phi_class_val = phi_classical_diffusion(x, c)
        phi_as_diff_val = phi_asymptotic_diffusion(x, c, method=method)
        
        # Plot
        ax.plot(x, phi_ex_val, label=r'$\phi(x)$ (Exact Transport)', color='#2c3e50', linewidth=2.0)
        
        if c > 0.0:
            ax.plot(x, phi_as_val, label=r'$\phi_{as}(x)$ (Asymptotic)', color='#e74c3c', linestyle='--', linewidth=1.5)
            ax.plot(x, phi_as_diff_val, label=r'$\phi_{as,diff}(x)$ (Asymptotic Diff)', color='#2ecc71', linestyle='-.', linewidth=1.8)
        
        ax.plot(x, phi_class_val, label=r'$\phi_{class}(x)$ (Classical Diff)', color='#e67e22', linestyle=':', linewidth=1.8)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        if i % 2 == 0:
            ax.set_ylabel('Scalar Flux (log scale)', fontsize=10)
        if i >= 6:
            ax.set_xlabel('$x$ [mean free paths]', fontsize=10)
            
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', shadow=True)
    
    # Hide the 8th subplot
    fig.delaxes(axes[7])
    
    plt.suptitle(f'Comparison of Transport and Diffusion Solutions ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_diffusion_errors(x, c_values, method='numerical', save_path=None):
    """
    Plots the relative error (%) of classical and asymptotic diffusion approximations
    relative to the exact transport solution.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    fig, axes = plt.subplots(4, 2, figsize=(12, 16), sharex=True)
    axes = axes.flatten()
    
    for i, c in enumerate(c_values):
        ax = axes[i]
        
        # Calculate solutions
        phi_as_val = phi_asymptotic(x, c, method=method)
        phi_tr_val = phi_transient(x, c)
        phi_ex_val = phi_as_val + phi_tr_val
        
        phi_class_val = phi_classical_diffusion(x, c)
        phi_as_diff_val = phi_asymptotic_diffusion(x, c, method=method)
        
        # Compute relative error (%)
        err_class = np.abs(phi_class_val - phi_ex_val) / phi_ex_val * 100.0
        
        ax.plot(x, err_class, label='Classical Diffusion', color='#e67e22', linewidth=2.0)
        
        if c > 0.0:
            err_as_diff = np.abs(phi_as_diff_val - phi_ex_val) / phi_ex_val * 100.0
            ax.plot(x, err_as_diff, label='Asymptotic Diffusion', color='#2ecc71', linewidth=2.0)
        
        # Formatting
        ax.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax.set_yscale('log')
        ax.grid(True, which="both", ls="--", alpha=0.5)
        
        if i % 2 == 0:
            ax.set_ylabel('Relative Error (%) (log scale)', fontsize=10)
        if i >= 6:
            ax.set_xlabel('$x$ [mean free paths]', fontsize=10)
            
        ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', shadow=True)
        
    # Hide the 8th subplot
    fig.delaxes(axes[7])
    
    plt.suptitle(f'Relative Error of Diffusion Approximations ({method.capitalize()} Method)', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path}")
    plt.close()

def plot_numerical_vs_analytic_diffusion(c_values=[0.5, 0.7, 0.9], D=1.0/3.0, save_path_comp=None, save_path_err=None):
    """
    Compares the numerical diffusion solution against the analytical classical diffusion solution
    for a subset of c values on [-a, 0]. Generates:
    1. A plot of the analytic vs numerical solutions.
    2. A plot of the relative error between them.
    """
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # 1. Comparison Plot
    fig1, axes1 = plt.subplots(len(c_values), 1, figsize=(8, 12), sharex=False)
    if len(c_values) == 1:
        axes1 = [axes1]
        
    # 2. Error Plot
    fig2, axes2 = plt.subplots(len(c_values), 1, figsize=(8, 12), sharex=False)
    if len(c_values) == 1:
        axes2 = [axes2]
        
    for i, c in enumerate(c_values):
        # Solve numerically
        x_num, phi_num = solve_diffusion_numerical(c, D=D)
        
        # Compute analytical solution on the same grid
        phi_ana = phi_classical_diffusion(x_num, c, D=D)
        
        # Calculate relative error (%)
        rel_err = np.abs(phi_num - phi_ana) / phi_ana * 100.0
        
        # Plot Comparison
        ax1 = axes1[i]
        ax1.plot(x_num, phi_ana, label='Analytical Classical Diffusion', color='#3498db', linewidth=2.0)
        ax1.plot(x_num, phi_num, label='Numerical (Shooting Method)', color='#e74c3c', linestyle='--', linewidth=2.0)
        ax1.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax1.set_ylabel('Scalar Flux', fontsize=10)
        ax1.set_xlabel('$x$ [mean free paths]', fontsize=10)
        ax1.legend(loc='lower left', frameon=True)
        ax1.grid(True, which="both", ls="--", alpha=0.5)
        
        # Plot Error
        ax2 = axes2[i]
        ax2.plot(x_num, rel_err, color='#9b59b6', linewidth=2.0)
        ax2.set_title(f'$c = {c}$', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Relative Error (%)', fontsize=10)
        ax2.set_yscale('log')
        ax2.set_xlabel('$x$ [mean free paths]', fontsize=10)
        ax2.grid(True, which="both", ls="--", alpha=0.5)
        
    fig1.suptitle('Numerical vs. Analytical Diffusion Solutions', fontsize=14, fontweight='bold', y=0.98)
    fig1.tight_layout()
    if save_path_comp:
        fig1.savefig(save_path_comp, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path_comp}")
    plt.close(fig1)
    
    fig2.suptitle('Relative Error of Numerical Diffusion Solution', fontsize=14, fontweight='bold', y=0.98)
    fig2.tight_layout()
    if save_path_err:
        fig2.savefig(save_path_err, bbox_inches='tight', dpi=300)
        logger.info(f"Saved figure to: {save_path_err}")
    plt.close(fig2)
