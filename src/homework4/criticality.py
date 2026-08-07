import logging
from tqdm import tqdm
from .config import SimulationConfig
from .simulation import SimulationHistory

logger = logging.getLogger(__name__)

def check_criticality(config: SimulationConfig, radius: float, num_histories: int = 100) -> bool:
    """
    Checks if a sphere of given radius is critical.
    A sphere is critical if at least one of the histories goes critical.
    """
    is_critical = False
    pbar = tqdm(total=num_histories, desc=f"Checking R={radius:.4f} cm", leave=False)
    for h in range(1, num_histories + 1):
        history = SimulationHistory(config, radius)
        for count, is_hist_critical in history.run_generations():
            pbar.set_postfix(particles=count)
            pbar.n = h
            pbar.refresh()
            if is_hist_critical:
                is_critical = True
                break
        if is_critical:
            break
    pbar.close()
            
    status_str = "critical sphere" if is_critical else "noncritical sphere"
    logger.info(f"Checking R = {radius:.4f} cm: {status_str}")
    return is_critical

def find_initial_bounds(config: SimulationConfig, initial_radius: float, num_histories: int = 100) -> tuple[float, float]:
    """
    Finds the initial [r_low, r_high] bounds where:
      - r_low is a noncritical sphere
      - r_high is a critical sphere
    """
    logger.info(f"Establishing search bounds. Initial estimate R = {initial_radius:.4f} cm...")
    is_initial_critical = check_criticality(config, initial_radius, num_histories)
    
    if is_initial_critical:
        r_high = initial_radius
        r_low = r_high / 2.0
        while check_criticality(config, r_low, num_histories):
            logger.info(f"  R={r_low:.4f} cm is a critical sphere. Halving again.")
            r_high = r_low
            r_low = r_low / 2.0
    else:
        r_low = initial_radius
        r_high = r_low * 2.0
        while not check_criticality(config, r_high, num_histories):
            logger.info(f"  R={r_high:.4f} cm is a noncritical sphere. Doubling again.")
            r_low = r_high
            r_high = r_high * 2.0
            
    logger.info(f"Bounds established: r_low = {r_low:.4f} cm (noncritical), r_high = {r_high:.4f} cm (critical)")
    return r_low, r_high

def bisection_search(config: SimulationConfig, r_low: float, r_high: float, num_histories: int = 100) -> float:
    """
    Performs bisection search between r_low and r_high until dr < 0.05 * r_mid.
    Returns the critical radius.
    """
    logger.info("Starting bisection search...")
    while True:
        r_mid = (r_low + r_high) / 2.0
        dr = r_high - r_low
        
        logger.info(f"Bisection: interval = [{r_low:.4f}, {r_high:.4f}] cm (dr/r = {dr/r_mid:.3f})")
        
        if dr < 0.05 * r_mid:
            break
            
        is_mid_critical = check_criticality(config, r_mid, num_histories)
        if is_mid_critical:
            r_high = r_mid
        else:
            r_low = r_mid
            
    return (r_low + r_high) / 2.0

def find_critical_radius(config: SimulationConfig, initial_radius: float = None, num_histories: int = 100) -> float:
    """
    Helper function that combines bounds initialization and bisection search.
    """
    if initial_radius is None:
        initial_radius = config.mfp
    r_low, r_high = find_initial_bounds(config, initial_radius, num_histories)
    return bisection_search(config, r_low, r_high, num_histories)
