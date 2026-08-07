import numpy as np
from .particle import Studention
from .config import SimulationConfig

class SimulationHistory:
    """Manages a single history (1 source particle and its descendants)."""
    def __init__(self, config: SimulationConfig, radius: float, max_particles: int = 10000):
        self.config = config
        self.radius = radius
        self.max_particles = max_particles
        
    def run_generations(self):
        """
        Generator that yields (particle_count, is_critical) at each generation.
        Terminates when 0 particles remain or the limit is exceeded.
        """
        current_generation = [Studention()]
        
        while len(current_generation) > 0:
            count = len(current_generation)
            if count > self.max_particles:
                yield count, True
                return
            
            yield count, False
            
            next_generation = []
            for p in current_generation:
                # 1. Draw distance to next event
                xi = np.random.random()
                s = -np.log(xi) * self.config.mfp
                
                # 2. Advance particle
                p.advance(s)
                
                # 3. Check boundary escape
                if p.radius > self.radius:
                    continue  # Escaped, deleted
                
                # 4. Determine event type and collect resulting particles
                survivors = p.determine_event(self.config.p_absorb, self.config.p_scatter)
                next_generation.extend(survivors)
                    
            current_generation = next_generation
            
        yield 0, False
    
    # deprecated to use tqdm.
    def run(self) -> bool:
        """
        Runs the history until 0 particles remain or active particles > max_particles.
        Returns True if critical (exceeded max_particles), False if noncritical (0 particles).
        """
        is_critical = False
        for _, is_crit in self.run_generations():
            if is_crit:
                is_critical = True
                break
        return is_critical
