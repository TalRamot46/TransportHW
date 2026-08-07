import numpy as np

class Studention:
    def __init__(self, position=None, direction=None):
        # Default position is origin [0.0, 0.0, 0.0]
        self.position = np.array(position if position is not None else [0.0, 0.0, 0.0], dtype=float)
        # Default direction is random unit vector
        if direction is not None:
            self.direction = np.array(direction, dtype=float)
        else:
            self.randomize_direction()

    def randomize_direction(self):
        """Generates a random unit vector uniformly distributed on a sphere."""
        phi = 2.0 * np.pi * np.random.random()
        cos_theta = 2.0 * np.random.random() - 1.0
        sin_theta = np.sqrt(1.0 - cos_theta**2)
        self.direction = np.array([
            sin_theta * np.cos(phi),
            sin_theta * np.sin(phi),
            cos_theta
        ])

    def advance(self, distance):
        """Advances the particle along its current direction."""
        self.position += distance * self.direction

    @property
    def radius(self):
        """Calculates distance from origin."""
        return np.linalg.norm(self.position)

    def determine_event(self, p_absorb: float, p_scatter: float):
        """
        Draws a random variable to decide between absorption, scattering, and fission.
        Returns a list of resulting Studention particles:
          - [] if absorbed (0 particles)
          - [self] if scattered (1 particle, direction randomized)
          - [p1, p2] if fission (2 new particles starting from current position)
        """
        event_roll = np.random.random()
        if event_roll < p_absorb:
            # Absorption: Studention is absorbed and disappears
            return []
        elif event_roll < p_absorb + p_scatter:
            # Scattering: Studention scatters in a random direction
            self.randomize_direction()
            return [self]
        else:
            # Fission: Two new Studentions are emitted in random directions
            return [
                Studention(position=self.position),
                Studention(position=self.position)
            ]
