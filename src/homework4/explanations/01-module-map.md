# 01 — The Module Map

**Six modules in a straight line from one particle to a critical mass, and two entry points
that do different things.**

## The files

| file | owns |
|---|---|
| `config.py` | `SimulationConfig`: the material, the three collision probabilities, and the derived `rho_A`, `Sigma_t`, `mfp`. |
| `particle.py` | `Studention`: position, direction, and the collision outcome draw. |
| `simulation.py` | `SimulationHistory`: one source particle and all its descendants, advanced generation by generation. |
| `criticality.py` | Is a given radius critical, bracket it, bisect it. |
| `analysis.py` | The three parameter studies of Questions 2–4, their CSV cache and their figures. |
| `main.py` | The Question 1 base case only. |
| `__init__.py` | Re-exports the public names, so callers write `from homework4 import ...`. |

## Two entry points

- **`python -m homework4.main`** — Question 1 alone. One `find_critical_radius` at the default
  config, printing the radius, volume and mass.
- **`python -m homework4.analysis`** — Questions 2, 3 and 4. Three studies, 57 critical-radius
  searches in total if the cache is cold ([04](04-studies-and-caching.md)).

`main.py` does not call `analysis.py` and vice versa; the base case is deliberately repeated
inside the density study as the `rho = 30` point.

## The call path

    Studention.determine_event        absorb / scatter / fission -> list of particles
      <- SimulationHistory.run_generations   one collision cycle per iteration
      <- criticality.check_criticality       is this radius critical? (bool)
      <- criticality.find_initial_bounds     halve or double until bracketed
      <- criticality.bisection_search        to 5% relative width
      <- criticality.find_critical_radius    the two above, in order
      <- main.main / analysis.run_*_study

Every layer above `run_generations` deals only in booleans and radii — the particle objects
never escape `SimulationHistory`.

## `SimulationConfig` derives, it does not store

Only `A`, `rho`, `NA`, `Sigma_t_barn` and the three probabilities are inputs. Everything the
simulation actually uses is computed in `__init__`:

    rho_A   = NA * rho / A                    atoms/cm^3
    Sigma_t = rho_A * Sigma_t_barn * 1e-24    cm^-1   (the barn conversion)
    mfp     = 1 / Sigma_t                     cm

so changing `rho` in a study changes the mean free path, which is the whole mechanism behind
Question 2. The constructor asserts `p_absorb + p_scatter + p_fission == 1` — the one input
validation in the package.

Note `NA = 0.6e24`, not Avogadro's number to full precision; that is the assignment's value.
