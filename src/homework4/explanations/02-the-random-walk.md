# 02 — The Random Walk

**One iteration of the loop in `run_generations` is a *collision cycle*, not a fission
generation. Three sampling details in it are easy to get wrong and are all currently right.**

## What one "generation" is

    while current_generation:
        for p in current_generation:
            s = -log(xi) * mfp          # 1. sample the flight
            p.advance(s)                # 2. move
            if p.radius > R: continue   # 3. escaped -> dropped
            next_generation += p.determine_event(...)   # 4. collide

Every surviving particle flies and collides **once** per iteration, whatever the outcome. So a
"generation" is one collision per particle, and a scatter advances the generation counter
exactly as a fission does. The expected population growth per generation is therefore

    0.2 * 0 + 0.5 * 1 + 0.3 * 2 = 1.1

in the base case — not `k`, and not the fission multiplicity. This is why the report's
multiplication parameter is `c = P_1 + 2 P_2`: it is the mean number of particles leaving a
collision, which is exactly the per-generation growth factor of this loop.

## Three sampling details

**The flight is exponential, not uniform.** `s = -log(xi) * mfp` with `xi` uniform on `(0,1]`
is the inverse-CDF sample of `exp(-Sigma_t s)`. Using `1 - xi` would be equally correct;
`np.random.random()` returns `[0,1)`, so `xi = 0` would give `log(0)`. It has probability zero
in floating point but is not impossible — if a run ever dies with a `-inf` flight, this is why.

**Directions are uniform on the sphere, not on the angles.** `randomize_direction` samples
`cos_theta` uniformly on `[-1,1]` and `phi` uniformly on `[0, 2pi)`. Sampling `theta`
uniformly instead would cluster particles at the poles and is the classic error here. This is
correct as written.

**Scattering is isotropic and memoryless.** On a scatter, `determine_event` calls
`randomize_direction()` on the existing particle and returns `[self]` — the same object
continues with a fresh direction and no memory of the old one. On fission it returns two
*new* particles and drops the original.

## The one copy that matters

    return [Studention(position=self.position), Studention(position=self.position)]

Both daughters are constructed from the *same* array. `Studention.__init__` does
`np.array(position, dtype=float)`, which **copies** — so the two daughters own independent
position vectors. If that constructor were ever changed to `np.asarray`, which does not copy
when the input is already a float array, both daughters and the parent would share one buffer
and the entire population would move as a single point. The simulation would still run and
still produce a plausible-looking critical radius, so this would not announce itself.

## Termination

`run_generations` is a generator yielding `(count, is_critical)`, and it stops on one of two
conditions:

- **`count > max_particles` (10000)** — yields `(count, True)` and returns. The history is
  declared critical.
- **the population reaches zero** — the `while` exits and it yields `(0, False)`.

`max_particles` is a constructor default on `SimulationHistory` and is never overridden by
`criticality.py`, so 10000 is the operative threshold everywhere.

`SimulationHistory.run()` is a thin wrapper over the generator, marked deprecated in the source
because `check_criticality` consumes the generator directly in order to drive its progress bar.
Nothing calls `run()`.
