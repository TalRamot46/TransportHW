# 03 — The Criticality Search

**"Critical" here is a population-survival threshold, not `k = 1`, and the predicate the
bisection is searching on is random — which changes what the answer means.**

## The operational definition

Three nested definitions, each in a different function:

1. **A history is critical** if its population ever exceeds `max_particles = 10000`
   (`run_generations`).
2. **A radius is critical** if **at least one** of `num_histories = 100` histories is critical
   (`check_criticality`).
3. **The critical radius** is where that flips (`bisection_search`).

Step 2 is the one to keep in mind. It is a statement about the *survival probability* of a
single source particle, not about a mean multiplication factor. A supercritical system still
extinguishes most of its histories early — the chain dies out before it can grow — so the
probability that at least one of 100 survives is what is actually being thresholded.

Two consequences follow directly:

- **The answer depends on `num_histories`.** More histories means more chances for one to run
  away, so a larger `N` will call smaller spheres critical. `num_histories` is a parameter of
  the definition, not just of the precision. Every study uses 100.
- **`check_criticality` breaks early**, on the first critical history, so the work done varies
  from 1 history (clearly supercritical) to 100 (clearly subcritical). The subcritical side is
  the expensive one. That is also why the `tqdm` bar sets `pbar.n = h` by hand rather than
  letting the loop drive it.

## Bracketing, then bisecting

`find_initial_bounds` starts from `initial_radius` (defaulting to one mean free path) and
either halves until it finds a noncritical radius or doubles until it finds a critical one,
returning `(r_low, r_high)` = (noncritical, critical).

`bisection_search` then halves the interval until

    dr < 0.05 * r_mid

and returns the midpoint. So the **stated tolerance is 5% of the interval width**, and since
the answer is the midpoint, the radius carries about **±2.5%** — before any Monte Carlo
uncertainty. Roughly 4–5 evaluations, each up to 100 histories.

## Why this is not an ordinary bisection

Bisection assumes a *deterministic monotone* predicate: below the root always false, above
always true. Here the predicate is random. Near the critical radius the survival probability
passes smoothly through the value where "at least one of 100" becomes a coin flip, so the same
radius can test critical on one call and noncritical on the next.

The consequences are worth being explicit about, because they are not visible in the output:

- The search **cannot converge below the width of that stochastic transition region**, no
  matter how small the tolerance is set. Tightening `0.05` would buy precision that is not
  there.
- A single unlucky draw is **permanent**: once `bisection_search` moves a bound it never
  revisits it, so one spurious "critical" at a subcritical radius drags the whole remaining
  search down with it. There is no re-test and no confidence interval.
- Repeated runs of the same configuration will differ by more than the 2.5% the tolerance
  suggests. The scatter visible between neighbouring points in the Question 3 and 4 figures is
  this, not physics.

Nothing in the code is wrong here — this is the standard way the assignment is posed. It is
recorded so that the scatter in the figures is read as method noise rather than as a real
dependence on `P_0, P_1, P_2` separately, which is precisely what the report's Question 3 and 4
Discussion sections have to decide. See [05](05-limitations.md) for what would tighten it.
