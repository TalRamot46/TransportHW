# 05 — Known Limitations

**What the current code cannot tell you, roughly in order of how much it affects the numbers
in the report. Nothing here is a bug; all of it is scope.**

## 1. No seed, so nothing is reproducible

Every draw goes through the global `np.random` and nothing is seeded. Two runs of the same
configuration give different critical radii, and there is no way to reproduce a figure from the
code alone — only from the cached CSV ([04](04-studies-and-caching.md)).

The fix is small if it is ever wanted: give `SimulationConfig` a `seed`, build a
`np.random.default_rng(seed)` once, and thread it through `Studention` and `SimulationHistory`
in place of the three `np.random.random()` calls. It would also make the studies comparable to
each other, since each would start from a known state.

## 2. No uncertainty on any quoted number

`find_critical_radius` returns a single float with no error bar. The bisection's 5% stopping
rule bounds the *interval*, not the uncertainty of the result, and the stochastic predicate
([03](03-criticality-search.md)) means repeat runs scatter by more than that.

This matters directly for the report's Question 3 and 4 Discussion sections, which ask whether
`c` alone determines the critical mass or whether `P_0, P_1, P_2` have separate effects.
**Answering that requires knowing the run-to-run spread**, or the scatter between
equal-`c` configurations cannot be attributed. The cheapest way to get it: pick two or three
configurations, run `find_critical_radius` five times each, and quote the spread.

## 3. The threshold definition is a parameter, not a constant

`max_particles = 10000` and `num_histories = 100` both enter the *definition* of critical, not
just its precision, and neither is swept. The assignment fixes both, so the numbers are the
requested ones — but "the critical radius" here means "the radius at which one of 100 source
particles has a reasonable chance of producing 10000 descendants", and that is not identical to
`k = 1`.

Assignment 1 Question 5 and Assignment 3 Question 5 compute critical radii for a genuine
`k = 1` condition, so a comparison against them is a comparison of two different definitions,
not a validation.

## 4. No variance reduction, and the subcritical side is the expensive one

Every history is analogue: no splitting, no Russian roulette, no weights. `check_criticality`
exits on the first critical history, so a clearly supercritical radius costs one history while
a subcritical one costs the full 100 — which is the wrong way round for the bisection, whose
later steps sit close to the transition.

## 5. Nothing is checked against a known answer

There is no verification of any kind in this assignment — no analytic limit, no benchmark, no
internal identity, unlike the other three assignments. Two cheap checks would cover most of the
sampling code in [02](02-the-random-walk.md):

- **The flight distribution.** Sample `s` a few hundred thousand times and compare the mean
  against `mfp` and the histogram against `exp(-s/mfp)/mfp`.
- **Direction isotropy.** Sample directions and check that each component's mean is ~0 and that
  `cos_theta` is uniform on `[-1,1]` — this is the one that would catch the pole-clustering
  error, and it would catch the shared-buffer failure mode described in
  [02](02-the-random-walk.md) as well.

A third, for the physics: with `P_2 = 0` the material cannot multiply, so no radius should ever
be critical. `check_criticality` should return `False` for arbitrarily large `R`. That is a
one-line test and it exercises the whole chain.
