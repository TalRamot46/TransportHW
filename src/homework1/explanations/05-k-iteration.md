# 05 — The `k` Iteration

**Two lines in `k_eigenvalue` look like housekeeping and are not: dividing the source by `k`,
and rescaling the flux every sweep.**

Report §4 derives the algorithm and the `c/(4c-3)` convergence rate. This is what the code does
with them.

## The loop

`k_eigenvalue(R, medium, n_cells)`: assemble the operator once, then per sweep —

1. `phi_new = solve_banded(ab, fission / k)` — one `O(N)` banded solve;
2. `k_new = k * sum(fission_new) / sum(fission)`;
3. rescale `phi, fission` by `1/max(phi_new)`;
4. stop when `|k_new - k| < 1e-10 |k_new|`.

**Dividing by `k`** (step 1) makes `M phi = F/k` an ordinary inhomogeneous problem with `M`
positive definite, solvable at *any* radius rather than only the critical one. That is what
turns criticality into `brentq` on `k(R) - 1` in `critical_radius`, and it is why the singular
eigenvalue problem is never confronted directly.

**Rescaling** (step 3) is not cosmetic. Each sweep multiplies the amplitude by roughly `k`, and
this problem needs hundreds of sweeps; without the rescale the iterate drifts by orders of
magnitude and `|k_new - k|` becomes a difference of numbers that have lost their significance.
The ratio in step 2 is homogeneous in `phi`, so the rescale changes nothing physically — note
that it is applied to `fission` as well, or the two would fall out of step.

The convergence test is evaluated **before** the rescale is applied but stored after; that
ordering is why `converged` is computed on its own line rather than inside the `if`.

## Why the sweep cap is 20000

As `c -> 1` the dominance ratio tends to 1 and convergence becomes arbitrarily slow — 280
sweeps already at `c = 1.02`, and radii away from critical are worse. Each sweep is one `O(N)`
banded solve, so a high cap costs nothing, and `_bracket` keeps `k` from ever being evaluated
far from critical. The whole of Question 4 runs in about a second, which is why none of the
acceleration schemes Bell & Glasstone discuss are implemented.

## The bracket starts narrow on purpose

`_bracket` widens by 5% per step around the analytic radius. The narrowness is the point: an
evaluation of `k` far from criticality is the expensive kind, because the dominance ratio is
worse there. Starting wide would cost more than the search saves.

`critical_radius` then calls `brentq` with `xtol=1e-12` on a quantity of order 1–12 mean free
paths, which is at the useful floor given the iteration converges to `1e-10` relative.
