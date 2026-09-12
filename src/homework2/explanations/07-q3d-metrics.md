# 07 — The Q3(d) Metrics

**`metrics.py` turns a profile comparison into three numbers; the only subtle one is that
`total_variation` integrates on its own grid, not on the solver's.**

Report §3, part 3(d) derives all three — eq. (46) for the misplaced fraction, eq. (49) for the
leak past the front, eq. (52) for the floor. What follows is only what the code does about them.

## The grid `total_variation` integrates on

`phi_exact` has a jump at `|x| = vt`, from the uncollided plateau `e^{-t}/(2t)` straight to
zero. The solver's grid is sized for the diffusion tail, not for that jump, so at `t = 1` in a
run that marches to `t = 15` only about 45 nodes sit inside the front and the one cell
straddling the discontinuity carries an `O(h)` error of its own.

So `total_variation` interpolates the numeric profile onto `N_FINE = 8001` points and evaluates
`phi_exact` there analytically. The interpolation costs nothing — the numeric profile is a
smooth Gaussian — and the jump error drops with the finer `h`. What is left is visible as the
residual node-count dependence, worst at the earliest time:

| `N` | 400 | 600 | 1000 | 2000 |
|---|---|---|---|---|
| `delta` at `t = 1`, `c = 0.6` | 22.77% | 22.40% | 22.20% | 22.12% |
| `delta` at `t = 15` | 10.065% | 10.095% | 10.110% | 10.116% |

The report quotes two decimals at `N = 2000`; the `t = 1` column is the one that would move
first, by about `0.1` in the second decimal.

## `gaussian_gap` is where the `t` cancels

`late_time_floor` is not measured, it is evaluated: the floor is the gap between two Gaussians
whose variances hold a fixed ratio, and report eq. (52) shows every `t` cancelling out of it.
The code therefore takes two `D` values and no time at all. `sorted()` on the pair is what lets
one expression serve both `D > D_inf` and `D < D_inf`, which is the difference between `c < 1`
and `c > 1` — the ratio `a^2` must be the larger over the smaller or the `erf` arguments go
negative under the root.

Checked against a direct 400001-point `L1` integral of the two closed-form Gaussians, the
formula agrees to seven digits at `t = 2` and at `t = 50` — the same value, which is the point.

## Two leakage functions, deliberately

`front_leakage` is the closed form `erfc(sqrt(vt/4D))`; `measured_front_leakage` is the mass the
solver actually has beyond `x = vt`. They exist as a pair so the check can compare them, and
they agree to `3%` (`2.207e-1` against `2.147e-1` at `t = 1`, `2.101e-6` against `2.032e-6` at
`t = 15`).

**That `3%` is not a bug.** It is the `O(h^2)` truncation measured in the worst place for it —
a tail six decades below the peak, past the point where `bulk()` stops calling a relative error
meaningful. Check 9 measures the scheme at `6.5e-5` *inside three standard deviations*; nothing
claims that accuracy out here, and the leak is still resolved to two digits.

## Why `1/(3c)` is not in `diffusion.py`

`late_time_coefficient` returns `D_CLASSICAL / c` and sits here rather than beside the other two
coefficients, because it is not an approximation anyone solves with — it is the coefficient the
*exact* flux spreads with, derived from the Question 1 identity, and it exists only to be the
reference that `late_time_floor` measures against. Putting it in `diffusion.py` would invite
`diffusion_coefficient(c, 'exact')`, which is not a thing.
