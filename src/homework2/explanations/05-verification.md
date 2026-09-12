# 05 — Verification

**Fourteen checks in `main.py`, run on every invocation. There is no separate test file, so
these are the whole safety net — the numbers below are from the current code.**

## 1. Normalisation — `check_normalisation`

`int phi(x,t;c) dx` must equal `e^{-(1-c)t}`: one particle emitted, decaying at the absorption
rate. Run with the **exact series** `G`, so it tests the flux and not the interpolation.

Relative error is `4e-15` to `9e-15` at `t = 1` and `t = 4` for every `c` in `0.6 … 1.5`,
loosening to `1.3e-9` at `c = 0.6, t = 15`. That worst case is the quadrature, not the flux:
at `t = 15` the integrand spans `e^{-15}` over a domain of half-width 15.

This is the strongest check in the file. It exercises `G`, `collided_integral`, `phi_c1` and
the Q1 scaling of `phi_exact` at once, on both sides of `c = 1`, and it would fail on any
error in the constant factors that the shape-based checks would miss.

## 2. Cost of the interpolation — `check_interpolation_error`

The same integral with the default interpolated `G`, whose departure from 1 is the price of
Paasschens' `G(w) ~ e^w sqrt(1+b/w)`:

| `t` | 0.3 | 1.0 | 3.0 | 10.0 | 30.0 |
|---|---|---|---|---|---|
| departure | +0.076% | +0.537% | +1.470% | +1.092% | +0.429% |

It peaks near `t = 3` at **1.5%** and decays either side, comfortably inside the ~2% Paasschens
quotes. That is the accuracy of every figure in the report, since the plots use the default.

## 3. The causal front — `check_front`

As `|x| -> vt` the collided bracket must vanish, leaving the uncollided plateau `e^{-t}/(2t)`.
At `t = 1`: `0.1839397971` against `0.1839397206`. At `t = 15` both are `1.02e-8`. This is the
check that the three front guards in `phi_c1` ([02](02-evaluating-the-closed-form.md)) agree
with each other.

## 4. The diffusion coefficient — `check_diffusion_coefficients`

`D0(c) = (1-c) nu0^2` must stay positive on both sides of `c = 1` and tend to `1/3` there:

| `c` | 0.6 | 0.8 | 1.0 | 1.2 | 1.5 |
|---|---|---|---|---|---|
| `D0` | 0.48588 | 0.39629 | 0.33333 | 0.28717 | 0.23745 |

and at `c = 1 ± 1e-6` it gives `0.33333360` and `0.33333307` against `1/3 = 0.33333333` — the
two-sided limit, approached from opposite directions. That is what justifies the short circuit
at exactly `c = 1` ([03](03-diffusion-module.md)).

## 5. The steady identity — `check_steady_identity`

`int_0^inf phi_diff(x,t) dt` must equal the closed-form steady solution
`e^{-kappa|x|}/(2 D kappa)`. Over `c = 0.6, 0.8` and `x = 0.5, 2.0`, both approximations, the
relative error runs from `2e-16` to `7e-14`.

This ties the time-dependent Green's function to the steady one Assignment 1 used, so the two
assignments are demonstrably solving the same equation. It is restricted to `c < 1` because
`_phi_steady` raises above it.

## 6. The diffusion limit — `check_diffusion_limit`

At late times the exact transport solution must relax onto the classical diffusion peak
`(4 pi t/3)^{-1/2}`:

| `t` | 10 | 100 | 300 |
|---|---|---|---|
| departure | +4.88% | +0.51% | +0.17% |

Falling roughly as `1/t`, which is the expected approach. This is the only check that tests
`exact.py` and `diffusion.py` **against each other** rather than each against its own identity,
and it is the quantitative form of the report's central claim — that diffusion is the
late-time limit of transport, and the `t = 1 … 15` figures sit in the range where it is not yet
reached.

## 7. Solver conservation — `check_solver_conservation`

`2 * trapezoid(u, x)` against `e^{-(1-c)t}`, over `c = 0.6, 1.0, 1.5` and `t = 1, 4, 15`. The
relative error is between `0` and `1.1e-15` in all nine cases — round-off, not discretisation.

That is not an accuracy result but a structural one: report eq. (45) says the only leak is
`-r h u[N-1]` at the far boundary, and the truncation is placed where `u[N-1]` is denormal. A
wrong factor in the symmetry row would break this immediately while leaving every plot looking
right, which is what the check is for.

## 8. Solver against the closed form — `check_solver_against_closed_form`

Max relative error within three standard deviations:

| `t` | 1 | 4 | 15 |
|---|---|---|---|
| error | 1.52e-3 | 3.81e-4 | 1.02e-4 |

identical for every `c`, which is the expected answer rather than a suspicious one: with the
classical `D = 1/3` the grid does not depend on `c` at all, and the absorption prefactor
cancels in a *relative* error.

The error falls as `1/t` — the signature of a fixed offset in the variance of the computed
Gaussian, whose relative effect decays as the true variance `2Dvt` grows. Check 11 shows that
offset is not the smeared source, so what is left is the scheme's own truncation.

## 9. Solver order — `check_solver_order`

Bulk error at `t = 4` against the node count, at fixed `r`:

| `N` | 500 | 1000 | 2000 |
|---|---|---|---|
| error | 1.007e-3 | 2.586e-4 | 6.466e-5 |
| ratio | — | 3.89 | 4.00 |

Second order, as the `O(dt) + O(h^2)` truncation with `dt ~ h^2` requires. This is the check
that would fail if the `j = 0` row were only first-order accurate — a mistake that leaves
conservation intact and so slips past check 7.

## 10. The stability edge — `check_solver_stability_edge`

400 sweeps of an isolated spike on 101 nodes, either side of `r = 1/2`:

| `r` | 0.49 | 0.51 |
|---|---|---|
| peak amplitude | 2.014e-2 | 1.310e5 |

Seven orders of magnitude apart across a 4% change in `r`. Report eq. (41) predicts exactly
this: past the limit the `theta = pi` mode grows as `|1-4r|^n`, so the scheme does not lose
accuracy, it explodes.

## 11. The source treatment — `check_solver_source_treatment`

Starting from the smeared delta versus from the analytic Gaussian at `WARM_T0`, compared at
`t = 4`: the two agree to `8.08e-7`, against a discretisation error of `6.47e-5` at the same
point.

**This is the empirical form of report eq. (43).** The source treatment sits nearly two orders
of magnitude below what limits the solver, so smoothing the initial spike — into a half-Gaussian
or anything else — would buy nothing. It also rules out the opposite worry: that the first cell
being the entire source leaves a defect the march never recovers from.

## 12. The misplaced fraction — `check_total_variation`

The part 3(d) headline, and the table the report quotes. Half the `L1` distance between the
numeric diffusion density and the exact one, for both approximations over `C_VALUES` and
`METRIC_TIMES`; [07](07-q3d-metrics.md) covers the grid it integrates on.

The one result worth pulling out is that the classical column is **not** monotone in `t` for
`c != 1`. At `c = 0.6` it reads `22.09, 8.32, 5.59, 7.71, 10.12` percent — it turns at about
`t = 4` and climbs. Check 13 is what explains that, and it is the reason this check exists
rather than a single late-time number.

## 13. The late-time floor — `check_late_time_floor`

Measured `delta` at `t = 40` against `late_time_floor`, which is evaluated rather than measured:

| `c` | classical, measured / floor | asymptotic, measured / floor |
|---|---|---|
| 0.6 | 11.48% / 12.29% | 2.41% / 3.24% |
| 1.0 | 0.51% / 0 | 0.51% / 0 |
| 1.5 | 10.12% / 9.78% | 1.94% / 1.60% |

The gap between measured and floor is the transport correction that has not died yet, and the
`c = 1` row sizes it: `0.51%` at `t = 40`, falling as `t^-1/2`. Read the other rows with that in
mind — `11.48` against `12.29` is approaching from below, `10.12` against `9.78` from above, and
both are inside that correction.

**This check is the one that would catch a wrong `D_inf`.** Setting it to `1/3` instead of
`1/(3c)` makes every classical floor zero, and the measured `c = 0.6` and `c = 1.5` columns
would then have nothing to converge to.

## 14. The leak past the front — `check_front_leakage`

`erfc(sqrt(vt/4D))` against the mass the solver puts beyond `x = vt`, with the front's distance
in standard deviations alongside:

| `t` | 1 | 4 | 15 |
|---|---|---|---|
| front, in sigmas | 1.2247 | 2.4495 | 4.7434 |
| closed form | 2.207e-1 | 1.431e-2 | 2.101e-6 |
| solver | 2.147e-1 | 1.410e-2 | 2.032e-6 |

The closed form is `c`-independent for the classical `D`, which is why the three `c` blocks
print identical numbers — that repetition is the check, not an oversight. The `3%` gap to the
solver is discussed in [07](07-q3d-metrics.md).
