# 05 — Verification

**Twelve checks in `main.py`, run on every invocation. There is no separate test file, so
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

## 12. Diffusion against exact transport — `check_diffusion_error`

The part 3(d) measurement, and the only check whose output the report quotes as a table rather
than as a bound. Signed relative error of each numeric diffusion flux against `phi_exact`, at
the origin and at `x = 0.9 vt`; the report's Table 2 is the `x = 0` half of it.

Two choices inside it are worth recording:

- **`form="series"`, not the default interpolation.** Check 2 puts the interpolation's own cost
  at up to `1.5%`, and the asymptotic error at `c = 0.8, t = 15` is `-0.59%` — the reference
  would otherwise be less accurate than the quantity being measured.
- **`0.9 vt`, not the front itself.** `phi_exact` is identically zero beyond `vt`, so the ratio
  has no limit there; `0.9` is far enough out to show the failure (`+144415%` at `c = 1.5`,
  `t = 15`) without dividing by zero.

The five `c` values are solved twice each, once per approximation, and each solve marches
through all of `SOLVER_TIMES` in one pass — ten marches, not thirty. That matters: this check
and `plot_diffusion_error` together are most of the runtime of `main.py`.
