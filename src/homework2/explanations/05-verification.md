# 05 — Verification

**Six checks in `main.py`, run on every invocation. There is no separate test file, so these
are the whole safety net — the numbers below are from the current code.**

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
