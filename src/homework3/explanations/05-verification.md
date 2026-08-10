# 05 — Verification

**The six checks in `main.py` and the numbers they produce.**

Run with `python -m homework3.main` from the repository root. Every number below is measured,
not asserted.

## 1. Particle balance, `int phi(x,t;c) dx = e^{-(1-c)t}`

With the exact series `G`, over all five `c` and `t = 1, 4, 15`, the relative error is
`4e-15` to `1e-9` (the larger values are the quadrature struggling with a flux spread over
`|x| < 15`, not the formula).

This is the strong check. It pins down three things at once that no other test touches: the
overall prefactor `sqrt(3/pi)/(2vt)`, the cancellation of `Sigma_t`, and the Jacobian factor
`c` in the Q1 scaling — get any of them wrong and the number is off by a clean constant.

## 2. Cost of the interpolated `G`

The same integral with Paasschens' interpolation instead of the exact series:

| `t` | 0.3 | 1 | 3 | 10 | 30 |
|---|---|---|---|---|---|
| `int phi dx` | 1.00076 | 1.00537 | 1.01470 | 1.01092 | 1.00429 |

So the interpolation costs at most `+1.5 %` in particle number, peaking near `t = 3` and
falling away in both directions. This is Paasschens' own quoted accuracy for his Eq. (36), and
it is the *modelling* error of the solution the assignment specifies — not an error in this
code. It bounds how much of any discrepancy with diffusion can be blamed on `G`.

## 3. Value at the causal front

As `|x| -> vt`, `w0 -> 0`, and the collided bracket cancels to zero, leaving the uncollided
plateau `e^{-t}/(2t)` alone:

| `t` | `phi(|x| -> vt)` | `e^{-t}/(2t)` |
|---|---|---|
| 1 | 0.1839398 | 0.1839397 |
| 4 | 0.00228946 | 0.00228945 |
| 15 | 1.02e-8 | 1.02e-8 |

Physically obvious — a particle arriving exactly at the front cannot have scattered, or it
would have fallen behind — and therefore a good structural test of the bracket.

## 4. `D0(c)` continuity through `c = 1`

`D0(1 - 1e-6) = 0.33333360` and `D0(1 + 1e-6) = 0.33333307`, both against `1/3 = 0.33333333`.
The two branches of the eigenvalue equation meet, so the piecewise definition in
[03](03-diffusion-coefficients.md) is a limit rather than a patch.

## 5. The steady identity

`int_0^inf phi_diff(x,t;c) dt` against the closed-form steady solution, for `c = 0.6, 0.8`,
`x = 0.5, 2.0`, both approximations: agreement to `2e-16`–`7e-14` relative.

This is the check that connects the time-dependent Green's functions the figures use to the
steady formula the assignment quotes (see [04](04-time-dependent-vs-steady.md)). Evaluated at
`x` very close to 0 the quadrature degrades to `~1e-4`, because the integrand acquires an
integrable `t^{-1/2}` singularity at `t -> 0`; that is `scipy.integrate.quad`, not the
formulas.

## 6. Diffusion limit at late times

`phi(0,t;1)` against the classical diffusion peak `(4 pi t/3)^{-1/2}`:

| `t` | exact | diffusion | error |
|---|---|---|---|
| 10 | 0.162054 | 0.154510 | `+4.9 %` |
| 100 | 0.049110 | 0.048860 | `+0.5 %` |
| 300 | 0.028258 | 0.028209 | `+0.2 %` |

The exact solution relaxes onto the diffusion solution, and the error falls roughly as `1/t` —
the expected rate, since the leading transport correction to diffusion is one collision time
out of `t`. This is the end state the `t = 15` panels of the figures are heading towards, and
it validates the two solutions *against each other* rather than each against itself.
