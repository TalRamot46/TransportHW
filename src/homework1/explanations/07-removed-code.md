# 07 — Removed Code and False Trails

**Three solvers and one boundary option were written, measured and deleted. Their numbers are
kept because they are the evidence for the choices the surviving code makes.**

## Three Q2 solvers, one survivor

Only `solve_diffusion_shooting` remains; the others are recoverable from branch
`homework1-q2-diffusion`.

| solver | max relative error | order |
|---|---|---|
| shooting, half-domain, radiation BC | `2.0e-8 %` | tolerance-limited |
| finite volume, half-domain, 500 cells | `1.49e-2 %` | `2.000` |
| finite volume, smeared source, 501 cells | `3.95e-2 %` | `2.000` |

The two finite-volume variants differ in how they treat the delta: the half-domain one converts
it to `J(0+) = 1/2`, the other keeps `[-a, a]` and smears it over the origin cell as
`S_i = 1/dx`. Both are second order, but the smeared source carries a **~14× larger error
constant** at equal `dx` (`6.49e-2` against `4.72e-3` relative `L2` on the coarsest mesh),
because it spreads the kink at the origin over a cell instead of putting it exactly on a cell
boundary. That measurement is the empirical half of the assignment's "how do you model a delta
source" question; report §2 is the analytic half.

## The `'robin'` boundary option

`build_medium` and `_banded_operator` once took a `boundary` argument selecting between an
extrapolated zero and a Robin condition at the physical surface. Both were implemented and
measured; neither could be shown more accurate than the other against any reference available
here, so the option was deleted in favour of one unambiguous condition. **`k_eigenvalue` now
always uses the extrapolated zero**, and the Robin numbers below are no longer reproducible
from the current code.

The whole distinction had been one denominator, `l0 + 0.5*h`:

| treatment | mesh covers | `l0` | condition realised |
|---|---|---|---|
| extrapolated (kept) | `[0, R + z0]` | `0` | `phi(R + z0) = 0` |
| robin (removed) | `[0, R]` | `z0` | `phi(R) + z0 phi'(R) = 0` |

so `l0 = 0` was never a special case — it degenerates to the zero-flux form.

It is worth recording *why* the spread mattered, since it is a real property of the model and
not solver error: each treatment reproduced its own analytic counterpart to `1.5e-6`, yet they
disagreed with each other by 0.26% at `c = 1.02` rising to 4.5% at `c = 2`, which reached the
critical masses as an **8–12% spread**. Extrapolated is the default because the analytic
relation everything is checked against, `R_c = pi/B - z0`, *is* the extrapolated statement — so
the comparison tests the numerics rather than the boundary model.

## Two dead ends, recorded so they are not re-run

**The `brentq` tolerance was never the limit.** The original shooting solver root-found on the
starting slope, and its `xtol = 2e-12` looked coarse against a root of `s = 1.362e-4`. It was
measured: `brentq` recovered the root to every printed digit. The accuracy was limited by the
zero-flux outer boundary alone — and since the equation is linear, the root-finder was
unnecessary in the first place, which is how the current one-pass-and-rescale form arrived.

**A 1% balance deficit was the check, not the solver.** `absorption_balance` once read
`0.99005` on finite-volume output. A trapezoid rule over *cell-averaged* data omits the
half-cell slivers at each end, discarding a fraction `kappa dx / 2` of the integral — exactly
1% at that resolution. The solver was right and the diagnostic was wrong. The shooting solver
returns point values, so the rule is correct there and `absorption_balance` needs no quadrature
argument at all ([03](03-shooting-solver.md)).
