# 04 — Why Shooting, and What the Finite-Volume Solvers Measured

**Finite-volume solvers were written, measured and then removed: shooting is `2.0e-8 %`
accurate against their `1.5e-2 %`. Their numbers are kept here because the comparison is
what answers the assignment's question about modelling a delta source.**

Three solvers existed at one point for the Question 2 slab problem. Only
`solve_diffusion_shooting` survives; the others are recoverable from the history of branch
`homework1-q2-diffusion`.

| solver | max relative error | order |
|---|---|---|
| shooting, half-domain, radiation BC | `2.0e-8 %` | tolerance-limited, see [03](03-delta-source-and-boundary.md) |
| finite volume, half-domain, 500 cells | `1.49e-2 %` | `2.000` |
| finite volume, smeared source, 501 cells | `3.95e-2 %` | `2.000` |

## The two delta-source treatments

The half-domain solver converts the delta into the boundary current `J(0+) = 1/2`; the
alternative keeps the full domain `[-a, a]` and smears the source over the cell containing
the origin, `S_i = 1/dx`. Both converge at second order, but the smeared source carries a
`~14x` larger error constant at equal `dx` (`6.49e-2` against `4.72e-3` relative `L2` on the
coarsest mesh), because it spreads the kink at the origin over a cell instead of placing it
exactly on a cell boundary. That measurement is the empirical half of the answer to "how
should a delta-function source be modelled?" — the symmetry argument in
[03](03-delta-source-and-boundary.md) is the analytic half.

## Two false trails, recorded so they are not re-run

**The `brentq` tolerance was not the limit.** The original shooting solver root-found on the
starting slope, and its `xtol = 2e-12` looked coarse against a root of `s = 1.362e-4`. It was
measured: `brentq` recovers the root to every printed digit. The accuracy was limited by the
zero-flux outer boundary alone — and the equation is linear, so the root-finder was
unnecessary in the first place.

**A `1 %` balance deficit was the check, not the solver.** `absorption_balance` initially
read `0.99005` on the finite-volume output. A trapezoid rule over *cell-averaged* data omits
the half-cell slivers at each end, discarding a fraction `kappa dx / 2` of the integral —
exactly `1 %` at that resolution. The shooting solver returns point values on a grid, so the
rule is correct there and `absorption_balance` needs no quadrature argument.

## One boundary-value problem, not six

The domain is scaled as `a = 10/kappa`, so in units of `kappa x` every `c` and both
approximations are the *same* problem. That is why the error and balance figures are
identical across all six cases: `c` and the approximation enter only through the scaling.
