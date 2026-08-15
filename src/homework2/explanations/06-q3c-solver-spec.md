# 06 — The Q3(c) Solver — Specification

**Nothing in this file is built yet.** Parts 3(a) and 3(b) are answered by closed forms, so
`src/homework2` currently contains no solver. Part 3(c) asks for a time-dependent numerical
diffusion code, and this is the specification for it — enough that writing `solver.py` is
mechanical.

## What it solves

Remove the absorption analytically first, leaving the plain heat equation with the pulse still
in it:

    du/dt = D v d2u/dx2 + v delta(x) delta(t),      phi(x,t) = e^{-(1-c) Sigma_t v t} u(x,t)

The `c`-dependence then enters only through `D` and through the prefactor reapplied on the way
out, so one solver covers both approximations and every `c`.

## Module layout

`src/homework2/solver.py`, five private helpers and one public entry point:

```python
def _grid(D, t_max, x_max, n_nodes)  -> x, h    # node-centred half domain [0, L]
def _laplacian(n_nodes, h, D)        -> A       # tridiagonal, BCs baked into two rows
def _pulse_source(x, h)              -> u0      # the delta, in node 0
def _warm_start(x, t0, D)            -> u0      # analytic Gaussian at t0 (diagnostic only)
def _step_matrices(A, dt)            -> banded (I - dt/2 A), (I + dt/2 A)
def solve_diffusion(c, approximation, times, start="pulse", ...) -> (x, {t: phi})
```

`solve_diffusion` takes `D` from `diffusion.diffusion_coefficient`, so it inherits the `c = 1`
handling for free ([03](03-diffusion-module.md)). The `start` switch is what makes the
convergence diagnostic possible: `"pulse"` is the answer to 3(c), `"warm"` starts from the
analytic Gaussian at a small `t0` and isolates the time stepping from the source treatment.

**Grid.** Node-centred, `x_j = j h` on the half domain `[0, L]` — all three solutions are even.
Node-centred rather than cell-centred so the symmetry condition sits exactly on a node and the
output is point values, which is what the plots and a trapezoid balance check want. `L` is set
by the far tail rather than guessed: `L = x_max + 8 sqrt(2 D v t_max)`, which at `t_max = 15`
and the largest `D` in play gives `L ~ 50`, where the Gaussian is `~1e-38` of its peak. With
`N = 4000`, `h = 0.0125`.

## The delta: smear it into cell 0, and know what that costs

Keep the standard treatment, `S_0 = 1/h`. The point is not that it is acceptable but that its
error can be written down exactly. A top-hat of width `h` has variance `h^2/12`, the diffusion
Green's function has variance `2Dvt`, and variances add under convolution, so the smeared
solution has second moment

    2 D v t + h^2/12  =  2 D v ( t + h^2/(24 D v) )

**The smeared source reproduces the exact solution at a shifted time**, `dt_eff = h^2/(24Dv)`.
That is more useful than "the error is `O(h^2)`" because it says where the error lives —
entirely in early time, decaying in relative terms as `dt_eff/t`:

| `h` | `dt_eff` | relative error at `t = 1` |
|---|---|---|
| 0.5 | 3.1e-2 | ~3% |
| 0.1 | 1.3e-3 | ~0.1% |
| 0.0125 | 2.0e-5 | ~2e-3% |

So on the mesh above the delta contributes essentially nothing, and it is smallest exactly
where diffusion is most accurate. **The time stepping, not the source, is the real problem.**

## Time: Crank–Nicolson, hand-written

`solve_ivp` is not an alternative to discretising time — BDF *is* a time discretisation. The
real distinction is who picks the step size and order. Discretise space first to get
`du/dt = A u`, then with `r = D v dt / h^2`:

| scheme | `g(theta)` | stable when | order |
|---|---|---|---|
| FTCS | `1 - 4r sin^2(theta/2)` | `r <= 1/2` | `O(dt) + O(h^2)` |
| BE | `1/(1 + 4r sin^2(theta/2))` | all `r` | `O(dt) + O(h^2)` |
| CN | `(1 - 2r sin^2)/(1 + 2r sin^2)` | all `r` | `O(dt^2) + O(h^2)` |

Three conclusions, each of which decides something:

- **FTCS is limited by stability, not accuracy.** Violating `r <= 1/2` blows up rather than
  degrading, and `dt <= h^2/(2Dv)` means halving `h` quarters `dt` — cost grows as `h^-3`.
- **Unconditional stability alone buys nothing.** BE is stable at any `dt` but only first order
  in time, so matching an `O(h^2)` spatial error still needs `dt ~ h^2`. BE alone is strictly
  worse than FTCS here — same step count, higher cost per step.
- **CN is the one that pays.** Second order in both, unconditionally stable, so `dt` is set by
  balancing `O(dt^2)` against `O(h^2)`, i.e. `dt ~ h` instead of `dt ~ h^2`.

At `h = 0.0125`, `D = 1/3`, to `t = 15`: FTCS and BE both need ~64,000 steps; CN needs a few
hundred. Every scheme is tridiagonal and the implicit ones cost one `solve_banded` per step,
which is `O(N)` — the same order as the explicit matrix–vector product.

**Startup.** CN's amplification factor tends to `-1` at `theta = pi`, so a discontinuous
initial condition — which a first-cell delta is — rings. Take the first one or two steps with
backward Euler, whose factor tends to `0` and damps those modes, then switch to CN. This is the
standard Rannacher startup and it costs nothing in order.
