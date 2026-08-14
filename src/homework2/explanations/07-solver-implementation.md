# 07 — Implementing the numerical solver

**A concrete specification for the Question 3(c) time-dependent diffusion code.**

The scheme and its justification are in [06](06-delta-source-numerics.md); this file is the
implementation detail — enough that writing `solver.py` is mechanical. Nothing here is built
yet.

## What is being solved

Only the absorption is removed analytically (Step 2 of the derivation, report §3.2.2), leaving
the solver the **heat equation** with the pulse still in it:

```
du/dt = D v d2u/dx2 + v delta(x) delta(t),      phi(x,t) = e^{-(1-c) Sigma_t v t} u(x,t)
```

The `c`-dependence now enters only through `D` and through the prefactor reapplied at the end.
The delta stays, and is discretised — that is the part Question 3(c) is actually asking about.

## Module layout

`src/homework2/solver.py`, five short functions:

```python
def _grid(D, t_max, x_max, n_nodes)       -> x, h        # node-centred half domain [0, L]
def _laplacian(n_nodes, h, D)             -> A           # sparse tridiagonal, BCs baked in
def _pulse_source(x, h)                   -> u0          # the delta, in node 0
def _warm_start(x, t0, D)                 -> u0          # analytic Gaussian at t0 (diagnostic)
def _step_matrices(A, dt)                 -> banded (I - dt/2 A), (I + dt/2 A)
def solve_diffusion(c, approximation, times, start="pulse", ...) -> (x, {t: phi})
```

`solve_diffusion` is the only public one. It picks `D = diffusion_coefficient(c, approximation)`
from `diffusion.py`, builds the grid and operator, sets the initial state from `start`,
integrates, and reapplies the absorption factor on the way out. The `start` switch is what makes
verification item 5 possible: `"pulse"` is the answer to 3(c), `"warm"` is the diagnostic run.

## Grid

Node-centred, `x_j = j h`, `j = 0..N`, `h = L/N`, on the half domain `[0, L]` — all three
solutions are even in `x`. Node-centred rather than cell-centred for two reasons: the symmetry
condition sits exactly on a node, and the output is point values, which is what the plots and
the trapezoid balance check want (see Pitfalls).

`L` is set by the far tail, not guessed. The widest Gaussian in the run has
`sigma_max = sqrt(2 D v t_max)`; taking eight of them past the plotted range,

```
L = x_max + 8 sqrt(2 D v t_max)
```

At `t_max = 15`, `x_max = 1.25 t_max = 18.75` and the largest `D` in play (`c = 0.6`
asymptotic, `D0 = 0.4859`), that is `L ~ 50`, where the Gaussian is `~1e-38` of its peak. With
`N = 4000`, `h = 0.0125`.

## The operator and its two boundary rows

Centred second difference, `(u_{j-1} - 2u_j + u_{j+1}) / h^2`, times `Dv`. Only the two end
rows need thought; both come from eliminating a ghost node.

**At `j = 0` — symmetry.** The solution is even, so `du/dx(0,t) = 0`. The ghost node is
`u_{-1} = u_1`, and the row becomes

```
(2 u_1 - 2 u_0) / h^2
```

This is a condition on the *shape*, not on the source; the pulse enters through the initial
state below, not through this row.

**At `j = N` — zero flux, far away.** `du/dx(L,t) = 0`, ghost `u_{N+1} = u_{N-1}`, row

```
(2 u_{N-1} - 2 u_N) / h^2
```

> **Why not the radiation condition here.** Assignment 1 replaced `phi(L) = 0` with
> `phi' = -kappa phi`, `kappa = sqrt(Sigma_a/D)`, and that fixed a 100 % boundary artifact. It
> does not transfer to this problem. That condition is exact because the steady solution *is*
> `e^{-kappa x}`, so imposing it annihilates the reflected mode. Here the solution is a
> spreading Gaussian, not an exponential, and `kappa` is not its local decay rate — worse,
> having factored the absorption out, `u` has no `kappa` at all, and for `c >= 1` the original
> `Sigma_a <= 0` makes `kappa` imaginary anyway. The honest choice is a reflecting boundary
> placed where nothing has arrived: the error is then bounded by the Gaussian tail at `L`,
> which the sizing above drives to `1e-38`. This is a domain-truncation error made negligible
> by construction rather than by a clever condition.

`A` is tridiagonal — build it with `scipy.sparse.diags` and keep it sparse.

## The initial state

### `start="pulse"` — the answer to 3(c)

Integrating across `t = 0` turns the pulse into initial data, `u(x, 0^+) = v delta(x)`, which
on the mesh becomes a spike in node 0. On the **half** domain the node at the origin represents
a half cell of width `h/2`, and it carries half the particles, so for a unit total source

```
u_0 = v / h,      u_j = 0  for j > 0
```

(The two factors of two cancel: half the source, spread over half a cell.) Check it by the same
trapezoid rule used for the balance test — `2 * trapezoid(u, x)` must come out `v`, not `2v` or
`v/2`. Getting this wrong is the single most likely bug in the whole solver, and it shows up as
a clean factor of two everywhere, at every time.

Integration then starts at `t = 0`.

### `start="warm"` — the diagnostic

Start instead at a small `t0` from the analytic Gaussian, chosen so its width spans several
cells. Five gives

```
t0 = (5h)^2 / (2 D v),    u(x, t0) = v / sqrt(4 pi D v t0) * exp(-x^2 / (4 D v t0))
```

which at `h = 0.0125`, `D = 1/3` is `t0 ~ 5.9e-3` — two orders of magnitude before the first
output time. Note this is the analytic solution **for `u`**, without the absorption factor.

Record `t0` in the returned metadata: it is a parameter of the method, and the convergence study
has to hold it fixed while `h` varies, or the two effects mix.

## Time integration

The scheme comparison and the reasoning behind this choice are in
[08](08-time-discretisation.md). What it settles: write **Crank–Nicolson with a backward-Euler
startup**, fixed step, hand-discretised.

### The scheme

With `r = D v dt / h^2` and the same `A` as above, each step is one banded solve:

```
Rannacher startup, steps 1..4:   (I - (dt/4) A) u^{n+1} = u^n
Crank-Nicolson, steps 5.. :      (I - (dt/2) A) u^{n+1} = (I + (dt/2) A) u^n
```

Written out, the CN interior row is

```
-r/2 u_{j-1}^{n+1} + (1+r) u_j^{n+1} - r/2 u_{j+1}^{n+1}
  =  r/2 u_{j-1}^n + (1-r) u_j^n + r/2 u_{j+1}^n
```

and the two end rows carry the same ghost-node substitutions as `A` — at `j = 0`,
`-r u_1^{n+1} + (1+r) u_0^{n+1} = r u_1^n + (1-r) u_0^n`, and likewise at `j = N` with
`u_{N-1}`. Build the three bands once and call `scipy.linalg.solve_banded` each step; the matrix
is constant, so `scipy.linalg.lu_factor` on the banded form once is better still.

**The startup is not optional here.** CN rings on a delta initial condition — the sign flip in
its amplification factor sets in at exactly `r > 1/2`, and the delta excites every mode equally.
Measured at `h = 0.0125`, `dt = h`: ten CN steps without the startup give a **negative** flux of
`-13.4` near the origin and a relative error of `2.2e+1`; with four backward-Euler quarter-steps
first, `+6.8e-2` and `1.3e-3`. See [08](08-time-discretisation.md).

Choose `dt ~ h` (not `h^2`), which balances the `O(dt^2)` and `O(h^2)` errors: ~1,200 steps to
`t = 15` at `h = 0.0125`, against the ~64,000 an explicit march would be forced into.

### Two library cross-checks, optional

**`solve_ivp(..., method="BDF", jac=lambda t, u: A, t_eval=times, rtol=1e-10)`.** Adaptive in
both step and order. Agreement with the hand-written CN is strong evidence that the *spatial*
operator and its boundary rows are right, since that is the only thing the two share. Pass `jac`
explicitly — without it BDF finite-differences a 4000x4000 Jacobian — and scale `atol` to the
initial peak rather than leaving it at `1e-6`.

**`scipy.sparse.linalg.expm_multiply(A * (t - t0), u0)`.** The *exact* solution of the
semi-discrete system — zero time-integration error by construction, so any discrepancy against
the analytic solution is purely spatial. This is what makes verification item 2 below a clean
measurement of the spatial order alone.

## Verification

1. **Against the analytic solution** at every output time: relative `L2` and max-norm over the
   plotted range. This is the headline number.
2. **Spatial order**, run with `start="warm"` and `expm_multiply` so neither the source
   smearing nor the time integration contaminates it: halve `h` at fixed `t0` and expect the
   error to fall by 4. Anything other than order 2 means a boundary row is wrong — that is what
   this test is really probing.
3. **Temporal order**, halving `dt` at fixed `h` with the Crank–Nicolson march: expect the error
   to fall by 4 until it hits the floor set by the spatial error. This is the study a
   hand-discretised scheme makes possible and a library integrator hides behind `rtol` — with
   `solve_ivp` the equivalent is tightening `rtol`, as `convergence_study` does for the shooting
   solver in `homework1/diffusion.py`, which measures tolerance-following rather than order.
4. **Particle balance**: `2 * trapezoid(phi, x) == exp(-(1-c)t)`, the factor 2 for the half
   domain.
5. **The claim in [06](06-delta-source-numerics.md)**: the difference between the `"pulse"` and
   `"warm"` runs on the same mesh isolates the source smearing. Check that `"pulse"` agrees
   with `phi_exact(x, t + dt_eff)`, `dt_eff = h^2/(24 D v)`, far better than with
   `phi_exact(x, t)`. That turns the claim into a measurement, and it is the one place the two
   `start` modes have to be run side by side.

## Pitfalls

- **The half-domain factor of 2** in every integral. Easy to lose, and it shows up as a clean
  factor-2 miss in the balance check rather than as anything subtle.
- **Trapezoid on point values is fine here**; the 1 % balance deficit recorded in
  `src/homework1/explanations/04-solver-history.md` came from applying it to *cell-averaged*
  data, which drops a half cell at each end. The node-centred grid avoids that entirely.
- **The grid depends on `D`, hence on both `c` and the approximation.** `L`, `h` and `t0` all
  scale with `D`, so build them per run rather than once. Sharing a grid between `c = 0.6`
  asymptotic (`D0 = 0.486`) and `c = 1.5` asymptotic (`D0 = 0.237`) silently changes the
  resolution by a factor of two in `h/sigma`.
- **Never integrate `phi` directly for `c > 1`.** At `c = 1.5, t = 15` it has grown by
  `e^{7.5} ~ 1808`; keep that in the analytic prefactor where it costs nothing.
- **`t0` is a method parameter, not a result.** Hold it fixed across a convergence study, and
  report it. It applies only to `start="warm"`.
- **The `u_0 = v/h` normalisation on a half domain** is where a factor of two hides. Verify it
  against the trapezoid rule before trusting anything downstream.
