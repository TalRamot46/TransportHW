# 07 — Implementing the numerical solver

**A concrete specification for the Question 3(c) time-dependent diffusion code.**

The scheme and its justification are in [06](06-delta-source-numerics.md); this file is the
implementation detail — enough that writing `solver.py` is mechanical. Nothing here is built
yet.

## What is being solved

After the two exact reductions of the derivation (report §3.2.2), the equation handed to the
solver is the **pure heat equation**:

```
du/dt = D v d2u/dx2,        phi(x,t) = e^{-(1-c) Sigma_t v t} u(x,t)
```

with the delta gone from the source and the absorption gone from the operator. Both the
`c`-dependence and the `delta` have been removed analytically; what remains for the numerics is
the one part that genuinely needs discretising.

## Module layout

`src/homework2/solver.py`, four short functions:

```python
def _grid(D, t_max, x_max, n_nodes)       -> x, h        # node-centred half domain [0, L]
def _laplacian(n_nodes, h, D)             -> A           # sparse tridiagonal, BCs baked in
def _warm_start(x, t0, D)                 -> u0          # analytic Gaussian at t0
def solve_diffusion(c, approximation, times, ...) -> (x, {t: phi})
```

`solve_diffusion` is the only public one. It picks `D = diffusion_coefficient(c, approximation)`
from `diffusion.py`, builds the grid and operator, warm-starts, integrates, and reapplies the
absorption factor on the way out.

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

**At `j = 0` — symmetry.** The source has already been absorbed into the initial condition, so
for `t > 0` there is no source at the origin and the condition is pure evenness,
`du/dx(0,t) = 0`. The ghost node is `u_{-1} = u_1`, and the row becomes

```
(2 u_1 - 2 u_0) / h^2
```

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

## Warm start

Choose `t0` so the initial Gaussian is resolved by the mesh: its width `sqrt(2 D v t0)` should
span several cells. Five gives

```
t0 = (5h)^2 / (2 D v)
```

which at `h = 0.0125`, `D = 1/3` is `t0 ~ 5.9e-3` — three orders of magnitude before the first
output time. The initial data is the analytic solution **for `u`**, i.e. without the absorption
factor:

```
u(x, t0) = v / sqrt(4 pi D v t0) * exp(-x^2 / (4 D v t0))
```

Record `t0` in the returned metadata: it is a parameter of the method, and the convergence
study has to hold it fixed while `h` varies, or the two effects mix.

## Time integration

Method of lines: `du/dt = A u` with `A` constant, linear and sparse. Two options, and it is
worth implementing both because they answer different questions.

**`solve_ivp(..., method="BDF", jac=lambda t, u: A, t_eval=times, rtol=1e-10)`.** Implicit, so
the step size is set by accuracy rather than by the explicit stability limit
`dt <= h^2/(2Dv) = 2.3e-4` — which would need about **64,000 steps** to reach `t = 15`. Passing
`jac` explicitly matters: without it BDF finite-differences a 4000x4000 Jacobian. `atol` should
be scaled to the initial peak (`u0.max() ~ 10` here), not left at its `1e-6` default.

**`scipy.sparse.linalg.expm_multiply(A * (t - t0), u0)`.** The *exact* solution of the
semi-discrete system — zero time-integration error by construction. This is the more useful of
the two as a diagnostic: any discrepancy against the analytic solution is then purely the
spatial discretisation, so the two error sources can be separated instead of being reported as
one number.

## Verification

1. **Against the analytic solution** at every output time: relative `L2` and max-norm over the
   plotted range. This is the headline number.
2. **Spatial order**, using `expm_multiply` so no time error contaminates it: halve `h` at
   fixed `t0` and expect the error to fall by 4. Anything other than order 2 means a boundary
   row is wrong — that is what this test is really probing.
3. **Temporal convergence**, using BDF at fixed `h`: tighten `rtol` and watch the error fall to
   the floor set by the spatial error. Same shape of study as `convergence_study` in
   `homework1/diffusion.py`.
4. **Particle balance**: `2 * trapezoid(phi, x) == exp(-(1-c)t)`, the factor 2 for the half
   domain.
5. **The claim in [06](06-delta-source-numerics.md)**: run the alternative — full domain, delta
   smeared into cell 0 at `S_0 = 1/h`, started at `t = 0` — and check its error against the
   analytic solution matches a time offset `dt_eff = h^2/(24 D v)`, i.e. that
   `phi_smeared(x,t)` agrees with `phi_exact(x, t + dt_eff)` far better than with
   `phi_exact(x,t)`. That turns the claim into a measurement.

## Pitfalls

- **The half-domain factor of 2** in every integral. Easy to lose, and it shows up as a clean
  factor-2 miss in the balance check rather than as anything subtle.
- **Trapezoid on point values is fine here**; the 1 % balance deficit recorded in `plan.md`
  §5.2 came from applying it to *cell-averaged* data, which drops a half cell at each end. The
  node-centred grid avoids that entirely.
- **The grid depends on `D`, hence on both `c` and the approximation.** `L`, `h` and `t0` all
  scale with `D`, so build them per run rather than once. Sharing a grid between `c = 0.6`
  asymptotic (`D0 = 0.486`) and `c = 1.5` asymptotic (`D0 = 0.237`) silently changes the
  resolution by a factor of two in `h/sigma`.
- **Never integrate `phi` directly for `c > 1`.** At `c = 1.5, t = 15` it has grown by
  `e^{7.5} ~ 1808`; keep that in the analytic prefactor where it costs nothing.
- **`t0` is a method parameter, not a result.** Hold it fixed across a convergence study, and
  report it.
