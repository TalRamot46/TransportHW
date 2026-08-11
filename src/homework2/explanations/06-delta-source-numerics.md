# 06 — Discretising a delta source

**What seeding the delta on the mesh actually costs, and why the time stepping is the real problem.**

## Scope

Question 3(c) asks for a *time-dependent numerical code for the diffusion equation*. That is a
separate deliverable from 3(a) and 3(b): those are answered by the analytic Green's function
derived in the report (§3.2), and having a closed form does not discharge 3(c).

Nothing here proposes avoiding finite differences. The spatial discretisation stays the standard
centred second difference — see [07](07-solver-implementation.md). What is worth changing is
**how the delta enters** and **how time is stepped**:

| | usual choice | recommended |
|---|---|---|
| space | centred second difference | **unchanged** |
| source | delta seeded in node 0, `u_0 = 1/h` | **keep** — `O(h^2)`, measured `2.1e-5` on the production mesh |
| time | explicit FTCS, `dt <= h^2/(2Dv)` | **Crank–Nicolson, implicit** — this is the real win; see [08](08-time-discretisation.md) |

So seeding the delta directly is fine, and is what 3(c) should submit. What follows establishes *why*
it is fine rather than assuming it, and describes a warm-start variant that is more accurate
still but is better used as a diagnostic than as the answer.

## What seeding the delta costs

### The time-offset argument, and where it applies

If `delta(x)` is replaced by a top-hat of width `h` and height `1/h`, the effect can be written
down exactly. A top-hat has variance `h^2/12`; the Green's function at time `t` has variance
`2Dvt`; variances add under convolution. So the smeared solution has second moment

```
2 D v t + h^2/12  =  2 D v ( t + h^2/(24 D v) )
```

i.e. it reproduces the exact solution **at a shifted time**, `dt_eff = h^2/(24 D v)`.

That argument is sound, but it describes a **cell-centred** discretisation, where the delta
genuinely becomes a top-hat spanning a cell. The scheme specified in
[07](07-solver-implementation.md) is **node-centred**, with the spike sitting on the node at
`x = 0`, and there the argument does not apply — as the measurement below shows.

### Measured: the node-centred scheme has no time offset at all

Its discrete second moment is `sum h x_j^2 u_j = 0`, because the spike sits at `x_j = 0`. And
the centred second difference is *exact* on quadratics, so summation by parts gives

```
d/dt <x^2> = 2 D v     exactly, for any h
```

Measured at `t = 1`, `D = 1/3`, for `n = 200, 400, 800, 1600` nodes: `<x^2> = 0.666667` in every
case, against `2Dt = 0.666667`. **There is no `dt_eff` to measure.** The earlier claim that this
scheme carries a time offset was wrong.

### What it does cost

Seeding the delta instead costs a *shape* error near the origin, from high modes the mesh cannot
resolve — same order, larger constant. Against a warm start on the identical mesh (`c = 1`,
classical, `t = 1` and `4`):

| `n` | `h` | warm start | delta seeded | ratio |
|---|---|---|---|---|
| 200 | 0.0908 | 4.59e-4 | 1.82e-3 | 4.0 |
| 400 | 0.0453 | 1.17e-4 | 5.04e-4 | 4.3 |
| 800 | 0.0226 | 2.98e-5 | 1.06e-4 | 3.6 |
| 1600 | 0.0113 | 7.52e-6 | 2.22e-5 | 3.0 |

Both columns fall by about 4 per refinement, so both are `O(h^2)`; seeding the delta multiplies
the error constant by roughly 3.5 and does not degrade the order. On the production mesh the
absolute error is `2.1e-5`.

So the conclusion stands even though the reasoning first offered for it was wrong: seeding the
delta is not the weak point of the scheme. The stability limit is.

## The real cost: the explicit time step

An explicit (FTCS) march is stable only for `dt <= h^2 / (2 D v)`. At `h = 0.0125` and
`D = 1/3` that is `dt <= 2.3e-4`, so reaching `t = 15` takes about **64,000 steps** — and
halving `h` quadruples it. This, not the delta, is what makes the naive scheme painful.

## The improvements

Three for the submitted code, then one diagnostic. Only the first is essential.

### 1. An implicit march, not an explicit one

Replace FTCS with **Crank–Nicolson**, one tridiagonal solve per step. It is unconditionally
stable, so `dt` is set by accuracy rather than by `dt <= h^2/(2Dv)`, and it is second order in
time as well as space, so the balance is `dt ~ h` instead of `dt ~ h^2`: the 64,000 steps above
become about 1,200, at essentially the same cost per step.

One caveat that matters precisely because the source is a delta: CN rings on rough initial data,
and a delta is the roughest there is. Four backward-Euler quarter-steps at the start fix it for
free. Both the analysis and the recipe are in [08](08-time-discretisation.md), together with why
`solve_ivp` avoids the problem automatically and why it is still not the better answer here.

### 2. Half domain, with the symmetry condition at the origin

All three solutions are even in `x`, so solve on `[0, L]` with `du/dx = 0` at `x = 0` and
`u_0 = 1/h` at the origin — half the source, over the half cell the boundary node owns, so the
two factors of two cancel. At the far end use a reflecting condition with `L` placed
beyond anything the pulse has reached, rather than `phi(L) = 0` — the latter is what produced
the spurious 100 % boundary error in Assignment 1. Half the unknowns, and no artificial
reflection.

### 3. Factor the absorption out before integrating

Substituting `phi = e^{-(1-c) Sigma_t v t} u` (Step 2 of the derivation) leaves the integrator
solving the pure heat equation for `u`. This matters most for `c > 1`: at `c = 1.5, t = 15` the
flux has grown by a factor of `e^{7.5} ~ 1808`, and an integrator working on `phi` directly
must track that growth to the same relative tolerance as everything else. Working on `u`, the
growth is applied analytically at the end and costs nothing.

### 4. The warm-start variant — a diagnostic, not the answer

Integrating the equation across `t = 0` (Step 1 of the derivation) gives

```
phi(x, 0+) = v delta(x)
```

and for `t > 0` the problem is source-free. The delta is therefore not a source term at all —
it is *initial data*, and initial data can be supplied at any starting time. Begin the march at
a small `t0` from the analytic Gaussian, which is smooth and well resolved:

```
u(x, t0) = v / sqrt(4 pi D v t0) * exp(-x^2 / (4 D v t0))
```

with `t0` chosen so the width `sqrt(2 D v t0)` spans several cells. The delta never touches the
mesh, so the high modes it cannot resolve are never excited.

**Why this is a diagnostic and not the submission.** Handing the solver the analytic solution
at `t0` is a weaker answer to "write a numerical code for the diffusion equation" — the delta
was never actually solved for, and a reader could fairly say so. Its real value is that it
*separates the two error sources*: run it and the remaining discrepancy is purely spatial
discretisation, so the difference between the two runs measures the source smearing alone. That
is what produced the measured table above, and it is the only way to separate the two
(verification item 5 in [07](07-solver-implementation.md)).

The instinct behind it is still the right one, and it is the same as Assignment 1, Q2, where
integrating across the source turned the delta into a boundary current (`plan.md` §2.4: "the
delta is never discretised"). The difference is that there it *was* the answer, because the
steady problem has no other way in.

## What about a spectral method?

Because the medium is homogeneous and `D` constant, a Fourier method diagonalises the operator
exactly, and the delta transforms *exactly* — `delta_hat(k) = 1`, no discretisation error at
all. It would be spectrally accurate in space and exact in time.

It is also not really a numerical solution of the diffusion equation: on this problem it
reduces to evaluating the analytic answer with an FFT in the middle, and it stops working the
moment `D` or `c` varies with position, which is the case any solver would exist for. Worth
knowing, not worth submitting.

## Summary

| difficulty | naive scheme | recommended |
|---|---|---|
| spatial operator | centred second difference | unchanged |
| delta source | seeded in node 0 | keep; `O(h^2)`, ~3.5x the warm-start error constant |
| time step | explicit, `dt <= h^2/(2Dv)`, ~6e4 steps | Crank–Nicolson + Rannacher startup, ~1.2e3 steps |
| far boundary | `phi(L) = 0`, reflects | reflecting, at an `L` the pulse has not reached |
| growth at `c > 1` | integrated numerically | factored out analytically |
| separating error sources | — | warm-start run, as a diagnostic |
