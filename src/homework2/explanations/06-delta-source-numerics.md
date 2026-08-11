# 06 — Discretising a delta source

**What the first-cell delta actually costs, and why the time stepping is the real problem.**

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
| source | delta smeared into cell 0, `S_0 = 1/h` | **keep** — the error is `dt_eff = h^2/(24Dv)`, negligible on any sane mesh |
| time | explicit FTCS, `dt <= h^2/(2Dv)` | **Crank–Nicolson, implicit** — this is the real win; see [08](08-time-discretisation.md) |

So the first-cell source is fine and is what 3(c) should submit. What follows establishes *why*
it is fine rather than assuming it, and describes a warm-start variant that is more accurate
still but is better used as a diagnostic than as the answer.

## What the first-cell source costs

Replacing `delta(x)` by a top-hat of width `h` and height `1/h` is a change of the initial
condition, and its effect can be written down exactly. A top-hat of width `h` has variance
`h^2/12`; the diffusion Green's function at time `t` has variance `2Dvt`; and variances add
under convolution. So the smeared solution has second moment

```
2 D v t + h^2/12  =  2 D v ( t + h^2/(24 D v) )
```

i.e. **the smeared source reproduces the exact solution at a shifted time**,
`dt_eff = h^2/(24 D v)`. That is a much more useful statement than "the error is `O(h^2)`",
because it says *where* the error lives: entirely in the early-time behaviour, decaying in
relative terms as `dt_eff / t`.

At `D = 1/3`, `v = 1`:

| `h` | `dt_eff` | relative error at `t = 1` |
|---|---|---|
| 0.5 | 3.1e-2 | ~3 % |
| 0.1 | 1.3e-3 | ~0.1 % |
| 0.0125 | 2.0e-5 | ~2e-3 % |

So on any reasonable mesh the smeared delta is a small, quantified error, and it is smallest
exactly where the diffusion approximation is *most* accurate (large `t`). The honest conclusion
is that the first-cell source is not the weak point of the scheme. The stability limit is.

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

All three solutions are even in `x`, so solve on `[0, L]` with `du/dx = 0` at `x = 0` and put
half the source in the first cell. At the far end use a reflecting condition with `L` placed
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
mesh: no smearing error, no `dt_eff`, no first cell to normalise.

**Why this is a diagnostic and not the submission.** Handing the solver the analytic solution
at `t0` is a weaker answer to "write a numerical code for the diffusion equation" — the delta
was never actually solved for, and a reader could fairly say so. Its real value is that it
*separates the two error sources*: run it and the remaining discrepancy is purely spatial
discretisation, so the difference between the two runs measures the source smearing alone. That
turns `dt_eff = h^2/(24Dv)` from a claim into a measurement (verification item 5 in
[07](07-solver-implementation.md)).

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
| delta source | smeared over cell 0 | keep; cost is a time offset `dt_eff = h^2/(24Dv)` |
| time step | explicit, `dt <= h^2/(2Dv)`, ~6e4 steps | Crank–Nicolson + Rannacher startup, ~1.2e3 steps |
| far boundary | `phi(L) = 0`, reflects | reflecting, at an `L` the pulse has not reached |
| growth at `c > 1` | integrated numerically | factored out analytically |
| separating error sources | — | warm-start run, as a diagnostic |
