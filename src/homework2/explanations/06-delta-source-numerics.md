# 06 — Discretising a delta source

**Why the delta should never reach the mesh, and what the first-cell approach actually costs.**

For Question 3(c) — a time-dependent numerical diffusion code. The obvious scheme is finite
differences with the source placed in the first cell at the correct normalisation,
`S_0 = 1/dx`. It works, and its error is quantifiable. But there is a strictly better option
that costs nothing, and the delta is not where the real difficulty lies anyway.

## What the first-cell source costs

Replacing `delta(x)` by a top-hat of width `h = dx` and height `1/h` is a change of the
initial condition, and its effect can be written down exactly. A top-hat of width `h` has
variance `h^2/12`; the diffusion Green's function at time `t` has variance `2Dvt`; and
variances add under convolution. So the smeared solution has second moment

```
2 D v t + h^2/12  =  2 D v ( t + h^2/(24 D v) )
```

i.e. **the smeared source reproduces the exact solution at a shifted time**,
`dt_eff = h^2/(24 D v)`. That is a much more useful statement than "the error is `O(h^2)`",
because it says *where* the error lives: entirely in the early-time behaviour, decaying in
relative terms as `dt_eff / t`.

At `D = 1/3`, `v = 1`:

| `dx` | `dt_eff` | relative error at `t = 1` |
|---|---|---|
| 0.5 | 3.1e-2 | ~3 % |
| 0.1 | 1.3e-3 | ~0.1 % |
| 0.0125 | 2.0e-5 | ~2e-3 % |

So on any reasonable mesh the smeared delta is a small error, and it is smallest exactly where
the diffusion approximation is *most* accurate (large `t`). The honest conclusion is that the
first-cell source is not the weak point of that scheme. The stability limit is.

## The real cost: the explicit time step

An explicit (FTCS) march is stable only for `dt <= dx^2 / (2 D v)`. At `dx = 0.0125` and
`D = 1/3` that is `dt <= 2.3e-4`, so reaching `t = 15` takes about **64,000 steps** — and
halving `dx` quadruples it. This, not the delta, is what makes the naive scheme painful.

## The better scheme

Four changes, in decreasing order of importance. Each one removes a difficulty rather than
resolving it more finely.

### 1. Move the delta into the initial condition, then warm-start

Integrating the equation across `t = 0` (Step 1 of the derivation in the report) gives

```
phi(x, 0+) = v delta(x)
```

and for `t > 0` the problem is source-free. The delta is therefore not a source term at all —
it is initial data, and initial data can be supplied at *any* starting time. Begin the march
at a small `t0` from the analytic Gaussian, which is smooth and perfectly resolved:

```
phi(x, t0) = v e^{-(1-c) Sigma_t v t0} / sqrt(4 pi D v t0) * exp(-x^2 / (4 D v t0))
```

pick `t0` so the width `sqrt(2 D v t0)` spans several cells (`t0 ~ 25 dx^2 / (2Dv)` gives five),
and march to `t = 1 ... 15`. **The delta never touches the mesh.** No smearing error, no
`dt_eff`, no first cell to normalise.

This is not circular. `t0` is an initial condition, not the answer; the scheme is still being
asked to propagate it correctly over the remaining `15 - t0`, and it is still checked against
the analytic solution at every output time. It is the same instinct as Assignment 1, Q2, where
integrating across the source turned the delta into a boundary current
(`plan.md` §2.4: "the delta is never discretised").

### 2. Method of lines, not a hand-rolled stepper

Discretise `x` only, leaving a linear ODE system `dphi/dt = A phi`, and hand it to
`scipy.integrate.solve_ivp` with `BDF` or `LSODA`. `A` is tridiagonal, so pass `jac_sparsity`
and the implicit solve stays `O(N)` per step. This removes the `dt <= dx^2/(2Dv)` restriction
entirely — an implicit method takes time steps set by *accuracy*, not stability, so the 64,000
steps above become a few hundred. Accuracy is then controlled by one knob, `rtol`, and
convergence is demonstrated by tightening it, exactly as `convergence_study` does for the
shooting solver in `homework1/diffusion.py`.

### 3. Half domain, with the symmetry condition at the origin

All three solutions are even in `x`. Solve on `[0, L]` with `dphi/dx = 0` at `x = 0` — for a
*pulse* source this is pure symmetry, with no source term in it, since the source has already
been absorbed into the initial condition by change 1. At the far end use the radiation
condition `dphi/dx = -kappa phi` with `kappa = sqrt(Sigma_a/D)` rather than `phi(L) = 0`; the
latter is what produced the spurious 100 % boundary error in Assignment 1 before it was
replaced. Half the unknowns, and no artificial reflection.

### 4. Factor the absorption out before integrating

Substituting `phi = e^{-(1-c) Sigma_t v t} u` (Step 2 of the derivation) leaves the integrator
solving the pure heat equation for `u`. This matters most for `c > 1`: at `c = 1.5, t = 15` the
flux has grown by a factor of `e^{7.5} ~ 1808`, and an integrator working on `phi` directly
must track that growth to the same relative tolerance as everything else. Working on `u`, the
growth is applied analytically at the end and costs nothing.

## What about a spectral method?

Because the medium is homogeneous and `D` constant, a Fourier method diagonalises the operator
exactly, and the delta transforms *exactly* — `delta_hat(k) = 1`, no discretisation error at
all. It would be spectrally accurate in space and exact in time.

It is also not really a numerical solution of the diffusion equation: on this problem it
reduces to evaluating the analytic answer with an FFT in the middle, and it stops working the
moment `D` or `c` varies with position, which is the case any solver would exist for. Worth
knowing, not worth submitting.

## Summary

| difficulty | naive scheme | better |
|---|---|---|
| delta source | smeared over cell 0, cost `dt_eff = dx^2/(24Dv)` | absorbed into the initial condition; warm-start at `t0` |
| time step | explicit, `dt <= dx^2/(2Dv)`, ~6e4 steps | method of lines + BDF, accuracy-limited |
| far boundary | `phi(L) = 0`, reflects | radiation condition, `dphi/dx = -kappa phi` |
| growth at `c > 1` | integrated numerically | factored out analytically |
