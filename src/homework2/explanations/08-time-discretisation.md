# 08 — Discretising time

**Why `solve_ivp` is a time discretisation too, and which fully-discrete scheme to write instead.**

`solve_ivp` is not an alternative to discretising time — BDF *is* a time discretisation, a
variable-order backward differentiation formula. The real distinction is **who chooses the step
size and the order**: a hand-written scheme uses the `dt` and order you fix in advance, a
library integrator adapts both to meet a tolerance. So the question "faster convergence or more
accuracy?" has a third answer: neither, it is about **stability and step-size selection**.

This file settles which scheme to write. Short version: **Crank–Nicolson with a
backward-Euler startup**, hand-written, each step a tridiagonal solve.

## The four candidates

Discretise space first (centred second difference, [07](07-solver-implementation.md)) to get
`du/dt = A u`, then discretise time. With `r = D v dt / h^2` the schemes are

```
FTCS (explicit):    u^{n+1} = u^n + dt A u^n
BE   (implicit):    (I - dt A) u^{n+1} = u^n
CN   (trapezoid):   (I - (dt/2) A) u^{n+1} = (I + (dt/2) A) u^n
```

All three are tridiagonal; the implicit ones cost one banded solve per step
(`scipy.linalg.solve_banded`, or the Thomas algorithm), which is `O(N)` — the same order as the
explicit matrix–vector product, and only a small constant more expensive.

## Von Neumann analysis

Substituting a Fourier mode `u_j^n = g^n e^{i j theta}` gives the amplification factor per step.
The eigenvalue of `A` for that mode is `lambda = -(4Dv/h^2) sin^2(theta/2)`, so

| scheme | `g(theta)` | stable when | `g` at `theta = pi`, large `r` |
|---|---|---|---|
| FTCS | `1 - 4r sin^2(theta/2)` | `r <= 1/2` | diverges |
| BE | `1 / (1 + 4r sin^2(theta/2))` | all `r` | `-> 0` |
| CN | `(1 - 2r sin^2(theta/2)) / (1 + 2r sin^2(theta/2))` | all `r` | `-> -1` |

Three consequences, and each one decides something:

**FTCS's limit is stability, not accuracy.** Violating `r <= 1/2` does not make the answer
worse, it makes it blow up. And the constraint is `dt <= h^2/(2Dv)`, so halving `h` quarters
`dt` — the cost grows as `h^-3`.

**Unconditional stability alone buys nothing.** BE is stable at any `dt`, but it is only
first-order in time, so matching a spatial error of `O(h^2)` still needs `dt ~ h^2` — the same
step count as FTCS, at a higher cost per step. BE alone is strictly worse than FTCS here.

**CN is the one that pays.** Second order in both, unconditionally stable, so `dt` is set by
balancing `O(dt^2)` against `O(h^2)`, i.e. `dt ~ h` rather than `dt ~ h^2`.

## What that is worth

At `h = 0.0125`, `D = 1/3`, `v = 1`, integrating to `t = 15`:

| scheme | order | `dt` set by | steps |
|---|---|---|---|
| FTCS | `O(dt) + O(h^2)` | stability, `dt <= 2.3e-4` | ~64,000 |
| BE | `O(dt) + O(h^2)` | accuracy, `dt ~ h^2` | ~64,000 |
| CN | `O(dt^2) + O(h^2)` | accuracy, `dt ~ h` | ~1,200 |
| BDF, adaptive | variable order | tolerance | a few hundred |

So the gain from CN is **cost at fixed accuracy — a factor of ~50 in step count** — not accuracy
itself. FTCS run at its stability limit is perfectly accurate; it is just forced to take steps
far smaller than accuracy requires. That is the honest answer to "faster convergence or more
accuracy": it is neither, it is that stability stops dictating the step size.

## The trap: Crank–Nicolson rings on a delta

CN is unconditionally *stable* but not unconditionally *smooth*. Look again at its `g`: the
numerator changes sign when `2r sin^2(theta/2) > 1`, so the highest modes come back with `g`
close to `-1` — they alternate in sign every step and barely decay. The sign flip happens at
`r > 1/2`, i.e. **exactly beyond the FTCS stability limit**, which is the whole range CN exists
to be used in.

Normally this is harmless because smooth initial data has little high-wavenumber content. Here
it is the worst case: the initial condition is a delta, which **excites every mode at equal
amplitude**. The result is a persistent sawtooth near the origin that decays only as `|g|^n`
with `|g| ~ 1`.

This is a known failure mode, not a bug to hunt — see Rannacher, *Finite element solution of
diffusion problems with irregular data*, Numer. Math. 43 (1984) 309.

### Measured

Ten CN steps from a delta at `h = 0.0125`, `dt = h`, `D = 1/3`, so `r = 26.7` and
`g(pi) = -0.963`:

| run | `min(u)` near the origin | sign flips | max relative error |
|---|---|---|---|
| CN alone | **`-13.4`** | 2 | **`2.2e+1`** |
| CN + 4 backward-Euler quarter-steps | `+6.8e-2` | 0 | `1.3e-3` |

The scheme does not merely lose accuracy — it returns a **negative flux**, and is wrong by a
factor of twenty. Four cheap startup steps take it to three good digits. Note also that `g(pi)`
decays as `0.963^n`, so the oscillation is still at 69 % of its initial amplitude after ten
steps and does not quietly go away with more marching.

### The fix: Rannacher startup

Take the first few steps with **backward Euler at a reduced step**, then switch to CN. BE has
`g -> 0` for the high modes (it is L-stable), so a handful of BE steps annihilates the rough
content that CN would otherwise ring on; the rest of the march is CN and keeps second-order
accuracy. The standard recipe is four BE steps of size `dt/4`, or two of `dt/2`:

```
steps 1..4 :  (I - (dt/4) A) u^{n+1} = u^n        # backward Euler, quarter step
steps 5..  :  (I - (dt/2) A) u^{n+1} = (I + (dt/2) A) u^n   # Crank-Nicolson
```

The cost is negligible and the global order stays 2. **Verify it rather than assume it**: plot
`u` near the origin after ten steps with and without the startup — the sawtooth is unmistakable.

This also explains why the warm-start mode of [07](07-solver-implementation.md) is useful for
more than diagnostics: started from a resolved Gaussian instead of a delta, CN needs no
Rannacher phase at all, because the rough modes were never excited.

### Why BDF sidesteps it for free

BDF1 and BDF2 are both L-stable, and `solve_ivp` starts at order 1 and raises the order as the
solution smooths. That is Rannacher startup, arrived at automatically by the order-selection
logic rather than by hand. It also shortens `dt` near `t = 0` where the solution is sharp and
lengthens it later where the Gaussian is nearly static — over `t = 0` to `15` the natural time
scale changes by orders of magnitude, and a fixed-`dt` scheme has to use the smallest step
everywhere.

That convenience is the honest case for `solve_ivp`, and it is not the same as being a better
method.

## Recommendation

**Write Crank–Nicolson with the Rannacher startup.** It is a genuine fully-discrete scheme —
which is what Question 3(c) asks for — it is about fifteen lines given a tridiagonal solve, it
is second order in both variables, and it makes the two discretisation errors independently
visible: refine `h` at fixed `dt` for the spatial order, refine `dt` at fixed `h` for the
temporal order. A library integrator hides the second of those behind a tolerance.

Keep `solve_ivp(method="BDF")` as a cross-check: agreement between an adaptive library
integrator and a hand-written fixed-step scheme is strong evidence that the *spatial* operator
and its boundary rows are right, since that is the only thing the two share.
