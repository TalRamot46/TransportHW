# Explanations — Assignment 2, Question 3

Background notes for the code in `src/homework2/`, whose report is §3 of
`docs/homework2/homework2.tex`. One idea per file; each file opens with its title and a
one-line statement of what it settles.

| # | Title | What it settles |
|---|---|---|
| [01](01-planar-flux-c1.md) | The planar flux at `c = 1` | Why the Q2 integral is closed-form, and why the code uses Dawson instead of `erfi` |
| [02](02-scaling-to-general-c.md) | Scaling to general `c` | How the Q1 identity carries the `c = 1` solution to every other `c`, including `c > 1` |
| [03](03-diffusion-coefficients.md) | The two diffusion coefficients | Where `D0(c) = (1-c) nu0^2` comes from, and why it stays positive above `c = 1` |
| [04](04-time-dependent-vs-steady.md) | Time-dependent vs. steady | Why Q3(a)(b) needs the Green's function, how the steady formula falls out of it, and where the `1/v` belongs |
| [05](05-verification.md) | Verification | The six checks in `main.py` and the numbers they produce |
| [06](06-delta-source-numerics.md) | Discretising a delta source | Why the delta should never reach the mesh, and what the first-cell approach actually costs |
| [07](07-solver-implementation.md) | Implementing the numerical solver | A concrete specification for the Question 3(c) code — not built yet |

Units throughout: `Sigma_t = v = 1`, as Assignment 2 specifies, so lengths are mean free
paths, times are mean free times, and the scalar flux equals the number density.
