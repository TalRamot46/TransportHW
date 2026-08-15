# 00 — Index

**These files explain the code in `src/homework2/`. The physics is in
`docs/homework2/homework2.tex` — read that first.**

Assignment 2 proves a scaling identity (Q1), derives the planar-source flux from Paasschens'
point-source solution (Q2), and compares it against classical and asymptotic diffusion (Q3).
The report derives all of it. What follows is the map from those formulas to the modules, the
numerical care the closed forms need, and what was verified.

Units throughout: `Sigma_t = v = 1`, so lengths are mean free paths, times are mean free times,
and the scalar flux equals the number density.

| # | Title | What it settles |
|---|---|---|
| [01](01-module-map.md) | The Module Map | Which file owns what, and why `main.py` is a test suite rather than a driver. |
| [02](02-evaluating-the-closed-form.md) | Evaluating the Closed Form | Why `exact.py` uses Dawson instead of `erfi`, and works in log space. |
| [03](03-diffusion-module.md) | The Diffusion Module | The one sign trick that makes `D0(c)` a single expression on both sides of `c = 1`. |
| [04](04-figures.md) | The Figures | The autoscaling heuristic, and why the exact curve is masked past the front. |
| [05](05-verification.md) | Verification | The six checks in `main.py` and the numbers they produce. |
| [06](06-q3c-solver-spec.md) | The Q3(c) Solver — Specification | The scheme chosen for code that is not written yet, and why. |

Nothing here re-derives the report. Where a physical result is needed it is cited by its
section or equation number in `homework2.tex`.
