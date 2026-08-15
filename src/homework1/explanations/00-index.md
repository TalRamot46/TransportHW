# 00 — Index

**These files explain the code in `src/homework1/`. The physics is in
`docs/homework1/homework1.tex` — read that first.**

Assignment 1 runs from Case's exact plane-source flux (Q1), through a planar diffusion solver
(Q2) and the criticality relations (Q3), to a spherical finite-volume `k` solver (Q4) and the
bare critical masses (Q5). The report derives all of it. What follows is the map from those
formulas to the modules, the internals that are not obvious from reading them, the code that
was removed and why, and what was verified.

| # | Title | What it settles |
|---|---|---|
| [01](01-module-map.md) | The Module Map | Which file owns what, and the call path from `main.py` to one table. |
| [02](02-eigenvalue-routines.md) | The Eigenvalue Routines | Which of the six `nu0` entry points to call, and the default that will catch you out. |
| [03](03-shooting-solver.md) | The Shooting Solver | Why there is no mesh and no root-find, and why all six Q2 cases are one problem. |
| [04](04-spherical-finite-volume.md) | The Spherical Finite Volume | How both boundary conditions become coefficients rather than rows. |
| [05](05-k-iteration.md) | The `k` Iteration | Why the flux is renormalised every sweep, and why the sweep cap is 20000. |
| [06](06-materials-and-data.md) | The Benchmark Data | The two U-235 rows, and the self-check that guards the table. |
| [07](07-removed-code.md) | Removed Code and False Trails | Three solvers and one boundary option that no longer exist, and two dead ends not worth re-running. |
| [08](08-verification.md) | Verification | Every check, in one place, with the number it produced. |

Nothing here re-derives the report. Where a physical result is needed it is cited by its
section or equation number in `homework1.tex`.
