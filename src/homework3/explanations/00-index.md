# 00 — Index

**These files explain the code in `src/homework3/`. The physics is in
`docs/homework3/homework3.tex` — read that first.**

Assignment 3's questions are 1 (reflected sphere, by diffusion), 2 (critical slab by P_N, twice
over) and 3–5 (critical slab, critical sphere and critical mass, all by S_N). The report derives
every formula and reports every number. What follows is the map from those formulas to the
modules, the internals that are not obvious from reading them, and what was done to trust the
results.

| # | Title | What it settles |
|---|---|---|
| [01](01-module-map.md) | The Module Map | Which file owns what, and the call path from `main.py` down to one table. |
| [02](02-sn-solver.md) | The S_N Solver | How one `cell_flux` serves both geometries, and the three traps in the sweep. |
| [03](03-k-iteration.md) | The Two-Level `k` Iteration | What each loop costs, and which tolerance is the cheap one. |
| [04](04-reflected-solver.md) | The Reflected-Sphere Solver | How three approximations collapse to one `Region` and one residual. |
| [05](05-assignment-1-reuse.md) | What Comes From Assignment 1 | Which imports carry physics and which are plumbing. |
| [06](06-verification.md) | Verification | Every check, in one place, with the number it produced. |
| [07](07-pn-box-solver.md) | The P_N Box Solver | How the banded system is laid out, and why the fission source cannot go in it. |
| [08](08-modal-benchmark.md) | The Modal Benchmark | The mesh-free second opinion, and the three ways it goes wrong if written naively. |

Nothing here re-derives the report. Where a physical result is needed it is cited by its
equation or table number in `homework3.tex`.
