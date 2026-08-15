# 00 — Index

Explanations for the Assignment 3 code in `src/homework3/`, which covers Questions 3, 4 and 5:
the critical slab half-thickness, the critical sphere radius, and the critical mass of U-235 and
Pu-239, all by S_N. Results and figures are in `docs/homework3/homework3.tex`.

| # | Title | What it settles |
|---|---|---|
| [01](01-sn-discretisation.md) | The S_N Sweep in a Slab and in a Sphere | One cell balance for both geometries, and where the `alpha` recursion comes from. |
| [02](02-splitting-c.md) | Splitting `c` into Scattering and Fission | Why the critical size does not depend on the split, so `Sigma_s = 0` is free. |
| [03](03-two-level-iteration.md) | The Two-Level Bell & Glasstone Iteration | What the inner and outer loops each cost, and the two dominance ratios that set it. |
| [04](04-convergence-and-checks.md) | What Was Checked, and What the Mesh Costs | That the mesh is free, so every tabulated departure is angular truncation. |
| [05](05-slab-sphere-and-pn.md) | The Slab/Sphere Relation, and `P_N` vs `S_{N+1}` | Why the two geometries converge from opposite sides, and where `P_N` and `S_{N+1}` part. |
| [06](06-negative-flux-fixup.md) | The Negative-Flux Fixup, and Why It Never Fires Here | That it is correct, and that these problems never reach it. |

Question 3 is covered by 01–06, Question 4 by 01 and 05, and Question 5 by 02 and 03.
Assignment 1 is reused throughout: `homework1.criticality` supplies the exact-transport
reference sizes, `homework1.materials` the benchmark cross sections and the mass formula, and
`homework1.figures` and `homework1.tables` the output helpers.
