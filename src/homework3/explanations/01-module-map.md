# 01 — The Module Map

**Two unrelated solvers share one set of output helpers: `reflected.py` answers Question 1 in
closed form, and `sn.py` with `slab.py`/`sphere.py` answers Questions 3–5 by iteration.**

## The files

| file | owns |
|---|---|
| `main.py` | Entry point. Creates the figure directory, then calls `report(figs)` on each question in turn. |
| `q1.py`, `q3.py`, `q4.py`, `q5.py` | One question each: the tables it prints and the figures it writes, and nothing else. |
| `reflected.py` | Question 1: `Region` parameters, the criticality residual, the flux profile. |
| `sn.py` | Geometry-independent S_N: quadrature, the diamond cell solve with its fixup, the two-level `k` iteration, the size root search. |
| `slab.py` | The slab sweep and its mesh. |
| `sphere.py` | The spherical sweep, which adds the angular redistribution term. |
| `plots.py` | The three-panel order scan shared by Questions 3 and 4. |
| `figures.py` | Assignment 1's matplotlib helpers, re-pointed at `docs/homework3/figs/`. |

Question 1 has no `plots.py` counterpart — its two figures are unlike the order scans, so they
live in `q1.py`.

## The call path, Questions 3–5

From the bottom up:

    sn.cell_flux             one cell, one link per outgoing face
      <- Solver._direction   one ordinate, swept across the mesh
      <- Solver.sweep        all ordinates -> scalar flux
      <- sn.inner_iteration  sweeps until the scattering source settles
      <- sn.k_eigenvalue     outer loop; returns KResult(k, x, phi, outers)
      <- sn.critical_size    brentq on k(size) - 1
      <- qN.report           table + figure

**`Solver` is a duck type, not a base class.** `sn.k_eigenvalue` touches only `.medium`,
`.n_cells`, `.volumes`, `.centres` and `.sweep(source)`; `SlabSolver` and `SphereSolver` each
supply those independently, with no shared parent and no registration. That is why adding the
sphere required no change at all to `sn.py`, and it is the seam to use for any further geometry.

## The call path, Question 1

Shallower, because there is nothing to iterate:

    reflected.region              (material, theory) -> Region
      <- reflected._residual      the interface balance, report eq. (11)
      <- reflected.critical_radius  brentq on (0, pi/B)
      <- q1.report                tables + figures

`q1.py` holds no physics beyond the choice of what to tabulate: the material and thickness
loops, the `BARE` translation table onto Assignment 1's naming, and the two plot builders.
