# 01 — The Module Map

**Three unrelated solvers share one set of output helpers: `reflected.py` answers Question 1 in
closed form, `pn/box.py`/`pn/modal.py` answer Question 2 twice over, and `sn/core.py` with
`sn/slab.py`/`sn/sphere.py` answers Questions 3–5 by iteration.**

## The files

The two iterative methods are subpackages; everything at the top level is either a question
driver or shared output plumbing.

| file | owns |
|---|---|
| `main.py` | Entry point. Creates the figure directory, then calls `report(figs)` on each question in turn. |
| `q1.py` … `q5.py` | One question each: the tables it prints and the figures it writes, and nothing else. |
| `reflected.py` | Question 1: `Region` parameters, the criticality residual, the flux profile. |
| `pn/algebra.py` | The P_N algebra shared by Question 2's two solutions: the streaming matrix, its parity blocks, the Marshak rows. |
| `pn/box.py` | Question 2's default solution: the banded box system and the `k` power iteration over it. |
| `pn/modal.py` | Question 2's analytic check: the modal elimination and `det H(a/2) = 0`. |
| `sn/core.py` | Geometry-independent S_N: quadrature, the two-level `k` iteration, the size root search. The cell solve and the sweeps live with their geometry. |
| `sn/slab.py` | The slab cell solve, its two sweeps, and its mesh. |
| `sn/sphere.py` | The spherical cell solve with its radial and angular faces, and its two sweeps. |
| `plots.py` | The three-panel order scan shared by Questions 2, 3 and 4. |
| `figures.py` | Assignment 1's matplotlib helpers, re-pointed at `docs/homework3/figs/`. |

**The split is by method, not by question.** `sn/` serves Questions 3, 4 and 5 and `pn/` serves
only Question 2, so the `qN.py` drivers stay outside both: `q5.py` reaches into `sn/` for its
solver and into Assignment 1 for its materials, and belongs to neither. The one cross-package
import is `pn/box.py` taking `critical_size` and `KResult` from `sn/core.py` — the reuse noted
in [07](07-pn-box-solver.md), and the reason `sn/` does not import from `pn/` in either
direction.

Question 2's own call path and layout are in [07](07-pn-box-solver.md) and
[08](08-modal-benchmark.md); the two entries it borrows from `sn/core.py` — `critical_size` and
`KResult` — are noted there.

Question 1 has no `plots.py` counterpart — its two figures are unlike the order scans, so they
live in `q1.py`.

## The call path, Questions 3–5

From the bottom up:

    slab.cell_flux /              one cell of the march, per geometry
    sphere.cell_flux
      <- Solver.sweep_backward /  one ordinate, marched across the mesh
         Solver.sweep_forward
      <- Solver.sweep_all_angles  every backward ordinate, then every forward one
      <- sn.core.run_sn_for_source  repeats that until the scattering source settles
      <- sn.core.k_eigenvalue       outer loop; returns KResult(k, x, phi, outers)
      <- sn.core.critical_size      brentq on k(size) - 1
      <- qN.report                table + figure

The middle names are deliberately graded, because the distinction is easy to lose:
**one `sweep_backward` or `sweep_forward` is one ordinate**, **one `sweep_all_angles` is every
ordinate once**, and
**`run_sn_for_source` is the S_N method itself** — it repeats `sweep_all_angles` to convergence.
It is named for what it is rather than for where it sits: it is a complete S_N solve at a fixed
fission source and would still be one if `k_eigenvalue` did not exist.

**`Solver` is a duck type, not a base class.** `sn.core.k_eigenvalue` and `sn.core.run_sn_for_source`
between them touch only `.medium`, `.n_cells`, `.volumes`, `.centres` and
`.sweep_all_angles(source)`; `SlabSolver` and `SphereSolver` each supply those independently,
with no shared parent and no registration. That is why adding the sphere required no change at
all to `sn/core.py`, and it is the seam to use for any further geometry.

## The call path, Question 1

Shallower, because there is nothing to iterate:

    reflected.region              (material, theory) -> Region
      <- reflected._residual      the interface balance, report eq. (4)
      <- reflected.critical_radius  brentq on (0, pi/k0)
      <- q1.report                tables + figures

`q1.py` holds no physics beyond the choice of what to tabulate: the material and thickness
loops, the `BARE` translation table onto Assignment 1's naming, and the two plot builders.
