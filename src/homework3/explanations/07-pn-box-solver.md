# 07 — The P_N Box Solver

**Method 1 of Question 2, in three files: `pn.py` builds the two matrices the report derives,
`pn_box.py` lays them out as one banded system and iterates `k` over it, `q2.py` tabulates.**

## The files, and the call path

| file | owns |
|---|---|
| `pn.py` | The algebra both methods share: `streaming_matrix`, `parity_blocks`, `marshak_matrix`. No solving. |
| `pn_box.py` | `BoxSystem` (assembly + factorisation), the power iteration, the size root search. |
| `pn_modal.py` | Method 2; see [08](08-modal-benchmark.md). |
| `q2.py` | The table and the figure, nothing else. |

    pn.streaming_matrix / pn.marshak_matrix
      <- BoxSystem._matrix          one sparse L, report eq. (20)
      <- BoxSystem.solve            one back-substitution against a midpoint source
      <- pn_box.pn_k_eigenvalue     the power iteration; returns sn.KResult
      <- pn_box.critical_half_thickness -> sn.critical_size   brentq on k(a/2) - 1
      <- q2.report

The last two lines are the reuse that matters: `sn.critical_size` and `sn.KResult` are taken
unchanged from the S_N side, so Questions 2 and 3 share one root finder, one bracket-widening
policy and one result type. Nothing in `sn.py` had to change to admit a method that has no
sweep and no ordinates.

## Why the matrix is laid out the way it is

Row order is **the (N+1)/2 symmetry rows, then one (N+1)-row block per cell, then the (N+1)/2
Marshak rows**. Any other order works algebraically; this one keeps the matrix banded, because
each cell block touches only nodes `j` and `j+1` and the two boundary groups sit at the two
ends. `splu` on the banded pattern is what makes the question cheap: the largest system here is
2010 by 2010 (`N = 9`, 200 cells), and `splu` factorises it in **3.0 ms** against **3.2 s** for
a dense `lu_factor` of the same matrix — a factor of a thousand, and the difference between one
(c, N) entry costing **0.22 s** and costing four minutes.

Two assembly details that look like tricks and are not:

- The cell block is written as one `hstack([collision - streaming, collision + streaming])` into
  a `width`-by-`2*width` slice. The left half multiplies `Phi_j`, the right half `Phi_{j+1}`;
  writing them together is what makes the loop body a single assignment.
- `rhs[n_conditions::width][:n_cells] = source` puts the fission source on the `n = 0` row of
  every cell block and nowhere else. The stride `width` steps one cell block; the offset
  `n_conditions` skips the symmetry rows; the `[:n_cells]` trim drops the Marshak rows, which
  the stride would otherwise reach.

`BoxSystem` factorises in its constructor, so one instance is one half-thickness. The outer
`brentq` therefore builds a new system per trial `a/2` — about eight per (c, N) — and each one
amortises its factorisation over the ten to twenty power iterations that follow.

## What was tried and rejected

**Absorbing `c/k` into the matrix.** The obvious simplification, and it destroys the method:
with the fission term on the left the system is homogeneous, and a homogeneous banded system
has only the zero solution except at the critical `k`, which is the unknown. The source has to
stay on the right, lagged one iteration. Report §2 states this; it is repeated here because the
code looks like it is missing a term until you know.

**A staggered mesh**, odd moments at cell centres and even moments at nodes, is the other
standard P_N discretisation and is also second order. It was not used because the Marshak rows
need every moment at the same point `x = a/2`, which a staggered mesh does not have and would
have to extrapolate — a first-order boundary closure bolted onto a second-order interior.

**A negative-flux fixup** has no counterpart here. `phi_n` for `n >= 1` is a current-like
quantity and is negative over half the slab by construction, so there is no sign to police;
[02](02-sn-solver.md)'s fixup is an S_N idea only.

## The one shared plot

`plots.py` grew two fields on `OrderScan` — `orders` and `family` — so Question 2 can reuse the
three-panel figure of Questions 3 and 4 with `N = 1, 3, 5, 9` and `P` labels. The defaults are
the old S_N values, and the Question 3 and 4 tables reproduce digit for digit after the change.
