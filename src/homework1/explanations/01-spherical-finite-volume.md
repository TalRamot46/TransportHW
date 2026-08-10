# 01 — Finite Volume in Spherical Coordinates

**Writing the discretisation as a neutron balance over spherical shells makes the
symmetry condition at the origin automatic: the face at `r = 0` has zero area.**

The equation to discretise is

    -D (1/r^2) d/dr (r^2 dphi/dr) + Sigma_a phi = (1/k) nu Sigma_f phi

Rather than discretising the derivatives directly, integrate over a spherical
shell cell `i` bounded by faces at `r_{i-1/2}` and `r_{i+1/2}`. The divergence
term becomes a difference of surface currents, and the balance reads

    A_{i+1/2} J_{i+1/2} - A_{i-1/2} J_{i-1/2} + Sigma_a V_i phi_i = S_i V_i

with `A = 4 pi r^2`, `V_i = (4pi/3)(r_{i+1/2}^3 - r_{i-1/2}^3)`, and Fick's law
`J = -D dphi/dr` approximated by a central difference between cell centres,
`J_{i+1/2} = -D (phi_{i+1} - phi_i) / h`.

## Why there is no special case at the origin

A direct discretisation of `(1/r^2) d/dr (r^2 dphi/dr)` divides by `r^2`, which
is singular at `r = 0` and forces an explicit L'Hôpital limit or a one-sided
symmetry condition. The balance form never divides by anything: the inner face
of the first cell is at `r = 0`, so `A_{1/2} = 0` and the term drops out on its
own. That *is* the condition `dphi/dr = 0` at the centre, imposed exactly rather
than approximated. In `_banded_operator` this is why the face conductance array
can be used without a branch for `i = 0`.

## Mesh layout

Cells are centred, not vertex-centred: `N` uniform cells on `[0, r_outer]` with
centres at `(i + 1/2) h`. The unknowns are cell averages, and no unknown sits on
a boundary, which keeps both the origin and the outer surface as face conditions
rather than as equations of their own.

The resulting matrix is tridiagonal and symmetric-definite (the fission source is
on the right-hand side, so the operator is pure removal plus leakage), and is
solved with `scipy.linalg.solve_banded` in `O(N)`.

## Measured order of accuracy

The scheme is second order. Refining the mesh at `c = 1.5` gives a relative error
of the critical radius falling from `5.28e-04` at `N = 25` to `5.16e-07` at
`N = 800`, a ratio of `4.00` per doubling at every step — see
`mesh_convergence` and figure `q4_mesh_convergence.pdf`. At the production
setting of `N = 400` the radius is within `2e-4 %` of analytic; see
[03](03-boundary-and-analytic.md).
