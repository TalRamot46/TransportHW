# 04 — The Spherical Finite Volume

**Both boundary conditions are coefficients inside ordinary cell balances — no ghost cells, no
replaced rows, no branch at the origin.**

Report §"The Spherical Finite-Volume Scheme in Detail" derives the scheme and its banded
storage index by index. This file is what a reader of `spherical.py` needs on top of it.

## Where the boundary conditions actually live

`_banded_operator` assembles `N` rows, and every one of them is a plain cell balance. The two
conditions are already inside them:

| condition | how it is realised | line |
|---|---|---|
| `phi'(0) = 0` at the centre | `conductance[0] = D * areas[0] / h`, and `areas[0] = 0` | falls out of the geometry |
| leakage at the outer face | `conductance[-1] = D * areas[-1] / (0.5 * h)` | one overwritten entry |

The origin needs **no code path of its own**: face 0 has zero area, so row 0 has no inward
coupling and the symmetry condition is imposed exactly by `A_0 = 0`. This is the payoff of
writing the discretisation as a balance rather than differencing `(1/r^2) d/dr (r^2 d/dr)`
directly — the latter divides by `r^2` and needs an explicit limit at the origin.

The operator is assembled **once**, before the first sweep, so every iterate satisfies both
conditions exactly and no sweep can drift off them.

## The mesh runs to `R + z0`, and part of it is fictitious

`_mesh` is called with `R + medium.z0`, not `R`. The shell `R < r < R + z0` is a numerical
device: material properties are continued into it and the flux plotted there is extrapolation,
not physical flux. Worth remembering when reading `q4_flux_profiles.pdf`.

Cells are centred, so no unknown ever sits on a boundary, and `h` is uniform — no cell is
`h + z0` wide.

## Two structural properties worth checking assembly against

- **`M` is symmetric and positive definite.** `ab[0]` and `ab[2]` are slices of the same
  conductance array, and the diagonal dominates strictly once `sigma_a > 0`. So the LU inside
  `solve_banded` needs no pivoting and the solve is `O(N)`.
- **Row sums are the neutron balance.** Summing a row makes the conductances cancel in pairs,
  leaving `sigma_a V_j` for interior rows and `sigma_a V_{N-1} + C_N` for the last. That
  identity is what `neutron_balance` measures, and it is the only independent check on the
  boundary coefficient — see [08](08-verification.md).

## `dominance_ratio` computes, it does not measure

`dominance_ratio(medium, R)` evaluates the *modal* prediction `k_2/k_1` from two analytic
diffusion modes. It never runs the iteration. It is a predictor to compare the measured sweep
count against, not an instrument — see [05](05-k-iteration.md).
