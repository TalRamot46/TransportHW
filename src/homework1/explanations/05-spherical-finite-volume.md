# 05 — Finite Volume in the Sphere

**Writing the discretisation as a neutron balance makes the symmetry condition at `r = 0`
automatic — the innermost face has zero area — and one leakage coefficient covers both
outer boundary treatments.**

Integrating `-D grad^2 phi + Sigma_a phi = S` over shell `i` gives

    A_{i+1/2} J_{i+1/2} - A_{i-1/2} J_{i-1/2} + Sigma_a V_i phi_i = S_i V_i

with `A = 4 pi r^2`, `V_i = (4pi/3)(r_{i+1/2}^3 - r_{i-1/2}^3)` and
`J_{i+1/2} = -D (phi_{i+1} - phi_i)/h`. Discretising the derivatives directly would divide
by `r^2` and need an explicit limit at the origin; the balance form never divides, and
`A_{1/2} = 0` imposes `phi'(0) = 0` exactly with no special case. Cells are centred, so no
unknown sits on a boundary. The matrix is tridiagonal, solved in `O(N)` by `solve_banded`.

## The outer face

Both standard treatments are the single condition `phi + l0 phi' = 0`. Eliminating the
surface flux gives one leakage coefficient with no `l0` in the numerator,

    J_s = D phi_{N-1} / (l0 + h/2)

so `l0 = 0` is not a special case — it degenerates to the zero-flux form.

| `boundary` | mesh runs to | `l0` | condition |
|---|---|---|---|
| `'extrapolated'` | `R + z0` | `0` | `phi(R + z0) = 0` |
| `'robin'` | `R` | `z0` | `phi(R) + z0 phi'(R) = 0` |

`'extrapolated'` is the default, because it is what the analytic relation assumes, making
the comparison a test of the numerics rather than of the boundary model.

## What it is compared against

At `k = 1` the equation is Helmholtz with `B^2 = (nu Sigma_f - Sigma_a)/D`, whose regular
solution `sin(Br)/r` vanishes at the extrapolated radius when `R_c = pi/B - z0`. The radius
depends on the cross sections only through `c`, since `Sigma_a` and `nu Sigma_f` enter only
as their difference. Substituting the two parameter sets recovers exactly the Question 3
relations ([04](04-criticality-relations.md)).

## Measured

At `N = 400` over `c = 1.02 … 2.0` the numerical radius is within `1.6e-4 %` to `2.4e-4 %`
of analytic, for both approximations. Refining at `c = 1.5` takes the relative error from
`5.28e-04` at `N = 25` to `5.16e-07` at `N = 800` — a factor `4.00` per doubling, i.e.
second order. The sign is consistently negative: the discrete operator is slightly too leaky
on a coarse mesh, so it reaches `k = 1` at a marginally smaller radius. `k` at the returned
radius is 1 to ten digits and the flux matches `sin(Br)/Br` to `3.2e-06`.
