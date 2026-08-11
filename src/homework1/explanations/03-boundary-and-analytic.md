# 03 — Boundary Conditions and the Analytic Radius

**One leakage coefficient, `D A / (l0 + h/2)`, covers both the extrapolated-zero
and the Robin boundary treatments; with `l0 = 0` on a mesh extended to `R + z0`
it reproduces the analytic `R_c = pi/B - z0` to `2e-4 %`.**

## The outer face

Both standard treatments of a bare surface are the single condition

    phi + l0 phi' = 0

applied at the outer mesh boundary, with `l0` the linear extrapolation length.
Eliminating the (unknown) surface flux between that condition and the
finite-difference current from the last cell centre gives the leakage through
the outer face:

    phi_s = phi_{N-1} * 2 l0 / (h + 2 l0)
    J_s   = D phi_s / l0 = D phi_{N-1} / (l0 + h/2)

The last form has no `l0` in the numerator, so `l0 = 0` is not a special case —
it degenerates cleanly to `D phi_{N-1} / (h/2)`, the zero-flux condition. That
is why `_banded_operator` takes a single `l0` and writes one expression.

The two treatments are then selected in `_outer_boundary`:

| `boundary` | mesh runs to | `l0` | condition realised |
|---|---|---|---|
| `'extrapolated'` | `R + z0` | `0` | `phi(R + z0) = 0` |
| `'robin'` | `R` | `z0` | `phi(R) + z0 phi'(R) = 0` |

`'extrapolated'` is the default, because it is what the analytic relations of
Question 3 assume and therefore what makes the comparison a like-for-like test
of the *numerics* rather than of the boundary model.

## The Robin condition as an independent check

`'robin'` is the analogue of `critical_dimensions_applied_bc` in
`criticality.py`, which solves `u cot u = 1 - u/(B l0)` on the flux shape rather
than placing an extrapolated zero. That gives a second, completely independent
analytic route to compare the solver against — one that exercises the boundary
coefficient itself rather than only the interior operator. Measured for
classical diffusion:

| `c` | extrapolated `R_c` | Robin `R_c` | Robin analytic | Robin vs. extrapolated |
|---|---|---|---|---|
| 1.02 | 12.158813 | 12.126930 | 12.126949 | `-0.26 %` |
| 1.10 | 5.069062 | 5.007865 | 5.007874 | `-1.21 %` |
| 1.50 | 1.898429 | 1.821650 | 1.821653 | `-4.04 %` |
| 2.00 | 1.147130 | 1.095638 | 1.095640 | `-4.49 %` |

The numerical Robin radius reproduces its own analytic counterpart to about
`1.5e-6` relative — the same discretisation error as the extrapolated case, so
the boundary coefficient `D A / (l0 + h/2)` is correct in both of its limits.

The two treatments agree to first order in `B l0` and separate as the system
shrinks, from `0.26 %` apart at `c = 1.02` to `4.5 %` at `c = 2`. That is not
solver error: it is the genuine ambiguity in what "the surface" means once the
system is only a few mean free paths across. It propagates to the critical mass
as an 8–12 % spread — see [05](05-critical-mass-results.md).

## What the numerical radius is compared against

Setting `k = 1` leaves `-D grad^2 phi + (Sigma_a - nu Sigma_f) phi = 0`, a
Helmholtz equation with buckling

    B^2 = (nu Sigma_f - Sigma_a) / D = Sigma_t (c - 1) / D

whose regular spherical solution is `sin(Br)/r`. Placing its first zero at the
extrapolated radius gives

    R_c = pi / B - z0

which is `analytic_critical_radius`. Note that the critical radius depends on
the cross sections only through `c`: `Sigma_a` and `nu Sigma_f` enter only as
their difference. `k` away from criticality does depend on both.

Substituting the two approximations recovers exactly the Question 3 relations,
in mean free paths:

- classical, `D = 1/(3 Sigma_t)`, `z0 = 2D`:
  `Sigma_t R_c = pi / sqrt(3(c-1)) - 2/3` — the `'marshak'` method;
- asymptotic, `D = (c-1)|nu0|^2/Sigma_t`, `z0 = z0(c)`:
  `Sigma_t R_c = pi |nu0| - z0(c)` — the `'transport'` method.

The second is the content of part 3(d): asymptotic diffusion with the exact
`z0(c)` *is* the exact-transport criticality relation, because its whole purpose
is to reproduce the transport eigenvalue.

## Measured agreement

At `N = 400` cells, over `c = 1.02 … 2.0`:

| approximation | `R_c` numerical vs. analytic |
|---|---|
| classical | `-1.6e-04 %` to `-2.4e-04 %` |
| asymptotic | `-1.6e-04 %` to `-2.1e-04 %` |

The difference is discretisation error alone and falls as `N^{-2}`
([01](01-spherical-finite-volume.md)). The sign is consistently negative: the
finite-volume operator is slightly too leaky on a coarse mesh, so it reaches
`k = 1` at a marginally smaller radius.

Two further checks, both in the run log: `k` at the returned radius is `1.0000000000`
to ten digits, and the converged flux matches `sin(Br)/(Br)` to `3.2e-06` across
the whole mesh.

One caveat when comparing against `criticality.py` directly: `build_medium` uses
the transcendental root `compute_nu0_magnitude_numerical`, whereas the
`'transport'` method there defaults to the analytic *fit*. The two differ by
about `0.1 %` at `c = 2`, which is why the asymptotic column matches Question 3
exactly at small `c` and drifts slightly at large `c`. That is a difference in
`|nu0|`, not in the solver.
