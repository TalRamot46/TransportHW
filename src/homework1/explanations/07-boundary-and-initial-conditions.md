# 07 — Boundary and Initial Conditions of the `k` Problem

**The steady problem takes exactly two conditions, `phi'(0) = 0` and `phi + l0 phi' = 0`, and
no initial condition at all; the flat starting flux is a property of the solver, and the two
boundary treatments spend the *same* `z0` in two different places.**

At fixed `k` this is a second-order ODE in `r`, so it needs two conditions and has exactly
two:

1. **`phi'(0) = 0`** at the centre. Symmetry — no direction is preferred at the origin —
   or equivalently regularity, which is what selects `sin(Br)/r` over the `cos(Br)/r` that
   diverges there.
2. **`phi + l0 phi' = 0`** at the outer mesh radius, with `l0` set by the treatment below.

There is **no initial condition**, because there is no time derivative. The flat
`phi = 1`, `k = 1` that `k_eigenvalue` starts from is the first iterate of the power
iteration, not physics: any positive starting vector converges to the same `k` and the same
shape, and only the sweep count changes ([08](08-k-iteration.md)).

## What the outer condition is standing in for

The physically correct condition on a bare surface is a **vacuum condition on the angular
flux**: no neutron streams inward, `psi(R, mu) = 0` for `mu < 0`. Diffusion theory carries
only `phi` and `J`, so it cannot state a condition on half the angular range. Everything it
cannot represent is compressed into the single number `l0`:

- `l0 = 2D` (Marshak) zeroes the incoming *partial current* instead of the incoming angular
  flux;
- `l0 = z0(c) ~ 0.7104/c` (Milne) makes the asymptotic flux extrapolate to zero where the
  exact transport solution says it does.

## The two treatments: same `z0`, spent differently

> **Removed from the code.** Both treatments were implemented and measured, but neither
> Robin variant could be shown to be more accurate than the extrapolated zero against any
> reference available here, so the choice was deleted in favour of one unambiguous condition
> at `r = R`. `k_eigenvalue` now always uses the extrapolated zero. The comparison below is
> kept because it is what settles what the discarded option *was*, and it is where the
> `8–12 %` mass spread once quoted for Question 5 came from; the Robin numbers in it are no
> longer reproducible from the current code.

Both use `N` **uniform** cells of width `h = r_outer/N` — no cell is ever `h + z0` wide.
What differs is where the mesh stops and what `l0` is:

| treatment | mesh covers | `h` | `l0` | condition realised |
|---|---|---|---|---|
| `'extrapolated'` | `[0, R + z0]` | `(R+z0)/N` | `0` | `phi(R + z0) = 0`, an extrapolated zero |
| `'robin'` | `[0, R]` | `R/N` | `z0` | `phi(R) + z0 phi'(R) = 0`, a Marshak-type surface condition |

So `'robin'` is not the treatment *without* the Milne correction — both carry it.
`'extrapolated'` spends `z0` stretching the domain and then imposes plain zero flux at the
far end; `'robin'` keeps the domain physical and spends `z0` inside the boundary condition.
In the code the whole distinction is one denominator, `l0 + 0.5*h`: the distance from the
last cell centre to wherever the flux is being driven to zero.

One consequence worth remembering when reading `q4_flux_profiles.pdf`: with
`'extrapolated'` the shell `R < r < R + z0` is **fictitious**. Material properties are
continued into it and the flux plotted there is the linear extrapolation, not real flux.

## What the Robin condition actually says

`z0` is defined by the Milne problem as the distance outside the surface at which the
*asymptotic* flux extrapolates to zero, so the extrapolated treatment is that definition
applied literally. The Robin form uses the same number for a different object:

| | statement | what vanishes at `R + z0` |
|---|---|---|
| `'extrapolated'` | `phi(R + z0) = 0` | the **solution itself**, continued past the surface |
| `'robin'` | `phi(R) + z0 phi'(R) = 0` | the **tangent** to the solution at `R` |

Rearranged, the Robin condition reads `-phi(R)/phi'(R) = z0`: extend the surface tangent
until it crosses zero and require that it does so `z0` out. Since `sin(Br)/r` is concave it
lies below its tangent, so the two are genuinely different constraints. They coincide only
when the flux is straight over a distance `z0`, which is the Taylor statement

    phi(R + z0) ~ phi(R) + z0 phi'(R)

— the Robin form is the first-order truncation of the extrapolated one, valid for `B z0 << 1`.

**Where it comes from.** Not from extrapolation at all, but from a current argument. With the
`P1` partial currents `J± = phi/4 -+ D phi'/2`, a bare surface has nothing coming back, so
`J-(R) = 0`:

    phi/4 + D phi'/2 = 0    =>    phi + 2D phi' = 0    =>    J(R) = phi(R)/2

which is Marshak's condition: all the flux at the surface is outgoing, and an isotropic
outgoing half-range distribution carries a current of half the flux. Stated generally it
fixes the ratio of leakage to surface flux, `J/phi = D/l0`, i.e. zero albedo.

So the constant in it is a **linear extrapolation length** `l0`, which is not the same
quantity as the **extrapolation distance** `z0` — they agree only where the flux is linear.
Diffusion's own answer is `l0 = 2D = 0.667` mfp; transport's Milne answer is
`z0 = 0.7104` mfp. Setting `l0 = z0` in `_outer_boundary` is best read as Marshak's boundary
condition with the transport-exact constant substituted for the `P1` one.

**Which constant the code actually puts in `l0`.** `SphericalMedium.z0` is not the same
object in the two branches, and it is worth naming the distinction rather than hiding it
behind one field:

| branch | `medium.z0` holds | that branch's Marshak `2D` | gap |
|---|---|---|---|
| classical, any `c` | `2D = 0.66667` | `0.66667` | none — it *is* Marshak |
| asymptotic, `c = 1.02` | `0.69651` (Milne) | `0.65616` | `5.8 %` |
| asymptotic, `c = 1.50` | `0.47127` (Milne) | `0.47490` | `0.8 %` |
| asymptotic, `c = 2.00` | `0.34815` (Milne) | `0.36804` | `5.7 %` |

For the classical branch the `'robin'` option is therefore exactly "diffusion with the
Marshak condition applied at the surface" — which is all the Question 4 boundary comparison
uses, since `_boundary_table` builds classical media only. The asymptotic branch instead puts
Milne's `z0(c)` where a `P1` derivation would put `2D`: not `P1`-consistent, but deliberate,
since the whole point of that branch is to carry transport constants, and `z0` is the
transport answer to precisely the question `l0` is asking. It is exact as `B z0 -> 0` and
drifts with curvature, like everything else here. Only Question 5's asymptotic Robin masses
depend on the choice.

**The two are inequivalent, and not merely by a sign.** The difference changes direction with
geometry: at `c = 2` the Robin sphere is *smaller* (`1.0956` against `1.1471`) while the
Robin slab is *larger* (`a/2 = 0.412` against `0.240`). The spherical shape carries an extra
`-1/r` in its logarithmic derivative, `phi'/phi = B cot(Br) - 1/r`, so its surface tangent is
steeper and reaches zero at a smaller radius than the slab's.

`'extrapolated'` is the default because the analytic relation everything is checked against,
`R_c = pi/B - z0`, was derived by placing the first zero of the asymptotic `sin(Br)/r` at the
extrapolated surface — that *is* the extrapolated statement, so the comparison tests the
numerics rather than the boundary model ([06](06-spherical-finite-volume.md)).

## Where the conditions live in the discrete system

In finite volume the boundary conditions **modify coefficients rather than add equations** —
there are no ghost cells and no rows replaced. All `N` rows are cell balances, and both
conditions are already inside them ([06](06-spherical-finite-volume.md)):

- `phi'(0) = 0` is `conductance[0] = D * areas[0] / h = 0`, because `A_0 = 0`. Row 0 has no
  inward coupling; the geometry imposes the condition, not the code.
- the outer condition is `conductance[-1] = D * areas[-1] / (l0 + h/2)`, on the diagonal of
  the last row only — a sink with no neighbour, which is what leakage is.

Both are assembled once, before the first sweep, so every iterate satisfies them exactly.

## Two homogeneous conditions, and why the problem is still solvable

Both conditions are homogeneous and so is the equation, so `phi = 0` always satisfies them
and any solution is fixed only up to a constant. That freedom is the normalisation, not a
third condition — which is why rescaling by `1/max(phi)` each sweep changes nothing.

A non-trivial solution then exists only for particular `(R, k)`, and that eigenvalue
condition *is* the criticality relation `R_c = pi/B - z0`. Dividing the fission source by `k`
is what keeps each individual sweep an ordinary inhomogeneous two-point boundary-value
problem — `M phi = F/k` with `M` positive definite, hence solvable at any radius — so the
singular eigenvalue problem is never confronted directly.
