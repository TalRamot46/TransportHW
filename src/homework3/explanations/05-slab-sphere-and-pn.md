# 05 — The Slab/Sphere Relation, and Why `P_N` and `S_{N+1}` Differ

**The slab converges onto the exact size from above and the sphere from below, and the gap
between them is exactly one extrapolation distance. In a slab `P_N` and `S_{N+1}` are the same
equations with different boundary conditions; in a sphere they are not even that.**

## `R_c = 2 (a/2) + z0`

Both geometries drive the asymptotic flux to zero at the same extrapolated surface, so
`a/2 + z0 = (pi/2) nu0` and `R_c + z0 = pi nu0` — subtracting one from twice the other leaves

    z0 = R_c - 2 (a/2)

which lets the Question 3 and Question 4 tables check each other. Case's tabulated `z0`, and
the value implied by the two S_N results:

| `c` | `z0` from `S_2` | from `S_4` | from `S_6` | from `S_10` | Case's table |
|---|---|---|---|---|---|
| 1.2 | 0.091 | 0.509 | 0.563 | 0.583 | 0.592 |
| 1.5 | 0.047 | 0.368 | 0.438 | 0.468 | 0.475 |
| 1.8 | 0.032 | 0.287 | 0.354 | 0.391 | 0.397 |

At `S_10` the implied `z0` is 1.5 % below Case's at every `c` — the same shortfall three times,
which is the relation's own asymptotic error rather than anything geometry-specific. At `S_2` it is nearly zero —
the low-order slab is far too thick and the low-order sphere far too small, so the difference
that should be a boundary layer is swallowed. That single row is the sharpest statement of what
low-order S_N gets wrong: **it is the boundary that is misrepresented, not the interior.**

## The opposite signs

`S_N` is a collocation method: it satisfies the transport equation exactly along `N`
directions and says nothing between them. In the slab the leaked current at the vacuum face is
`sum_{mu_m > 0} w_m mu_m psi_m`, and a coarse quadrature *underestimates* the leakage of a flux
that is strongly peaked towards `mu = 1` there, so the slab must be made thicker than it should
be to stay critical — the tables come down from above. In the sphere the same coarse quadrature
also mishandles the angular redistribution term, which sweeps flux towards `mu = +1` as `r`
grows; that term is a loss from the interior directions, is *over*-represented at low `N`, and
pushes the critical radius the other way. Hence `+38 %` at `S_2` in the slab against `-5.5 %`
in the sphere, at the same `c = 1.8`.

## `P_N` versus `S_{N+1}`

In *slab* geometry the two are algebraically equivalent: expanding `psi` in `N+1` Legendre
moments and truncating gives the same `N+1` coupled ODEs as collocating at the `N+1`
Gauss–Legendre points, because Gauss–Legendre quadrature of order `N+1` integrates the
Legendre polynomials up to degree `2N+1` exactly. So `P_N` and `S_{N+1}` have the same
*interior* solution, and any difference between them is a difference of boundary condition:

- `S_{N+1}` imposes `psi(a/2, mu_m) = 0` on each of the `(N+1)/2` incoming ordinates. That is
  the Mark condition, and it is what this code does.
- `P_N` needs `(N+1)/2` conditions on moments instead. Marshak's are
  `integral_0^1 P_{2j+1}(mu) psi dmu = 0`; Mark's are the same zeros as `S_{N+1}` uses.

With **Mark** conditions the two methods agree to machine precision; with **Marshak** they
differ, most at low `N`, because the moment conditions weight the whole incoming half-range
whereas the ordinate conditions only pin isolated directions. The `S_2` result of Question 3
makes this concrete: `a/2 = 0.78001` mean free paths at `c = 1.5`, and `P_1` with the Mark
extrapolation `l0 = 1/sqrt(3)` gives `arctan(1/(B l0))/B = 0.78001` — the same five digits,
from `homework1.criticality.critical_dimensions_applied_bc`. The same routine with the Marshak
`l0 = 2/3` gives `0.72348`, a 7 % difference between two methods that share every interior
equation.

In the *sphere* there is no such equivalence to begin with. The angular derivative couples the
Legendre moments differently from the way the `alpha` recursion couples the ordinates, so
`P_1` and `S_2` are genuinely different approximations there: `S_2` gives `R_c = 1.607` mean
free paths at `c = 1.5` where `P_1`-Mark gives `1.919`.
