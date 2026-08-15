# 08 — The Neutron Balance Check

**Production equals absorption plus leakage to `1e-13` relative at every critical radius
tested; the non-trivial part is the leakage, the only independent check on the boundary
discretisation.**

`neutron_balance` integrates the converged critical flux three ways,

    production = sum(nu Sigma_f phi_i V_i) / k
    absorption = sum(Sigma_a phi_i V_i)
    leakage    = A_outer * D phi_{N-1} / (l0 + h/2)

and returns `|production - absorption - leakage| / production`. Over `1.02 < c < 2` and both
approximations the residual runs from `6e-16` to `1.3e-13` — machine precision.

## What it does and does not test

**Not a test:** the ratio of production to absorption. Both are `sum(phi_i V_i)` times a
constant, so `production / absorption = nu Sigma_f / (k Sigma_a)` holds whatever the flux is.
At `c = 2` that shows up as absorption and leakage coming out exactly equal (`7.597619`),
which is only `nu Sigma_f / Sigma_a = 2`.

**A real test:** that the leakage closes the gap. It is computed from the boundary formula
using `phi_{N-1}` alone, while the other two are interior volume integrals; nothing forces
them to agree unless the discrete operator conserves neutrons cell by cell and the boundary
coefficient is the exact discrete counterpart of its condition. It stays at machine precision
for both boundary treatments and both approximations.

This is the spherical analogue of `absorption_balance` in `diffusion.py`, and a stronger
check: the slab version compares against an analytic tail the truncated domain only
approximates, whereas this is an exact identity of the discrete system.

At `c = 1.5` classical, with `Sigma_a = 1`: production `32.233965`, absorption `21.489310`,
leakage `10.744655`, residual `2.3e-14`.

## The leakage fraction is fixed by the cross sections alone

At `k = 1` the tautology `A/P = Sigma_a/(nu Sigma_f)` and the balance together give

    L/P = 1 - Sigma_a/(nu Sigma_f)

with no reference to `D`, to the geometry, or to the boundary treatment. For the test media
(`Sigma_t = Sigma_a = 1`, `nu Sigma_f = c`, hence no scattering) that is `(c-1)/c`, and the
measured fractions match it to six digits — identically for both approximations, despite
their different radii.

The `c` in that expression is doing no work of its own, which matters when carrying the
number across to Question 5: the benchmark materials have the same `c` but plenty of
scattering, so their fractions come out quite differently — `61.7 %` for Pu-239 (`c = 1.50`)
and `55.6 %` for U-235, against the `33 %` of the scattering-free test medium at the same
`c = 1.5`. Two thirds of the neutrons born in a critical Pu-239 sphere never have a chance to
be absorbed in it, which is the quantitative reason the surface condition matters as much as
it does ([06](06-spherical-finite-volume.md)).
