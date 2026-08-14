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
leakage `10.744655`, residual `2.3e-14`. Roughly a third of the neutrons leak out of a
critical sphere this small — the quantitative reason the boundary treatment matters as much
as it does ([06](06-spherical-finite-volume.md)).
