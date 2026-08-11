# 07 — The Neutron Balance Check

**Production equals absorption plus leakage to `1e-13` relative at every
critical radius tested; the non-trivial part is the leakage term, which is the
only independent check on the boundary discretisation.**

`neutron_balance` integrates the converged critical flux three ways:

    production = sum(nu Sigma_f phi_i V_i) / k
    absorption = sum(Sigma_a phi_i V_i)
    leakage    = A_outer * D phi_{N-1} / (l0 + h/2)

and returns the relative residual `|production - absorption - leakage| /
production`. Measured over `1.02 < c < 2` for both approximations, the residual
ranges from `6e-16` to `1.3e-13` — machine precision.

## What this actually tests, and what it does not

It is worth being precise, because the check looks stronger than it is in one
respect and is genuinely strong in another.

**Not a test:** the ratio of production to absorption. Both are the same integral
`sum(phi_i V_i)` multiplied by a constant, so
`production / absorption = nu Sigma_f / (k Sigma_a)` holds identically, whatever
the flux is. At `c = 2` this shows up as absorption and leakage coming out
exactly equal (`7.597619` each), which is just `nu Sigma_f / Sigma_a = 2`.

**A real test:** that the *leakage* closes the gap. Leakage is computed from the
boundary formula at the outer face, using only `phi_{N-1}`, while the other two
terms are volume integrals over the interior. Nothing forces them to agree
unless the discrete operator conserves neutrons cell by cell and the boundary
coefficient is the exact discrete counterpart of the condition it represents.
That the residual is at machine precision, and stays there for both boundary
treatments and both approximations, is what confirms it.

This is the spherical analogue of `absorption_balance` in
`homework1/diffusion.py`, which checks the Question 2 slab solver against its
unit source, and it is a stronger check here: the slab version compares against
an analytic tail that the truncated domain only approximates, whereas this one
is an exact identity of the discrete system.

## Sample values

At `c = 1.5`, classical, in mean free paths with `Sigma_a = 1`:

| quantity | value |
|---|---|
| production | 32.233965 |
| absorption | 21.489310 |
| leakage | 10.744655 |
| residual | `2.3e-14` |

Roughly a third of the neutrons leak out of a critical sphere this small, which
is the quantitative statement of why the boundary condition matters so much to
the answer ([03](03-boundary-and-analytic.md)).
