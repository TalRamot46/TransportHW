# 04 — The Asymptotic `D` Above `c = 1`

**`D = (1-c) nu0^2` needs no modification on the multiplying branch: the
eigenvalue turning imaginary and the sign of `1-c` flipping cancel exactly.**

Question 2 fixed the asymptotic diffusion coefficient below `c = 1` by requiring
the diffusion Green's function to reproduce Case's asymptotic mode, giving
`D = (1-c) nu0^2` (see the header of `homework1/diffusion.py`). Question 4 needs
the same coefficient for `c > 1`, where `nu0` is imaginary.

Write `nu0 = i |nu0|`. Then `nu0^2 = -|nu0|^2`, and

    D = (1 - c) nu0^2 = (1 - c)(-|nu0|^2) = (c - 1) |nu0|^2

Both factors changed sign, so `D` stays positive — as it must. The formula is
the analytic continuation of the `c < 1` one, not a new definition, which is
worth stating explicitly because CLAUDE.md's warning applies here: a
`c`-dependent formula derived below 1 usually does *not* survive above it, and
this one does only because of the double sign flip.

`build_medium` writes it in physical units as

    D = (c - 1) |nu0(c)|^2 / Sigma_t

with `|nu0|` dimensionless (in mean free paths) from
`compute_nu0_magnitude_numerical`, the already-validated root of
`c arctan(k0) = k0`.

## What this buys

The point of the asymptotic coefficient is that the diffusion buckling then
equals the transport one:

    B^2 = Sigma_t (c - 1) / D = Sigma_t^2 / |nu0|^2   =>   1/B = |nu0| / Sigma_t

So the asymptotic diffusion equation has exactly the relaxation length of the
true transport solution, whereas classical diffusion has `1/B = 1/sqrt(3(c-1))`,
which is only its `c -> 1` limit. The gap is the reason the two approximations
give different critical radii, and it widens with `c`: `1/B` classical is `0.8 %`
above `|nu0|` at `c = 1.02` but `34 %` above it at `c = 2`.

Measured in the Question 4 sweep, the critical radii differ accordingly:

| `c` | classical `Sigma_t R_c` | asymptotic `Sigma_t R_c` | difference |
|---|---|---|---|
| 1.02 | 12.1588 | 12.0275 | `+1.1 %` |
| 1.10 | 5.0691 | 4.8729 | `+4.0 %` |
| 1.50 | 1.8984 | 1.6937 | `+12.1 %` |
| 2.00 | 1.1471 | 0.9995 | `+14.8 %` |

Classical diffusion always gives the *larger* radius, so it over-predicts the
critical mass — by a factor of about `1.5` in volume at `c = 1.5`. See
[05](05-critical-mass-results.md).
