# 07 — Critical Masses, and the Two U-235 Rows

**The bare critical masses are 9.19 kg (Pu-239) and 32.70 kg (U-235) under asymptotic
diffusion, about 40 % more under classical; the U-235 row from the task prompt is not
self-consistent and is reported separately.**

The mass follows from the radius alone, `M = rho (4/3) pi R_c^3`, so Question 5 is
Question 4 with realistic cross sections. Both materials have `Sigma_t = 0.32640 cm^-1`, a
mean free path of `3.0637 cm`, so both spheres are only about two mean free paths in radius.

| material | `c` | approximation | `R_c` [cm] | `M_c` [kg] |
|---|---|---|---|---|
| Pu-239, `rho = 15.7` | 1.50 | classical | 5.8163 | 12.940 |
| Pu-239 | 1.50 | asymptotic | 5.1890 | 9.188 |
| U-235, `rho = 19.0` | 1.30 | classical | 8.1031 | 42.345 |
| U-235 | 1.30 | asymptotic | 7.4339 | 32.696 |

The gap is entirely the relaxation length ([02](02-diffusion-coefficients.md)): classical
`1/sqrt(3(c-1))` is too long by `12 %` at `c = 1.5` and `7 %` at `c = 1.3`, and the mass goes
as the cube — mass ratios `1.41` (Pu-239) and `1.30` (U-235). The asymptotic values are the
ones to quote, since asymptotic diffusion with the exact `z0(c)` reproduces the
exact-transport criticality relation.

These are one-group, isotropic-scattering numbers, not predictions of real critical masses:
at two mean free paths the boundary layer is a large fraction of the system and diffusion of
either flavour is outside its regime. The transport benchmark radii of Sood, Forster &
Parsons (2003) are the right reference to check them against.

## The two U-235 rows

| source | `Sigma_c` | `Sigma_s` | `Sigma_t` | `c` quoted | `c` implied |
|---|---|---|---|---|---|
| Assignment PDF | 0.013056 | 0.248064 | 0.32640 | 1.30 | 1.3000 |
| Task prompt | 0.015672 | 0.180448 | 0.26140 | 1.50 | 1.3646 |

Both sum correctly to their `Sigma_t`, but only the PDF row reproduces its own quoted `c`,
so it is what `BENCHMARK['U-235']` holds. The prompt row is kept as `PROMPT_U235` and
reported alongside, giving 56.89 kg classical and 42.55 kg asymptotic — roughly `30 %` more
mass, too large a difference to leave implicit. Note that it is *heavier* despite the higher
`c`: its `Sigma_t` is lower, so the sphere is physically bigger though fewer mean free paths
across. `q5.report` checks the cross-section sum of every row it prints, so a future typo is
caught rather than propagated.
