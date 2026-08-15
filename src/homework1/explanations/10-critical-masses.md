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

## Sensitivity to the boundary treatment

*Historical: the Robin option these numbers came from has since been removed, and the code
now implements only the extrapolated zero — see the note in
[07](07-boundary-and-initial-conditions.md). The spread is recorded because it is a real
property of the model, even though only the left-hand column is still reproducible.*

The masses above use the extrapolated zero. Applying the same condition at the physical
surface instead ([06](06-spherical-finite-volume.md)) moved them by far more than any
numerical error, again through the cube:

| material | approximation | `M_c` extrapolated | `M_c` Robin | difference |
|---|---|---|---|---|
| Pu-239 | classical | 12.940 kg | 11.432 kg | `-11.7 %` |
| Pu-239 | asymptotic | 9.188 kg | 8.307 kg | `-9.6 %` |
| U-235 | classical | 42.345 kg | 38.679 kg | `-8.7 %` |
| U-235 | asymptotic | 32.696 kg | 30.212 kg | `-7.6 %` |

Both treatments reproduce their own analytic result to `1.5e-6`, so this 8–12 % spread is
model uncertainty, not solver error — and it is comparable to the difference between the two
diffusion approximations themselves. Any single quoted mass should be read with it in mind.

These are one-group, isotropic-scattering numbers, not predictions of real critical masses:
at two mean free paths the boundary layer is a large fraction of the system — `61.7 %` of the
neutrons produced in the Pu-239 sphere leak out of it, and `55.6 %` in U-235
([09](09-neutron-balance.md)) — and diffusion of either flavour is outside its regime. The transport benchmark radii of Sood, Forster & Parsons (2003) are the
right reference to check them against. The other three benchmark rows (H2O, Fe, Na) have
`c <= 1`, so no bare critical sphere exists for them at any radius.

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
