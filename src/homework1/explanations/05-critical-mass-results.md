# 05 — Critical Masses of U-235 and Pu-239

**With the assignment's benchmark cross sections the bare critical masses come
out at 9.19 kg (Pu-239) and 32.70 kg (U-235) under asymptotic diffusion, and
about 40 % higher under classical diffusion.**

The mass follows from the radius alone, `M = rho (4/3) pi R_c^3`, so the whole
question reduces to Question 4 run with realistic cross sections. Both materials
have `Sigma_t = 0.32640 cm^-1`, a mean free path of `3.0637 cm`.

## Results

| material | `c` | approximation | `R_c` [cm] | `Sigma_t R_c` [mfp] | `M_c` [kg] |
|---|---|---|---|---|---|
| Pu-239, `rho = 15.7` | 1.50 | classical | 5.8163 | 1.8984 | 12.940 |
| Pu-239 | 1.50 | asymptotic | 5.1890 | 1.6937 | 9.188 |
| U-235, `rho = 19.0` | 1.30 | classical | 8.1031 | 2.6449 | 42.345 |
| U-235 | 1.30 | asymptotic | 7.4339 | 2.4264 | 32.696 |

Numerical and analytic radii agree to four decimal places in cm; the numbers
above are from the `k = 1` search on a 400-cell mesh.

## Why the two approximations differ so much

Entirely through the relaxation length, as set out in
[04](04-asymptotic-diffusion-coefficient.md). Classical diffusion uses
`1/B = 1/sqrt(3(c-1))`, which is the `c -> 1` limit of the transport `|nu0|` and
is too long by `12 %` at `c = 1.5` and `7 %` at `c = 1.3`. A longer relaxation
length means a larger critical sphere, and the mass goes as the cube:

- Pu-239: radius ratio `1.121`, mass ratio `1.41`;
- U-235: radius ratio `1.090`, mass ratio `1.30`.

So the classical approximation over-predicts the critical mass by 30–40 % on
these materials, and the error grows with `c` — Pu-239 at `c = 1.5` is worse
than U-235 at `c = 1.3`. The asymptotic values are the ones to quote, since
asymptotic diffusion with the exact `z0(c)` reproduces the exact-transport
criticality relation ([03](03-boundary-and-analytic.md)).

## Caveats on the physical realism

These are one-group, isotropic-scattering numbers, and should not be read as
predictions of real critical masses. Two effects dominate the residual error:

- **Transport, not diffusion, near the surface.** Both spheres are only about
  2 mean free paths in radius, so the boundary layer is a sizeable fraction of
  the system and the asymptotic flux shape is not the whole story. Diffusion
  theory of any flavour is being used far outside the regime where it is
  accurate.
- **One energy group.** The benchmark collapses the whole spectrum onto a single
  group by construction; that is what makes it an analytic benchmark, not a
  physical model.

The published transport benchmark radii in Sood, Forster & Parsons (2003) are
the right reference to check these against, and that comparison has deliberately
been left for the reader with the paper to hand rather than filled in from
memory.

The U-235 numbers above use the assignment PDF's cross sections; see
[06](06-u235-data-discrepancy.md) for the row given in the task prompt, which
gives 42.55 kg instead.
