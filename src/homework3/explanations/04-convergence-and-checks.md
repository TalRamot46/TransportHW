# 04 — What Was Checked, and What the Mesh Costs

**The spatial mesh is free — the critical sizes are mesh converged at 50 cells — so every
departure in the tables is angular. Three independent checks fix the sweep itself.**

## The mesh does not matter

Critical radius of the `c = 1.5` sphere at `S_10`, over four mesh refinements:

| cells | `R_c` [mfp] |
|---|---|
| 50 | 1.685945 |
| 100 | 1.685954 |
| 200 | 1.685957 |
| 800 | 1.685957 |

Seven microns of mean free path across a factor of sixteen in mesh. Diamond difference is
second order and the critical flux is a smooth cosine-like shape with no material interfaces,
so there is nothing for the mesh to resolve badly. `sphere.N_CELLS = 100` and
`slab.N_CELLS = 200` are both far inside the plateau; the sphere is the smaller of the two only
because it is the expensive geometry.

**The consequence for the tables: every number in the `departure` column of Questions 3 to 5 is
angular truncation, not spatial.**

## Check 1 — uncollided transport in the sphere

A pure absorber with a uniform unit source has the exact centre flux
`phi(0) = (1 - e^{-Sigma R})/Sigma`. At `Sigma = 1`, `R = 2` the exact value is `0.864665`, and
the code gives `0.864733` (`S_2`), `0.864645` (`S_4`), `0.864602` (`S_32`), all at 800 cells.
This exercises the areas, the volumes, the `mu = -1` starting direction and the reflection at
`r = 0`, with no scattering to hide behind.

## Check 2 — the exact slab benchmark

The one-speed critical slab at `c = 1.5` has half-thickness `0.605055` mean free paths — the
Sood et al. Pu-239 slab, `1.853722` cm at `Sigma_t = 0.32640` cm^-1, which is the same
cross-section table Assignment 1 Question 5 uses. The S_N sequence walks down onto it:

| order | `a/2` [mfp] | departure |
|---|---|---|
| `S_10` | 0.609042 | +0.66 % |
| `S_16` | 0.606407 | +0.22 % |
| `S_24` | 0.605631 | +0.10 % |
| `S_32` | 0.605373 | +0.05 % |
| `S_48` | 0.605195 | +0.02 % |

Monotone, and roughly second order in `1/N`.

## Check 3 — Assignment 1's exact-transport relation

`critical_dimensions(c, 'transport-ref')` gives `a/2 = (pi/2) nu0 - z0` and
`R_c = pi nu0 - z0` from Case's tabulated `nu0` and `z0`. At `S_10` the S_N results sit within
0.4 % of it in the slab and 0.3 % in the sphere, at every `c` — and they approach it from
*opposite* sides, which is why the two geometries are worth running together. See `05`.
