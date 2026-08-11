# 00 — Index

Explanations for the Assignment 1 spherical criticality code, Questions 4 and 5.
Questions 1–3 are documented in `docs/homework1/STATUS.md` and `plan.md`.

| # | Title | What it settles |
|---|---|---|
| [01](01-spherical-finite-volume.md) | Finite Volume in Spherical Coordinates | Why the balance form needs no special case at `r = 0`, and how the mesh is built. |
| [02](02-bell-glasstone-k-iteration.md) | The Bell & Glasstone `k` Iteration | Why `k` is updated by a ratio of fission-source integrals, and why the flux must be renormalised each sweep. |
| [03](03-boundary-and-analytic.md) | Boundary Conditions and the Analytic Radius | How one leakage coefficient covers both the extrapolated zero and the Robin condition, and what `R_c = pi/B - z0` is being compared against. |
| [04](04-asymptotic-diffusion-coefficient.md) | The Asymptotic `D` Above `c = 1` | Why `D = (1-c) nu0^2` carries over unchanged to the multiplying branch. |
| [05](05-critical-mass-results.md) | Critical Masses of U-235 and Pu-239 | The measured radii and masses, and why the two approximations differ by the amount they do. |
| [06](06-u235-data-discrepancy.md) | The Two U-235 Cross-Section Rows | Which U-235 data is used and why the prompt's row is not self-consistent. |
