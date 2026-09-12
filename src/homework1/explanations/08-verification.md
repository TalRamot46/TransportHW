# 08 — Verification

**Every check that was run, with the number it produced. Each is something that would have
failed loudly if the code were wrong.**

## The eigenvalue

`compute_nu0_magnitude_numerical` reproduces every entry of Case, de Hoffmann & Placzek's
Table 8 Part II to **`4.7e-06`** — `q3.report` computes that maximum and logs it on every run,
so it is a live check rather than a recorded one. The fit `compute_nu0_approx` sits within
`5e-3 %` of the root over `c = 0.5 … 0.95`, and its continuation above `c = 1` within
`1e-4 %` at `c = 1.02` rising to `0.10 %` at `c = 2`.

## The planar solver

**Tolerance convergence.** The relative `L2` error against the closed form falls from `9e-06`
at `rtol = 1e-4` to `1.3e-12` at `rtol = 1e-11` — one for one with the request, which is what
a shooting method with no mesh should do.

**Neutron balance.** `Sigma_a * integral(phi dx) = 0.99998807` against a unit source. The
residual is the physical tail beyond `a = 10/kappa`, **not** solver error: it does not move
when `rtol` is tightened, and it does move when `n_diffusion_lengths` is raised.

**Flat error profile.** The relative error is flat across the domain rather than growing toward
`x = a`, which is the signature that the radiation condition is doing its job. A zero-flux
truncation instead gives 100% error at `x = a` that no tolerance reduces.

## The spherical solver

**Against the analytic radius.** At `N = 400` over `c = 1.02 … 2.0` the numerical critical
radius is within `1.6e-4 %` to `2.4e-4 %` of `pi/B - z0`, for both approximations.

**Mesh order.** At `c = 1.5` the relative error goes from `5.28e-04` at `N = 25` to `5.16e-07`
at `N = 800` — a factor `4.00` per doubling, i.e. clean second order. The sign is consistently
negative: the discrete operator is slightly too leaky on a coarse mesh, so it reaches `k = 1`
at a marginally smaller radius.

**Flux shape.** `k` at the returned radius is 1 to ten digits, and the flux matches
`sin(Br)/Br` to `3.2e-06`.

**Neutron balance.** `neutron_balance` returns
`|production - absorption - leakage| / production` between `6e-16` and `1.3e-13` over
`1.02 < c < 2` and both approximations.

What that does and does not test is worth being precise about. The production-to-absorption
*ratio* is **not** a test — both are `sum(phi_i V_i)` times a constant, so
`production/absorption = nu Sigma_f / (k Sigma_a)` holds whatever the flux is. (At `c = 2` this
surfaces as absorption and leakage coming out exactly equal, which is only
`nu Sigma_f / Sigma_a = 2`.) The **leakage closing the gap** is the real test: it is computed
from the boundary formula using `phi_{N-1}` alone, while the other two are interior volume
integrals, and nothing forces agreement unless the discrete operator conserves neutrons cell by
cell *and* the boundary coefficient is the exact discrete counterpart of its condition. It is
the only independent check on `conductance[-1]`.

**Dominance ratio.** The measured convergence rate matches the modal prediction `c/(4c-3)` to
six digits at every `c`, and — because `D` cancels out of it — the classical and asymptotic
runs take an *identical* sweep count despite their different radii:

| `c` | 1.02 | 1.05 | 1.10 | 1.20 | 1.50 | 2.00 |
|---|---|---|---|---|---|---|
| `c/(4c-3)` | 0.944444 | 0.875000 | 0.785714 | 0.666667 | 0.500000 | 0.400000 |
| sweeps | 280 | 133 | 79 | 50 | 32 | 25 |

The residual six-digit disagreement is that the radius used is the *numerical* critical radius,
which is not exactly critical.

## The data

`q5._mass_table` checks `Sigma_f + Sigma_c + Sigma_s == Sigma_t` to `1e-12` on all six rows it
prints, including `PROMPT_U235`. All pass — the prompt row's inconsistency is in its quoted
`c`, not in its cross-section sum ([06](06-materials-and-data.md)).
