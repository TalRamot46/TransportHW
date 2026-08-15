# 03 — The Two-Level Bell & Glasstone Iteration

**The outer loop rescales the fission source by the ratio of successive fission integrals; the
inner loop converges the scattering source at a fixed `k`. Their costs are set by two different
dominance ratios, and both come out where the modal estimate says they should.**

`k_eigenvalue` in `sn.py`:

1. start from a flat flux and `k = 1`;
2. **outer** — form the fission source density `nu Sigma_f phi / k`;
3. **inner** — sweep `phi <- L^-1 (Sigma_s phi + nu Sigma_f phi_outer / k)` until `phi` settles;
4. `k_new = k * integral(nu Sigma_f phi_new dV) / integral(nu Sigma_f phi dV)`;
5. renormalise by `1 / max(phi)` and repeat until `|k_new - k| < 1e-9 k_new`.

Steps 4 and 5 are Assignment 1's, and for the same reasons — see
`src/homework1/explanations/08-k-iteration.md`. Dividing the fission source by `k` makes a
solution exist at *every* size, so criticality is a root search on the size rather than an
eigensolve; renormalising keeps the iterate `O(1)` so that `|k_new - k|` does not become a
difference of numbers that have lost their significance.

## What each loop costs

Measured at the tolerances the code uses (inner `1e-8`, outer `1e-9`, `N = 8`, at a critical
size):

| problem | outers | sweeps | sweeps per outer |
|---|---|---|---|
| slab, `c = 1.5` | 16 | 16 | 1.0 |
| sphere, `c = 1.5` | 35 | 35 | 1.0 |
| sphere, Pu-239 | 24 | 246 | 10.2 |

One sweep per outer in the `c`-only problems is not an accident: they are written with
`Sigma_s = 0` (see `02`), so a single sweep inverts the transport operator exactly and the
inner loop has nothing to do. Pu-239 carries `Sigma_s / Sigma_t = 0.69` and pays about ten
sweeps per outer for it.

## The two dominance ratios

Assignment 1 derives that at a critical size the outer error falls by `c/(4c-3)` per sweep,
from the two lowest *spherical* diffusion modes `B_n = n pi / r_outer`. In a slab the modes are
`B_n = (2n-1) pi / (2 a_outer)`, so `B_2 = 3 B_1` instead of `2 B_1` and the same algebra gives
`c/(9c-8)`. At `c = 1.5` that is 0.500 for the sphere and 0.273 for the slab, predicting
`log(1e-9)/log(0.5) = 30` and `log(1e-9)/log(0.273) = 16` outers. Measured: 35 and 16. The slab
lands exactly; the sphere is a few outers long because `S_8` is not diffusion.

## Tolerances

The inner tolerance is the cheap knob. Tightening it from `1e-6` to `1e-11` triples the sweep
count (149 to 446) and does not move `k` in its first eight digits, because an inexact inner
solve is absorbed by the next outer. `1e-8` is comfortably inside that plateau.

The root search costs 8 evaluations of `k`, `brentq` on `k(size) - 1` bracketed within 3 % of
Assignment 1's exact-transport size. The bracket starts narrow deliberately: a `k` evaluation
far from criticality is the expensive kind, since the dominance ratio there is worse.
