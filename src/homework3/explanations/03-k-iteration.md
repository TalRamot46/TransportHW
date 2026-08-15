# 03 — The Two-Level `k` Iteration

**The outer loop rescales the fission source, the inner loop converges the scattering source at
fixed `k`, and two lines that look cosmetic are not.**

## The loop, in `sn.k_eigenvalue`

1. flat flux, `k = 1`;
2. **outer** — form the fission source density `nu Sigma_f phi / k`;
3. **inner** (`run_sn`) — sweep `phi <- L^-1 (Sigma_s phi + fission)` until `phi` settles;
4. `k_new = k * P_new / P_old`, with `P = sum(nu Sigma_f phi V)`;
5. renormalise by `1/max(phi)`; repeat until `|k_new - k| < 1e-9 |k_new|`.

Two of those are worth not removing:

- **Dividing the fission source by `k`** (step 2) makes a solution exist at *every* system size.
  That is what turns criticality into a root search on the size — `sn.critical_size` — instead
  of an eigensolve. Report §"Eigenvalue" makes the same point physically.
- **The renormalisation** (step 5) is not tidiness. Without it the iterate's amplitude drifts by
  a factor `k` per outer, and `|k_new - k|` becomes a difference of numbers that have lost their
  significant digits long before the tolerance is met.

Steps 4 and 5 are Assignment 1's, for the same reasons; see
`src/homework1/explanations/08-k-iteration.md`.

## What each loop costs

At the tolerances the code uses (inner `1e-8`, outer `1e-9`, `N = 8`, at a critical size):

| problem | outers | sweeps | sweeps per outer |
|---|---|---|---|
| slab, `c = 1.5` | 16 | 16 | 1.0 |
| sphere, `c = 1.5` | 35 | 35 | 1.0 |
| sphere, Pu-239 | 24 | 246 | 10.2 |

One sweep per outer is the `Sigma_s == 0` short circuit in `run_sn`
(see [02](02-sn-solver.md)). Pu-239 carries real scattering and pays about ten.

The outer count is predictable in advance, which is useful when deciding whether a slow run is a
bug. At a critical size the outer error falls by the ratio of the two lowest diffusion modes:
`c/(4c-3)` in the sphere, and `c/(9c-8)` in the slab, where `B_n = (2n-1)pi/(2 a_outer)` makes
`B_2 = 3 B_1` rather than `2 B_1`. At `c = 1.5` those are `0.500` and `0.273`, predicting
`log(1e-9)/log(ratio)` = 30 and 16 outers. Measured: 35 and 16 — the slab lands exactly, the
sphere runs a few long because `S_8` is not diffusion.

## The root search

`sn.critical_size` costs 8 evaluations of `k`. `_bracket` widens by only 3% per step
deliberately: an evaluation far from criticality is the expensive kind, since the dominance
ratio is worse there. Its tolerance is *relative* because the sizes are mean free paths in
Questions 3 and 4 and centimetres in Question 5, so no absolute `xtol` would suit both.

Which tolerance matters: the inner one is the cheap knob. Tightening it from `1e-6` to `1e-11`
triples the sweep count and does not move `k` in its first eight digits, because an inexact
inner solve is absorbed by the next outer — see [06](06-verification.md).
