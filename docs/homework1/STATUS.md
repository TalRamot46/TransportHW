# Homework 1 — Implementation Status

Status of `src/homework1` against `instruction_files/Assignment1.pdf`, as of 2026-08-07.

## Summary

| Question | Status |
|---|---|
| 1a — asymptotic / transient components + relative contributions | Done |
| 1b — classical diffusion | Done |
| 1c — asymptotic diffusion | Done |
| 1d — relative error of each diffusion approximation | Done |
| 1 — all seven required `c` values plotted and regenerated | Done |
| 2 — numerical diffusion code with delta-source | Done (branch `homework1-q2-diffusion`) |
| 3 — critical `a/2` and `R_c` for `1.02 < c < 2` (5 methods) | Not started |
| 4 — spherical diffusion code, `k = 1` critical radius | Not started |
| 5 — bare critical mass of U-235 and Pu-239 | Not started |

## Question 1 — exact solution (Case's method)

Implemented in `src/homework1/exact_solution.py` and `src/homework1/diffusion.py`:

- Discrete eigenvalue `nu0` via root of `arctanh(k)/k = 1/c`, plus the closed-form
  approximation, with a printed comparison table (`main.py:36-46`).
- `phi_as`, `phi_tr` (quadrature over `nu` in `(0,1)`), and `phi = phi_as + phi_tr`.
- Classical diffusion (`diffusion.py:6-26`) and asymptotic diffusion (`diffusion.py:28-50`).
- All four plots: flux components, relative contributions, four-way comparison,
  relative errors (`plots.py`).

Question 1 is complete. Three fixes were applied on 2026-08-07:

- **`c = 0` restored** in `main.py:28` (was `1e-1`), so the case list now matches the
  assignment exactly. No solver change was needed — `phi_as` already returns 0 and
  `phi_tr` already reduces to the pure-absorber form `(1/2) E_1(|x|)` at `c = 0`.
- **The four Q1 plot calls were uncommented** in `main.py`, so `main()` regenerates all
  eight Q1 figures (four plots x two `nu0` methods) instead of only the Q2 pair.
- **Axis-label bug in the panel grid fixed.** The grid is built with `sharex=True`, which
  hides x tick labels on every panel that is not at the bottom of its column. With seven
  `c` values the eighth panel is deleted, which exposed the `c = 0.9` panel at the bottom
  of the right column with its tick labels still hidden and no x label. The grid setup and
  teardown now live in the `_make_grid` / `_finalize_grid` helpers in `plots.py`, shared by
  all four Q1 figures, which re-enable `labelbottom` on the exposed panels and derive the
  row count and panel indices from `len(c_values)` instead of hard-coding 4x2 and index 7.

## Question 2 — numerical diffusion code

`solve_diffusion_numerical` (`diffusion.py:52-103`) reduces `phi'' = ((1-c)/D) phi` to a
first-order system and shoots on `phi'(-a)`. The delta source is handled by solving the
half-domain `x` in `[-a, 0]` and imposing the current boundary condition
`phi'(0) = 1/(2D)` (that is, `-D phi' = 1/2` of the unit source) rather than smearing the
delta over a cell. Compared against the analytic solution with an error plot
(`plots.py:225`).

All three gaps are now closed on branch `homework1-q2-diffusion` — see `plan.md` for the
derivation and measured results. Summary:

- **One solver**, `solve_diffusion_shooting`: half-domain, using the symmetry-derived
  current condition `J(0+) = 1/2` at the source and a radiation condition at the outer
  boundary, with a single scaled integration instead of root-finding. Maximum relative
  error `2.0e-8 %`.
- **Both approximations** are covered by that one solver, since asymptotic diffusion is the
  same equation with `D = (1-c) nu0^2`.
- **Convergence study** added over the integrator tolerance — the meaningful analogue of
  mesh refinement for a shooting method. Error falls from `9e-6` at `rtol = 1e-4` to
  `1.3e-12` at `rtol = 1e-11`.
- The zero-flux outer boundary was the real defect in the original solver, forcing a `100 %`
  relative error at `x = -a`. The radiation condition removes it; the error is now flat
  across the domain.
- Finite-volume solvers were implemented, measured, and then removed in favour of keeping
  only the shooting method. Their numbers are preserved in `plan.md` section 5.2, and the
  code is recoverable from this branch's history — likely wanted for Q4, whose Bell &
  Glasstone `k` iteration operates on a discretised mesh operator.
- `solve_diffusion_numerical` is retained, with a docstring explaining its limitation, so
  the original Q2 figure still builds.
- The stale "3 decay lengths" comment is fixed.

## Questions 3, 4, 5 — not started

Nothing in `src/homework1` covers these; `docs/homework1/homework1.tex:237-251` holds
placeholder stubs.

- **Q3.** Needs a supercritical branch of the eigenvalue that the current code cannot
  produce: `compute_nu0_numerical` raises for `c >= 1`, and for `c > 1` the root moves off
  the real axis (`1/c = arctan(eta)/eta` with `nu0 = i/eta`). Also missing: the
  extrapolation distance `z0(c)` / `l0(c)`, Marshak and Mark boundary conditions, and the
  planar / spherical critical-size solves.
- **Q4.** No spherical diffusion solver and no Bell & Glasstone power iteration for `k`.
  `src/homework4/criticality.py` finds a critical radius, but by Monte Carlo bisection for
  Assignment 4 — a different method, not reusable beyond the bracketing idea.
- **Q5.** No cross-section table in the repo and no mass calculation. Gated on Q3/Q4
  producing `R_c` first.

## Packaging note

`pyproject.toml` still declares `package-dir = {"" = "code"}`, but the sources now live in
`src/`. Until that is updated, `homework1` cannot be imported from an install of the
project; the figures above were regenerated with `PYTHONPATH=src`.

## Suggested order of work

1. ~~Restore `c = 0` and uncomment the Q1 plot calls.~~ Done 2026-08-07.
2. ~~Add the Q2 grid-convergence sweep and the asymptotic-diffusion variant.~~ Done
   2026-08-07 on branch `homework1-q2-diffusion`.
3. Extend the eigenvalue solver to `c > 1`, then build Q3 on top of it. The finite-volume
   solver added for Q2 is the intended foundation for Q4.
4. Q4 spherical `k`-eigenvalue code, validated against Q3's analytic `R_c`.
5. Q5 critical masses, using the Q3/Q4 radii and the Sood benchmark cross sections.
