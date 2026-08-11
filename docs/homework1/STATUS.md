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
| 3 — critical `a/2` and `R_c` for `1.02 < c < 2` (5 methods) | 3a–3c done; 3d–3e not started |
| 4 — spherical diffusion code, `k = 1` critical radius | Done |
| 5 — bare critical mass of U-235 and Pu-239 | Done |

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

## Question 3 — critical dimensions

Parts 3(a)–3(c) are implemented in `src/homework1/criticality.py`. All three evaluate the
same pair of relations and differ only in where the two inputs come from, so they are
selected by a `method` key (`METHODS`): `'transport'` (3a), `'marshak'` (3b), `'mark'`
(3c), plus `'transport-ref'` and `'transport-q+'` for the reference evaluations of 3a.

### 3(a) — exact transport

    a/2 = (pi/2)|nu0(c)| - z0(c),    Sigma_t R_c = pi |nu0(c)| - z0(c)

both in mean free paths, with `|nu0(c)|` and `z0(c)` from the approximate formulas.

- **Multiplying branch of the eigenvalue** added to `exact_solution.py`. For `c > 1` the
  root moves onto the imaginary axis; with `nu0 = i|nu0|` and `k0 = 1/|nu0|` the
  transcendental equation becomes `c arctan(k0) = k0`. Both a root-finder
  (`compute_nu0_magnitude_numerical`) and the continuation of the Question 1 fit,
  `|nu0| = 1/sqrt(c^p(c) - 1)` (`compute_nu0_magnitude_approx`), are provided. The fit is
  within `1e-4 %` at `c = 1.02` and `0.10 %` at `c = 2`; the root reproduces Case's
  Table 8 Part II to `4.7e-6`.
- **Extrapolation distance** from the two-term expansion about `c = 1`. The course notes
  print the quadratic coefficient as `q = -0.0199`, but Case's Table 23 (`c z0` rising on
  *both* sides of `c = 1`) requires `q = +0.0199`, which then matches the table to its
  four printed digits at `c = 0.9` and `c = 1.1`. Both signs are computed and plotted;
  the difference is `< 0.07 %` in `a/2` below `c = 1.2` but `3.2 %` at `c = 2`.
- **Results** span `a/2 = 5.67 -> 0.33` and `Sigma_t R_c = 12.03 -> 1.00` mean free paths
  over `1.02 < c < 2`. The error of the approximate route is dominated by `z0`, not by
  `nu0`, and stays under `0.1 %` up to `c = 1.25`.
- Figures `q3_critical_dimensions.pdf` and `q3_extrapolation_distance.pdf`.

### 3(b), 3(c) — classical diffusion with Marshak and Mark conditions

Classical diffusion in a multiplying medium gives `phi'' + B^2 phi = 0` with
`B^2 = (c-1)/D = 3(c-1)`, i.e. the *same* flux shapes as transport. The criticality
relations are therefore unchanged in form; only the inputs move:

    |nu0(c)| -> 1/B = 1/sqrt(3(c-1)),    z0 -> 2/3 (Marshak) or 1/sqrt(3) (Mark)

- Neither diffusion input carries the true `c`-dependence, which is the whole story of the
  comparison: `1/B` is the leading term of `|nu0|` as `c -> 1+` (0.8 % high at `c = 1.02`,
  34 % high at `c = 2`), and the two extrapolation distances are constants while the
  transport `z0(c)` falls as `~0.7104/c`.
- The two errors act in opposite directions and carry different weights in the two
  geometries (`pi/2` vs `pi` on the relaxation length, the same `z0` subtracted in both).
  Marshak's `a/2` is therefore non-monotonic in error — `+5.5 %` at `c ~ 1.23`, zero near
  `c = 1.52`, `-26 %` at `c = 2` — while both spherical radii are over-estimated
  throughout, `+15 %` (Marshak) and `+24 %` (Mark) at `c = 2`. Mark's `a/2` ends at `+1 %`
  at `c = 2`, which is a cancellation, not accuracy.
- `critical_dimensions_applied_bc` solves `phi + l0 phi' = 0` on the flux shape instead of
  using an extrapolated zero, giving `B a/2 = arctan(1/(B l0))` and
  `u cot u = 1 - u/(B l0)`. The two agree to first order in `B l0` (0.1 % at `c = 1.02`)
  and diverge as the system shrinks (factor 1.7 in the Marshak slab at `c = 2`), which is
  what shows the good large-`c` agreement above to be a cancellation.
- Figure `q3_method_comparison.pdf`.

Still open: parts **3(d)–3(e)** — the asymptotic diffusion approximation with the exact
`l0(c)`/`z0(c)`, and with `l0(c)` from a modified Marshak-like boundary condition. The
linear extrapolation length `l0(c)` is not implemented, and the spherical relation used
here omits the Winslow curvature correction.

## Question 4 — spherical `k`-eigenvalue code

Implemented in `src/homework1/spherical.py`, with figures in
`src/homework1/plots_spherical.py`. Derivations are in `src/homework1/explanations/`.

- **Finite volume on cell-centred spherical shells.** Writing the discretisation as a
  neutron balance rather than as a difference of derivatives removes the `1/r^2`
  singularity at the origin: the innermost face has zero area, so `phi'(0) = 0` is imposed
  exactly with no special case. Tridiagonal, solved in `O(N)` by `solve_banded`.
- **One leakage coefficient, `D A / (l0 + h/2)`,** covers both the extrapolated zero
  (`l0 = 0` on a mesh run out to `R + z0`) and the Robin condition (`l0 = z0` at `R`)
  without a branch. The extrapolated form is the default, since it is what the Question 3
  relations assume.
- **Bell & Glasstone source iteration** (pp. 189–192) for `k`, updated by the ratio of
  successive fission-source integrals, with the flux renormalised each sweep. The critical
  radius follows from `brentq` on `k(R) - 1`.
- **Both approximations** come from the same solver: classical `D = 1/(3 Sigma_t)`,
  `z0 = 2D`; asymptotic `D = (c-1)|nu0|^2/Sigma_t`, `z0 = z0(c)`. The asymptotic
  coefficient is the analytic continuation of Q2's `D = (1-c) nu0^2` — both `nu0^2` and
  `1-c` flip sign above `c = 1`, so `D` stays positive.
- **Agreement with analytic** `R_c = pi/B - z0` is `1.6e-4 %` to `2.4e-4 %` over
  `1.02 < c < 2` at `N = 400`, for both approximations. Second-order convergence measured
  at exactly `4.00` per doubling, `5.28e-04` at `N = 25` to `5.16e-07` at `N = 800`.
  `k` at the returned radius is 1 to ten digits and the flux matches `sin(Br)/Br` to
  `3.2e-6`.
- Figures `q4_critical_radius.pdf`, `q4_mesh_convergence.pdf`, `q4_flux_profiles.pdf`.

## Question 5 — bare critical masses

`src/homework1/materials.py` holds the Sood benchmark table and the mass calculation.

| Material | Approximation | `R_c` [cm] | `M_c` [kg] |
|---|---|---|---|
| Pu-239 (`c = 1.50`) | classical | 5.8163 | 12.940 |
| Pu-239 | asymptotic | 5.1890 | 9.188 |
| U-235 (`c = 1.30`) | classical | 8.1031 | 42.345 |
| U-235 | asymptotic | 7.4339 | 32.696 |

Classical diffusion over-predicts the mass by 30–40 %, because its relaxation length
`1/sqrt(3(c-1))` is only the `c -> 1` limit of `|nu0(c)|` and the mass goes as the cube of
the radius. Figure `q5_criticality.pdf`.

The U-235 row supplied in the task prompt (`Sigma_t = 0.26140`, `c` quoted as 1.50) is not
self-consistent — its own cross sections give `c = 1.365` — and is not the assignment's
row. The assignment values are used; the prompt row is retained as `PROMPT_U235` and
reported alongside, giving 56.89 / 42.55 kg instead. See
`explanations/06-u235-data-discrepancy.md`.

## Packaging note

`pyproject.toml` now declares `package-dir = {"" = "src"}` with `homework1` among its
packages, so `.\.venv\Scripts\python.exe -m homework1.main` runs from the repository root
as documented in `CLAUDE.md`. The earlier `PYTHONPATH=src` workaround is no longer needed.

## Building the report

Use `.\docs\build.ps1 homework1`, not `latexmk`.

On this machine a file that has just been written cannot immediately be overwritten; the
write fails with a Windows EINVAL, reported by pdflatex as ``! I can't write on file
`homework1.pdf'.`` `latexmk` runs its pdflatex passes back to back, so its second pass
always collides with the output of its first. Invoking `latexmk.pl` through `perl` directly
fails the same way — `latexmk` is not the cause. A single `pdflatex` pass succeeds, which is
why this only shows up when a document is *rebuilt*.

`docs/build.ps1` deletes the output before each pass, which avoids the collision. The same
workaround is applied to the figures in `_savefig` (`src/homework1/plots.py`). Ruled out as
causes: OneDrive (not on this path), Windows Defender Controlled Folder Access (disabled),
third-party antivirus (none registered), and stale file locks (an exclusive open succeeds).

## Suggested order of work

1. ~~Restore `c = 0` and uncomment the Q1 plot calls.~~ Done 2026-08-07.
2. ~~Add the Q2 grid-convergence sweep and the asymptotic-diffusion variant.~~ Done
   2026-08-07 on branch `homework1-q2-diffusion`.
3. ~~Extend the eigenvalue solver to `c > 1`, then build Q3 on top of it.~~ Done for parts
   3(a)–3(c) on 2026-08-08; parts 3(d)–3(e) need `l0(c)`. The finite-volume solver added
   for Q2 is the intended foundation for Q4.
4. ~~Q4 spherical `k`-eigenvalue code, validated against Q3's analytic `R_c`.~~ Done.
5. ~~Q5 critical masses, using the Q3/Q4 radii and the Sood benchmark cross sections.~~
   Done.
6. Remaining: parts **3(d)–3(e)**. Note that 3(d) is now partly covered — Q4's asymptotic
   mode is exactly asymptotic diffusion with the exact `z0(c)` — but in spherical geometry
   only, and `l0(c)` is still not implemented.
