# 01 — The Module Map

**Five question modules that only format output, sitting on four solver modules that carry all
the physics and never print anything.**

## The files

| file | owns |
|---|---|
| `main.py` | Entry point. Creates the figure directory, then calls `report(figs)` on `q1`…`q5` in order. |
| `q1.py` … `q5.py` | One question each: which cases to run, what to tabulate, what to plot. No physics. |
| `exact_solution.py` | Case's discrete eigenvalue and the exact plane-source flux (asymptotic + transient). |
| `diffusion.py` | Planar diffusion: the two `D` choices, the Green's function, the shooting solver. |
| `criticality.py` | The Q3 relations, `z0(c)`, and Case's tabulated data as arrays. |
| `spherical.py` | The spherical finite-volume operator, the `k` iteration, the critical-radius search. |
| `materials.py` | The Sood benchmark cross-section rows and the mass formula. |
| `figures.py`, `tables.py` | Matplotlib helpers and logged tables, so no other module carries a format string. |

The `q*.py` / solver split is strict, and it is what the later assignments depend on:
homework3 imports `criticality`, `exact_solution`, `spherical` and `materials` and never
touches a `q*` module. Anything added to a `q*` file is, by that convention, unreusable.

## The call paths

Questions 1 and 2, planar:

    exact_solution.phi_exact          asymptotic + transient components
    diffusion.solve_diffusion_shooting  one solve_ivp pass, rescaled
      <- q1.report / q2.report        tables + figures

Questions 4 and 5, spherical:

    spherical.build_medium            (sigma_t, sigma_a, nu_sigma_f, approximation) -> SphericalMedium
      <- spherical._banded_operator   the tridiagonal removal matrix, assembled once
      <- spherical.k_eigenvalue       source iteration; returns KResult(k, r, phi, sweeps)
      <- spherical.critical_radius    brentq on k(R) - 1
      <- q4.report / q5.report        tables + figures

Question 3 has no solver at all — `criticality.critical_dimensions(c, method)` is a closed-form
evaluation, and `METHODS` is the whole of it: a dict mapping each method name onto a
`(relaxation-length source, z0 source)` pair. Adding a sixth method means adding one row there,
not a code path.

## Two things to know before editing

**`build_medium` is the only place an approximation name is interpreted.** It takes
`'classical'` or `'asymptotic'` and returns a frozen `SphericalMedium` holding `D`,
`sigma_a`, `nu_sigma_f`, `z0` and `c`. Everything downstream reads those five numbers and
never asks which branch produced them.

**`critical_dimensions` returns a positional tuple**, `(nu, z0, a/2, Sigma_t R_c)`, and callers
index it by integer — `q3.HALF, q3.RADIUS = 2, 3`, and homework3 uses the same indices.
Reordering that tuple would silently change results in two assignments.
