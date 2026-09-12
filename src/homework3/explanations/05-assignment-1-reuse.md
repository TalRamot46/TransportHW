# 05 — What Comes From Assignment 1

**Six imports. The first four carry physics and are worth checking when a number looks wrong;
the rest is plumbing.**

| import | supplies | used by |
|---|---|---|
| `exact_solution.compute_nu0_numerical` (`c<1`), `compute_nu0_magnitude_numerical` (`c>1`) | Case's discrete eigenvalue, both branches | `reflected.relaxation_rate` |
| `criticality.extrapolation_distance` | `z0(c)`, the Milne fit | `reflected.region` |
| `criticality.critical_dimensions(c, 'transport-ref')` | the exact-transport reference sizes the S_N tables are measured against | `plots.reference_size`, `q5.radius_guess` |
| `spherical.build_medium`, `analytic_critical_radius` | the bare critical radii Question 1 compares against | `q1.bare_radius` |
| `materials.BENCHMARK`, `FISSILE`, `critical_mass` | the Sood cross-section rows and the mass formula | `q1`, `q5` |
| `tables.log_section`, `log_table`; `figures.*` | output helpers, so no module here carries a format string | all |

Three of these have a sharp edge:

**The exact root, not the fit.** `exact_solution` offers both
`compute_nu0_magnitude_numerical` (the transcendental root) and `compute_nu0_magnitude_approx`
(an analytic fit), and its dispatcher `compute_nu0_magnitude` defaults to **`'approx'`**.
`reflected.relaxation_rate` bypasses the dispatcher and calls the numerical routines by name,
so Question 1 is on the exact root throughout. At `c = 1.5` that is `k0 = 1.45110`, Case's
Table 8 value to all six printed digits.

**Two different naming schemes for the approximations.** `spherical.build_medium` takes
`'classical'` / `'asymptotic'`; this assignment uses `'classic'` / `'asymptotic'` /
`'zimmerman'`. `q1.BARE` is the translation table, and it deliberately maps *both* asymptotic
theories onto the same bare sphere — with no reflector there is no interface, so there is no
jump to distinguish them.

**`critical_dimensions` is indexed positionally.** It returns `(nu, z0, a/2, Sigma_t R_c)`, and
callers pick an element by integer: `q3.HALF_THICKNESS = 2` and index `3` for the radius, via
`plots.reference_size(c, index)`. Reordering that tuple in Assignment 1 would silently change
every reference column here.
