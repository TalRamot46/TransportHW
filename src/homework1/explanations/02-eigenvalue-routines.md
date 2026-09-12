# 02 — The Eigenvalue Routines

**`exact_solution.py` exposes six ways to ask for `nu0`. Picking the wrong one is silent, and
the two dispatchers do not default the same way.**

## The six entry points

Report §1 derives the eigenvalue; this is only which function computes it.

| function | branch | route |
|---|---|---|
| `compute_nu0_numerical(c)` | `c < 1` | `brentq` on `arctanh(k0)/k0 = 1/c` |
| `compute_nu0_approx(c)` | `c < 1` | the course fit `1/sqrt(1 - c^p(c))` |
| `compute_nu0(c, method=...)` | `c < 1` | dispatcher, **defaults to `'numerical'`** |
| `compute_nu0_magnitude_numerical(c)` | `c > 1` | `brentq` on `c arctan(k0) = k0` |
| `compute_nu0_magnitude_approx(c)` | `c > 1` | the same fit, whose radicand has gone negative |
| `compute_nu0_magnitude(c, method=...)` | `c > 1` | dispatcher, **defaults to `'approx'`** |

**The two dispatchers default oppositely.** `compute_nu0` gives the exact root;
`compute_nu0_magnitude` gives the fit. That asymmetry is deliberate — Q1 compares fit against
root below `c = 1`, so the *sub*critical default is the reference — but it means
`compute_nu0_magnitude(c)` with no argument silently returns fitted values. Callers that need
the root above `c = 1` ask for it explicitly: `spherical.build_medium` passes
`method='numerical'`, and `homework3.reflected.relaxation_rate` bypasses the dispatcher and
calls `compute_nu0_magnitude_numerical` by name.

## Why the `c > 1` bracket is written the way it is

`compute_nu0_magnitude_numerical` brackets `f(k0) = c arctan(k0) - k0` on `[1e-12, c pi/2]`
rather than `[0, ...]`. `f(0) = 0` is a *trivial* root, so a bracket starting at zero can be
returned instead of the physical one. The lower end is lifted off it deliberately, and the
upper end works because `f'(0) = c - 1 > 0` and `f(c pi/2) < 0`.

The tolerances (`xtol=1e-15, rtol=8.9e-16`) are at the floor because this root feeds critical
radii through `pi |nu0|`, where the error is multiplied by about three.

## `phi_transient` is the expensive call

`phi_transient` runs a `scipy.integrate.quad` per evaluation point, with the integrand cut off
by hand below `nu = 1e-12` and at `nu >= 1` where `N_nu` is zero or singular. It dominates the
Q1 runtime — `phi_asymptotic` is one exponential. If a Q1 figure feels slow, that is where it
is going, and the grid `q1.X` (500 points, seven `c` values, two methods) is the knob.

## The `c = 0` special case

`phi_asymptotic` returns zero at `c = 0` rather than calling `compute_nu0`, because there is no
discrete eigenvalue for a pure absorber. `q1._diffusion_solutions` and
`diffusion.diffusion_coefficients` each carry the matching guard — the asymptotic
approximation is skipped at `c = 0` rather than raising. Three separate places encode the same
fact; if `c = 0` ever misbehaves, expect to fix all three.
