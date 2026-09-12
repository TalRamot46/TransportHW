# 03 — The Diffusion Module

**One signed square keeps `D0(c) = (1-c) nu0^2` a single expression on both sides of `c = 1`,
where the eigenvalue it is built from does not exist.**

Report §3.1 derives `D0(c)`; this is how `diffusion.py` evaluates it without branching on the
physics.

## `nu0_squared` carries the sign, not the caller

Above `c = 1` the discrete eigenvalue is imaginary, so there is no real `nu0` to square.
`nu0_squared` returns the **signed** square instead:

| `c` | returns | from |
|---|---|---|
| `< 1` | `+nu0^2` | `homework1.compute_nu0_numerical` |
| `> 1` | `-|nu0|^2` | `homework1.compute_nu0_magnitude_numerical` |
| `= 1` | `np.inf` | no discrete eigenvalue exists |

Then `(1.0 - c) * nu0_squared(c)` is positive on both sides without the caller knowing which
branch it is on — both factors flip together. This is why `diffusion_coefficient` has one
expression rather than an `if c > 1`.

**The `c = 1` case is handled before that product, not by it.** `inf * 0` is `nan`, so
`diffusion_coefficient` short-circuits to `D_CLASSICAL` at exactly `c = 1`. The `np.inf`
return is therefore never multiplied — it exists so that a caller reaching for `nu0_squared(1)`
gets an obviously wrong number rather than a plausible one. The limit from either side really
is `1/3`, which `check_diffusion_coefficients` measures at `c = 1 ± 1e-6`.

Both eigenvalue routines are called by name rather than through
`homework1.exact_solution`'s dispatchers, so this module is always on the exact root — the
`compute_nu0_magnitude` dispatcher would have defaulted to the analytic fit.

## `_phi_steady` raises above `c = 1`

The steady Green's function only exists for `c < 1`; above it the flux grows without bound and
`kappa = sqrt((1-c)/D)` would be imaginary. `_phi_steady` raises rather than returning `nan`,
which is why `check_steady_identity` in `main.py` runs over `c = 0.6, 0.8` only while every
other check runs the full `C_VALUES` including `1.2` and `1.5`.

## The `1/v` is absent on purpose

The module docstring records it: the `1/v` belongs on the time derivative whenever the unknown
is the flux `phi = v n` rather than the density `n` (report §3.2.1). With `Sigma_t = v = 1` the
two forms coincide, so it does not appear in the code. Anyone restoring physical units must put
it back — its absence here is a units choice, not a modelling one.
