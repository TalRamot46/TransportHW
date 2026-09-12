# 04 — The Reflected-Sphere Solver

**Two criticality equations and three theories: `Region` carries the theory, `_continuous` and
`_fixed` carry the equations, and nothing else in the module knows which is running.**

## The map to the report

| report | code |
|---|---|
| eq. (3), `k0(c_C)` and `1/nu0(c_R)` | `relaxation_rate` |
| eq. (14), `mu0(c)` | `partial_current_factor` |
| Table 1, the `D0`, `z0`, `rho_2/1` of one theory | `region`, `jump_ratio` |
| eq. (8), theories (a) and (b) | `_continuous` |
| eq. (13), theory (c) | `_fixed` |
| `a` → `R_c` in cm | `critical_radius` |
| eq. (2), both regions | `flux_profile` |

Each equation function is one line of algebra: report eq. (8) or (13) with its right-hand side
moved across, so the critical `a` is a zero. Their shared pieces are named after the terms they
are — `_far_face_return` is `coth([d+z0]/nu0)/nu0`, `_leakage` is `D_R/(D_C k0)`,
`_interface_in_reflector_mfp` is `b - d`. `_continuous` differs from `_fixed` by exactly the two
`1/r` terms of report eq. (42), and by nothing else.

Everything works in the optical radius `a` of report eq. (1); `critical_radius` divides by
`sigma_t` once, at the end, to hand back cm. The reflector's side of the interface is counted in
*reflector* mean free paths, `b - d`, not in `a` — none of these pairs shares a mean free path.

## Why the split is not symmetric

The obvious refactor is one formula for all three theories, with the curvature terms switched
off for Zimmerman. It would be wrong, and report Appendix A is why: (a) and (b) impose the exact
matching conditions of the diffusion operator, while (c) transcribes a *plane* amplitude ratio
onto a sphere and so imposes it on `r phi`. **(c) therefore does not reduce to (b) as
`mu0_C -> mu0_R`.** That is a property of the pair of theories, not a bug to be fixed.

The cost, also in Appendix A: `_fixed` does not conserve `J` at the interface. The finite-volume
cross-check of [06](06-verification.md) is conservative by construction, so it validates
`_continuous` and can say nothing about `_fixed`.

## The `c = 1` branch is real, not defensive

Sodium is a pure scatterer, so `c_R = 1` **exactly** — not nearly. `nu0` diverges, `rate` is
`0`, and four expressions become `0/0`. All four limits are elementary and all four are taken:

| function | `rate != 0` | `rate == 0` |
|---|---|---|
| `_far_face_return` | `coth(L/nu0)/nu0` | `1/L` |
| `_decay` | `nu0 sinh(s/nu0)` | `s` |
| `partial_current_factor` | the two log forms | `1/2` |
| `region` → `D0` | `abs(c-1)/rate^2` | `1/3` |

Letting a near-zero `rate` divide out numerically would survive iron (`rate = 0.077`) and fail
outright on sodium.

## Two things worth not rediscovering

**`np.sinc`.** `flux_profile` writes the core shape as `np.sinc(k0 sigma_t r / pi)`, because
`np.sinc(x) = sin(pi x)/(pi x)` supplies the value `1` at `r = 0` that `sin(k0 a)/(k0 a)` cannot.

**The bracket is exact, not a guess.** Both residuals run from `+inf` at `a -> 0` to `-inf` at
`a = pi/k0`, and `pi/k0` is precisely the unreflected limit, so `(0, pi/k0)` is guaranteed to
bracket the fundamental mode — `brentq` gets it with no widening search, unlike
`sn.core._bracket`.
