# 04 — The Reflected-Sphere Solver

**Three approximations, one code path: everything that distinguishes them is packed into a
`Region` and a single scalar `g`.**

## `Region` absorbs the whole theory choice

`reflected.region(material, theory)` returns the six numbers report Table 2 tabulates —
`sigma_t, c, D0, rate, mu0, z0`. After that call, nothing downstream knows or asks which
approximation it is running:

- `'classic'` returns `mu0 = 0.5` for *every* material, so `jump_ratio` gives `g = 1` and flux
  continuity falls out with no branch anywhere.
- `'asymptotic'` and `'zimmerman'` build **identical** `Region`s. They differ only in that
  `jump_ratio` returns `mu0_C/mu0_R` for the second.

`_residual` is report equation (4) transcribed literally. `_setup` exists so that
`critical_radius` and `flux_profile` cannot disagree about the pair they are solving — the
flux profile must be evaluated with exactly the `g` that produced the radius, or the interface
value is inconsistent.

## The `c = 1` branch is real, not defensive

Sodium is a pure scatterer, so `c_R = 1` **exactly** — not nearly. `nu0` diverges, `rate` is
`0`, and four expressions become `0/0`. All four limits are elementary and all four are taken
explicitly:

| function | `rate != 0` | `rate == 0` |
|---|---|---|
| `_coth_over_length` | `kappa coth(kappa L)` | `1/L` |
| `_decay` | `sinh(kappa s)/kappa` | `s` |
| `partial_current_factor` | the two log forms | `1/2` |
| `region` → `D0` | `abs(c-1)/rate^2` | `1/3` |

Letting a near-zero `rate` divide out numerically would survive iron (`rate = 0.077`) and fail
outright on sodium. `relaxation_rate` is the single place the `c > 1` / `c < 1` split appears,
and it dispatches to Assignment 1's two validated solvers — see
[05](05-assignment-1-reuse.md).

## Two small things worth not rediscovering

**`np.sinc`.** `flux_profile` writes the core shape as `np.sinc(B r / pi)`, because
`np.sinc(x) = sin(pi x)/(pi x)` supplies the value `1` at `r = 0` that `sin(Br)/(Br)` cannot.

**The bracket is exact, not a guess.** `_residual` runs from `+inf` at `R -> 0` to `-inf` at
`R = pi/B`, and `pi/B` is precisely the unreflected limit, so `(0, pi/B)` is guaranteed to
bracket the fundamental mode. `critical_radius` hands that interval straight to `brentq` with
no widening search — unlike `sn._bracket`, which has no such analytic endpoint available.

## The alternative that was rejected

`_residual` keeps the curvature term `(D_C - g D_R)/R`. The tempting alternative is to match
`-D du/dr` instead of the true current `-D dphi/dr`: `u = r phi` obeys a planar equation, and
Zimmerman's derivation is planar, so the substitution looks natural and cancels the term
outright.

It is not conservative. `-D du/dr` equals `r J + D phi`, so making *that* continuous conserves
`r J + D phi` rather than `J`, and quietly leaks particles at the interface. Zimmerman's
Eq. (9) is a statement about the physical net current.

The choice is worth up to a third of the answer, so it is recorded rather than left to be
rediscovered — Pu-239 core:

| reflector | `d` [mfp] | `J` continuous | curvature dropped |
|---|---|---|---|
| water  | 3  | 4.3321 cm | 4.1906 cm (−3.3%) |
| iron   | 10 | 4.1070 cm | 3.5037 cm (−14.7%) |
| sodium | 10 | 5.3080 cm | 3.4701 cm (−34.6%) |

The finite-volume cross-check in [06](06-verification.md) reproduces the left-hand column, as
any conservative discretisation must.
