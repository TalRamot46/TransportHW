# 01 — The planar flux at `c = 1`

**Why the Q2 integral is closed-form, and why the code uses Dawson instead of `erfi`.**

Implemented in `exact.py`. The derivation itself is in
[homework2.tex](../../../docs/homework2/homework2.tex) §2 and, in more detail, in
[explanation.md](../../../explanation.md) §2.2–2.4.

## What the source is

Paasschens' solution is the Green's function of an **instantaneous isotropic point source**,
his Eq. (4):

```
S(r, t, s_hat) = delta(r) delta(t)
```

— point in space, instantaneous in time, isotropic in direction. Three checks that it cannot be
a steady `delta(r)`:

- His abstract and §I describe "the spreading of a **pulse** of particles", with `t` "the time
  after **pulse generation**".
- His Eq. (13b) gives the uncollided term `P_0 = e^{-ct/l} delta(r - ct) / (Omega_d r^{d-1})`,
  which in `d = 3` is the first term of the assignment's formula. A delta *at* `r = vt` is the
  signature of an instantaneous release: every uncollided particle is at radius exactly `vt`. A
  steady source would give `e^{-Sigma_t r}/(4 pi r^2)`, with no `t` in it.
- The assignment asks for the solution at `t = 1, 2, 3, 4, 7, 15`, which is only a question for
  a pulse.

> **Notation trap.** Paasschens' `c` is the particle *velocity*, not the scattering ratio; `l`
> is the mean free path, `l_a` the absorption length, and `D = cl/d`. He sets `l_a -> infinity`
> — purely scattering — which is why his result is the `c = 1` case in the assignment's
> notation. His `D = cl/3 = v/(3 Sigma_t)` carries a velocity, where the `D = 1/(3 Sigma_t)`
> used here is a pure length; the difference is exactly the `1/v` on the time derivative
> discussed in [04](04-time-dependent-vs-steady.md), and `D_Paasschens = v * D_here`.

## The result being coded

Superposing the Paasschens point-source solution over the source plane gives

```
phi(x,t;1) = e^{-t}/(2t) * [ 1 + sqrt(3/pi) * int_0^{w0} sqrt(w) G(w) dw ] * Theta(t - |x|)
w0 = t (1 - x^2/t^2)^{3/4}
```

Two structural features are worth keeping in mind while reading the code:

- **The uncollided term is the `1`.** For a *plane* pulse the uncollided flux is not a delta
  but a plateau, `e^{-t}/(2t)` for `|x| < t`, because the delta in `r` is smeared over the
  plane. So the collided part is measured in units of that plateau, and the bracket is
  dimensionless.
- **`Sigma_t` cancels.** `a^{3/2} = (v Sigma_t t)^{3/2}` cancels the `Sigma_t^{3/2}` from the
  `(4 pi vt/(3 Sigma_t))^{3/2}` denominator exactly. This is the check that the substitutions
  were done right.

## Why there is no quadrature

The awkward exponents `1/8` and `3/4` in Paasschens' formula were engineered to interpolate
between `d = 2` and `d = 4`. Under `u = 1 - xi^2` followed by `w = a u^{3/4}` they collapse as
`1/6 + 1/3 = 1/2`, leaving a plain `sqrt(w)`. That `sqrt(w)` is exactly what rationalises the
interpolation `G(w) = e^w sqrt(1 + b/w)`, `b = 2.026`:

```
sqrt(w) G(w) = e^w sqrt(w + b)
```

which is elementary. `int_0^{w0} e^w sqrt(w+b) dw` integrates by parts into `erfi`.

This also disposes of the endpoint singularity flagged in
[planar_spherical_relation.md](../../../docs/homework2/planar_spherical_relation.md): the
`(t-r)^{-1/4}` behaviour is real, but it is an artifact of the variable `r`, not of the
integral. In `w` it is gone.

## Why Dawson and not `erfi`

The textbook form,

```
int_0^{w0} e^w sqrt(w+b) dw
    = e^{w0} sqrt(w0+b) - sqrt(b) - (sqrt(pi)/2) e^{-b} [ erfi(sqrt(w0+b)) - erfi(sqrt(b)) ]
```

overflows: `erfi(z) ~ e^{z^2}`, and `z^2 = w0 + b` reaches `t + b`. At `t = 15` that is already
`e^17`, and it is being multiplied by an `e^{-t}` prefactor — a large number times a small one,
which is exactly the arrangement that loses precision.

Substituting `erfi(z) = (2/sqrt(pi)) e^{z^2} D(z)`, with `D` the Dawson function
(`scipy.special.dawsn`), lets the exponential factor out:

```
int_0^{w0} e^w sqrt(w+b) dw = e^{w0} [ sqrt(w0+b) - D(sqrt(w0+b)) ] - [ sqrt(b) - D(sqrt(b)) ]
```

`collided_integral` returns this **already divided by `e^{w0}`**, so every quantity it handles
is `O(1)`. `phi_c1` then multiplies by `e^{-(t - w0)}`, which is bounded by 1 because `w0 <= t`
always. Nothing in the evaluation ever exceeds order unity, at any `t`.

## The two forms of `G`

`G(w, form=...)` and `collided_integral(w, form=...)` accept either:

- `"interp"` — Paasschens' Eq. (36b), `e^w sqrt(1+b/w)`. The default: closed form, and accurate
  to the few percent Paasschens himself quotes.
- `"series"` — the exact series `8(3w)^{-3/2} sum_N [Gamma(3N/4+3/2)/Gamma(3N/4)] w^N/N!`,
  summed in log-space via `gammaln` to avoid overflow. The same `sqrt(w)` reduction makes it
  integrable term by term, `int_0^{w0} sqrt(w) G = (8/3^{3/2}) sum_N coeff_N w0^N / N`.

The series is what the normalisation check uses, because it conserves particles *exactly*
(see [05](05-verification.md)); the interpolation is what the figures use, since it is the
solution the assignment quotes.

> The argument of the Gamma functions is `3N/4`, not `N/4` — easy to misread in the scanned
> PDF. It is forced by the exponent `3N/4 - 1` in Paasschens' Eq. (35).
