# 03 — The two diffusion coefficients

**Where `D0(c) = (1-c) nu0^2` comes from, and why it stays positive above `c = 1`.**

Implemented in `diffusion.py`. Classical and asymptotic diffusion solve the *same* equation

```
dn/dt - D n'' + (1-c) n = delta(x) delta(t)
```

and differ only in `D`. So the code has one solver, `_phi_diffusion`, and two coefficients.

## Classical: `D = 1/3`

`D = 1/(3 Sigma_t)`, the `1/3` being the slab-geometry correction to the pure-1D case: it is
`<mu^2> = (1/2) int_{-1}^{1} mu^2 dmu` — the mean square of the direction cosine projected on
`x`, for an isotropic flux. Purely geometric, with no reference to `c`.

## Asymptotic: `D = D0(c)/Sigma_t`, `D0(c) = (1-c) nu0^2`

The point of asymptotic diffusion is to make the diffusion solution decay at the *transport*
rate rather than the diffusion one. The exact transport equation has a discrete eigenmode
`exp(-Sigma_t |x| / nu0)`, with `nu0` the root of

```
c nu0 arctanh(1/nu0) = 1
```

The steady diffusion equation decays as `exp(-sqrt((1-c) Sigma_t / D) |x|) = exp(-sqrt((1-c)/D0) Sigma_t |x|)`.
Matching the two exponents,

```
(1-c)/D0 = 1/nu0^2      =>      D0(c) = (1-c) nu0^2
```

That is the whole derivation. Note it is fixed by the *decay length*, which is why the
asymptotic solution beats the classical one most visibly in the tails of the figures.

## The `c > 1` branch

Above `c = 1` the transcendental equation has no real root: the eigenvalue moves onto the
imaginary axis, `nu0 = i|nu0|`. Using `arctanh(i k) = i arctan(k)` the equation becomes real
again, `c |nu0| arctan(1/|nu0|) = 1` — the form tabulated by Case, de Hoffmann & Placzek.

`nu0_squared(c)` therefore returns a **signed** square: `+nu0^2` below `c = 1` and `-|nu0|^2`
above it. The sign is not cosmetic — it is what makes the single expression `D0 = (1-c) nu0^2`
correct on both sides:

```
c < 1:   D0 = (1-c)(+nu0^2)  > 0    since 1-c > 0
c > 1:   D0 = (1-c)(-|nu0|^2) > 0   since both factors flip sign
```

An imaginary `nu0` turns `exp(-x/nu0)` into `cos(x/|nu0|)` — the flux shape of a critical
system. In the *steady* problem of Assignment 1 that is the whole story. Here, in the
time-dependent problem, the medium is infinite and the multiplication shows up instead as the
growing factor `e^{(c-1)t}`; `D0` remains a perfectly ordinary positive diffusion coefficient.

## The limit at `c = 1`

There is no discrete eigenvalue at `c = 1` (`nu0 -> infinity`), so the code returns `1/3`
directly. It is a genuine limit, not a patch. Expanding the transcendental equation for large
`nu0`, `1 + 1/(3 nu0^2) + ... = 1/c`, gives `nu0^2 = c/(3(1-c))` and hence

```
D0 = (1-c) nu0^2 -> c/3 -> 1/3
```

The same expansion of the `c > 1` equation gives `|nu0|^2 = c/(3(c-1))` and again `D0 -> c/3`.
So `D0` is continuous through `c = 1` and the two approximations coincide there — which is
visible in the `c = 1` figure, where the orange and green curves lie on top of each other.

Measured, from `main.py`:

| `c` | 0.6 | 0.8 | 1.0 | 1.2 | 1.5 |
|---|---|---|---|---|---|
| `D0` | 0.48588 | 0.39629 | 0.33333 | 0.28717 | 0.23745 |

`D0` decreases monotonically through `c = 1`: more scattering per absorption means a shorter
decay length, hence a smaller effective diffusion coefficient.

## Reuse

The two root-finders are **not** reimplemented here. `homework1.exact_solution` already has
`compute_nu0_numerical` (`c < 1`) and `compute_nu0_magnitude_numerical` (`c > 1`), both
`brentq` on a properly bracketed interval, and both already checked against Case's Table 8 in
Assignment 1. They import only numpy and scipy, so nothing of homework1's plotting stack comes
with them.
