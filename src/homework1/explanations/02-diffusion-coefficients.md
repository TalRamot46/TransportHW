# 02 — The Two Diffusion Coefficients

**Both approximations solve the same equation and differ only in `D`; the asymptotic
`D = (1-c) nu0^2` carries over unchanged above `c = 1`.**

In mean free paths (`Sigma_t = 1`) the plane-source problem is

    -D phi'' + (1 - c) phi = delta(x)

with the closed form `phi = exp(-kappa|x|) / (2 D kappa)`, `kappa = sqrt(Sigma_a/D)`.

- **Classical:** `D = 1/3`, the `P1` value.
- **Asymptotic:** `D = (1 - c) nu0^2`, fixed by requiring the closed form to reproduce
  Case's asymptotic mode `exp(-|x|/nu0) / (2(1-c)nu0)`. Matching the decay rate gives
  `kappa = 1/nu0`, and matching the amplitude then gives `D`.

## Above `c = 1`

With `nu0 = i|nu0|` both `nu0^2` and `1 - c` change sign, so

    D = (1 - c) nu0^2 = (c - 1) |nu0|^2

stays positive. This is analytic continuation, not a new definition — worth stating because
a `c`-dependent formula derived below 1 usually does *not* survive above it.

What it buys is the true relaxation length: the diffusion buckling becomes
`B = Sigma_t/|nu0|`, whereas classical diffusion has `1/B = 1/sqrt(3(c-1))`, which is only
the `c -> 1` limit of `|nu0|` — `0.8 %` high at `c = 1.02` and `34 %` high at `c = 2`. That
gap is the whole difference between the two critical radii in
[09](09-critical-masses.md).
