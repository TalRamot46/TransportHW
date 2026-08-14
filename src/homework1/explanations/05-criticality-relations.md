# 04 — The Criticality Relations and the Sign of `q`

**All five methods evaluate `a/2 = (pi/2) nu - z0` and `Sigma_t R_c = pi nu - z0`; they
differ only in where `nu` and `z0` come from. Case's Table 23 says the printed `q = -0.0199`
has the wrong sign.**

For `c > 1` the asymptotic flux is trigonometric ([01](01-case-eigenvalue.md)):
`cos(x/|nu0|)` in the slab (the sine is killed by symmetry) and `sin(r/|nu0|)/r` in the
sphere (the cosine by regularity). A bare system is critical when that shape reaches zero
at the extrapolation distance `z0` outside the surface, which gives the two relations
above — the quarter-wave zero of the cosine and the first zero of the sine.

| method | `nu` | `z0` | part |
|---|---|---|---|
| `transport` | fitted `\|nu0\|` | fitted `z0(c)` | 3(a) |
| `transport-ref` | transcendental root | Case's Table 23 | reference for 3(a) |
| `marshak` | `1/sqrt(3(c-1))` | `2/3` | 3(b) |
| `mark` | `1/sqrt(3(c-1))` | `1/sqrt(3)` | 3(c) |

`transport` uses no diffusion coefficient at all, but it still neglects the transient flux
near the surface, so it is asymptotic rather than exact — least defensible at `c = 2`,
where the critical sphere is barely one mean free path in radius.

Marshak sets the incoming partial current to zero, `phi + 2D phi' = 0`, so `z0 = 2D = 2/3`;
Mark sets the incoming angular flux to zero along `mu = 1/sqrt(3)`, giving `z0 = 1/sqrt(3)`.
Both are constants, whereas the transport `z0(c) ~ 0.7104/c` halves across the range — that,
plus `1/B` being only the `c -> 1` limit of `|nu0|`, is the whole story of the comparison.

## The sign of the quadratic correction

The expansion about `c = 1` is `z0(c) = 0.710446 [1 + q (1-c)^2] / c`, whose leading term is
exact in the sense that `0.710446/c = 0.710446 nu0 arctanh(1/nu0)`. The notes print
`q = -0.0199`, but the tabulated product `c z0` rises away from `0.710446` on **both** sides
of `c = 1`, which a negative `q` cannot do. With `q = +0.0199` the fit matches Table 23 to
its four printed digits at `c = 0.9` and `c = 1.1`, so the printed minus sign looks like a
transcription error. Measured maximum error against Table 23 over `1.02 < c < 2`:
`2.63 %` with `q = -0.0199`, `1.33 %` with `q = +0.0199`. Results are quoted with the
printed sign and both are plotted; below `c = 1.2` the choice moves `a/2` by under `0.07 %`.

## Extrapolated zero vs. applying the condition

Setting the flux to zero at `s + z0` is not the same as imposing `phi + l0 phi' = 0` on the
shape, which gives `B a/2 = arctan(1/(B l0))` and `u cot u = 1 - u/(B l0)`. The two agree to
first order in `B l0` (`0.1 %` at `c = 1.02`) and separate as the system shrinks — a factor
`1.7` in the Marshak slab at `c = 2`. That is what shows the good large-`c` agreement of the
extrapolated Marshak slab to be a cancellation of two errors, not accuracy.
