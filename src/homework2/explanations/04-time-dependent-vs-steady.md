# 04 — Time-dependent vs. steady

**Why Q3(a)(b) needs the Green's function, and how the steady formula falls out of it.**

## The two objects

Question 3 asks for the diffusion solution to be added to the Q2 plots. Those plots are
`phi(x,t)` at `t = 1, 2, 3, 4, 7, 15`, so the curve added to them has to carry a `t`. The
familiar planar diffusion result,

```
-D phi'' + (1-c) Sigma_t phi = delta(x)
phi(x) = sqrt(3)/(2 sqrt(1-c)) exp(-sqrt(3(1-c)) Sigma_t |x|)
```

is the solution for a source that has been *on forever* — a steady state. It has no `t`, and
so cannot be overlaid on a `t`-resolved plot.

What matches the pulse source `delta(x) delta(t)` of Q1 and Q2 is the time-dependent Green's
function, which for `(1/v) dphi/dt - D phi'' + (1-c) Sigma_t phi = delta(x) delta(t)` is a
spreading Gaussian (derived in the report, `docs/homework2/homework2.tex` §3.2):

```
n_diff(x,t;c) = e^{-(1-c)t} / sqrt(4 pi D t) * exp(-x^2 / (4 D t))
```

`sqrt(2 D t)` is the diffusive width and `e^{-(1-c)t}` the absorption (or, above `c = 1`,
multiplication) of the population as a whole. The derivation runs: integrate the pulse into the
initial condition, factor the absorption out, Fourier transform, invert.

## A note on the `1/v`

The time term is `(1/v) dphi/dt` whenever the unknown is the **flux** `phi = v n`, and plain
`dn/dt` when it is the **density**. Dimensionally: `D` is a length and `Sigma_a` an inverse
length, so `D phi''` and `Sigma_a phi` both carry `phi/length`, while `dphi/dt` carries
`phi/time` — only `1/v` reconciles them. It is the same pairing of `d/dx` with `(1/v) d/dt`
that appears in the transport equation this descends from.

With `Sigma_t = v = 1` the two forms coincide, so nothing in the code changes; the general-`v`
Green's function is `phi = v e^{-(1-c) Sigma_t v t} / sqrt(4 pi D v t) * exp(-x^2/(4 D v t))`,
and setting `v = Sigma_t = 1` gives what `_phi_diffusion` computes.

## They are the same solution

The steady solution is the time-integral of the pulse solution — physically, a steady source
is a pulse source fired repeatedly, so superposing the pulse response over all elapsed times
gives the steady response. Doing the integral, with `tau = vt` and `Sigma_t = v = 1`:

```
int_0^inf  e^{-(1-c)tau} / sqrt(4 pi D tau) * exp(-x^2/(4 D tau)) dtau
    = 1/(2 sqrt(D (1-c))) * exp(-sqrt((1-c)/D) |x|)
```

At `D = 1/3` that is `sqrt(3)/(2 sqrt(1-c)) exp(-sqrt(3(1-c))|x|)` — exactly the quoted steady
formula. At `D = D0 = (1-c) nu0^2` it is `exp(-|x|/nu0) / (2(1-c) nu0)`, the asymptotic one.

So neither form is more correct than the other; they answer different questions. Both are in
`diffusion.py`:

| function | what it is |
|---|---|
| `phi_classical_diffusion(x,t,c)`, `phi_asymptotic_diffusion(x,t,c)` | the pulse response — what the figures plot |
| `phi_steady_classical(x,c)`, `phi_steady_asymptotic(x,c)` | the steady response |

and the identity above is run as a numerical check in `main.py`, agreeing to `~1e-14`
relative. That check is doing real work: it ties the coefficient in front of the Gaussian to
the coefficient in front of the exponential, so a normalisation slip in either one shows up
immediately.

The steady pair raises on `c >= 1`. That is not a limitation of the code — for a
non-absorbing or multiplying infinite medium the steady flux is genuinely infinite, since
nothing removes the particles a steady source keeps adding. Only the time-dependent form
exists there, which is another reason the figures use it.

## What the comparison then shows

Two failures of diffusion are visible in the figures, and they are separate:

- **Infinite propagation speed.** The Gaussian is non-zero for all `x`, while the exact
  solution is strictly zero beyond `|x| = vt`. The diffusion curves leak past the vertical
  line in every panel; this never goes away, it only becomes irrelevant as the bulk of the
  particles retreat from the front.
- **The wrong decay length.** Classical diffusion gets the tail slope wrong by construction;
  asymptotic diffusion is built to get it right. The gap between the orange and green curves
  in the tails *is* `D0(c)` versus `1/3`, and it is largest where `|1-c|` is largest.

Both errors shrink with `t` — at `t = 15` all three curves nearly coincide over the bulk —
because after many collisions the angular flux is nearly isotropic, which is the assumption
diffusion is built on.
