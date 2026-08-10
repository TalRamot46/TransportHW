# 02 — Scaling to general `c`

**How the Q1 identity carries the `c = 1` solution to every other `c`, including `c > 1`.**

Implemented as `phi_exact` in `exact.py` — three lines, because all the work was done in Q1.

## The identity

Question 1 proves, for the angular flux,

```
psi(x,mu,t;c) = c e^{-(1-c) v Sigma_t t} psi(cx, mu, ct; 1)
```

Integrating both sides over `mu` leaves the relation untouched — the factor is independent of
`mu` — so the scalar flux obeys the same thing:

```
phi(x,t;c) = c e^{-(1-c)t} phi(cx, ct; 1)
```

So the planar Green's function is only ever constructed once, at `c = 1`, and every other `c`
is a rescaling of it. That is why `exact.py` has one real function and one wrapper.

## Reading the three factors

Each factor does one job, and knowing which is which makes the code obvious:

- `e^{-(1-c)t}` — the removal deficit. `psi(cx,ct;1)` is removed at rate `c Sigma_t` rather than
  `Sigma_t`; the exponential makes up the difference `(1-c) Sigma_t`.
- `x -> cx`, `t -> ct` — the rescaling that produces that altered removal rate in the first
  place.
- the prefactor `c` — **not** from the PDE. The equation is linear and homogeneous, so it can
  never fix an overall constant; `c` is the Jacobian of `x -> cx` acting on the source delta,
  `delta(cx) = delta(x)/c`. It is fixed by the initial condition, not the operator.

## Consequences the code relies on

**The support scales correctly.** `phi(cx, ct; 1)` is non-zero for `|cx| < ct`, i.e. `|x| < t`.
The causal front sits at `|x| = vt` for every `c`, which is right: `c` changes what happens to
a particle at a collision, not how fast it flies.

**`c > 1` needs no special handling.** For a multiplying medium `1 - c < 0`, so `e^{-(1-c)t}`
becomes exponential *growth* and `x -> cx` compresses rather than stretches. Both are just
arithmetic; nothing in the formula breaks. This is worth stating because the *diffusion*
coefficient does need a separate branch above `c = 1` (see [03](03-diffusion-coefficients.md)) —
the exact solution does not.

**Particle balance follows for free.** Since `int phi(x,t;1) dx = 1`,

```
int phi(x,t;c) dx = c e^{-(1-c)t} * (1/c) * 1 = e^{-(1-c)t}
```

the `1/c` coming from the change of variable in the integral. It cancels the prefactor exactly.
This is check 1 in [05](05-verification.md), and it is a genuinely independent test of the
prefactor: get the Jacobian wrong and the normalisation is off by a factor of `c`.
