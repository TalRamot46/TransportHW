# 09 — The `c` Split Is Free

**`multiplying_medium` discards a free parameter — $\Sigma_s$ — and the critical size does not
notice. The argument is one line of operator algebra, and it is what licenses both
$\Sigma_s = 0$ in `sn.py` and the single $c/k$ source in `pn_box.py`.**

Questions 2–4 are given $c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t$ and nothing else. Every split
$\nu\Sigma_f = c - \Sigma_s$ realises the stated $c$, and those are genuinely different media
with genuinely different $k$. `sn.multiplying_medium` picks $\Sigma_s = 0$, $\nu\Sigma_f = c$
without comment. This file is why that is allowed.

## The half-thickness does not move while $k$ converges

The first thing to be clear about, because it is where the worry usually starts: the two loops
are *nested, not interleaved*. `sn.k_eigenvalue` runs at a **fixed** half-thickness — it is
baked into the solver at construction — and returns one converged number $k(\ell)$.
`slab.critical_half_thickness` only then hands that scalar function to `brentq`. So the
question is not about iteration paths at all. It is the static question: **is the zero set of
$k(\ell) - 1$ independent of the split?**

## The operator that does not know about the split

Fix a half-thickness $\ell$, let $T_\ell\psi = \mu\,\partial_x\psi + \Sigma_t\psi$ be streaming
plus total collision under the reflective/vacuum conditions, and let

$$\mathcal{K}_\ell\,\phi \;=\; \int_{-1}^{1}\Big(T_\ell^{-1}\,\tfrac{1}{2}\phi\Big)\,d\mu$$

— emit isotropically, transport once, collect the scalar flux. In the code $\mathcal{K}_\ell$
**is** `solver.sn_iteration`: it takes a source density and returns a scalar flux, and it reads
only `sigma_t` and the mesh. Same for $P_N$, where it is `BoxSystem.solve` behind a matrix built
at `pn_box.py:30` from `sigma_t` and `dx` alone. Neither $\Sigma_s$ nor $\nu\Sigma_f$ nor $k$
reaches it.

The converged fixed point of the two loops, for any split, is

$$\phi \;=\; \Sigma_s\,\mathcal{K}_\ell\phi \;+\; \frac{1}{k}\,\nu\Sigma_f\,\mathcal{K}_\ell\phi
\;=\; \Big(\Sigma_s + \frac{\nu\Sigma_f}{k}\Big)\,\mathcal{K}_\ell\phi$$

Both terms pass through the *same* $\mathcal{K}_\ell$, which is exactly what `run_sn` does at
`sn.py:85`. So the bracket is a scalar and must be the reciprocal of an eigenvalue of
$\mathcal{K}_\ell$; power iteration converges to the dominant one, $\lambda_0(\ell)$, which is
simple with a positive eigenvector (Perron–Frobenius on the sweep matrix — nonnegative and
irreducible, the fixup never firing in Questions 2–5, see [02](02-sn-solver.md)). Hence

$$\Sigma_s + \frac{\nu\Sigma_f}{k} \;=\; \frac{1}{\lambda_0(\ell)}
\qquad\Longrightarrow\qquad
k(\ell) \;=\; \frac{\nu\Sigma_f\,\lambda_0(\ell)}{1 - \Sigma_s\,\lambda_0(\ell)}$$

which is split-dependent, as it must be. But at criticality,

$$k = 1
\quad\Longleftrightarrow\quad \Sigma_s + \nu\Sigma_f = \frac{1}{\lambda_0(\ell)}
\quad\Longleftrightarrow\quad \lambda_0(\ell) = \frac{1}{c\,\Sigma_t}$$

The split survives in $k(\ell)$ and cancels at $k = 1$, because there it enters only through the
sum, and the sum *is* $c\,\Sigma_t$. The critical half-thickness is the root of an equation
$\Sigma_s$ does not appear in. Since $\lambda_0$ rises monotonically from $0$ to $1/\Sigma_t$,
that root exists and is unique exactly when $c > 1$.

Nothing in this used a continuum limit — only that $\mathcal{K}_\ell$ is linear and that both
terms go through it. So it holds for the 200-cell discrete operator too: the *computed*
half-thickness, discretisation error and all, is split-independent, which is what the report's
tables need.

## Measured

$c = 1.5$, $S_8$, 200 cells. The critical half-thickness is constant to $1.1\times10^{-11}$ —
`brentq`'s own tolerance — while $k$ at a *non-critical* half-thickness moves by a factor of two:

| $\Sigma_s$ | $\nu\Sigma_f$ | critical $a/2$ | $k$ at $a/2 = 0.5$ |
|---|---|---|---|
| 0.0 | 1.5 | 0.612235408166 | 0.922016 |
| 0.5 | 1.0 | 0.612235397229 | 0.887414 |
| 1.0 | 0.5 | 0.612235406845 | 0.797614 |
| 1.4 | 0.1 | 0.612235403122 | 0.440782 |

Stronger: taking the single value $\lambda_0(0.5) = 0.614677468$ from the $\Sigma_s = 0$ row and
predicting the other three from the formula above reproduces the measured $k$ to $10^{-9}$. So
$\mathcal{K}_\ell$ really is the same operator in all four runs, and the split really does enter
only through that scalar bracket.

## What would break it

- **A source term that bypasses $\mathcal{K}_\ell$.** The argument needs both terms to go
  through one operator. Folding $\Sigma_s$ into the removal cross section instead of the source
  — a perfectly ordinary thing to do — would leave $\mathcal{K}_\ell$ split-dependent and the
  cancellation would not happen.
- **A fixup that fires.** The clamp at `sn.py:78` is nonlinear. It never triggers in Questions
  2–5, but a problem where it did would not have a Perron eigenvalue to converge to.
- **Anisotropic scattering.** $\Sigma_s$ would then enter through more than the $n = 0$ moment
  and would no longer be summable into $c$.

## What the choice buys

Only speed and one less symbol. $\Sigma_s = 0$ puts the whole source under $k$, which is what
lets report eq. (6) be written in $c$ alone, and it makes the inner iteration exact in a single
sweep — `run_sn` returns immediately at `sn.py:89`. Questions 3–4 cost one sweep per outer;
Question 5, with real cross sections, costs about ten. See [03](03-k-iteration.md).
