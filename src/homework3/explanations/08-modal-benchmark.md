# 08 — The Modal Benchmark

**`pn_modal.py` solves the same P_N system with no mesh at all, which is the only reason the
Δx² error of [07](07-pn-box-solver.md) can be quoted as a number rather than estimated.**

Method 2 is 60 lines and runs in microseconds, so `q2.py` calls it inside the table loop with
nothing cached. Its value is entirely as a second opinion: it shares `pn.py` with Method 1 and
nothing else — no discretisation, no iteration, no root search in `k`.

## Three places it goes wrong if written naively

**The sign of `K^2`.** The elimination produces `AB Phi'' = Sigma_0 Phi`, so the Helmholtz form
carries a minus: `K^2 = -(AB)^-1 Sigma_0`, report eq. (23). Dropping it turns every `cos` into a
`cosh`, and `det H(a/2)` then has no zero at all — the failure is loud, but only if you know the
sign was the suspect.

**Complex round-off from `eig`.** `K^2` is not symmetric as stored, though it is similar to a
symmetric matrix, so `np.linalg.eig` can return eigenvalues with imaginary parts of order
`1e-17`. Those are taken as round-off and discarded with `.real` before anything compares them
to zero — including `values.max()`, which raises on a complex array.

**Overflow in the boundary-layer columns.** `kappa_j a/2` reaches only about 6 for the cases in
the report, but the same expression at large `c` or high `N` overflows `cosh`. Every hyperbolic
column is divided by its own `cosh(kappa_j a/2)`, leaving `v_j - (B v_j) kappa_j tanh(kappa_j a/2)`,
which is bounded for all `a/2`. Scaling a column multiplies the determinant by a non-zero factor
and cannot move its zeros, so the root is untouched.

## Why the root is found by scanning first

`_first_root` samples 200 points on `(0, pi/(2 B_0))` and hands the first sign change to
`brentq`, rather than calling `brentq` on the whole interval. Two reasons: the determinant is a
product of trigonometric and hyperbolic factors and can change sign more than once, and only the
**first** crossing is the fundamental; and the endpoint `pi/(2 B_0)` is a cutoff argued from the
physics — the critical slab is always thinner than the first zero of its own fundamental cosine
— not a place the determinant is known to have any particular sign, so it cannot be used as a
bracket on its own.

The interval is closed at both ends on purpose. At `a = 0` every `cos` is 1 and every `tanh` is
0, so `det H(0) = det(M_even V)`, finite and generically non-zero; at the cutoff the cosine
column survives through its `sin` term. Neither endpoint needs an epsilon.

## What it is not

It is not an independent derivation — both methods take their matrices from the same three
functions in `pn.py`, so an error in `marshak_matrix` would corrupt both identically and the
agreement in report §2's gap column would say nothing. That is what the analytic `P_1` check in
[06](06-verification.md) is for: it goes around `pn.py` entirely.
