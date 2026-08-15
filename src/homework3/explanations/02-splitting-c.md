# 02 — Splitting `c` into Scattering and Fission

**Questions 3 and 4 give only `c`, not the split between `Sigma_s` and `nu Sigma_f`; the
critical size does not depend on the split, so the code counts every secondary as a fission
neutron.**

The `k`-eigenvalue problem divides only the fission part of the source:

    mu dpsi/dx + Sigma_t psi = (1/2) [ Sigma_s phi + (1/k) nu Sigma_f phi ]

At `k = 1` the bracket collapses to `(Sigma_s + nu Sigma_f) phi = c Sigma_t phi`, which depends
on `c` alone. Any split with `Sigma_s + nu Sigma_f = c Sigma_t` therefore has **the same
critical size**, and `multiplying_medium(c)` may take the simplest one,

    Sigma_t = 1,    Sigma_s = 0,    nu Sigma_f = c

What the split does change is the *path* to that size. With `Sigma_s = 0` the scattering
iteration is empty — one sweep inverts the transport operator exactly, which is why
`inner_iteration` returns immediately in that case — and `k` away from criticality means
something different: here `k_inf = nu Sigma_f / Sigma_a = c`, so the `k(size)` curves of
Questions 3 and 4 saturate at `c`, whereas a split with scattering would saturate elsewhere.
Only the crossing `k = 1` is common to all splits.

Question 5 has no such freedom: the benchmark rows carry their own `Sigma_s`, and there the
inner iteration does real work — about five sweeps per outer for Pu-239, against one for the
`c`-only problems.

## Why `c > 1` is not a special case here

Nothing in the sweep or in the iteration is restricted to `c < 1`. The discrete eigenvalue
`nu0` of Assignment 1 turns imaginary above `c = 1` and its formulae need care, but the S_N
equations are algebraic in `c` and multiplying media are simply the case where a critical size
exists at all. The only place the distinction survives is in the *reference* values, which come
from `homework1.criticality.critical_dimensions` and use `compute_nu0_magnitude` — Assignment
1's `c > 1` branch.
