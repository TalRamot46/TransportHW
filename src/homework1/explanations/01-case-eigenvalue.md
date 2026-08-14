# 01 — The Discrete Eigenvalue

**`nu0` is the root of `c nu0 arctanh(1/nu0) = 1`; above `c = 1` it turns imaginary and
the same equation becomes `c arctan(k0) = k0`.**

For `c < 1`, with `k0 = 1/nu0 in (0,1)`, the equation is `arctanh(k0)/k0 = 1/c`, solved by
`brentq` in `compute_nu0_numerical`. The course fit

    nu0 = 1 / sqrt(1 - c^p(c)),   p(c) = 2.47412 + 0.00363081/c^2 - 0.0352458 c + 0.557498/c

is `compute_nu0_approx`, within `5e-3 %` of the root over `c = 0.5 … 0.95`.

## Above `c = 1`

No real root survives: the eigenvalue moves onto the imaginary axis, `nu0 = i|nu0|`.
Using `arctanh(ik) = i arctan(k)` the equation is real again,

    c arctan(k0) = k0,   k0 = 1/|nu0|

the form tabulated by Case, de Hoffmann & Placzek (Table 8, Part II). The root reproduces
every entry of that table to `4.7e-06`.

The fit continues onto this branch by itself: its radicand `1 - c^p(c)` turns negative,
which *is* the statement that `nu0` is imaginary, so `|nu0| = 1/sqrt(c^p(c) - 1)`. Measured
against the root, that is `1e-4 %` at `c = 1.02` and `0.10 %` at `c = 2`.

An imaginary `nu0` turns the infinite-medium mode `exp(-x/nu0)` into `cos(x/|nu0|)` — the
flux shape of a finite critical system, which is what [05](05-criticality-relations.md)
builds on.
