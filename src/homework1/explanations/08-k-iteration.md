# 06 — The Bell & Glasstone `k` Iteration

**`k` is updated by the ratio of successive fission-source integrals, and the flux must be
renormalised each sweep or the convergence test loses its digits.**

The algorithm (Bell & Glasstone, pp. 189–192) as implemented in `k_eigenvalue`:

1. start from a flat flux and `k = 1`;
2. form the fission source `S = (1/k) nu Sigma_f phi`;
3. solve the removal problem `(-D grad^2 + Sigma_a) phi_new = S`, one banded solve;
4. set `k_new = k * integral(nu Sigma_f phi_new dV) / integral(nu Sigma_f phi dV)`;
5. repeat until `|k_new - k| < 1e-10 k_new`.

Step 4 is the trick. Dividing the source by `k` means a solution exists at *every* radius,
not only the critical one, so criticality becomes a root search on `R` — `brentq` on
`k(R) - 1`, which is monotonically increasing — instead of an eigensolve.

**Renormalisation.** The converged flux is defined only up to a constant, and each sweep
multiplies its amplitude by roughly `k`. Over the several hundred sweeps this problem needs,
the iterate would drift by many orders of magnitude and `|k_new - k|` would be a difference
of numbers whose significance has been lost. Rescaling by `1/max(phi)` keeps it `O(1)` and
changes nothing physically, because the ratio in step 4 is homogeneous in `phi`.

## The convergence rate depends on `c` alone

The error falls by the dominance ratio per sweep, which with the spherical modes
`B_n = n pi / r_outer` is `(Sigma_a + D B_1^2)/(Sigma_a + D B_2^2)` — `dominance_ratio` in
the code. At a *critical* radius it collapses: criticality means
`Sigma_a + D B_1^2 = nu Sigma_f`, and `B_2 = 2 B_1`, so

    k_1 / k_0 = nu Sigma_f / (4 nu Sigma_f - 3 Sigma_a) = c / (4c - 3)

`D` has cancelled. That is why the classical and asymptotic runs take an *identical* number
of sweeps at every `c` despite their different critical radii, and it is the sharpest check
that the iteration behaves as the theory says:

| `c` | `c/(4c-3)` | measured | sweeps (both approximations) |
|---|---|---|---|
| 1.02 | 0.944444 | 0.944444 | 280 |
| 1.05 | 0.875000 | 0.875000 | 133 |
| 1.10 | 0.785714 | 0.785714 | 79 |
| 1.20 | 0.666667 | 0.666667 | 50 |
| 1.50 | 0.500000 | 0.500000 | 32 |
| 2.00 | 0.400000 | 0.400000 | 25 |

The measured ratio agrees to six digits; the residual is that the radius used is the
*numerical* critical radius, which is not exactly critical.

**Why the sweep cap is 20000.** As `c -> 1` the ratio tends to 1 and convergence becomes
arbitrarily slow, and radii away from critical are worse still. The cap costs nothing (each
sweep is one `O(N)` solve), and the root search brackets within `+/- 5 %` of the analytic
radius so `k` is never evaluated far from critical. The whole Question 4 sweep runs in about
a second, so the acceleration schemes Bell & Glasstone discuss are not needed here.
