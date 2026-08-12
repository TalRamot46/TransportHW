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

**Why the sweep cap is 20000.** Convergence is at the dominance ratio
`k_1/k_0 = (Sigma_a + D B_0^2)/(Sigma_a + D B_1^2)`, which tends to 1 for a large weakly
multiplying sphere — about `0.94` at `c = 1.02`. The cap costs nothing (each sweep is one
`O(N)` solve), and the root search brackets within `+/- 5 %` of the analytic radius so that
`k` is never evaluated far from critical. The whole Question 4 sweep runs in about 0.2 s, so
the acceleration schemes Bell & Glasstone discuss are not needed here.
