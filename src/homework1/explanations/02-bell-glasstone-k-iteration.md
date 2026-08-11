# 02 — The Bell & Glasstone `k` Iteration

**`k` is updated by the ratio of successive fission-source integrals, not by
solving an eigenvalue problem; the flux must be renormalised each sweep or the
iteration loses precision before it converges.**

The algorithm of Bell & Glasstone, pp. 189–192, is power iteration on the
fission source. Written in the form implemented in `k_eigenvalue`:

1. Start from a flat flux `phi^(0)` and `k^(0) = 1`.
2. Form the fission source `S^(n) = (1/k^(n)) nu Sigma_f phi^(n)`.
3. Solve the *removal* problem `(-D grad^2 + Sigma_a) phi^(n+1) = S^(n)`.
   This is a plain linear solve — the operator is fixed, so it is one banded
   back-substitution per sweep.
4. Update

       k^(n+1) = k^(n) * [ integral nu Sigma_f phi^(n+1) dV ]
                       / [ integral nu Sigma_f phi^(n)   dV ]

5. Repeat until `k` stops moving.

Step 4 is the whole trick. Dividing the source by `k` means a solution exists at
*every* radius, not only at the critical one, so the code never has to search for
an eigenvalue directly. The ratio of successive fission integrals measures how
much the medium multiplied over one generation; `k` is the fixed point at which
that ratio is one. This is why criticality can then be found by a root search on
the radius rather than by an eigensolver — see [03](03-boundary-and-analytic.md).

## Why the flux is renormalised

The converged flux is defined only up to a constant, and each sweep multiplies
its amplitude by roughly `k`. Left alone, `phi` drifts geometrically — over the
several hundred sweeps this problem needs, a `k` even a few per cent from unity
moves the amplitude by many orders of magnitude, and the convergence test on
`|k^(n+1) - k^(n)|` starts comparing numbers whose difference has been lost to
rounding. Rescaling by `1/max(phi)` after each sweep keeps the iterate at
`O(1)`. It changes nothing physically, because both integrals in step 4 are
linear in the flux and the scale cancels in the ratio.

## Why the sweep cap is high

Source iteration converges at the dominance ratio, the ratio of the first
harmonic to the fundamental:

    k_1 / k_0 = (Sigma_a + D B_0^2) / (Sigma_a + D B_1^2)

For a large, weakly multiplying sphere both bucklings are small, so this
approaches 1 and convergence is slow. At `c = 1.02` the critical sphere is
about 12 mean free paths across and the ratio is near `0.94`, needing several
hundred sweeps for a tolerance of `1e-10`; radii well away from critical are
worse still. The cap is therefore set at 20000 sweeps, which costs nothing —
each sweep is one `O(N)` solve — but two things are done to avoid paying it:

- the root search brackets tightly around the analytic radius (`+/- 5 %`), so
  `k` is never evaluated at a radius far from critical;
- the whole `c` sweep of Question 4 runs in about 0.2 s.

Accelerating the iteration (Chebyshev or Aitken extrapolation, both discussed by
Bell & Glasstone) would be the next step if this ever became a bottleneck. It
has not, so the plain iteration is kept because it is the one the assignment
asks for.
