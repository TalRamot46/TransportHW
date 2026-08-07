# Plan — Homework 1, Question 2 (numerical diffusion with a delta source)

Branch: `homework1-q2-diffusion`

## 1. Assessment of the current implementation

`solve_diffusion_numerical` (`src/homework1/diffusion.py:52`) is **correct**. Measured
against the analytic Green's function it gives a median relative error of `4.5e-5 %` for
`c = 0.5, 0.7, 0.9`, and it is exact at `x = 0`. Two things are worth making explicit.

### 1.1 The symmetry assumption is already in the code

The equation being solved is

```
-D phi''(x) + Sigma_a phi(x) = delta(x),     Sigma_a = 1 - c,   D = 1/3
```

Integrating across the source over `(-eps, +eps)` gives the jump condition

```
phi'(0-) - phi'(0+) = 1/D
```

The Green's function is even, `phi(x) = phi(-x)`, so `phi'(0+) = -phi'(0-)` and therefore

```
phi'(0-) = +1/(2D)          equivalently   J(0+) = -D phi'(0+) = 1/2
```

which is exactly the `target = 1.0 / (2.0 * D)` the code shoots at, on the half-domain
`[-a, 0]`. So the answer to "can I assume symmetry to get the derivative at zero" is that
**this is already what the code does** — the half-domain plus the current boundary
condition *is* the symmetry reduction. The physical reading of `J(0+) = 1/2` is that half
the source neutrons stream right and half stream left. This should be derived explicitly in
the report, since it is the answer to the assignment's "how should a delta-function source
be modeled in a numerical scheme?".

### 1.2 The real defect: the far boundary, not the source

The left boundary uses `phi(-a) = 0` with `a = 10/kappa`. The true Green's function at that
point is small but **non-zero**, so the relative error there is exactly `100 %`, and this
dominates `figs/diffusion_numerical_errors.pdf` — the plot currently reads as if the solver
were failing, when in fact the interior is accurate to 5 decimal places.

Two further consequences:

- This is a **modeling** error (domain truncation), not a discretization error. It does not
  shrink under mesh refinement, only under increasing `a`. A convergence study run on the
  current code would therefore show an error floor and look broken.
- The amplitude bias from the Dirichlet condition is only `-2.1e-9` relative, so the
  interior is fine; the artifact is purely local to `x = -a`.

Ruled out as error sources (measured, not assumed): the `brentq` tolerance — it recovers
the shooting slope `s = 1.362e-4` to all printed digits — and the `solve_ivp` tolerances,
which contribute the `~4.5e-5 %` interior floor.

Minor: the comment on line 75 says "3 decay lengths" but the code uses `a = 10/kappa`, i.e.
10 diffusion lengths.

## 2. Preferable approach

### 2.1 Keep the symmetry reduction, fix the far boundary

Replace `phi(-a) = 0` with the **radiation (Robin) condition**

```
phi'(-a) = kappa * phi(-a),        kappa = sqrt(Sigma_a / D)
```

For `x < 0` the decaying infinite-medium solution is `phi = A e^(kappa x)`, which satisfies
this identically. Imposing it annihilates the reflected mode `e^(-kappa x)`, so the
truncated domain reproduces the infinite-medium Green's function **exactly, for any `a`**.

This is legitimate rather than circular: `kappa` is a property of the differential operator
(`kappa^2 = Sigma_a / D`), not of the solution being validated against. It is the standard
non-reflecting boundary for this problem. With it, the remaining error is purely numerical,
so a convergence study becomes meaningful.

### 2.2 Drop the root-finder — the problem is linear

The ODE is linear and homogeneous away from the source, so the shooting function is
*exactly* proportional to `s`: `y(x; s) = s * y(x; 1)`. `brentq` is currently performing
~40 integrations to find the root of a straight line. One integration suffices:

```
integrate once from -a to 0 with (phi, phi') = (1, kappa)      # Robin BC, arbitrary scale
scale = target / phi'(0)                                        # match the source condition
phi   = scale * phi
```

Same answer, no bracketing heuristic, no tolerance to tune.

### 2.3 Add a finite-volume solver as the primary method

Shooting is a fine check but a poor foundation. A cell-centred finite-volume discretisation
solved as a tridiagonal system is:

- the natural place to answer the delta-source question, because it permits a direct
  comparison of the two modelling choices (below);
- second-order accurate, so it produces a clean `O(dx^2)` convergence plot;
- **the code Q4 needs anyway** (spherical geometry, `k`-eigenvalue power iteration). Writing
  it here means Q4 is a geometry change plus an outer iteration, not a rewrite.

### 2.4 Model the delta source two ways and compare

This is the assignment's actual question, so answer it empirically:

- **(a) Half-domain, current boundary condition.** `x in [0, a]`, with `-D phi'(0) = 1/2`
  from the symmetry argument. The kink at the origin sits exactly on a boundary, so it is
  represented exactly.
- **(b) Full domain, source smeared over one cell.** `x in [-a, a]` with `S_i = 1/dx` in the
  cell containing the origin, zero elsewhere. This is the approach that generalises to
  sources that are not delta functions, but it smears the kink over one cell.

Expected result, to be confirmed: both converge at `O(dx^2)`, with (b) carrying a larger
constant and a visible local error at the origin. That comparison is the report's answer.

### 2.5 Cover both approximations at no extra cost

The asymptotic diffusion solution is the *same* equation with a different diffusion
coefficient. Matching amplitude and decay rate to `phi = 1/(2(1-c) nu0) exp(-|x|/nu0)`:

```
Sigma_a = 1 - c        (unchanged)
D_asy   = (1 - c) nu0^2        =>   kappa = 1/nu0
```

(Check: `kappa^2 = Sigma_a / D_asy = 1/nu0^2`.) So one solver covers both approximations by
parameter choice, and the assignment's "either ... or" can be answered with "both".

## 3. Work items

1. `solve_diffusion_shooting` — half-domain, Robin BC, single scaled integration (2.1–2.2).
2. `solve_diffusion_fv` — cell-centred tridiagonal solve, half-domain, current BC (2.3).
3. `solve_diffusion_fv_full` — full domain, one-cell smeared source (2.4b).
4. `diffusion_coefficients(c, approximation)` — returns `(D, Sigma_a)` for `"classical"` /
   `"asymptotic"` (2.5); refactor `phi_classical_diffusion` / `phi_asymptotic_diffusion`
   onto the shared closed form.
5. Convergence study over `dx`, reporting observed order for methods 2 and 3.
6. Plots: solution comparison, error vs `x` (log), and a convergence plot with an `O(dx^2)`
   reference slope.
7. Fix the stale "3 decay lengths" comment.
8. Keep `solve_diffusion_numerical` working so the existing Q2 figure path does not break.

## 4. Verification

- Every solver reproduces the closed form to the tolerance claimed, for
  `c = 0.5, 0.7, 0.9` and both approximations.
- Observed convergence order is `2.0 +/- 0.1` for the finite-volume methods.
- Neutron balance: `Sigma_a * integral(phi dx) = 1` (the unit source) to solver tolerance.
- The Robin BC removes the `100 %` endpoint artifact from the error plot.

## 5. Results (implemented, measured)

All items in section 3 are implemented. Measured on `c = 0.5, 0.7, 0.9` for both
approximations:

| Solver | Max relative error | Neutron balance |
|---|---|---|
| Shooting, radiation BC | `2.0e-8 %` | `0.99998807` |
| Finite volume, half-domain, 500 cells | `1.49e-2 %` | `0.99995459` |
| Finite volume, smeared source, 501 cells | `3.95e-2 %` | `0.99995457` |

- **Convergence order is `2.000` for both finite-volume treatments**, across every `c` and
  both approximations. The smeared source has a `~14x` larger error constant at equal `dx`,
  exactly as anticipated in 2.4: `6.49e-2` vs `4.72e-3` relative L2 error at the coarsest
  mesh. Same order, worse constant — that comparison is the report's answer to the
  assignment's delta-source question.
- **The endpoint artifact is gone.** The error is now flat across the domain at
  `<= 1.5e-2 %`, instead of reaching `100 %` at the outer boundary.
- The residual `4.5e-5` balance deficit is not solver error: it is the physical tail of the
  Green's function beyond 10 diffusion lengths, which the truncated domain does not carry.
- `D_asy = (1-c) nu0^2` reproduces `phi_asymptotic_diffusion` to `1e-13 %`, confirming the
  identification in 2.5.
- Every case shares identical error and balance figures because the domain is scaled as
  `a = 10/kappa`; in units of `kappa x` all cases are the *same* boundary-value problem, so
  `c` and the approximation enter only through the scaling. Worth stating in the report.

### Note on the trapezoid trap

The first balance check reported `0.99005` for the half-domain solver — a 1% deficit that
looks like a broken solver. It is not: applying the trapezoid rule to *cell-averaged* data
omits the half-cell slivers at each end, discarding a fraction `kappa dx / 2` of the
integral, which is exactly 1% at the default resolution. `absorption_balance` therefore
takes a `quadrature` argument, and the finite-volume callers pass `'midpoint'`.

### Discarded hypothesis

The `brentq` tolerance in the original solver was initially suspected of limiting accuracy,
on the grounds that an absolute `xtol=2e-12` would be coarse relative to a small root. It
was measured and is not a factor: the root is `s = 1.362e-4`, and `brentq` recovers it to
all printed digits. The original solver's accuracy was limited only by the outer boundary
condition.
