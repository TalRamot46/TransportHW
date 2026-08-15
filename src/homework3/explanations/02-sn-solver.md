# 02 — The S_N Solver

**`cell_flux` is the entire spatial discretisation in one function; everything
geometry-specific is pushed out into the list of links its caller passes it.**

## One cell solve for both geometries

`sn.cell_flux(removal, source, links)` solves

    sum over links of (c_out psi_out - c_in psi_in)  +  removal * psi  =  source

closing every link with the same diamond relation `psi_out = 2 psi - psi_in`. A link is a
`(c_out, c_in, psi_in)` triple describing one outgoing face. That indirection is why the sphere
needed no new cell solve:

| caller | links passed |
|---|---|
| `SlabSolver._direction` | one — `(abs(mu), abs(mu), psi_in)` |
| `SphereSolver._starting_direction` | one — `(1, 1, psi_in)`, the `mu = -1` slab sweep |
| `SphereSolver._direction` | two — the spatial face, plus the angular face `(alpha_{m+1/2}, alpha_{m-1/2}, psi_low)` |

Report equations (3) and (5) are the two cell balances being solved; the code never sees them
separately, only as link lists of length one and two.

## The fixup loop provably terminates

The negative-flux fixup clamps an offending outgoing flux to zero and re-solves the balance
with it held there. It has to be a loop, because clamping one face can drive another negative.
It is bounded by `len(links) + 1` passes: each pass clamps at least one more face, and once
every face is clamped there is no diamond closure left to produce a negative value, so the pass
after that returns. **The trailing `max(value, 0.0)` after the loop is therefore unreachable**
— it is there as a guard, not as a code path anything relies on.

## Three traps in the sweep

**Ordinate order is load-bearing.** `sn.ordinates` returns ascending `mu`, and both `sweep`
methods depend on it twice: every inward ordinate (`mu < 0`) must be swept before the outward
ordinate it feeds through the reflective boundary, and in the sphere the `alpha` recursion must
run upward from `mu = -1`. Reordering the ordinates produces a wrong answer, not an error.

**The `alpha` factor of two.** `sphere.angular_coefficients` pairs the recursion
`alpha_{m+1/2} = alpha_{m-1/2} - w_m mu_m` with the cell coefficient `(A_out - A_in)/w_m`. The
pair `{-2 w_m mu_m, (A_out - A_in)/(2 w_m)}` is the *same* scheme written differently — but
mixing one half of each **moved the `c = 1.5` critical radius from 1.686 to 1.607 mfp**, a 5%
error that no mesh refinement removes. If a spherical result is a few percent off while the
mesh convergence looks clean, look here first.

**The centre of the sphere is not a special case.** `areas[0] = 0`, so the innermost face
carries no current and the cell balance does not determine the flux there; the diamond closure
does, and that value is what the `r = 0` reflection hands to the outward sweep. No branch
implements this — it falls out of `a_in = areas[i] = 0` in `_direction`.

## Why `multiplying_medium` sets `Sigma_s = 0`

Report §"The meaning of `c`" establishes that the critical size depends on `c` alone, so
`multiplying_medium(c)` may take the simplest split, `Sigma_t = 1, Sigma_s = 0, nu Sigma_f = c`.
The consequence lives in `inner_iteration`, which returns after a single sweep when
`sigma_s == 0`: with no scattering, one sweep inverts the transport operator exactly. Questions
3 and 4 therefore cost one sweep per outer, while Question 5 — real cross sections,
`Sigma_s/Sigma_t = 0.69` — costs about ten. See [03](03-k-iteration.md).

## Mesh sizes

`slab.N_CELLS = 200` and `sphere.N_CELLS = 100` are both far inside the converged plateau
measured in [06](06-verification.md). The sphere gets the smaller of the two only because it is
the expensive geometry.
