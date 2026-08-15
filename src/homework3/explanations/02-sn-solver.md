# 02 — The S_N Solver

**`cell_flux` is the entire spatial discretisation in one function; everything
geometry-specific is pushed out into the list of faces its caller passes it.**

## One cell solve for both geometries

`sn.cell_flux(removal, source, faces)` solves the cell balance

    sum over faces of (a_out psi_out - a_in psi_in)  +  removal * psi  =  source

closing every face with the same diamond relation `psi_out = 2 psi - psi_in`, which collapses
to the single expression the loop evaluates:

    psi = [ source + sum_f (a_out + a_in) psi_in ] / [ removal + 2 sum_f a_out ]

Those two are report equations (7) and (8), derived there once for both geometries.

An `sn.Face(a_out, a_in, psi_in)` is **one outgoing face** of the cell: its two balance
coefficients, and the flux arriving through it. The coefficients differ only where the face is
weighted differently on the way in and on the way out, which is the whole reason there are two
of them:

| caller | faces passed | why the two coefficients differ |
|---|---|---|
| `SlabSolver._sweep` | one — `Face(abs(mu), abs(mu), psi_in)` | they don't; a slab face has no area factor |
| `SphereSolver._starting_direction` | one — `Face(1, 1, psi_in)`, the `mu = -1` slab sweep | they don't; `|mu| = 1` and the areas cancel |
| `SphereSolver._sweep` | two — spatial, then angular | the two areas of a shell; the two alphas of a bin |

Substituting the slab's row into the collapsed form returns report eq. (3) verbatim, which is
the quickest way to convince yourself the abstraction is not hiding anything.

**Reading the two coefficients off a balance.** Given a discrete balance, `a_out` is whatever
multiplies the flux on the face the sweep is travelling *towards*, and `a_in` whatever
multiplies the face it came *from*. In `_sweep` that is the line

    a_out, a_in = (inner, outer) if inward else (outer, inner)

— an inward ordinate leaves through the inner face, so the two swap. Everything else about the
sweep direction is already handled by the loop order.

## The fixup loop provably terminates

The negative-flux fixup clamps an offending outgoing flux to zero and re-solves the balance
with it held there. It has to be a loop, because clamping one face can drive another negative.
It is bounded by `len(faces) + 1` passes: each pass clamps at least one more face, and once
every face is clamped there is no diamond closure left to produce a negative value, so the pass
after that returns. **The trailing `max(value, 0.0)` after the loop is therefore unreachable**
— it is there as a guard, not as a code path anything relies on.

## Three traps in the sweep

**Ordinate order is load-bearing.** `sn.ordinates` returns ascending `mu`, and both `sn_iteration`
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
implements this — it falls out of `a_in = areas[i] = 0` in `_sweep`.

## Why `multiplying_medium` sets `Sigma_s = 0`

Report §"The meaning of `c`" establishes that the critical size depends on `c` alone, so
`multiplying_medium(c)` may take the simplest split, `Sigma_t = 1, Sigma_s = 0, nu Sigma_f = c`.
The consequence lives in `run_sn`, which returns after a single iteration when
`sigma_s == 0`: with no scattering, one sweep inverts the transport operator exactly. Questions
3 and 4 therefore cost one sweep per outer, while Question 5 — real cross sections,
`Sigma_s/Sigma_t = 0.69` — costs about ten. See [03](03-k-iteration.md).

## Mesh sizes

`slab.N_CELLS = 200` and `sphere.N_CELLS = 100` are both far inside the converged plateau
measured in [06](06-verification.md). The sphere gets the smaller of the two only because it is
the expensive geometry.
