# 02 — The S_N Solver

**Each geometry owns its own `cell_flux` and its own pair of sweeps, so the direction of a march
is visible in the function you are reading rather than in a flag it was passed.**

## Two cell solves, one per geometry

`slab.cell_flux(removal, source, mu_abs, psi_in)` is report eq. (22) forward and (23) backward —
the same expression, since the two differ only in which face is the inflow, and the caller has
already decided that by choosing a sweep. It has one closure, so its fixup is a single `if`.

`sphere.cell_flux(removal, source, radial_face, angular_face)` is report eq. (30), reached in the
report by substituting the two closures (27) and (28) into the balance (25). Its two faces are
not alike and are named for what they feed:

| face | weighted by | its outflow feeds |
|---|---|---|
| `radial_face` | the two shell areas `A_{i±1/2}` | the **next radial cell** of the same ordinate |
| `angular_face` | the two coefficients `alpha_{m±1/2}` | the **next ordinate** in the same cell |

A `Face(a_out, a_in, psi_in)` carries two weights because a shell has two areas and an angular
bin two alphas. `a_out` multiplies the flux on the face the march is travelling *towards*,
`a_in` the face it came *from*.

**Only the radial face swaps.** `sweep_forward` builds
`Face(mu*areas[i+1], mu*areas[i], psi_in)` and `sweep_backward` builds
`Face(mu*areas[i], mu*areas[i+1], psi_in)` — inner and outer exchanged. The angular face is
built once, by `_angular_face`, and is identical in both: `d(mu)/ds >= 0` means angle is
one-way, so its inflow is always at `m-1/2` and its outflow always at `m+1/2`.

`sphere._starting_direction` imports `slab.cell_flux` rather than reimplementing it. That is not
code sharing for its own sake: at `mu = -1` the factor `1 - mu^2` vanishes, the angular term
leaves the equation, and what remains *is* the plane balance at `|mu| = 1`.

## Backward, then forward

`sweep_all_angles` runs every `mu_m < 0` ordinate and then every `mu_m > 0` one, each half in
ascending `mu`. That ordering is load-bearing twice over and satisfying both at once is why it
is written this way:

- every backward ordinate must be swept before the forward ordinate `-mu_m` that the reflective
  boundary feeds from it — `sweep_backward` returns that inflow as its last value, and
  `sweep_all_angles` parks it in `reflected[mirror - m]`;
- in the sphere the `alpha` recursion must run upward from `mu = -1`, so `psi_low` threads
  through the backward loop and straight on into the forward one.

Since `core.ordinates` returns ascending `mu`, splitting at the sign is a *reordering of
nothing* — the negatives already come first. Reordering the ordinates produces a wrong answer,
not an error.

| report | code |
|---|---|
| the starting column `psi_{i,1/2}` at `mu = -1` | `SphereSolver._starting_direction` |
| eq. (23), the backward march | `sweep_backward` |
| eq. (22), the forward march | `sweep_forward` |
| eq. (30), one cell of it | `sphere.cell_flux` with its two faces |
| the `r = 0` reflection | `reflected[mirror - m]` in `sweep_all_angles` |

The one thing the code shows that the report does not is where the two outgoing fluxes *go*.
They come back positionally,

    psi[i], psi_in, psi_high[i] = self._cell(...)

and are consumed on different axes: `psi_in` is rebound in place and is the next cell's inflow
on the *same* pass of the `for i` loop, while `psi_high[i]` is stored into a column handed back
as `psi_low` on the *next* ordinate.

`_starting_direction` is the extra march the recursion needs, so the sphere costs `N + 1`
marches per S_N iteration where the slab costs `N`.

## The fixup terminates

The slab's fixup is one branch: with a single closure there is nothing left to flip once the
face is clamped. The sphere's is a loop, because clamping one face lowers `psi` and can drive
the other negative — the report's "Fixup" paragraph argues both why the test must be per face
and why clamping one lowers `psi`. It is bounded by two extra passes: each pass clamps at least
one more of the two faces, and with both clamped there is no closure left to produce a negative
value. **The trailing `max(value, 0.0)` is therefore unreachable** — a guard, not a code path.

## Two traps in the sweep

**The `alpha` factor of two.** `sphere.angular_coefficients` pairs the recursion
`alpha_{m+1/2} = alpha_{m-1/2} - w_m mu_m` with the cell coefficient `(A_out - A_in)/w_m`. The
pair `{-2 w_m mu_m, (A_out - A_in)/(2 w_m)}` is the *same* scheme written differently — but
mixing one half of each **moved the `c = 1.5` critical radius from 1.686 to 1.607 mfp**, a 5%
error that no mesh refinement removes. If a spherical result is a few percent off while the
mesh convergence looks clean, look here first.

**The centre of the sphere is not a special case.** `areas[0] = 0`, so the innermost face
carries no current and the cell balance does not determine the flux there; the diamond closure
does, and that value is what the `r = 0` reflection hands to the forward sweep. No branch
implements this — it falls out of `areas[0] = 0` in `sweep_backward`.

## Why `multiplying_medium` sets `Sigma_s = 0`

Report paragraph "The meaning of `c`" establishes that the critical size depends on `c` alone, so
`multiplying_medium(c)` may take the simplest split, `Sigma_t = 1, Sigma_s = 0, nu Sigma_f = c`.
[09](09-the-c-split.md) is the argument, and the measurement that the four splits agree.
The consequence lives in `run_sn_for_source`, which returns after a single iteration when
`sigma_s == 0`: with no scattering, one sweep inverts the transport operator exactly. Measured,
Questions 3 and 4 take **exactly one** sweep per inner, while Question 5 — real cross sections,
`Sigma_s/Sigma_t = 0.69` — takes up to **21**. See [03](03-k-iteration.md).

## Mesh sizes

`slab.N_CELLS = 200` and `sphere.N_CELLS = 100` are both far inside the converged plateau
measured in [06](06-verification.md). The sphere gets the smaller of the two only because it is
the expensive geometry.
