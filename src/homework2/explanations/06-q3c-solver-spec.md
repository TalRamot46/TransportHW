# 06 — The Q3(c) Solver

**`solver.py` is 70 lines because the absorption is factored out before anything is
discretised: what is stepped is the bare heat equation, and `c` never enters the loop.**

Report §3, *Part 3(c)* derives the scheme, its stability limit and its conservation identity.
This file is the map from those formulas to the module; the numbers they produce are in
[05](05-verification.md).

## The shape of the module

| function | owns |
|---|---|
| `_grid` | the half domain `[0, L]`, `L` set by the tail rather than guessed |
| `_time_step` | the largest `dt` satisfying `r <= 1/2`, times `SAFETY` |
| `_pulse` | the delta as `V/h` in node 0 — report eq. (35) |
| `_warm` | the analytic Gaussian at `WARM_T0`, a diagnostic start |
| `step` | one FTCS sweep, public so the stability check can drive it directly |
| `_march` | `n_steps` sweeps at fixed `r` |
| `solve` | `D` from `diffusion.py`, then one march per requested time |
| `mass` | `2 * trapezoid(u, x)`, the half domain doubled |

`solve` takes `D` from `diffusion.diffusion_coefficient`, so it inherits the `c = 1` limit and
the `c > 1` sign handling for free ([03](03-diffusion-module.md)). Everything `c`-dependent
happens outside the loop: `D` picks the step size, and the prefactor
`exp(-(1-c) V t)` is reapplied on the way out. That is why one 70-line module answers both
approximations, every `c`, and both sides of `c = 1`.

## Two decisions that are not obvious from the code

**The step is shrunk to land on each requested time.** `solve` takes
`ceil((target - t) / dt)` steps of size `(target - t) / n_steps` rather than stepping by `dt`
and stopping past the target. Shrinking only lowers `r`, so stability is untouched, and it
removes a time error of up to one `dt` that would otherwise be indistinguishable from a
discretisation error in the comparison against the closed form.

**`step` is public, `_march` is not.** The stability check has to run the scheme *past* its own
limit, which `solve` will not do — `_time_step` caps `r` and the target snapping lowers it
further. Driving `step` directly on a 101-node spike is the honest way to show the blow-up.

## Two traps

**The symmetry row is not the interior row.** `new[0] += 2*r*(u[1] - u[0])` — the factor of two
comes from the ghost value `u_{-1} = u_1`. Dropping it leaves the scheme stable and
plausible-looking while quietly leaking particles at the origin. Nothing in a plot would show
it; the conservation check would fail immediately.

**`r` is defined with `D V`, not `D`.** Under this assignment's units `V = 1` and the two
coincide, so an error here costs nothing today and everything the moment the units change.
`V` is kept explicit in `r`, in `_pulse` and in the prefactor for that reason.

## Cost

`dt` falls as `h^2` and the array grows as `h^-1`, so **the total work goes as `h^-3`**.
At the default `n_nodes = 2000` a march to `t = 15` is ~2×10^4 sweeps of a 2001-element array.
Halving the node count costs one bit of spatial accuracy and buys a factor of eight — which is
the trade to make first if this ever becomes slow. Crank–Nicolson would remove the `r <= 1/2`
limit altogether and let `dt ~ h`, at the price of a banded solve per step; it is the right
answer if the mesh ever needs to be fine, and the wrong complication at this one.
