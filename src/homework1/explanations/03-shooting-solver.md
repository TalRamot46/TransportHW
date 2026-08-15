# 03 — The Shooting Solver

**`solve_diffusion_shooting` has no mesh, no root-find and no source term. All three are
consequences of the equation being linear and homogeneous away from the origin.**

Report §2 derives the delta-to-current conversion and the radiation condition; this is how the
function is built on them.

## One integration, then a rescale

The solver never searches for a starting amplitude. It integrates once from `a` to `0` from
`[phi, phi'] = [1, -kappa]` — the radiation condition with unit amplitude — and then divides
by whatever `phi'(0)` came out:

    phi = sol.y[0][::-1] * (-SOURCE_CURRENT / D_eff) / sol.y[1][-1]

Because the ODE is linear and homogeneous, that single rescale imposes `J(0+) = 1/2` exactly.
Integrating **inwards** is the other half of it: the decaying solution grows in that
direction, which is the numerically stable one to follow.

Two consequences for anyone reading the code:

- `SOURCE_CURRENT = 0.5` is the entire source term. There is no delta anywhere in the numerics.
- The accuracy knob is `rtol`, not a cell count. `convergence_study` therefore sweeps tolerance
  where a mesh-based solver would sweep resolution, and the error tracks `rtol` one for one
  (see [08](08-verification.md)).

## The `x_eval` path

`x_eval` lets `q1` evaluate the solver on its own plotting grid, so the exact-transport and
diffusion curves are compared point for point with no interpolation. Two details make it work:

- **`x = 0` is always integrated**, even when the caller's grid omits it, because the rescale
  reads `phi'(0)`. The extra point is prepended and then stripped — that is what the
  `at_origin` flag is doing.
- **The grid is validated against `a`**, since `a = 10/kappa` depends on `c`, and a grid valid
  for one `c` may exceed the domain for another.

## All six Q2 cases are one problem

The domain is `a = 10/kappa`, so measured in units of `kappa x` every `c` and both
approximations are the **same** boundary-value problem. `c` and the approximation enter only
through the scaling. That is why the Q2 error and balance figures look identical across all six
panels — not a plotting bug, and not a check that has been passed six times.

## `absorption_balance` needs no quadrature argument

It applies `np.trapezoid` to point values on the returned grid and doubles the half-line
integral. That is only correct because the shooting solver returns **point** values; a
cell-averaged solver would need the end slivers, which is exactly what once produced a
spurious 1% deficit — see [07](07-removed-code.md).
