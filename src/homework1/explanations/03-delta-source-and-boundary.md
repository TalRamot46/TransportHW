# 03 — The Delta Source and the Outer Boundary

**The delta source is never discretised: symmetry turns it into the current condition
`J(0+) = 1/2`, and the truncated outer boundary carries a radiation condition.**

Integrating `-D phi'' + Sigma_a phi = delta(x)` across the origin gives the jump
`phi'(0-) - phi'(0+) = 1/D`. The infinite-medium Green's function is even, so
`phi'(0+) = -1/(2D)`, that is

    J(0+) = -D phi'(0+) = 1/2

— half the source neutrons stream each way. The problem is then solved on `x >= 0` and the
singularity never enters the numerics. This is the assignment's question of how to model a
delta source: you don't, you convert it.

## The outer boundary

`x = a` is a truncation of an infinite medium, not a physical surface. Zero flux there
forces `phi` to vanish where the true solution is small but non-zero, giving a `100 %`
relative error at `x = a` that no tolerance reduces. The radiation condition

    phi'(a) = -kappa phi(a)

is satisfied identically by `exp(-kappa x)`, kills the growing mode, and makes the
truncation exact for any `a`. The measured error is then flat across the domain, and
`Sigma_a integral(phi dx) = 0.99998807` — the residual is the physical tail beyond
`a = 10/kappa`, not solver error.

## Why one integration suffices

Away from the source the equation is linear and homogeneous, so its solution is exactly
proportional to its starting amplitude. `solve_diffusion_shooting` integrates once from `a`
to `0` (the direction in which the solution grows, which is the stable one) starting from
`[1, -kappa]`, then rescales so that `-D phi'(0) = 1/2`. No root-finding.

Accuracy is set by the integrator tolerance, not by a mesh, so the convergence study
tightens `rtol`: the relative `L2` error falls from `9e-06` at `rtol = 1e-4` to `1.3e-12`
at `rtol = 1e-11`, one for one with the request.
