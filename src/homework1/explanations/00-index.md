# 00 — Index

Explanations for the Assignment 1 code in `src/homework1/`, in question order. Results and
figures are in `docs/homework1/homework1.tex`; the status of each question is in
`docs/homework1/STATUS.md`.

| # | Title | What it settles |
|---|---|---|
| [01](01-case-eigenvalue.md) | The Discrete Eigenvalue | How `nu0` is found, and why it turns imaginary above `c = 1`. |
| [02](02-diffusion-coefficients.md) | The Two Diffusion Coefficients | Where `D = (1-c) nu0^2` comes from, and why it survives above `c = 1`. |
| [03](03-delta-source-and-boundary.md) | The Delta Source and the Outer Boundary | Why the source is a current condition and the truncation a radiation condition. |
| [04](04-solver-history.md) | Why Shooting, and What the Finite-Volume Solvers Measured | What the removed solvers cost, and two false trails not worth re-running. |
| [05](05-criticality-relations.md) | The Criticality Relations and the Sign of `q` | The one pair of relations behind all five methods, and why the printed `q` has the wrong sign. |
| [06](06-spherical-finite-volume.md) | Finite Volume in the Sphere | Why the balance form needs no special case at `r = 0`, and what the two boundary treatments cost. |
| [07](07-boundary-and-initial-conditions.md) | Boundary and Initial Conditions of the `k` Problem | Which two conditions the problem takes, why there is no initial one, and what `z0` is spent on. |
| [08](08-k-iteration.md) | The Bell & Glasstone `k` Iteration | Why `k` is a ratio of fission integrals, and why its convergence rate is `c/(4c-3)`. |
| [09](09-neutron-balance.md) | The Neutron Balance Check | What production = absorption + leakage does and does not verify. |
| [10](10-critical-masses.md) | Critical Masses, and the Two U-235 Rows | The measured masses, their 8–12 % boundary spread, and which U-235 data is used. |

Questions 1 and 2 are covered by 01–04, Question 3 by 05, Question 4 by 06–09, and
Question 5 by 10.
