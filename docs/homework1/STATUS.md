# Homework 1 — Status

`src/homework1` against `instruction_files/Assignment1.pdf`. Derivations are in
`homework1.tex`; the code is documented in `src/homework1/explanations/`.

| Question | Status |
|---|---|
| 1a–1d — flux components, both diffusion approximations, their errors | Done |
| 2 — numerical diffusion code with a delta source | Done |
| 3a–3c — critical `a/2` and `R_c` for `1.02 < c < 2` | Done |
| 3d–3e — asymptotic diffusion with the exact and modified-Marshak `l0(c)` | Not started |
| 4 — spherical `k = 1` critical radius | Done |
| 5 — bare critical mass of U-235 and Pu-239 | Done |

## Layout

| File | Holds |
|---|---|
| `exact_solution.py` | `nu0` on both branches, `phi_as`, `phi_tr` |
| `diffusion.py` | diffusion coefficients, Green's function, shooting solver |
| `criticality.py` | the five Question 3 methods |
| `spherical.py` | finite volume, the `k` iteration, dominance ratio and neutron balance |
| `materials.py` | the Sood benchmark table and the mass |
| `figures.py`, `tables.py` | shared plotting and logging helpers |
| `q1.py` … `q5.py` | one report plus its figures per question |

Run with `.\.venv\Scripts\python.exe -m homework1.main` from the repository root.

## Measured

- `nu0` fit vs. root: `5e-3 %` below `c = 1`, `0.10 %` at `c = 2`; the root matches
  Case's Table 8 Part II to `4.7e-06`.
- Q2 solver vs. closed form: `2e-8 %`; balance `Sigma_a int(phi) = 0.99998807`; `L2` error
  tracks `rtol` from `9e-06` down to `1.3e-12`.
- Q3 approximate inputs: error dominated by `z0`, under `0.1 %` in `a/2` to `c ~ 1.25`,
  `3.2 %` at `c = 2`. The printed `q = -0.0199` has the wrong sign — see
  the report, §3 "The Sign of the Quadratic Correction".
- Q4 vs. analytic `R_c = pi/B - z0`: `1.6e-4 %` to `2.4e-4 %` at `N = 400`, second order
  at exactly `4.00` per doubling.
- Q4 verification: sweep counts follow `c/(4c-3)` to six digits and are identical for both
  approximations; neutron balance closes to `1e-13`.
- Q5: Pu-239 `9.19 kg`, U-235 `32.70 kg` asymptotic; 30–40 % more classically.

The solver has one boundary condition at the outer surface: the extrapolated zero at
`R + z0`. A second, Robin-type treatment was implemented, measured and then removed —
see `explanations/07-removed-code.md`.

## Building the report

`.\docs\build.ps1 homework1`, never `latexmk` — see the header of `docs/build.ps1`. On this
machine a just-written file cannot immediately be overwritten (Windows EINVAL), so both the
build script and `figures.savefig` delete the target first. The log is affected as much as
the pdf — pdflatex stops on `I can't write on file homework1.log` before it reaches the pdf
— so the script now deletes both before every pass. Ruled out as causes: OneDrive,
Controlled Folder Access, antivirus, stale locks.
