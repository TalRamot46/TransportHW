# 06 — The Analytic Comparison

**Question 5 adds no simulation. It reads the Question 3 and 4 caches and overlays two
closed-form curves, both of which are one call each into Assignment 1.**

## Two new modules

| file | owns |
|---|---|
| `analytic.py` | The two bare-sphere critical radii of report §5, in mean free paths. |
| `q5.py` | Loads both cached sweeps, prints the comparison table, writes the figure. |

`analysis.py` calls `run_analytic_comparison()` last, after the three studies, so
`python -m homework4.analysis` produces the Question 5 figure too. `q5.py` also runs standalone
(`python -m homework4.q5`), which is the fast path: it never touches the simulation, so it
returns in about a second whether or not the caches were just rebuilt.

## What is reused

`analytic.py` contains no formula of its own. Both theories are one call into Assignment 1:

    classical_radius -> homework1.criticality.critical_dimensions(c, 'marshak')
    milne_radius     -> homework1.criticality.critical_dimensions(c, 'transport-ref')

Those are parts (b) and (a) of Assignment 1 Question 3, already validated there against Case's
Table 8 and Table 23. `RADIUS = 3` picks `Sigma_t R_c` out of the four-tuple that function
returns. The figure and table helpers are `homework1.figures` and `homework1.tables`, so the
panel styling and the delete-before-write save are inherited rather than repeated.

## The variant that was tried and dropped

Both theories are `Sigma_t R_c = pi/B - z0` and differ in the pair `(B, z0)`, as report §5 sets
out. That invites a third, intermediate curve: the asymptotic diffusion coefficient
`D0 = (c-1)|nu0|^2` paired with the *diffusion* surface condition `z0 = 2 D0` rather than the
Milne one. It was implemented and then removed, for two reasons.

The first is that it is not what this assignment means by asymptotic diffusion, which is the
Milne-corrected solution and nothing else. The second is that it barely separates from it:
substituting `D0` into `B = sqrt((c-1)/D0)` gives `B = 1/|nu0|` exactly, so the intermediate
curve **shares a buckling with asymptotic Milne** and can differ only through `z0`, by a few
hundredths of a mean free path. The two were indistinguishable on the radius and mass panels and
crossed near `c = 1.5`, where `2 D0` and the Milne `z0(c)` agree to four decimals.

Adding it back would be a fifth entry in `homework1.criticality.METHODS`, which has no such
pairing, or a four-line function here. It is not worth either.

## Units, and the one conversion

`analytic.py` works entirely in mean free paths, matching Assignment 1. `q5.py` multiplies by
`config.mfp` once, in `build_theories` and again in `_comparison_table`, and converts to grams
with `critical_mass`. Nothing else in the module carries a unit, so a density change would flow
through correctly — though the caches are all at `rho = 30`, so this is untested.

## Grouping the c values

The cached `c` column is `P1 + 2*P2` in floating point, so the two rows that should both read
`1.05` differ in the sixteenth digit and `np.unique` would split them. `load_study` rounds the
column to six decimals; the tabulated probabilities are exact multiples of `0.01`, so this is a
repair rather than a tolerance.

## What the comparison shows

The table printed by `_comparison_table` is reproduced as Table 1 of the report. The classical
departure runs from `-0.38%` at `c = 1.05` to `-14.19%` at `c = 2.00` and is monotone in `c`
apart from the low-`c` scatter. Asymptotic Milne stays inside `2%` across the range but for two
low-`c` points.

This is not a verification of the Monte Carlo. The bisection stops at `dr < 0.05r`
([03](03-criticality-search.md)), so each radius carries about `±2.5%` of its own, which is the
same size as the asymptotic residual. It bounds the Monte Carlo from outside instead: agreement
this close over a factor of seven in radius rules out a systematic error much above the stopping
tolerance, which is the closest thing to a check this assignment has
([05](05-limitations.md) §5).
