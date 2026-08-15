# 06 — Verification

**Every check that was run, with the number it produced. None of these asserts that the code is
correct; each is something that would have failed loudly if it were not.**

## Question 1 — the reflected sphere

**Against an independent discretisation.** The closed-form root of report equation (13) was
checked against a two-region spherical finite-volume `k`-eigenvalue solve — 3000 cells,
harmonic-mean face diffusion coefficients across the material jump, the same extrapolated-zero
outer condition — for Pu-239 behind all three reflectors at `d = 1, 3, 10` mfp in both
continuous theories. Agreement is better than **0.4%** everywhere, and better than **0.1%**
except sodium at `d = 10`, whose 124 cm extrapolation layer is the least well resolved by a
uniform mesh.

That check is what makes the sodium anomaly of report §1 reportable: two independent
discretisations of the same equation agree, so the anomaly is in the model and not in the
solver. It also fixes the interface convention — a conservative finite-volume scheme conserves
`J`, and it reproduces the current-continuous column of [04](04-reflected-solver.md), not the
curvature-dropped one.

**`partial_current_factor`.** Against direct half-range integration of the discrete mode
`psi ∝ (nu0 - mu)^-1`: agrees to **eight digits** at `c = 0.9` and `c = 0.998`. Its `c > 1`
branch was checked to be the exact analytic continuation by substituting `nu0 = i/k0` into the
`c < 1` form — the result has zero imaginary part and matches the coded branch to **ten digits**
at `c = 1.3` and `c = 1.5`.

**`relaxation_rate`.** `k0(1.5) = 1.45110`, Case's Table 8 value to six digits.

## Questions 3–5 — the S_N solver

**The mesh is free.** Critical radius of the `c = 1.5` sphere at `S_10`:

| cells | 50 | 100 | 200 | 800 |
|---|---|---|---|---|
| `R_c` [mfp] | 1.685945 | 1.685954 | 1.685957 | 1.685957 |

Seven microns of mean free path across a factor of sixteen in mesh. **Every departure in the
report's Questions 3–5 tables is therefore angular truncation, not spatial** — which is the
claim those tables are built to make, so it is the load-bearing check of the three.

**Uncollided transport in the sphere.** A pure absorber with a uniform unit source has
`phi(0) = (1 - e^{-Sigma R})/Sigma`. At `Sigma = 1, R = 2` the exact value is `0.864665`; the
code gives `0.864733` (`S_2`), `0.864645` (`S_4`), `0.864602` (`S_32`) at 800 cells. This
exercises the areas, the volumes, the `mu = -1` starting direction and the `r = 0` reflection
with no scattering to hide behind — the parts of `sphere.py` the critical-size results would
not isolate.

**The exact slab benchmark.** At `c = 1.5` the one-speed critical half-thickness is `0.605055`
mfp. The S_N sequence descends onto it monotonically — `0.609042` (`S_10`), `0.606407`
(`S_16`), `0.605631` (`S_24`), `0.605195` (`S_48`) — roughly second order in `1/N`.

**The fixup fires when it should.** A pure absorber 20 mfp thick on **4 cells**, `S_10`, source
confined to the first cell: the fixup fires on 5 of 40 cell solves and the scalar flux stays
non-negative, `[0.968, 0.032, 0, 0]`. Without it the last two cells oscillate in sign.

**And never fires in Questions 3–5.** Counted over whole `k` calculations at `c = 1.5`, `S_10`,
critical size: 0 of 640 and 0 of 17000 cell solves in the slab (4 and 100 cells), 0 of 1364 and
0 of 38500 in the sphere. The fission source is proportional to the flux, so no cell is ever
starved, and a critical system is thin enough that `Sigma_t dx/|mu| < 2` even on 4 cells.
**The report's S_N tables can therefore be read as pure diamond-difference results**, with no
fixup-induced first-order error mixed in.

**The inner tolerance is the cheap knob.** Tightening `run_sn`'s tolerance from `1e-6`
to `1e-11` triples the sweep count (149 to 446) and does not move `k` in its first eight
digits: an inexact inner solve is absorbed by the next outer. `1e-8` sits comfortably inside
that plateau.
