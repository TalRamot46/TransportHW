# 06 — The Benchmark Data

**`materials.py` and the two tables in `criticality.py` are the only hard-coded data in the
assignment. Two of the four have a trap in them.**

## `BENCHMARK`, and the second U-235 row

`BENCHMARK` holds the five Sood, Forster & Parsons rows; `FISSILE = ('Pu-239', 'U-235')` names
the two with `c > 1`, and the other three have `c <= 1` so no bare critical sphere exists for
them at any radius.

`PROMPT_U235` is a **sixth** row that is deliberately not in `BENCHMARK`. It is the U-235 data
given in the task prompt rather than the assignment PDF, and it is not self-consistent: its own
cross sections give `c = 1.3646`, not the `c = 1.50` quoted with them (report §5). It is kept
as a module-level constant and reported alongside so the discrepancy is visible rather than
silently resolved.

`q5._mass_table` iterates `[BENCHMARK[n] for n in FISSILE] + [PROMPT_U235]` — which is why
adding a row to `BENCHMARK` will *not* put it in the Q5 table, and why `PROMPT_U235` appears
there despite not being in the dict.

**The guard.** `Material.sigma_t_sum` exists only so `q5._mass_table` can check
`Sigma_f + Sigma_c + Sigma_s == Sigma_t` on every row it prints and log a warning otherwise.
All six rows currently pass. That check is the reason a future transcription typo gets caught
instead of propagating into a critical mass.

## `Material.c` includes fission

`c = (nu Sigma_f + Sigma_s) / Sigma_t`, not `Sigma_s / Sigma_t`. Every secondary counts,
which is what makes `c > 1` possible at all. The reflector rows have `nu = Sigma_f = 0`, so
for them the two definitions coincide.

## Case's tables are arrays, not fits

`criticality.py` carries two transcribed tables:

| constant | what it holds | why it is stored this way |
|---|---|---|
| `CASE_TABLE_23_C` / `_CZ0` | the **product** `c z0(c)`, `c = 0 … 3.0` | the product barely moves — 0.7104 to 0.7199 across all of `c >= 0.8` — so interpolating it is far more accurate than interpolating `z0`, which halves over the same range |
| `CASE_TABLE_8_C` / `_K0` | the root `k0 = 1/|nu0|`, `c = 1.0 … 2.0` | reference the eigenvalue solver is checked against in `q3.report` |

`extrapolation_distance_table` interpolates the product and *then* divides by `c` — reversing
those two steps would lose most of the accuracy the table is stored for. It raises outside
`c in [0, 3]` rather than extrapolating.

Only Table 23 is on the `'transport-ref'` path (`METHODS['transport-ref'] = ('exact',
'table')`), so that method is usable for `1 < c <= 3`; Table 8 is used solely as a check in
`q3.report` and constrains nothing.
