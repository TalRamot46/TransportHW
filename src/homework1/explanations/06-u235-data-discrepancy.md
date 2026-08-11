# 06 — The Two U-235 Cross-Section Rows

**The code uses the U-235 row from the assignment PDF; the row supplied in the
task prompt is not self-consistent, and both are reported so the difference is
visible rather than silently resolved.**

Two different U-235 rows were in play. Pu-239 is identical in both and is not
affected.

| source | `nu` | `Sigma_f` | `Sigma_c` | `Sigma_s` | `Sigma_t` | `c` quoted |
|---|---|---|---|---|---|---|
| Assignment PDF (Sood et al.) | 2.70 | 0.065280 | 0.013056 | 0.248064 | 0.32640 | 1.30 |
| Task prompt | 2.70 | 0.065280 | 0.015672 | 0.180448 | 0.26140 | 1.50 |

Both rows are consistent in the sense that `Sigma_f + Sigma_c + Sigma_s` equals
the tabulated `Sigma_t`. They differ in the scattering ratio actually implied:

    c = (nu Sigma_f + Sigma_s) / Sigma_t

- PDF row: `(2.70 x 0.065280 + 0.248064) / 0.32640 = 1.3000` — matches its
  quoted `c` exactly.
- Prompt row: `(2.70 x 0.065280 + 0.180448) / 0.26140 = 1.3646` — does **not**
  match the `c = 1.50` quoted alongside it.

The prompt row therefore cannot be reconciled with itself, while the PDF row
reproduces both its own `c` and the published benchmark. `BENCHMARK['U-235']`
holds the PDF values and is what the results use.

The prompt row is kept as `PROMPT_U235` and reported on its own line, so the
consequence of the discrepancy is measured rather than argued:

| row | `c` | classical `M_c` | asymptotic `M_c` |
|---|---|---|---|
| PDF | 1.3000 | 42.35 kg | 32.70 kg |
| prompt | 1.3646 | 56.89 kg | 42.55 kg |

The choice moves the answer by roughly `30 %` in mass, which is far too large to
leave implicit. Note the prompt row's larger mass despite its higher `c`: its
`Sigma_t` is lower, so its mean free path is longer and the sphere is physically
bigger even though it is fewer mean free paths across.

`report_q5` also checks `Sigma_f + Sigma_c + Sigma_s` against the tabulated
`Sigma_t` for every row it prints and warns on a mismatch, so a future typo in
the table is caught rather than propagated.
