# 06 — The Negative-Flux Fixup, and Why It Never Fires Here

**It is implemented, it works, and in all three questions it activates on exactly zero cell
solves. That is a statement about the problems, not about the code.**

Diamond difference produces `psi_out = 2 psi_i - psi_in`, which turns negative once a cell is
optically thick along the ordinate: with no source, the sign flips at `Sigma_t dx / |mu| > 2`.
The fixup in `cell_flux` is the set-to-zero one — clamp the offending outgoing flux to zero and
re-solve the cell balance with it held there — applied to the spatial and the angular outgoing
face alike, and repeated until nothing is negative. Since clamping every face terminates it,
the loop runs at most once per link.

## It fires when it should

A pure absorber, half-thickness 20 mean free paths on **4 cells** (`dx = 5` mfp), `S_10`, with
the source confined to the first cell: the fixup fires on 5 of 40 cell solves and the scalar
flux stays non-negative, `[0.968, 0.032, 0, 0]`. Without it the third and fourth cells would
oscillate in sign.

## It does not fire in Questions 3 to 5

Counted over whole `k` calculations at `c = 1.5`, `S_10`, at the critical size:

| cells | slab | sphere |
|---|---|---|
| 4 | 0 / 640 | 0 / 1364 |
| 16 | 0 / 2720 | 0 / 6160 |
| 100 | 0 / 17000 | 0 / 38500 |

Two reasons, and both are properties of a criticality eigenvalue problem:

1. **The source is everywhere.** The fission source is proportional to the flux itself, so no
   cell is ever starved the way the shielding test above starves its third cell. The numerator
   of the sweep never gets small enough for `2 psi_i - psi_in` to go under.
2. **The cells are thin.** A critical system is a few mean free paths across, and even 4 cells
   put `dx = 0.15` mfp in the slab — `Sigma_t dx / |mu| = 1.0` at the most grazing `S_10`
   ordinate, still under 2.

The fixup is required by the assignment and is the right thing to carry, but the results of
Questions 3 to 5 would be identical without it. **Where it does matter is that the tables can
be read as pure diamond-difference results**, with no fixup-induced first-order error mixed in.
