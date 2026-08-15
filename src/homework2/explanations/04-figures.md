# 04 — The Figures

**Three constants in `plots.py` decide what the Q3 comparison actually shows, and the y-limit
heuristic is the one that would break first.**

## Framing the panel

The figures are log-scale over six times spanning `t = 1` to `t = 15`, so nothing can be
autoscaled by matplotlib without either flattening the early panels or clipping the late ones.
`_y_limits` frames **the exact curve only**:

    top    = nanmax(exact) * Y_PAD_ABOVE          # 3.0
    bottom = min(nanmin(exact) * Y_PAD_BELOW,     # 0.25
                 top * 10^-Y_MIN_DECADES)         # at least 2 decades

The `Y_MIN_DECADES` floor is what stops an early-time panel — where the exact solution barely
varies across the front — from being blown up into a band of noise. The diffusion curves are
deliberately *not* consulted: they leak past the causal front and can run orders of magnitude
below the exact curve, and letting them set the limits would compress the region of interest to
nothing. They are allowed to run off the bottom of the panel.

`X_REACH = 1.25` extends the grid 25% past the front for the same reason — the point of the
figure is partly that diffusion puts flux where transport cannot, so that region has to be
visible.

## The exact curve is masked, the diffusion curves are not

`_panel_curves` returns `np.where(x < t, exact, np.nan)` for the exact flux and the raw arrays
for the other two. The `nan` makes matplotlib break the line at the causal front rather than
drawing it to zero, so the front reads as a boundary of support instead of a steep decay. The
grey `axvline` at `x = t` marks it.

Note the mask here is `x < t`, on a grid that starts at `0` — the figures show only the right
half, since the solution is even.

## Practical notes

- **`close(fig)` is called explicitly** at the end of `plot_comparison_for_c`. Five figures per
  run is enough to matter, and `savefig` deliberately does not close, so a caller can keep
  drawing.
- **`savefig` deletes before writing.** Overwriting a PDF in place fails intermittently on this
  machine with `OSError: [Errno 22]`; any new writer must do the same (see the repo's
  `CLAUDE.md`).
- **`make_grid` is three-column** and `label_grid` removes the unused panels afterwards, so
  `TIMES` need not be a multiple of three.
