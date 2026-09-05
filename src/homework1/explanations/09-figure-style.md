# 09 — The Figure Style

**Every panel is `case_label` + `panel` over a `subplots`/`make_grid` frame; the panel sizes
read backwards and the titles are gone on purpose.**

`figures.py` is Assignment 3's `homework3/figures.py`, ported so both reports look like one
set of figures: serif text, ticks inwards on all four sides, a light solid grid behind the
data, a transparent save, and the palette constants `NAVY, ORANGE, GREEN, RED, BLUE, PURPLE`
that `q1`…`q5` import instead of writing hex.

## The three calls

    fig, axes = subplots(rows, cols)     # or make_grid(n) for one panel per case
    ax.plot(...)                          # colours from the palette constants
    case_label(ax, '$c = 0.7$')           # only where the panel needs naming
    panel(ax, xlabel, ylabel, log=, legend=)
    finish(fig); savefig(fig, path)

## Three things that will catch you out

**The panel sizes are inverted.** `PANEL_WIDTH`/`GRID_WIDTH` are figure inches, and the
report scales the result to `\textwidth`; a *smaller* panel is scaled down less, so its
13 pt label lands larger on the page. Raising these numbers makes the text smaller. The
current values are the previous ones times ~0.7, which is why every figure kept its aspect
ratio and its place on the page while the text grew.

**There is no title, anywhere.** `panel` takes no title and `finish` takes no suptitle — the
report caption carries the description. What a title used to identify moves to one of three
places: the y-label (the quantity), the legend (the curves), or `case_label` (the case —
this panel's `c`, approximation or material), which annotates inside the axes on a
translucent white box and takes an axes-fraction `xy` chosen per figure to clear the data.

**The y-label is rotated, so its length competes with the panel *height*.** About twenty
characters is the ceiling at 13 pt in a 3.2 in panel. Longer labels do not wrap; they run
into the row above, which is what `'Error vs. exact transport (%), $c = 0.5$'` did before it
became `'Relative error (%)'` plus a `case_label`. Keep the y-label to the quantity alone —
`$L_2$ error`, `$a/2$ [mfp]`, `Difference (%)` — and put the qualifier elsewhere.

## Verified

All eighteen figures were regenerated and read at 55 dpi after the port: no label crosses a
panel boundary, no legend hides a curve, and `q4_mesh_convergence` has explicit `set_xticks`
because under two decades matplotlib labels the log minor ticks too and they collide.
