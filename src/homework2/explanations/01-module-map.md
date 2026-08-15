# 01 — The Module Map

**Five modules, no solver: every curve in Assignment 2 is a closed form, so `main.py` spends
its time checking identities rather than integrating anything.**

## The files

| file | owns |
|---|---|
| `exact.py` | The Q2 planar flux: Paasschens' `G`, the collided integral, and the Q1 scaling to general `c`. |
| `diffusion.py` | The two diffusion Green's functions, time-dependent and steady, and `D0(c)`. |
| `plots.py` | The Q3 comparison figure: one figure per `c`, one panel per `t`. |
| `figures.py` | Matplotlib helpers — headless backend, three-column grid, safe save. |
| `main.py` | Six verification checks, then the figures. |

## `main.py` is a test suite

Unlike Assignments 1 and 3, there is no `qN.report(figs)` structure, because parts 3(a) and
3(b) have nothing to solve — all three curves being compared are closed forms. So `main.py` is
six `check_*` functions ([05](05-verification.md)) followed by one `generate_figures()`. Each
check is an independent identity; none of them feeds the figures.

The consequence is that **the checks are the only regression tests in the assignment**. There
is no separate test file, and nothing else would notice if `exact.py` broke.

## The call path

    exact.G(w, form)                Paasschens' G: interpolation or exact series
      <- exact.collided_integral    int_0^w0 sqrt(w) G(w) dw, scaled by e^{-w0}
      <- exact.phi_c1               the c = 1 planar flux, report eq. (7)
      <- exact.phi_exact            scaled to any c by report eq. (9)
      <- plots._panel_curves        alongside the two diffusion curves
      <- plots.plot_comparison_for_c

`diffusion.py` sits beside that path, not in it: `_phi_diffusion` is a two-line Gaussian and
everything else in the module is a thin wrapper choosing `D`.

## The `form` argument threads all the way through

`G`, `collided_integral`, `phi_c1` and `phi_exact` each take `form="interp"` or `"series"`,
defaulting to `"interp"`. It selects Paasschens' interpolation `G(w) ~ e^w sqrt(1 + b/w)`
against the exact series. The figures use the default; `check_normalisation` passes
`"series"` so that the normalisation test is not measuring the interpolation's own error, and
`check_interpolation_error` then measures exactly that difference. Keep the argument threaded
if you add a caller — dropping it silently makes any new check untestable against the exact
form.
