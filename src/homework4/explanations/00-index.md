# 00 — Index

**These files explain the code in `src/homework4/`. The physics and the results are in
`docs/homework4/homework4.tex` — read that first.**

Assignment 4 is an explicit Monte Carlo simulation of Studention transport in a sphere of
Exercisium: sample a flight, advance, test escape, draw a collision outcome, repeat. There is
no equation being discretised, so almost everything worth writing down about this code is a
sampling or bookkeeping decision rather than a derivation.

| # | Title | What it settles |
|---|---|---|
| [01](01-module-map.md) | The Module Map | Which file owns what, and the two entry points. |
| [02](02-the-random-walk.md) | The Random Walk | How a history advances, what a "generation" actually is, and the one copy that matters. |
| [03](03-criticality-search.md) | The Criticality Search | What "critical" means operationally here, and why the bisection is not an ordinary bisection. |
| [04](04-studies-and-caching.md) | The Studies and Their Cache | Why re-running `analysis.py` does not re-simulate, and how to force it to. |
| [05](05-limitations.md) | Known Limitations | What the current code cannot tell you, in order of how much it matters. |
| [06](06-the-analytic-comparison.md) | The Analytic Comparison | How Question 5 gets both closed-form curves out of Assignment 1, and the third one that was dropped. |

Nothing here re-derives the report.

**Status.** `docs/homework4/homework4.tex` still has the Question 1 results as an
`[Insert ...]` placeholder. Questions 2 to 5 are written, and the code produces the four figures
the report references.
