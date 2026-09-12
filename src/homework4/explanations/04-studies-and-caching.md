# 04 — The Studies and Their Cache

**`analysis.py` caches every result to CSV, so a second run re-plots without re-simulating.
That is a feature until you change the physics.**

## The three studies

| function | report | configurations | swept |
|---|---|---|---|
| `run_density_study` | Q2 | 11 | `rho` linear on `[10, 100]` g/cm³, probabilities fixed |
| `run_probability_study` | Q3 | 13 | `(P_0, P_1, P_2)` triples, `rho = 30` fixed |
| `run_deep_probability_study` | Q4 | 33 | a wider grid of the same, transcribed from the assignment table |

57 critical-radius searches on a cold cache, each up to 100 histories per bisection step.

Both probability studies plot against `c = P_1 + 2 P_2`, the mean number of particles leaving a
collision ([02](02-the-random-walk.md)). The grids deliberately include triples that **share a
`c` but differ in `P_0`** — that is what makes the Q3 and Q4 Discussion questions answerable at
all, so preserve that property if you edit the lists.

## The cache

Each study writes `docs/homework4/data/<name>.csv` and, on the next run, loads it instead of
simulating:

```python
if os.path.exists(data_path):
    data = np.loadtxt(data_path, delimiter=',', skiprows=1)
else:
    ...simulate...
    np.savetxt(data_path, ..., header=header, comments='')
```

**The check is existence only.** It does not compare parameters, so if you change
`SimulationConfig`, the probability grids, `num_histories`, or anything in `particle.py` or
`simulation.py`, the stale CSV is still loaded and the figures silently keep the old numbers.

**To force a re-run, delete the CSV** — `docs/homework4/data/density_study.csv` and the two
probability files. That is the whole invalidation mechanism.

The CSVs are the only durable record of a run, since nothing is seeded
([05](05-limitations.md)); re-simulating gives different numbers. Treat them as data, not as
build artefacts.

## Figure conventions

`analysis.py` does **not** use `homework1.figures` the way Assignments 2 and 3 do — it sets up
matplotlib itself at module import:

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'cm'

to match the LaTeX body font, and each study builds its own two-panel figure directly. Points
are coloured by the swept variable (`cmap='rainbow'`) with black edges, and the mass panel is
log-scaled while the radius panel is linear — mass spans decades where radius does not.

`BASE_DIR` is resolved by walking three `dirname` calls up from `__file__`, so the output paths
are absolute and `analysis.py` works from any working directory. This is unlike the rest of the
repo, where figure paths are relative to the repository root and the working directory
matters — worth knowing if you copy code between assignments.

Figures go to `docs/homework4/figs/` as PDF with `transparent=True`. Note that these writers do
**not** delete before writing, which the repo's `CLAUDE.md` requires as a workaround for
intermittent `OSError: [Errno 22]` on this machine; if a save ever fails there, that is the fix.
