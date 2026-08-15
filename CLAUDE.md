# TransportHW — working rules

Coursework repository for a neutron transport course. Assignment *N* has its code in
`src/homeworkN/` and its report in `docs/homeworkN/homeworkN.tex`. **`homeworkN` always means
Assignment N, never "question N"** — Assignment 2's Question 3 lives in `src/homework2/` and
in §3 of `docs/homework2/homework2.tex`, not in a `homework3` of its own.

The reports are the user's own write-ups. Additions to them go in
`\textcolor{blue}{...}` / `{\color{blue} ... }`, the convention already established in
`homework2.tex`, so they can be reviewed before being accepted. Do not restructure or
reword the surrounding text.

## Code style

These are the user's standing preferences. They apply to every file written here.

- **Keep it simple and short.** Short files, short functions. If a function is getting long
  or heavily loaded, split the load into helper functions rather than adding nesting.
- **One or two lines of comment on each function** — a docstring saying what it is, not how it
  works. Save the *how* and the *why* for the explanations directory (below).
- Reserve longer in-code comments for things that would otherwise look like mistakes: a
  workaround, a sign convention, a numerically-motivated rewrite.
- Reuse what already exists in the repo before writing a new implementation. In particular
  `homework1.exact_solution` holds the discrete-eigenvalue solvers `compute_nu0_numerical`
  (`c < 1`) and `compute_nu0_magnitude_numerical` (`c > 1`), both already validated.

## Explanations

`src/homeworkN/explanations/` documents **the code**, not the physics. The intended reading
flow is the report first, then these: `docs/homeworkN/homeworkN.tex` carries the idea — every
derivation, every formula, every result — and the explanations carry the map from those
formulas to the modules that implement them.

So these files hold: the module map and the call path, what each non-obvious internal function
does and why it is shaped that way, numerical branches and traps, alternatives that were tried
and rejected, and what was verified and how. They do **not** re-derive anything that belongs in
the report. If a derivation is missing, it goes in the `.tex`, not here.

- One idea per file, named `NN-short-kebab-title.md`, numbered in reading order.
- Each file opens with `# NN — Short Title` and a one-line bold statement of what it settles.
- `00-index.md` is the index: a table of every file with its number, title, and that one-line
  hook. Keep it in sync when adding a file.
- Cite the report by equation or table number instead of restating its content.
- Reference measured numbers rather than asserting correctness.
- Keep them short. A file that repeats the report is worse than no file at all.

## Physics conventions

- Assignment 2 is dimensionless: `Sigma_t = v = 1`. Lengths are mean free paths, times are
  mean free times, and the scalar flux equals the number density.
- `c > 1` is a real case in this course (multiplying media). The discrete eigenvalue turns
  imaginary there; check that any `c`-dependent formula still holds above 1 rather than
  assuming `c < 1`.

## Build and run

- Run a homework from the repository root: `.\.venv\Scripts\python.exe -m homeworkN.main`.
  Figures are written to `docs/homeworkN/figs/` on paths relative to the root, so the working
  directory matters.
- **Build LaTeX with `.\docs\build.ps1 homeworkN`, never latexmk.** latexmk cannot build in
  this repository — see the comment header of `docs/build.ps1` for why. A second argument
  builds a document whose name differs from its directory, e.g. `.\docs\build.ps1 homework2 extra`.
- **Delete a figure or PDF before overwriting it.** Writing over an existing file fails
  intermittently on this machine with `OSError: [Errno 22] Invalid argument`. Both
  `docs/build.ps1` and `homework2/figures.py::savefig` already do this; any new writer must
  too.

## Report style

Reports are LaTeX documents under `docs/`. Match the preamble of the existing ones. When the
user asks for a "minimal" report, that means the formulas and the figures, with no
walk-through prose — but the derivations still belong here rather than in `explanations/`.
