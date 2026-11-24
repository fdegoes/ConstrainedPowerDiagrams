# Figures of the paper

Status of every figure, its generating script and its data. "Exists" means
the PDF is in `latex/figures/` (or `latex/`) and referenced in the paper.
"Planned" means it is not made yet. Paths are relative to the repository root.

## 1. Figures in the paper

| figure | label | source | script | data | status |
|---|---|---|---|---|---|
| Fig. 1, introduction | `fig:intro` | `latex/figures/intro.pdf` | `benchmarks/figures/plot_intro.py` (`introFigure.cpp`) | generated | exists |
| alpha ablation | `fig:alpha-ablation` | `latex/figures/alpha-ablation-*.pdf` | `benchmarks/figures/plot_alpha_ablation.py` | `latex/figures/alpha-ablation.json` | exists |
| synthetic gallery (2D) | `fig:synthetic-gallery` | `latex/figures/synthetic/synthetic-noise*-*.pdf` | `benchmarks/synthetic/gallery.py` | `benchmarks/synthetic/results/` | exists |
| time vs n (2D) | `fig:scaling-time` | `latex/figures/synthetic/scaling-time.pdf` | `benchmarks/synthetic/plot_scaling.py` | `benchmarks/overnight/2026-11-20/log-synthetic.csv` | exists; regenerate after the Linux run |
| iterations vs n (2D) | `fig:scaling-iterations` | `latex/figures/synthetic/scaling-iterations.pdf` | `benchmarks/synthetic/plot_scaling.py` | same | exists; regenerate after the Linux run |
| convergence, Bourne cases | `fig:convergence` | `latex/figures/convergence.pdf` | `benchmarks/figures/plot_convergence.py` | `benchmarks/figures/cache/convergence/` | exists; subgradient curve added, regenerate after the Linux run |
| EBSD grain map | `fig:ebsd` | `latex/figures/ebsd-*.pdf` | `benchmarks/figures/plot_ebsd.py` | Bourne EBSD data | exists |
| foam | `fig:foam` | `latex/figures/foam-*.pdf` | `benchmarks/figures/plot_foam.py` | foam photo | exists |
| flow blocks | `fig:flow-blocks` | `latex/figures/flow-blocks.tex` | `latex/figures/flow_blocks.py` | generated | exists |
| nullspace, reducible | `fig:appendix-reducible` and nullspace figures | `latex/figures/nullspace/*.pdf` | `latex/figures/nullspace_*.py` | generated | exists |

## 2. Planned figures

### 2.1 Timing and iterations, 3D synthetic (similar to Figs. 5 and 6)

- Source: `benchmarks/synthetic/results3d/ours.csv` (our max H) and
  `benchmarks/synthetic/results3d/oursG.csv` (method G, run by
  `benchmarks/OVERNIGHT_3D_PLAN.md`).
- Script: not written yet. It should reuse the layout of `plot_scaling.py`
  (one panel per noise level, medians over seeds, 25–75 % bands).
- Blocked on: the G run on the Linux machine, and on the load problem in
  `ours.csv` (see the notes on the outliers; re-time on the idle machine).

### 2.2 Statistics table, our synthetic set vs Bourne's data (2D)

- A table with, for each data set, the number of cases, the sizes, the
  volume types, the noise levels, the iteration and time statistics, and the
  outcome of each method (converged, time cap, iteration cap).
- Sources: `benchmarks/synthetic/results/ours.csv`,
  `benchmarks/synthetic/results/theirs.csv`, and the realizable cases in
  `benchmarks/bourne/results/theirs-realizable/theirs.csv`.
- Not written yet.

### 2.3 Power Particles fluid examples (2D)

Data: the original Power Particles inputs in
`~/Documents/Projects/Before2015/pp/data/2D/`. The runner is
`benchmarks/pp2d/run_examples.py` (`runPP2D`, fluid-only cases).

- **2D, four vortices, three volume distributions.** Only
  `yesweights/band2_5k_4vort.pp` exists. The uniform and radial versions
  of the four-vortex case are not in the data folder, so they would have to
  be generated. Open question: are the three distributions the ones of the
  Power Particles paper (uniform, band, radial), and do we generate the
  missing inputs or drop the uniform and radial cases?
- **2D splash.** `2D/splash/` has `initial_{low,med,high}.pp` and `rect.line`.
  Needs a case entry in `run_examples.py` and a runner call.
- **3D, four examples as in the Power Particles paper.** Candidates in
  `3D/`: `smoke` (Fig. "Smoke", 20k particles), `arc` (Splash de Triomphe,
  600k particles), `blade` (moving obstacle), `corridor` (dam break in a
  U-shaped corridor). Open question: which four.
- **Gap for 3D fluids:** there is no 3D runner. `runPP2D` is 2D only, and
  `src/core/integrator.h` is dimension-generic but has no 3D driver,
  output or renderer. A `runPP3D` target and a 3D renderer are needed
  before these figures can be made.
- **New 3D jet example.** No input exists. It needs a new initial state
  (particles and domain) and a case description. Open question: the jet
  geometry, the particle count, and the boundary conditions.

## 3. Open decisions

1. Which four 3D Power Particles examples go in the paper.
2. Whether the uniform and radial four-vortex cases are generated.
3. The geometry and resolution of the new 3D jet.
4. Whether the 2D statistics table includes the approximate Bourne cases
   (EBSD) or only the synthetic and realizable ones.
