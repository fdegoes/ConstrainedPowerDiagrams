# Figures of the paper

Status of every figure, its generating script and its data. Paths are
relative to the repository root.

Status values:

- **exists**: the PDF is in the paper's figures folder and referenced.
- **stale**: the PDF exists but must be regenerated (new data or code).
- **script ready**: the generating script exists and runs, but the figure
  is not made yet.
- **not started**: nothing written yet.
- **blocked**: needs an input or a runner that does not exist yet (the
  blocker is stated).

## 1. Figures in the paper

| figure | label | source | script | data | status |
|---|---|---|---|---|---|
| Fig. 1, introduction | `fig:intro` | `latex/figures/intro.pdf` | `benchmarks/figures/plot_intro.py` (`introFigure.cpp`) | generated | exists |
| alpha ablation | `fig:alpha-ablation` | `latex/figures/alpha-ablation-*.pdf` | `benchmarks/figures/plot_alpha_ablation.py` | `latex/figures/alpha-ablation.json` | exists |
| synthetic gallery (2D) | `fig:synthetic-gallery` | `latex/figures/synthetic/synthetic-noise*-*.pdf` | `benchmarks/synthetic/gallery.py` | `benchmarks/synthetic/results/` | exists |
| time vs n (2D) | `fig:scaling-time` | `latex/figures/synthetic/scaling-time.pdf` | `benchmarks/synthetic/plot_scaling.py` | `benchmarks/overnight/2026-11-20/log-synthetic.csv` | stale: our timings are under load (see the outlier notes); re-time on the idle Linux machine |
| iterations vs n (2D) | `fig:scaling-iterations` | `latex/figures/synthetic/scaling-iterations.pdf` | `benchmarks/synthetic/plot_scaling.py` | same | exists (iteration counts are not affected by load); regenerate with the new log |
| convergence, Bourne cases | `fig:convergence` | `latex/figures/convergence.pdf` | `benchmarks/figures/plot_convergence.py` | `benchmarks/figures/cache/convergence/` | stale: subgradient curve added in the script; regenerate on Linux (`--rerun`) |
| EBSD grain map | `fig:ebsd` | `latex/figures/ebsd-*.pdf` | `benchmarks/figures/plot_ebsd.py` | Bourne EBSD data | exists |
| foam | `fig:foam` | `latex/figures/foam-*.pdf` | `benchmarks/figures/plot_foam.py` | foam photo | exists |
| flow blocks | `fig:flow-blocks` | `latex/figures/flow-blocks.tex` | `latex/figures/flow_blocks.py` | generated | exists |
| nullspace, reducible | `fig:appendix-reducible` and nullspace figures | `latex/figures/nullspace/*.pdf` | `latex/figures/nullspace_*.py` | generated | exists |

## 2. Figures planned

| figure | source | script | status |
|---|---|---|---|
| 2.1 timing and iterations, 3D synthetic (like Figs. 5 and 6) | `benchmarks/synthetic/results3d/ours.csv` (our H, exists), `oursG.csv` (method G) | not written (reuse `plot_scaling.py` layout) | blocked: G 3D run on Linux (`benchmarks/OVERNIGHT_3D_PLAN.md`), and re-timing of our 3D runs on the idle machine |
| 2.2 statistics table, our synthetic set vs Bourne's 2D data | `benchmarks/synthetic/results/ours.csv`, `theirs.csv`, `benchmarks/bourne/results/theirs-realizable/theirs.csv` | not written | not started |
| 2.3a fluid comparison, 2D (PP, GM, constrained, softened alpha > 0) | `benchmarks/pp2d/run_examples.py` output (`energy.csv`, frames) | not written | blocked: softened and constrained fail at initialization on `yesweights-band2` ("Lloyd iterations did not converge"); PP and GM run. Figure to be redone on Linux |
| 2.3b four vortices, three volume distributions (2D) | `yesweights/band2_5k_4vort.pp` exists | not written | blocked: uniform and radial four-vortex inputs do not exist; band2 has not been run |
| 2.3c splash (2D) | `2D/splash/initial_{low,med,high}.pp`, `rect.line` | needs a case entry in `run_examples.py` | not started (inputs exist) |
| 2.3d four 3D examples as in the Power Particles paper | `3D/smoke`, `3D/arc`, `3D/blade`, `3D/corridor` (inputs exist) | none | blocked: no 3D fluid runner (`runPP2D` is 2D only); the choice of four is open |
| 2.3e new 3D jet example | none | none | blocked: no input; geometry, particle count and boundary conditions to be chosen |

## 3. Notes on the fluid figures

- The runner `benchmarks/pp2d/run_examples.py` already supports the four
  integrators used in the comparison (`--integrator powerparticles softened
  gm constrained`). The softened scheme needs `--alpha` different from the GM
  limit; the default is `0.01`, and a second value (for example `0.1`) is
  not yet chosen.
- GM is `gm` with `--stiffness 1`, the alpha -> 0 limit of the softened
  scheme.
- The Mac test on `yesweights-band2` (20 frames) gave: PP ok (volume error
  0.016), GM ok (6.6e-5), softened and constrained failed at initialization.
  The Lloyd settings (`--lloydtol`, `--lloyditers`) are the first thing to
  check.

## 4. Open decisions

1. Which four 3D Power Particles examples go in the paper.
2. Whether the uniform and radial four-vortex cases are generated.
3. The geometry, resolution and boundary conditions of the new 3D jet.
4. Which alpha values for the softened scheme in the fluid comparison.
5. Whether the 2D statistics table includes the approximate Bourne cases
   (EBSD) or only the synthetic and realizable ones.
