# 2D Power Particles examples

Runs the 2D examples of the original Power Particles code
(`~/Documents/Projects/Before2015/pp/data/2D`) with any of our integrators
(`src/core/integrator.h`), to reproduce Power Particles and compare it with
the constrained, softened and GM schemes.

For now only the fluid-only cases, in which the fluid fills the domain, are
supported: no free surface (air particles), no solids, no viscosity
(`notes/viscosity.md`).

## Files

- `runPP2D.cpp`: one run. It reads PP's `.line` domain, `.pp`/`.dat`
  particles and `param.txt`, steps with CFL substeps as `Sim::advance` does,
  and writes per frame:
  - `frame_NNNN.csv`: the full state, x, y, vx, vy, ux, uy, w, p, volume,
    and the cell's cx, cy, cellVolume;
  - one row of `energy.csv` (`FlowMeasures`: energies, momentum, angular
    momentum, plus substeps, Newton iterations and errors).
  `--no-csv` skips the frame CSVs and `--no-energy` skips `energy.csv`; both
  are written by default, and `run.txt` always is. `render.py` needs the
  frame CSVs. Target `runPP2D` in `CMakeLists.txt`.
- `run_examples.py`: runs the 4 fluid-only cases (`--list`) with the chosen
  integrators.

## Usage

```sh
cmake --build build --target runPP2D
python3 benchmarks/pp2d/run_examples.py --frames 100            # PP, softened, GM
python3 benchmarks/pp2d/run_examples.py --case yesweights-band2 --integrator powerparticles
build/runPP2D <dir>/square.line <dir>/uniform_5k.pp <dir>/param.txt out/ --integrator gm --frames 50
```

## Cases

| case | parameters |
|---|---|
| `yesweights-{uniform,band2,radial}`, `yesweights-band2-4vort` | dt 0.01, cap_tol 1e-6, cap_max_iters 5 |

All cases have 5000 particles in [−0.5, 0.5]², no gravity and unit density.
The `noweights` examples are skipped. They use the same particles with
cap_max_iters 0, i.e. Power Particles without the weight solve, so their
volumes drift, and none of our schemes reproduces that. The stored target volumes are rounded (they sum to
1 + 1.2e-8), so they are rescaled to the domain area.

## How each integrator starts

- `powerparticles`: from the loaded state, as the 2015 code does.
- `softened`: its own `Initialize()`, which first moves the particles to a centroidal
  diagram (Lloyd, `--lloydtol`, 1e-6 in `run_examples.py`), so its frame 0 no
  longer matches pp2.
- `gm`: its own `Initialize()`, at the loaded positions.
- `constrained`: not run by default. Its Lloyd start does not reach 1e-9 on
  5000 particles, and with a looser start (`--lloydtol 1e-6`) its first step
  reaches the iteration cap.

## Comparison with the original pp2

`compare.py` compares a `runPP2D` run with the original pp2 frame by frame,
particle by particle: positions, velocities, kinetic energy and cell-volume
errors. It writes `compare.csv` and `compare.png`, and with `--frames` it
renders side-by-side diagrams.

```sh
external/bourne/.venv/bin/python benchmarks/pp2d/compare.py \
    benchmarks/pp2d/results/yesweights-radial/powerparticles <pp2_out_dir> --plot --frames 0 300
```

### Building pp2 (headless, voro++)

`pp2/` holds:
- `build_pp2.sh`: builds `bin/pp2` from the original sources;
- `params/yesweights.txt`: the param file switched to `library: Voro++`.

The pp2 runs are not kept (about 1.2 GB per 300-frame case). Regenerate
them as needed, which takes about 20 s per case.

```sh
sh benchmarks/pp2d/pp2/build_pp2.sh
benchmarks/pp2d/pp2/bin/pp2 -r -i <data>/yesweights/radial_5k.pp -d <data>/yesweights/square.line \
    -p benchmarks/pp2d/pp2/params/yesweights.txt -n 300 -f 1 -o <out_dir>
```

pp2 has no GUI here, uses the bundled voro++, and needs CGAL only for its
2D kernel types (header-only; `external/cgal-6.0.1`) plus Boost headers,
taken from `/Applications/RV64.app/Contents/include`. Compile the sources
in a scratch copy:

    c++ -std=c++17 -O3 -ffp-contract=off -DNDEBUG -DCGAL_NDEBUG \
        -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION -DCGAL_DISABLE_GMP -DCGAL_NO_CORE \
        -I. -isystem <boost> -I<cgal>/*/include ...

then link all objects with `-framework Accelerate`.

`-ffp-contract=off` is required. With FMA contraction (clang's default on
arm64), voro++'s radical-cell bounds test is evaluated with an extra
rounding. It then stops the cell search too early whenever the container has
more than one block, so the cells overlap (areas summing to 428 instead of
1), neighbors are one-sided, and pp2's linear solves diverge.

The param files say `library: CGAL`. Run them with `library: Voro++`
instead: the same algorithm, with cells from voro++.

### Results (300 frames, dt 0.01)

| case | max \|Δq\| / h at frame 150 / 300 | kinetic energy (rel. diff.) | max volume error, ours / pp2 |
|---|---|---|---|
| band2 | 1.6e-3 / 0.20 | ≤ 9e-7 | 2.3e-2 / 2.3e-2 |
| radial | 5.6e-3 / 0.84 | ≤ 1.2e-6 | 2.0e-2 / 2.0e-2 |
| uniform | 4.5e-3 / 0.33 | ≤ 9e-6 | 3.8e-3 / 3.8e-3 |

The trajectories agree to 1e-4 h for the first ~100 frames, then separate
slowly. The two codes compute cells with different libraries (geogram,
voro++) and solve with different linear solvers (Cholesky, PCG), and the
round-off grows with the flow's mixing. The kinetic-energy curves and the
volume errors match, and the frame-300 images are indistinguishable.

