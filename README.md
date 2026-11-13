# PowerParticles (PP)

Houdini SOPs for power diagrams clipped to a domain, built on
[geogram](https://github.com/BrunoLevy/geogram), which lives in `external/geogram`.
This is research code. geogram's bundled TetGen is AGPL-licensed, and that is
fine for this project.

| SOP | Operator | Does |
|-----|----------|------|
| PP Triangulate 2D | `pp_triangulate2d` | closed polygons → constrained Delaunay triangles |
| PP Tetrahedralize | `pp_tetrahedralize` | closed triangulated surface → tetrahedra (TetGen) |
| PP Power Diagram | `pp_powerdiagram` | sites + triangle/tet domain → clipped power cells |
| PP Power Diagram Solver | `pp_powerdiagramsolver` | sites + target volume fractions → weights (Newton) |
| PP Feasible Volumes | `pp_feasiblevolumes` | sites as target centroids → feasible volume fractions (max-entropy) |

Typical networks:

    2D:  polygons --> PP Triangulate 2D --> PP Power Diagram (input 2)
    3D:  surface  --> PP Tetrahedralize --> PP Power Diagram (input 2)

    volumes: sites --> PP Feasible Volumes --> PP Power Diagram Solver --> PP Power Diagram
             (each node takes the triangle/tet domain on input 2)

## Build and install

Run these in a regular terminal. The Claude Code sandbox can't write to
`~/Library/Preferences`.

    sh config.sh                 # or: sh config.sh <houdini_user_dir>
    cmake --build build -j8
    ./build/testPowerDiagram     # standalone checks, no Houdini needed
    ./build/testVolumeSolver
    ./build/testJointDual
    ./build/testDualSolvers
    ./build/testIntegrators

The install layout follows `../geomt`, under the Houdini user directory (by
default `~/Library/Preferences/houdini/20.5`):

| File | Installed to |
|------|--------------|
| `PowerParticlesHDK.dylib` (all three SOPs) | `<user_dir>/dso/` |
| `src/houdini/help/nodes/sop/*.txt` | `<user_dir>/help/nodes/sop/` |

CMake asks `hython` for the user directory, the same way
`houdini_get_default_install_dir()` does. If `hython` can't run (for example,
no license), it uses the standard location instead. Set `-DPP_HOUDINI_HOME=...`
to override it.

geogram is compiled as a static library inside this build. The sparse
Cholesky of the Newton systems is SuiteSparse's CHOLMOD on macOS and Linux
(`-DPP_SPARSE_SOLVER=Cholmod`, with `-DSUITESPARSE_DIR=<prefix>`; config.sh
uses `~/MyLibs/SuiteSparse-5.4-arm64` on macOS). On Linux, config.sh looks
for SuiteSparse in `/pixar/d2/sets/vfx2025.002`, then `~/mylibs/SuiteSparse`,
and for Eigen ≥ 3.3 in the same set, then `~/mylibs/eigen337`; override them
with `SUITESPARSE_PATH` and `EIGEN3_PATH`. It also runs that set's cmake
(override with `CMAKE`), so build with `make -C build -j8` if the cmake on
`PATH` is older than 3.12. Set
`PP_SPARSE_SOLVER=Accelerate` (macOS) or `Eigen` to switch backends.

## Datasets and benchmarks

All benchmarks invert power diagrams: given target volumes V̄ and centroids
c̄, find generators and weights whose diagram has them. They need the
`benchBourne` and `makeSynthetic` targets and the Python environment with
pysdot in `external/bourne/.venv`:

```sh
cmake --build build --target benchBourne makeSynthetic
export MPLCONFIGDIR=/tmp/mplcfg
PY=external/bourne/.venv/bin/python
```

### List of datasets

| dataset | source | dim | N | exact solution? | details |
|---|---|---|---|---|---|
| synthetic 2D | ours | 2 | 10–10000 (600 cases) | yes | `benchmarks/synthetic/README.md` |
| synthetic 3D | ours | 3 | 10–10000 (600 cases) | yes | `benchmarks/synthetic/README.md` |
| `2025-invert` | Bourne 2025 §5.1 | 2 | 20 | yes | `benchmarks/bourne/README.md` |
| `2025-perturbed-small` | Bourne 2025 §5.2.1 | 2 | 20 | no (fit only) | `benchmarks/bourne/README.md` |
| `2025-perturbed-large` | Bourne 2025 §5.2.2 | 2 | 20 | no (fit only) | `benchmarks/bourne/README.md` |
| `2025-ebsd` | Bourne 2025 §5.3 | 2 | 243 | no (fit only) | `benchmarks/bourne/README.md` |
| `2026-voronoi` | Bourne 2026 §6.3.1 | 2 | 20 | yes | `benchmarks/bourne/README.md` |
| `2026-aligned` | Bourne 2026 §6.3.1 | 2 | 728 | yes | `benchmarks/bourne/README.md` |
| `2026-ebsd` | Bourne 2026 §6.3.2 | 2 | 243 | no (fit only) | `benchmarks/bourne/README.md` |

Bourne 2025 is *Inverting Laguerre tessellations* (Bourne, Pearce & Roper);
Bourne 2026 is *Semi-discrete convex order and Laguerre tessellation fitting*
(Bourne, Gallouët, Mérigot & Natale). Their code is in `external/bourne/`.

### Synthetic (ours)

Inverse problems in the unit square (2D) or cube (3D), built from a known
power diagram. A case is named `noise<k>-<volumes>-<n>-<seed>`:

- noise 0, 1, 2 (dart-throwing radius 0.8, 0.5, 0.25 of the cell size);
- volumes `uniform` or `voronoi`;
- n one of 20 log-spaced sizes from 10 to 10000;
- seeds 1–5.

Cases are generated the first time a driver needs them, in
`benchmarks/synthetic/results/` (2D) or `results3d/` (3D), and reused
afterwards. All cases have an exact solution, so they use our H solver;
`run_synthetic_g.py` solves them with G instead (sites q = c̄, fixed
weights w = w_Y − |Y − c̄|²).
`run_theirs.py` can be slow: FW runs up to `--max-iter` (100000) or
`--max-seconds` (1800 s) per case, e.g. about 290 s on
`noise2-uniform-10-1` without meeting their gap rule.

```sh
# Our solver on every 2D case (writes results/ours.csv); --dim 3 for 3D
$PY benchmarks/synthetic/run_synthetic.py
# A subset
$PY benchmarks/synthetic/run_synthetic.py --n 10 100 --noise 2 --volumes uniform --seeds 1
# Max G (shifts only) with fixed weights chosen so that lambda = w_Y (writes results/oursG.csv)
$PY benchmarks/synthetic/run_synthetic_g.py
# Tables of the H vs. G report (benchmarks/synthetic/REPORT.md)
$PY benchmarks/synthetic/report.py
# Bourne et al.'s methods (FW, SLSQP; 2D only, seed 1 by default; writes results/theirs.csv)
$PY benchmarks/synthetic/run_theirs.py --method fw slsqp
# Generate every case and print the data statistics
$PY benchmarks/synthetic/stats.py [--dim 3]
# Why one case fails
$PY benchmarks/synthetic/diagnose.py noise2-uniform-1129-3
# Figures: the n = 10 gallery for the paper, or every 2D case for browsing
$PY benchmarks/synthetic/gallery.py
$PY benchmarks/synthetic/render_all.py
```

### Bourne et al. (2025, 2026)

The seven examples of the two papers, generated as in their notebooks. The
realizable ones (exact solution) run with `run_realizable.py`, the others
with `run_approximate.py`. Each runs their method and ours, under their
termination rule and a common residual, and writes JSON and `summary.txt`
to `benchmarks/bourne/results/realizable/` or `.../approximate/`.

```sh
# All cases of a group
$PY benchmarks/bourne/run_realizable.py      # 2025-invert, 2026-voronoi, 2026-aligned
$PY benchmarks/bourne/run_approximate.py     # 2025-perturbed-small, -large, 2025-ebsd, 2026-ebsd
# One dataset
$PY benchmarks/bourne/run_realizable.py --case 2026-aligned
$PY benchmarks/bourne/run_approximate.py --case 2025-ebsd
# List the cases of a driver
$PY benchmarks/bourne/run_realizable.py --list
```

The realizable cases use our H solver. The others have no exact solution
(perturbed centroids, or measured EBSD grains, which are not polygonal), so
they use the softened solver H_alpha over a sweep of stiffnesses α. Their
2026 Frank–Wolfe solves the same problem with α = t − 1.

How each dataset is built, the choice of α, all the options, and the known
issues (currently: our solver does not start on `2025-ebsd`) are in
`benchmarks/bourne/README.md`.

### Overnight comparison (ours vs. Bourne et al.)

`benchmarks/run_overnight.sh` runs our H and their methods on all the
realizable data (the Bourne cases and the synthetic 2D grid), each stopped
on its own rule at 1e-8, and writes the tables and timing plots for the
paper. Run it on an idle machine:

```sh
OVERNIGHT_QUICK=1 bash benchmarks/run_overnight.sh   # smoke test, ~4 min, /tmp/overnight-quick
caffeinate -i bash benchmarks/run_overnight.sh       # full run, worst case ~5 h
```

What it runs, the rules, and the output are in `benchmarks/OVERNIGHT.md`.

## Why the domain must be simplicial

geogram's restricted Voronoi diagram never clips against a polygon or a
surface directly. It walks a simplicial mesh of the domain: triangles in 2D,
tetrahedra in 3D. Each simplex is a convex starting cell that gets cut by the
power planes of nearby sites. The cut pieces of a cell are then glued together
across shared edges or faces, using the simplex adjacency. That is why a
non-convex domain (or one with holes or cavities) works: each simplex is
convex even when the domain isn't. Building that mesh is a separate step,
handled by the two helper SOPs.

## PP Power Diagram

Inputs:

1. **Sites**: points `P`, plus an optional float attribute `w` and an optional
   vector attribute `u`. The diagram uses site `P + u` with weight `w + |u|^2`.
2. **Clipping domain**: chosen by **Mode** (2D or 3D).

### Requirements on the domain

These are **assumed, not checked**. A domain that breaks them gives a wrong
output.

**2D**: triangle polygons in the XY, YZ or ZX plane.
- Conforming: two triangles either don't touch, or share a whole edge or point.
  No T-junctions and no overlaps.
- Welded: neighboring triangles share points.
- No zero-area triangles.
- Winding doesn't matter. The domain may be non-convex, have holes, or have
  several pieces.

**3D**: tetrahedron primitives. Other primitives are ignored.
- Conforming: two tets either don't touch, or share a whole face, edge or
  point. No T-junctions and no overlaps.
- Welded: neighboring tets share points.
- No zero-volume tets.
- Tet orientation doesn't matter. The output winding follows the tet
  orientation, so inverted tets are flipped internally. The domain may be
  non-convex, have cavities, or have several pieces.

### Output

- Closed polygons. Each has a prim int `cell` holding the index of its site.
  Points aren't shared between cells.
- In 2D, one polygon per boundary loop of each cell. Outer boundaries are
  wound counter-clockwise, and boundaries around a hole inside a cell are
  wound clockwise. In 3D, the boundary faces of each cell, with outward
  normals.
- A detail int `nb_empty_cells`: the number of sites with no cell. This
  happens when a site is outside the domain or hidden by heavier neighbors.
- **Boundary Merge Angle** (3D): at the default of 0 the domain boundary is
  exact. A positive value merges nearly coplanar boundary faces of each cell
  into one polygon. That gives fewer polygons but flattens curved domains, and
  on domains with cavities it was measured to give wrong volumes.

## PP Triangulate 2D

- Input: closed polygons in the chosen plane. Each polygon is one boundary
  loop.
- The region follows the even-odd rule, so a loop inside another is a hole.
- Loops must not cross each other or themselves, have no zero-length edges,
  and must not repeat their first point at the end.
- Output: welded, counter-clockwise triangles that use only the polygon points.

## PP Tetrahedralize

- Input: triangles only.
- With **Preprocess Surface** off (the default), the surface must be welded,
  closed, manifold, and free of self-intersections. Turning it on makes geogram
  try to repair the surface first.
- **Refine** adds interior points (TetGen quality meshing). The target is set
  by **Quality** (radius-edge ratio, about 1 to 2).
- Output: tetrahedra with positive volume. The surface triangles are faces of
  the output tets.

## Power-diagram duals (`src/core/jointDual.h`, `src/core/dualSolver.h`)

The algorithms of `latex/appendix-algorithms.tex` are built from four
pieces:

- `ShiftedDiagram`: the power diagram of the shifted sites
  `p = q + u`, `lambda = w + |u|² + 2 uᵀ(q − c̄)` (eq:shifted-site), with the
  cell volumes and centroids and the facet coefficients of the Hessian blocks.
  The per-cell quantities come from `ComputePowerCells` (`src/core/powerCells.h`),
  built with geogram's restricted Voronoi callbacks over the `Domain` mesh.
- `JointDual`: the concave objectives as functions of the multipliers:
  `H` (shifts and weights), `H_α` (softened centroids), `G` (shifts, fixed
  weights) and `F` (weights). It gives the value, the gradient, the relative
  residual (every `|V_i − V̄_i| / V̄_i` and centroid error over `L`) and
  `N = −∇² + η M`, where the mass `M` is `V̄_i r_i²` on the shifts and `V̄_i`
  on the weights (`r_i = V̄_i^{1/d}`), so that `η` has the same units on both
  blocks.
- `DualSolver`: Newton's method on these objectives. `VolumeSolver` (F) takes
  the Newton step with a backtracking line search; the other solvers
  (`CentroidSolver`, `JointSolver`, `SoftenedSolver`, `ConstrainedJointSolver`)
  use Levenberg–Marquardt, `(N + μ diag(N)) p = g`, with Nielsen's update of
  `μ`. Steps are accepted when no cell is empty and the objective increases;
  the solve converges only when the residual is at most the tolerance.
  `ConstrainedJointSolver` borders the step with the linear constraints of the
  constrained integrator (`latex/appendix-bordered.tex`).
- `MakeJointStart`: the start of an `H` solve, an `F` solve at the generators
  `(1 − s) q + s c̄` written as multipliers of `H`. `s = 0.5` is the most robust
  for static solves; the integrators use `s = 0`.

`CentroidalDiagram` builds a centroidal power diagram with prescribed volumes
(Lloyd iterations on `F`), the initial state of `ConstrainedIntegrator` and
`PenaltyIntegrator` (`latex/flow.tex`).

## Feasible volumes (`src/core/feasibleVolumes.h`)

`FeasibleVolumeSolver::Solve(targets, domainCentroid, dim, fractions)` follows
`latex/contributions.tex`, "Remarks / Feasibility". It returns the
maximum-entropy fractions `f_i > 0` with `Σ f_i = 1` and
`Σ f_i c̄_i = c_D`. `MaxEntropyEnergy` is its (d+1)-dimensional convex dual:

- value: `Σ exp(−y·l_i) + y_d`, with `l_i = [c̄_i, 1]`;
- gradient: `[0, 1] − Σ f_i l_i`;
- Hessian: `Σ f_i l_i l_iᵀ`.

The targets are centered on the domain centroid and scaled before solving.
Fractions exist exactly when `c_D` is inside the convex hull of the targets;
otherwise the result is reported as infeasible.

Eigen is header-only. CMake looks for it in `EIGEN3_INCLUDE_DIR`, which
defaults to `~/MyLibs/eigen-3.4.0`, the same copy `../geomt` uses.

## Coding conventions

- **Files** start with a lowercase letter: `dualSolver.h`,
  `powerDiagramSOP.cpp`.
- **Classes and structs** start with a capital letter: `DualSolver`,
  `PowerDiagramSOP`.
- **Methods and free functions** start with a capital-letter verb: `Solve`,
  `Build`, `GetGradient`, `GetNumVariables`, `ComputePowerDiagram2D`,
  `Tetrahedralize`.
- **Private methods and file-local helper functions** (in a `.cpp`'s
  anonymous namespace) start with an underscore: `_Solve`,
  `_ComputeDampedStep`, `_Factorize`, `_ChainLoops`. This includes the private
  parameter accessors of the SOPs (`_MODE()`, `_PLANE()`, ...).
- **Data members** of classes and structs: `_name` (underscore, then
  lowerCamelCase: `_params`, `_facetPts`, `Parameters::_maxIterations`). The
  lowercase letter after the underscore tells them apart from `_Helper`
  functions.
- **Comments** use `//`, including documentation and trailing member
  comments (no `///`, `///<` or Doxygen tags).
- **Braces** around every `if`, `else`, `for` and `while` body, even a
  single line.
- **Logical operators** are spelled `and`, `or` and `not`, never `&&`, `||`
  and `!`. Conditions chained over several lines are split into separate
  `if`s (or a named helper).
- **No exceptions**: no `try`/`catch`. Functions return `bool` plus an error
  string, and inputs are validated (sizes, index ranges, finite values) before
  they reach geogram.
- **Small named helpers**: long functions are split into steps with
  descriptive verb names (`ReadDomain`, `AppendPolygons`, `ComputeCauchyStep`,
  ...). Shared code lives in `core/geometry.h`, `core/geogramMesh.h` and
  `houdini/geometryIO.h`.
- **Overrides of HDK virtuals** keep the HDK's names (`cookMySop`,
  `updateParmsFlags`, `inputLabel`, `myConstructor`). Operator and parameter
  names (`pp_powerdiagram`, `weightattrib`, ...) stay lowercase, as Houdini
  expects.
