// Standalone checks of the Newton solvers of latex/appendix-algorithms.tex:
// F (volumes, checked against exploragram's semi-discrete optimal
// transport), H_alpha, H, and H with the linear constraints of the
// constrained integrator.

#include "core/domain.h"
#include "core/dualSolver.h"
#include "core/geogramContext.h"
#include "testUtils.h"

#include <exploragram/optimal_transport/optimal_transport_2d.h>
#include <exploragram/optimal_transport/optimal_transport_3d.h>
#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace
{

using PP::Domain;
using PP::DualSolver;
using PPTest::Check;
using PPTest::FormatValue;

std::string
_Describe(const DualSolver::Result& r)
{
    return std::to_string(r._iterations) + " iterations, residual " + FormatValue(r._residual) + ", " +
           DualSolver::GetStatusString(r._status) + (r._error.empty() ? "" : ", " + r._error);
}

// max_i |V_i - Vbar_i| / |domain|.
double
_GetVolumeError(const DualSolver::Result& r, const std::vector<double>& volumes, double D)
{
    double e = 0.0;
    for (size_t i = 0; i < volumes.size(); ++i) {
        e = std::max(e, std::fabs(r._diagram.GetVolume(int(i)) - volumes[i]) / D);
    }
    return e;
}

// max_i |c_i - cbar_i|.
double
_GetCentroidError(const DualSolver::Result& r, const std::vector<double>& targets, int dim)
{
    double e = 0.0;
    for (int i = 0; i < r._diagram.GetNumSites(); ++i) {
        for (int c = 0; c < dim; ++c) {
            e = std::max(e, std::fabs(r._diagram.GetCentroid(i)[c] - targets[dim * i + c]));
        }
    }
    return e;
}

// Laguerre centroids of the uniform-volume semi-discrete transport, by
// exploragram (Newton, as in its Euler simulator).
std::vector<double>
_ComputeExploragramCentroids(Domain& D, const std::vector<double>& sites)
{
    const int dim = D.GetDimension();
    const GEO::index_t n = GEO::index_t(sites.size() / dim);
    std::vector<double> centroids(sites.size(), 0.0);
    PP::GeogramLock lock;
    GEO::Mesh& M = D.GetMesh();
    if (dim == 2) {
        // exploragram expects 3 coordinates in 2D (z = 0).
        GEO::compute_Laguerre_centroids_2d(&M, n, sites.data(), centroids.data(), nullptr, false, 0, nullptr,
                                           0, 0.0, nullptr, nullptr, 1000);
    } else {
        GEO::compute_Laguerre_centroids_3d(&M, n, sites.data(), centroids.data(), nullptr, false, 1000);
        // exploragram leaves the mesh in dimension 3, but Domain keeps it
        // lifted to dimension 4 (extra coordinate 0).
        M.vertices.set_dimension(4);
    }
    return centroids;
}

void
_TestVolume(Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    const std::vector<double> sites = PPTest::RandomSites(n, dim, 0.05, 0.95, 21);
    const std::vector<double> volumes = PPTest::ScaleFractions(PPTest::RandomFractions(n, 3), D.GetVolume());
    std::vector<double> w;
    PP::VolumeSolver solver;
    const DualSolver::Result r = solver.Solve(D, sites, volumes, w);
    const double err = _GetVolumeError(r, volumes, D.GetVolume());
    Check(r.HasConverged() and err < 1e-9,
          name + " F: prescribed volumes reached (" + _Describe(r) + ", max rel err " + FormatValue(err) + ")");

    // Uniform volumes: same Laguerre centroids as exploragram.
    const std::vector<double> uniform(n, D.GetVolume() / n);
    std::vector<double> wu;
    const DualSolver::Result ru = solver.Solve(D, sites, uniform, wu);
    const std::vector<double> ref = _ComputeExploragramCentroids(D, sites);
    double cErr = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            cErr = std::max(cErr, std::fabs(ru._diagram.GetCentroid(i)[c] - ref[dim * i + c]));
        }
    }
    // exploragram stops at a 1% relative cell-measure error (epsilon 0.01).
    Check(ru.HasConverged() and cErr < 5e-3,
          name + " F: uniform-volume centroids match exploragram (max diff " + FormatValue(cErr) + ")");
}

// Half the bound of eq:regularized-distinct, min_ij |q_i - q_j| / (2 diam),
// which keeps the generators of the maximizer distinct.
double
_GetDistinctAlpha(const Domain& D, const std::vector<double>& sites)
{
    const int dim = D.GetDimension();
    const int n = int(sites.size() / dim);
    double dmin2 = HUGE_VAL;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d2 = 0.0;
            for (int c = 0; c < dim; ++c) {
                d2 += (sites[dim * i + c] - sites[dim * j + c]) * (sites[dim * i + c] - sites[dim * j + c]);
            }
            dmin2 = std::min(dmin2, d2);
        }
    }
    double diam2 = 0.0;
    for (int c = 0; c < dim; ++c) {
        diam2 += 1.0;  // the test domains fit in the unit box or [0, 2]^2
    }
    const double diam = std::sqrt(diam2) * (D.GetVolume() > 1.5 ? 2.0 : 1.0);
    return 0.25 * std::sqrt(dmin2) / diam;
}

void
_TestSoftened(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    const std::vector<double> sites = PPTest::RandomSites(n, dim, 0.05, 0.95, 31);
    const std::vector<double> volumes(n, D.GetVolume() / n);
    const double alpha = _GetDistinctAlpha(D, sites);
    std::vector<double> u, w;
    PP::SoftenedSolver solver;
    const DualSolver::Result r = solver.Solve(D, sites, sites, volumes, alpha, u, w);
    const double vErr = _GetVolumeError(r, volumes, D.GetVolume());
    Check(r.HasConverged() and vErr < 1e-9, name + " H_alpha: converges with exact volumes (alpha " +
                                                FormatValue(alpha) + ", " + _Describe(r) + ")");

    // Stationary shifts u_i = alpha (cbar_i - c_i) (eq:regularized-gradient-shift).
    double uErr = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            const double expected = alpha * (sites[dim * i + c] - r._diagram.GetCentroid(i)[c]);
            uErr = std::max(uErr, std::fabs(u[dim * i + c] - expected));
        }
    }
    Check(uErr < 1e-8, name + " H_alpha: shifts are alpha (cbar - c) (max err " + FormatValue(uErr) + ")");
}

// Targets realized by a power diagram: the volumes and centroids of the
// cells of random weighted sites.
void
_MakeRealizableTargets(
    const Domain& D,
    int n,
    std::vector<double>& sites,
    std::vector<double>& targets,
    std::vector<double>& volumes,
    std::vector<double>& w)
{
    const int dim = D.GetDimension();
    std::vector<double> generators = PPTest::RandomSites(n, dim, 0.1, 0.9, 41);
    w.assign(n, 0.0);
    PP::VolumeSolver vs;
    const std::vector<double> fractions = PPTest::RandomFractions(n, 42);
    vs.Solve(D, generators, PPTest::ScaleFractions(fractions, D.GetVolume()), w);
    PP::ShiftedDiagram d;
    std::string err;
    d.Build(D, generators, generators, std::vector<double>(generators.size(), 0.0), w, err);
    targets = d.GetCells()._centroids;
    volumes = d.GetCells()._volumes;
    // Sites near, but not at, the generators, so that the start (u = 0,
    // w = the weights of the generators) has no empty cell.
    sites = generators;
    const std::vector<double> noise = PPTest::RandomSites(n, dim, -1e-3, 1e-3, 43);
    for (size_t k = 0; k < sites.size(); ++k) {
        sites[k] += noise[k];
    }
}

void
_TestJoint(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    std::vector<double> sites, targets, volumes;
    std::vector<double> wRef;
    _MakeRealizableTargets(D, n, sites, targets, volumes, wRef);

    // Start of H: F at the generators halfway between the sites and the
    // targets.
    std::vector<double> u, w;
    const DualSolver::Result rF = PP::MakeJointStart(D, sites, targets, volumes, 0.5, u, w);
    Check(rF.HasConverged(), name + " H start: F at the midpoints converges (" + _Describe(rF) + ")");
    const std::vector<double> uStart = u, wStart = w;  // Solve() overwrites u, w

    PP::JointSolver solver;
    const DualSolver::Result r = solver.Solve(D, sites, targets, volumes, u, w);
    const double vErr = _GetVolumeError(r, volumes, D.GetVolume());
    const double cErr = _GetCentroidError(r, targets, dim);
    Check(r.HasConverged() and vErr < 1e-9 and cErr < 1e-8,
          name + " H from F: realizable volumes and centroids are recovered (" + _Describe(r) + ", centroid err " +
              FormatValue(cErr) + ")");

    // Without the stretch limit (it only shortens steps), the same exact
    // solution.
    PP::JointSolver unlimited;
    unlimited.GetParameters()._stretchLimit = 0.0;
    std::vector<double> ub = uStart, wb = wStart;
    const DualSolver::Result rb = unlimited.Solve(D, sites, targets, volumes, ub, wb);
    const double vErrB = _GetVolumeError(rb, volumes, D.GetVolume());
    const double cErrB = _GetCentroidError(rb, targets, dim);
    Check(rb.HasConverged() and vErrB < 1e-9 and cErrB < 1e-8,
          name + " H without stretch limit: realizable volumes and centroids are recovered (" + _Describe(rb) +
              ", centroid err " + FormatValue(cErrB) + ")");

    // MakeJointStart: its multipliers give the power diagram of the
    // generators (1 - s) q + s cbar with the prescribed volumes.
    std::vector<double> u1, w1;
    PP::MakeJointStart(D, sites, targets, volumes, 1.0, u1, w1);
    PP::ShiftedDiagram d1;
    std::string err;
    d1.Build(D, sites, targets, u1, w1, err);
    double pErr = 0.0;
    for (size_t k = 0; k < targets.size(); ++k) {
        pErr = std::max(pErr, std::fabs(d1.GetGenerators()[k] - targets[k]));
    }
    double v1 = 0.0;
    for (int i = 0; i < n; ++i) {
        v1 = std::max(v1, std::fabs(d1.GetVolume(i) - volumes[i]) / volumes[i]);
    }
    Check(pErr < 1e-14 and v1 < 1e-9,
          name + " MakeJointStart(s = 1): generators at the targets with the prescribed volumes");

    // Duplicate sites with zero multipliers give coincident generators,
    // which cannot both have a cell: rejected at the start.
    std::vector<double> dup = sites;
    std::copy(sites.begin() + dim, sites.begin() + 2 * dim, dup.begin());
    std::vector<double> u0(sites.size(), 0.0), w0(volumes.size(), 0.0);
    const DualSolver::Result rd = PP::JointSolver().Solve(D, dup, targets, volumes, u0, w0);
    Check(rd._status == DualSolver::Status::InvalidInput and rd._error.find("coincide") != std::string::npos,
          name + " H: duplicate sites are rejected at the start (" + _Describe(rd) + ")");
}

void
_TestConstrained(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    std::vector<double> sites, targets, volumes;
    std::vector<double> w;
    _MakeRealizableTargets(D, n, sites, targets, volumes, w);
    // Direction: the sites themselves (a centroidal-like start, p^0 = q^0).
    const std::vector<double> direction = sites;
    std::vector<double> u;
    PP::ConstrainedJointSolver solver;
    PP::ConstrainedJointSolver::Gauge gauge;
    const DualSolver::Result r = solver.Solve(D, sites, targets, volumes, direction, u, w, gauge);
    const double vErr = _GetVolumeError(r, volumes, D.GetVolume());
    Check(r.HasConverged() and vErr < 1e-9, name + " constrained H: converges with exact volumes (" + _Describe(r) + ")");

    // C z = 0 at the solution.
    double dil = 0.0;
    std::vector<double> trans(dim, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            dil += volumes[i] * u[dim * i + c] * direction[dim * i + c];
            trans[c] += volumes[i] * u[dim * i + c];
        }
    }
    double cz = std::fabs(dil);
    for (double t : trans) {
        cz = std::max(cz, std::fabs(t));
    }
    Check(cz < 1e-10, name + " constrained H: shifts satisfy the dilation and translation constraints (" +
                          FormatValue(cz) + ")");

    if (not r.HasConverged()) {
        return;
    }
    // Stationarity: the centroids are the targets corrected along the
    // homothety, c_i = cbar_i - s p_i - t (eq:flow-shake-step-corr).
    double sErr = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            const double expected = targets[dim * i + c] - gauge._dilation * direction[dim * i + c] -
                                    gauge._translation[c];
            sErr = std::max(sErr, std::fabs(r._diagram.GetCentroid(i)[c] - expected));
        }
    }
    Check(sErr < 1e-8, name + " constrained H: centroids are cbar - s p - t (max err " + FormatValue(sErr) +
                           ", s = " + FormatValue(gauge._dilation) + ")");

    // Shifts violating the constraints are rejected.
    std::vector<double> bad(sites.size(), 0.01), wb = w;
    const DualSolver::Result rb = solver.Solve(D, sites, targets, volumes, direction, bad, wb, gauge);
    Check(rb._status == DualSolver::Status::InvalidInput,
          name + " constrained H: infeasible initial shifts are rejected");
}

// A strip of n cells in [0,n] x [0,1], whose dual graph is a path: n - 2 cut
// cells, hence n - 2 localized homotheties. From targets slightly off a
// realizable configuration, the correction along the homotheties alone
// cannot reach a realizable one, but it can with the branch gauges of the cut
// cells (latex/flow.tex, eq:flow-line-constraint-arm).
void
_TestConstrainedStrip(int n)
{
    const int dim = 2;
    Domain D;
    PPTest::MakeDomain2D({0, 0, double(n), 0, double(n), 1, 0, 1}, {4}, D);
    std::vector<double> generators(dim * n), w(n, 0.0);
    for (int i = 0; i < n; ++i) {
        generators[dim * i] = 0.5 + i;
        generators[dim * i + 1] = 0.5 + (i % 2 ? 0.08 : -0.06);
    }
    PP::ShiftedDiagram d;
    std::string err;
    d.Build(D, generators, generators, std::vector<double>(generators.size(), 0.0), w, err);
    const std::vector<double> volumes = d.GetCells()._volumes;
    std::vector<std::pair<int, int>> edges;
    for (size_t f = 0; f < d.GetNumFacets(); ++f) {
        edges.emplace_back(d.GetFacet(f)._i, d.GetFacet(f)._j);
    }
    const std::vector<PP::CutCellSide> sides = PP::FindCutCellSides(n, edges);
    Check(int(edges.size()) == n - 1 and int(sides.size()) == n - 2,
          "strip: the dual graph is a path with " + std::to_string(sides.size()) + " cut-cell sides");

    // Sites near the generators, and targets moved off the realizable centroids.
    std::vector<double> sites = generators, targets = d.GetCells()._centroids;
    const std::vector<double> noise = PPTest::RandomSites(n, dim, -1.0, 1.0, 51);
    for (size_t k = 0; k < sites.size(); ++k) {
        sites[k] += 1e-3 * noise[k];
        targets[k] += 1e-2 * noise[(k + 3) % noise.size()];
    }
    const std::vector<double>& direction = sites;
    PP::ConstrainedJointSolver solver;

    // Homothety gauges only: no realizable correction.
    {
        std::vector<double> u, wh = w;
        PP::ConstrainedJointSolver::Gauge gauge;
        const DualSolver::Result r = solver.Solve(D, sites, targets, volumes, direction, u, wh, gauge);
        Check(not r.HasConverged(),
              "strip: without the arm gauges the correction does not converge (" + _Describe(r) + ")");
    }

    // With the branch gauges.
    std::vector<double> u, wa = w;
    PP::ConstrainedJointSolver::Gauge gauge;
    const DualSolver::Result r = solver.Solve(D, sites, targets, volumes, direction, sides, u, wa, gauge);
    const double vErr = _GetVolumeError(r, volumes, D.GetVolume());
    Check(r.HasConverged() and vErr < 1e-9, "strip: with the arm gauges it converges with exact volumes (" +
                                                _Describe(r) + ")");
    if (not r.HasConverged()) {
        return;
    }
    // C z = 0 on every row, including the branches.
    double cz = 0.0, dil = 0.0;
    std::vector<double> trans(dim, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            dil += volumes[i] * u[dim * i + c] * direction[dim * i + c];
            trans[c] += volumes[i] * u[dim * i + c];
        }
    }
    cz = std::max({std::fabs(dil), std::fabs(trans[0]), std::fabs(trans[1])});
    for (const PP::CutCellSide& side : sides) {
        double arm = 0.0;
        for (int i : side._side) {
            for (int c = 0; c < dim; ++c) {
                arm += volumes[i] * u[dim * i + c] * (direction[dim * i + c] - direction[dim * side._cut + c]);
            }
        }
        cz = std::max(cz, std::fabs(arm));
    }
    Check(cz < 1e-10, "strip: shifts satisfy all the constraints (" + FormatValue(cz) + ")");
    // Stationarity: c_i = cbar_i - s p_i - t - sum_B s_B (p_i - p_{c_B}).
    std::vector<double> expected(targets);
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            expected[dim * i + c] -= gauge._dilation * direction[dim * i + c] + gauge._translation[c];
        }
    }
    for (size_t l = 0; l < sides.size(); ++l) {
        for (int i : sides[l]._side) {
            for (int c = 0; c < dim; ++c) {
                expected[dim * i + c] -=
                    gauge._arms[l] * (direction[dim * i + c] - direction[dim * sides[l]._cut + c]);
            }
        }
    }
    double sErr = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int c = 0; c < dim; ++c) {
            sErr = std::max(sErr, std::fabs(r._diagram.GetCentroid(i)[c] - expected[dim * i + c]));
        }
    }
    Check(sErr < 1e-8, "strip: centroids are cbar - s p - t - sum s_B (p - p_c) (max err " + FormatValue(sErr) + ")");
}

} // namespace

int
main()
{
    Domain square, L, cube;
    PPTest::MakeUnitSquare(square);
    PPTest::MakeLShape(L);
    PPTest::MakeUnitCube(cube);

    _TestVolume(square, "2D square", 100);
    _TestVolume(cube, "3D cube", 60);
    _TestSoftened(square, "2D square", 60);
    _TestSoftened(L, "2D L-shape", 60);
    _TestSoftened(cube, "3D cube", 40);
    _TestJoint(square, "2D square", 40);
    _TestJoint(cube, "3D cube", 30);
    _TestConstrained(square, "2D square", 40);
    _TestConstrained(cube, "3D cube", 30);
    _TestConstrainedStrip(4);
    _TestConstrainedStrip(6);
    return PPTest::Report();
}
