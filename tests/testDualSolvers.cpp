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

    PP::JointSolver solver;
    const DualSolver::Result r = solver.Solve(D, sites, targets, volumes, u, w);
    const double vErr = _GetVolumeError(r, volumes, D.GetVolume());
    const double cErr = _GetCentroidError(r, targets, dim);
    Check(r.HasConverged() and vErr < 1e-9 and cErr < 1e-8,
          name + " H from F: realizable volumes and centroids are recovered (" + _Describe(r) + ", centroid err " +
              FormatValue(cErr) + ")");

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
    return PPTest::Report();
}
