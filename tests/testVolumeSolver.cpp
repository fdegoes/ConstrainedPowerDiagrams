// Standalone checks of the volume-constrained power diagram solver and of
// the max-entropy feasibility solver.

#include "core/domain.h"
#include "core/feasibleVolumes.h"
#include "core/powerCells.h"
#include "core/powerDiagram.h"
#include "core/tetrahedralize.h"
#include "core/triangulate2D.h"
#include "core/dualSolver.h"
#include "core/jointDual.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace
{

using PP::Domain;
using PP::DualSolver;
using PP::JointDual;
using Vector = JointDual::Vector;

int failures = 0;

void
_Check(bool ok, const std::string& what)
{
    std::printf("[%s] %s\n", ok ? " OK " : "FAIL", what.c_str());
    if (not ok) {
        ++failures;
    }
}

bool
_MakeDomain2D(const std::vector<double>& loops, const std::vector<int>& sizes, Domain& D)
{
    std::vector<double> pts;
    std::vector<int> tris;
    std::string err;
    return PP::Triangulate2D(loops, sizes, pts, tris, err) and Domain::Make2D(pts, tris, D, err);
}

bool
_MakeUnitCube(Domain& D)
{
    std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                               0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    std::vector<int> tris = {0, 2, 1, 0, 3, 2,  4, 5, 6, 4, 6, 7,
                             0, 1, 5, 0, 5, 4,  2, 3, 7, 2, 7, 6,
                             1, 2, 6, 1, 6, 5,  0, 4, 7, 0, 7, 3};
    std::vector<double> tetPts;
    std::vector<int> tets;
    std::string err;
    return PP::Tetrahedralize(pts, tris, PP::TetrahedralizeOptions(), tetPts, tets, err) and
           Domain::Make3D(tetPts, tets, D, err);
}

std::vector<double>
_RandomSites(int n, int dim, double lo, double hi, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(lo, hi);
    std::vector<double> s(n * dim);
    for (double& x : s) {
        x = U(rng);
    }
    return s;
}

// Random positive fractions summing to 1 (between 0.5/n and 1.5/n).
std::vector<double>
_RandomFractions(int n, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0.5, 1.5);
    std::vector<double> f(n);
    double sum = 0.0;
    for (double& x : f) {
        sum += (x = U(rng));
    }
    for (double& x : f) {
        x /= sum;
    }
    return f;
}

void
_TestDomain()
{
    Domain sq, L, cube;
    _Check(_MakeDomain2D({0, 0, 1, 0, 1, 1, 0, 1}, {4}, sq) and
              std::fabs(sq.GetVolume() - 1.0) < 1e-14 and
              std::fabs(sq.GetCentroid()[0] - 0.5) < 1e-14 and
              std::fabs(sq.GetCentroid()[1] - 0.5) < 1e-14,
          "Domain unit square: area 1, centroid (0.5, 0.5)");
    // L = [0,2]x[0,1] U [0,1]x[1,2]: area 3, centroid (5/6, 5/6).
    _Check(_MakeDomain2D({0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2}, {6}, L) and
              std::fabs(L.GetVolume() - 3.0) < 1e-14 and
              std::fabs(L.GetCentroid()[0] - 5.0 / 6.0) < 1e-14 and
              std::fabs(L.GetCentroid()[1] - 5.0 / 6.0) < 1e-14,
          "Domain L-shape: area 3, centroid (5/6, 5/6)");
    _Check(_MakeUnitCube(cube) and std::fabs(cube.GetVolume() - 1.0) < 1e-14 and
              std::fabs(cube.GetCentroid()[2] - 0.5) < 1e-14,
          "Domain unit cube: volume 1, centroid (0.5, 0.5, 0.5)");
}

// Finite-difference checks of the derivatives of F on one domain.
void
_TestDerivatives(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    std::vector<double> sites = _RandomSites(n, dim, 0.05, 0.95, 11);
    std::vector<double> targets(n, D.GetVolume() / n);
    const JointDual F(D, sites, sites, targets, HUGE_VAL, JointDual::Unknowns::Weights);

    std::mt19937 rng(5);
    std::uniform_real_distribution<double> U(-0.01, 0.01);
    Vector w(n);
    for (int i = 0; i < n; ++i) {
        w[i] = U(rng);
    }
    auto build = [&](const Vector& x, PP::ShiftedDiagram& d) {
        std::string err;
        F.Build(x, d, err);
    };

    PP::ShiftedDiagram d;
    build(w, d);
    const PP::PowerCells& cells = d.GetCells();
    double total = 0.0;
    for (double v : cells._volumes) {
        total += v;
    }
    _Check(std::fabs(total - D.GetVolume()) < 1e-10 * D.GetVolume(),
          name + ": cell volumes sum to the domain volume");

    // Centroids vs. a direct check: sum V_i c_i = V_D c_D.
    double centroidErr = 0.0;
    for (int c = 0; c < dim; ++c) {
        double m = 0.0;
        for (int i = 0; i < n; ++i) {
            m += cells._volumes[i] * cells._centroids[dim * i + c];
        }
        centroidErr = std::max(centroidErr, std::fabs(m - D.GetVolume() * D.GetCentroid()[c]));
    }
    _Check(centroidErr < 1e-10, name + ": cell centroids satisfy sum V_i c_i = V_D c_D");

    // Gradient vs. central differences of the value.
    Vector g;
    F.GetGradient(w, d, g);
    double gradErr = 0.0, gradScale = 0.0;
    const double h = 1e-6;
    for (int i = 0; i < n; ++i) {
        Vector wp = w, wm = w;
        wp[i] += h;
        wm[i] -= h;
        PP::ShiftedDiagram dp, dm;
        build(wp, dp);
        build(wm, dm);
        const double fd = (F.GetValue(wp, dp) - F.GetValue(wm, dm)) / (2 * h);
        gradErr = std::max(gradErr, std::fabs(fd - g[i]));
        gradScale = std::max(gradScale, std::fabs(g[i]));
    }
    _Check(gradErr < 1e-6 * std::max(1.0, gradScale) * D.GetVolume(),
          name + ": gradient Vbar_i - V_i matches finite differences of F (err " + std::to_string(gradErr) + ")");

    // -Hessian vs. central differences of the gradient.
    JointDual::SparseMatrix N;
    F.GetNegatedHessian(d, 0.0, N);
    const Eigen::MatrixXd Nd = Eigen::MatrixXd(N);
    double hessErr = 0.0, hessScale = Nd.cwiseAbs().maxCoeff();
    for (int j = 0; j < n; ++j) {
        Vector wp = w, wm = w;
        wp[j] += h;
        wm[j] -= h;
        PP::ShiftedDiagram dp, dm;
        build(wp, dp);
        build(wm, dm);
        Vector gp, gm;
        F.GetGradient(wp, dp, gp);
        F.GetGradient(wm, dm, gm);
        const Vector col = (gp - gm) / (2 * h);
        hessErr = std::max(hessErr, (col + Nd.col(j)).lpNorm<Eigen::Infinity>());
    }
    _Check(hessErr < 1e-4 * hessScale,
          name + ": Hessian (1/2 A_ij / |q_i - q_j|) matches finite differences of the volumes (rel err " +
              std::to_string(hessErr / hessScale) + ")");
}

// Two sites split the unit square / cube by the plane x = 0.5: the shared
// facet is the unit segment / square centered at the domain center.
void
_TestFacetMoments(const Domain& D, const std::string& name)
{
    const int dim = D.GetDimension();
    std::vector<double> sites(2 * dim, 0.5);
    sites[0] = 0.25;
    sites[dim] = 0.75;
    PP::PowerCells cells;
    std::string err;
    const bool ok = PP::ComputePowerCells(D, sites, {0.0, 0.0}, cells, err);
    _Check(ok and cells._facets.size() == 1, name + ": two cells share one facet");
    if (not ok or cells._facets.size() != 1) {
        return;
    }
    const PP::PowerCells::Facet& f = cells._facets[0];

    double centroidErr = 0.0;
    for (int c = 0; c < dim; ++c) {
        centroidErr = std::max(centroidErr, std::fabs(f._centroid[c] - 0.5));
    }
    _Check(std::fabs(f._area - 1.0) < 1e-12 and centroidErr < 1e-12,
          name + ": facet measure 1 and centroid at the domain center");

    // Central second moment of the unit segment / square in the plane x = 0.5:
    // 0 along x, 1/12 along each in-facet axis, no cross terms.
    double momentErr = 0.0;
    for (int r = 0; r < dim; ++r) {
        for (int c = 0; c < dim; ++c) {
            const double expected = (r == c and r > 0) ? 1.0 / 12.0 : 0.0;
            momentErr = std::max(momentErr, std::fabs(f._secondMoment[dim * r + c] - expected));
        }
    }
    _Check(momentErr < 1e-12, name + ": facet second moment is diag(0, 1/12" +
                                 std::string(dim == 3 ? ", 1/12" : "") + ")");
}

void
_TestSolve(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    std::vector<double> sites = _RandomSites(n, dim, 0.05, 0.95, 21);
    std::vector<double> fractions = _RandomFractions(n, 3);
    std::vector<double> targets(n);
    for (int i = 0; i < n; ++i) {
        targets[i] = fractions[i] * D.GetVolume();
    }

    DualSolver::Parameters P;
    P._volumeTolerance = 1e-8;
    PP::VolumeSolver solver(P);
    std::vector<double> w;
    const DualSolver::Result r = solver.Solve(D, sites, targets, w);
    // The tolerance is relative to each target volume.
    double err = 0.0;
    for (int i = 0; i < n; ++i) {
        err = std::max(err, std::fabs(r._diagram.GetVolume(i) - targets[i]) / targets[i]);
    }
    _Check(r.HasConverged() and err <= 1e-8,
          name + ": prescribed volumes reached (" + std::to_string(r._iterations) +
              " iterations, max |V_i - Vbar_i| / Vbar_i " + std::to_string(err) + ", " +
              DualSolver::GetStatusString(r._status) + (r._error.empty() ? "" : ", " + r._error) + ")");
}

// The solved weights reproduce the target areas through ComputePowerDiagram2D.
void
_TestRoundTrip()
{
    std::vector<double> loops = {0, 0, 3, 0, 3, 3, 2, 3, 2, 1, 1, 1, 1, 3, 0, 3};  // U-shape
    std::vector<double> triPts;
    std::vector<int> tris;
    std::string err;
    PP::Triangulate2D(loops, {8}, triPts, tris, err);
    Domain D;
    Domain::Make2D(triPts, tris, D, err);

    const int n = 30;
    std::vector<double> sites = _RandomSites(n, 2, 0.05, 2.95, 31);
    // Keep only sites inside the U (outside the notch [1,2]x[1,3]).
    for (int i = 0; i < n; ++i) {
        if (sites[2 * i] > 1.0 and sites[2 * i] < 2.0 and sites[2 * i + 1] > 1.0) {
            sites[2 * i + 1] = 0.5 * sites[2 * i + 1] - 0.45;
        }
    }
    std::vector<double> fractions = _RandomFractions(n, 4);
    std::vector<double> targets(n);
    for (int i = 0; i < n; ++i) {
        targets[i] = fractions[i] * D.GetVolume();
    }

    // Shifted sites, as the SOP does: p = P + u, omega = w + |u|^2.
    std::vector<double> u = _RandomSites(n, 2, -0.05, 0.05, 32);
    std::vector<double> p(2 * n);
    for (int k = 0; k < 2 * n; ++k) {
        p[k] = sites[k] + u[k];
    }

    PP::VolumeSolver solver;
    std::vector<double> omega;
    const DualSolver::Result r = solver.Solve(D, p, targets, omega);
    _Check(r.HasConverged(), "U-shape with shifts: solver converges (" +
                                std::to_string(r._iterations) + " iterations)");

    std::vector<double> w(n), omega2(n);
    for (int i = 0; i < n; ++i) {
        w[i] = omega[i] - (u[2 * i] * u[2 * i] + u[2 * i + 1] * u[2 * i + 1]);
        omega2[i] = w[i] + (u[2 * i] * u[2 * i] + u[2 * i + 1] * u[2 * i + 1]);
    }
    PP::PolyMesh cells;
    PP::ComputePowerDiagram2D(p, omega2, triPts, tris, cells, err);
    std::vector<double> areas(n, 0.0);
    size_t first = 0;
    for (size_t poly = 0; poly < cells.GetNumPolygons(); ++poly) {
        const int nv = cells._polySizes[poly];
        double a = 0.0;
        for (int k = 0; k < nv; ++k) {
            const double* P0 = &cells._points[2 * cells._polyVerts[first + k]];
            const double* P1 = &cells._points[2 * cells._polyVerts[first + (k + 1) % nv]];
            a += P0[0] * P1[1] - P1[0] * P0[1];
        }
        areas[cells._polyCell[poly]] += 0.5 * a;
        first += nv;
    }
    double err2 = 0.0;
    for (int i = 0; i < n; ++i) {
        err2 = std::max(err2, std::fabs(areas[i] - targets[i]) / D.GetVolume());
    }
    _Check(err2 < 1e-8, "U-shape with shifts: PP Power Diagram with w = omega - |u|^2 "
                       "reproduces the target areas");
}

void
_TestFeasible(const Domain& D, const std::string& name, int n)
{
    const int dim = D.GetDimension();
    std::vector<double> sites = _RandomSites(n, dim, 0.05, 0.95, 41);
    PP::FeasibleVolumeSolver solver;
    std::vector<double> f;
    PP::FeasibleVolumeSolver::Result r = solver.Solve(sites, D.GetCentroid(), dim, f);
    double sum = 0.0;
    std::vector<double> m(dim, 0.0);
    bool positive = true;
    for (int i = 0; i < n; ++i) {
        sum += f[i];
        positive = positive and f[i] > 0.0;
        for (int c = 0; c < dim; ++c) {
            m[c] += f[i] * sites[dim * i + c];
        }
    }
    double cErr = 0.0;
    for (int c = 0; c < dim; ++c) {
        cErr = std::max(cErr, std::fabs(m[c] - D.GetCentroid()[c]));
    }
    _Check(r._feasible and positive and std::fabs(sum - 1.0) < 1e-9 and cErr < 1e-9,
          name + ": feasible fractions are positive, sum to 1 and match the centroid (" +
              std::to_string(r._iterations) + " iterations)");

    // Feasible fractions are valid volume targets: the solver reaches them.
    std::vector<double> targets(n);
    for (int i = 0; i < n; ++i) {
        targets[i] = f[i] * D.GetVolume();
    }
    PP::VolumeSolver vs;
    std::vector<double> w;
    const DualSolver::Result vr = vs.Solve(D, sites, targets, w);
    _Check(vr.HasConverged(), name + ": feasible fractions are reached by the volume solver");

    // All sites on one side of the centroid: infeasible.
    std::vector<double> side = sites;
    for (int i = 0; i < n; ++i) {
        side[dim * i] = D.GetCentroid()[0] + 0.1 + 0.3 * (sites[dim * i] - 0.05);
    }
    PP::FeasibleVolumeSolver::Result r2 = solver.Solve(side, D.GetCentroid(), dim, f);
    _Check(not r2._feasible, name + ": sites on one side of the centroid are reported infeasible");
}

} // namespace

int
main()
{
    _TestDomain();

    Domain square, cube;
    _MakeDomain2D({0, 0, 1, 0, 1, 1, 0, 1}, {4}, square);
    _MakeUnitCube(cube);

    _TestDerivatives(square, "2D square", 25);
    _TestDerivatives(cube, "3D cube", 20);
    _TestFacetMoments(square, "2D square");
    _TestFacetMoments(cube, "3D cube");

    _TestSolve(square, "2D square", 100);
    _TestSolve(cube, "3D cube", 60);
    _TestRoundTrip();

    _TestFeasible(square, "2D square", 40);
    _TestFeasible(cube, "3D cube", 40);

    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
