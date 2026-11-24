// The instance of the showcase figure of the introduction (latex/intro.tex,
// fig:intro): a power diagram of the unit square with prescribed volumes and
// centroids, solved by H from given sites or from sites at the centroids.
//
// From generators Y, the power diagram of Y with the volumes (VolumeSolver,
// F, at Y) gives realizable targets (Vbar, cbar). JointSolver (H) then solves
// for the shifts u and weights w of the sites q, and the tool writes the
// generators p = q + u and the shifted weights lambda of the solution.
//
// The generators of a diagram with prescribed volumes and centroids are
// unique only up to a homothety p -> s p + t (s > 0), under which the shifted
// weights become |s p_i + t|^2 - s (|p_i|^2 - lambda_i) and the cells do not
// change. With "closest 1" the tool writes the representative closest to the
// sites, min_{s,t} sum_i Vbar_i |s p_i + t - q_i|^2, instead of the one the
// solver reaches.
//
//   introFigure <in.txt> <out.txt>
//
// in.txt   "n", "generators", "sites" records: the number of cells, then 2 n
//          values each for Y and q; optionally "targets" (2 n values, cbar
//          given directly, instead of "generators"), "volumes" (n values, scaled
//          to sum to the area; uniform by default), "sitesAtTargets 1" (q =
//          cbar, instead of the "sites" record), "closest 1" and "alpha a"
//          (H_alpha instead of H, see below), "domain" (a polygon instead of
//          the unit square: 2 values per vertex of one loop), "centroidOnly m" (G instead of
//          H, see below)
// out.txt  "n", "volumes", "targets", "sites", "generators",
//          "shifts", "weights" (the shifted weights lambda) records

#include "core/domain.h"
#include "core/dualSolver.h"
#include "core/jointDual.h"
#include "core/shiftedDiagram.h"
#include "core/triangulate2D.h"

#include <Eigen/Dense>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../common/boxDomain.h"

namespace
{

// Volumes and centroids of the cells of the generators with weights (as a
// diagram of H with zero shifts).
bool
_GetCells(
    const PP::Domain& domain,
    const std::vector<double>& generators,
    const std::vector<double>& weights,
    std::vector<double>& volumes,
    std::vector<double>& centroids,
    std::string& err)
{
    const size_t n = generators.size() / 2;
    const std::vector<double> ones(n, 1.0 / double(n));
    PP::JointDual dual(domain, generators, generators, ones, HUGE_VAL, PP::JointDual::Unknowns::Weights);
    PP::JointDual::Vector z = PP::JointDual::Vector::Zero(Eigen::Index(n));
    for (size_t i = 0; i < n; ++i) {
        z[Eigen::Index(i)] = weights.empty() ? 0.0 : weights[i];
    }
    PP::ShiftedDiagram diagram;
    if (not dual.Build(z, diagram, err)) {
        return false;
    }
    if (diagram.GetNumEmptyCells() > 0) {
        err = "a cell is empty";
        return false;
    }
    volumes.resize(n);
    centroids.resize(2 * n);
    for (size_t i = 0; i < n; ++i) {
        volumes[i] = diagram.GetVolume(int(i));
        centroids[2 * i] = diagram.GetCentroid(int(i))[0];
        centroids[2 * i + 1] = diagram.GetCentroid(int(i))[1];
    }
    return true;
}

void
_WriteRecord(std::ostream& out, const char* key, const std::vector<double>& v)
{
    out << key;
    for (double x : v) {
        out << " " << x;
    }
    out << "\n";
}

} // namespace

int
main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: introFigure <in.txt> <out.txt>" << std::endl;
        return 2;
    }
    std::ifstream in(argv[1]);
    std::map<std::string, std::vector<double>> records;
    for (std::string line; std::getline(in, line);) {
        std::istringstream ls(line);
        std::string key;
        ls >> key;
        for (double x; ls >> x;) {
            records[key].push_back(x);
        }
    }
    const int n = records["n"].empty() ? 0 : int(records["n"][0]);
    const std::vector<double> given = records["targets"];
    const std::vector<double> Y = given.empty() ? records["generators"] : given;
    const bool sitesAtTargets = not records["sitesAtTargets"].empty() and records["sitesAtTargets"][0] != 0.0;
    std::vector<double> q = records["sites"];
    if (sitesAtTargets) {
        q.assign(Y.size(), 0.0);
    }
    if (n <= 0 or Y.size() != size_t(2 * n) or q.size() != size_t(2 * n)) {
        std::cerr << "cannot read " << argv[1] << std::endl;
        return 2;
    }

    // The unit square, or the polygon of the "domain" record (2 values per
    // vertex, one loop, possibly non-convex).
    std::string err;
    PP::Domain domain;
    double area = 1.0;
    const std::vector<double>& loop = records["domain"];
    if (loop.empty()) {
        if (not PPBench::MakeBoxDomain(2, {0.0, 0.0, 1.0, 1.0}, domain, err)) {
            std::cerr << err << std::endl;
            return 2;
        }
    } else {
        std::vector<double> pts;
        std::vector<int> tris;
        if (loop.size() < 6 or loop.size() % 2 != 0 or
            not PP::Triangulate2D(loop, {int(loop.size() / 2)}, pts, tris, err) or
            not PP::Domain::Make(2, pts, tris, domain, err)) {
            std::cerr << "domain: " << err << std::endl;
            return 2;
        }
        const size_t m = loop.size() / 2;
        area = 0.0;
        for (size_t k = 0; k < m; ++k) {
            const size_t l = (k + 1) % m;
            area += 0.5 * (loop[2 * k] * loop[2 * l + 1] - loop[2 * l] * loop[2 * k + 1]);
        }
        area = std::abs(area);
    }
    std::vector<double> prescribed(size_t(n), area / double(n));
    if (not records["volumes"].empty()) {
        prescribed = records["volumes"];
        double sum = 0.0;
        for (double v : prescribed) {
            sum += v;
        }
        if (prescribed.size() != size_t(n) or not (sum > 0.0)) {
            std::cerr << "cannot read " << n << " volumes from " << argv[1] << std::endl;
            return 2;
        }
        for (double& v : prescribed) {
            v *= area / sum;
        }
    }

    // The targets: given, or those of the diagram of Y with the prescribed
    // volumes.
    std::vector<double> volumes, targets;
    if (not given.empty()) {
        volumes = prescribed;
        targets = given;
    } else {
        PP::DualSolver::Parameters fParams;
        fParams._volumeTolerance = 1e-13;
        fParams._maxIterations = 200;
        std::vector<double> yWeights;
        const PP::DualSolver::Result fRes = PP::VolumeSolver(fParams).Solve(domain, Y, prescribed, yWeights);
        if (not fRes.HasConverged() and not (fRes._volumeResidual <= 1e-10)) {
            std::cerr << "volume solve at Y: " << PP::DualSolver::GetStatusString(fRes._status) << " " << fRes._error
                      << std::endl;
            return 1;
        }
        if (not _GetCells(domain, Y, yWeights, volumes, targets, err)) {
            std::cerr << "diagram of Y: " << err << std::endl;
            return 1;
        }
    }

    if (sitesAtTargets) {
        q = targets;
    }

    // The H solve from the sites q, or the H_alpha solve with "alpha a": the
    // stiffness a if a > 0, else -a times the practical bound
    // (eq:regularized-practical) at the centroids c0 of F at q, or with
    // "guaranteed 1" -a times the bound eq:regularized-distinct.
    double alpha = records["alpha"].empty() ? 0.0 : records["alpha"][0];
    if (alpha < 0.0) {
        PP::DualSolver::Parameters fParams;
        fParams._volumeTolerance = 1e-13;
        fParams._maxIterations = 200;
        std::vector<double> w0, V0, c0;
        const PP::DualSolver::Result r0 = PP::VolumeSolver(fParams).Solve(domain, q, volumes, w0);
        if ((not r0.HasConverged() and not (r0._volumeResidual <= 1e-10)) or
            not _GetCells(domain, q, w0, V0, c0, err)) {
            std::cerr << "volume solve at q: " << PP::DualSolver::GetStatusString(r0._status) << " " << err
                      << std::endl;
            return 1;
        }
        // and the guaranteed bound (eq:regularized-distinct), with the residual
        // resid = sum_i Vbar_i |c0_i - cbar_i|^2; "guaranteed 1" scales that one
        double resid = 0.0;
        for (int i = 0; i < n; ++i) {
            resid += volumes[size_t(i)] * (std::pow(c0[size_t(2 * i)] - targets[size_t(2 * i)], 2) +
                                           std::pow(c0[size_t(2 * i + 1)] - targets[size_t(2 * i + 1)], 2));
        }
        double guaranteed = HUGE_VAL;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double d = std::hypot(q[size_t(2 * i)] - q[size_t(2 * j)], q[size_t(2 * i + 1)] - q[size_t(2 * j + 1)]);
                guaranteed = std::min(
                    guaranteed, d / std::sqrt(resid * (1.0 / volumes[size_t(i)] + 1.0 / volumes[size_t(j)])));
            }
        }
        double bound = HUGE_VAL;
        for (int i = 0; i < n; ++i) {
            const double ri = std::hypot(c0[size_t(2 * i)] - targets[size_t(2 * i)],
                                         c0[size_t(2 * i + 1)] - targets[size_t(2 * i + 1)]);
            for (int j = i + 1; j < n; ++j) {
                const double rj = std::hypot(c0[size_t(2 * j)] - targets[size_t(2 * j)],
                                             c0[size_t(2 * j + 1)] - targets[size_t(2 * j + 1)]);
                const double d = std::hypot(q[size_t(2 * i)] - q[size_t(2 * j)], q[size_t(2 * i + 1)] - q[size_t(2 * j + 1)]);
                if (ri + rj > 0.0) {
                    bound = std::min(bound, d / (ri + rj));
                }
            }
        }
        const bool useGuaranteed = not records["guaranteed"].empty() and records["guaranteed"][0] != 0.0;
        alpha = -alpha * (useGuaranteed ? guaranteed : bound);
        std::cout << "guaranteed bound " << guaranteed << ", practical bound " << bound << ", alpha " << alpha << "; ";
    }
    PP::DualSolver::Parameters params;
    params._volumeTolerance = 1e-10;
    params._shiftTolerance = 1e-10;
    if (not records["maxIterations"].empty()) {
        params._maxIterations = int(records["maxIterations"][0]);
    }
    if (not records["log"].empty() and records["log"][0] != 0.0) {
        // Per iterate: the objective, the residuals, the cell volumes relative
        // to the mean, and the closest generator pair relative to the cell size.
        params._onIterate = [n, area](const PP::JointDual&, const PP::DualSolver::Vector&,
                                      const PP::DualSolver::Result& r) {
            const PP::ShiftedDiagram& D = r._diagram;
            double vmin = HUGE_VAL, vmax = 0.0, dmin = HUGE_VAL;
            for (int i = 0; i < D.GetNumSites(); ++i) {
                vmin = std::min(vmin, D.GetVolume(i));
                vmax = std::max(vmax, D.GetVolume(i));
            }
            for (size_t f = 0; f < D.GetNumFacets(); ++f) {
                dmin = std::min(dmin, D.GetGeneratorDistance(f));
            }
            const double mean = area / double(n);
            std::cerr << "  it " << r._iterations << "  value " << r._value << "  shift res " << r._shiftResidual
                      << "  vol res " << r._volumeResidual << "  V/mean [" << vmin / mean << ", " << vmax / mean
                      << "]  empty " << D.GetNumEmptyCells() << "  min |p_i - p_j| / h " << dmin / std::sqrt(mean)
                      << std::endl;
        };
    }
    std::vector<double> shifts, weights;
    const auto solveStart = std::chrono::steady_clock::now();
    PP::DualSolver::Result res;
    // "centroidOnly m": G with fixed weights instead, w = 0 (m = 1) or the
    // weights of F at q with the prescribed volumes (m = 2); the volumes are
    // then free, and the realized ones are written.
    const int centroidOnly = records["centroidOnly"].empty() ? 0 : int(records["centroidOnly"][0]);
    if (centroidOnly == 1 or centroidOnly == 2) {
        std::vector<double> w(size_t(n), 0.0);
        if (centroidOnly == 2) {
            PP::DualSolver::Parameters fParams;
            fParams._volumeTolerance = 1e-13;
            fParams._maxIterations = 200;
            const PP::DualSolver::Result r0 = PP::VolumeSolver(fParams).Solve(domain, q, volumes, w);
            if (not r0.HasConverged() and not (r0._volumeResidual <= 1e-10)) {
                std::cerr << "volume solve at q: " << PP::DualSolver::GetStatusString(r0._status) << std::endl;
                return 1;
            }
        }
        res = PP::CentroidSolver(params).Solve(domain, q, targets, w, volumes, shifts);
        weights = w;
    } else if (alpha > 0.0) {
        // from the diagram of the sites with the prescribed volumes (F at q)
        const PP::DualSolver::Result start = PP::MakeJointStart(domain, q, targets, volumes, 0.0, shifts, weights);
        if (not start.HasConverged()) {
            std::cerr << "start: " << PP::DualSolver::GetStatusString(start._status) << " " << start._error
                      << std::endl;
            return 1;
        }
        res = PP::SoftenedSolver(params).Solve(domain, q, targets, volumes, alpha, shifts, weights);
    } else {
        const PP::DualSolver::Result start = PP::MakeJointStart(domain, q, targets, volumes, 0.5, shifts, weights);
        if (not start.HasConverged()) {
            std::cerr << "start: " << PP::DualSolver::GetStatusString(start._status) << " " << start._error
                      << std::endl;
            return 1;
        }
        res = PP::JointSolver(params).Solve(domain, q, targets, volumes, shifts, weights);
    }
    const double solveSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - solveStart).count();
    if (not res.HasConverged()) {
        std::cerr << (centroidOnly ? "G" : alpha > 0.0 ? "H_alpha" : "H") << " solve: " << PP::DualSolver::GetStatusString(res._status)
                  << " (" << res._iterations << " iterations) " << res._error << std::endl;
        return 1;
    }

    // Generators and shifted weights (eq:shifted-site), and a check of the
    // cells they give against the targets.
    std::vector<double> p(2 * static_cast<size_t>(n)), lambda(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        double l = weights[size_t(i)];
        for (int c = 0; c < 2; ++c) {
            const size_t k = size_t(2 * i + c);
            p[k] = q[k] + shifts[k];
            l += shifts[k] * shifts[k] + 2.0 * shifts[k] * (q[k] - targets[k]);
        }
        lambda[size_t(i)] = l;
    }
    // H_alpha has a unique maximizer, so its generators are kept as they are.
    const bool closest = alpha <= 0.0 and not centroidOnly and not records["closest"].empty() and records["closest"][0] != 0.0;
    double s = 1.0, t[2] = {0.0, 0.0};
    if (closest) {
        // Normal equations of min sum_i Vbar_i |s p_i + t - q_i|^2 in (s, t).
        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
        Eigen::Vector3d b = Eigen::Vector3d::Zero();
        for (int i = 0; i < n; ++i) {
            const double vi = volumes[size_t(i)];
            for (int c = 0; c < 2; ++c) {
                const size_t k = size_t(2 * i + c);
                const Eigen::Vector3d row(p[k], c == 0 ? 1.0 : 0.0, c == 1 ? 1.0 : 0.0);
                A += vi * row * row.transpose();
                b += vi * row * q[k];
            }
        }
        const Eigen::Vector3d x = A.ldlt().solve(b);
        s = x[0];
        t[0] = x[1];
        t[1] = x[2];
        if (not (s > 0.0)) {
            std::cerr << "the closest homothety has s = " << s << std::endl;
            return 1;
        }
        for (int i = 0; i < n; ++i) {
            double pp = 0.0, pq = 0.0;
            for (int c = 0; c < 2; ++c) {
                const size_t k = size_t(2 * i + c);
                pp += p[k] * p[k];
                p[k] = s * p[k] + t[c];
                pq += p[k] * p[k];
                shifts[k] = p[k] - q[k];
            }
            lambda[size_t(i)] = pq - s * (pp - lambda[size_t(i)]);
        }
    }
    std::vector<double> V, C;
    if (not _GetCells(domain, p, lambda, V, C, err)) {
        std::cerr << "diagram of p: " << err << std::endl;
        return 1;
    }
    double volumeError = 0.0, centroidError = 0.0;
    const double size = std::sqrt(area / double(n));
    for (int i = 0; i < n; ++i) {
        volumeError = std::max(volumeError, std::abs(V[size_t(i)] - volumes[size_t(i)]) / volumes[size_t(i)]);
        centroidError = std::max(centroidError, std::hypot(C[size_t(2 * i)] - targets[size_t(2 * i)],
                                                           C[size_t(2 * i + 1)] - targets[size_t(2 * i + 1)]) /
                                                    size);
    }

    if (centroidOnly) {
        // the volumes are an output of G
        double vmin = HUGE_VAL, vmax = 0.0;
        for (double v : V) {
            vmin = std::min(vmin, v);
            vmax = std::max(vmax, v);
        }
        std::cout << "volumes / (area / n): " << vmin * n / area << " to " << vmax * n / area << "; ";
        volumes = V;
        volumeError = 0.0;
    }
    std::ofstream out(argv[2]);
    out.precision(17);
    out << "n " << n << "\n";
    _WriteRecord(out, "volumes", volumes);
    _WriteRecord(out, "targets", targets);
    _WriteRecord(out, "sites", q);
    _WriteRecord(out, "generators", p);
    _WriteRecord(out, "shifts", shifts);
    _WriteRecord(out, "weights", lambda);
    std::cout << (centroidOnly ? "G: " : alpha > 0.0 ? "H_alpha: " : "H: ") << res._iterations << " iterations, max relative volume error " << volumeError
              << ", max centroid error / cell size " << centroidError << ", " << solveSeconds
              << " s (start and solve)";
    if (closest) {
        std::cout << "; closest homothety s = " << s << ", t = (" << t[0] << ", " << t[1] << ")";
    }
    std::cout << std::endl;
    return 0;
}
