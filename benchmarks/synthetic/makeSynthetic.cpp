// Realizable data of the synthetic benchmark (benchmarks/synthetic): from
// generators Y in the unit square (2D) or cube (3D) and a volume distribution, the power
// diagram of Y with those volumes (VolumeSolver, F, at Y). Its volumes and
// centroids are the targets (Vbar, cbar), and Y with the weights of the
// diagram is the known solution.
//
//   makeSynthetic <generators.txt> <volumes> <out.txt>
//
// generators.txt  "n" or "n d" (d = 2 or 3, default 2), then d n values (Y)
// volumes         uniform      Vbar_i = 1 / n
//                 voronoi      the volumes of the Voronoi cells of Y
//                 <file.txt>   n, then the n volumes (normalized to sum 1)
// out.txt         "n", "volumes", "targets", "weights" records (the format of
//                 benchBourne's case files, plus the weights of the solution)

#include "core/domain.h"
#include "core/dualSolver.h"
#include "core/jointDual.h"
#include "core/shiftedDiagram.h"

#include "../common/boxDomain.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{

// Volumes of the cells of the generators with weights (as a diagram of H with
// zero shifts).
bool
_GetVolumes(
    const PP::Domain& domain,
    const std::vector<double>& generators,
    const std::vector<double>& weights,
    std::vector<double>& volumes,
    std::vector<double>& centroids,
    std::string& err)
{
    const int d = domain.GetDimension();
    const size_t n = generators.size() / size_t(d);
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
    volumes.resize(n);
    centroids.resize(size_t(d) * n);
    for (size_t i = 0; i < n; ++i) {
        volumes[i] = diagram.GetVolume(int(i));
        for (int k = 0; k < d; ++k) {
            centroids[size_t(d) * i + size_t(k)] = diagram.GetCentroid(int(i))[k];
        }
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
    if (argc != 4) {
        std::cerr << "usage: makeSynthetic <generators.txt> <uniform|voronoi|volumes.txt> <out.txt>" << std::endl;
        return 2;
    }
    std::ifstream in(argv[1]);
    std::string header;
    std::getline(in, header);
    std::istringstream hs(header);
    int n = 0, d = 2;
    hs >> n;
    if (not (hs >> d)) {
        d = 2;
    }
    if (n <= 0 or (d != 2 and d != 3)) {
        std::cerr << "cannot read the header of " << argv[1] << std::endl;
        return 2;
    }
    std::vector<double> Y(size_t(d) * size_t(n));
    for (double& y : Y) {
        in >> y;
    }
    if (not in) {
        std::cerr << "cannot read " << argv[1] << std::endl;
        return 2;
    }

    std::string err;
    PP::Domain domain;
    std::vector<double> box(size_t(2 * d), 0.0);
    std::fill(box.begin() + d, box.end(), 1.0);  // the unit square or cube
    const bool ok = PPBench::MakeBoxDomain(d, box, domain, err);
    if (not ok) {
        std::cerr << err << std::endl;
        return 2;
    }

    const std::string mode = argv[2];
    std::vector<double> volumes(size_t(n), 1.0 / double(n)), centroids;
    if (mode == "voronoi") {
        if (not _GetVolumes(domain, Y, {}, volumes, centroids, err)) {
            std::cerr << err << std::endl;
            return 2;
        }
    } else if (mode != "uniform") {
        std::ifstream vin(mode);
        int m = 0;
        vin >> m;
        double sum = 0.0;
        for (double& v : volumes) {
            vin >> v;
            sum += v;
        }
        if (not vin or m != n or not (sum > 0.0)) {
            std::cerr << "cannot read " << n << " volumes from " << mode << std::endl;
            return 2;
        }
        for (double& v : volumes) {
            v /= sum;
        }
    }

    // The diagram of Y with these volumes.
    std::vector<double> weights;
    PP::DualSolver::Parameters params;
    params._volumeTolerance = 1e-13;
    params._maxIterations = 200;
    const PP::DualSolver::Result res = PP::VolumeSolver(params).Solve(domain, Y, volumes, weights);
    if (not res.HasConverged() and not (res._volumeResidual <= 1e-10)) {
        std::cerr << "volume solve at Y: " << PP::DualSolver::GetStatusString(res._status) << " (residual "
                  << res._volumeResidual << ") " << res._error << std::endl;
        return 1;
    }
    std::vector<double> realized;
    if (not _GetVolumes(domain, Y, weights, realized, centroids, err)) {
        std::cerr << err << std::endl;
        return 2;
    }

    std::ofstream out(argv[3]);
    out.precision(17);
    out << "n " << n << "\n";
    if (d != 2) {
        out << "dim " << d << "\n";
    }
    // The realized volumes: equal to the prescribed ones up to the solve's
    // tolerance, and exactly those of the diagram whose centroids are cbar.
    _WriteRecord(out, "volumes", realized);
    _WriteRecord(out, "targets", centroids);
    _WriteRecord(out, "weights", weights);
    std::cout << "volume residual " << res._volumeResidual << " after " << res._iterations << " iterations"
              << std::endl;
    return 0;
}
