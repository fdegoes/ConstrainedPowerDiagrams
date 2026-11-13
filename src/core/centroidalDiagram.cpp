#include "core/centroidalDiagram.h"
#include "core/domain.h"
#include "core/geometry.h"

#include <algorithm>
#include <cmath>

namespace PP
{

namespace
{

// max_i |c_i - q_i|.
double
_ComputeMaxCentroidDistance(const ShiftedDiagram& diagram, const std::vector<double>& sites)
{
    const int dim = diagram.GetDimension();
    double d2 = 0.0;
    for (int i = 0; i < diagram.GetNumSites(); ++i) {
        d2 = std::max(d2, ComputeDistance2(diagram.GetCentroid(i), &sites[dim * i], dim));
    }
    return std::sqrt(d2);
}

} // namespace

CentroidalDiagram::Result
CentroidalDiagram::Solve(
    const Domain& domain,
    const std::vector<double>& volumes,
    std::vector<double>& sites,
    std::vector<double>& weights) const
{
    Result res;
    if (not domain.IsValid()) {
        res._error = "Invalid clipping domain.";
        return res;
    }
    BoundingBox box;
    box.Add(domain.GetPoints(), domain.GetDimension());
    const double diam = box.GetDiagonal(domain.GetDimension());

    VolumeSolver solver(_params._volumeSolver);
    for (;;) {
        const DualSolver::Result vr = solver.Solve(domain, sites, volumes, weights);
        if (not vr.HasConverged()) {
            res._error = std::string("Volume solve failed: ") + DualSolver::GetStatusString(vr._status) +
                         (vr._error.empty() ? "." : ": " + vr._error);
            return res;
        }
        res._residual = _ComputeMaxCentroidDistance(vr._diagram, sites) / diam;
        if (res._residual <= _params._tolerance) {
            res._converged = true;
            return res;
        }
        if (res._iterations >= _params._maxIterations) {
            return res;
        }
        sites = vr._diagram.GetCells()._centroids;
        ++res._iterations;
    }
}

ParticleState
CentroidalDiagram::MakeInitialState(const std::vector<double>& sites, const std::vector<double>& weights)
{
    ParticleState s;
    s._positions = sites;
    s._velocities.assign(sites.size(), 0.0);
    s._shifts.assign(sites.size(), 0.0);
    s._weights = weights;
    return s;
}

} // namespace PP
