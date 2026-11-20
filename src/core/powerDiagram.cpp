#include "core/powerDiagram.h"
#include "core/domain.h"
#include "core/geogramContext.h"
#include "core/liftSites.h"
#include "core/powerCells2D.h"
#include "core/powerCells3D.h"

namespace PP
{


bool
ComputePowerDiagram(
    const Domain& domain,
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    double boundaryAngle,
    PolyMesh& out,
    std::string& err)
{
    out = PolyMesh();
    if (not domain.IsValid()) {
        err = "Invalid clipping domain.";
        return false;
    }
    out._dim = domain.GetDimension();
    if (not CheckSites(sites, weights, out._dim, err)) {
        return false;
    }

    // Cells read from the regular triangulation (powerCells2D.h,
    // powerCells3D.h).
    GeogramLock lock;
    const PowerTriangulation pd = MakePowerTriangulation(domain, sites, weights);
    if (out._dim == 2) {
        if (not ComputePowerPolygons2D(domain, int(weights.size()), pd, out)) {
            err = "The power cells could not be built (unbounded cell).";
            return false;
        }
    } else {
        ComputePowerPolygons3D(domain, int(weights.size()), pd, boundaryAngle, out);
    }
    return true;
}

bool
ComputePowerDiagram2D(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& triPts,
    const std::vector<int>& tris,
    PolyMesh& out,
    std::string& err)
{
    out = PolyMesh();
    out._dim = 2;
    Domain domain;
    if (not CheckSites(sites, weights, 2, err) or not Domain::Make2D(triPts, tris, domain, err)) {
        return false;
    }
    return ComputePowerDiagram(domain, sites, weights, 0.0, out, err);
}

bool
ComputePowerDiagram3D(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& tetPts,
    const std::vector<int>& tets,
    double boundaryAngle,
    PolyMesh& out,
    std::string& err)
{
    out = PolyMesh();
    out._dim = 3;
    Domain domain;
    if (not CheckSites(sites, weights, 3, err) or not Domain::Make3D(tetPts, tets, domain, err)) {
        return false;
    }
    return ComputePowerDiagram(domain, sites, weights, boundaryAngle, out, err);
}

} // namespace PP
