#include "core/powerDiagram.h"
#include "core/cellBackend.h"
#include "core/domain.h"
#include "core/geogramContext.h"
#include "core/liftSites.h"

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

    GeogramLock lock;
    ClippedPowerDiagram pd = MakeClippedPowerDiagram(domain, sites, weights);
    if (pd._hiddenSeed) {
        err = "Some power cells are empty and the diagram cannot be clipped.";
        return false;
    }
    domain.GetBackend().ExtractPolygons(*pd._rvd, GEO::index_t(weights.size()), boundaryAngle, out);
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
