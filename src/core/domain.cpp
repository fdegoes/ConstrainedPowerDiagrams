#include "core/domain.h"
#include "core/cellBackend.h"
#include "core/geogramContext.h"
#include "core/geogramMesh.h"
#include "core/geometry.h"

#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <utility>

namespace PP
{

void
Domain::MeshDeleter::operator()(GEO::Mesh* mesh) const
{
    delete mesh;
}

bool
Domain::Make(
    int dim,
    const std::vector<double>& pts,
    const std::vector<int>& elems,
    Domain& out,
    std::string& err)
{
    out = Domain();

    const CellBackend* backend = GetCellBackend(dim);
    if (backend == nullptr) {
        err = "The clipping domain must be 2D or 3D.";
        return false;
    }
    const int perElem = backend->GetElementSize();
    if (not CheckElementIndices(elems, pts.size() / size_t(dim), perElem, perElem == 3 ? "triangle" : "tetrahedron",
                                err)) {
        return false;
    }

    if (not AreAllFinite(pts)) {
        err = "The clipping domain has non-finite coordinates.";
        return false;
    }

    double measure = 0.0;
    std::vector<double> centroid;
    backend->ComputeMeasureAndCentroid(pts, elems, measure, centroid);

    if (not (measure > 0.0)) {
        err = dim == 2 ? "The clipping domain has zero area." : "The clipping domain has zero volume.";
        return false;
    }

    GeogramLock lock;
    out._mesh.reset(backend->MakeMesh(pts, elems).release());
    out._backend = backend;
    out._dim = dim;
    out._volume = measure;
    out._centroid = std::move(centroid);
    out._points = pts;
    return true;
}

bool
Domain::Make2D(
    const std::vector<double>& triPts,
    const std::vector<int>& tris,
    Domain& out,
    std::string& err)
{
    return Make(2, triPts, tris, out, err);
}

bool
Domain::Make3D(
    const std::vector<double>& tetPts,
    const std::vector<int>& tets,
    Domain& out,
    std::string& err)
{
    return Make(3, tetPts, tets, out, err);
}

} // namespace PP
