#include "core/tetrahedralize.h"
#include "core/geogramContext.h"
#include "core/geogramMesh.h"
#include "core/geometry.h"

#include <geogram/mesh/mesh_tetrahedralize.h>

#include <utility>

namespace PP
{

namespace
{

GEO::MeshTetrahedralizeParameters
_MakeTetGenParameters(const TetrahedralizeOptions& options)
{
    GEO::MeshTetrahedralizeParameters params;
    params.preprocess = options._preprocess;
    params.refine = options._refine;
    params.refine_quality = options._quality;
    params.verbose = false;
    return params;
}

// Tets of M, reoriented to positive volume.
std::vector<int>
_GetPositiveTets(const GEO::Mesh& M, const std::vector<double>& pts)
{
    std::vector<int> tets;
    tets.reserve(4 * M.cells.nb());
    for (GEO::index_t t = 0; t < M.cells.nb(); ++t) 
    {
        int v[4];
        for (GEO::index_t lv = 0; lv < 4; ++lv) {
            v[lv] = int(M.cells.vertex(t, lv));
        }
        
        if (ComputeSixSignedVolume(
            &pts[3 * v[0]], 
            &pts[3 * v[1]], 
            &pts[3 * v[2]], 
            &pts[3 * v[3]]) < 0.0) {
            std::swap(v[0], v[1]);
        }

        tets.insert(tets.end(), v, v + 4);
    }
    return tets;
}

} // namespace

bool
Tetrahedralize(
    const std::vector<double>& surfPts,
    const std::vector<int>& surfTris,
    const TetrahedralizeOptions& options,
    std::vector<double>& tetPts,
    std::vector<int>& tets,
    std::string& err)
{
    tetPts.clear();
    tets.clear();
    
    if (surfTris.size() < 12) {
        err = "The surface needs at least 4 triangles.";
        return false;
    }

    if (not CheckElementIndices(surfTris, surfPts.size() / 3, 3, "triangle", err)) {
        return false;
    }
    
    if (not AreAllFinite(surfPts)) {
        err = "The surface has non-finite coordinates.";
        return false;
    }

    GeogramLock lock;
    GEO::Mesh M;
    SetMeshVertices(M, surfPts, 3, 3);
    AddMeshTriangles(M, surfTris);
    M.facets.connect();

    if (not GEO::mesh_tetrahedralize(M, _MakeTetGenParameters(options)) or M.cells.nb() == 0) {
        err = options._preprocess
            ? "TetGen failed, even after preprocessing."
            : "TetGen failed: the surface must be closed, manifold, welded "
              "and free of self-intersections (or enable Preprocess).";
        return false;
    }

    tetPts = GetMeshVertices(M, 3);
    tets = _GetPositiveTets(M, tetPts);
    return true;
}

} // namespace PP
