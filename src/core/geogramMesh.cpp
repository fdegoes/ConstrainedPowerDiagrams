#include "core/geogramMesh.h"
#include "core/geometry.h"

#include <utility>

namespace PP
{

bool
CheckElementIndices(
    const std::vector<int>& elems,
    size_t nbPts,
    int perElem,
    const char* what,
    std::string& err)
{
    if (elems.empty() or elems.size() % perElem != 0) {
        err = std::string("The mesh needs at least one ") + what + ".";
        return false;
    }

    for (int v : elems) 
    {
        if (v < 0 or size_t(v) >= nbPts) 
        {
            err = std::string(what) + " vertex index out of range.";
            return false;
        }
    }
    return true;
}

void
SetMeshVertices(
    GEO::Mesh& M, 
    const std::vector<double>& pts, 
    int dim, 
    int meshDim)
{
    const GEO::index_t nbPts = GEO::index_t(pts.size() / dim);
    M.vertices.set_dimension(GEO::index_t(meshDim));
    M.vertices.create_vertices(nbPts);
    for (GEO::index_t v = 0; v < nbPts; ++v) 
    {
        double* p = M.vertices.point_ptr(v);
        for (int c = 0; c < meshDim; ++c) 
        {
            p[c] = (c < dim) ? pts[dim * v + c] : 0.0;
        }
    }
}

void
AddMeshTriangles(
    GEO::Mesh& M, 
    const std::vector<int>& tris)
{
    for (size_t t = 0; t < tris.size(); t += 3) 
    {
        M.facets.create_triangle(
            GEO::index_t(tris[t]), 
            GEO::index_t(tris[t + 1]), 
            GEO::index_t(tris[t + 2]));
    }
}

void
AddMeshTets(
    GEO::Mesh& M, 
    const std::vector<double>& pts, 
    const std::vector<int>& tets)
{
    for (size_t t = 0; t < tets.size(); t += 4) 
    {
        GEO::index_t v[4] = 
        {
            GEO::index_t(tets[t]), 
            GEO::index_t(tets[t + 1]),
            GEO::index_t(tets[t + 2]), 
            GEO::index_t(tets[t + 3])
        };
        
        const double sixVol = ComputeSixSignedVolume(
            &pts[3 * v[0]], 
            &pts[3 * v[1]], 
            &pts[3 * v[2]], 
            &pts[3 * v[3]]);
        
        if (sixVol < 0.0) {
            std::swap(v[0], v[1]);
        }
        
        M.cells.create_tet(v[0], v[1], v[2], v[3]);
    }
}

std::vector<double>
GetMeshVertices(const GEO::Mesh& M, int dim)
{
    std::vector<double> pts(size_t(dim) * M.vertices.nb());
    for (GEO::index_t v = 0; v < M.vertices.nb(); ++v) {
        for (int c = 0; c < dim; ++c) {
            pts[dim * v + c] = M.vertices.point_ptr(v)[c];
        }
    }
    return pts;
}

} // namespace PP
