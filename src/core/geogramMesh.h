#ifndef PP_CORE_GEOGRAMMESH_H
#define PP_CORE_GEOGRAMMESH_H

// Conversions between raw arrays and geogram meshes (internal to the core
// library). Callers must hold a GeogramLock.

#include <geogram/mesh/mesh.h>

#include <string>
#include <vector>

namespace PP
{

// Checks that elems has a positive multiple of perElem indices, all
// in [0, nbPts).
bool CheckElementIndices(
    const std::vector<int>& elems,
    size_t nbPts,
    int perElem,
    const char* what,
    std::string& err);

// Copies pts (dim coordinates per point) into the vertices of M,
// whose dimension becomes meshDim >= dim (extra coordinates are 0).
void SetMeshVertices(
    GEO::Mesh& M, 
    const std::vector<double>& pts, 
    int dim, 
    int meshDim);

// Adds the triangles of tris (3 indices each) to M.
void AddMeshTriangles(
    GEO::Mesh& M, 
    const std::vector<int>& tris);

// Adds the tets of tets (4 indices each) to M, reoriented to positive
// volume using pts (3 coordinates per point).
void AddMeshTets(
    GEO::Mesh& M, 
    const std::vector<double>& pts, 
    const std::vector<int>& tets);

// Copies the first dim coordinates of the vertices of M.
std::vector<double> GetMeshVertices(
    const GEO::Mesh& M, 
    int dim);

} // namespace PP

#endif // PP_CORE_GEOGRAMMESH_H
