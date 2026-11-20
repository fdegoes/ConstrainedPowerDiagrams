#ifndef PP_CORE_CELLBACKEND_H
#define PP_CORE_CELLBACKEND_H

// The dimension-specific parts of the domain of clipped power diagrams
// (internal to the core library): the measure of the domain mesh, its
// geogram storage and adjacency, and the regular triangulation to use. One
// backend per supported dimension (2: triangle domains, 3: tetrahedral
// domains); the cells themselves are clipped by powerCells2D.h and
// powerCells3D.h. Callers must hold a GeogramLock where geogram is used.

#include <geogram/mesh/mesh.h>

#include <memory>
#include <vector>

namespace PP
{

class CellBackend
{
public:
    virtual ~CellBackend() = default;

    virtual int GetDimension() const = 0;

    // Number of vertices of a domain element (3: triangle, 4: tetrahedron).
    virtual int GetElementSize() const = 0;

    // geogram's Delaunay implementation of the regular triangulation of the
    // lifted sites (dimension d + 1).
    virtual const char* GetDelaunayName() const = 0;

    // Total measure (area or volume) and centroid of the elements of a
    // simplicial mesh (d coordinates per point).
    virtual void ComputeMeasureAndCentroid(
        const std::vector<double>& pts,
        const std::vector<int>& elems,
        double& measure,
        std::vector<double>& centroid) const = 0;

    // The domain mesh, with element adjacency. Its points have d + 1
    // coordinates (the last one 0), the dimension of the lifted sites.
    virtual std::unique_ptr<GEO::Mesh> MakeMesh(
        const std::vector<double>& pts,
        const std::vector<int>& elems) const = 0;

    // Elements of a domain mesh (facets in 2D, cells in 3D) and their
    // adjacency.
    virtual GEO::index_t GetNumElements(const GEO::Mesh& M) const = 0;
    virtual GEO::index_t GetElementVertex(const GEO::Mesh& M, GEO::index_t e, GEO::index_t lv) const = 0;
    virtual GEO::index_t GetNumAdjacent(const GEO::Mesh& M, GEO::index_t e) const = 0;
    virtual GEO::index_t GetAdjacent(const GEO::Mesh& M, GEO::index_t e, GEO::index_t k) const = 0;
};

// The backend of dimension dim, or nullptr if dim is not supported (2 or 3).
const CellBackend* GetCellBackend(int dim);

} // namespace PP

#endif // PP_CORE_CELLBACKEND_H
