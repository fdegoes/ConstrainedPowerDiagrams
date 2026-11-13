#ifndef PP_CORE_CELLBACKEND_H
#define PP_CORE_CELLBACKEND_H

// The dimension-specific parts of clipped power diagrams (internal to the
// core library): how the domain mesh is stored and walked, how the
// fragments of the restricted Voronoi diagram are integrated, and how the
// cells are extracted as polygons. One backend per supported dimension (2:
// triangle domains, 3: tetrahedral domains); everything above it (Domain,
// ComputePowerCells, ShiftedDiagram, the solvers) is written once for both.
// Callers must hold a GeogramLock where geogram is used.

#include <geogram/mesh/mesh.h>
#include <geogram/voronoi/RVD.h>

#include <memory>
#include <vector>

namespace PP
{

class Accumulator;
struct PolyMesh;

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

    // Whether the restricted Voronoi diagram clips volumes (tetrahedra) or
    // surfaces (triangles).
    virtual bool IsVolumetric() const = 0;

    // Total measure (area or volume) and centroid of the elements of a
    // simplicial mesh (d coordinates per point).
    virtual void ComputeMeasureAndCentroid(
        const std::vector<double>& pts,
        const std::vector<int>& elems,
        double& measure,
        std::vector<double>& centroid) const = 0;

    // The domain mesh in the lifted space of the weighted sites (dimension
    // d + 1, extra coordinate 0), with element adjacency.
    virtual std::unique_ptr<GEO::Mesh> MakeMesh(
        const std::vector<double>& pts,
        const std::vector<int>& elems) const = 0;

    // The mesh of the parallel walk: a copy of mesh, refined (2D: midpoint
    // subdivision) until it has at least minElements elements, in Hilbert
    // order and with adjacency.
    virtual std::unique_ptr<GEO::Mesh> MakeWalkMesh(const GEO::Mesh& mesh, size_t minElements) const = 0;

    // Elements of a domain mesh (facets in 2D, cells in 3D) and their
    // adjacency.
    virtual GEO::index_t GetNumElements(const GEO::Mesh& M) const = 0;
    virtual GEO::index_t GetElementVertex(const GEO::Mesh& M, GEO::index_t e, GEO::index_t lv) const = 0;
    virtual GEO::index_t GetNumAdjacent(const GEO::Mesh& M, GEO::index_t e) const = 0;
    virtual GEO::index_t GetAdjacent(const GEO::Mesh& M, GEO::index_t e, GEO::index_t k) const = 0;

    // Restricts rvd to the elements [first, last) of its mesh.
    virtual void SetElementRange(GEO::RestrictedVoronoiDiagram& rvd, GEO::index_t first, GEO::index_t last) const = 0;

    // Integrates the cells and facets of rvd into acc.
    virtual void Walk(GEO::RestrictedVoronoiDiagram& rvd, Accumulator& acc) const = 0;

    // Boundary polygons of the cells of rvd (see PolyMesh).
    virtual void ExtractPolygons(
        GEO::RestrictedVoronoiDiagram& rvd,
        GEO::index_t nbSites,
        double boundaryAngle,
        PolyMesh& out) const = 0;
};

// The backend of dimension dim, or nullptr if dim is not supported (2 or 3).
const CellBackend* GetCellBackend(int dim);

} // namespace PP

#endif // PP_CORE_CELLBACKEND_H
