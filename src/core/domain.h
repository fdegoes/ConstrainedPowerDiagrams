#ifndef PP_CORE_DOMAIN_H
#define PP_CORE_DOMAIN_H

#include <memory>
#include <string>
#include <vector>

namespace GEO
{
class Mesh;
}

namespace PP
{

class CellBackend;

// Simplicial mesh of a clipping domain (triangles in 2D, tetrahedra
// in 3D), ready to be intersected with power cells.
//
// Building it once and reusing it avoids redoing the mesh setup (adjacency,
// orientation) at every evaluation of a power diagram, e.g. inside Newton
// iterations.
//
// Requirements on the mesh (not checked): conforming (no T-junctions, no
// overlaps), welded, and no degenerate elements. Element orientation does
// not matter.
class Domain
{
public:
    Domain() = default;
    ~Domain() = default;

    Domain(Domain and ) noexcept = default;
    Domain& operator=(Domain and ) noexcept = default;

    Domain(const Domain&) = delete;
    Domain& operator=(const Domain&) = delete;

    // A domain of dimension dim (2: triangles, 3: tetrahedra).
    // pts    dim doubles per vertex
    // elems  dim + 1 vertex indices per element
    static bool Make(
        int dim,
        const std::vector<double>& pts,
        const std::vector<int>& elems,
        Domain& out,
        std::string& err);

    // Make(2, ...).
    // triPts 2 doubles per vertex
    // tris   3 vertex indices per triangle
    static bool Make2D(
        const std::vector<double>& triPts,
        const std::vector<int>& tris,
        Domain& out,
        std::string& err);

    // Make(3, ...).
    // tetPts 3 doubles per vertex
    // tets   4 vertex indices per tetrahedron
    static bool Make3D(
        const std::vector<double>& tetPts,
        const std::vector<int>& tets,
        Domain& out,
        std::string& err);

    bool IsValid() const
    { return _mesh != nullptr; }

    int GetDimension() const
    { return _dim; }

    double GetVolume() const
    { return _volume; }

    const std::vector<double>& GetCentroid() const
    { return _centroid; }

    const std::vector<double>& GetPoints() const
    { return _points; }

    // Simplicial mesh with adjacency, embedded in the lifted space of the
    // weighted sites (dimension d + 1, extra coordinate 0). For the core
    // library only (callers need the geogram headers).
    GEO::Mesh& GetMesh() const
    { return *_mesh; }

    // The dimension-specific geometry (cellBackend.h). For the core
    // library only.
    const CellBackend& GetBackend() const
    { return *_backend; }

    // The mesh of the parallel integration of power cells (ComputePowerCells):
    // the same domain, in Hilbert order (and in 2D refined by midpoint
    // subdivision until each part has a few triangles), split into parts of
    // contiguous elements, part p being [ranges[p], ranges[p + 1]). For the
    // core library only.
    GEO::Mesh& GetWalkMesh() const
    { return *_walkMesh; }
    const std::vector<unsigned>& GetWalkRanges() const
    { return _walkRanges; }

private:
    // Defined in domain.cpp, where GEO::Mesh is complete, so that this
    // header does not need the geogram headers.
    struct MeshDeleter
    {
        void operator()(GEO::Mesh* mesh) const;
    };

    // Builds _walkMesh and _walkRanges from _mesh.
    void _MakeWalkMesh();

    std::unique_ptr<GEO::Mesh, MeshDeleter> _mesh;
    std::unique_ptr<GEO::Mesh, MeshDeleter> _walkMesh;
    std::vector<unsigned>                   _walkRanges;
    const CellBackend*                      _backend = nullptr;
    int                                     _dim = 0;
    double                                  _volume = 0.0;
    std::vector<double>                     _centroid;
    std::vector<double>                     _points;
};

} // namespace PP

#endif // PP_CORE_DOMAIN_H
