#include "core/cellBackend.h"
#include "core/cellMoments.h"
#include "core/geogramMesh.h"
#include "core/geometry.h"

#include <cmath>

namespace PP
{

namespace
{

// ---------------------------------------------------------------- measures

// Total area and centroid of a 2D triangle mesh.
void
_ComputeAreaAndCentroid(
    const std::vector<double>& pts,
    const std::vector<int>& tris,
    double& area,
    std::vector<double>& centroid)
{
    area = 0.0;
    centroid.assign(2, 0.0);
    for (size_t t = 0; t < tris.size(); t += 3) {
        const double* a = &pts[2 * tris[t]];
        const double* b = &pts[2 * tris[t + 1]];
        const double* c = &pts[2 * tris[t + 2]];
        const double A = 0.5 * std::fabs(ComputeDoubleSignedArea(a, b, c));
        area += A;
        for (int k = 0; k < 2; ++k) {
            centroid[k] += A * (a[k] + b[k] + c[k]) / 3.0;
        }
    }
    if (area > 0.0) {
        for (double& x : centroid) {
            x /= area;
        }
    }
}

// Total volume and centroid of a tet mesh.
void
_ComputeVolumeAndCentroid(
    const std::vector<double>& pts,
    const std::vector<int>& tets,
    double& volume,
    std::vector<double>& centroid)
{
    volume = 0.0;
    centroid.assign(3, 0.0);
    for (size_t t = 0; t < tets.size(); t += 4) {
        const double* p[4] = {&pts[3 * tets[t]], &pts[3 * tets[t + 1]],
                              &pts[3 * tets[t + 2]], &pts[3 * tets[t + 3]]};
        const double V = std::fabs(ComputeSixSignedVolume(p[0], p[1], p[2], p[3])) / 6.0;
        volume += V;
        for (int k = 0; k < 3; ++k) {
            centroid[k] += V * 0.25 * (p[0][k] + p[1][k] + p[2][k] + p[3][k]);
        }
    }
    if (volume > 0.0) {
        for (double& x : centroid) {
            x /= volume;
        }
    }
}

// ---------------------------------------------------------------- backends

// 2D: triangle domains, embedded in the plane z = 0 of the lifted space.
class CellBackend2D final : public CellBackend
{
public:
    int GetDimension() const override { return 2; }
    int GetElementSize() const override { return 3; }
    const char* GetDelaunayName() const override { return "BPOW2d"; }

    void ComputeMeasureAndCentroid(
        const std::vector<double>& pts,
        const std::vector<int>& elems,
        double& measure,
        std::vector<double>& centroid) const override
    {
        _ComputeAreaAndCentroid(pts, elems, measure, centroid);
    }

    std::unique_ptr<GEO::Mesh> MakeMesh(const std::vector<double>& pts, const std::vector<int>& elems) const override
    {
        // Lifted space (x, y, 0), with triangle-triangle adjacency.
        std::unique_ptr<GEO::Mesh> mesh(new GEO::Mesh);
        SetMeshVertices(*mesh, pts, 2, 3);
        AddMeshTriangles(*mesh, elems);
        mesh->facets.connect();
        return mesh;
    }

    GEO::index_t GetNumElements(const GEO::Mesh& M) const override { return M.facets.nb(); }
    GEO::index_t GetElementVertex(const GEO::Mesh& M, GEO::index_t e, GEO::index_t lv) const override
    { return M.facets.vertex(e, lv); }
    GEO::index_t GetNumAdjacent(const GEO::Mesh& M, GEO::index_t e) const override { return M.facets.nb_vertices(e); }
    GEO::index_t GetAdjacent(const GEO::Mesh& M, GEO::index_t e, GEO::index_t k) const override
    { return M.facets.adjacent(e, k); }
};

// 3D: tetrahedral domains, lifted to (x, y, z, 0).
class CellBackend3D final : public CellBackend
{
public:
    int GetDimension() const override { return 3; }
    int GetElementSize() const override { return 4; }
    const char* GetDelaunayName() const override { return "BPOW"; }

    void ComputeMeasureAndCentroid(
        const std::vector<double>& pts,
        const std::vector<int>& elems,
        double& measure,
        std::vector<double>& centroid) const override
    {
        _ComputeVolumeAndCentroid(pts, elems, measure, centroid);
    }

    std::unique_ptr<GEO::Mesh> MakeMesh(const std::vector<double>& pts, const std::vector<int>& elems) const override
    {
        // Positive tets, with tet-tet adjacency, then lifted to (x, y, z, 0).
        std::unique_ptr<GEO::Mesh> mesh(new GEO::Mesh);
        SetMeshVertices(*mesh, pts, 3, 3);
        AddMeshTets(*mesh, pts, elems);
        mesh->cells.connect();
        mesh->vertices.set_dimension(4);
        return mesh;
    }

    GEO::index_t GetNumElements(const GEO::Mesh& M) const override { return M.cells.nb(); }
    GEO::index_t GetElementVertex(const GEO::Mesh& M, GEO::index_t e, GEO::index_t lv) const override
    { return M.cells.vertex(e, lv); }
    GEO::index_t GetNumAdjacent(const GEO::Mesh& M, GEO::index_t e) const override { return M.cells.nb_facets(e); }
    GEO::index_t GetAdjacent(const GEO::Mesh& M, GEO::index_t e, GEO::index_t k) const override
    { return M.cells.adjacent(e, k); }
};

static_assert(kMaxDimension >= 3, "the moment arrays must hold the 3D backend");

} // namespace

const CellBackend*
GetCellBackend(int dim)
{
    static const CellBackend2D backend2D;
    static const CellBackend3D backend3D;
    switch (dim) {
    case 2:
        return &backend2D;
    case 3:
        return &backend3D;
    default:
        return nullptr;
    }
}

} // namespace PP
