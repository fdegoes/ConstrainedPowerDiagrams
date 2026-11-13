#include "core/cellBackend.h"
#include "core/cellMoments.h"
#include "core/geogramMesh.h"
#include "core/geometry.h"
#include "core/powerDiagram.h"

#include <geogram/mesh/mesh_reorder.h>
#include <geogram/voronoi/RVD_callback.h>
#include <geogram/voronoi/generic_RVD_polygon.h>
#include <geogram/voronoi/generic_RVD_vertex.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <utility>

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

// Splits each triangle of tris into 4 (midpoint subdivision), sharing the
// edge midpoints so that the result stays welded.
void
_SubdivideTriangles(std::vector<double>& pts, std::vector<int>& tris)
{
    std::map<std::pair<int, int>, int> midpoints;
    auto midpoint = [&](int a, int b) {
        const std::pair<int, int> key(std::min(a, b), std::max(a, b));
        const auto it = midpoints.find(key);
        if (it != midpoints.end()) {
            return it->second;
        }
        const int m = int(pts.size() / 2);
        pts.push_back(0.5 * (pts[2 * a] + pts[2 * b]));
        pts.push_back(0.5 * (pts[2 * a + 1] + pts[2 * b + 1]));
        midpoints.emplace(key, m);
        return m;
    };
    std::vector<int> out;
    out.reserve(4 * tris.size());
    for (size_t t = 0; t < tris.size(); t += 3) {
        const int a = tris[t], b = tris[t + 1], c = tris[t + 2];
        const int ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
        out.insert(out.end(), {a, ab, ca, ab, b, bc, ca, bc, c, ab, bc, ca});
    }
    tris = std::move(out);
}

// ---------------------------------------------------------------- integration

// Index of a bisector shared by the symbolic vertices a and b, or -1:
// the edge (a, b) of a 2D fragment then lies on that bisector.
int
_FindSharedBisector(const GEOGen::SymbolicVertex& a, const GEOGen::SymbolicVertex& b)
{
    for (GEO::index_t k = 0; k < a.nb_bisectors(); ++k) {
        const GEO::index_t j = a.bisector(GEO::signed_index_t(k));
        if (b.has_bisector(j)) {
            return int(j);
        }
    }
    return -1;
}

// 2D: one polygon per (cell, domain triangle) intersection.
class PolygonAccumulator : public GEO::RVDPolygonCallback
{
public:
    explicit PolygonAccumulator(Accumulator& acc) : _acc(acc) { }

    void operator()(GEO::index_t v, GEO::index_t t, const GEOGen::Polygon& P) const override
    {
        GEO::RVDPolygonCallback::operator()(v, t, P);
        const int i = int(v);
        if (i >= _acc.GetNumSites() or P.nb_vertices() < 3) {
            return;
        }  // sentinel site or degenerate fragment
        _acc.AddCellFragment(i, _ComputeFragmentMoments(i, P));
        _AddBisectorEdges(i, P);
    }

private:
    // Moments of the fragment about site i (fan triangles).
    Moments _ComputeFragmentMoments(int i, const GEOGen::Polygon& P) const
    {
        const double* q = _acc.GetSite(i);
        Moments m(2);
        double a[2], b[2], c[2];
        MakeRelative(P.vertex(0).point(), q, 2, a);
        for (GEO::index_t k = 1; k + 1 < P.nb_vertices(); ++k) {
            MakeRelative(P.vertex(k).point(), q, 2, b);
            MakeRelative(P.vertex(k + 1).point(), q, 2, c);
            const double* tri[3] = {a, b, c};
            m.AddSimplex(tri, 3, 0.5 * ComputeDoubleSignedArea(a, b, c));
        }
        return m;
    }

    // Edges whose extremities lie on the same bisector are facet pieces.
    void _AddBisectorEdges(int i, const GEOGen::Polygon& P) const
    {
        const double* q = _acc.GetSite(i);
        const GEO::index_t nv = P.nb_vertices();
        for (GEO::index_t k = 0; k < nv; ++k) {
            const GEOGen::Vertex& va = P.vertex(k);
            const GEOGen::Vertex& vb = P.vertex((k + 1) % nv);
            const int j = _FindSharedBisector(va.sym(), vb.sym());
            if (j < 0) {
                continue;
            }
            double a[2], b[2];
            MakeRelative(va.point(), q, 2, a);
            MakeRelative(vb.point(), q, 2, b);
            const double* seg[2] = {a, b};
            Moments m(2);
            m.AddSimplex(seg, 2, std::sqrt(ComputeDistance2(a, b, 2)));
            _acc.AddFacetPiece(i, j, m);
        }
    }

    Accumulator& _acc;
};

// 3D: one closed convex polyhedron per (cell, domain tet) intersection,
// reported facet by facet.
class PolyhedronAccumulator : public GEO::RVDPolyhedronCallback
{
public:
    explicit PolyhedronAccumulator(Accumulator& acc) : _acc(acc) { }

    void begin_polyhedron(GEO::index_t seed, GEO::index_t) override
    {
        _seed = int(seed);
        _facetPts.clear();
        _facetSizes.clear();
        _facetSeeds.clear();
    }

    void begin_facet(GEO::index_t facetSeed, GEO::index_t) override
    {
        _facetSeeds.push_back(facetSeed == GEO::index_t(-1) ? -1 : int(facetSeed));
        _facetSizes.push_back(0);
    }

    void vertex(const double* geometry, const GEOGen::SymbolicVertex&) override
    {
        // geometry is in the lifted space; keep x, y, z.
        _facetPts.insert(_facetPts.end(), geometry, geometry + 3);
        ++_facetSizes.back();
    }

    void end_polyhedron() override
    {
        if (_seed < 0 or _seed >= _acc.GetNumSites() or _facetPts.empty()) {
            return;
        }  // sentinel site
        Moments cell(3);
        size_t first = 0;
        for (size_t f = 0; f < _facetSizes.size(); ++f) {
            _AddFacet(first, _facetSizes[f], _facetSeeds[f], cell);
            first += size_t(_facetSizes[f]);
        }
        _acc.AddCellFragment(_seed, cell);
    }

private:
    // Point v of the current polyhedron, relative to the current site.
    void _GetRelativePoint(size_t v, double* r) const
    {
        MakeRelative(&_facetPts[3 * v], _acc.GetSite(_seed), 3, r);
    }

    // Adds the fan tets (apex: first point of the polyhedron) of facet
    // [first, first + size) to cell, and the facet itself to the shared
    // facet with site facetSeed (if any).
    void _AddFacet(size_t first, int size, int facetSeed, Moments& cell) const
    {
        double apex[3], a[3], b[3], c[3];
        _GetRelativePoint(0, apex);
        _GetRelativePoint(first, a);
        Moments facet(3);
        for (int k = 1; k + 1 < size; ++k) {
            _GetRelativePoint(first + k, b);
            _GetRelativePoint(first + k + 1, c);
            const double* tet[4] = {apex, a, b, c};
            cell.AddSimplex(tet, 4, ComputeSixSignedVolume(apex, a, b, c) / 6.0);

            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const double w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            double n[3];
            ComputeCross(u, w, n);
            const double* tri[3] = {a, b, c};
            facet.AddSimplex(tri, 3, 0.5 * std::sqrt(ComputeDot(n, n, 3)));
        }
        if (facetSeed >= 0) {
            _acc.AddFacetPiece(_seed, facetSeed, facet);
        }
    }

    Accumulator&        _acc;
    int                 _seed = -1;
    std::vector<double> _facetPts;
    std::vector<int>    _facetSizes;
    std::vector<int>    _facetSeeds;
};

// ---------------------------------------------------------------- polygons

using HalfEdgeMap = std::unordered_multimap<GEO::index_t, GEO::index_t>;

// Directed edge (a, b) packed in 64 bits.
uint64_t
_MakeEdgeKey(GEO::index_t a, GEO::index_t b)
{
    return (uint64_t(a) << 32) | uint64_t(b);
}

// Returns true if b is collinear with its neighbors a and c, and between
// them, or if an incident edge has zero length.
bool
_IsStraightVertex(const double* a, const double* b, const double* c)
{
    const double e1[2] = {b[0] - a[0], b[1] - a[1]};
    const double e2[2] = {c[0] - b[0], c[1] - b[1]};
    const double l1 = std::hypot(e1[0], e1[1]);
    const double l2 = std::hypot(e2[0], e2[1]);
    if (l1 == 0.0 or l2 == 0.0) {
        return true;
    }
    const double cross = e1[0] * e2[1] - e1[1] * e2[0];
    const double dot = e1[0] * e2[0] + e1[1] * e2[1];
    return dot > 0.0 and std::fabs(cross) <= 1e-10 * l1 * l2;
}

// Drops vertices where the boundary goes straight through (they come from
// the internal triangulation edges of the clipping domain).
void
_RemoveCollinear(std::vector<int>& loop, const std::vector<double>& pts)
{
    bool changed = true;
    while (changed and loop.size() > 3) {
        changed = false;
        for (size_t k = 0; k < loop.size() and loop.size() > 3; ++k) {
            const size_t prev = (k + loop.size() - 1) % loop.size();
            const size_t next = (k + 1) % loop.size();
            if (_IsStraightVertex(&pts[2 * loop[prev]], &pts[2 * loop[k]], &pts[2 * loop[next]])) {
                loop.erase(loop.begin() + k);
                changed = true;
                --k;
            }
        }
    }
}

// Facets of the restricted Voronoi diagram mesh, grouped by site.
std::vector<std::vector<GEO::index_t>>
_GroupFacetsBySite(GEO::Mesh& rvdMesh, GEO::index_t nbSites)
{
    GEO::Attribute<GEO::index_t> region(rvdMesh.facets.attributes(), "region");
    std::vector<std::vector<GEO::index_t>> cellFacets(nbSites);
    for (GEO::index_t f = 0; f < rvdMesh.facets.nb(); ++f) {
        if (region[f] < nbSites) {
            cellFacets[region[f]].push_back(f);
        }
    }
    return cellFacets;
}

// Vertices of facet f of M, oriented counter-clockwise.
std::vector<int>
_GetCCWFacet(const GEO::Mesh& M, GEO::index_t f, const std::vector<double>& pts)
{
    const GEO::index_t n = M.facets.nb_vertices(f);
    std::vector<int> frag(n);
    for (GEO::index_t lv = 0; lv < n; ++lv) {
        frag[lv] = int(M.facets.vertex(f, lv));
    }
    if (ComputeLoopSignedArea(frag, pts) < 0.0) {
        std::reverse(frag.begin(), frag.end());
    }
    return frag;
}

// Boundary half-edges of the union of the given facets: edges that are not
// shared by two of them. Facets are oriented counter-clockwise first, so
// shared edges cancel whatever the orientation of the input triangles.
HalfEdgeMap
_GetBoundaryHalfEdges(
    const GEO::Mesh& rvdMesh,
    const std::vector<GEO::index_t>& facets,
    const std::vector<double>& pts)
{
    std::unordered_map<uint64_t, int> count;  // (a, b) -> multiplicity
    for (GEO::index_t f : facets) {
        const std::vector<int> frag = _GetCCWFacet(rvdMesh, f, pts);
        for (size_t lv = 0; lv < frag.size(); ++lv) {
            const GEO::index_t a = GEO::index_t(frag[lv]);
            const GEO::index_t b = GEO::index_t(frag[(lv + 1) % frag.size()]);
            if (a == b) {
                continue;
            }
            auto opposite = count.find(_MakeEdgeKey(b, a));
            if (opposite == count.end()) {
                ++count[_MakeEdgeKey(a, b)];
            }
            else if (--opposite->second == 0) {
                count.erase(opposite);
            }
        }
    }

    HalfEdgeMap next;
    for (const auto& he : count) {
        for (int c = 0; c < he.second; ++c) {
            next.emplace(GEO::index_t(he.first >> 32), GEO::index_t(he.first & 0xffffffffu));
        }
    }
    return next;
}

// Chains boundary half-edges (consumed) into closed loops.
std::vector<std::vector<int>>
_ChainLoops(HalfEdgeMap& next)
{
    std::vector<std::vector<int>> loops;
    while (not next.empty()) {
        auto start = next.begin();
        const GEO::index_t first = start->first;
        GEO::index_t cur = start->second;
        next.erase(start);
        std::vector<int> loop{int(first)};
        while (cur != first) {
            loop.push_back(int(cur));
            auto it = next.find(cur);
            if (it == next.end()) {
                break;
            }  // open chain (should not happen)
            cur = it->second;
            next.erase(it);
        }
        loops.push_back(std::move(loop));
    }
    return loops;
}

void
_AppendPolygon2D(PolyMesh& out, const std::vector<int>& loop, const std::vector<double>& pts, int cell)
{
    out._polySizes.push_back(int(loop.size()));
    out._polyCell.push_back(cell);
    for (int v : loop) {
        out._polyVerts.push_back(int(out.GetNumPoints()));
        out._points.insert(out._points.end(), &pts[2 * v], &pts[2 * v] + 2);
    }
}

// Builds the boundary loops of each cell of a 2D restricted Voronoi diagram
// mesh (one fragment per cell/triangle intersection).
void
_ExtractCellBoundaries2D(GEO::Mesh& rvdMesh, GEO::index_t nbSites, PolyMesh& out)
{
    const std::vector<double> pts = GetMeshVertices(rvdMesh, 2);
    const std::vector<std::vector<GEO::index_t>> cellFacets = _GroupFacetsBySite(rvdMesh, nbSites);
    for (GEO::index_t s = 0; s < nbSites; ++s) {
        if (cellFacets[s].empty()) {
            continue;
        }
        HalfEdgeMap next = _GetBoundaryHalfEdges(rvdMesh, cellFacets[s], pts);
        // Outer boundaries come out counter-clockwise, boundaries of holes
        // inside the cell clockwise.
        for (std::vector<int>& loop : _ChainLoops(next)) {
            _RemoveCollinear(loop, pts);
            if (loop.size() >= 3) {
                _AppendPolygon2D(out, loop, pts, int(s));
            }
        }
    }
}

// Copies the facets of a 3D cell mesh built by BuildRVDMesh, skipping the
// sentinel sites.
void
_ExtractCellFacets3D(GEO::Mesh& cells, GEO::index_t nbSites, PolyMesh& out)
{
    GEO::Attribute<int> seedId(cells.facets.attributes(), "seed_id");
    out._points = GetMeshVertices(cells, 3);
    for (GEO::index_t f = 0; f < cells.facets.nb(); ++f) {
        const GEO::index_t n = cells.facets.nb_vertices(f);
        if (n < 3 or seedId[f] < 0 or GEO::index_t(seedId[f]) >= nbSites) {
            continue;
        }
        out._polySizes.push_back(int(n));
        out._polyCell.push_back(seedId[f]);
        for (GEO::index_t lv = 0; lv < n; ++lv) {
            out._polyVerts.push_back(int(cells.facets.vertex(f, lv)));
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
    bool IsVolumetric() const override { return false; }

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

    std::unique_ptr<GEO::Mesh> MakeWalkMesh(const GEO::Mesh& mesh, size_t minElements) const override
    {
        // A few triangles per part, so that the parts balance.
        std::vector<double> pts = GetMeshVertices(mesh, 2);
        std::vector<int> tris;
        for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f) {
            for (GEO::index_t lv = 0; lv < 3; ++lv) {
                tris.push_back(int(mesh.facets.vertex(f, lv)));
            }
        }
        while (tris.size() / 3 < minElements) {
            _SubdivideTriangles(pts, tris);
        }
        std::unique_ptr<GEO::Mesh> walk(new GEO::Mesh);
        SetMeshVertices(*walk, pts, 2, 3);
        AddMeshTriangles(*walk, tris);
        // Hilbert order, so that contiguous element ranges are compact.
        GEO::mesh_reorder(*walk, GEO::MESH_ORDER_HILBERT);
        walk->facets.connect();
        return walk;
    }

    GEO::index_t GetNumElements(const GEO::Mesh& M) const override { return M.facets.nb(); }
    GEO::index_t GetElementVertex(const GEO::Mesh& M, GEO::index_t e, GEO::index_t lv) const override
    { return M.facets.vertex(e, lv); }
    GEO::index_t GetNumAdjacent(const GEO::Mesh& M, GEO::index_t e) const override { return M.facets.nb_vertices(e); }
    GEO::index_t GetAdjacent(const GEO::Mesh& M, GEO::index_t e, GEO::index_t k) const override
    { return M.facets.adjacent(e, k); }

    void SetElementRange(GEO::RestrictedVoronoiDiagram& rvd, GEO::index_t first, GEO::index_t last) const override
    {
        rvd.set_facets_range(first, last);
    }

    void Walk(GEO::RestrictedVoronoiDiagram& rvd, Accumulator& acc) const override
    {
        PolygonAccumulator callback(acc);
        rvd.for_each_polygon(callback, true, false, false);
    }

    void ExtractPolygons(GEO::RestrictedVoronoiDiagram& rvd, GEO::index_t nbSites, double, PolyMesh& out) const override
    {
        GEO::Mesh rvdMesh;
        rvd.compute_RVD(rvdMesh);
        _ExtractCellBoundaries2D(rvdMesh, nbSites, out);
    }
};

// 3D: tetrahedral domains, lifted to (x, y, z, 0).
class CellBackend3D final : public CellBackend
{
public:
    int GetDimension() const override { return 3; }
    int GetElementSize() const override { return 4; }
    const char* GetDelaunayName() const override { return "BPOW"; }
    bool IsVolumetric() const override { return true; }

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
        // Positive tets (the RVD orients the output facets after them), with
        // tet-tet adjacency, then lifted to (x, y, z, 0).
        std::unique_ptr<GEO::Mesh> mesh(new GEO::Mesh);
        SetMeshVertices(*mesh, pts, 3, 3);
        AddMeshTets(*mesh, pts, elems);
        mesh->cells.connect();
        mesh->vertices.set_dimension(4);
        return mesh;
    }

    std::unique_ptr<GEO::Mesh> MakeWalkMesh(const GEO::Mesh& mesh, size_t) const override
    {
        std::unique_ptr<GEO::Mesh> walk(new GEO::Mesh);
        walk->copy(mesh);
        // Hilbert order, so that contiguous element ranges are compact.
        GEO::mesh_reorder(*walk, GEO::MESH_ORDER_HILBERT);
        walk->cells.connect();
        return walk;
    }

    GEO::index_t GetNumElements(const GEO::Mesh& M) const override { return M.cells.nb(); }
    GEO::index_t GetElementVertex(const GEO::Mesh& M, GEO::index_t e, GEO::index_t lv) const override
    { return M.cells.vertex(e, lv); }
    GEO::index_t GetNumAdjacent(const GEO::Mesh& M, GEO::index_t e) const override { return M.cells.nb_facets(e); }
    GEO::index_t GetAdjacent(const GEO::Mesh& M, GEO::index_t e, GEO::index_t k) const override
    { return M.cells.adjacent(e, k); }

    void SetElementRange(GEO::RestrictedVoronoiDiagram& rvd, GEO::index_t first, GEO::index_t last) const override
    {
        rvd.set_tetrahedra_range(first, last);
    }

    void Walk(GEO::RestrictedVoronoiDiagram& rvd, Accumulator& acc) const override
    {
        PolyhedronAccumulator callback(acc);
        rvd.for_each_polyhedron(callback, true, false, false);
    }

    void ExtractPolygons(
        GEO::RestrictedVoronoiDiagram& rvd,
        GEO::index_t nbSites,
        double boundaryAngle,
        PolyMesh& out) const override
    {
        GEO::Mesh cells;
        {
            GEO::BuildRVDMesh callback(cells);
            callback.set_simplify_boundary_facets(true, boundaryAngle);
            callback.set_generate_ids(true);
            rvd.for_each_polyhedron(callback);
        }
        _ExtractCellFacets3D(cells, nbSites, out);
    }
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
