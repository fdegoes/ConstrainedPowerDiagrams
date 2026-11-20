#include "core/powerCells3D.h"
#include "core/cellClip.h"
#include "core/cellMoments.h"
#include "core/domain.h"
#include "core/geometry.h"
#include "core/liftSites.h"
#include "core/powerDiagram.h"

#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <utility>

namespace PP
{

namespace
{

// Stride of the lifted sites (x, y, z, h) of LiftSites().
constexpr int _stride = 4;

// Relative tolerance under which a vertex is taken to lie on a cutting
// plane, so that it is neither duplicated nor split off as a sliver.
constexpr double _onPlaneTolerance = 1e-14;

// Plane sign * (n . (x - o) + c) = 0; the kept side is >= 0.
struct _Plane
{
    double _n[3], _o[3], _c, _sign;

    double Evaluate(const double* x) const
    {
        return _sign * (_n[0] * (x[0] - _o[0]) + _n[1] * (x[1] - _o[1]) + _n[2] * (x[2] - _o[2]) + _c);
    }
};

// Identifier of the plane of box face k (power half-spaces are identified by
// their neighbor site j >= 0 and domain faces f by -2 - f).
constexpr int _boxPlane = -(1 << 30);

// A convex polyhedron: vertices and faces, each a cycle of vertices oriented
// outwards (right-hand rule) with a label: the neighbor site of a power
// half-space, -1 on the starting box, or -2 - f on the plane of the domain
// face f. Each vertex records the planes it lies on, which identifies the
// same point in the fragments of a cell (to weld the cell's polygons).
// Faces are stored flat, and cuts reuse thread-local buffers.
struct _Polyhedron
{
    static constexpr int kMaxPlanes = 8;

    struct Vertex
    {
        double _x[3];
        int    _planes[kMaxPlanes];  // sorted
        int    _numPlanes;
    };

    std::vector<Vertex> _vertices;
    std::vector<int>    _faceStart{0};  // face f: [_faceStart[f], _faceStart[f + 1]) of _faceVerts
    std::vector<int>    _faceVerts;
    std::vector<int>    _faceLabels;

    int GetNumFaces() const { return int(_faceLabels.size()); }
    bool IsEmpty() const { return GetNumFaces() < 4; }
    int GetFaceSize(int f) const { return _faceStart[size_t(f) + 1] - _faceStart[size_t(f)]; }
    int GetFaceVertex(int f, int k) const { return _faceVerts[size_t(_faceStart[size_t(f)] + k)]; }
    int GetFaceLabel(int f) const { return _faceLabels[size_t(f)]; }
    const double* GetPoint(int v) const { return _vertices[size_t(v)]._x; }

    void MakeBox(const double* lo, const double* hi)
    {
        _vertices.clear();
        for (int k = 0; k < 8; ++k) {
            Vertex V;
            V._x[0] = (k & 1) ? hi[0] : lo[0];
            V._x[1] = (k & 2) ? hi[1] : lo[1];
            V._x[2] = (k & 4) ? hi[2] : lo[2];
            V._planes[0] = _boxPlane + ((k & 1) ? 1 : 0);
            V._planes[1] = _boxPlane + ((k & 2) ? 3 : 2);
            V._planes[2] = _boxPlane + ((k & 4) ? 5 : 4);
            V._numPlanes = 3;
            _vertices.push_back(V);
        }
        _faceVerts = {0, 4, 6, 2, 1, 3, 7, 5, 0, 1, 5, 4, 2, 6, 7, 3, 0, 2, 3, 1, 4, 5, 7, 6};
        _faceStart = {0, 4, 8, 12, 16, 20, 24};
        _faceLabels.assign(6, -1);
    }

    // Bounding box of the vertices of the faces.
    void GetBounds(double* lo, double* hi) const
    {
        for (int c = 0; c < 3; ++c) {
            lo[c] = HUGE_VAL;
            hi[c] = -HUGE_VAL;
        }
        for (int v : _faceVerts) {
            for (int c = 0; c < 3; ++c) {
                lo[c] = std::min(lo[c], GetPoint(v)[c]);
                hi[c] = std::max(hi[c], GetPoint(v)[c]);
            }
        }
    }

    // +1 if all the vertices are strictly on the side P > 0, -1 if none is,
    // 0 otherwise (the plane may cut).
    int GetSide(const _Plane& P) const
    {
        bool anyIn = false, anyOut = false;
        for (int v : _faceVerts) {
            const double d = P.Evaluate(GetPoint(v));
            anyIn = anyIn or d > 0.0;
            anyOut = anyOut or d <= 0.0;
            if (anyIn and anyOut) {
                return 0;
            }
        }
        return anyIn ? 1 : -1;
    }

    // Keeps the part on the side P >= 0 of the plane with identifier id,
    // closing it with a face of the given label. Vertices within tol of the
    // plane are taken on it.
    void Clip(int id, const _Plane& P, int label, double tol)
    {
        thread_local std::vector<double> g;
        thread_local std::vector<char> used;
        thread_local std::vector<std::pair<int64_t, int>> cut;
        thread_local std::vector<int> start, verts, labels, cap;
        thread_local std::vector<std::pair<double, int>> angles;
        const size_t nv = _vertices.size();
        g.assign(nv, 0.0);
        used.assign(nv, 0);
        bool anyIn = false, anyOut = false;
        for (int v : _faceVerts) {
            if (used[size_t(v)]) {
                continue;
            }
            used[size_t(v)] = 1;
            double d = P.Evaluate(GetPoint(v));
            if (std::fabs(d) <= tol) {
                d = 0.0;
            }
            g[size_t(v)] = d;
            anyIn = anyIn or d > 0.0;
            anyOut = anyOut or d < 0.0;
        }
        if (not anyOut) {
            return;
        }
        if (not anyIn) {
            _faceStart.assign(1, 0);
            _faceVerts.clear();
            _faceLabels.clear();
            return;
        }
        for (size_t v = 0; v < nv; ++v) {
            if (used[v] and g[v] == 0.0) {
                _AddPlane(_vertices[v], id);
            }
        }
        // One new vertex per cut edge, shared by its two faces.
        cut.clear();
        auto split = [&](int a, int b) {
            const int lo = std::min(a, b), hi = std::max(a, b);
            const int64_t key = (int64_t(lo) << 32) | int64_t(hi);
            for (const auto& c : cut) {
                if (c.first == key) {
                    return c.second;
                }
            }
            cut.emplace_back(key, int(_vertices.size()));
            _vertices.push_back(_MakeCut(lo, hi, g[size_t(lo)], g[size_t(hi)], id));
            g.push_back(0.0);
            return int(_vertices.size() - 1);
        };
        start.assign(1, 0);
        verts.clear();
        labels.clear();
        bool onPlaneFace = false;
        for (int f = 0; f < GetNumFaces(); ++f) {
            const int m = GetFaceSize(f);
            const size_t first = verts.size();
            bool allOn = true;
            for (int k = 0; k < m; ++k) {
                const int a = GetFaceVertex(f, k), b = GetFaceVertex(f, (k + 1) % m);
                const double ga = g[size_t(a)], gb = g[size_t(b)];
                allOn = allOn and ga == 0.0;
                if (ga >= 0.0) {
                    verts.push_back(a);
                }
                if ((ga > 0.0 and gb < 0.0) or (ga < 0.0 and gb > 0.0)) {
                    verts.push_back(split(a, b));
                }
            }
            onPlaneFace = onPlaneFace or allOn;
            if (verts.size() - first >= 3) {
                labels.push_back(GetFaceLabel(f));
                start.push_back(int(verts.size()));
            } else {
                verts.resize(first);
            }
        }
        std::swap(_faceStart, start);
        std::swap(_faceVerts, verts);
        std::swap(_faceLabels, labels);
        if (onPlaneFace) {
            return;  // a face already lies on the plane: it is the closing face
        }
        // The closing face: the vertices on the plane, in counter-clockwise
        // order seen from outside (-P's gradient).
        cap.clear();
        used.assign(_vertices.size(), 0);
        for (int v : _faceVerts) {
            if (g[size_t(v)] == 0.0 and not used[size_t(v)]) {
                used[size_t(v)] = 1;
                cap.push_back(v);
            }
        }
        if (cap.size() < 3) {
            return;
        }
        const double m[3] = {-P._sign * P._n[0], -P._sign * P._n[1], -P._sign * P._n[2]};
        int axis = 0;
        for (int c = 1; c < 3; ++c) {
            if (std::fabs(m[c]) < std::fabs(m[axis])) {
                axis = c;
            }
        }
        double e[3] = {0.0, 0.0, 0.0}, u[3], w[3];
        e[axis] = 1.0;
        ComputeCross(m, e, u);
        ComputeCross(m, u, w);
        double center[3] = {0.0, 0.0, 0.0};
        for (int v : cap) {
            for (int c = 0; c < 3; ++c) {
                center[c] += GetPoint(v)[c] / double(cap.size());
            }
        }
        angles.clear();
        for (int v : cap) {
            double d[3];
            for (int c = 0; c < 3; ++c) {
                d[c] = GetPoint(v)[c] - center[c];
            }
            angles.emplace_back(std::atan2(ComputeDot(d, w, 3), ComputeDot(d, u, 3)), v);
        }
        std::sort(angles.begin(), angles.end());
        for (const auto& a : angles) {
            _faceVerts.push_back(a.second);
        }
        _faceLabels.push_back(label);
        _faceStart.push_back(int(_faceVerts.size()));
    }

private:
    static void _AddPlane(Vertex& V, int id)
    {
        int k = 0;
        while (k < V._numPlanes and V._planes[k] < id) {
            ++k;
        }
        if ((k < V._numPlanes and V._planes[k] == id) or V._numPlanes == kMaxPlanes) {
            return;
        }
        for (int l = V._numPlanes; l > k; --l) {
            V._planes[l] = V._planes[l - 1];
        }
        V._planes[k] = id;
        ++V._numPlanes;
    }

    // The point where the edge (a, b), with a < b, crosses the plane id,
    // interpolated along the edge from the smaller index (solving the planes
    // of the edge and id is ill-conditioned when a face is nearly parallel
    // to the plane); it lies on the planes of the edge and on id.
    Vertex _MakeCut(int a, int b, double ga, double gb, int id) const
    {
        const Vertex& A = _vertices[size_t(a)];
        const Vertex& B = _vertices[size_t(b)];
        Vertex V;
        V._numPlanes = int(std::set_intersection(A._planes, A._planes + A._numPlanes, B._planes,
                                                 B._planes + B._numPlanes, V._planes) - V._planes);
        const double t = ga / (ga - gb);
        for (int c = 0; c < 3; ++c) {
            V._x[c] = A._x[c] + t * (B._x[c] - A._x[c]);
        }
        _AddPlane(V, id);
        return V;
    }
};

// The domain tetrahedra: each face plane once (from its sorted vertices), so
// that the two tets sharing a face cut along the same plane, with the side
// of each tet, and a grid of their bounding boxes.
class _DomainTets
{
public:
    explicit _DomainTets(const Domain& domain) : _grid(_Build(domain)) { }

    size_t GetNumTets() const { return _sides.size(); }
    bool IsBoundaryFace(int f) const { return _boundary[size_t(f)] != 0; }
    // The 3 mesh vertices of face f, sorted.
    const std::array<GEO::index_t, 3>& GetFaceVertices(int f) const { return _faceVertices[size_t(f)]; }

    // The 4 planes of tet t (kept side: inside) and their face indices.
    _Plane GetPlane(int t, int k, int& face) const
    {
        face = _faces[size_t(t)][size_t(k)];
        _Plane P = _planes[size_t(face)];
        P._sign = _sides[size_t(t)][size_t(k)];
        return P;
    }

    void Find(const double* lo, const double* hi, int tag, std::vector<int>& stamp, std::vector<int>& out) const
    { _grid.Find(lo, hi, tag, stamp, out); }

    const double* GetLow() const { return _lo; }
    const double* GetHigh() const { return _hi; }

private:
    BoxGrid<3> _Build(const Domain& domain)
    {
        const GEO::Mesh& M = domain.GetMesh();
        std::map<std::array<GEO::index_t, 3>, int> index;
        std::vector<double> lo, hi;
        for (int c = 0; c < 3; ++c) {
            _lo[c] = HUGE_VAL;
            _hi[c] = -HUGE_VAL;
        }
        for (GEO::index_t t = 0; t < M.cells.nb(); ++t) {
            std::array<int, 4> faces;
            std::array<double, 4> sides;
            for (GEO::index_t lf = 0; lf < 4; ++lf) {
                // Facet lf of a tet is opposite its vertex lf.
                std::array<GEO::index_t, 3> key = {M.cells.facet_vertex(t, lf, 0), M.cells.facet_vertex(t, lf, 1),
                                                   M.cells.facet_vertex(t, lf, 2)};
                std::sort(key.begin(), key.end());
                auto it = index.find(key);
                if (it == index.end()) {
                    const double* p0 = M.vertices.point_ptr(key[0]);
                    const double* p1 = M.vertices.point_ptr(key[1]);
                    const double* p2 = M.vertices.point_ptr(key[2]);
                    _Plane P;
                    const double a[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
                    const double b[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
                    ComputeCross(a, b, P._n);
                    std::copy(p0, p0 + 3, P._o);
                    P._c = 0.0;
                    P._sign = 1.0;
                    it = index.emplace(key, int(_planes.size())).first;
                    _planes.push_back(P);
                    _boundary.push_back(M.cells.adjacent(t, lf) == GEO::NO_INDEX);
                    _faceVertices.push_back(key);
                }
                faces[lf] = it->second;
                const _Plane& P = _planes[size_t(it->second)];
                sides[lf] = P.Evaluate(M.vertices.point_ptr(M.cells.vertex(t, lf))) > 0.0 ? 1.0 : -1.0;
            }
            _faces.push_back(faces);
            _sides.push_back(sides);
            for (int c = 0; c < 3; ++c) {
                double l = HUGE_VAL, h = -HUGE_VAL;
                for (GEO::index_t lv = 0; lv < 4; ++lv) {
                    const double x = M.vertices.point_ptr(M.cells.vertex(t, lv))[c];
                    l = std::min(l, x);
                    h = std::max(h, x);
                }
                lo.push_back(l);
                hi.push_back(h);
                _lo[c] = std::min(_lo[c], l);
                _hi[c] = std::max(_hi[c], h);
            }
        }
        return BoxGrid<3>(lo, hi);
    }

    std::vector<_Plane>                _planes;
    std::vector<char>                  _boundary;
    std::vector<std::array<GEO::index_t, 3>> _faceVertices;
    std::vector<std::array<int, 4>>    _faces;
    std::vector<std::array<double, 4>> _sides;
    double                             _lo[3], _hi[3];
    BoxGrid<3>                         _grid;
};

// Visits, for each site i in [first, last) with a non-empty cell,
// visit(i, fragments) with the cell clipped to each domain tet it overlaps.
template <class Visit>
void
_ForEachCell(const PowerTriangulation& pd, const _DomainTets& domain, int first, int last, const Visit& visit)
{
    const GEO::Delaunay& T = *pd._delaunay;
    const double* dlo = domain.GetLow();
    const double* dhi = domain.GetHigh();
    double box[2][3], diag2 = 0.0;
    for (int c = 0; c < 3; ++c) {
        const double pad = 0.01 * (dhi[c] - dlo[c]) + 1e-12;
        box[0][c] = dlo[c] - pad;
        box[1][c] = dhi[c] + pad;
        diag2 += (box[1][c] - box[0][c]) * (box[1][c] - box[0][c]);
    }
    const double scale = std::sqrt(diag2);
    GEO::vector<GEO::index_t> neighbors;
    std::vector<std::pair<double, GEO::index_t>> order;
    std::vector<int> stamp(domain.GetNumTets(), -1), candidates;
    std::vector<_Polyhedron> fragments;
    _Polyhedron cell;
    for (int i = first; i < last; ++i) {
        if (T.vertex_cell(GEO::index_t(i)) == GEO::NO_INDEX) {
            continue;  // hidden site: empty cell
        }
        // The power half-spaces of the neighbors, the nearest first.
        T.get_neighbors(GEO::index_t(i), neighbors);
        const double* si = &pd._lifted[_stride * size_t(i)];
        order.clear();
        for (GEO::index_t j : neighbors) {
            order.emplace_back(ComputeDistance2(si, &pd._lifted[_stride * size_t(j)], 3), j);
        }
        std::sort(order.begin(), order.end());
        cell.MakeBox(box[0], box[1]);
        for (const auto& nb : order) {
            const int j = int(nb.second);
            // Canonical plane of the pair (a, b), a < b: power_b - power_a,
            // relative to s_a, with power_k(x) = |x - s_k|^2 + h_k^2.
            const int a = std::min(i, j), b = std::max(i, j);
            const double* sa = &pd._lifted[_stride * size_t(a)];
            const double* sb = &pd._lifted[_stride * size_t(b)];
            _Plane P;
            double d2 = 0.0;
            for (int c = 0; c < 3; ++c) {
                P._n[c] = 2.0 * (sa[c] - sb[c]);
                P._o[c] = sa[c];
                d2 += (sb[c] - sa[c]) * (sb[c] - sa[c]);
            }
            P._c = d2 + sb[3] * sb[3] - sa[3] * sa[3];
            P._sign = i == a ? 1.0 : -1.0;
            cell.Clip(j, P, j, _onPlaneTolerance * 2.0 * std::sqrt(d2) * scale);
            if (cell.IsEmpty()) {
                break;
            }
        }
        if (cell.IsEmpty()) {
            continue;
        }
        double lo[3], hi[3];
        cell.GetBounds(lo, hi);
        domain.Find(lo, hi, i, stamp, candidates);
        fragments.clear();
        for (int t : candidates) {
            // The planes of tet t that cut the cell; none if the cell is
            // inside the tet (then no other tet meets it), skip the tet if
            // the cell is outside one of them.
            int cutting[4], numCutting = 0;
            bool outside = false;
            for (int k = 0; k < 4 and not outside; ++k) {
                int face;
                const _Plane P = domain.GetPlane(t, k, face);
                const int side = cell.GetSide(P);
                outside = side < 0;
                if (side == 0) {
                    cutting[numCutting++] = k;
                }
            }
            if (outside) {
                continue;
            }
            if (numCutting == 0) {
                fragments.assign(1, cell);
                break;
            }
            _Polyhedron f = cell;
            for (int c = 0; c < numCutting and not f.IsEmpty(); ++c) {
                int face;
                const _Plane P = domain.GetPlane(t, cutting[c], face);
                const double norm = std::sqrt(ComputeDot(P._n, P._n, 3));
                f.Clip(-2 - face, P, -2 - face, _onPlaneTolerance * norm * scale);
            }
            if (not f.IsEmpty()) {
                fragments.push_back(std::move(f));
            }
        }
        if (not fragments.empty()) {
            visit(i, fragments);
        }
    }
}

bool
_IntegrateCells(const PowerTriangulation& pd, const _DomainTets& domain, int first, int last, Accumulator& acc)
{
    const int n = acc.GetNumSites();
    _ForEachCell(pd, domain, first, last, [&](int i, const std::vector<_Polyhedron>& fragments) {
        const double* q = acc.GetSite(i);
        const double zero[3] = {0.0, 0.0, 0.0};
        for (const _Polyhedron& f : fragments) {
            // Fan tets from site i over the faces (outward), and the faces
            // on bisectors of real neighbors.
            Moments m(3);
            for (int F = 0; F < f.GetNumFaces(); ++F) {
                const int label = f.GetFaceLabel(F);
                Moments facet(3);
                double a[3], b[3], c[3];
                MakeRelative(f.GetPoint(f.GetFaceVertex(F, 0)), q, 3, a);
                for (int k = 1; k + 1 < f.GetFaceSize(F); ++k) {
                    MakeRelative(f.GetPoint(f.GetFaceVertex(F, k)), q, 3, b);
                    MakeRelative(f.GetPoint(f.GetFaceVertex(F, k + 1)), q, 3, c);
                    const double* tet[4] = {zero, a, b, c};
                    m.AddSimplex(tet, 4, ComputeDeterminant3(a, b, c) / 6.0);
                    if (label > i and label < n) {
                        const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
                        const double w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
                        double x[3];
                        ComputeCross(u, w, x);
                        const double* tri[3] = {a, b, c};
                        facet.AddSimplex(tri, 3, 0.5 * std::sqrt(ComputeDot(x, x, 3)));
                    }
                }
                if (label > i and label < n) {
                    acc.AddFacetPiece(i, label, facet);
                }
            }
            acc.AddCellFragment(i, m);
        }
    });
    return true;
}

// ---------------------------------------------------------------- polygons

// Loops of the union of oriented polygons (indices into the welded points of
// a cell): shared edges cancel and the remaining ones are chained.
std::vector<std::vector<int>>
_MergeLoops(const std::vector<std::vector<int>>& pieces)
{
    std::map<std::pair<int, int>, int> count;
    for (const std::vector<int>& p : pieces) {
        for (size_t k = 0; k < p.size(); ++k) {
            const int a = p[k], b = p[(k + 1) % p.size()];
            if (a == b) {
                continue;
            }
            const auto opposite = count.find({b, a});
            if (opposite != count.end()) {
                if (--opposite->second == 0) {
                    count.erase(opposite);
                }
            } else {
                ++count[{a, b}];
            }
        }
    }
    std::multimap<int, int> next;
    for (const auto& e : count) {
        for (int c = 0; c < e.second; ++c) {
            next.emplace(e.first.first, e.first.second);
        }
    }
    std::vector<std::vector<int>> loops;
    while (not next.empty()) {
        auto start = next.begin();
        const int head = start->first;
        int cur = start->second;
        next.erase(start);
        std::vector<int> loop{head};
        while (cur != head) {
            loop.push_back(cur);
            const auto it = next.find(cur);
            if (it == next.end()) {
                break;  // open chain (not expected)
            }
            cur = it->second;
            next.erase(it);
        }
        loops.push_back(std::move(loop));
    }
    return loops;
}

// Unit normal of a polygon (Newell's method).
void
_ComputeNormal(const std::vector<int>& loop, const std::vector<double>& pts, double* nrm)
{
    nrm[0] = nrm[1] = nrm[2] = 0.0;
    for (size_t k = 0; k < loop.size(); ++k) {
        const double* p = &pts[3 * size_t(loop[k])];
        const double* q = &pts[3 * size_t(loop[(k + 1) % loop.size()])];
        nrm[0] += (p[1] - q[1]) * (p[2] + q[2]);
        nrm[1] += (p[2] - q[2]) * (p[0] + q[0]);
        nrm[2] += (p[0] - q[0]) * (p[1] + q[1]);
    }
    const double l = std::sqrt(ComputeDot(nrm, nrm, 3));
    for (int c = 0; c < 3; ++c) {
        nrm[c] = l > 0.0 ? nrm[c] / l : 0.0;
    }
}

void
_AppendCellPolygons(
    int i,
    int n,
    const _DomainTets& domain,
    double cosAngle,
    const std::vector<_Polyhedron>& fragments,
    PolyMesh& out)
{
    // Weld the vertices of the fragments, relative to the cell size.
    double lo[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL}, hi[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};
    for (const _Polyhedron& f : fragments) {
        double l[3], h[3];
        f.GetBounds(l, h);
        for (int c = 0; c < 3; ++c) {
            lo[c] = std::min(lo[c], l[c]);
            hi[c] = std::max(hi[c], h[c]);
        }
    }
    const double tol = 1e-10 * std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
    std::vector<double> pts;
    // A vertex is identified by the power half-spaces it lies on and by the
    // mesh simplex (face, edge or vertex) common to its domain face planes,
    // the same in every tet around it; vertices that this does not determine
    // are welded within tol.
    std::map<std::vector<int>, int> byPlanes;
    auto identify = [&](const _Polyhedron::Vertex& V, std::vector<int>& key) {
        key.clear();
        std::vector<GEO::index_t> simplex;
        bool onMesh = false;
        for (int k = 0; k < V._numPlanes; ++k) {
            const int p = V._planes[k];
            if (p >= 0 or p < -(1 << 29)) {
                key.push_back(p);  // power half-space or box
                continue;
            }
            const std::array<GEO::index_t, 3>& f = domain.GetFaceVertices(-2 - p);
            if (not onMesh) {
                simplex.assign(f.begin(), f.end());
                onMesh = true;
            } else {
                std::vector<GEO::index_t> common;
                std::set_intersection(simplex.begin(), simplex.end(), f.begin(), f.end(), std::back_inserter(common));
                simplex = std::move(common);
            }
        }
        // Constraints: one per half-space, and 1, 2 or 3 for a mesh face,
        // edge or vertex.
        const size_t constraints = key.size() + (onMesh and not simplex.empty() ? 4 - simplex.size() : 0);
        if (constraints < 3 or (onMesh and simplex.empty())) {
            return false;
        }
        key.push_back(-1);
        for (GEO::index_t v : simplex) {
            key.push_back(int(v));
        }
        return true;
    };
    std::vector<int> key;
    auto weld = [&](const _Polyhedron::Vertex& V) {
        const bool known = identify(V, key);
        if (known) {
            const auto it = byPlanes.find(key);
            if (it != byPlanes.end()) {
                return it->second;
            }
        }
        const double* x = V._x;
        int k = -1;
        for (size_t q = 0; q < pts.size() / 3 and k < 0; ++q) {
            if (std::fabs(pts[3 * q] - x[0]) <= tol and std::fabs(pts[3 * q + 1] - x[1]) <= tol and
                std::fabs(pts[3 * q + 2] - x[2]) <= tol) {
                k = int(q);
            }
        }
        if (k < 0) {
            pts.insert(pts.end(), x, x + 3);
            k = int(pts.size() / 3 - 1);
        }
        if (known) {
            byPlanes.emplace(key, k);
        }
        return k;
    };
    // Pieces on the facets (by neighbor) and on the domain boundary; the
    // pieces on internal tet faces and on the box are dropped.
    std::map<int, std::vector<std::vector<int>>> facets;
    std::vector<std::vector<int>> boundary;
    for (const _Polyhedron& f : fragments) {
        for (int F = 0; F < f.GetNumFaces(); ++F) {
            const int label = f.GetFaceLabel(F);
            const bool onFacet = label >= 0 and label < n;
            const bool onBoundary = label <= -2 and label > -(1 << 29) and domain.IsBoundaryFace(-2 - label);
            if (not onFacet and not onBoundary) {
                continue;
            }
            std::vector<int> loop;
            for (int k = 0; k < f.GetFaceSize(F); ++k) {
                const int w = weld(f._vertices[size_t(f.GetFaceVertex(F, k))]);
                if (loop.empty() or loop.back() != w) {
                    loop.push_back(w);
                }
            }
            while (loop.size() > 1 and loop.front() == loop.back()) {
                loop.pop_back();
            }
            if (loop.size() < 3) {
                continue;
            }
            (onFacet ? facets[label] : boundary).push_back(std::move(loop));
        }
    }
    std::map<int, int> outIndex;  // welded point -> point of out
    auto emit = [&](std::vector<int> loop) {
        RemoveCollinearVertices3D(loop, pts);
        if (loop.size() < 3) {
            return;
        }
        out._polySizes.push_back(int(loop.size()));
        out._polyCell.push_back(i);
        for (int v : loop) {
            auto it = outIndex.find(v);
            if (it == outIndex.end()) {
                it = outIndex.emplace(v, int(out.GetNumPoints())).first;
                out._points.insert(out._points.end(), &pts[3 * size_t(v)], &pts[3 * size_t(v)] + 3);
            }
            out._polyVerts.push_back(it->second);
        }
    };
    auto emitGroup = [&](const std::vector<std::vector<int>>& pieces) {
        const std::vector<std::vector<int>> loops = _MergeLoops(pieces);
        if (loops.size() == 1) {
            emit(loops[0]);
        } else {
            for (const std::vector<int>& p : pieces) {
                emit(p);
            }
        }
    };
    for (const auto& f : facets) {
        emitGroup(f.second);
    }
    // Boundary patches: pieces joined across shared edges whose normals are
    // within the merge angle.
    const size_t m = boundary.size();
    std::vector<std::array<double, 3>> normals(m);
    for (size_t k = 0; k < m; ++k) {
        _ComputeNormal(boundary[k], pts, normals[k].data());
    }
    std::map<std::pair<int, int>, size_t> edgeOwner;
    for (size_t k = 0; k < m; ++k) {
        for (size_t e = 0; e < boundary[k].size(); ++e) {
            edgeOwner[{boundary[k][e], boundary[k][(e + 1) % boundary[k].size()]}] = k;
        }
    }
    std::vector<size_t> parent(m);
    std::iota(parent.begin(), parent.end(), size_t(0));
    std::function<size_t(size_t)> root = [&](size_t k) { return parent[k] == k ? k : parent[k] = root(parent[k]); };
    for (size_t k = 0; k < m; ++k) {
        for (size_t e = 0; e < boundary[k].size(); ++e) {
            const auto it = edgeOwner.find({boundary[k][(e + 1) % boundary[k].size()], boundary[k][e]});
            if (it != edgeOwner.end() and ComputeDot(normals[k].data(), normals[it->second].data(), 3) >= cosAngle) {
                parent[root(k)] = root(it->second);
            }
        }
    }
    std::map<size_t, std::vector<std::vector<int>>> patches;
    for (size_t k = 0; k < m; ++k) {
        patches[root(k)].push_back(boundary[k]);
    }
    for (const auto& p : patches) {
        emitGroup(p.second);
    }
}

} // namespace

void
ComputePowerCells3D(
    const Domain& domain,
    const std::vector<double>& sites,
    const PowerTriangulation& pd,
    bool parallel,
    Accumulator& acc)
{
    const _DomainTets tets(domain);
    ForEachSiteChunk(acc.GetNumSites(), parallel, sites, 3, acc, [&](int first, int last, Accumulator& part) {
        return _IntegrateCells(pd, tets, first, last, part);
    });
}

void
ComputePowerPolygons3D(
    const Domain& domain,
    int numSites,
    const PowerTriangulation& pd,
    double boundaryAngle,
    PolyMesh& out)
{
    const _DomainTets tets(domain);
    // Coplanar within round-off at angle 0.
    const double cosAngle = std::min(std::cos(boundaryAngle * M_PI / 180.0), 1.0 - 1e-12);
    _ForEachCell(pd, tets, 0, numSites, [&](int i, const std::vector<_Polyhedron>& fragments) {
        _AppendCellPolygons(i, numSites, tets, cosAngle, fragments, out);
    });
}

} // namespace PP
