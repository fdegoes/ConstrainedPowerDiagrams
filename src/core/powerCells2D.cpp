#include "core/powerCells2D.h"
#include "core/cellClip.h"
#include "core/cellMoments.h"
#include "core/domain.h"
#include "core/geometry.h"
#include "core/liftSites.h"
#include "core/powerDiagram.h"

#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <thread>
#include <utility>

namespace PP
{

namespace
{

// Stride of the lifted sites (x, y, h) of LiftSites().
constexpr int _stride = 3;

// A vertex of a cell polygon and the label of its outgoing edge: the
// neighbor site across it (>= the number of sites: a sentinel), or -1 on a
// domain triangle edge.
struct _Vertex
{
    double _x[2];
    int    _label;
};

// Supporting line of a domain triangle edge, with a canonical orientation
// (from the lexicographically smaller endpoint), so that the two triangles
// sharing an edge evaluate it, and its intersections, identically.
struct _Line
{
    double _o[2];   // origin (smaller endpoint)
    double _d[2];   // direction (to the larger endpoint)
    double _sign;   // +1 or -1: the triangle is on the side sign * f >= 0

    double Evaluate(const double* x) const
    { return _d[0] * (x[1] - _o[1]) - _d[1] * (x[0] - _o[0]); }
};

struct _Triangle
{
    _Line _edges[3];
};

// The domain triangles, counter-clockwise, and a grid of their bounding
// boxes.
class _DomainTriangles
{
public:
    explicit _DomainTriangles(const Domain& domain) : _grid(_Build(domain, _triangles)) { }

    const _Triangle& Get(int t) const { return _triangles[size_t(t)]; }
    size_t GetNumTriangles() const { return _triangles.size(); }

    // The triangles whose bounding box overlaps [lo, hi] (see BoxGrid).
    void Find(const double* lo, const double* hi, int tag, std::vector<int>& stamp, std::vector<int>& out) const
    { _grid.Find(lo, hi, tag, stamp, out); }

private:
    static BoxGrid<2> _Build(const Domain& domain, std::vector<_Triangle>& triangles)
    {
        const GEO::Mesh& M = domain.GetMesh();
        std::vector<double> lo, hi;
        for (GEO::index_t f = 0; f < M.facets.nb(); ++f) {
            double p[3][2];
            for (GEO::index_t lv = 0; lv < 3; ++lv) {
                const double* x = M.vertices.point_ptr(M.facets.vertex(f, lv));
                p[lv][0] = x[0];
                p[lv][1] = x[1];
            }
            if (ComputeDoubleSignedArea(p[0], p[1], p[2]) < 0.0) {
                std::swap(p[1], p[2]);
            }
            _Triangle T;
            for (int k = 0; k < 3; ++k) {
                const double* a = p[k];
                const double* b = p[(k + 1) % 3];
                const bool forward = a[0] < b[0] or (a[0] == b[0] and a[1] < b[1]);
                const double* o = forward ? a : b;
                const double* e = forward ? b : a;
                T._edges[k] = {{o[0], o[1]}, {e[0] - o[0], e[1] - o[1]}, forward ? 1.0 : -1.0};
            }
            for (int c = 0; c < 2; ++c) {
                lo.push_back(std::min({p[0][c], p[1][c], p[2][c]}));
                hi.push_back(std::max({p[0][c], p[1][c], p[2][c]}));
            }
            triangles.push_back(T);
        }
        return BoxGrid<2>(lo, hi);
    }

    std::vector<_Triangle> _triangles;
    BoxGrid<2>             _grid;
};

// Power vertex (radical center) of the regular triangle (a, b, c) of the
// lifted sites: the point with equal power distance to the three sites;
// false if they are collinear.
bool
_ComputePowerVertex(const double* a, const double* b, const double* c, double* x)
{
    // 2 (s_k - s_a) . y = |s_k - s_a|^2 + h_k^2 - h_a^2, with x = s_a + y.
    const double u[2] = {b[0] - a[0], b[1] - a[1]};
    const double v[2] = {c[0] - a[0], c[1] - a[1]};
    const double ru = 0.5 * (u[0] * u[0] + u[1] * u[1] + b[2] * b[2] - a[2] * a[2]);
    const double rv = 0.5 * (v[0] * v[0] + v[1] * v[1] + c[2] * c[2] - a[2] * a[2]);
    const double det = u[0] * v[1] - u[1] * v[0];
    if (not (det != 0.0)) {
        return false;
    }
    x[0] = a[0] + (ru * v[1] - rv * u[1]) / det;
    x[1] = a[1] + (u[0] * rv - v[0] * ru) / det;
    return true;
}

// Clips the convex polygon in against the half-plane of L (Sutherland-Hodgman),
// keeping the edge labels; the new edge along L gets label -1.
void
_Clip(const std::vector<_Vertex>& in, const _Line& L, std::vector<_Vertex>& out)
{
    out.clear();
    const size_t m = in.size();
    for (size_t k = 0; k < m; ++k) {
        const _Vertex& A = in[k];
        const _Vertex& B = in[(k + 1) % m];
        const double fa = L.Evaluate(A._x), fb = L.Evaluate(B._x);
        const bool ain = L._sign * fa >= 0.0, bin = L._sign * fb >= 0.0;
        if (ain) {
            out.push_back(A);
        }
        if (ain != bin) {
            const double t = fa / (fa - fb);
            _Vertex I = {{A._x[0] + t * (B._x[0] - A._x[0]), A._x[1] + t * (B._x[1] - A._x[1])}, ain ? -1 : A._label};
            out.push_back(I);
        }
    }
}

// Visits the fragments of the cells of the sites [first, last): for each
// site i with a non-empty cell, visit(i, fragments) with the cell clipped to
// each domain triangle it overlaps (convex, counter-clockwise, with edge
// labels); false if a cell is unbounded.
template <class Visit>
bool
_ForEachCell(
    const PowerTriangulation& pd,
    const std::vector<double>& powerVertices,
    const _DomainTriangles& domain,
    int first,
    int last,
    const Visit& visit)
{
    const GEO::Delaunay& T = *pd._delaunay;
    std::vector<_Vertex> cell, a, b;
    std::vector<std::vector<_Vertex>> fragments;
    std::vector<int> stamp(domain.GetNumTriangles(), -1), candidates;
    for (int i = first; i < last; ++i) {
        const GEO::index_t t0 = T.vertex_cell(GEO::index_t(i));
        if (t0 == GEO::NO_INDEX) {
            continue;  // hidden site: empty cell
        }
        // The triangles around i, counter-clockwise: in the triangle
        // (i, b, c), the next one is across the edge (i, c), opposite b, and
        // the cell edge between their power vertices lies on the bisector of
        // i and c.
        cell.clear();
        GEO::index_t t = t0;
        do {
            const GEO::index_t lv = T.index(t, GEO::index_t(i));
            const _Vertex V = {{powerVertices[2 * t], powerVertices[2 * t + 1]},
                               int(T.cell_vertex(t, (lv + 2) % 3))};
            cell.push_back(V);
            t = T.cell_adjacent(t, (lv + 1) % 3);
            if (t == GEO::NO_INDEX or cell.size() > T.nb_cells()) {
                return false;  // on the hull: unbounded cell (not expected)
            }
        } while (t != t0);

        double lo[2] = {HUGE_VAL, HUGE_VAL}, hi[2] = {-HUGE_VAL, -HUGE_VAL};
        for (const _Vertex& V : cell) {
            for (int c = 0; c < 2; ++c) {
                lo[c] = std::min(lo[c], V._x[c]);
                hi[c] = std::max(hi[c], V._x[c]);
            }
        }
        domain.Find(lo, hi, i, stamp, candidates);
        fragments.clear();
        for (int k : candidates) {
            const _Triangle& D = domain.Get(k);
            a = cell;
            for (int e = 0; e < 3 and a.size() >= 3; ++e) {
                _Clip(a, D._edges[e], b);
                std::swap(a, b);
            }
            if (a.size() >= 3) {
                fragments.push_back(a);
            }
        }
        if (not fragments.empty()) {
            visit(i, fragments);
        }
    }
    return true;
}

// Moments of the cells of the sites [first, last) and of their facets.
bool
_IntegrateCells(
    const PowerTriangulation& pd,
    const std::vector<double>& powerVertices,
    const _DomainTriangles& domain,
    int first,
    int last,
    Accumulator& acc)
{
    const int n = acc.GetNumSites();
    return _ForEachCell(pd, powerVertices, domain, first, last, [&](int i, const std::vector<std::vector<_Vertex>>& fragments) {
        const double* q = acc.GetSite(i);
        for (const std::vector<_Vertex>& a : fragments) {
            // Moments of the fragment about site i (fan triangles), and its
            // edges on bisectors of real neighbors.
            Moments m(2);
            double r0[2], r1[2], r2[2];
            MakeRelative(a[0]._x, q, 2, r0);
            for (size_t v = 1; v + 1 < a.size(); ++v) {
                MakeRelative(a[v]._x, q, 2, r1);
                MakeRelative(a[v + 1]._x, q, 2, r2);
                const double* tri[3] = {r0, r1, r2};
                m.AddSimplex(tri, 3, 0.5 * ComputeDoubleSignedArea(r0, r1, r2));
            }
            acc.AddCellFragment(i, m);
            for (size_t v = 0; v < a.size(); ++v) {
                const int j = a[v]._label;
                if (j <= i or j >= n) {
                    continue;  // domain edge, sentinel, or recorded by j
                }
                MakeRelative(a[v]._x, q, 2, r1);
                MakeRelative(a[(v + 1) % a.size()]._x, q, 2, r2);
                const double* seg[2] = {r1, r2};
                Moments f(2);
                f.AddSimplex(seg, 2, std::sqrt(ComputeDistance2(r1, r2, 2)));
                acc.AddFacetPiece(i, j, f);
            }
        }
    });
}

// Power vertices of the triangles of pd (once, so that the cells sharing one
// use identical coordinates); false if a triangle is degenerate.
bool
_ComputePowerVertices(const PowerTriangulation& pd, std::vector<double>& powerVertices)
{
    const GEO::Delaunay& T = *pd._delaunay;
    powerVertices.resize(2 * size_t(T.nb_cells()));
    for (GEO::index_t t = 0; t < T.nb_cells(); ++t) {
        const double* p[3];
        for (GEO::index_t lv = 0; lv < 3; ++lv) {
            p[lv] = &pd._lifted[_stride * size_t(T.cell_vertex(t, lv))];
        }
        if (not _ComputePowerVertex(p[0], p[1], p[2], &powerVertices[2 * size_t(t)])) {
            return false;
        }
    }
    return true;
}

// Boundary loops of the union of the fragments of one cell, appended to out:
// fragments are welded at their shared vertices (which the clipping
// computes in each fragment, equal up to round-off), the edges shared by two
// fragments cancel, and the remaining ones are chained into loops (outer
// loops counter-clockwise, holes clockwise). Vertices where a loop goes
// straight through (from the domain triangle edges) are dropped.
void
_AppendCellLoops(int i, const std::vector<std::vector<_Vertex>>& fragments, PolyMesh& out)
{
    // Weld tolerance relative to the size of the cell.
    double lo[2] = {HUGE_VAL, HUGE_VAL}, hi[2] = {-HUGE_VAL, -HUGE_VAL};
    for (const std::vector<_Vertex>& a : fragments) {
        for (const _Vertex& V : a) {
            for (int c = 0; c < 2; ++c) {
                lo[c] = std::min(lo[c], V._x[c]);
                hi[c] = std::max(hi[c], V._x[c]);
            }
        }
    }
    const double tol = 1e-10 * std::max(hi[0] - lo[0], hi[1] - lo[1]);
    std::vector<double> pts;
    auto weld = [&](const double* x) {
        for (size_t k = 0; k < pts.size() / 2; ++k) {
            if (std::fabs(pts[2 * k] - x[0]) <= tol and std::fabs(pts[2 * k + 1] - x[1]) <= tol) {
                return int(k);
            }
        }
        pts.insert(pts.end(), x, x + 2);
        return int(pts.size() / 2 - 1);
    };
    std::map<std::pair<int, int>, int> count;  // directed edge -> multiplicity
    for (const std::vector<_Vertex>& a : fragments) {
        for (size_t v = 0; v < a.size(); ++v) {
            const int p = weld(a[v]._x), q = weld(a[(v + 1) % a.size()]._x);
            if (p == q) {
                continue;
            }
            const auto opposite = count.find({q, p});
            if (opposite != count.end()) {
                if (--opposite->second == 0) {
                    count.erase(opposite);
                }
            } else {
                ++count[{p, q}];
            }
        }
    }
    std::multimap<int, int> next;
    for (const auto& e : count) {
        for (int c = 0; c < e.second; ++c) {
            next.emplace(e.first.first, e.first.second);
        }
    }
    std::map<int, int> outIndex;  // welded point -> point of out
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
        RemoveCollinearVertices(loop, pts);
        if (loop.size() < 3) {
            continue;
        }
        // The points of the cell are shared by its loops (see PolyMesh).
        out._polySizes.push_back(int(loop.size()));
        out._polyCell.push_back(i);
        for (int v : loop) {
            auto it = outIndex.find(v);
            if (it == outIndex.end()) {
                it = outIndex.emplace(v, int(out.GetNumPoints())).first;
                out._points.insert(out._points.end(), &pts[2 * size_t(v)], &pts[2 * size_t(v)] + 2);
            }
            out._polyVerts.push_back(it->second);
        }
    }
}

} // namespace

bool
ComputePowerCells2D(
    const Domain& domain,
    const std::vector<double>& sites,
    const PowerTriangulation& pd,
    bool parallel,
    Accumulator& acc)
{
    std::vector<double> powerVertices;
    if (not _ComputePowerVertices(pd, powerVertices)) {
        return false;
    }
    const _DomainTriangles triangles(domain);
    return ForEachSiteChunk(acc.GetNumSites(), parallel, sites, 2, acc, [&](int first, int last, Accumulator& part) {
        return _IntegrateCells(pd, powerVertices, triangles, first, last, part);
    });
}

bool
ComputePowerPolygons2D(const Domain& domain, int numSites, const PowerTriangulation& pd, PolyMesh& out)
{
    std::vector<double> powerVertices;
    if (not _ComputePowerVertices(pd, powerVertices)) {
        return false;
    }
    const _DomainTriangles triangles(domain);
    return _ForEachCell(pd, powerVertices, triangles, 0, numSites,
                        [&](int i, const std::vector<std::vector<_Vertex>>& fragments) {
                            _AppendCellLoops(i, fragments, out);
                        });
}

} // namespace PP
