#include "core/domainBoundary.h"
#include "core/cellBackend.h"
#include "core/domain.h"
#include "core/geogramContext.h"

#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace PP
{

namespace
{

double
_Dot(const double* a, const double* b, int dim)
{
    double s = 0.0;
    for (int c = 0; c < dim; ++c) {
        s += a[c] * b[c];
    }
    return s;
}

// Closest point to x on segment ab (dim coordinates).
void
_ClosestOnSegment(const double* x, const double* a, const double* b, int dim, double* out)
{
    double ab[3], ax[3];
    for (int c = 0; c < dim; ++c) {
        ab[c] = b[c] - a[c];
        ax[c] = x[c] - a[c];
    }
    const double l2 = _Dot(ab, ab, dim);
    const double t = l2 > 0.0 ? std::clamp(_Dot(ax, ab, dim) / l2, 0.0, 1.0) : 0.0;
    for (int c = 0; c < dim; ++c) {
        out[c] = a[c] + t * ab[c];
    }
}

// Closest point to x on triangle abc in 3D (Ericson, Real-Time Collision
// Detection, 5.1.5).
void
_ClosestOnTriangle(const double* p, const double* a, const double* b, const double* c, double* out)
{
    double ab[3], ac[3], ap[3];
    for (int k = 0; k < 3; ++k) {
        ab[k] = b[k] - a[k];
        ac[k] = c[k] - a[k];
        ap[k] = p[k] - a[k];
    }
    const double d1 = _Dot(ab, ap, 3), d2 = _Dot(ac, ap, 3);
    auto set = [&](const double* q) { std::copy(q, q + 3, out); };
    auto lerp = [&](const double* q, const double* e, double t) {
        for (int k = 0; k < 3; ++k) {
            out[k] = q[k] + t * e[k];
        }
    };
    if (d1 <= 0.0 and d2 <= 0.0) {
        return set(a);
    }
    double bp[3];
    for (int k = 0; k < 3; ++k) {
        bp[k] = p[k] - b[k];
    }
    const double d3 = _Dot(ab, bp, 3), d4 = _Dot(ac, bp, 3);
    if (d3 >= 0.0 and d4 <= d3) {
        return set(b);
    }
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 and d1 >= 0.0 and d3 <= 0.0) {
        return lerp(a, ab, d1 / (d1 - d3));
    }
    double cp[3];
    for (int k = 0; k < 3; ++k) {
        cp[k] = p[k] - c[k];
    }
    const double d5 = _Dot(ab, cp, 3), d6 = _Dot(ac, cp, 3);
    if (d6 >= 0.0 and d5 <= d6) {
        return set(c);
    }
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 and d2 >= 0.0 and d6 <= 0.0) {
        return lerp(a, ac, d2 / (d2 - d6));
    }
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 and (d4 - d3) >= 0.0 and (d5 - d6) >= 0.0) {
        double bc[3];
        for (int k = 0; k < 3; ++k) {
            bc[k] = c[k] - b[k];
        }
        return lerp(b, bc, (d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    for (int k = 0; k < 3; ++k) {
        out[k] = a[k] + ab[k] * v + ac[k] * w;
    }
}

} // namespace

DomainBoundary::DomainBoundary(const Domain& domain)
{
    if (not domain.IsValid()) {
        return;
    }
    _dim = domain.GetDimension();
    const CellBackend& backend = domain.GetBackend();
    GeogramLock lock;
    const GEO::Mesh& M = domain.GetMesh();
    const GEO::index_t nv = M.vertices.nb();
    _points.resize(size_t(_dim) * nv);
    for (GEO::index_t v = 0; v < nv; ++v) {
        std::copy(M.vertices.point_ptr(v), M.vertices.point_ptr(v) + _dim, &_points[size_t(_dim) * v]);
    }
    const int perElem = backend.GetElementSize();
    const GEO::index_t ne = backend.GetNumElements(M);
    for (GEO::index_t e = 0; e < ne; ++e) {
        for (int lv = 0; lv < perElem; ++lv) {
            _elements.push_back(int(backend.GetElementVertex(M, e, GEO::index_t(lv))));
        }
        // Faces with no neighbor. 2D: edge k of a facet is (k, k + 1);
        // 3D: geogram's local facet k of a tet.
        for (GEO::index_t k = 0; k < backend.GetNumAdjacent(M, e); ++k) {
            if (backend.GetAdjacent(M, e, k) != GEO::NO_INDEX) {
                continue;
            }
            if (_dim == 2) {
                _faces.push_back(int(M.facets.vertex(e, k)));
                _faces.push_back(int(M.facets.vertex(e, (k + 1) % 3)));
            } else {
                for (GEO::index_t j = 0; j < 3; ++j) {
                    _faces.push_back(int(M.cells.facet_vertex(e, k, j)));
                }
            }
        }
    }
}

bool
DomainBoundary::_IsInside(const double* x) const
{
    const size_t perElem = size_t(_dim + 1);
    for (size_t e = 0; e < _elements.size(); e += perElem) {
        // Barycentric coordinates by the signed measures of the sub-simplices.
        const double* v0 = &_points[size_t(_dim) * size_t(_elements[e])];
        if (_dim == 2) {
            const double* v1 = &_points[2 * size_t(_elements[e + 1])];
            const double* v2 = &_points[2 * size_t(_elements[e + 2])];
            auto area = [](const double* a, const double* b, const double* c) {
                return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            };
            const double A = area(v0, v1, v2);
            const double a0 = area(x, v1, v2) / A, a1 = area(v0, x, v2) / A, a2 = area(v0, v1, x) / A;
            if (a0 >= 0.0 and a1 >= 0.0 and a2 >= 0.0) {
                return true;
            }
        } else {
            const double* v[4] = {v0, &_points[3 * size_t(_elements[e + 1])], &_points[3 * size_t(_elements[e + 2])],
                                  &_points[3 * size_t(_elements[e + 3])]};
            auto vol = [](const double* a, const double* b, const double* c, const double* d) {
                const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
                const double w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
                const double z[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
                return u[0] * (w[1] * z[2] - w[2] * z[1]) - u[1] * (w[0] * z[2] - w[2] * z[0]) +
                       u[2] * (w[0] * z[1] - w[1] * z[0]);
            };
            const double V = vol(v[0], v[1], v[2], v[3]);
            bool in = true;
            for (int k = 0; k < 4 and in; ++k) {
                const double* s[4] = {v[0], v[1], v[2], v[3]};
                s[k] = x;
                in = vol(s[0], s[1], s[2], s[3]) / V >= 0.0;
            }
            if (in) {
                return true;
            }
        }
    }
    return false;
}

double
DomainBoundary::_FaceDistance2(size_t f, const double* x, double* closest) const
{
    const int* v = &_faces[size_t(_dim) * f];
    if (_dim == 2) {
        _ClosestOnSegment(x, &_points[2 * size_t(v[0])], &_points[2 * size_t(v[1])], 2, closest);
    } else {
        _ClosestOnTriangle(x, &_points[3 * size_t(v[0])], &_points[3 * size_t(v[1])], &_points[3 * size_t(v[2])],
                           closest);
    }
    double d2 = 0.0;
    for (int c = 0; c < _dim; ++c) {
        d2 += (x[c] - closest[c]) * (x[c] - closest[c]);
    }
    return d2;
}

void
DomainBoundary::Project(const double* x, double* closest) const
{
    double best = std::numeric_limits<double>::infinity();
    double y[3];
    for (size_t f = 0; f < GetNumFaces(); ++f) {
        const double d2 = _FaceDistance2(f, x, y);
        if (d2 < best) {
            best = d2;
            std::copy(y, y + _dim, closest);
        }
    }
}

double
DomainBoundary::Evaluate(const double* x) const
{
    double closest[3];
    Project(x, closest);
    double d2 = 0.0;
    for (int c = 0; c < _dim; ++c) {
        d2 += (x[c] - closest[c]) * (x[c] - closest[c]);
    }
    const double d = std::sqrt(d2);
    return _IsInside(x) ? -d : d;
}

} // namespace PP
