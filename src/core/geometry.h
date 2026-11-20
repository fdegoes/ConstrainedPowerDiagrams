#ifndef PP_CORE_GEOMETRY_H
#define PP_CORE_GEOMETRY_H

// Small geometric helpers on raw coordinate arrays (internal to the core
// library).

#include <algorithm>
#include <cmath>
#include <vector>

namespace PP
{

// Twice the signed area of the 2D triangle (a, b, c) (> 0 if CCW).
inline double
ComputeDoubleSignedArea(const double* a, const double* b, const double* c)
{
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

// Determinant of the 3x3 matrix with rows u, v, w.
inline double
ComputeDeterminant3(const double* u, const double* v, const double* w)
{
    return u[0] * (v[1] * w[2] - v[2] * w[1]) -
           u[1] * (v[0] * w[2] - v[2] * w[0]) +
           u[2] * (v[0] * w[1] - v[1] * w[0]);
}

// Six times the signed volume of the tet (a, b, c, d) (> 0 if positive).
inline double
ComputeSixSignedVolume(const double* a, const double* b, const double* c, const double* d)
{
    const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    const double w[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
    return ComputeDeterminant3(u, v, w);
}

// Cross product r = u x v.
inline void
ComputeCross(const double* u, const double* v, double* r)
{
    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];
}

inline double
ComputeDot(const double* u, const double* v, int dim)
{
    double s = 0.0;
    for (int c = 0; c < dim; ++c) {
        s += u[c] * v[c];
    }
    return s;
}

// Squared distance between two points.
inline double
ComputeDistance2(const double* a, const double* b, int dim)
{
    double s = 0.0;
    for (int c = 0; c < dim; ++c) {
        s += (a[c] - b[c]) * (a[c] - b[c]);
    }
    return s;
}

// Returns true if b is collinear with its neighbors a and c, and between
// them, or if an incident edge has zero length.
inline bool
IsStraightVertex(const double* a, const double* b, const double* c)
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

// Drops the vertices of a closed 2D loop (indices into pts) where it goes
// straight through, e.g. on the internal edges of a triangulated domain.
inline void
RemoveCollinearVertices(std::vector<int>& loop, const std::vector<double>& pts)
{
    bool changed = true;
    while (changed and loop.size() > 3) {
        changed = false;
        for (size_t k = 0; k < loop.size() and loop.size() > 3; ++k) {
            const size_t prev = (k + loop.size() - 1) % loop.size();
            const size_t next = (k + 1) % loop.size();
            if (IsStraightVertex(&pts[2 * loop[prev]], &pts[2 * loop[k]], &pts[2 * loop[next]])) {
                loop.erase(loop.begin() + k);
                changed = true;
                --k;
            }
        }
    }
}

// Drops the vertices of a closed 3D loop (indices into pts, 3 coordinates
// each) where it goes straight through, or that repeat their predecessor.
inline void
RemoveCollinearVertices3D(std::vector<int>& loop, const std::vector<double>& pts)
{
    auto straight = [&](int a, int b, int c) {
        const double* p = &pts[3 * size_t(a)];
        const double* q = &pts[3 * size_t(b)];
        const double* r = &pts[3 * size_t(c)];
        const double e1[3] = {q[0] - p[0], q[1] - p[1], q[2] - p[2]};
        const double e2[3] = {r[0] - q[0], r[1] - q[1], r[2] - q[2]};
        const double l1 = std::sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        const double l2 = std::sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
        if (l1 == 0.0 or l2 == 0.0) {
            return true;
        }
        const double x[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]};
        const double dot = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
        return dot > 0.0 and std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) <= 1e-10 * l1 * l2;
    };
    bool changed = true;
    while (changed and loop.size() > 3) {
        changed = false;
        for (size_t k = 0; k < loop.size() and loop.size() > 3; ++k) {
            const size_t prev = (k + loop.size() - 1) % loop.size();
            const size_t next = (k + 1) % loop.size();
            if (straight(loop[prev], loop[k], loop[next])) {
                loop.erase(loop.begin() + long(k));
                changed = true;
                --k;
            }
        }
    }
}

// Axis-aligned bounding box of dim-dimensional points.
struct BoundingBox
{
    double _lo[3] = {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    double _hi[3] = {-HUGE_VAL, -HUGE_VAL, -HUGE_VAL};

    void Add(const std::vector<double>& pts, int dim)
    {
        for (size_t i = 0; i + dim <= pts.size(); i += dim) 
        {
            for (int c = 0; c < dim; ++c) 
            {
                _lo[c] = std::min(_lo[c], pts[i + c]);
                _hi[c] = std::max(_hi[c], pts[i + c]);
            }
        }
    }

    double GetCenter(int c) const 
    { return 0.5 * (_lo[c] + _hi[c]); }
    
    double GetExtent(int c) const 
    { return _hi[c] - _lo[c]; }

    double GetDiagonal(int dim) const
    {
        double d2 = 0.0;
        for (int c = 0; c < dim; ++c) {
            d2 += GetExtent(c) * GetExtent(c);
        }
        return std::sqrt(d2);
    }
};

// True if every value is a finite number.
inline bool
AreAllFinite(const std::vector<double>& values)
{
    return std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); });
}

} // namespace PP

#endif // PP_CORE_GEOMETRY_H
