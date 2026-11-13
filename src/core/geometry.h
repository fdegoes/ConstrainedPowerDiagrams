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

// Signed area of a closed 2D polygon given by indices into pts.
inline double
ComputeLoopSignedArea(
    const std::vector<int>& loop, 
    const std::vector<double>& pts)
{
    double a = 0.0;
    for (size_t k = 0; k < loop.size(); ++k) {
        const double* p = &pts[2 * loop[k]];
        const double* q = &pts[2 * loop[(k + 1) % loop.size()]];
        a += p[0] * q[1] - p[1] * q[0];
    }
    return 0.5 * a;
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
