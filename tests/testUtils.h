#ifndef PP_TESTS_TESTUTILS_H
#define PP_TESTS_TESTUTILS_H

// Shared helpers of the standalone tests.

#include "core/domain.h"
#include "core/tetrahedralize.h"
#include "core/triangulate2D.h"

#include <cstdio>
#include <random>
#include <string>
#include <vector>

namespace PPTest
{

inline int&
GetFailures()
{
    static int failures = 0;
    return failures;
}

inline void
Check(bool ok, const std::string& what)
{
    std::printf("[%s] %s\n", ok ? " OK " : "FAIL", what.c_str());
    std::fflush(stdout);
    if (not ok) {
        ++GetFailures();
    }
}

inline int
Report()
{
    const int failures = GetFailures();
    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}

inline bool
MakeDomain2D(const std::vector<double>& loops, const std::vector<int>& sizes, PP::Domain& D)
{
    std::vector<double> pts;
    std::vector<int> tris;
    std::string err;
    return PP::Triangulate2D(loops, sizes, pts, tris, err) and PP::Domain::Make2D(pts, tris, D, err);
}

inline bool
MakeUnitSquare(PP::Domain& D)
{
    return MakeDomain2D({0, 0, 1, 0, 1, 1, 0, 1}, {4}, D);
}

// L = [0,2]x[0,1] U [0,1]x[1,2], area 3.
inline bool
MakeLShape(PP::Domain& D)
{
    return MakeDomain2D({0, 0, 2, 0, 2, 1, 1, 1, 1, 2, 0, 2}, {6}, D);
}

inline bool
MakeUnitCube(PP::Domain& D)
{
    std::vector<double> pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
                               0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    std::vector<int> tris = {0, 2, 1, 0, 3, 2,  4, 5, 6, 4, 6, 7,
                             0, 1, 5, 0, 5, 4,  2, 3, 7, 2, 7, 6,
                             1, 2, 6, 1, 6, 5,  0, 4, 7, 0, 7, 3};
    std::vector<double> tetPts;
    std::vector<int> tets;
    std::string err;
    return PP::Tetrahedralize(pts, tris, PP::TetrahedralizeOptions(), tetPts, tets, err) and
           PP::Domain::Make3D(tetPts, tets, D, err);
}

// Short scientific notation, e.g. 1.2e-08.
inline std::string
FormatValue(double v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.2g", v);
    return buf;
}

inline std::vector<double>
RandomValues(size_t count, double lo, double hi, unsigned seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(lo, hi);
    std::vector<double> s(count);
    for (double& x : s) {
        x = U(rng);
    }
    return s;
}

inline std::vector<double>
RandomSites(int n, int dim, double lo, double hi, unsigned seed)
{
    return RandomValues(size_t(n) * dim, lo, hi, seed);
}

// Random positive fractions summing to 1 (between 0.5/n and 1.5/n).
inline std::vector<double>
RandomFractions(int n, unsigned seed)
{
    std::vector<double> f = RandomValues(n, 0.5, 1.5, seed);
    double sum = 0.0;
    for (double x : f) {
        sum += x;
    }
    for (double& x : f) {
        x /= sum;
    }
    return f;
}

inline std::vector<double>
ScaleFractions(const std::vector<double>& f, double total)
{
    std::vector<double> v(f);
    for (double& x : v) {
        x *= total;
    }
    return v;
}

} // namespace PPTest

#endif // PP_TESTS_TESTUTILS_H
