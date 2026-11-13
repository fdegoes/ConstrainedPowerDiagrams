// Standalone checks of the clipped power diagram (no Houdini needed).

#include "core/powerDiagram.h"
#include "core/tetrahedralize.h"
#include "core/triangulate2D.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace
{

int failures = 0;

void
_Check(bool ok, const char* what)
{
    std::printf("[%s] %s\n", ok ? " OK " : "FAIL", what);
    std::fflush(stdout);
    if (not ok) {
        ++failures;
    }
}

double
_ComputePolygonArea(const PP::PolyMesh& m, size_t poly, size_t first)
{
    double a = 0.0;
    int n = m._polySizes[poly];
    for (int k = 0; k < n; ++k) {
        int i = m._polyVerts[first + k];
        int j = m._polyVerts[first + (k + 1) % n];
        a += m._points[2 * i] * m._points[2 * j + 1] - m._points[2 * j] * m._points[2 * i + 1];
    }
    return 0.5 * a;
}

std::vector<double>
_ComputeCellAreas2D(const PP::PolyMesh& m, size_t nbSites)
{
    std::vector<double> areas(nbSites, 0.0);
    size_t first = 0;
    for (size_t p = 0; p < m.GetNumPolygons(); ++p) {
        areas[m._polyCell[p]] += _ComputePolygonArea(m, p, first);
        first += m._polySizes[p];
    }
    return areas;
}

// Signed volume of closed polygonal shells (fan triangulation, divergence).
std::vector<double>
_ComputeCellVolumes3D(const PP::PolyMesh& m, size_t nbSites)
{
    std::vector<double> vol(nbSites, 0.0);
    size_t first = 0;
    for (size_t p = 0; p < m.GetNumPolygons(); ++p) {
        int n = m._polySizes[p];
        const double* a = &m._points[3 * m._polyVerts[first]];
        for (int k = 1; k + 1 < n; ++k) {
            const double* b = &m._points[3 * m._polyVerts[first + k]];
            const double* c = &m._points[3 * m._polyVerts[first + k + 1]];
            double det = a[0] * (b[1] * c[2] - b[2] * c[1]) -
                         a[1] * (b[0] * c[2] - b[2] * c[0]) +
                         a[2] * (b[0] * c[1] - b[1] * c[0]);
            vol[m._polyCell[p]] += det / 6.0;
        }
        first += n;
    }
    return vol;
}

void
_MakeUnitCube(std::vector<double>& pts, std::vector<int>& tris)
{
    pts = {0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0,
           0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1};
    // Outward-oriented triangles.
    tris = {0, 2, 1, 0, 3, 2,  4, 5, 6, 4, 6, 7,
            0, 1, 5, 0, 5, 4,  2, 3, 7, 2, 7, 6,
            1, 2, 6, 1, 6, 5,  0, 4, 7, 0, 7, 3};
}

// Triangulates the polygon loops and computes the 2D power diagram.
bool
_ComputeDiagram2D(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& loopPts,
    const std::vector<int>& loopSizes,
    PP::PolyMesh& m,
    std::string& err)
{
    std::vector<double> triPts;
    std::vector<int> tris;
    return PP::Triangulate2D(loopPts, loopSizes, triPts, tris, err) and
           PP::ComputePowerDiagram2D(sites, weights, triPts, tris, m, err);
}

// Tetrahedralizes the surface and computes the 3D power diagram.
bool
_ComputeDiagram3D(
    const std::vector<double>& sites,
    const std::vector<double>& weights,
    const std::vector<double>& surfPts,
    const std::vector<int>& surfTris,
    PP::PolyMesh& m,
    std::string& err)
{
    std::vector<double> tetPts;
    std::vector<int> tets;
    return PP::Tetrahedralize(surfPts, surfTris, PP::TetrahedralizeOptions(),
                              tetPts, tets, err) and
           PP::ComputePowerDiagram3D(sites, weights, tetPts, tets, 0.0, m, err);
}

void
_TestTriangulate2D()
{
    // Square with a square hole: area 1 - 0.25.
    std::vector<double> loops = {0, 0, 1, 0, 1, 1, 0, 1,
                                 0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75};
    std::vector<double> triPts;
    std::vector<int> tris;
    std::string err;
    bool ok = PP::Triangulate2D(loops, {4, 4}, triPts, tris, err);
    _Check(ok, ("triangulate2D square with hole " + err).c_str());
    if (not ok) {
        return;
    }
    double area = 0.0;
    bool ccw = true;
    for (size_t t = 0; t < tris.size(); t += 3) {
        const double* a = &triPts[2 * tris[t]];
        const double* b = &triPts[2 * tris[t + 1]];
        const double* c = &triPts[2 * tris[t + 2]];
        double a2 = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
        ccw = ccw and a2 > 0.0;
        area += 0.5 * a2;
    }
    _Check(std::fabs(area - 0.75) < 1e-12, "triangulate2D even-odd rule leaves the hole out");
    _Check(ccw, "triangulate2D outputs counter-clockwise triangles");
    _Check(triPts.size() == 16, "triangulate2D adds no points");
}

void
_Test2D()
{
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const size_t n = 50;
    std::vector<double> sites, weights;
    for (size_t i = 0; i < n; ++i) {
        sites.push_back(U(rng));
        sites.push_back(U(rng));
        weights.push_back(0.0);
    }
    std::vector<double> square = {0, 0, 1, 0, 1, 1, 0, 1};

    PP::PolyMesh m;
    std::string err;
    bool ok = _ComputeDiagram2D(sites, weights, square, {4}, m, err);
    _Check(ok, ("2D Voronoi in unit square " + err).c_str());
    if (not ok) {
        return;
    }
    auto areas = _ComputeCellAreas2D(m, n);
    double total = 0.0;
    bool allPositive = true;
    for (double a : areas) {
        total += a;
        allPositive = allPositive and a > 0.0;
    }
    _Check(std::fabs(total - 1.0) < 1e-9, "2D areas sum to the square area");
    _Check(allPositive, "2D every unweighted site has a cell");
    _Check(m.GetNumPolygons() == n, "2D one polygon per cell in a convex domain");

    // Increasing one weight grows its cell.
    weights[0] = 0.05;
    PP::PolyMesh m2;
    _ComputeDiagram2D(sites, weights, square, {4}, m2, err);
    auto areas2 = _ComputeCellAreas2D(m2, n);
    _Check(areas2[0] > areas[0], "2D heavier weight gives a larger cell");

    // Clockwise triangles give the same result.
    std::vector<double> triPts;
    std::vector<int> tris;
    PP::Triangulate2D(square, {4}, triPts, tris, err);
    for (size_t t = 0; t < tris.size(); t += 3) {
        std::swap(tris[t + 1], tris[t + 2]);
    }
    PP::PolyMesh m3;
    PP::ComputePowerDiagram2D(sites, weights, triPts, tris, m3, err);
    auto areas3 = _ComputeCellAreas2D(m3, n);
    double diff = 0.0;
    for (size_t i = 0; i < n; ++i) {
        diff += std::fabs(areas3[i] - areas2[i]);
    }
    _Check(diff < 1e-9, "2D result is independent of triangle orientation");
}

void
_Test2DNonConvex()
{
    // U-shaped domain: the top site's cell is split by the notch.
    std::vector<double> ushape = {0, 0, 3, 0, 3, 3, 2, 3, 2, 1, 1, 1, 1, 3, 0, 3};
    std::vector<double> sites = {1.5, 2.5,   1.5, 0.5};
    std::vector<double> weights = {0.0, 0.0};

    PP::PolyMesh m;
    std::string err;
    bool ok = _ComputeDiagram2D(sites, weights, ushape, {8}, m, err);
    _Check(ok, ("2D U-shaped domain " + err).c_str());
    if (not ok) {
        return;
    }
    int piecesOfSite0 = 0;
    for (int c : m._polyCell) {
        piecesOfSite0 += (c == 0);
    }
    _Check(piecesOfSite0 == 2, "2D non-convex domain splits a cell in two pieces");
    auto areas = _ComputeCellAreas2D(m, 2);
    _Check(std::fabs(areas[0] + areas[1] - 7.0) < 1e-9, "2D U-shape areas sum to 7");

    // A single site in a square with a hole: one outer loop (CCW) and one
    // hole loop (CW), whose signed areas add up to the domain area.
    std::vector<double> holed = {0, 0, 1, 0, 1, 1, 0, 1,
                                 0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75};
    PP::PolyMesh mh;
    ok = _ComputeDiagram2D({0.1, 0.1}, {0.0}, holed, {4, 4}, mh, err);
    _Check(ok and mh.GetNumPolygons() == 2, "2D cell around a hole has two boundary loops");
    if (ok) {
        _Check(std::fabs(_ComputeCellAreas2D(mh, 1)[0] - 0.75) < 1e-9, "2D holed cell area is 0.75");
    }
}

void
_Test3D()
{
    std::vector<double> pts;
    std::vector<int> tris;
    _MakeUnitCube(pts, tris);

    std::mt19937 rng(7);
    std::uniform_real_distribution<double> U(0.05, 0.95);
    const size_t n = 30;
    std::vector<double> sites, weights;
    for (size_t i = 0; i < n; ++i) {
        for (int c = 0; c < 3; ++c) {
            sites.push_back(U(rng));
        }
        weights.push_back(0.0);
    }

    PP::PolyMesh m;
    std::string err;
    bool ok = _ComputeDiagram3D(sites, weights, pts, tris, m, err);
    _Check(ok, ("3D Voronoi in unit cube " + err).c_str());
    if (not ok) {
        return;
    }
    auto vol = _ComputeCellVolumes3D(m, n);
    double total = 0.0;
    bool allPositive = true;
    for (double v : vol) {
        total += v;
        allPositive = allPositive and v > 0.0;
    }
    std::printf("       3D total volume = %.12f, polygons = %zu\n", total, m.GetNumPolygons());
    _Check(std::fabs(total - 1.0) < 1e-6, "3D volumes sum to the cube volume");
    _Check(allPositive, "3D every unweighted site has a cell (outward orientation)");

    weights[0] = 0.05;
    PP::PolyMesh m2;
    _ComputeDiagram3D(sites, weights, pts, tris, m2, err);
    auto vol2 = _ComputeCellVolumes3D(m2, n);
    _Check(vol2[0] > vol[0], "3D heavier weight gives a larger cell");

    // A very light site gets hidden.
    weights[0] = 0.0;
    weights[1] = -10.0;
    PP::PolyMesh m3;
    _ComputeDiagram3D(sites, weights, pts, tris, m3, err);
    std::set<int> used(m3._polyCell.begin(), m3._polyCell.end());
    _Check(used.count(1) == 0, "3D a site with very low weight has an empty cell");
    weights[1] = 0.0;

    // Inverted tets give the same, outward-oriented cells.
    std::vector<double> tetPts;
    std::vector<int> tets;
    PP::Tetrahedralize(pts, tris, PP::TetrahedralizeOptions(), tetPts, tets, err);
    for (size_t t = 0; t < tets.size(); t += 4) {
        std::swap(tets[t], tets[t + 1]);
    }
    PP::PolyMesh m4;
    PP::ComputePowerDiagram3D(sites, weights, tetPts, tets, 0.0, m4, err);
    auto vol4 = _ComputeCellVolumes3D(m4, n);
    double diff = 0.0;
    for (size_t i = 0; i < n; ++i) {
        diff += std::fabs(vol4[i] - vol[i]);
    }
    _Check(diff < 1e-9, "3D result is independent of tet orientation");
}

void
_Test3DCavity()
{
    // Unit cube with a cubic cavity [0.4, 0.6]^3: volume 1 - 0.008.
    std::vector<double> pts, inner;
    std::vector<int> tris, innerTris;
    _MakeUnitCube(pts, tris);
    _MakeUnitCube(inner, innerTris);
    for (double& x : inner) {
        x = 0.4 + 0.2 * x;
    }
    for (int& v : innerTris) {
        v += 8;
    }
    pts.insert(pts.end(), inner.begin(), inner.end());
    tris.insert(tris.end(), innerTris.begin(), innerTris.end());

    std::vector<double> sites = {0.2, 0.2, 0.2, 0.8, 0.8, 0.8, 0.2, 0.8, 0.5,
                                 0.8, 0.2, 0.5, 0.5, 0.5, 0.1};
    std::vector<double> weights(5, 0.0);
    PP::PolyMesh m;
    std::string err;
    bool ok = _ComputeDiagram3D(sites, weights, pts, tris, m, err);
    _Check(ok, ("3D cube with cavity " + err).c_str());
    if (not ok) {
        return;
    }
    double total = 0.0;
    for (double v : _ComputeCellVolumes3D(m, 5)) {
        total += v;
    }
    _Check(std::fabs(total - 0.992) < 1e-9, "3D cavity volumes sum to 0.992");
}

} // namespace

int
main()
{
    _TestTriangulate2D();
    _Test2D();
    _Test2DNonConvex();
    _Test3D();
    _Test3DCavity();
    std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures,
                failures == 1 ? "" : "s");
    return failures ? 1 : 0;
}
