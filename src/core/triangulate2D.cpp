#include "core/triangulate2D.h"
#include "core/geogramContext.h"
#include "core/geometry.h"

#include <geogram/delaunay/CDT_2d.h>

#include <algorithm>
#include <utility>

namespace PP
{

namespace
{

bool
_CheckLoops(
    const std::vector<double>& loopPts, 
    const std::vector<int>& loopSizes, 
    std::string& err)
{
    size_t total = 0;
    for (int n : loopSizes) 
    {
        if (n < 3) 
        {
            err = "Every loop needs at least 3 vertices.";
            return false;
        }
        total += size_t(n);
    }

    if (loopSizes.empty() or total != loopPts.size() / 2) 
    {
        err = "Loop sizes do not match the number of loop vertices.";
        return false;
    }
    
    if (not AreAllFinite(loopPts)) 
    {
        err = "The loops have non-finite coordinates.";
        return false;
    }
    
    return true;
}

void
_InsertLoops(
    GEO::CDT2d& cdt, 
    const std::vector<double>& loopPts, 
    const std::vector<int>& loopSizes)
{
    BoundingBox box;
    box.Add(loopPts, 2);
    const double pad = 0.1 * std::max(box.GetExtent(0), box.GetExtent(1)) + 1e-6;
    cdt.create_enclosing_rectangle(
        box._lo[0] - pad, 
        box._lo[1] - pad, 
        box._hi[0] + pad, 
        box._hi[1] + pad);

    const GEO::index_t nbPts = GEO::index_t(loopPts.size() / 2);
    std::vector<GEO::index_t> ids(nbPts);
    cdt.insert(nbPts, loopPts.data(), ids.data());

    GEO::index_t first = 0;
    for (int n : loopSizes) 
    {
        for (int k = 0; k < n; ++k) 
        {
            const GEO::index_t a = ids[first + GEO::index_t(k)];
            const GEO::index_t b = ids[first + GEO::index_t((k + 1) % n)];
            if (a != b) {
                cdt.insert_constraint(a, b);
            }
        }
        first += GEO::index_t(n);
    }
}

void
_ExtractTriangles(
    const GEO::CDT2d& cdt, 
    std::vector<double>& triPts, std::vector<int>& tris)
{
    std::vector<int> cdtToOut(cdt.nv(), -1);
    auto GetOutVertex = [&](GEO::index_t v) 
    {
        if (cdtToOut[v] < 0) 
        {
            cdtToOut[v] = int(triPts.size() / 2);
            const GEO::vec2 p = cdt.point(v);
            triPts.push_back(p.x);
            triPts.push_back(p.y);
        }
        return cdtToOut[v];
    };

    for (GEO::index_t t = 0; t < cdt.nT(); ++t) 
    {
        int tv[3] = 
        {
            GetOutVertex(cdt.Tv(t, 0)), 
            GetOutVertex(cdt.Tv(t, 1)), 
            GetOutVertex(cdt.Tv(t, 2))
        };
        
        if (ComputeDoubleSignedArea(
            &triPts[2 * tv[0]], 
            &triPts[2 * tv[1]], 
            &triPts[2 * tv[2]]) < 0.0) 
        {
            std::swap(tv[1], tv[2]);
        }

        tris.insert(tris.end(), tv, tv + 3);
    }
}

} // namespace

bool
Triangulate2D(
    const std::vector<double>& loopPts,
    const std::vector<int>& loopSizes,
    std::vector<double>& triPts,
    std::vector<int>& tris,
    std::string& err)
{
    tris.clear();
    triPts.clear();

    if (not _CheckLoops(loopPts, loopSizes, err)) {
        return false;
    }

    GeogramLock lock;

    GEO::CDT2d cdt;
    _InsertLoops(cdt, loopPts, loopSizes);
    
    cdt.remove_external_triangles(true);
    if (cdt.nT() == 0) {
        err = "The loops do not enclose any area.";
        return false;
    }
    
    _ExtractTriangles(cdt, triPts, tris);
    return true;
}

} // namespace PP
