#include "houdini/geometryIO.h"

#include <GA/GA_Handle.h>
#include <GA/GA_Iterator.h>
#include <GA/GA_PolyCounts.h>
#include <GEO/GEO_PrimPoly.h>
#include <GEO/GEO_PrimTetrahedron.h>
#include <GU/GU_Detail.h>
#include <UT/UT_Vector3.h>

namespace PP
{

std::vector<double>
ReadPositions(const GU_Detail* gdp, const int* axes, int dim)
{
    std::vector<double> pts;
    pts.reserve(size_t(dim) * gdp->getNumPoints());
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        const UT_Vector3D p(gdp->getPos3(ptoff));
        for (int c = 0; c < dim; ++c) {
            pts.push_back(p[axes[c]]);
        }
    }
    return pts;
}

bool
ReadTriangles(const GU_Detail* gdp, std::vector<int>& tris, const char* hint, std::string& err)
{
    tris.clear();
    const GEO_Primitive* prim;
    GA_FOR_ALL_PRIMITIVES(gdp, prim)
    {
        if (prim->getTypeId() != GA_PRIMPOLY) {
            continue;
        }
        if (prim->getVertexCount() != 3) {
            err = std::string("The mesh must be triangulated (") + hint + ").";
            return false;
        }
        for (GA_Size i = 0; i < 3; ++i) {
            tris.push_back(int(gdp->pointIndex(prim->getPointOffset(i))));
        }
    }
    if (tris.empty()) {
        err = std::string("The input has no triangles (") + hint + ").";
        return false;
    }
    return true;
}

std::vector<int>
ReadTetrahedra(const GU_Detail* gdp)
{
    std::vector<int> tets;
    const GEO_Primitive* prim;
    GA_FOR_ALL_PRIMITIVES(gdp, prim)
    {
        if (prim->getTypeId() != GA_PRIMTETRAHEDRON) {
            continue;
        }
        for (GA_Size i = 0; i < 4; ++i) {
            tets.push_back(int(gdp->pointIndex(prim->getPointOffset(i))));
        }
    }
    return tets;
}

double
ComputeAverageHeight(const GU_Detail* gdp, int axis)
{
    double sum = 0.0;
    GA_Size count = 0;
    const GEO_Primitive* prim;
    GA_FOR_ALL_PRIMITIVES(gdp, prim)
    {
        if (prim->getTypeId() != GA_PRIMPOLY) {
            continue;
        }
        for (GA_Size i = 0; i < prim->getVertexCount(); ++i) {
            sum += prim->getPos3(i)[axis];
        }
        count += prim->getVertexCount();
    }
    return count ? sum / double(count) : 0.0;
}

void
ReadLoops(
    const GU_Detail* gdp,
    const int* axes,
    std::vector<double>& loopPts,
    std::vector<int>& loopSizes,
    double& height)
{
    loopPts.clear();
    loopSizes.clear();
    double sum = 0.0;
    const GEO_Primitive* prim;
    GA_FOR_ALL_PRIMITIVES(gdp, prim)
    {
        const GA_Size n = prim->getVertexCount();
        if (prim->getTypeId() != GA_PRIMPOLY or n < 3) {
            continue;
        }
        for (GA_Size i = 0; i < n; ++i) {
            const UT_Vector3D p(prim->getPos3(i));
            loopPts.push_back(p[axes[0]]);
            loopPts.push_back(p[axes[1]]);
            sum += p[axes[2]];
        }
        loopSizes.push_back(int(n));
    }
    height = loopPts.empty() ? 0.0 : sum / (0.5 * double(loopPts.size()));
}

std::vector<UT_Vector3>
MakePositions3D(const std::vector<double>& pts, int dim, const int* axes, double height)
{
    const size_t n = pts.size() / dim;
    std::vector<UT_Vector3> pos(n);
    for (size_t i = 0; i < n; ++i) {
        UT_Vector3D p(0.0, 0.0, 0.0);
        for (int c = 0; c < dim; ++c) {
            p[axes[c]] = pts[dim * i + c];
        }
        if (dim == 2) {
            p[axes[2]] = height;
        }
        pos[i] = UT_Vector3(p);
    }
    return pos;
}

GA_Offset
AppendPolygons(
    GU_Detail* gdp,
    const std::vector<UT_Vector3>& pos,
    const std::vector<int>& polySizes,
    const std::vector<int>& polyVerts)
{
    if (pos.empty() or polySizes.empty()) {
        return GA_INVALID_OFFSET;
    }
    GA_PolyCounts sizes;
    for (int n : polySizes) {
        sizes.append(n);
    }
    return GEO_PrimPoly::buildBlock(gdp, pos.data(), GA_Size(pos.size()), sizes, polyVerts.data(), true);
}

void
AppendTetrahedra(GU_Detail* gdp, const std::vector<double>& pts, const std::vector<int>& tets)
{
    const GA_Size nbPts = GA_Size(pts.size() / 3);
    const GA_Offset startPt = gdp->appendPointBlock(nbPts);
    for (GA_Size i = 0; i < nbPts; ++i) {
        gdp->setPos3(startPt + GA_Offset(i), UT_Vector3D(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]));
    }
    GEO_PrimTetrahedron::buildBlock(gdp, startPt, nbPts, GA_Size(tets.size() / 4), tets.data());
}

void
SetDetailInt(GU_Detail* gdp, const char* name, int value)
{
    GA_RWHandleI h(gdp->addIntTuple(GA_ATTRIB_DETAIL, name, 1));
    h.set(GA_Offset(0), value);
}

void
SetDetailFloat(GU_Detail* gdp, const char* name, double value)
{
    GA_RWHandleF h(gdp->addFloatTuple(GA_ATTRIB_DETAIL, name, 1));
    h.set(GA_Offset(0), float(value));
}

void
SetDetailString(GU_Detail* gdp, const char* name, const char* value)
{
    GA_RWHandleS h(gdp->addStringTuple(GA_ATTRIB_DETAIL, name, 1));
    h.set(GA_Offset(0), value);
}

} // namespace PP
