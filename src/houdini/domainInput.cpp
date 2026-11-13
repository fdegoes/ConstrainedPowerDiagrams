#include "houdini/domainInput.h"
#include "houdini/geometryIO.h"

#include "core/domain.h"

#include <GA/GA_Handle.h>
#include <GA/GA_Iterator.h>
#include <GEO/GEO_Primitive.h>
#include <GU/GU_Detail.h>
#include <PRM/PRM_ChoiceList.h>
#include <PRM/PRM_Name.h>
#include <UT/UT_Vector3.h>

namespace PP
{

namespace
{

PRM_Name modeChoices[] = {
    PRM_Name("2d", "2D (Triangle Mesh)"),
    PRM_Name("3d", "3D (Tet Mesh)"),
    PRM_Name(),
};

} // namespace

PRM_Name       modeParmName("mode", "Mode");
PRM_ChoiceList modeParmMenu(PRM_CHOICELIST_SINGLE, modeChoices);

void
ReadSites(
    const GU_Detail* gdp,
    bool is2D,
    const int* axes,
    const char* weightAttr,
    const char* shiftAttr,
    SiteInput& out)
{
    out = SiteInput();
    out._dim = is2D ? 2 : 3;

    GA_ROHandleF weightH;
    if (weightAttr and *weightAttr) {
        weightH = GA_ROHandleF(gdp, GA_ATTRIB_POINT, weightAttr);
        if (not weightH.isValid()) {
            out._warning += std::string("No float point attribute '") + weightAttr + "', using 0. ";
        }
    }
    GA_ROHandleV3 shiftH;
    if (shiftAttr and *shiftAttr) {
        shiftH = GA_ROHandleV3(gdp, GA_ATTRIB_POINT, shiftAttr);
        if (not shiftH.isValid()) {
            out._warning += std::string("No vector point attribute '") + shiftAttr + "', using 0. ";
        }
    }

    const GA_Size n = gdp->getNumPoints();
    out._offsets.reserve(n);
    out._sites.reserve(out._dim * n);
    out._weights.reserve(n);
    out._shiftSq.reserve(n);
    GA_Offset ptoff;
    GA_FOR_ALL_PTOFF(gdp, ptoff)
    {
        UT_Vector3D p(gdp->getPos3(ptoff));
        UT_Vector3D u(0.0, 0.0, 0.0);
        if (shiftH.isValid()) {
            u = UT_Vector3D(shiftH.get(ptoff));
        }
        if (is2D) {
            u[axes[2]] = 0.0;
        }  // only the in-plane part of the shift matters
        const double w = weightH.isValid() ? double(weightH.get(ptoff)) : 0.0;
        const UT_Vector3D s = p + u;
        for (int c = 0; c < out._dim; ++c) {
            out._sites.push_back(s[axes[c]]);
        }
        out._shiftSq.push_back(u.length2());
        out._weights.push_back(w + u.length2());
        out._offsets.push_back(ptoff);
    }
}

bool
ReadDomain(
    const GU_Detail* gdp,
    bool is2D,
    const int* axes,
    std::vector<double>& pts,
    std::vector<int>& elems,
    double& planeOffset,
    std::string& err)
{
    pts = ReadPositions(gdp, axes, is2D ? 2 : 3);
    planeOffset = 0.0;
    if (is2D) {
        planeOffset = ComputeAverageHeight(gdp, axes[2]);
        return ReadTriangles(gdp, elems, "append a PP Triangulate 2D SOP", err);
    }
    elems = ReadTetrahedra(gdp);
    if (elems.empty()) {
        err = "The second input has no tetrahedra (append a PP Tetrahedralize SOP).";
        return false;
    }
    return true;
}

bool
ReadDomain(const GU_Detail* gdp, bool is2D, const int* axes, Domain& domain, std::string& err)
{
    std::vector<double> pts;
    std::vector<int> elems;
    double planeOffset = 0.0;
    if (not ReadDomain(gdp, is2D, axes, pts, elems, planeOffset, err)) {
        return false;
    }
    return is2D ? Domain::Make2D(pts, elems, domain, err) : Domain::Make3D(pts, elems, domain, err);
}

} // namespace PP
