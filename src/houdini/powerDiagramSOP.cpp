#include "houdini/powerDiagramSOP.h"
#include "houdini/domainInput.h"
#include "houdini/geometryIO.h"
#include "houdini/plane.h"

#include "core/powerDiagram.h"

#include <GA/GA_Handle.h>
#include <GU/GU_Detail.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <PRM/PRM_Include.h>

#include <algorithm>
#include <string>
#include <vector>

namespace PP
{

namespace
{

PRM_Name weightAttribName("weightattrib", "Weight Attribute");
PRM_Default weightAttribDefault(0, "w");

PRM_Name shiftAttribName("shiftattrib", "Shift Attribute");
PRM_Default shiftAttribDefault(0, "u");

PRM_Name boundaryAngleName("boundaryangle", "Boundary Merge Angle");
PRM_Default boundaryAngleDefault(0.0);
PRM_Range boundaryAngleRange(PRM_RANGE_RESTRICTED, 0.0, PRM_RANGE_RESTRICTED, 180.0);

// Primitive attribute "cell": the site of each polygon.
void
_WriteCellAttribute(GU_Detail* gdp, GA_Offset firstPrim, const std::vector<int>& polyCell)
{
    if (polyCell.empty()) {
        return;
    }
    GA_RWHandleI cellH(gdp->addIntTuple(GA_ATTRIB_PRIMITIVE, "cell", 1));
    for (size_t i = 0; i < polyCell.size(); ++i) {
        cellH.set(firstPrim + GA_Offset(i), polyCell[i]);
    }
}

// Sites without any polygon (hidden by weights or outside the domain).
int
_CountEmptyCells(const std::vector<int>& polyCell, size_t nbSites)
{
    std::vector<char> hasCell(nbSites, 0);
    for (int c : polyCell) {
        if (c >= 0 and size_t(c) < nbSites) {
            hasCell[c] = 1;
        }
    }
    return int(std::count(hasCell.begin(), hasCell.end(), 0));
}

} // namespace

PRM_Template PowerDiagramSOP::myTemplateList[] = {
    PRM_Template(PRM_ORD, 1, &modeParmName, PRMzeroDefaults, &modeParmMenu),
    PRM_Template(PRM_ORD, 1, &planeParmName, PRMzeroDefaults, &planeParmMenu),
    PRM_Template(PRM_STRING, 1, &weightAttribName, &weightAttribDefault),
    PRM_Template(PRM_STRING, 1, &shiftAttribName, &shiftAttribDefault),
    PRM_Template(PRM_FLT, 1, &boundaryAngleName, &boundaryAngleDefault,
                 nullptr, &boundaryAngleRange),
    PRM_Template(),
};

OP_Node*
PowerDiagramSOP::myConstructor(OP_Network* net, const char* name, OP_Operator* op)
{
    return new PowerDiagramSOP(net, name, op);
}

PowerDiagramSOP::PowerDiagramSOP(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{ }

PowerDiagramSOP::~PowerDiagramSOP() { }

const char*
PowerDiagramSOP::inputLabel(OP_InputIdx idx) const
{
    switch (idx) {
    case 0: return "Sites (P, w, u)";
    case 1: return "Clipping Polytope";
    default: return "Invalid";
    }
}

bool
PowerDiagramSOP::updateParmsFlags()
{
    bool changed = SOP_Node::updateParmsFlags();
    const bool is2D = (_MODE() == MODE_2D);
    changed |= setVisibleState("plane", is2D);
    changed |= setVisibleState("boundaryangle", not is2D);
    return changed;
}

OP_ERROR
PowerDiagramSOP::cookMySop(OP_Context& context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    const fpreal t = context.getTime();
    gdp->clearAndDestroy();

    const GU_Detail* siteGdp = inputGeo(0, context);
    const GU_Detail* clipGdp = inputGeo(1, context);
    const bool is2D = (_MODE() == MODE_2D);
    const int* axes = GetPlaneAxes(is2D ? _PLANE() : PLANE_XY);

    // Sites: s = P + u, weight = w + |u|^2.
    UT_String weightName, shiftName;
    _WEIGHTATTRIB(weightName, t);
    _SHIFTATTRIB(shiftName, t);
    SiteInput in;
    ReadSites(siteGdp, is2D, axes, weightName.c_str(), shiftName.c_str(), in);
    if (not in._warning.empty()) {
        addWarning(SOP_MESSAGE, in._warning.c_str());
    }
    const std::vector<double>& sites = in._sites;
    const std::vector<double>& weights = in._weights;
    if (weights.empty()) {
        addError(SOP_MESSAGE, "The first input has no points.");
        return error();
    }

    std::vector<double> domPts;
    std::vector<int> elems;
    double planeOffset = 0.0;
    std::string err;
    if (not ReadDomain(clipGdp, is2D, axes, domPts, elems, planeOffset, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    PolyMesh cells;
    const bool ok = is2D
        ? ComputePowerDiagram2D(sites, weights, domPts, elems, cells, err)
        : ComputePowerDiagram3D(sites, weights, domPts, elems, _BOUNDARYANGLE(t), cells, err);
    if (not ok) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    const GA_Offset firstPrim =
        AppendPolygons(gdp, MakePositions3D(cells._points, cells._dim, axes, planeOffset),
                       cells._polySizes, cells._polyVerts);
    _WriteCellAttribute(gdp, firstPrim, cells._polyCell);
    SetDetailInt(gdp, "nb_empty_cells", _CountEmptyCells(cells._polyCell, weights.size()));

    return error();
}

} // namespace PP
