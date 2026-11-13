#include "houdini/triangulate2DSOP.h"
#include "houdini/geometryIO.h"
#include "houdini/plane.h"

#include "core/triangulate2D.h"

#include <GU/GU_Detail.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <PRM/PRM_Include.h>

#include <string>
#include <vector>

namespace PP
{

PRM_Template Triangulate2DSOP::myTemplateList[] = {
    PRM_Template(PRM_ORD, 1, &planeParmName, PRMzeroDefaults, &planeParmMenu),
    PRM_Template(),
};

OP_Node*
Triangulate2DSOP::myConstructor(OP_Network* net, const char* name, OP_Operator* op)
{
    return new Triangulate2DSOP(net, name, op);
}

Triangulate2DSOP::Triangulate2DSOP(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{ }

Triangulate2DSOP::~Triangulate2DSOP() { }

const char*
Triangulate2DSOP::inputLabel(OP_InputIdx idx) const
{
    return idx == 0 ? "Boundary Polygons" : "Invalid";
}

OP_ERROR
Triangulate2DSOP::cookMySop(OP_Context& context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    gdp->clearAndDestroy();
    const GU_Detail* inGdp = inputGeo(0, context);
    const int* axes = GetPlaneAxes(_PLANE());

    // Every polygon is one boundary loop.
    std::vector<double> loopPts;
    std::vector<int> loopSizes;
    double height = 0.0;
    ReadLoops(inGdp, axes, loopPts, loopSizes, height);
    if (loopSizes.empty()) {
        addError(SOP_MESSAGE, "The input has no polygon with at least 3 vertices.");
        return error();
    }

    std::vector<double> triPts;
    std::vector<int> tris;
    std::string err;
    if (not Triangulate2D(loopPts, loopSizes, triPts, tris, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    // Welded triangles at the average height of the input.
    AppendPolygons(gdp, MakePositions3D(triPts, 2, axes, height),
                   std::vector<int>(tris.size() / 3, 3), tris);
    return error();
}

} // namespace PP
