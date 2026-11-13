#include "houdini/tetrahedralizeSOP.h"

#include "houdini/geometryIO.h"

#include "core/tetrahedralize.h"

#include <GU/GU_Detail.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <PRM/PRM_Include.h>

#include <string>
#include <vector>

namespace PP
{

namespace
{

PRM_Name preprocessName("preprocess", "Preprocess Surface");
PRM_Name refineName("refine", "Refine");
PRM_Name qualityName("quality", "Quality");
PRM_Default qualityDefault(2.0);
PRM_Range qualityRange(PRM_RANGE_RESTRICTED, 1.0, PRM_RANGE_UI, 2.0);

} // namespace

PRM_Template TetrahedralizeSOP::myTemplateList[] = {
    PRM_Template(PRM_TOGGLE, 1, &preprocessName, PRMzeroDefaults),
    PRM_Template(PRM_TOGGLE, 1, &refineName, PRMzeroDefaults),
    PRM_Template(PRM_FLT, 1, &qualityName, &qualityDefault, nullptr, &qualityRange),
    PRM_Template(),
};

OP_Node*
TetrahedralizeSOP::myConstructor(OP_Network* net, const char* name, OP_Operator* op)
{
    return new TetrahedralizeSOP(net, name, op);
}

TetrahedralizeSOP::TetrahedralizeSOP(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{ }

TetrahedralizeSOP::~TetrahedralizeSOP() { }

const char*
TetrahedralizeSOP::inputLabel(OP_InputIdx idx) const
{
    return idx == 0 ? "Closed Triangulated Surface" : "Invalid";
}

bool
TetrahedralizeSOP::updateParmsFlags()
{
    bool changed = SOP_Node::updateParmsFlags();
    changed |= enableParm("quality", _REFINE() != 0);
    return changed;
}

TetrahedralizeOptions
TetrahedralizeSOP::_ReadOptions(fpreal t)
{
    TetrahedralizeOptions options;
    options._preprocess = _PREPROCESS() != 0;
    options._refine = _REFINE() != 0;
    options._quality = _QUALITY(t);
    return options;
}

OP_ERROR
TetrahedralizeSOP::cookMySop(OP_Context& context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    const fpreal t = context.getTime();
    gdp->clearAndDestroy();
    const GU_Detail* inGdp = inputGeo(0, context);

    static const int xyz[3] = {0, 1, 2};
    const std::vector<double> surfPts = ReadPositions(inGdp, xyz, 3);
    std::vector<int> surfTris;
    std::string err;
    if (not ReadTriangles(inGdp, surfTris, "append a Divide SOP", err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    std::vector<double> tetPts;
    std::vector<int> tets;
    if (not Tetrahedralize(surfPts, surfTris, _ReadOptions(t), tetPts, tets, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }
    AppendTetrahedra(gdp, tetPts, tets);
    return error();
}

} // namespace PP
