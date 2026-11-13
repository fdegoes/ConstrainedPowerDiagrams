#include "houdini/feasibleVolumesSOP.h"
#include "houdini/domainInput.h"
#include "houdini/geometryIO.h"
#include "houdini/newtonParms.h"
#include "houdini/plane.h"

#include "core/domain.h"
#include "core/feasibleVolumes.h"

#include <GA/GA_Handle.h>
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

PRM_Name outAttribName("outattrib", "Output Attribute");
PRM_Default outAttribDefault(0, "volume");

PRM_Name toleranceName("tolerance", "Tolerance");
PRM_Default toleranceDefault(1e-10);

PRM_Name switcherName("tabs", "Tabs");

} // namespace

PRM_Template*
FeasibleVolumesSOP::GetTemplates()
{
    static std::vector<PRM_Template> templates;
    static PRM_Default tabs[] = {
        PRM_Default(3, "Main"),
        PRM_Default(kNumNewtonParms + 1, "Newton Solver"),
    };
    if (templates.empty()) {
        templates.push_back(PRM_Template(PRM_SWITCHER, 2, &switcherName, tabs));
        // Main (3)
        templates.push_back(PRM_Template(PRM_ORD, 1, &modeParmName, PRMzeroDefaults, &modeParmMenu));
        templates.push_back(PRM_Template(PRM_ORD, 1, &planeParmName, PRMzeroDefaults, &planeParmMenu));
        templates.push_back(PRM_Template(PRM_STRING, 1, &outAttribName, &outAttribDefault));
        // Newton Solver (kNumNewtonParms + 1)
        templates.push_back(PRM_Template(PRM_FLT_LOG, 1, &toleranceName, &toleranceDefault));
        AppendNewtonParmTemplates(templates);
        templates.push_back(PRM_Template());
    }
    return templates.data();
}

OP_Node*
FeasibleVolumesSOP::myConstructor(OP_Network* net, const char* name, OP_Operator* op)
{
    return new FeasibleVolumesSOP(net, name, op);
}

FeasibleVolumesSOP::FeasibleVolumesSOP(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{ }

FeasibleVolumesSOP::~FeasibleVolumesSOP() { }

const char*
FeasibleVolumesSOP::inputLabel(OP_InputIdx idx) const
{
    switch (idx) {
    case 0: return "Sites (target centroids)";
    case 1: return "Clipping Domain";
    default: return "Invalid";
    }
}

bool
FeasibleVolumesSOP::updateParmsFlags()
{
    bool changed = SOP_Node::updateParmsFlags();
    changed |= setVisibleState("plane", _MODE() == MODE_2D);
    return changed;
}

bool
FeasibleVolumesSOP::_WriteFractions(
    const SiteInput& in,
    const std::vector<double>& fractions,
    const char* attrib)
{
    GA_RWHandleF outH(gdp->addFloatTuple(GA_ATTRIB_POINT, attrib, 1));
    if (not outH.isValid()) {
        addError(SOP_MESSAGE, "Cannot create the output attribute.");
        return false;
    }
    for (size_t i = 0; i < in._offsets.size(); ++i) {
        outH.set(in._offsets[i], float(fractions[i]));
    }
    outH.bumpDataId();
    return true;
}

OP_ERROR
FeasibleVolumesSOP::cookMySop(OP_Context& context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    const fpreal t = context.getTime();
    duplicateSource(0, context);
    const bool is2D = (_MODE() == MODE_2D);
    const int* axes = GetPlaneAxes(is2D ? _PLANE() : PLANE_XY);

    // Targets: point positions, projected on the plane in 2D.
    SiteInput in;
    ReadSites(gdp, is2D, axes, nullptr, nullptr, in);
    const size_t n = in._offsets.size();
    if (n == 0) {
        addError(SOP_MESSAGE, "The first input has no points.");
        return error();
    }

    Domain domain;
    std::string err;
    if (not ReadDomain(inputGeo(1, context), is2D, axes, domain, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    const NewtonParms newton = ReadNewtonParameters(this, t);
    FeasibleVolumeSolver::Parameters params;
    params._maxIterations = newton._maxIterations;
    params._regularization = newton._regularization;
    params._maxHalvings = newton._maxHalvings;
    params._tolerance = _TOLERANCE(t);
    FeasibleVolumeSolver solver(params);
    std::vector<double> fractions;
    const FeasibleVolumeSolver::Result res =
        solver.Solve(in._sites, domain.GetCentroid(), in._dim, fractions);
    if (not res._feasible) {
        addError(SOP_MESSAGE, res._error.c_str());
        return error();
    }

    UT_String outName;
    _OUTATTRIB(outName, t);
    if (not _WriteFractions(in, fractions, outName.c_str())) {
        return error();
    }
    SetDetailFloat(gdp, "domain_volume", domain.GetVolume());
    return error();
}

} // namespace PP
