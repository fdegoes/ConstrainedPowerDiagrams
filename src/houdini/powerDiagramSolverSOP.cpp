#include "houdini/powerDiagramSolverSOP.h"
#include "houdini/domainInput.h"
#include "houdini/geometryIO.h"
#include "houdini/newtonParms.h"
#include "houdini/plane.h"

#include "core/domain.h"
#include "core/dualSolver.h"

#include <GA/GA_Handle.h>
#include <GU/GU_Detail.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <PRM/PRM_Include.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

namespace PP
{

namespace
{

PRM_Name constraintName("constraint", "Constraint");
PRM_Name constraintChoices[] = {
    PRM_Name("volume", "Volume"),
    PRM_Name(),
};
PRM_ChoiceList constraintMenu(PRM_CHOICELIST_SINGLE, constraintChoices);

PRM_Name targetAttribName("targetattrib", "Target Volume Attribute");
PRM_Default targetAttribDefault(0, "volume");

PRM_Name weightAttribName("weightattrib", "Initial Weight Attribute");
PRM_Default weightAttribDefault(0, "w");

PRM_Name shiftAttribName("shiftattrib", "Shift Attribute");
PRM_Default shiftAttribDefault(0, "u");

PRM_Name outWeightAttribName("outweightattrib", "Output Weight Attribute");
PRM_Default outWeightAttribDefault(0, "w");

PRM_Name volumeTolName("volumetolerance", "Volume Tolerance");
PRM_Default volumeTolDefault(1e-8);


PRM_Name switcherName("tabs", "Tabs");

} // namespace

PRM_Template*
PowerDiagramSolverSOP::GetTemplates()
{
    static std::vector<PRM_Template> templates;
    static PRM_Default tabs[] = {
        PRM_Default(9, "Main"),
        PRM_Default(kNumNewtonParms + 1, "Newton Solver"),
    };
    if (templates.empty()) {
        templates.push_back(PRM_Template(PRM_SWITCHER, 2, &switcherName, tabs));
        // Main (9)
        templates.push_back(PRM_Template(PRM_ORD, 1, &modeParmName, PRMzeroDefaults, &modeParmMenu));
        templates.push_back(PRM_Template(PRM_ORD, 1, &planeParmName, PRMzeroDefaults, &planeParmMenu));
        templates.push_back(PRM_Template(PRM_ORD, 1, &constraintName, PRMzeroDefaults, &constraintMenu));
        templates.push_back(PRM_Template(PRM_SEPARATOR));
        templates.push_back(PRM_Template(PRM_STRING, 1, &targetAttribName, &targetAttribDefault));
        templates.push_back(PRM_Template(PRM_STRING, 1, &weightAttribName, &weightAttribDefault));
        templates.push_back(PRM_Template(PRM_STRING, 1, &shiftAttribName, &shiftAttribDefault));
        templates.push_back(PRM_Template(PRM_SEPARATOR));
        templates.push_back(PRM_Template(PRM_STRING, 1, &outWeightAttribName, &outWeightAttribDefault));
        // Newton Solver (kNumNewtonParms + 1)
        templates.push_back(PRM_Template(PRM_FLT_LOG, 1, &volumeTolName, &volumeTolDefault));
        AppendNewtonParmTemplates(templates);
        templates.push_back(PRM_Template());
    }
    return templates.data();
}

OP_Node*
PowerDiagramSolverSOP::myConstructor(OP_Network* net, const char* name, OP_Operator* op)
{
    return new PowerDiagramSolverSOP(net, name, op);
}

PowerDiagramSolverSOP::PowerDiagramSolverSOP(OP_Network* net, const char* name, OP_Operator* op)
    : SOP_Node(net, name, op)
{ }

PowerDiagramSolverSOP::~PowerDiagramSolverSOP() { }

const char*
PowerDiagramSolverSOP::inputLabel(OP_InputIdx idx) const
{
    switch (idx) {
    case 0: return "Sites (P, targets, w, u)";
    case 1: return "Clipping Domain";
    default: return "Invalid";
    }
}

bool
PowerDiagramSolverSOP::updateParmsFlags()
{
    bool changed = SOP_Node::updateParmsFlags();
    changed |= setVisibleState("plane", _MODE() == MODE_2D);
    return changed;
}

bool
PowerDiagramSolverSOP::_ReadTargetVolumes(
    const SiteInput& in,
    const char* attrib,
    double domainVolume,
    std::vector<double>& targets)
{
    GA_ROHandleF targetH(gdp, GA_ATTRIB_POINT, attrib);
    if (not targetH.isValid()) {
        addError(SOP_MESSAGE, ("No float point attribute '" + std::string(attrib) +
                               "' with the target volume fractions (e.g. append a PP Feasible "
                               "Volumes SOP).").c_str());
        return false;
    }
    targets.resize(in._offsets.size());
    for (size_t i = 0; i < targets.size(); ++i) {
        targets[i] = targetH.get(in._offsets[i]);
    }
    if (not std::all_of(targets.begin(), targets.end(), [](double f) { return f > 0.0; })) {
        addError(SOP_MESSAGE, "Target volume fractions must be positive.");
        return false;
    }

    const double sum = std::accumulate(targets.begin(), targets.end(), 0.0);
    if (std::fabs(sum - 1.0) > 1e-6) {
        addWarning(SOP_MESSAGE, ("Target volume fractions sum to " + std::to_string(sum) +
                                 "; normalizing them to 1.").c_str());
    }
    for (double& v : targets) {
        v *= domainVolume / sum;
    }
    return true;
}

DualSolver::Parameters
PowerDiagramSolverSOP::_ReadSolverParameters(fpreal t)
{
    const NewtonParms newton = ReadNewtonParameters(this, t);
    DualSolver::Parameters params;
    params._maxIterations = newton._maxIterations;
    params._regularization = newton._regularization;
    params._lineSearch._maxHalvings = newton._maxHalvings;
    params._log = newton._log;
    params._volumeTolerance = _VOLUMETOL(t);
    return params;
}

void
PowerDiagramSolverSOP::_ReportSolve(const DualSolver::Result& res)
{
    if (not res._warning.empty()) {
        addWarning(SOP_MESSAGE, res._warning.c_str());
    }
    if (not res.HasConverged()) {
        addWarning(SOP_MESSAGE, (std::string("The solver did not converge: ") +
                                 DualSolver::GetStatusString(res._status) + ".").c_str());
    }
}

bool
PowerDiagramSolverSOP::_WriteWeights(
    const SiteInput& in,
    const std::vector<double>& omega,
    const DualSolver::Result& res,
    const char* attrib,
    double domainVolume)
{
    // w = omega - |u|^2, and the achieved volume fractions.
    GA_RWHandleF outW(gdp->addFloatTuple(GA_ATTRIB_POINT, attrib, 1));
    GA_RWHandleF cellVol(gdp->addFloatTuple(GA_ATTRIB_POINT, "cell_volume", 1));
    if (not outW.isValid() or not cellVol.isValid()) {
        addError(SOP_MESSAGE, "Cannot create the output weight attribute.");
        return false;
    }
    for (size_t i = 0; i < in._offsets.size(); ++i) {
        outW.set(in._offsets[i], float(omega[i] - in._shiftSq[i]));
        cellVol.set(in._offsets[i], float(res._diagram.GetVolume(int(i)) / domainVolume));
    }
    outW.bumpDataId();
    cellVol.bumpDataId();

    SetDetailString(gdp, "solver_status", DualSolver::GetStatusString(res._status));
    SetDetailInt(gdp, "solver_iterations", res._iterations);
    SetDetailInt(gdp, "solver_converged", res.HasConverged() ? 1 : 0);
    SetDetailFloat(gdp, "solver_max_volume_error", res._volumeResidual);
    return true;
}

OP_ERROR
PowerDiagramSolverSOP::cookMySop(OP_Context& context)
{
    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT) {
        return error();
    }

    const fpreal t = context.getTime();
    duplicateSource(0, context);
    const bool is2D = (_MODE() == MODE_2D);
    const int* axes = GetPlaneAxes(is2D ? _PLANE() : PLANE_XY);

    UT_String targetName, weightName, shiftName, outWeightName;
    _TARGETATTRIB(targetName, t);
    _WEIGHTATTRIB(weightName, t);
    _SHIFTATTRIB(shiftName, t);
    _OUTWEIGHTATTRIB(outWeightName, t);

    // Sites: s = P + u, initial effective weight w + |u|^2.
    SiteInput in;
    ReadSites(gdp, is2D, axes, weightName.c_str(), shiftName.c_str(), in);
    if (not in._warning.empty()) {
        addWarning(SOP_MESSAGE, in._warning.c_str());
    }
    if (in._offsets.empty()) {
        addError(SOP_MESSAGE, "The first input has no points.");
        return error();
    }

    Domain domain;
    std::string err;
    if (not ReadDomain(inputGeo(1, context), is2D, axes, domain, err)) {
        addError(SOP_MESSAGE, err.c_str());
        return error();
    }

    std::vector<double> targets;
    if (not _ReadTargetVolumes(in, targetName.c_str(), domain.GetVolume(), targets)) {
        return error();
    }

    std::vector<double> omega = in._weights;
    const DualSolver::Result res = VolumeSolver(_ReadSolverParameters(t)).Solve(domain, in._sites, targets, omega);
    if (not res._error.empty()) {
        addError(SOP_MESSAGE, res._error.c_str());
        return error();
    }
    _ReportSolve(res);
    _WriteWeights(in, omega, res, outWeightName.c_str(), domain.GetVolume());
    return error();
}

} // namespace PP
